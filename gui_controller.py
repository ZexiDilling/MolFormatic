import configparser
from rdkit.Chem import MolToSmiles
import multiprocessing as mp
from multiprocessing import Queue
import PySimpleGUI as sg


from config_writer import ConfigWriter
from gui_function_info_bio import update_bio_info_values, heatmap_low_colour_update, heatmap_mid_colour_update, \
    heatmap_high_colour_update, hit_low_colour_update, hit_mid_colour_update, hit_high_colour_update
from gui_function_info_calculations import calculate_dose
from gui_function_info_lcms import sample_selection_mode_update, lcms_calculation, lcms_drawing
from upstarts_values import database_guard, gui_data_fetcher, all_table_data, plate_table_table_heading_mp, \
    window_tables, window_1_lcms, window_1_plate_layout

from gui_function_setup_extra import method_do_update, add_source_wells_update, execute_button_pressed, database_tab_pressed, \
    database_responsible_import, database_customers_import, database_vendors_import, database_academia_company_import, \
    database_place_type_import, database_location_import

from gui_function_general_startup import start_up_database
from gui_function_setup_lcms import lcms_importer, lcms_info_overview, lcms_reporting
from gui_function_table_lcms import table_group_lcms, date_set_update, batch_list_box_update
from gui_function_general_menu import menu_open, menu_save, menu_about, help_info_controller, sorting_the_tables
from gui_function_setup_plate_layout import plate_layout_draw_groups, colour_target_update, plate_archive, \
    sample_type_update, dose_sample_amount, dose_dilution_replicates, dose_colouring
from gui_function_table_plate import plate_chooser_update, barcode_list_box_update, table_limiter_update, \
    table_group_tables, clear_plate_table_update
from gui_function_setup_search import search_compound, sub_search_method_update_values, \
    search_daughter_plates_update_values, search_mother_plates_update_values, search_plate_layout
from gui_function_setup_bio import bio_report_update, bio_report_hits_update, bio_report_smiles_update, \
    bio_report_threshold_update, bio_report_amount_update, bio_plate_layout, bio_compound_data_update, \
    bio_report_compound_id_update, bio_experiment_add_to_database_update, bio_combined_report_update, bio_settings, \
    bio_analyse_type, new_assay, assay_name_update, bio_fina_report_multi_conc, bio_calculate, bio_blank_run, \
    send_to_info_window
from gui_function_setup_database import update_compound, update_plates
from gui_function_setup_simulation import simulation_input_update, simulation_run
from gui_function_table_compound import tree_database_update, compound_table_refreshed, compound_table_export
from gui_function_table_bio import table_tab_group_pressed_update, experiment_table_assay_list_update, \
    experiment_table_assay_run_update, \
    experiment_table_plate_update, assay_table_double_click, compound_table_double_click, plate_table_double_click, \
    update_overview_compound, bio_exp_compound_list
from gui_function_setup_worklist import worklist_tab_clicked, worklist_control_layout_update, worklist_generator
from gui_layout import GUILayout
from gui_function_general_plate_drawing import on_up, save_layout, delete_layout, rename_layout, draw_layout, on_move, \
    bio_canvas, export_layout


def main(config, queue_gui, queue_mol):
    """
    The main control modul for the GUI.

    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :return: This is a gui control modul. Returns depending on what is being done in the GUI
    """

    # importing config writer, to write data to the config file
    global worklist_mp_plates_list
    config_writer = ConfigWriter(config)
    db_active = database_guard(config, config_writer)
    # The archive could be removed from here, but then there would be a call to the DB everytime a layout is used...
    # Grabs the updated data from the database

    plate_list, assay_list, archive_plates_dict = gui_data_fetcher(db_active, config) # TODO move archive_plates_dict to the database

    gui_layout = GUILayout(config, plate_list)

    window = gui_layout.full_layout()
    window.maximize()

    window["-BIO_ASSAY_NAME-"].update(values=assay_list)

    # temp_wave_test = True

    #   WINDOW 1 - PLATE LAYOUT #
    dose_colour_dict = {}
    well_dict = {}
    well_dict_bio_info = {}

    #   WINDOW 2 - BIO EXP  #
    graph_bio_exp = window["-BIO_INFO_CANVAS-"]
    window_1_plate_layout["graph_plate"] = window["-RECT_BIO_CANVAS-"]

    #   WINDOW 2 - PURITY INFO  #
    lc_graph_showing = [keys for keys in list(config["lc_mapping"].keys())]
    colour_select = {}
    for keys in list(config["plate_colouring"].keys()):
        colour_select[keys] = config["plate_colouring"][keys]

    # COMPOUND TABLE CONSTANTS #
    all_data = None
    compound_info_tables = {"-COMPOUND_INFO_INFO_MP-": "-COMPOUND_INFO_INFO_MP_TABLE-",
                            "-COMPOUND_INFO_INFO_DP-": "-COMPOUND_INFO_INFO_DP_TABLE-",
                            "-COMPOUND_INFO_INFO_ASSAY-": "-COMPOUND_INFO_INFO_ASSAY_TABLE-",
                            "-COMPOUND_INFO_INFO_HITS-": "-COMPOUND_INFO_INFO_HITS_TABLE-",
                            "-COMPOUND_INFO_INFO_TRANSFERS-": "-COMPOUND_INFO_INFO_TRANSFERS_TABLE-",
                            "-COMPOUND_INFO_INFO_PURITY-": "-COMPOUND_INFO_INFO_PURITY_USED_TABLE-"}


    # Makes it possible to double-click on the table
    window["-BIO_EXP_COMPOUND_TABLE-"].bind('<Double-Button-1>', "+-double click-")

    # Table stuff
    window.Element("-BIO_INFO_MATRIX_TABLE-").Widget.configure(displaycolumns=[])
    window.Element("-PLATE_TABLE_TABLE-").Widget.configure(displaycolumns=plate_table_table_heading_mp)
    search_reverse = {}

    while True:
        event, values = window.read(timeout=100)
        if event == sg.WIN_CLOSED or event == "Exit":
            queue_mol.put(("close", None))
            break

        while not queue_gui.empty():
            msg_type, data = queue_gui.get()
            if msg_type == 'smiles':
                mol = data
                smiles = MolToSmiles(mol)
                window["-SUB_SEARCH_SMILES-"].update(smiles)

        if event == "-START_UP_DB-":
            start_up_database(config, db_active, window)

        #   WINDOW MENU         ###
        if event == "Open    Ctrl-O":
            window = menu_open(config, config_writer, window, gui_layout)

        if event == "Save    Ctrl-S":
            menu_save()

        if event == "About...":
            menu_about()

        if event == "Info":
            help_info_controller(config)

        #   WINDOW 1 - SEARCH     ###
        if event == "-SEARCH_AC-":
            search_compound(window, values, config)

        if event == "-SUB_SEARCH_METHOD-":
            sub_search_method_update_values(window, values)

        if event == "-SEARCH_PLATE_PRODUCTION-" and values["-SEARCH_PLATE_PRODUCTION-"] == "Daughter Plates":
            search_daughter_plates_update_values(window)

        if event == "-SEARCH_PLATE_PRODUCTION-" and values["-SEARCH_PLATE_PRODUCTION-"] == "Mother Plates":
            search_mother_plates_update_values(window)

        if event == "-SEARCH_PLATE_LAYOUT-":
            search_plate_layout(window, values, archive_plates_dict)

        if event == "-SUB_SEARCH_DRAW_MOL-":
            queue_mol.put(("start_draw_tool", values["-SUB_SEARCH_SMILES-"]))

        #     WINDOW 1 - BIO DATA         ###
        if event == "-BIO_COMBINED_REPORT-":
            bio_report_update(window, values)

        if event == "-BIO_FINAL_REPORT_INCLUDE_HITS-":
            bio_report_hits_update(window, values)

        if event == "-BIO_FINAL_REPORT_INCLUDE_SMILES-":
            bio_report_smiles_update(window, values)

        if event == "-BIO_FINAL_REPORT_USE_THRESHOLD-":
            bio_report_threshold_update(window, values)

        if event == "-BIO_FINAL_REPORT_USE_AMOUNT-":
            bio_report_amount_update(window, values)

        if event == "-BIO_COMPOUND_DATA-":
            bio_compound_data_update(window, values)

        if event == "-BIO_REPORT_ADD_COMPOUND_IDS-":
            bio_report_compound_id_update(window, values)

        if event == "-BIO_EXPERIMENT_ADD_TO_DATABASE-":
            bio_experiment_add_to_database_update(window, values)

        if event == "-BIO_PLATE_LAYOUT-":
            well_dict = bio_plate_layout(config, window, values, well_dict, archive_plates_dict)

        if event == "-BIO_SAMPLE_TYPE-":
            well_dict = bio_plate_layout(config, window, values, well_dict, archive_plates_dict)

        if event == "-BIO_COMBINED_REPORT-":
            bio_combined_report_update(window, values)

        if event == "-BIO_REPORT_SETTINGS-" or event == "-PURITY_ADVANCED_SETTINGS-":
            bio_settings(config, window)     #TODO fix

        if event == "-BIO_ANALYSE_TYPE-":
            bio_analyse_type(window, values)

        # Add a new assay to the database
        if event == "-BIO_NEW_ASSAY-":
            new_assay(config, window, plate_list)

        if event == "-BIO_ASSAY_NAME-":
            assay_name_update(config, window, values)

        if event == "-BIO_FINAL_REPORT_CONCENTRATION_MULTIPLE-":
            bio_fina_report_multi_conc(window, values)

        if event == "-BIO_CALCULATE-":
            bio_calculate(config, window, values, plate_list, archive_plates_dict)

        if event == "-BIO_BLANK_RUN-":
            bio_blank_run(config, window, values)

        if event == "-BIO_SEND_TO_INFO-":
            send_to_info_window(config, window, values)

        #   WINDOW 1 - LCMS           ###
        if event == "-PURITY_DATA_IMPORT-":
            lcms_importer(config, window, values)

        if event == "-PURITY_INFO_PURITY_OVERVIEW_IMPORT-":
            lcms_info_overview(config, window, values)

        if event == "-PURITY_DATA_REPORT-":
            lcms_reporting(config, window, values)

        #     WINDOW 1 - PLATE LAYOUT     ###
        if event == "-PLATE_LAYOUT_DRAW_GROUPS-":
            plate_layout_draw_groups(window, values, well_dict)

        # Dose Response Tab
        if event == "-DOSE_SAMPLE_AMOUNT-" and values["-EQUAL_SPLIT-"]:
            dose_colour_dict = dose_sample_amount(window, event, values, dose_colour_dict)

        if event == "-DOSE_DILUTIONS-" or event == "-DOSE_REPLICATES-":
            dose_colour_dict = dose_dilution_replicates(window, event, values, dose_colour_dict)

        if event == "-DOSE_COLOUR_LOW-" or event == "-DOSE_COLOUR_HIGH-":
            dose_colour_dict = dose_colouring(window, event, values)

        # State mapping
        if event == "-PLATE_LAYOUT_COLOUR_CHOSE_TARGET-":
            colour_target_update(window, values)

        # Plate Layout functions:
        if event == "-ARCHIVE_PLATES-":
            plate_archive(window, values)

        if event == "-RECT_SAMPLE_TYPE-":
            sample_type_update(window, values)

        if event == "-DRAW-":
            well_dict = draw_layout(config, window, values, well_dict, archive_plates_dict)

        if event == "-EXPORT_LAYOUT-":
            export_layout(config, window, values, well_dict)

        if event == "-SAVE_LAYOUT-":
            save_layout(config, window, values, well_dict, archive_plates_dict)

        if event == "-DELETE_LAYOUT-":
            delete_layout(config, window, values)

        if event == "-RENAME_LAYOUT-":
            rename_layout(config, window, values)

        # Used both for Plate layout and Bio Info
        # prints coordinate and well under the plate layout
        try:
            event.endswith("+MOVE")

        except AttributeError:
            pass

        else:
            if event.endswith("+MOVE") and type(event) != tuple:
                on_move(window, values, graph_bio_exp, well_dict_bio_info, well_dict)

        if event == "-RECT_BIO_CANVAS-":
            bio_canvas(values)

        # it does not always detect this event:
        try:
            event.endswith("+UP")
        except AttributeError:
            pass
        else:
            if event.endswith("+UP"):
                well_dict = on_up(window, values, well_dict, dose_colour_dict, colour_select)

        #     WINDOW 1 - UPDATE Database      ###
        if event == "-UPDATE_COMPOUND-":
            window = update_compound(config, window, values, gui_layout)

        if event == "-UPDATE_MP-" or event == "-UPDATE_DP-":
            update_plates(config, window, values, event)

        if event == "-UPDATE_BIO-":
            pass
            # THIS SHOULD BE DELETED !!

        if event == "-UPDATE_AUTO-":
            pass
            # THIS SHOULD BE DELETED !!

        #     WINDOW 1 - Worklist     ###
        if event == "-TAB_GROUP_ONE-":
            temp_assay_list, worklist_mp_plates_list = worklist_tab_clicked(config, window, values)

        if event == "-WORKLIST_CONTROL_LAYOUT-":
            worklist_control_layout_update(window, values)

        if event == "-WORKLIST_USE_POSITIVE_CONTROL-":
            window["-WORKLIST_POSITIVE_CONTROL_ID-"].update(disabled=not values["-WORKLIST_USE_POSITIVE_CONTROL-"])

        if event == "-WORKLIST_USE_NEGATIVE_CONTROL-":
            window["-WORKLIST_NEGATIVE_CONTROL_ID-"].update(disabled=not values["-WORKLIST_USE_NEGATIVE_CONTROL-"])

        if event == "-WORKLIST_USE_BONUS_COMPOUND-":
            window["-WORKLIST_BONUS_COMPOUND_ID-"].update(disabled=not values["-WORKLIST_USE_BONUS_COMPOUND-"])

        if event == "-WORKLIST_GENERATE-":
            worklist_generator(config, window, values, worklist_mp_plates_list, archive_plates_dict)

        #       WINDOW 1 - EXTRA            ###
        if event == "-PD_METHOD_DD-":
            method_do_update(window, values)

        if event == "-PD_ADD_SOURCE_WELLS-":
            add_source_wells_update(window, values)

        if event == "-PD_EXECUTE_BUTTON-":
            execute_button_pressed(config, window, values, archive_plates_dict)

        if event == "-EXTRA_SUB_TABS-":
            print("WHAT IS THIS???????? ARG???????? NOT WORKING")

        if event == "-EXTRA_SUB_DATABASE_TABS-":
            database_tab_pressed(config, window, values, db_active)

        if event == "-EXTRA_DATABASE_RESPONSIBLE_IMPORT_DB-":
            database_responsible_import(config, window, values)

        if event == "-EXTRA_DATABASE_CUSTOMERS_IMPORT_DB-":
            database_customers_import(config, window, values)

        if event == "-EXTRA_DATABASE_VENDORS_IMPORT_DB-":
            database_vendors_import(config, window, values)

        if event == "-EXTRA_DATABASE_AC_IMPORT_DB-":
            database_academia_company_import(config, window, values)

        if event == "-EXTRA_DATABASE_PLACE_TYPE_IMPORT_DB-":
            database_place_type_import(config, window, values)

        if event == "-EXTRA_DATABASE_LOCATION_IMPORT_DB-":
            database_location_import(config, window, values)

        #       WINDOW 1 - SIMULATE         ###
        if event == "-SIM_INPUT_EQ-":
            simulation_input_update(window, values)

        if event == "-SIM_RUN-":
            simulation_run(window, values)

        #     WINDOW TABLES - COMPOUND TABLE      ###
        if event == "-TREE_DB-":
            compound_id = tree_database_update(config, window, values, compound_data)

        if event == "-C_TABLE_REFRESH-":
            treedata, all_data, compound_data, counter = compound_table_refreshed(config, window, values)

        if event == "-C_TABLE_EXPORT-":
            compound_table_export(config, window, values, all_data, archive_plates_dict)

        #   WINDOW TABLE - BIO EXPERIMENT TABLE     ###

        if event == "-TABLE_TAB_GRP-":
            table_tab_group_pressed_update(config, window, values)

        if event == "-BIO_EXP_TABLE_ASSAY_LIST_BOX-":
            bio_exp_assay_runs = experiment_table_assay_list_update(config, window, values)

        if event == "-BIO_EXP_ASSAY_RUN_TABLE-"  or event == "-BIO_EXP_APPROVED_PLATES_ONLY-":
            bio_exp_plate_data = experiment_table_assay_run_update(config, window, values, bio_exp_assay_runs)

        if event == "-BIO_EXP_PLATE_TABLE-" or event == "-BIO_EXP_APPROVED_COMPOUNDS_ONLY-" or \
                event == "-REFRESH_BIO_TABLE-":
            bio_exp_compound_data = experiment_table_plate_update(config, window, values, bio_exp_plate_data)

        if event == "-BIO_EXP_ASSAY_RUN_TABLE-+-double click-":
            assay_table_double_click(window, values)

        if event == "-BIO_EXP_PLATE_TABLE-+-double click-":
            plate_table_double_click(window, values)

        if event == "-BIO_EXP_COMPOUND_TABLE-+-double click-":
            compound_table_double_click(config, window, values, bio_exp_compound_data)

        if event == "-BIO_EXP_EXPORT_COMPOUNDS-" and values["-BIO_EXP_PLATE_TABLE-"]:
            bio_exp_compound_list(config, event, values, bio_exp_compound_data)

        #   WINDOW TABLE - LC EXPERIMENT    ###
        if event == "-TABLE_TAB_GRP-":
            table_group_lcms(config, window, values)

        if event == "-LC_MS_TABLE_DATE_START_TARGET-" or event == "-LC_MS_TABLE_DATE_END_TARGET-":
            date_set_update(config, window, values)

        if event == "-LC_MS_TABLE_BATCH_LIST_BOX-":
            batch_list_box_update(config, window, values, all_table_data)

        #   WINDOW TABLE - PLATE TABLE      ###
        if event == "-TABLE_TAB_GRP-":
            temp_mp_plates = table_group_tables(config, window, values)

        if event == "-PLATE_TABLE_CLEAR-":
            temp_mp_plates = clear_plate_table_update(config, window, values)

        if event == "-PLATE_TABLE_CHOOSER-" or event == "-PLATE_TABLE_START_DATE_TARGET-" or \
                event == "-PLATE_TABLE_END_DATE_TARGET-":
            plate_chooser_update(config, window, values)

        if event == "-PLATE_TABLE_BARCODE_LIST_BOX-":
            barcode_list_box_update(config, window, values)

        if event == "-PLATE_TABLE_BUTTON_LIMITER-":
            table_limiter_update(config, window, values, temp_mp_plates)

        #   WINDOW 2 - COMPOUND INFO    ###
        if event == "-COMPOUND_INFO_SEARCH_COMPOUND_ID-":
            update_overview_compound(config, window, values, None)

        if event in compound_info_tables:
            print(f"compound info table - {compound_info_tables[event]}")

        if event == "-COMPOUND_INFO_SEND_TO_SEARCH-":
            window["-SUB_SEARCH_SMILES-"].update(value=values["-COMPOUND_INFO_ID-"])

        #   WINDOW 2 - BIO INFO         ###
        if event == "-BIO_INFO_STATES-" or event == "-BIO_INFO_ANALYSE_METHOD-" or event == "-BIO_INFO_PLATES-":
            update_bio_info_values(values, window, window_tables["plate_bio_info"])

        # Updating Sub setting data
        if event == "-BIO_INFO_HEATMAP_LOW_COLOUR_TARGET-":
            heatmap_low_colour_update(window, values)

        if event == "-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-":
            heatmap_mid_colour_update(window, values)

        if event == "-BIO_INFO_HEATMAP_HIGH_COLOUR_TARGET-":
            heatmap_high_colour_update(window, values)

        if event == "-BIO_INFO_HIT_MAP_LOW_COLOUR_TARGET-":
            hit_low_colour_update(window, values)

        if event == "-BIO_INFO_HIT_MAP_MID_COLOUR_TARGET-":
            hit_mid_colour_update(window, values)

        if event == "-BIO_INFO_HIT_MAP_HIGH_COLOUR_TARGET-":
            hit_high_colour_update(window, values)

        if event == "-BIO_INFO_BOUNDS_BUTTON-":
            sg.PopupError("This is not working")
            # ToDo make this button work. Should get a small popup, to choose all the bins for the bio analysis.
        # ToDO rethink the whole BIO info stuff
        # if event == "-BIO_INFO_SUB_SETTINGS_TABS-" and values["-BIO_INFO_SUB_SETTINGS_TABS-"] == "Plate Overview" \
        #         and bio_info_sub_setting_tab_plate_overview_calc or event == "-BIO_INFO_PLATE_OVERVIEW_METHOD_LIST-" \
        #         or event == "-BIO_INFO_PLATE_OVERVIEW_STATE_LIST-" or event == "-BIO_INFO_PLATE_OVERVIEW_PLATE-":
        #     method = values["-BIO_INFO_PLATE_OVERVIEW_METHOD_LIST-"]
        #     state = values["-BIO_INFO_PLATE_OVERVIEW_STATE_LIST-"]
        #     plate = values["-BIO_INFO_PLATE_OVERVIEW_PLATE-"]
        #     all_table_data["-BIO_INFO_OVERVIEW_TABLE-"] = sub_settings_plate_overview(plate_bio_info, method, plate,
        #                                                                               state)
        #     window["-BIO_INFO_OVERVIEW_TABLE-"].update(values=all_table_data["-BIO_INFO_OVERVIEW_TABLE-"])
        #     bio_info_sub_setting_tab_plate_overview_calc = False
        #
        # if event == "-BIO_INFO_SUB_SETTINGS_TABS-" and values["-BIO_INFO_SUB_SETTINGS_TABS-"] == "Overview" \
        #         and bio_info_sub_setting_tab_overview_calc or event == "-BIO_INFO_OVERVIEW_METHOD-" \
        #         or event == "-BIO_INFO_OVERVIEW_STATE-":
        #     method = values["-BIO_INFO_OVERVIEW_METHOD-"]
        #     state = values["-BIO_INFO_OVERVIEW_STATE-"]
        #     sub_settings_overview_table_data = sub_settings_overview(plate_bio_info, method, state)
        #     all_table_data["-BIO_INFO_OVERVIEW_AVG_TABLE-"], all_table_data["-BIO_INFO_OVERVIEW_STDEV_TABLE-"], \
        #     all_table_data["-BIO_INFO_OVERVIEW_Z_PRIME_TABLE-"] = sub_settings_overview_table_data
        #     window["-BIO_INFO_OVERVIEW_AVG_TABLE-"].update(values=all_table_data["-BIO_INFO_OVERVIEW_AVG_TABLE-"])
        #     window["-BIO_INFO_OVERVIEW_STDEV_TABLE-"].update(values=all_table_data["-BIO_INFO_OVERVIEW_STDEV_TABLE-"])
        #     window["-BIO_INFO_OVERVIEW_Z_PRIME_TABLE-"].update(
        #         values=all_table_data["-BIO_INFO_OVERVIEW_Z_PRIME_TABLE-"])
        #     bio_info_sub_setting_tab_overview_calc = False
        #
        # # if event == "-BIO_INFO_SUB_SETTINGS_TABS-" and values["-BIO_INFO_SUB_SETTINGS_TABS-"] == "List" \
        # #         and bio_info_sub_setting_tab_list_calc or event == "-BIO_INFO_LIST_METHOD-" \
        # #         or event == "-BIO_INFO_LIST_STATE-" or event == "-BIO_INFO_LIST_CALC-":
        # #
        # #     method = values["-BIO_INFO_LIST_METHOD-"]
        # #     state = values["-BIO_INFO_LIST_STATE-"]
        # #     calc = values["-BIO_INFO_LIST_CALC-"]
        # #     sub_setting_list_table_data = sub_settings_list(plate_bio_info, method, state, calc)
        # #     window["-BIO_INFO_LIST_TABLE-"].update(values=sub_setting_list_table_data)
        # #     bio_info_sub_setting_tab_list_calc = False
        #
        # if event == "-BIO_INFO_SUB_SETTINGS_TABS-" and values["-BIO_INFO_SUB_SETTINGS_TABS-"] == "Z-Prime" \
        #         and bio_info_sub_setting_tab_z_prime_calc:
        #     z_prime_data = sub_settings_z_prime(plate_bio_info)
        #     all_table_data["-BIO_INFO_Z_PRIME_LIST_TABLE-"], z_prime_max_barcode, z_prime_max_value, z_prime_min_barcode \
        #         , z_prime_min_value = z_prime_data
        #     window["-BIO_INFO_Z_PRIME_LIST_TABLE-"].update(values=all_table_data["-BIO_INFO_Z_PRIME_LIST_TABLE-"])
        #     window["-BIO_INFO_Z_PRIME_MAX_BARCODE-"].update(value=z_prime_max_barcode)
        #     window["-BIO_INFO_Z_PRIME_MAX_VALUE-"].update(value=z_prime_max_value)
        #     window["-BIO_INFO_Z_PRIME_MIN_BARCODE-"].update(value=z_prime_min_barcode)
        #     window["-BIO_INFO_Z_PRIME_MIN_VALUE-"].update(value=z_prime_min_value)
        #
        #     bio_info_sub_setting_tab_z_prime_calc = False
        #
        # if event == "-BIO_INFO_SUB_SETTINGS_TABS-" and values["-BIO_INFO_SUB_SETTINGS_TABS-"] == "Hit List" \
        #         and bio_info_sub_setting_tab_hit_list_calc:
        #     plate = values["-BIO_INFO_HIT_LIST_PLATES-"]
        #     method = values["-BIO_INFO_HIT_LIST_METHOD-"]
        #     state = values["-BIO_INFO_HIT_LIST_STATE-"]
        #
        #     pora_thresholds = {
        #         "low": {"min": float(values["-BIO_INFO_PORA_LOW_MIN_HIT_THRESHOLD-"]),
        #                 "max": float(values["-BIO_INFO_PORA_LOW_MAX_HIT_THRESHOLD-"])},
        #         "mid": {"min": float(values["-BIO_INFO_PORA_MID_MIN_HIT_THRESHOLD-"]),
        #                 "max": float(values["-BIO_INFO_PORA_MID_MAX_HIT_THRESHOLD-"])},
        #         "high": {"min": float(values["-BIO_INFO_PORA_HIGH_MIN_HIT_THRESHOLD-"]),
        #                  "max": float(values["-BIO_INFO_PORA_HIGH_MAX_HIT_THRESHOLD-"])}
        #     }
        #
        #     sub_settings_hit_list_table_data = sub_settings_hit_list(plate_bio_info, plate, method, state,
        #                                                              bio_info_state_dict, pora_thresholds)
        #
        #     all_table_data["-BIO_INFO_HIT_LIST_LOW_TABLE-"], all_table_data["-BIO_INFO_HIT_LIST_MID_TABLE-"], \
        #     all_table_data["-BIO_INFO_HIT_LIST_HIGH_TABLE-"] = sub_settings_hit_list_table_data
        #     window["-BIO_INFO_HIT_LIST_LOW_TABLE-"].update(values=all_table_data["-BIO_INFO_HIT_LIST_LOW_TABLE-"])
        #     window["-BIO_INFO_HIT_LIST_MID_TABLE-"].update(values=all_table_data["-BIO_INFO_HIT_LIST_MID_TABLE-"])
        #     window["-BIO_INFO_HIT_LIST_HIGH_TABLE-"].update(values=all_table_data["-BIO_INFO_HIT_LIST_HIGH_TABLE-"])
        #
        #     bio_info_sub_setting_tab_hit_list_calc = False
        #
        # if event == "-BIO_INFO_PORA_LOW_MIN_HIT_THRESHOLD-" or event == "-BIO_INFO_PORA_LOW_MAX_HIT_THRESHOLD-" or \
        #         event == "-BIO_INFO_PORA_MID_MIN_HIT_THRESHOLD-" or event == "-BIO_INFO_PORA_MID_MAX_HIT_THRESHOLD-" or \
        #         event == "-BIO_INFO_PORA_HIGH_MIN_HIT_THRESHOLD-" or event == "-BIO_INFO_PORA_HIGH_MAX_HIT_THRESHOLD-":
        #
        #     if not bio_info_sub_setting_tab_hit_list_calc:
        #         bio_info_sub_setting_tab_hit_list_calc = True
        #
        # if event == "-BIO_INFO_MAPPING-" or event == "-BIO_INFO_ANALYSE_METHOD-" or event == "-BIO_INFO_PLATES-" \
        #         or event == "-BIO_INFO_RE_DRAW-" \
        #         or event == "-BIO_INFO_SUB_SETTINGS_TABS-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
        #         or event == "-BIO_INFO_HIT_MAP_LOW_COLOUR_TARGET-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
        #         or event == "-BIO_INFO_HIT_MAP_MID_COLOUR_TARGET-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
        #         or event == "-BIO_INFO_HIT_MAP_HIGH_COLOUR_TARGET-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
        #         or event == "-BIO_INFO_PORA_LOW_MIN_HIT_THRESHOLD-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
        #         and values["-BIO_INFO_PORA_LOW_MIN_HIT_THRESHOLD-"] \
        #         or event == "-BIO_INFO_PORA_LOW_MAX_HIT_THRESHOLD-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
        #         and values["-BIO_INFO_PORA_LOW_MAX_HIT_THRESHOLD-"] \
        #         or event == "-BIO_INFO_PORA_MID_MIN_HIT_THRESHOLD-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
        #         and values["-BIO_INFO_PORA_MID_MIN_HIT_THRESHOLD-"] \
        #         or event == "-BIO_INFO_PORA_MID_MAX_HIT_THRESHOLD-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
        #         and values["-BIO_INFO_PORA_MID_MAX_HIT_THRESHOLD-"] \
        #         or event == "-BIO_INFO_PORA_HIGH_MIN_HIT_THRESHOLD-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
        #         and values["-BIO_INFO_PORA_HIGH_MIN_HIT_THRESHOLD-"] \
        #         or event == "-BIO_INFO_PORA_HIGH_MAX_HIT_THRESHOLD-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
        #         and values["-BIO_INFO_PORA_HIGH_MAX_HIT_THRESHOLD-"] \
        #         or event == "-BIO_INFO_STATE_LIST_BOX-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
        #         or event == "-BIO_INFO_STATE_LIST_BOX-" and values["-BIO_INFO_MAPPING-"] == "Heatmap" \
        #         or event == "-BIO_INFO_HEATMAP_LOW_COLOUR_TARGET-" and values["-BIO_INFO_MAPPING-"] == "Heatmap" \
        #         or event == "-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-" and values["-BIO_INFO_MAPPING-"] == "Heatmap" \
        #         or event == "-BIO_INFO_HEATMAP_HIGH_COLOUR_TARGET-" and values["-BIO_INFO_MAPPING-"] == "Heatmap" \
        #         or event == "-BIO_INFO_HEAT_PERCENTILE_LOW-" and values["-BIO_INFO_MAPPING-"] == "Heatmap" \
        #         and values["-BIO_INFO_HEAT_PERCENTILE_LOW-"] \
        #         or event == "-BIO_INFO_HEAT_PERCENTILE_MID-" and values["-BIO_INFO_MAPPING-"] == "Heatmap" \
        #         and values["-BIO_INFO_HEAT_PERCENTILE_MID-"] \
        #         or event == "-BIO_INFO_HEAT_PERCENTILE_HIGH-" and values["-BIO_INFO_MAPPING-"] == "Heatmap" \
        #         and values["-BIO_INFO_HEAT_PERCENTILE_HIGH-"]:
        #
        #     if plate_bio_info:
        #         if values["-BIO_INFO_MAPPING-"] == "State Mapping":
        #             gui_tab = "bio_exp"
        #             archive = True
        #
        #             well_dict_bio_info, bio_info_min_x, bio_info_min_y, bio_info_max_x, bio_info_max_y, off_set \
        #                 = draw_plate(config, graph_bio_exp, bio_info_plate_size, bio_info_state_dict, gui_tab, archive)
        #
        #         if values["-BIO_INFO_MAPPING-"] == "Heatmap":
        #             mapping = {
        #                 "mapping": values["-BIO_INFO_MAPPING-"],
        #                 "colours": {"low": [values["-BIO_INFO_HEATMAP_LOW_COLOUR_TARGET-"],
        #                                     values["-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-"]],
        #                             "high": [values["-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-"],
        #                                      values["-BIO_INFO_HEATMAP_HIGH_COLOUR_TARGET-"]]},
        #                 "states": values["-BIO_INFO_STATE_LIST_BOX-"],
        #                 "percentile": {"low": float(values["-BIO_INFO_HEAT_PERCENTILE_LOW-"]),
        #                                "mid": float(values["-BIO_INFO_HEAT_PERCENTILE_MID-"]),
        #                                "high": float(values["-BIO_INFO_HEAT_PERCENTILE_HIGH-"])}
        #             }
        #
        #             gui_tab = "bio_exp"
        #             plate = values["-BIO_INFO_PLATES-"]
        #             analyse_method = values["-BIO_INFO_ANALYSE_METHOD-"]
        #
        #             temp_plate_bio_info = plate_bio_info[plate]["plates"][analyse_method]["wells"]
        #             well_dict_bio_info, bio_info_min_x, bio_info_min_y, bio_info_max_x, bio_info_max_y, off_set \
        #                 = draw_plate(config, graph_bio_exp, bio_info_plate_size, temp_plate_bio_info, gui_tab,
        #                              mapping=mapping, state_dict=bio_info_state_dict)
        #
        #         if values["-BIO_INFO_MAPPING-"] == "Hit Map":
        #             mapping = {
        #                 "mapping": values["-BIO_INFO_MAPPING-"],
        #                 "lower_bound_start": float(values["-BIO_INFO_PORA_LOW_MIN_HIT_THRESHOLD-"]),
        #                 "lower_bound_end": float(values["-BIO_INFO_PORA_LOW_MAX_HIT_THRESHOLD-"]),
        #                 "middle_bound_start": float(values["-BIO_INFO_PORA_MID_MIN_HIT_THRESHOLD-"]),
        #                 "middle_bound_end": float(values["-BIO_INFO_PORA_MID_MAX_HIT_THRESHOLD-"]),
        #                 "higher_bound_start": float(values["-BIO_INFO_PORA_HIGH_MIN_HIT_THRESHOLD-"]),
        #                 "higher_bound_end": float(values["-BIO_INFO_PORA_HIGH_MAX_HIT_THRESHOLD-"]),
        #                 "low_colour": values["-BIO_INFO_HIT_MAP_LOW_COLOUR_TARGET-"],
        #                 "mid_colour": values["-BIO_INFO_HIT_MAP_MID_COLOUR_TARGET-"],
        #                 "high_colour": values["-BIO_INFO_HIT_MAP_HIGH_COLOUR_TARGET-"],
        #                 "states": values["-BIO_INFO_STATE_LIST_BOX-"]
        #             }
        #
        #             plate = values["-BIO_INFO_PLATES-"]
        #             analyse_method = values["-BIO_INFO_ANALYSE_METHOD-"]
        #             gui_tab = "bio_exp"
        #
        #             temp_plate_bio_info = plate_bio_info[plate]["plates"][analyse_method]["wells"]
        #             well_dict_bio_info, bio_info_min_x, bio_info_min_y, bio_info_max_x, bio_info_max_y, off_set \
        #                 = draw_plate(config, graph_bio_exp, bio_info_plate_size, temp_plate_bio_info, gui_tab,
        #                              mapping=mapping, state_dict=bio_info_state_dict)
        #
        # if event == "-BIO_INFO_MATRIX_POPUP-" and plate_bio_info:
        #     method = values["-BIO_INFO_MATRIX_METHOD-"]
        #     state = values["-BIO_INFO_MATRIX_STATE-"]
        #     calc = values["-BIO_INFO_MATRIX_CALC-"]
        #
        #     matrix_popup(plate_bio_info, calc_values, state_values, method_values, calc, sub_settings_matrix, state,
        #                  method)
        #
        # if event == "-BIO_INFO_Z_PRIME_MATRIX_BUTTON-" and plate_bio_info:
        #     matrix_popup(plate_bio_info, calc_values, state_values, method_values, "z_prime", sub_settings_matrix)
        #
        # if event == "-BIO_INFO_MATRIX_BUTTON-" and plate_bio_info:
        #     data_dict = plate_bio_info
        #     state = values["-BIO_INFO_MATRIX_STATE-"]
        #     if state == "z_prime":
        #         calc = None
        #         method = None
        #     else:
        #         calc = values["-BIO_INFO_MATRIX_CALC-"]
        #         method = values["-BIO_INFO_MATRIX_METHOD-"]
        #
        #     try:
        #         all_table_data["-BIO_INFO_MATRIX_TABLE-"], display_columns = sub_settings_matrix(data_dict, calc,
        #                                                                                          method, state)
        #     except KeyError:
        #         sg.popup_error("Please select all information")
        #     else:
        #         window.Element("-BIO_INFO_MATRIX_TABLE-").Widget.configure(displaycolumns=display_columns)
        #
        #         window["-BIO_INFO_MATRIX_TABLE-"].update(values=all_table_data["-BIO_INFO_MATRIX_TABLE-"])
        #     # window["-BIO_INFO_MATRIX_TABLE-"].update(headings=headings)
        #
        # if event == "-BIO_INFO_EXPORT-":
        #     sg.popup("This functions does nothing ATM ")

        #   WINDOW 2 - PURITY INFO  ###
        if event == "-PURITY_INFO_SAMPLE_SELECTION-":
            sample_selection_mode_update(window, values)

        if event == "-PURITY_INFO_RE_CALC-":
            purity_peak_list_table_data, peak_table_data = lcms_calculation(config, window, values)

        if event == "-PURITY_INFO_DRAW_STUFF-":
            sg.Popup("DO NOT WORK ATM :D Still figuring out how to do stuff")

        if event == "-PURITY_INFO_CANVAS-":
            sg.PopupError(f"PURITY_INFO_CANVAS: {values}")

        #   WINDOW 2 - DRAWING THINGY!!!    ###
        if event == "-PURITY_INFO_SAMPLE_BOX-":
            lcms_drawing(values["-PURITY_INFO_SAMPLE_BOX-"], window, values, event, peak_table_data, lc_graph_showing,
                         purity_peak_list_table_data)

        if event == "-PURITY_INFO_GRAPH_SHOWING-":
            lcms_drawing(window_1_lcms["purity_info_samples"], window, values, event, peak_table_data, lc_graph_showing,
                         purity_peak_list_table_data)

        if event == "-PURITY_INFO_PURITY_PEAK_LIST_TABLE-" or event == "-PURITY_INFO_PURITY_OVERVIEW_TABLE-":
            lcms_drawing(window_1_lcms["purity_info_values"], window, values, event, peak_table_data, lc_graph_showing,
                         purity_peak_list_table_data)

        if event == "-PURITY_INFO_PEAK_TABLE-" or event == "-PURITY_INFO_DRAW_PEAKS-":
            lcms_drawing(window_1_lcms["plot_style"], window, values, event, peak_table_data, lc_graph_showing,
                         purity_peak_list_table_data)

        #   WINDOW 2 - CALCULATIONS DOSE  ###
        if event == "-CALCULATIONS_BUTTON_DOSE_CALCULATION-":
            calculate_dose(window, values)

        # Sorting when clicking on Table headers. All tables should be in here execpt compound table, as it is a tree
        if isinstance(event, tuple):
            sorting_the_tables(window, event, search_reverse)


if __name__ == "__main__":
    from draw_tool_handler import launch_draw_tool

    config = configparser.ConfigParser()
    config.read("config.ini")
    queue_gui = Queue()
    queue_mol = Queue()
    process_gui = mp.Process(target=main, args=(config, queue_gui, queue_mol))
    process_gui.start()

    while True:
        msg, smiles = queue_mol.get()
        print(smiles)
        if msg is None:
            break
        elif msg == "start_draw_tool":
            # Start the PySide6 process
            process_side6 = mp.Process(target=launch_draw_tool, args=(queue_gui, queue_mol, smiles))
            process_side6.start()
            process_side6.join()
        elif msg == "close":
            break


