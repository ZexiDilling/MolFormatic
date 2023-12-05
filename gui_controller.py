import configparser
from rdkit.Chem import MolToSmiles
import multiprocessing as mp
from multiprocessing import Queue
import PySimpleGUI as sg

from config_writer import ConfigWriter
from database_handler import DataBaseFunctions
from gui_function_info_bio import colour_chooser_update, bio_info_window_update, \
    bio_info_plate_list_update, bio_info_plate_update
from gui_function_info_calculations import calculate_dose
from gui_function_info_lcms import sample_selection_mode_update, lcms_calculation, lcms_drawing
from gui_popup import popup_table
from start_up_values import database_guard, all_table_data, compound_info_tables, window_1_lcms, \
    start_up_gui, search_reverse, colour_chooser_buttons, bio_info_tables, assay_updater_list
from gui_function_setup_extra import method_do_update, add_source_wells_update, execute_button_pressed, \
    database_tab_pressed, \
    database_responsible_import, database_customers_import, database_vendors_import, database_academia_company_import, \
    database_place_type_import, database_location_import
from gui_function_general_startup import start_up_database
from gui_function_setup_lcms import lcms_importer, lcms_info_overview, lcms_reporting
from gui_function_table_lcms import table_group_lcms, date_set_update, batch_list_box_update
from gui_function_general_menu import menu_open, menu_save, menu_about, help_info_controller, sorting_the_tables
from gui_function_setup_plate_layout import plate_layout_draw_groups, colour_target_update, plate_archive, \
    plate_list_updater, dose_sample_amount, dose_dilution_replicates, dose_colouring
from gui_function_table_plate import plate_chooser_update, barcode_list_box_update, table_limiter_update, \
    table_group_tables, clear_plate_table_update
from gui_function_setup_search import search_compound, sub_search_method_update_values, \
    search_daughter_plates_update_values, search_mother_plates_update_values, search_sample_counter_update
from gui_function_setup_bio import bio_report_update, bio_report_hits_update, bio_report_smiles_update, \
    bio_report_threshold_update, bio_report_amount_update, bio_plate_layout, bio_compound_data_update, \
    bio_report_compound_id_update, bio_experiment_add_to_database_update, bio_combined_report_update, bio_settings, \
    new_assay, assay_name_update, bio_finaL_report_multi_conc, bio_calculate, bio_blank_run, \
    send_to_info_window, assay_drop_down_updates
from gui_function_setup_database import update_compound, update_plates
from gui_function_setup_simulation import simulation_input_update, simulation_run
from gui_function_table_compound import tree_database_update, compound_table_refreshed, compound_table_export
from gui_function_table_bio import table_tab_group_pressed_update, experiment_table_assay_list_update, \
    experiment_table_assay_run_update, \
    experiment_table_plate_update, compound_table_double_click, \
    update_overview_compound, bio_exp_compound_list, bio_tables_double_clicked
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
    config_writer = ConfigWriter(config)
    db_active = database_guard(config, config_writer)
    dbf = DataBaseFunctions(config)
    gui_layout = GUILayout(config, dbf)
    window = gui_layout.full_layout()
    window.maximize()

    well_dict_bio_info, well_dict, dose_colour_dict, colour_select, graph_bio_exp, lc_graph_showing = start_up_gui(config, window)

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
            search_sample_counter_update(dbf, window, values)

        if event == "-SUB_SEARCH_DRAW_MOL-":
            queue_mol.put(("start_draw_tool", values["-SUB_SEARCH_SMILES-"]))

        #     WINDOW 1 - BIO DATA         ###
        if event == "-TAB_GROUP_ONE-" and values["-TAB_GROUP_ONE-"] == "Bio Data":
            assay_drop_down_updates(dbf, window)

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
            well_dict = bio_plate_layout(dbf, config, window, values, well_dict)

        if event == "-BIO_COMBINED_REPORT-":
            bio_combined_report_update(window, values)

        if event == "-BIO_REPORT_SETTINGS-" or event == "-PURITY_ADVANCED_SETTINGS-":
            bio_settings(config, window)     #TODO fix ??

        if event == "-BIO_ANALYSE_TYPE-":
            print("Needs to change plate-layout notes depending on the choose... ")
            plate_list_updater(dbf, window, values, event)

        if event == "-BIO_NEW_ASSAY-":
            new_assay(config, dbf, window)

        if event == "-BIO_ASSAY_NAME-":
            assay_name_update(config, window, values)

        if event == "-BIO_FINAL_REPORT_CONCENTRATION_MULTIPLE-":
            bio_finaL_report_multi_conc(window, values)

        if event == "-BIO_CALCULATE-":
            bio_calculate(dbf, config, values)

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
        if event == "-RECT_SAMPLE_TYPE-":
            plate_list_updater(dbf, window, values, event)

        if event == "-DRAW-":
            well_dict = draw_layout(dbf, config, window, values, well_dict)

        if event == "-EXPORT_LAYOUT-":
            export_layout(config, window, values, well_dict)

        if event == "-SAVE_LAYOUT-":
            save_layout(dbf, config, window, values, well_dict)

        if event == "-DELETE_LAYOUT-":
            delete_layout(dbf, window, values)

        if event == "-RENAME_LAYOUT-":
            rename_layout(dbf, window, values)

        # Used both for Plate layout and Bio Info
        # prints coordinate and well under the plate layout
        try:
            event.endswith("+MOVE")

        except AttributeError:
            pass

        else:
            if event.endswith("+MOVE") and type(event) != tuple:
                on_move(window, values, graph_bio_exp, well_dict, well_dict_bio_info)

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

        #     WINDOW 1 - Worklist     ###
        if event == "-TAB_GROUP_ONE-" and values["-TAB_GROUP_ONE-"] == "Worklist":
            temp_assay_list, worklist_mp_plates_list = worklist_tab_clicked(config, window)

        if event == "-WORKLIST_CONTROL_LAYOUT-":
            worklist_control_layout_update(window, values)

        if event == "-WORKLIST_USE_POSITIVE_CONTROL-":
            window["-WORKLIST_POSITIVE_CONTROL_ID-"].update(disabled=not values["-WORKLIST_USE_POSITIVE_CONTROL-"])

        if event == "-WORKLIST_USE_NEGATIVE_CONTROL-":
            window["-WORKLIST_NEGATIVE_CONTROL_ID-"].update(disabled=not values["-WORKLIST_USE_NEGATIVE_CONTROL-"])

        if event == "-WORKLIST_USE_BONUS_COMPOUND-":
            window["-WORKLIST_BONUS_COMPOUND_ID-"].update(disabled=not values["-WORKLIST_USE_BONUS_COMPOUND-"])

        if event == "-WORKLIST_GENERATE-":
            worklist_generator(dbf, config, window, values, worklist_mp_plates_list)

        #       WINDOW 1 - EXTRA            ###
        if event == "-PD_METHOD_DD-":
            method_do_update(window, values)

        if event == "-PD_ADD_SOURCE_WELLS-":
            add_source_wells_update(window, values)

        if event == "-PD_EXECUTE_BUTTON-":
            execute_button_pressed(dbf, config, window, values)

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
            compound_table_export(dbf, config, window, values, all_data)

        #   WINDOW TABLE - BIO EXPERIMENT TABLE     ###
        if event == "-TABLE_TAB_GRP-":
            table_tab_group_pressed_update(config, window, values)

        if event == "-BIO_EXP_TABLE_ASSAY_LIST_BOX-":
            experiment_table_assay_list_update(config, window, values)

        if event == "-BIO_EXP_ASSAY_RUN_TABLE-" or event == "-BIO_EXP_APPROVED_PLATES_ONLY-":
            experiment_table_assay_run_update(config, window, values)

        if event == "-BIO_EXP_PLATE_TABLE-" or event == "-BIO_EXP_APPROVED_COMPOUNDS_ONLY-" or \
                event == "-REFRESH_BIO_TABLE-":
            experiment_table_plate_update(config, window, values)

        if event in bio_info_tables:
            bio_tables_double_clicked(window, values, event)

        if event == "-BIO_EXP_COMPOUND_TABLE-+-double click-":
            compound_table_double_click(config, window, values, event)

        if event == "-BIO_EXP_EXPORT_COMPOUNDS-" and values["-BIO_EXP_PLATE_TABLE-"]:
            bio_exp_compound_list(config, event, values)

        #   WINDOW TABLE - LC EXPERIMENT    ###
        if event == "-TABLE_TAB_GRP-":
            table_group_lcms(config, window, values)

        if event == "-LC_MS_TABLE_DATE_START_TARGET-" or event == "-LC_MS_TABLE_DATE_END_TARGET-":
            date_set_update(config, window, values)

        if event == "-LC_MS_TABLE_BATCH_LIST_BOX-":
            batch_list_box_update(config, window, values, all_table_data)

        if event == "-LC_MS_SAMPLE_TABLE-+-double click-":
            print(f"DOUBLE CLICKING LC SAMPLE TABLE THINGY DATYA !!!!! - {event}")

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

        if event == "-PLATE_TABLE_UPDATE_VOLUME-":
            sg.PopupError("Missing function - should update the volume of the selected plate's compounds based on ECHO survey file")

        if event == "-PLATE_TABLE_TABLE-+-double click-":
            compound_table_double_click(config, window, values, event)

        #   WINDOW 2 - COMPOUND INFO    ###
        if event == "-COMPOUND_INFO_SEARCH_COMPOUND_ID-":
            update_overview_compound(config, window, values, None)

        if event in compound_info_tables:
            print(f"compound info table - {compound_info_tables[event]}")

        if event == "-COMPOUND_INFO_SEND_TO_SEARCH-":
            window["-SUB_SEARCH_SMILES-"].update(value=values["-COMPOUND_INFO_ID-"])

        if event == "-COMPOUND_INFO_MP-+-double click-":
            popup_table("-COMPOUND_INFO_INFO_MP_TABLE-")

        if event == "-COMPOUND_INFO_DP-+-double click-":
            popup_table("-COMPOUND_INFO_INFO_DP_TABLE-")

        if event == "-COMPOUND_INFO_ASSAY-+-double click-":
            popup_table("-COMPOUND_INFO_INFO_ASSAY_TABLE-")

        if event == "-COMPOUND_INFO_HITS-+-double click-":
            popup_table("-COMPOUND_INFO_INFO_HITS_TABLE-")

        if event == "-COMPOUND_INFO_PURITY-+-double click-":
            popup_table("-COMPOUND_INFO_INFO_PURITY_USED_TABLE-")

        #   WINDOW 2 - BIO INFO         ###
        # Updating Sub setting data
        if event in colour_chooser_buttons:
            colour_chooser_update(window, values, event)

        if event in assay_updater_list:
            bio_info_window_update(dbf, window, values)

        if event == "-BIO_INFO_RUN_DROPDOWN-":
            bio_info_plate_list_update(dbf, window, values, None)

        if event == "-BIO_INFO_PLATES_DROPDOWN-":
            well_dict_bio_info = bio_info_plate_update(dbf, config, window, values, event, well_dict_bio_info)

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


