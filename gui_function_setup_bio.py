import copy
from PySimpleGUI import popup, popup_get_text, popup_error, popup_get_file


from database_functions import _get_list_of_names_from_database
from draw_basic import draw_plate
from bio_functions import bio_compound_info_from_worklist, bio_import_handler_single_point, bio_dead_plate_handler, \
    plate_layout_setup, set_colours
from database_functions import grab_table_data
from gui_popup import assay_generator
from gui_settings_control import GUISettingsController
from compound_plate_formatting import plate_layout_re_formate


def assay_drop_down_updates(dbf, window):
    table_data = dbf.find_column_data("assay", "assay_name")
    assay_list = [assays for assays in table_data]
    window["-BIO_ASSAY_NAME-"].update(values=assay_list)


def bio_report_update(window, values):
    if values["-BIO_COMBINED_REPORT-"]:
        window["-BIO_FINAL_REPORT_INCLUDE_HITS-"].update(disabled=False)
        window["-BIO_FINAL_REPORT_INCLUDE_SMILES-"].update(disabled=False)
    else:
        window["-BIO_FINAL_REPORT_INCLUDE_HITS-"].update(disabled=True)
        window["-BIO_FINAL_REPORT_INCLUDE_SMILES-"].update(disabled=True)
        window["-BIO_FINAL_REPORT_USE_THRESHOLD-"].update(disabled=True)
        window["-BIO_FINAL_REPORT_USE_AMOUNT-"].update(disabled=True)
        window["-BIO_FINAL_REPORT_THRESHOLD-"].update(value="")
        window["-BIO_FINAL_REPORT_THRESHOLD-"].update(disabled=True)
        window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(value="")
        window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(disabled=True)


def bio_report_hits_update(window, values):
    if values["-BIO_FINAL_REPORT_INCLUDE_HITS-"] is True:
        window["-BIO_FINAL_REPORT_USE_THRESHOLD-"].update(disabled=False)
        window["-BIO_FINAL_REPORT_USE_THRESHOLD-"].update(value=False)
        window["-BIO_FINAL_REPORT_USE_AMOUNT-"].update(disabled=False)
        window["-BIO_FINAL_REPORT_USE_AMOUNT-"].update(value=False)
    else:
        window["-BIO_FINAL_REPORT_USE_THRESHOLD-"].update(disabled=True)
        window["-BIO_FINAL_REPORT_USE_AMOUNT-"].update(disabled=True)
        window["-BIO_FINAL_REPORT_THRESHOLD-"].update(value="")
        window["-BIO_FINAL_REPORT_THRESHOLD-"].update(disabled=True)
        window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(value="")
        window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(disabled=True)


def bio_report_smiles_update(window, values):
    if values["-BIO_FINAL_REPORT_INCLUDE_SMILES-"] is True:
        window["-BIO_FINAL_REPORT_INCLUDE_STRUCTURE-"].update(disabled=True)
    else:
        window["-BIO_FINAL_REPORT_INCLUDE_STRUCTURE-"].update(disabled=False)


def bio_report_threshold_update(window, values):
    if values["-BIO_FINAL_REPORT_USE_THRESHOLD-"] is True:
        window["-BIO_FINAL_REPORT_USE_AMOUNT-"].update(value=False)
        window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(value="")
        window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(disabled=True)
        window["-BIO_FINAL_REPORT_THRESHOLD-"].update(disabled=False)
    else:
        window["-BIO_FINAL_REPORT_THRESHOLD-"].update(value="")
        window["-BIO_FINAL_REPORT_THRESHOLD-"].update(disabled=True)


def bio_report_amount_update(window, values):
    if values["-BIO_FINAL_REPORT_USE_AMOUNT-"] is True:
        window["-BIO_FINAL_REPORT_USE_THRESHOLD-"].update(value=False)
        window["-BIO_FINAL_REPORT_THRESHOLD-"].update(value="")
        window["-BIO_FINAL_REPORT_THRESHOLD-"].update(disabled=True)
        window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(disabled=False)
    else:
        window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(value="")
        window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(disabled=True)


def bio_report_compound_id_update(window, values):
    if values["-BIO_REPORT_ADD_COMPOUND_IDS-"]:
        window["-BIO_COMPOUND_DATA-"].update(value=True)


def bio_experiment_add_to_database_update(window, values):
    if values["-BIO_EXPERIMENT_ADD_TO_DATABASE-"]:
        window["-BIO_COMPOUND_DATA-"].update(value=True)


def bio_compound_data_update(window, values):
    if values["-BIO_COMPOUND_DATA-"]:
        window["-BIO_EXPERIMENT_ADD_TO_DATABASE-"].update(value=True)
    else:
        window["-BIO_REPORT_ADD_COMPOUND_IDS-"].update(value=False)


def bio_combined_report_update(window, values):
    if values["-BIO_COMBINED_REPORT-"] and not values["-FINAL_BIO_NAME-"]:

        final_report_name = popup_get_text("Final Report Name?")
        if final_report_name:
            window["-FINAL_BIO_NAME-"].update(value=final_report_name)
        else:
            window["-BIO_COMBINED_REPORT-"].update(value=False)
            window["-BIO_FINAL_REPORT_INCLUDE_HITS-"].update(disabled=True)
            window["-BIO_FINAL_REPORT_INCLUDE_SMILES-"].update(disabled=True)
            window["-BIO_FINAL_REPORT_USE_THRESHOLD-"].update(disabled=True)
            window["-BIO_FINAL_REPORT_USE_AMOUNT-"].update(disabled=True)
            window["-BIO_FINAL_REPORT_THRESHOLD-"].update(value="")
            window["-BIO_FINAL_REPORT_THRESHOLD-"].update(disabled=True)
            window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(value="")
            window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(disabled=True)


def bio_settings(config, window):
    gsc = GUISettingsController(config)
    reports = gsc.main_settings_controller()
    if reports:
        set_colours(window, reports)


def bio_plate_layout(dbf, config, window, values, well_dict):
    if values["-BIO_PLATE_LAYOUT-"]:
        well_dict.clear()
        plate_data = dbf.find_data_single_lookup("plate_layout", values["-BIO_PLATE_LAYOUT-"], "layout_name")[0]
        plate_layout = eval(plate_data[5])
        well_dict = copy.deepcopy(plate_layout)
        well_dict = plate_layout_re_formate(config, well_dict)
        plate_size = plate_data[1]
        archive = True
        gui_tab = "bio"
        graph_bio = window["-BIO_CANVAS-"]
        sample_type = values["-BIO_SAMPLE_TYPE-"]  # ToDo figure out why this is needed

        if sample_type != "Custom":
            try:
                draw_plate(config, graph_bio, plate_size, well_dict, gui_tab, archive)
            except KeyError:
                print(f"There is no place layout for {sample_type}")
        return well_dict
    else:
        return well_dict


def assay_name_update(config, window, values):
    temp_assay_data = grab_table_data(config, "assay", single_row=True, data_value=values["-BIO_ASSAY_NAME-"],
                                      headline="assay_name")
    window["-BIO_PLATE_LAYOUT-"].update(value=temp_assay_data[0][4])


def new_assay(config, dbf, window):
    plate_list = _get_list_of_names_from_database(dbf, "plate_layout", "plate_name")
    assay_check = assay_generator(config, plate_list)

    if assay_check:
        assay_table_data = grab_table_data(config, "assay")
        assay_list = []
        for row in assay_table_data[0]:
            assay_list.append(row[1])

        window["-BIO_ASSAY_NAME-"].update(values=assay_list)


def bio_finaL_report_multi_conc(window, values):
    popup("This functions does nothing ATM ")


def bio_calculate(dbf, config, values):
    bio_breaker = False
    if not values["-BIO_PLATE_LAYOUT-"]:
        popup_error("Please choose a plate layout")
    elif not values["-BIO_IMPORT_FOLDER-"]:
        popup_error("Please choose an import folder")
    elif values["-BIO_COMBINED_REPORT-"] and not values["-BIO_EXPORT_FOLDER-"]:
        popup_error("Please choose an export folder")
    elif values["-BIO_EXPORT_TO_EXCEL-"] and not values["-BIO_EXPORT_FOLDER-"]:
        popup_error("Please choose an export folder")
    elif values["-BIO_COMBINED_REPORT-"] and not values["-FINAL_BIO_NAME-"]:
        popup_error("Please choose an Report name")
    elif values["-BIO_EXPERIMENT_ADD_TO_DATABASE-"] and not values["-BIO_ASSAY_NAME-"]:
        popup_error("Please choose an Assay name")
    elif values["-BIO_EXPERIMENT_ADD_TO_DATABASE-"] and not values["-BIO_RESPONSIBLE-"]:
        popup_error("Please choose an Responsible ")
    elif values["-BIO_EXPERIMENT_ADD_TO_DATABASE-"] and not values["-BIO_FINAL_REPORT_CONCENTRATION_NUMBER-"]:
        popup_error("Please write a concentration for the run")
    elif values["-BIO_FINAL_REPORT_INCLUDE_HITS-"] and not values["-BIO_FINAL_REPORT_HIT_AMOUNT-"] \
            and values["-BIO_FINAL_REPORT_INCLUDE_HITS-"] and not values["-BIO_FINAL_REPORT_THRESHOLD-"]:
        popup_error("Please either select amount of hits or the threshold for the score, if "
                       "Hits are to be included in the report")
    # Missing setting move moving files after analyse is done.
    # elif not values["-BIO_ANALYSE_TYPE-"]:
    #     popup_error("Please choose an analyse type")
    # ToDo add guards for wrong file-formate
    else:
        # Sets values for different parametors
        bio_import_folder = values["-BIO_IMPORT_FOLDER-"]
        # default_plate_layout = archive_plates_dict[values["-BIO_PLATE_LAYOUT-"]]
        default_plate_layout = values["-BIO_PLATE_LAYOUT-"]
        include_hits = values["-BIO_FINAL_REPORT_INCLUDE_HITS-"]
        try:
            threshold = float(values["-BIO_FINAL_REPORT_THRESHOLD-"])
        except ValueError:
            threshold = values["-BIO_FINAL_REPORT_THRESHOLD-"]
        try:
            hit_amount = int(values["-BIO_FINAL_REPORT_HIT_AMOUNT-"])
        except ValueError:
            hit_amount = values["-BIO_FINAL_REPORT_HIT_AMOUNT-"]
        include_smiles = values["-BIO_FINAL_REPORT_INCLUDE_SMILES-"]
        final_report_name = values["-FINAL_BIO_NAME-"]
        export_to_excel = values["-BIO_EXPORT_TO_EXCEL-"]
        same_layout = values["-BIO_PLATE_LAYOUT_CHECK-"]
        include_structure = values["-BIO_FINAL_REPORT_INCLUDE_STRUCTURE-"]

        bio_export_folder = values["-BIO_EXPORT_FOLDER-"] # TODO make a default output folder?

        if not same_layout:
            # If there are difference between what layout each plate is using, or if you know some data needs
            # to be dismissed, you can choose different plate layout for each plate.
            plate_to_layout = plate_layout_setup(dbf, bio_import_folder, values["-BIO_PLATE_LAYOUT-"])
            if plate_to_layout is None:
                bio_breaker = True
        else:
            # If all plate uses the same plate layout
            plate_to_layout = default_plate_layout
        if not bio_breaker:
            # If this is checked, it will ask for worklist, that can be converted to a sample dict,
            # that can be used for finding sample info in the database.
            if values["-BIO_COMPOUND_DATA-"]:
                bio_sample_list = popup_get_file("Please select worklist files", multiple_files=True)

                if bio_sample_list:
                    bio_sample_list = bio_sample_list.split(";")
                    bio_sample_dict, all_destination_plates = bio_compound_info_from_worklist(config, bio_sample_list)
                else:
                    bio_sample_dict = None
                    bio_breaker = True
                    all_destination_plates = []
            else:
                all_destination_plates = None
                bio_sample_dict = None
            if not bio_breaker:

                analyse_method = values["-BIO_ANALYSE_TYPE-"]
                add_compound_ids = values["-BIO_REPORT_ADD_COMPOUND_IDS-"]
                combined_report_check = values["-BIO_COMBINED_REPORT-"]
                import_to_database_check = values["-BIO_EXPERIMENT_ADD_TO_DATABASE-"]
                responsible = values["-BIO_RESPONSIBLE-"]
                assay_name = values["-BIO_ASSAY_NAME-"]

                if analyse_method.casefold() == "single":
                    # Get concentration for the samples
                    try:
                        float(values["-BIO_FINAL_REPORT_CONCENTRATION_NUMBER-"])
                    except ValueError:
                        temp_concentration = float(popup_get_text("Please provide a concentration in uM ?? "
                                                                     "\n numbers only"))
                    else:
                        temp_concentration = float(values["-BIO_FINAL_REPORT_CONCENTRATION_NUMBER-"])

                    bio_import_handler_single_point(dbf, config, bio_import_folder, plate_to_layout, analyse_method,
                                                    bio_sample_dict, bio_export_folder, add_compound_ids,
                                                    export_to_excel, all_destination_plates,
                                                    combined_report_check, import_to_database_check,
                                                    final_report_name, include_hits,
                                                    threshold, hit_amount, include_smiles, include_structure,
                                                    assay_name, responsible, temp_concentration)

                elif analyse_method.casefold() == "dose response":
                    print(f"HEY!!! - we are doing {analyse_method}")

                popup("Done")


def bio_blank_run(config, window, values):
    # Adds data to the database for runs that have died before data have been produced.
    bio_breaker = False
    if not values["-BIO_PLATE_LAYOUT-"]:
        popup_error("Please choose a plate layout")
    elif not values["-BIO_ASSAY_NAME-"]:
        popup_error("Please choose an Assay name")
    elif not values["-BIO_RESPONSIBLE-"]:
        popup_error("Please choose an Responsible ")

    # Missing setting move moving files after analyse is done.
    # elif not values["-BIO_ANALYSE_TYPE-"]:
    #     sg.popup_error("Please choose an analyse type")
    else:
        default_plate_layout = values["-BIO_PLATE_LAYOUT-"]
        bio_sample_list = popup_get_file("Please select worklist files", multiple_files=False)
        if bio_sample_list is not None:
            worklist = bio_sample_list
        else:
            worklist = None
            bio_breaker = True

        if not bio_breaker:
            responsible = values["-BIO_RESPONSIBLE-"]
            assay_name = values["-BIO_ASSAY_NAME-"]
            analyse_method = values["-BIO_ANALYSE_TYPE-"]
            bio_dead_plate_handler(config, assay_name, worklist, analyse_method, responsible)

            popup("Done")


def send_to_info_window(config, window, values):
    popup_error("Needs to be updated to do something else")
    # ToDO this needs to be updated to fit with the new formate of bio data!
    #
    # if not values["-BIO_PLATE_LAYOUT-"]:
    #     sg.popup_error("Please choose a plate layout")
    # elif not values["-BIO_IMPORT_FOLDER-"]:
    #     sg.popup_error("Please choose an import folder")
    #
    # else:
    #     bio_import_folder = values["-BIO_IMPORT_FOLDER-"]
    #     plate_layout = archive_plates_dict[values["-BIO_PLATE_LAYOUT-"]]
    #     analyse_method = values["-BIO_ANALYSE_TYPE-"]
    #     write_to_excel = False
    #     _, all_plates_data, date = bio_data(config, bio_import_folder, plate_layout,
    #                                         bio_plate_report_setup,
    #                                         analyse_method, write_to_excel)
    #
    #     gui_tab = "bio_exp"
    #     archive = True
    #
    #     file_name = "bio_experiments.txt"
    #     # plate_dict_name = bio_exp_table_data[values["-BIO_EXP_PLATE_TABLE-"][0]][2]
    #     plate_bio_info = all_plates_data
    #
    #     bio_info_plate_layout = plate_layout
    #     bio_info_plate_size = plate_layout["plate_type"]
    #     bio_info_state_dict = copy.deepcopy(plate_layout["well_layout"])
    #     bio_info_state_dict = plate_layout_re_formate(bio_info_state_dict)
    #
    #     well_dict_bio_info, bio_info_min_x, bio_info_min_y, bio_info_max_x, bio_info_max_y \
    #         = draw_plate(config, graph_bio_exp, bio_info_plate_size, bio_info_state_dict, gui_tab, archive)
    #
    #     bio_info_plates = []
    #     bio_info_states = []
    #     bio_info_analyse_method = []
    #     bio_info_calc = []
    #     for plates in plate_bio_info:
    #         bio_info_plates.append(plates)
    #         for method in plate_bio_info[plates]["calculations"]:
    #             if method != "other":
    #                 if method not in bio_info_analyse_method:
    #                     bio_info_analyse_method.append(method)
    #                 for state in plate_bio_info[plates]["calculations"][method]:
    #                     if state not in bio_info_states:
    #                         bio_info_states.append(state)
    #                     for calc in plate_bio_info[plates]["calculations"][method][state]:
    #                         if calc not in bio_info_calc:
    #                             bio_info_calc.append(calc)
    #             if method == "other":
    #                 for calc_other in plate_bio_info[plates]["calculations"][method]:
    #                     if calc_other not in bio_info_calc:
    #                         bio_info_calc.append(calc_other)
    #
    #     # Main settings info
    #     window["-BIO_INFO_ANALYSE_METHOD-"].Update(value="original")
    #     window["-BIO_INFO_MAPPING-"].Update(value="state mapping")
    #     window["-BIO_INFO_ANALYSE_METHOD-"].update(values=bio_info_analyse_method,
    #                                                value=bio_info_analyse_method[0])
    #     window["-BIO_INFO_PLATES-"].update(values=bio_info_plates, value=bio_info_plates[0])
    #     window["-BIO_INFO_STATES-"].update(values=bio_info_states)
    #
    #     # Map settings
    #     list_box_index = None
    #     for index_state, values in enumerate(bio_info_states):
    #         if values == "sample":
    #             list_box_index = index_state
    #         if not list_box_index:
    #             list_box_index = 0
    #
    #     window["-BIO_INFO_STATE_LIST_BOX-"].update(values=bio_info_states, set_to_index=list_box_index)
    #
    #     # Matrix Settings
    #     window["-BIO_INFO_MATRIX_METHOD-"].update(values=bio_info_analyse_method)
    #     window["-BIO_INFO_MATRIX_STATE-"].update(values=bio_info_states)
    #     window["-BIO_INFO_MATRIX_CALC-"].update(values=bio_info_calc)
    #
    #     # # List settings
    #     # window["-BIO_INFO_LIST_METHOD-"].update(values=bio_info_analyse_method, value=bio_info_analyse_method[0])
    #     # window["-BIO_INFO_LIST_STATE-"].update(values=bio_info_states, value=bio_info_states[0])
    #     # window["-BIO_INFO_LIST_CALC-"].update(values=bio_info_calc, value=bio_info_calc[0])
    #
    #     # Overview settings
    #     window["-BIO_INFO_PLATE_OVERVIEW_METHOD_LIST-"].update(values=bio_info_analyse_method,
    #                                                            set_to_index=len(bio_info_analyse_method) - 1)
    #     window["-BIO_INFO_PLATE_OVERVIEW_STATE_LIST-"].update(values=bio_info_states,
    #                                                           set_to_index=len(bio_info_states) - 1)
    #     window["-BIO_INFO_PLATE_OVERVIEW_PLATE-"].update(values=bio_info_plates, value=bio_info_plates[0])
    #     window["-BIO_INFO_OVERVIEW_METHOD-"].update(values=bio_info_analyse_method,
    #                                                 value=bio_info_analyse_method[0])
    #     window["-BIO_INFO_OVERVIEW_STATE-"].update(values=bio_info_states, value=bio_info_states[0])
    #
    #     # HIT List settings
    #     window["-BIO_INFO_HIT_LIST_PLATES-"].update(values=bio_info_plates, value=bio_info_plates[0])
    #     window["-BIO_INFO_HIT_LIST_METHOD-"].update(values=bio_info_analyse_method,
    #                                                 value=bio_info_analyse_method[-1])
    #     window["-BIO_INFO_HIT_LIST_STATE-"].update(values=bio_info_states, value="sample")
    #
    #     # Popup Matrix
    #     method_values = bio_info_analyse_method
    #     state_values = bio_info_states
    #     calc_values = bio_info_calc
    #
    #     # bio_info_sub_setting_tab_mapping_calc = True
    #     # bio_info_sub_setting_tab_matrix_calc = True
    #     # bio_info_sub_setting_tab_list_calc = True
    #     bio_info_sub_setting_tab_overview_calc = True
    #     bio_info_sub_setting_tab_plate_overview_calc = True
    #     bio_info_sub_setting_tab_z_prime_calc = True
    #     bio_info_sub_setting_tab_hit_list_calc = True