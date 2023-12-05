from PySimpleGUI import PopupError, Popup, PopupYesNo
from natsort import natsorted

from database_controller import FetchData
from helpter_functions import sort_table
from lcms_functions import get_peak_information, import_ms_data, lcms_data_to_db, add_start_end_time, lcms_ops, \
    grab_sample_data, lcms_to_compounds, name_changer
from gui_guards import guard_purity_data_wavelength
from gui_popup import sample_to_compound_name_controller, ms_raw_name_guard
from start_up_values import all_table_data, window_1_lcms, ms_mode_selector


def lcms_importer(config, window, values):
    if window_1_lcms["purit_info_values"]:
        # clearing purity info window:
        window["-PURITY_INFO_RT_START-"].update(value="")
        window["-PURITY_INFO_RT_END-"].update(value="")
        window["-PURITY_INFO_WAVELENGTH-"].update(value="")
        window["-PURITY_INFO_BIN-"].update(value="")
        window["-PURITY_INFO_MZ-"].update(value="")
        window["-PURITY_INFO_BATCH_BOX-"].update(values=[])
        window["-PURITY_INFO_SAMPLE_BOX-"].update(values=[])
        window["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"].update(values="")
        window["-PURITY_INFO_PEAK_TABLE-"].update(values="")
        window["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"].update(values="")
        all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"] = None
        all_table_data["-PURITY_INFO_PEAK_TABLE-"] = None
        all_table_data["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"] = None
        window["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"].update(select_rows=[])
        window["-PURITY_INFO_PEAK_TABLE-"].update(select_rows=[])
        window["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"].update(select_rows=[])
        window_1_lcms["purity_data_added_to_db"] = False
        window_1_lcms["purity_data"] = None

        window["-PURITY_DATA_IMPORT-"].update(text="Import Data")
        window_1_lcms["purity_info_values"] = False

    elif not window_1_lcms["purit_info_values"]:  # ToDo add threading
        temp_wavelength = values["-PURITY_DATA_UV_WAVE-"]
        temp_wave_test, wavelength_data = guard_purity_data_wavelength(temp_wavelength)

        if not values["-PURITY_DATA_IMPORT_FOLDER-"]:
            PopupError("Please select a folder")
        # Check for excel file if needed
        elif values["-PURITY_DATA_CALC_PURITY-"] and not values["-PURITY_DATA_COMPOUND_DATA-"] and not values[
            "-PURITY_DATA_USE_COMPOUNDS-"]:
            PopupError("Please select an import file")

        elif not temp_wave_test:
            PopupError(wavelength_data)

        else:

            folder = values["-PURITY_DATA_IMPORT_FOLDER-"]
            all_data = import_ms_data(folder)

            # GAURD# Checking if there are UV data
            if isinstance(all_data, str):
                Popup("Missing UV Data")
            else:
                _, purity_samples, window_1_lcms["purity_data"], missing_samples = all_data
                if missing_samples:
                    guard = PopupYesNo(
                        f"Missing following data: {missing_samples}. Do you want to continue (for now this will stop the process!!!)")  # ToDo Fix this, make it possible to continue
                else:
                    compound_data = values["-PURITY_DATA_USE_COMPOUNDS-"]
                    slope_threshold = int(values["-PURITY_DATA_SLOPE_THRESHOLD-"])
                    uv_threshold = int(values["-PURITY_DATA_UV_THRESHOLD-"])
                    rt_solvent_peak = float(values["-PURITY_DATA_RT_SOLVENT-"])
                    sample_data = values["-PURITY_DATA_COMPOUND_DATA-"]
                    wavelength_data = values["-PURITY_DATA_UV_WAVE-"]

                    # Make sure that the names are correct for the compound data. Will change the name in purity-data
                    if compound_data:
                        fd = FetchData(config)
                        # ToDo duplicate code later on, for checking data that is not compound data
                        new_names = sample_to_compound_name_controller(config, window_1_lcms["purity_data"], fd,
                                                                       window_1_lcms["purity_sample_layout_import"],
                                                                       window_1_lcms["purity_sample_layout_export"],
                                                                       sort_table)

                        if new_names:
                            for sample in new_names:
                                if new_names[sample] == "Delete":
                                    window_1_lcms["purity_data"].pop(sample)
                                else:
                                    window_1_lcms["purity_data"][new_names[sample]] = window_1_lcms["purity_data"].pop(sample)
                        purity_samples = []
                        for samples in window_1_lcms["purity_data"]:
                            purity_samples.append(samples)

                    if values["-PURITY_DATA_ADD_TO_DATABASE-"]:
                        _ = lcms_data_to_db(config, window_1_lcms["purity_data"])
                        window_1_lcms["purity_data_added_to_db"] = True

                    peak_information, peak_table_data, sample_peak_dict = get_peak_information(
                        window_1_lcms["purity_data"], slope_threshold, uv_threshold, rt_solvent_peak, sample_data,
                        wavelength_data)

                    if values["-PURITY_DATA_CALC_PURITY-"]:
                        if not values["-PURITY_DATA_USE_COMPOUNDS-"]:
                            window_1_lcms["sample_data_file"] = values["-PURITY_DATA_COMPOUND_DATA-"]
                        else:
                            window_1_lcms["sample_data_file"] = "compound_data"

                        ms_mode = ms_mode_selector[values["-PURITY_DATA_MS_MODE-"]]
                        delta_mass = float(values["-PURITY_DATA_MS_DELTA-"])
                        mz_threshold = int(values["-PURITY_DATA_MS_THRESHOLD-"])
                        peak_amounts = int(values["-PURITY_DATA_MS_PEAKS-"])

                        sample_data, db_data = grab_sample_data(window_1_lcms["sample_data_file"],
                                                                window_1_lcms["purity_data"], config)

                        if not compound_data:
                            # GUARD Check if file names the same. #ToDO this is duplicated earlier. FIX ! ! !
                            raw_data_samples = natsorted([keys for keys in window_1_lcms["purity_data"]])
                            excel_data_samples = natsorted([keys for keys in sample_data])

                            if raw_data_samples != excel_data_samples:
                                # If they are not the same, a popup will show, where you can set names for the data.
                                new_names = ms_raw_name_guard(raw_data_samples, excel_data_samples, db_data,
                                                              config)

                            name_changer(new_names, window_1_lcms["purity_data"], sample_data, peak_information,
                                         sample_peak_dict)

                        all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"], \
                        purity_peak_list_table_data = lcms_ops(sample_data, window_1_lcms["purity_data"],
                                                               peak_information, ms_mode, delta_mass, mz_threshold,
                                                               peak_amounts)

                        add_start_end_time(purity_peak_list_table_data, sample_peak_dict)
                        window["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"]. \
                            update(values=all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"])

                    window["-PURITY_INFO_SAMPLE_BOX-"].update(values=purity_samples)
                    window["-PURITY_DATA_IMPORT-"].update(text="Clear Purity Info")
                    window_1_lcms["purity_info_values"] = True
                    # window["-PURITY_INFO_OVERVIEW_TABLE-"].update(values=all_table_data["-PURITY_INFO_OVERVIEW_TABLE-"])


def lcms_info_overview(config, window, values):
    if window_1_lcms["purity_data"]:  # ToDo Check if this works
        if not values["-PURITY_DATA_USE_COMPOUNDS-"]:
            PopupError("Not compound data. Can't be added to the database. ")
        else:
            if not window_1_lcms["purity_data_added_to_db"]:
                batch_dict = lcms_data_to_db(config, window_1_lcms["purity_data"])

            lcms_to_compounds(config, all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"])


def lcms_reporting(config, window, values):
    Popup("Not working atm - should create a report over the data")
