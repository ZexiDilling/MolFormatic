from os.path import isfile
from pathlib import Path

import PySimpleGUI as sg

from bio_data_functions import txt_to_xlsx, original_data_dict
from bio_date_handler import BIOAnalyser
from bio_dose_controller import dose_response_controller
from bio_dose_excle_handler import dose_excel_controller
from bio_report_setup import bio_final_report_controller
from csv_handler import CSVReader
from database_functions import get_number_of_rows, grab_table_data
from database_handler import DataBaseFunctions
from draw_basic import draw_plate
from extra_functions import unit_converter
from file_handler import get_file_list
from helpter_functions import sort_table
from gui_popup import new_headlines_popup, assay_run_naming, bio_data_approval_table, bio_dose_response_set_up, \
    dead_run_naming, plate_layout_chooser
from info import matrix_header


def bio_compound_info_from_worklist(config, bio_sample_list):
    """
    Gets a dictionary, with each Destination plate in a worklist as a key, and each well in the plate have a dictionary
    with the source plate and source well. It is used to look up compounds in the database
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param sg: module
    :param bio_sample_list: A list of files
    :type bio_sample_list: list
    :return: a dict over destiantion plates with source plate and source well for each well in the destination plate
    :rtype dict
    """

    # right_headlines = ["source_plates", "destination_plates", "source_well", "destination_well"]
    right_headlines_v1 = [headlines for headlines in config["worklist_headlines_v1"]]
    temp_headlines = [headlines for headlines in config["worklist_headlines_v2"]]
    right_headlines_v2 = []
    for headline in temp_headlines:
        right_headlines_v2.append(config["worklist_headlines_v2"][headline])
    all_destination_plates = []
    sample_dict = {}
    for files in bio_sample_list:
        file_headlines = CSVReader.grab_headlines(files)
        if file_headlines[0] in right_headlines_v1:
            right_headlines = right_headlines_v1
            config_headline = "worklist_headlines_v1"
        elif file_headlines[0] in right_headlines_v2:
            right_headlines = right_headlines_v2
            config_headline = "worklist_headlines_v2"
        else:
            right_headlines = right_headlines_v1
            config_headline = "worklist_headlines_v1"

        new_headline = None
        msg, file_headlines, sample_dict, all_destination_plates = CSVReader.echo_worklist_to_dict(
            config, config_headline, files, right_headlines, new_headline, sample_dict, all_destination_plates)

        # if the file uses wrong names for the headlines, this will give a popup with the "wrong headlines" and an
        # option to change them.
        if msg == "Wrong headlines":
            print(f"{files} have the wrong format")
            new_headline = new_headlines_popup(sort_table, right_headlines, file_headlines)
            msg, file_headlines, sample_dict, all_destination_plates = CSVReader.echo_worklist_to_dict(
                config, config_headline, files, right_headlines, new_headline, sample_dict, all_destination_plates)
        elif msg == "Not a CSV file":
            sg.popup_error("Wrong file fomate!!!")              # ToDo sort out this guard so it dose not crash.
            continue

    return sample_dict, all_destination_plates


def bio_data(config, folder, plate_to_layout, archive_plates_dict, analysis_method,
             bio_sample_dict, save_location, add_compound_ids, write_to_excel=True):

    """
    Handles the Bio data.

    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param folder: The folder where the raw data is located
    :type folder: str
    :param plate_to_layout: Either the name of the default plate_layout that all plates are using or a dict over each
        plate and the layout for that place
    :type plate_to_layout: str or dict
    :param archive_plates_dict: The layout for the plate with values for each well, what state they are in
    :type archive_plates_dict: dict
    :param bio_sample_dict: None or a dict of sample ide, per plate analysed
    :type bio_sample_dict: dict
    :param analysis_method: The analysis method
    :type analysis_method: str
    :param save_location: where to save all the excel files
    :type save_location: str
    :param add_compound_ids: Will add the compound ID to each well on the hit-list
    :type add_compound_ids: bool
    :param write_to_excel: A bool statement to see if the results should be exported to excel
    :type write_to_excel: bool
    :return: All the data for the plates raw data, and their calculations
    :rtype: dict
    """
    # needs to reformat plate-layout to use well ID instead of numbers...
    bioa = BIOAnalyser(config)
    file_list = get_file_list(folder)
    all_plates_data = {}
    used_plates = []

    plate_layout_dict = {}
    file_amount = len(file_list)
    for file_index, files in enumerate(file_list):
        if isfile(files) and files.endswith(".txt"):
            files = txt_to_xlsx(files)
        if isfile(files) and files.endswith(".xlsx"):
            temp_barcode = files.split("/")[-1].split(".")[0]
            if type(plate_to_layout) is dict:
                plate_layout_dict = plate_to_layout
                temp_plate_layout_name = plate_to_layout[temp_barcode]
            else:
                plate_layout_dict[temp_barcode] = plate_to_layout
                temp_plate_layout_name = plate_to_layout
            if temp_plate_layout_name == "skip":
                continue
            else:
                temp_plate_layout = archive_plates_dict[temp_plate_layout_name]

                all_data, well_row_col, well_type, barcode, date = original_data_dict(files, temp_plate_layout)
                if not all_data:
                    return False
                used_plates.append(barcode)
                all_plates_data[barcode] = bioa.bio_data_controller(files, temp_plate_layout, all_data, well_row_col,
                                                                    well_type, analysis_method, write_to_excel,
                                                                    bio_sample_dict, save_location, add_compound_ids)

        else:
            print(f"{files} is not the right formate")
        print(f"file: {file_index + 1} / {file_amount}")
    return True, all_plates_data, date, used_plates, plate_layout_dict


def bio_full_report(config, analyse_method, all_plate_data, output_folder,
                    final_report_name, include_hits, threshold, hit_amount, include_smiles, bio_sample_dict,
                    plate_to_layout, archive_plates_dict, include_structure):
    """
    Writes the final report for the bio data
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param analyse_method: The analysed method used for the data
    :type analyse_method: str
    :param all_plate_data: All the data for all the plates, raw and calculations
    :type all_plate_data: dict
    :param output_folder: The output folder, where the final report ends up
    :type output_folder: str
    :param final_report_name: The name for the report
    :type final_report_name: str
    :param include_hits: If the final reports should include hits
    :type include_hits: bool
    :param threshold: If the report should include hits, then this is the threshold for what is consideret as hits
        It can either use threshold or hit_amount
    :type threshold: int or None
    :param hit_amount: If the report should include hits, then how many hits should be included
        It can either use threshold or hit_amount
    :type hit_amount: int or None
    :param include_smiles: If there should be id and smiles on the hits
    :type include_smiles: bool
    :param bio_sample_dict: A dict for ID's and smiles for each well on each plate
    :type bio_sample_dict: dict or None
    :param plate_to_layout: a dicts for the plate_layout
    :type plate_to_layout: dict
    :param archive_plates_dict: the dict over the layouys
    :type archive_plates_dict: dict
    :param include_structure: boolen to determen if the file report should include png of the structure
    :type include_structure: bool
    :return: A excel report file with all the data
    """

    output_file = f"{output_folder}/{final_report_name}.xlsx"
    bio_final_report_controller(config, analyse_method, all_plate_data, output_file,
                                include_hits, threshold, hit_amount, include_smiles, bio_sample_dict, plate_to_layout,
                                archive_plates_dict, include_structure)


def bio_experiment_to_database(config, assay_data, used_plates, all_plates_data, plate_analyse_methods, responsible,
                               bio_sample_dict, concentration, all_compound_data, sub_layouts, dismissed_plates,
                               dead_plates):
    """
    Controls adding Biological data, to two different databases: Biological_plate_data and Biological_compound_data

    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param assay_data: Data for the assay
    :type assay_data: dict
    :param all_plates_data: The data for all the plates
    :type all_plates_data: dict
    :param plate_analyse_methods: The method that have been used to analyse the data
    :type plate_analyse_methods: list
    :param responsible: The responsible person for the run
    :type responsible: dict or str
    :param bio_sample_dict: All the data for the compounds
    :type bio_sample_dict: dict
    :param concentration: Either a single value if the analys type is "single". It should be a dict for Dose-Reponse
    :type concentration: str or dict
    :param all_compound_data:
    :type all_compound_data: dict
    :param sub_layouts:
    :type sub_layouts: dict
    :return:
    """

    # set-up for the import:
    plate_table = "biological_plate_data"
    compound_table = "biological_compound_data"
    # Set-up plate counter to update the assay table
    compound_counter = 0
    plate_counter = 0

    # connects with the database to be able to check data up against it, to make sure that data is not duplicated
    dbf = DataBaseFunctions(config)

    # Add plates to the database
    total_plate_amount = len(all_plates_data)
    sub_layout_table = "plate_layout_sub"

    # add new sub_layouts to db:
    # for plates in sub_layouts:
    #     if sub_layouts[plates]["new"]:
    #         temp_plate_sub = sub_layouts[plates]["name"]
    #         temp_plate_layout = sub_layouts[plates]["layout"]
    #         temp_sub_layout_data = {
    #             "plate_sub": temp_plate_sub,
    #             "plate_main": assay_data["plate_layout"],
    #             "plate_layout": f"{temp_plate_layout['well_layout']}"
    #         }
    #         dbf.add_records_controller(sub_layout_table, temp_sub_layout_data)

    for plate_index, plates in enumerate(used_plates):

        # A guard to check if the plate is already in the database. If this is the case it skips to the next plate
        guard = dbf.find_data_single_lookup(plate_table, plates, "plate_name")
        if guard:
            print(plates)
            print("skipping plate")
            continue

        # Set up for the import of each plate:
        exp_id = get_number_of_rows(config, plate_table) + 1
        if plates in dismissed_plates:
            print(f"dissed: {all_plates_data[plates]}")
            plate_raw_data = f"{all_plates_data[plates]['plates'][plate_analyse_methods[0]]['wells']}"
            process_data = "Null"
            z_prime = float(all_plates_data[plates]["calculations"]["other"]["z_prime"])
            plate_approval = f"{all_plates_data[plates]['approved']}"
            plate_note = f"{all_plates_data[plates]['note']}"
            assay_run = f"{all_plates_data[plates]['run_name']}"
            analysis_method = f"{all_plates_data[plates]['analysis_method']}"
            skipped_wells = f"{all_plates_data[plates]['skipped_wells']}"
        elif plates in dead_plates:
            print(f"dead_plate: {all_plates_data[plates]}")
            plate_raw_data = "None"
            process_data = "Null"
            z_prime = "Null"
            plate_approval = f"{all_plates_data[plates]['approved']}"
            plate_note = "Dead Plate - See Run note"
            assay_run = f"{all_plates_data[plates]['run_name']}"
            analysis_method = f"{all_plates_data[plates]['analysis_method']}"
            skipped_wells = f"{all_plates_data[plates]['skipped_wells']}"
        else:
            print(f"Good: {all_plates_data[plates]}")
            plate_raw_data = f"{all_plates_data[plates]['plates'][plate_analyse_methods[0]]['wells']}"
            process_data = f"{all_plates_data[plates]['plates'][plate_analyse_methods[-1]]['wells']}"
            z_prime = float(all_plates_data[plates]["calculations"]["other"]["z_prime"])
            plate_approval = f"{all_plates_data[plates]['approved']}"
            plate_note = f"{all_plates_data[plates]['note']}"
            assay_run = f"{all_plates_data[plates]['run_name']}"
            analysis_method = f"{all_plates_data[plates]['analysis_method']}"
            skipped_wells = f"{all_plates_data[plates]['skipped_wells']}"

        if sub_layouts[plates]["new"]:
            temp_plate_layout = sub_layouts[plates]["name"]
        else:
            temp_plate_layout = assay_data["assay_name"]

        temp_plate_data = {
            "exp_id": exp_id,
            "assay_run": assay_run,
            "plate_name": plates,
            "raw_data": plate_raw_data,
            "process_data": process_data,
            "z_prime": z_prime,
            "responsible": responsible,
            "approval": plate_approval,
            "note": plate_note,
            "plate_layout": temp_plate_layout,
            "skipped_wells": skipped_wells,
            "analysis_method": analysis_method
        }

        # Adds the plate to the database
        dbf.add_records_controller(plate_table, temp_plate_data)

        plate_counter += 1

        # Loops through the wells on each plate.
        for wells in bio_sample_dict[plates]:

            # Gets the bio_data_id
            bio_data_id = get_number_of_rows(config, compound_table) + 1

            compound_id = bio_sample_dict[plates][wells]["compound_id"]

            if not compound_id:
                table = "compound_mp"
                well_headline = "mp_well"
                temp_well_data = bio_sample_dict[plates][wells]["source_well"]
                plate_headline = "mp_barcode"
                temp_plate_data = bio_sample_dict[plates][wells]["source_plate"]
                temp_row_data = dbf.find_data_double_lookup(table, temp_well_data, temp_plate_data, well_headline,
                                                            plate_headline)

                try:
                    temp_row_data[0][3]
                except IndexError:
                    print(f"Missing data for MP: {temp_plate_data} well: {temp_well_data}")
                else:
                    compound_id = temp_row_data[0][3]
                    bio_sample_dict[plates][wells]["compound_id"] = compound_id
                    try:
                        score = all_plates_data[plates]["plates"][plate_analyse_methods[-1]]["wells"][wells]
                    except KeyError:
                        score = ""
                    #     hit = False
                    # else:

                        # if score <= assay_data["hit_threshold"]:
                        #     hit = True
                        # else:
                        #     hit = False
                    # Hit should always be False. The Hits need to be chosen manually.
                    hit = False

                    # The raw data is the initial value before any calculations have been done.
                    try:
                        compound_raw_data = all_plates_data[plates]["plates"][plate_analyse_methods[0]]["wells"][wells]
                    except KeyError:
                        compound_raw_data = 0

                    # gets the concentration
                    if type(concentration) == str or type(concentration) == float:
                        temp_concentration = concentration
                    else:
                        temp_concentration = concentration[wells]
                    # check if the well have been transferred before adding it to the database.
                    if wells in all_plates_data[plates]["transferred_wells"]:
                        # Not sure if this is needed...
                        compound_note = all_compound_data[plates][wells]["note"]
                        compound_approval = all_compound_data[plates][wells]["approved"]
                        transferred = all_compound_data[plates][wells]["transfereed"]

                    else:
                        compound_note = "Not transferred"
                        compound_approval = "0"
                        transferred = "0"

                    temp_compound_data = {
                        "bio_data_id": bio_data_id,
                        "compound_id": compound_id,
                        "assay_plate": plates,
                        "assay_well": wells,
                        "score": score,
                        "hit": hit,
                        "concentration": temp_concentration,
                        "raw_data": compound_raw_data,
                        "approval": compound_approval,
                        "note": compound_note,
                        "transferred": transferred
                    }

                    dbf.add_records_controller(compound_table, temp_compound_data)
                    compound_counter += 1

        print(f"{plate_index + 1} / {total_plate_amount} have been uploaded to the database - last plate was: {plates}")


def bio_import_handler_single_point(config, bio_import_folder, plate_to_layout, archive_plates_dict,
                                    analyse_method, bio_sample_dict, bio_export_folder, add_compound_ids, export_to_excel,
                                    all_destination_plates, combined_report_check, import_to_database_check,
                                    final_report_name, include_hits, threshold, hit_amount,
                                    include_smiles, include_structure, assay_name, responsible, concentration):

    worked, all_plates_data, date, used_plates, plate_to_layout = \
        bio_data(config, bio_import_folder, plate_to_layout, archive_plates_dict,
                 analyse_method, bio_sample_dict, bio_export_folder, add_compound_ids, export_to_excel)

    # Check if there should be produced a combined report over all the data
    if combined_report_check:
        bio_full_report(config, analyse_method, all_plates_data,
                        bio_export_folder, final_report_name, include_hits,
                        threshold, hit_amount, include_smiles, bio_sample_dict, plate_to_layout,
                        archive_plates_dict, include_structure)

    # Check if the data should be added to the database
    if import_to_database_check:

        # Makes a popup to assign assay_run names to plates.
        if len(all_destination_plates) > len(used_plates):
            dead_run_check = sg.PopupYesNo("Not all plates from worklist have data."
                                           "Do you wish to add plates without data to the database?")
        else:
            dead_run_check = "No"

        # Adds dead plates to the all_plates_data dict for later use
        if dead_run_check.casefold() == "yes":
            dead_plates = _add_dead_plates(all_destination_plates, used_plates, all_plates_data, bio_sample_dict)
            used_plates = all_destination_plates
        else:
            dead_plates = []

        check, transfer_dict, dismissed_plates = assay_run_naming(config, all_plates_data, analyse_method,
                                                                  used_plates, assay_name, all_destination_plates,
                                                                  dead_run_check, bio_compound_info_from_worklist)
        if check:
            if len(dismissed_plates) == len(used_plates):
                all_plates_are_dismissed = True
            else:
                all_plates_are_dismissed = False
            # Set up values for the database.

            # # Grabs the value for the assay from the database
            assay_data_row = grab_table_data(config, table_name="assay", single_row=True,
                                             data_value=assay_name, headline="assay_name")
            plate_layout = assay_data_row[0][4]
            plate_data = grab_table_data(config, table_name="plate_layout", single_row=True, data_value=plate_layout,
                                         headline="plate_name")
            plate_size = plate_data[0][2]
            assay_data = {"assay_name": assay_data_row[0][2],
                          "plate_layout": plate_layout,
                          "z_prime_threshold": assay_data_row[0][5],
                          "hit_threshold": assay_data_row[0][6],
                          "plate_size": plate_size}

            # Open a popup window where data can be checked before being added to the database.
            all_plates_data, all_compound_data, plate_analyse_methods, sub_layouts = \
                bio_data_approval_table(draw_plate, config, all_plates_data, assay_data, plate_to_layout,
                                        archive_plates_dict, transfer_dict, dismissed_plates, all_plates_are_dismissed,
                                        dead_plates)

            if type(all_plates_data) != str:
                # Adds the approved data to the database
                bio_experiment_to_database(config, assay_data, used_plates, all_plates_data, plate_analyse_methods,
                                           responsible, bio_sample_dict, concentration, all_compound_data, sub_layouts,
                                           dismissed_plates, dead_plates)

    return


def dose_response_import_handler(config, assay_name, plate_reader_files, worklist, plate_layout, add_to_database,
                                 excel_report, save_location,  include_id, include_pic, dose_out="uM"):
    """
    The Dose response module for taking care of importing dose response data to the database.
    For now the reporting is also placed here... should maybe be moved somewhere else
    :param config:
    :param plate_reader_files: A list of the file name of the plate reader data
    :type plate_reader_files: list
    :param worklist: The worklist for the echo
    :type worklist: str
    :param plate_layout: The layout for the plate.
    :type plate_layout: dict
    :param add_to_database: If the data should be added to the database or not
    :type add_to_database: bool
    :param excel_report: If there should be created a report over the data
    :type excel_report: bool
    :param save_location: If there is a report created, where to save it
    :type save_location: str
    :param include_id: If The report should include the ID's of the compounds
    :type include_id: bool
    :param include_pic: If the report should include pictures of the compound
    :type include_pic: bool
    :param dose_out: What unit the Dose response should be in
    :type dose_out: str
    :return:
    """

    check = bio_dose_response_set_up(config, worklist, assay_name, plate_reader_files, bio_compound_info_from_worklist)

    if type(check) == str:
        return check
    else:
        all_data, dose_response_curve_shape, method_calc_reading_50, vol_needed_pure = check

        if include_id:
            well_to_id_dict = _dose_response_destination_plate_well_to_compound_id(worklist)
            plate_group_to_compound_id = _dose_response_destination_well_to_group(well_to_id_dict, plate_layout)
        else:
            plate_group_to_compound_id = None

        dose_layout = _dose_repsonse_list(vol_needed_pure)
        all_dose_readings = {}
        file_list = []
        for file in plate_reader_files:
            file_list.append(file)
            # Get raw data:
            single_plate_data, well_row_col, well_type, barcode, date = original_data_dict(file, plate_layout)

            temp_dose_readings = {}
            all_dose_readings[barcode] = _dose_response_data_to_dict(temp_dose_readings, plate_layout, single_plate_data,
                                                                     barcode, dose_layout, dose_out)

            all_dose_readings[barcode]["state_data"] = _dose_response_get_state_data(single_plate_data, well_type)

        all_dose_data = {}
        for plates in all_dose_readings:

            all_dose_data[plates] = dose_response_controller(config, dose_response_curve_shape,
                                                             all_dose_readings[plates], method_calc_reading_50)
        if excel_report:
            dose_excel_controller(config, file_list, all_dose_data, plate_group_to_compound_id, save_location,
                                  include_id, include_pic)

        if add_to_database:

            _add_dose_response_data_to_the_database(all_data)


def _add_dose_response_data_to_the_database(all_data):

    print(all_data)


def _bio_add_dead_plates_to_db(config, dbf, assay_data, assay_run_name, plate_list, analysis_method, responsible,
                               plate_table="biological_plate_data"):
    for plate_index, plates in enumerate(plate_list):

        # A guard to check if the plate is already in the database. If this is the case it skips to the next plate
        guard = dbf.find_data_single_lookup(plate_table, plates, "plate_name")
        if guard:
            print(plates)
            print("skipping plate")
            continue

        # Set up for the import of each plate:
        exp_id = get_number_of_rows(config, plate_table) + 1
        plate_raw_data = "None"
        process_data = "Null"
        z_prime = "Null"
        plate_approval = f"False"
        plate_note = "Dead Plate - See Run note"
        assay_run = assay_run_name
        analysis_method = analysis_method
        skipped_wells = f"All"
        temp_plate_layout = assay_data["plate_layout"]

        temp_plate_data = {
            "exp_id": exp_id,
            "assay_run": assay_run,
            "plate_name": plates,
            "raw_data": plate_raw_data,
            "process_data": process_data,
            "z_prime": z_prime,
            "responsible": responsible,
            "approval": plate_approval,
            "note": plate_note,
            "plate_layout": temp_plate_layout,
            "skipped_wells": skipped_wells,
            "analysis_method": analysis_method
        }

        # Adds the plate to the database
        dbf.add_records_controller(plate_table, temp_plate_data)


def bio_dead_plate_handler(config, assay_name, worklist, analysis_method, responsible):

    # Grab list of plates from the CSV file
    csv_reader = CSVReader
    plate_list = csv_reader.echo_grab_plate_names(worklist)

    # Popup to set-up the data for the plates. Date for the run, notes and so on.
    check, assay_run_name = dead_run_naming(config, assay_name, plate_list, worklist, bio_compound_info_from_worklist)

    if check:
        # # Grabs the value for the assay from the database
        assay_data_row = grab_table_data(config, table_name="assay", single_row=True,
                                         data_value=assay_name, headline="assay_name")
        plate_layout = assay_data_row[0][4]
        plate_data = grab_table_data(config, table_name="plate_layout", single_row=True, data_value=plate_layout,
                                     headline="plate_name")
        plate_size = plate_data[0][2]
        assay_data = {"assay_name": assay_data_row[0][2],
                      "plate_layout": plate_layout,
                      "z_prime_threshold": assay_data_row[0][5],
                      "hit_threshold": assay_data_row[0][6],
                      "plate_size": plate_size}

        dbf = DataBaseFunctions(config)

        _bio_add_dead_plates_to_db(config, dbf, assay_data, assay_run_name, plate_list, analysis_method, responsible)


def sub_settings_matrix(data_dict, calc, method, state):
    values = {}
    temp_plate_name = ["", ""]
    for plates in data_dict:
        temp_plate_name.append(plates)
        if state:
            values[plates] = data_dict[plates]["calculations"][method][state][calc]
        else:
            values[plates] = data_dict[plates]["calculations"]["other"][calc]

    # Writes name in headers
    table_data = [temp_plate_name]

    # Writes the first row of data
    temp_row_data = []
    for index, plates in enumerate(values):
        if index == 0:
            temp_row_data.append("")
            temp_row_data.append("")
        temp_row_data.append(values[plates])
    table_data.append(temp_row_data)

    # Writes the rest of the Matrix
    for index_row, plate_row in enumerate(values):
        temp_row_data = []
        for index_col, plate_col in enumerate(values):
            #  writes plate names in column 1
            if index_col == 0:
                temp_row_data.append(plate_row)
                temp_row_data.append(values[plate_row])
            try:
                temp_value = (values[plate_col] / values[plate_row]) * 100
                temp_value = round(temp_value, 2)
            except ZeroDivisionError:
                temp_value = ""
            temp_row_data.append(temp_value)

        table_data.append(temp_row_data)

    display_columns = []
    for x in range(len(table_data)):
        display_columns.append(matrix_header[x])

    return table_data, display_columns


def sub_settings_list(data_dict, method, state, calc):

    row_data = []
    for plates in data_dict:
        temp_row_data = []
        if calc != "z_prime":
            temp_row_data.append(plates)
            temp_row_data.append(data_dict[plates]["calculations"][method][state][calc])
        else:
            temp_row_data.append(plates)
            temp_row_data.append(data_dict[plates]["calculations"]["other"][calc])
        row_data.append(temp_row_data)

    return row_data


def sub_settings_plate_overview(data_dict, method, plate, state):
    row_data = []
    for methods in data_dict[plate]["calculations"]:
        include_method = True
        if methods in method:
            for states in data_dict[plate]["calculations"][methods]:
                include_state = True
                if states in state:
                    for calc in data_dict[plate]["calculations"][methods][states]:
                        if include_method:
                            temp_row_data = [methods]
                        else:
                            temp_row_data = [[]]

                        if include_state:
                            temp_row_data.append(states)
                        else:
                            temp_row_data.append([])
                        temp_row_data.append(calc)
                        temp_row_data.append(data_dict[plate]["calculations"][methods][states][calc])
                        row_data.append(temp_row_data)
                        include_method = False
                        include_state = False

    for calc in data_dict[plate]["calculations"]["other"]:
        temp_row_data = [[], calc, [], data_dict[plate]["calculations"]["other"][calc]]
        row_data.append(temp_row_data)

    return row_data


def sub_settings_overview(data_dict, method, state):
    row_data_avg = []
    row_data_stdev = []

    for plates in data_dict:
        for calc in data_dict[plates]["calculations"][method][state]:
            if calc == "avg":
                temp_row_avg = [plates, data_dict[plates]["calculations"][method][state][calc]]
                row_data_avg.append(temp_row_avg)
            elif calc == "stdev":
                temp_row_stdev = [plates, data_dict[plates]["calculations"][method][state][calc]]
                row_data_stdev.append(temp_row_stdev)

    row_data_z_prime, _, _ = listing_z_prime(data_dict)

    return row_data_avg, row_data_stdev, row_data_z_prime


def listing_z_prime(data_dict):
    row_data_z_prime = []
    z_prime_dict = {}
    z_prime_values = []
    for plates in data_dict:
        for calc in data_dict[plates]["calculations"]["other"]:
            temp_row_z_prime = [plates, data_dict[plates]["calculations"]["other"][calc]]
            z_prime_dict[data_dict[plates]["calculations"]["other"][calc]] = plates
            z_prime_values.append(data_dict[plates]["calculations"]["other"][calc])

        row_data_z_prime.append(temp_row_z_prime)

    return row_data_z_prime, z_prime_dict, z_prime_values


def sub_settings_z_prime(data_dict):
    row_data_z_prime, z_prime_dict, z_prime_values = listing_z_prime(data_dict)

    z_prime_max_value = max(z_prime_values)
    z_prime_min_value = min(z_prime_values)

    z_prime_max_barcode = z_prime_dict[z_prime_max_value]
    z_prime_min_barcode = z_prime_dict[z_prime_min_value]

    return row_data_z_prime, z_prime_max_barcode, z_prime_max_value, z_prime_min_barcode, z_prime_min_value


def sub_settings_hit_list(data_dict, plate, method, state, state_dict, pora_thresholds):

    row_data_low = []
    row_data_mid = []
    row_data_high = []

    for well in data_dict[plate]["plates"][method]["wells"]:
        if isinstance(well, str):
            if state == state_dict[well]["state"]:
                well_value = data_dict[plate]["plates"][method]["wells"][well]
                if pora_thresholds["low"]["min"] < well_value < pora_thresholds["low"]["max"]:
                    temp_row_data = [well, well_value]
                    row_data_low.append(temp_row_data)
                if pora_thresholds["mid"]["min"] < well_value < pora_thresholds["mid"]["max"]:
                    temp_row_data = [well, well_value]
                    row_data_mid.append(temp_row_data)
                if pora_thresholds["high"]["min"] < well_value < pora_thresholds["high"]["max"]:
                    temp_row_data = [well, well_value]
                    row_data_high.append(temp_row_data)

    return row_data_low, row_data_mid, row_data_high


def _add_dead_plates(all_destination_plates, used_plates, all_plates_data, bio_sample_dict):
    """
    Adds dead plates to "all_plates_data" to make it available for adding to the database
    :param all_destination_plates:
    :param used_plates:
    :param all_plates_data:
    :param bio_sample_dict:
    :return:
    """
    dead_plates = []
    for plates in all_destination_plates:
        if plates not in used_plates:
            all_plates_data[plates] = {}
            all_plates_data[plates]["transferes"] = []
            dead_plates.append(plates)

    return dead_plates


def _dose_response_destination_plate_well_to_compound_id(worklist):
    """
    Gets compound id for each well in the destination plate from a worklist in csv format
    :param file:
    :return:
    """

    well_to_id_dict = {}
    with open(worklist) as f:
        for row_index, line in enumerate(f):
            values = line.split(";")
            for value_index, value in enumerate(values):
                value = value.strip()
                if row_index == 0:
                    if value == "destination_plates":
                        destination_plate_index = value_index
                    elif value == "destination_well":
                        destination_well_index = value_index
                    elif value == "compound_id":
                        compound_id_index = value_index

                else:
                    if value_index == destination_plate_index:
                        temp_destination_plate = value
                    elif value_index == destination_well_index:
                        temp_destination_well = value
                    elif value_index == compound_id_index:
                        temp_compound_id = value

            if row_index != 0 and "filler" not in temp_compound_id:
                try:
                    well_to_id_dict[temp_destination_plate]
                except KeyError:
                    well_to_id_dict[temp_destination_plate] = {temp_destination_well: temp_compound_id}
                else:
                    well_to_id_dict[temp_destination_plate][temp_destination_well] = temp_compound_id

    return well_to_id_dict


def _dose_response_destination_well_to_group(well_to_id_dict, plate_layout):
    plate_group_to_compound_id = {}
    for plates in well_to_id_dict:
        for counter in plate_layout:
            if plate_layout[counter]["group"] != 0:
                temp_name = f"{plates}_{plate_layout[counter]['group']}_{plate_layout[counter]['replicate']}"
                temp_well = plate_layout[counter]["well_id"]
                try:
                    well_to_id_dict[plates][temp_well]
                except KeyError:
                    continue
                else:
                    plate_group_to_compound_id[temp_name] = well_to_id_dict[plates][temp_well]

    return plate_group_to_compound_id


def _dose_repsonse_list(vol_needed_pure):

    dose_layout = []
    for counter in vol_needed_pure:
        dose_layout.append(vol_needed_pure[counter]["new_conc"])

    return dose_layout


def _dose_response_data_to_dict(temp_dose_data, plate_layout, all_data, barcode, dose_layout, dose_out):

    for counter in plate_layout:
        if plate_layout[counter]["group"] != 0:
            temp_name = f"{barcode}_{plate_layout[counter]['group']}_{plate_layout[counter]['replicate']}"
            temp_well = plate_layout[counter]["well_id"]
            well_value = all_data["plates"]["original"]["wells"][temp_well]
            temp_concentration = plate_layout[counter]["concentration"] - 1
            temp_dose = dose_layout[temp_concentration]
            temp_dose = round(float(unit_converter(temp_dose, old_unit_out=False, new_unit_out=dose_out,
                                                   as_list=True)[0]), 4)
            try:
                temp_dose_data[temp_name]
            except KeyError:
                temp_dose_data[temp_name] = {"reading": {"raw": [well_value]},
                                             "dose": {"raw": [temp_dose], "unit": dose_out}}
            else:
                temp_dose_data[temp_name]["reading"]["raw"].append(well_value)
                temp_dose_data[temp_name]["dose"]["raw"].append(temp_dose)

    return temp_dose_data


def _dose_response_get_state_data(single_plate_data, well_type):
    temp_state_data = {}
    for states in well_type:
        if states != "sample" and states != "empty":
            for wells in well_type[states]:

                temp_well_value = single_plate_data["plates"]["original"]["wells"][wells]
                try:
                    temp_state_data[states]
                except KeyError:
                    temp_state_data[states] = [temp_well_value]
                else:
                    temp_state_data[states].append(temp_well_value)

    return temp_state_data


def _folder_digger_for_file_names(folder):
    """
    List all files in a folder and all sup folders as there name alone:
    :param folder: A folder
    :type folder: str
    :return: A list of file names without suffix
    :rtype: list
    """
    all_files = Path(folder).glob("**/*")
    files = [file for file in all_files if file.is_file()]
    file_list = []
    for file in files:
        temp_file = file.name.split(".")[0]
        file_list.append(temp_file)

    return file_list


def plate_layout_setup(folder, default_plate_layout, plate_layout_list):
    """
    Gets all files from a folder, and list them in a table, where the user can choose the layout for each plate
    :param folder:
    :type folder: str
    :param default_plate_layout: The default plate layout
    :type default_plate_layout: str
    :param plate_layout_list: A list with all the plate layouts
    :type plate_layout_list: list
    :return: a dicts with each plate having its own layout
    :rtype: dict
    """

    files = _folder_digger_for_file_names(folder)
    all_plate_layouts = plate_layout_chooser(files, default_plate_layout, plate_layout_list)

    return all_plate_layouts


def set_colours(window, reports):
    """
    Update all the input colour fields with new colours, after changes in the settings.
    :param window: The sg window
    :type window: PySimpleGUI.PySimpleGUI.Window
    :param reports: All the settings from the menu
    :type reports: dict
    :return:
    """
    _, bio_plate_report_setup, ms, simple_settings = reports

    window["-PLATE_LAYOUT_COLOUR_BOX_SAMPLE-"].\
        update(background_color=simple_settings["plate_colouring"]["sample"])
    window["-PLATE_LAYOUT_COLOUR_BOX_BLANK-"].\
        update(background_color=simple_settings["plate_colouring"]["blank"])
    window["-PLATE_LAYOUT_COLOUR_BOX_NAX-"].\
        update(background_color=simple_settings["plate_colouring"]["max"])
    window["-PLATE_LAYOUT_COLOUR_BOX_MINIMUM-"].\
        update(background_color=simple_settings["plate_colouring"]["minimum"])
    window["-PLATE_LAYOUT_COLOUR_BOX_POSITIVE-"].\
        update(background_color=simple_settings["plate_colouring"]["positive"])
    window["-PLATE_LAYOUT_COLOUR_BOX_NEGATIVE-"].\
        update(background_color=simple_settings["plate_colouring"]["negative"])
    window["-PLATE_LAYOUT_COLOUR_BOX_EMPTY-"].\
        update(background_color=simple_settings["plate_colouring"]["empty"])
    window["-BIO_PLATE_LAYOUT_COLOUR_BOX_SAMPLE-"].\
        update(background_color=simple_settings["plate_colouring"]["sample"])
    window["-BIO_PLATE_LAYOUT_COLOUR_BOX_BLANK-"].\
        update(background_color=simple_settings["plate_colouring"]["blank"])
    window["-BIO_PLATE_LAYOUT_COLOUR_BOX_NAX-"].\
        update(background_color=simple_settings["plate_colouring"]["max"])
    window["-BIO_PLATE_LAYOUT_COLOUR_BOX_MINIMUM-"].\
        update(background_color=simple_settings["plate_colouring"]["minimum"])
    window["-BIO_PLATE_LAYOUT_COLOUR_BOX_POSITIVE-"].\
        update(background_color=simple_settings["plate_colouring"]["positive"])
    window["-BIO_PLATE_LAYOUT_COLOUR_BOX_NEGATIVE-"].\
        update(background_color=simple_settings["plate_colouring"]["negative"])
    window["-BIO_PLATE_LAYOUT_COLOUR_BOX_EMPTY-"].\
        update(background_color=simple_settings["plate_colouring"]["empty"])

    # window["-MP_DP_UPDATE_PLATE_LAYOUT_COLOUR_BOX_SAMPLE-"].\
    #     update(background_color=simple_settings["plate_colouring"]["sample"])
    # window["-MP_DP_UPDATE_PLATE_LAYOUT_COLOUR_BOX_BLANK-"].\
    #     update(background_color=simple_settings["plate_colouring"]["blank"])
    # window["-MP_DP_UPDATE_PLATE_LAYOUT_COLOUR_BOX_NAX-"].\
    #     update(background_color=simple_settings["plate_colouring"]["max"])
    # window["-MP_DP_UPDATE_PLATE_LAYOUT_COLOUR_BOX_MINIMUM-"].\
    #     update(background_color=simple_settings["plate_colouring"]["minimum"])
    # window["-MP_DP_UPDATE_PLATE_LAYOUT_COLOUR_BOX_POSITIVE-"].\
    #     update(background_color=simple_settings["plate_colouring"]["positive"])
    # window["-MP_DP_UPDATE_PLATE_LAYOUT_COLOUR_BOX_NEGATIVE-"].\
    #     update(background_color=simple_settings["plate_colouring"]["negative"])
    # window["-MP_DP_UPDATE_PLATE_LAYOUT_COLOUR_BOX_EMPTY-"].\
    #     update(background_color=simple_settings["plate_colouring"]["empty"])

    window["-BIO_INFO_HEATMAP_LOW_COLOUR_BOX-"].\
        update(background_color=bio_plate_report_setup["heatmap_colours"]["low"])
    window["-BIO_INFO_HEATMAP_MID_COLOUR_BOX-"].\
        update(background_color=bio_plate_report_setup["heatmap_colours"]["mid"])
    window["-BIO_INFO_HEATMAP_HIGH_COLOUR_BOX-"].\
        update(background_color=bio_plate_report_setup["heatmap_colours"]["high"])
    window["-BIO_INFO_HIT_MAP_TH_1_COLOUR_BOX-"].\
        update(background_color=bio_plate_report_setup["pora_threshold"]["colour"]["th_1"])
    window["-BIO_INFO_HIT_MAP_TH_2_COLOUR_BOX-"].\
        update(background_color=bio_plate_report_setup["pora_threshold"]["colour"]["th_2"])
    window["-BIO_INFO_HIT_MAP_TH_3_COLOUR_BOX-"].\
        update(background_color=bio_plate_report_setup["pora_threshold"]["colour"]["th_3"])