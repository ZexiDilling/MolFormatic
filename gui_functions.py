from pathlib import Path

import PySimpleGUI as sg
import configparser
from datetime import date
from operator import itemgetter
from os.path import isfile, isdir

import natsort
from natsort import natsorted
from os import mkdir

from gui_popup import new_headlines_popup, plate_layout_chooser
from csv_handler import CSVWriter, CSVConverter, CSVReader
from lc_data_handler import LCMSHandler
from database_handler import DataBaseFunctions
from database_controller import FetchData, AddData
from file_handler import get_file_list
from bio_data_functions import original_data_dict, well_row_col_type, txt_to_xlsx
from bio_report_setup import bio_final_report_controller
from bio_date_handler import BIOAnalyser
from info import matrix_header, plate_384_row
from json_handler import dict_writer, dict_reader
from pickle_handler import df_writer
from heatmap import Heatmap
from config_writer import ConfigWriter
from plate_formatting import plate_layout_to_well_ditc, daughter_plate_generator, plate_layout_re_formate
from excel_handler import export_plate_layout, plate_dilution_write_vol_well_amount, plate_dilution_excel, \
    well_compound_list
from data_miner import dm_controller
from plate_dilution import PlateDilution
from visualization import *
from lcms_uv_integration import Integration
from lcms_ms_search import mass_search
from chem_operators import ChemOperators


def config_update(config):
    fd = FetchData(config)
    cw = ConfigWriter(config)
    # database_specific_commercial
    search_limiter = {
        "ac": {"value": "Commercial",
               "operator": "=",
               "target_column": "ac",
               "use": True}}
    rows = fd.data_search("origin", search_limiter)
    simple_settings = {"database_specific_commercial": {},
                       "database_specific_academic": {}}
    for row in rows:
        simple_settings["database_specific_commercial"][f"vendor_{rows[row]['ac_id']}"] = rows[row]["origin"]

    # database_specific_academia
    search_limiter = {
        "ac": {"value": "Academic",
               "operator": "=",
               "target_column": "ac",
               "use": True}}
    rows = fd.data_search("origin", search_limiter)
    for row in rows:
        simple_settings["database_specific_commercial"][f"academia_{rows[row]['ac_id']}"] = rows[row]["origin"]

    cw.run(simple_settings, "simple_settings", True)


def sort_table(table, cols, reverse):
    """ sort a table by multiple columns
        table: a list of lists (or tuple of tuples) where each inner list
               represents a row
        cols:  a list (or tuple) specifying the column numbers to sort by
               e.g. (1,0) would sort by column 1, then by column 0
    """
    for col in reversed(cols):
        try:
            table = natsorted(table, key=itemgetter(col), reverse=reverse)
        except Exception as e:
            sg.popup_error('Error in sort_table', 'Exception in sort_table', e)
    reverse = not reverse
    return table, reverse


def amount_samples(plate_amount, samples_per_plate=384):
    """
    calculate the amount of samples needed depending on amount of mother plates.

    :param plate_amount: Amount of plates needed
    :type plate_amount: int
    :param samples_per_plate: Amount of samples per plate
    :type samples_per_plate: int
    :return: Total amount of samples needed
    :rtype: int
    """
    if not isinstance(samples_per_plate, int):
        return None
    if not samples_per_plate > 0:
        return None

    if isinstance(plate_amount, int):
        if plate_amount < 1:
            return None
        return plate_amount * samples_per_plate
    else:
        return None


def _compound_list(config, mp_amount, min_mp, samples_per_plate, ignore_active, sub_search, smiles,
                   sub_search_methode, threshold, source_table, fd, search_limiter):
    """
    Generate list of compounds, based on number of motherplates only, or by sub_structure search.
    Generate comPOUND file for fecthing tubes from the comPOUND freezer.

    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param mp_amount: Amount of samples
    :type mp_amount: int
    :param samples_per_plate: amount of samples per plate
    :type samples_per_plate: int
    :param ignore_active: If the list needs to take into account compounds already in MotherPlates with more than
        1 uL volume left. - This might needs to be lowered...
    :type ignore_active: bool
    :param smiles: smiles code to compare compounds with for sub search
    :type smiles: str
    :param sub_search: true or false, If it will use structure search to find compounds or not
    :type sub_search: bool
    :param sub_search_methode: What method to use for making a substructure search
    :type sub_search_methode: str
    :param threshold: This is for sub_searchs. how alike the compounds should minimum be
    :type threshold: float
    :param source_table: The table from the database, where the samples are coming from. and the table where the
        structure search is used. should always be compound_main
    :type source_table: str
    :param search_limiter: A dict over values to search for in the db
    :type search_limiter: dict
    :return:
        - compound_list: A list of compounds
        - liquid_warning_list: A warning for compounds that are close to zero
    :rtype:
        - list
        - list
    """
    if mp_amount:
        sample_amount = amount_samples(mp_amount, samples_per_plate)
    else:
        sample_amount = None
    plated_compounds = []
    if not ignore_active:
        # plated_compounds = [compounds for compounds in fd.data_search(config["Tables"]["compound_mp_table"], None)]
        plated_compounds = [[test for row, test in data.items() if row == "compound_id"][0] for _, data in
                            fd.data_search(config["Tables"]["compound_mp_table"], None).items()]

    # Gets a list of compounds, based on search criteria
    items = fd.list_limiter(sample_amount, min_mp, samples_per_plate, source_table, sub_search, sub_search_methode, smiles,
                            threshold, ignore_active, plated_compounds, search_limiter)

    return items


def table_update_tree(mp_amount, min_mp, samples_per_plate, ignore_active, sub_search, smiles, sub_search_methode,
                      threshold, source_table, search_limiter, config):
    """
    Updates the compound table with compounds depending on search criteria

    :param mp_amount: amount of mother plates to find compounds from
    :type mp_amount: int
    :param min_mp: If all the samples needs  to be from as few plates as possible
    :type min_mp: bool
    :param samples_per_plate: amount of sample per plate
    :type samples_per_plate: int
    :param ignore_active: If it should take into account witch compounds are allready in MotherPlates.
        True = All compounds, False = only compounds not found in MotherPlates
    :type ignore_active: bool
    :param sub_search: Is it uses structure search for finding the compounds:
    :type sub_search: bool
    :param smiles: smiles code for the structure search
    :type smiles: str
    :param sub_search_methode: what method to use for the structure search
    :type sub_search_methode: str
    :param threshold: threshold value for how alike the compounds should be to the smiles code
    :type threshold: float
    :param source_table: what table to look for compounds in. (not sure if this one makes sense...
    :type source_table: str
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param search_limiter: A dict over values to search for in the db
    :type search_limiter: dict
    :return:
        - treedata: The data for the "tree" table
        - all_data: A dict over all the data
        - rows: A dict for each row in the database
        - counter: Number of compounds
    :rtype:
        - PySimpleGUI.PySimpleGUI.TreeData
        - dicts
        - dicts
        - int
    """

    fd = FetchData(config)
    all_data = {}
    all_data_headlines = ["compound_list", "liquid_warning_list", "row_data", "mp_data", "mp_mapping", "plate_count"]

    temp_all_data = _compound_list(config, mp_amount, min_mp, samples_per_plate, ignore_active, sub_search, smiles,
                                     sub_search_methode, threshold, source_table, fd, search_limiter)

    if not temp_all_data:
        return None
    else:
        for data_index, values in enumerate(temp_all_data):
            all_data[all_data_headlines[data_index]] = values

        # if source_table == "join_main_mp":
        #     rows = fd.list_to_rows(all_data["compound_list"], source_table)
        # elif source_table == config["Tables"]["compound_main"]:
        #     rows = fd.list_to_rows(all_data["compound_list"], source_table)

        if source_table == "join_main_mp":
            search_limiter["join_tables"][config["Tables"]["compound_mp_table"]]["compound_id"]["value"] = \
                all_data["compound_list"]
            search_limiter["join_tables"][config["Tables"]["compound_mp_table"]]["compound_id"]["use"] = True
            search_limiter_tree = search_limiter["join_tables"]

        elif source_table == config["Tables"]["compound_main"]:
            search_limiter_tree = {source_table: {"value": all_data["compound_list"],
                                                  "operator": "IN",
                                                  "target_column": "compound_id",
                                                  "use": True}}
        rows = {}

        temp_dict = fd.data_search(source_table, search_limiter_tree)
        for key, value in temp_dict.items():
            rows[key] = value

        counter = 0
        treedata = sg.TreeData()

        for compound_id in rows:

            temp_list = []
            for key in rows[compound_id]:
                if key == "png":
                    temp_png = rows[compound_id][key]
                else:
                    temp_list.append(rows[compound_id][key])
            counter += 1
            if counter < 100:
                treedata.Insert("", compound_id, "", temp_list, icon=temp_png)
            else:
                treedata.Insert("", compound_id, "", temp_list, icon="")

        return treedata, all_data, rows, counter


def compound_export(folder, compound_list):
    """
    Export the list of compounds to a CSV file

    :param folder: The destination folder for the data
    :type folder: str
    :param compound_list: A list of all the compounds that needs to be extrated from the freezer
    :type compound_list: list
    :return: A CSV file that can be used for the comPOUND freezer
    """

    csvw = CSVWriter()
    csvw.compound_freezer_writer(folder, compound_list)


def compound_counter(config, table):
    """
    Gets the amount of compounds for a specific table

    :param table: The table for the compounds
    :type table: str
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :return: The number of compounds in the table
    :rtype: int
    """
    fd = FetchData(config)
    return len(fd.data_search(table, None))


def update_database(config, table, data, file_type=None):
    """
    Update the database with data. Both for adding data to the different tables, but updating tables with new values

    :param data: Output file from the plate_butler system, tube files for the comPOUND freezer, or other data.
    :type data: str or dict
    :param file_type: If the file is a CSV file, it needs a specific file_type to know how to handle the file. as the
        CSV files are different depending on where they are coming from.
    :type file_type: str
    :param table: What table needs to be updated/where the data is going
    :type table: str
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :return: Updated database with values
    """

    ad_db = AddData(config)
    ad_db.add_controller(table, data, file_type)


def purity_handler(folder, uv_one, uv_same_wavelength, wavelength, uv_threshold, rt_solvent_peak, ms_delta, ms_mode,
                   ms_threshold):
    """
    Takes raw data from LC/MS (UV and MS data) and mass from the database, per compound. and see if they can find the
    mass in the raw data, and then find the purity of the peak (based on retention time) for each compound

    :param folder: Folder with the raw-data
    :type folder: str
    :param uv_one: If it uses a single wavelength per compound or uses the whole PDA-range
    :type uv_one: bool
    :param uv_same_wavelength: if it uses the same wavelength, what wavelength that is
    :type uv_same_wavelength: bool
    :param wavelength: set to all, if there is no wavelength for each individuel compound
    :type wavelength: float
    :param uv_threshold: minimum threshold for the UV signal. anything below will be ignored
    :type uv_threshold: float
    :param rt_solvent_peak: retention time for the solvent peak
    :type rt_solvent_peak: float
    :param ms_delta: When looking for the MS data, how precise the mass should fit with the data.
    :type ms_delta: float
    :param ms_mode: if you are looking at the positive or negative. Have not been set up to look at both yet.
    :type ms_mode: str
    :param ms_threshold: minimum threshold for the MS signal. anything below will be ignored
    :type ms_threshold: float
    :return: An updated database with MS data. I think ??
    """

    lc_h = LCMSHandler()
    file_list = get_file_list(folder)
    compound_info = lc_h.lc_controller(file_list, uv_one, uv_same_wavelength, wavelength, uv_threshold, rt_solvent_peak,
                                       ms_delta, ms_mode, ms_threshold)
    print("MISSING FILE TYPE!!! ")
    #update_database(compound_info, "purity_data", )


def draw_plate(config, graph, plate_type, well_data_dict, gui_tab, archive_plate=False, sample_layout=None, mapping=None
               , state_dict=None):
    """
    Draws different plate type in on a canvas/graph.

    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param graph: The canvas/graph that is setup in sg.window
    :type graph: PySimpleGUI.PySimpleGUI.Graph
    :param plate_type: what platetype that it needs to draw. there are 3 options 96, 384 or 1536.
    :type plate_type: str
    :param well_dict: A dict over wells, this is used for drawing saved plate layouts. The dict hold information for
        what type each well is (sample, blank or paint) and the colour of the well to draw, or the value of the sample
        if the dict is from experimental data.
    :type well_dict: dict
    :param archive_plate: bool to see if it needs to draw a blank plate or a saved plate
    :type archive_plate: bool
    :param gui_tab: what tab the plate is drawn on. differet tabs differe sizes:
    :type gui_tab: str
    :param sample_layout: This is for single point, or multiple samples with same ID.
    :type sample_layout: str
    :param mapping: Information to colour wells in specific colours, depending on what state mapping is used.
        There are 3 states - state Mapping, heatmap and hit mapping
    :type mapping: dict
    :param state_dict: A dict over what state each sample is in
    :type state_dict: dict
    :return:
        - well_dict: a dict over the wells, name, state, colour and number.
        - min_x: coordinate boundaries for the plate on the canvas
        - min_y: coordinate boundaries for the plate on the canvas
        - max_x: coordinate boundaries for the plate on the canvas
        - max_y: coordinate boundaries for the plate on the canvas
    :rtype:
        - dict
        - int
        - int
        - int
        - int
    """
    if mapping and mapping["mapping"] == "Heatmap":

        heatmap = Heatmap()

        heatmap_dict = heatmap.dict_convert(well_data_dict, state_dict, mapping["states"])

        colour_dict, well_percentile, max_values, min_values = heatmap.heatmap_colours(heatmap_dict, mapping["percentile"], mapping["colours"])

    well_dict = {}
    graph.erase()
    if gui_tab == "bio":
        well_size = 20
        start_x = 5
        start_y = 165
    elif gui_tab == "bio_exp":
        well_size = 20
        start_x = 5
        start_y = 165
    else:
        well_size = 40
        start_x = 10
        start_y = 335

    fill_colour = config["plate_colouring"]["empty"]
    line_colour = "black"
    well_state = "empty"
    size = {"plate_96": well_size, "plate_384": well_size / 2, "plate_1536": well_size / 4}
    start_x = start_x + size[plate_type]
    rows = {"plate_96": 12, "plate_384": 24, "plate_1536": 48}
    columns = {"plate_96": 8, "plate_384": 16, "plate_1536": 32}
    well_id_col = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T",
                   "U", "V", "W", "X", "Y", "Z", "AA", "AB", "AC", "AD", "AE", "AF"]
    sample_layout_dict = {
        "single point": 1,
        "duplicate": 2,
        "triplicate": 3
    }
    counter = 0
    sample_counter = 0
    group = 1

    for row in range(rows[plate_type]):
        for column in range(columns[plate_type]):
            bottom_left = (start_x + row * size[plate_type],
                           start_y - column * size[plate_type])
            top_right = (bottom_left[0] - size[plate_type],
                         bottom_left[1] - size[plate_type])
            well_id = f"{well_id_col[column]}{row+1}"
            if archive_plate:
                counter += 1
                # print(well_dict)
                well_state = well_data_dict[well_id]["state"]
                fill_colour = config["plate_colouring"][well_state]

            if sample_layout and sample_layout.casefold() != "single point":
                if well_state == "sample":
                    sample_counter += 1
                    temp_colour = group % 200
                    if group % 2 == 0:
                        fill_colour = f"#FFFF{format(temp_colour, '02x')}"
                    else:
                        fill_colour = f"#FF{format(temp_colour, '02x')}FF"
                    if sample_counter % sample_layout_dict[sample_layout] == 0:
                        group += 1
            else:
                group = counter

            if mapping:
                if mapping["mapping"] == "Heatmap":
                    if state_dict[well_id]["state"] in mapping["states"]:
                        try:
                            fill_colour = heatmap.get_well_colour(colour_dict, well_percentile, well_data_dict, well_id)
                        except ZeroDivisionError:
                            fill_colour = "#FFFFFF"
                    else:
                        fill_colour = "#FFFFFF"
                elif mapping["mapping"] == "Hit Map":

                    if state_dict[well_id]["state"] in mapping["states"]:
                        if mapping["lower_bound_start"] < well_data_dict[well_id] < mapping["lower_bound_end"]:
                            fill_colour = mapping["low_colour"]
                        elif mapping["middle_bound_start"] < well_data_dict[well_id] < mapping["middle_bound_end"]:
                            fill_colour = mapping["mid_colour"]
                        elif mapping["higher_bound_start"] < well_data_dict[well_id] < mapping["higher_bound_end"]:
                            fill_colour = mapping["high_colour"]
                        else:
                            fill_colour = "#FFFFFF"
                    else:
                        fill_colour = "#FFFFFF"

            temp_well = graph.DrawRectangle(bottom_left, top_right, line_color=line_colour, fill_color=fill_colour)
            well_dict[temp_well] = {}
            well_dict[temp_well]["group"] = group
            well_dict[temp_well]["well_id"] = well_id
            well_dict[temp_well]["state"] = well_state
            well_dict[temp_well]["colour"] = fill_colour

    min_x = start_x - size[plate_type]
    min_y = start_y - (columns[plate_type] * size[plate_type])
    max_x = start_x + (rows[plate_type] * size[plate_type])-size[plate_type]
    max_y = start_y

    return well_dict, min_x, min_y, max_x, max_y


def bio_compound_info_from_worklist(config, sg, bio_sample_list):
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
        msg, file_headlines, sample_dict = CSVReader.echo_worklist_to_dict(config, config_headline, files, right_headlines, new_headline,
                                                                           sample_dict)

        # if the file uses wrong names for the headlines, this will give a popup with the "wrong headlines" and an
        # option to change them.
        if msg == "Wrong headlines":
            print(files)
            new_headline = new_headlines_popup(right_headlines, file_headlines)
            msg, file_headlines, sample_dict = CSVReader.echo_worklist_to_dict(config, config_headline, files, right_headlines, new_headline,
                                                                               sample_dict)
        elif msg == "Not a CSV file":
            sg.popup_error("Wrong file fomate!!!")              # ToDo sort out this guard so it dose not crash.
            continue

    return sample_dict


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


def bio_data(config, folder, plate_layout, archive_plates_dict, bio_plate_report_setup, analysis, bio_sample_dict,
             save_location, add_compound_ids, write_to_excel=True):

    """
    Handles the Bio data.

    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param folder: The folder where the raw data is located
    :type folder: str
    :param plate_layout: Either the name of the default plate_layout that all plates are using or a dict over each
        plate and the layout for that place
    :type plate_layout: str or dict
    :param archive_plates_dict: The layout for the plate with values for each well, what state they are in
    :type archive_plates_dict: dict
    :param bio_plate_report_setup: The setup for what is included in the report
    :type bio_plate_report_setup: dict
    :param bio_sample_dict: None or a dict of sample ide, per plate analysed
    :type bio_sample_dict: dict
    :param analysis: The analysis method
    :type analysis: str
    :param save_location: where to save all the excel files
    :type save_location: str
    :param add_compound_ids: Will add the compound ID to each well on the hit-list
    :type add_compound_ids: bool
    :return: All the data for the plates raw data, and their calculations
    :rtype: dict
    """
    # needs to reformat plate-layout to use well ID instead of numbers...
    bioa = BIOAnalyser(config, bio_plate_report_setup)
    file_list = get_file_list(folder)
    all_plates_data = {}
    used_plates = []
    for files in file_list:
        if isfile(files) and files.endswith(".txt"):
            files = txt_to_xlsx(files)
        if isfile(files) and files.endswith(".xlsx"):
            if type(plate_layout) is dict:
                temp_barcode = files.split("/")[-1].split(".")[0]
                temp_plate_layout_name = plate_layout[temp_barcode]
            else:
                temp_plate_layout_name = plate_layout
            if temp_plate_layout_name == "skip":
                continue
            else:
                temp_plate_layout = archive_plates_dict[temp_plate_layout_name]
                all_data, well_row_col, well_type, barcode, date = original_data_dict(files, temp_plate_layout)
                if not all_data:
                    return False
                used_plates.append(barcode)
                all_plates_data[barcode] = bioa.bio_data_controller(files, temp_plate_layout, all_data, well_row_col, well_type,
                                                                    analysis, write_to_excel, bio_sample_dict,
                                                                    save_location, add_compound_ids)
        else:
            print(f"{files} is not the right formate")
    return True, all_plates_data, date, used_plates


def bio_full_report(config, analyse_method, all_plate_data, used_plates, final_report_setup, output_folder, final_report_name,
                    include_hits, threshold, hit_amount, include_smiles, bio_sample_dict):
    """
    Writes the final report for the bio data
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param analyse_method: The analysed method used for the data
    :type analyse_method: str
    :param all_plate_data: All the data for all the plates, raw and calculations
    :type all_plate_data: dict
    :param used_plates: A list of plates that are used
    :type used_plates: list
    :param final_report_setup: The settings for the report
    :type final_report_name: dict
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
    :return: A excel report file with all the data
    """

    output_file = f"{output_folder}/{final_report_name}.xlsx"
    bio_final_report_controller(config, analyse_method, all_plate_data, used_plates, output_file, final_report_setup,
                                include_hits, threshold, hit_amount, include_smiles, bio_sample_dict)


def mp_production_2d_to_pb_simulate(folder_output, barcodes_2d, mp_name, trans_vol):
    """
    A simulation modul, for simulating output data

    :param folder_output: Output folder
    :type folder_output: str
    :param barcodes_2d: The folder with the 2D barcodes
    :type barcodes_2d: str
    :param mp_name: The name used for the MotherPlates
    :type mp_name: str
    :param trans_vol: The amount of volume to transferee
    :type trans_vol: float
    :return: CSV file resembling the one produced by the PB
    """
    barcodes_2d_files = get_file_list(barcodes_2d)
    csvc = CSVConverter()
    csvc.mp_in_to_out(folder_output, barcodes_2d_files, mp_name, trans_vol)


def compound_freezer_to_2d_simulate(tube_file, output_folder):
    """
    A simulation modul, for simulating output data

    :param tube_file: The CSV file with all the tube ID's for the comPOUND freezer
    :type tube_file: str
    :param output_folder: The output folder for the CSV files
    :type output_folder: str
    :return: A lot of CSV files. 1 per 96 compounds
    """
    csv_w = CSVWriter()
    csv_r = CSVReader
    tube_dict = csv_r.tube_list_to_dict(tube_file)
    csv_w.compound_freezer_to_2d_csv_simulate(tube_dict, output_folder)


def bio_experiment_to_database(assay_name, plate_data, plate_layout, date, responsible, config, bio_files):


    table = "bio_experiment"

    raw_data = f"{assay_name}_{date}"
    exp_id = get_number_of_rows(config, table) + 1

    data_dict = {
        "exp_id": exp_id,
        "assay_name": assay_name,
        "raw_data": raw_data,
        "plate_layout": plate_layout,
        "responsible": responsible,
        "date": date
    }
    temp_dict = {raw_data: plate_data}
    dict_writer(bio_files, temp_dict)
    update_database(data_dict, table, None, config)


def grab_table_data(config, table_name, search_limiter=None):
    fd = FetchData(config)

    rows = fd.data_search(table_name, search_limiter)

    all_table_data = []
    headlines = []

    for row in rows:
        temp_data = []
        for data in rows[row]:
            headlines.append(data)
            temp_data.append(rows[row][data])

        all_table_data.append(temp_data)

    return all_table_data, headlines


def update_bio_info_values(values, window, plate_bio_info):
    temp_plate_name = values["-BIO_INFO_PLATES-"]
    temp_analyse_method = values["-BIO_INFO_ANALYSE_METHOD-"]
    temp_state = values["-BIO_INFO_STATES-"]

    window["-INFO_BIO_AVG-"].update(
        value=plate_bio_info[temp_plate_name]["calculations"][temp_analyse_method][temp_state]["avg"])
    window["-INFO_BIO_STDEV-"].update(
        value=plate_bio_info[temp_plate_name]["calculations"][temp_analyse_method][temp_state]["stdev"])
    window["-INFO_BIO_Z_PRIME-"].update(value=plate_bio_info[temp_plate_name]["calculations"]["other"]["z_prime"])


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


def dp_creator(config, plate_layout, sample_amount, mp_data, transferee_volume, dp_name, output_folder):
    csv_w = CSVWriter()

    dp_layout = plate_layout_to_well_ditc(plate_layout)
    # generate a dict for Daughter_Plates
    dp_dict = daughter_plate_generator(mp_data, sample_amount, dp_name, dp_layout, transferee_volume)

    # generate CSV-file for PlateButler
    csv_w.dp_writer(config, dp_dict, output_folder)

    # Generate list over mp needed
    csv_w.plate_list_writer(mp_data, output_folder)


def update_plate_table(compound_id, config):
    plate_tables = {"compound_mp": "MP", "compound_dp": "DP"}

    search_limiter_plates = {"compound_id": {"value": compound_id,
                                             "operator": "=",
                                             "target_column": "compound_id",
                                             "use": True}}
    plate_data = []
    for tables in plate_tables:
        temp_data = []
        all_data_plate, _ = grab_table_data(config, tables, search_limiter_plates)
        if all_data_plate:
            all_data_plate = all_data_plate[0]
            if tables == "compound_mp":
                temp_data.append(all_data_plate[1])
                temp_data.append(plate_tables[tables])
                temp_data.append(all_data_plate[2])
                temp_data.append(all_data_plate[3])
                temp_data.append(all_data_plate[4])
            elif tables == "compound_dp":
                temp_data.append(all_data_plate[3])
                temp_data.append(plate_tables[tables])
                temp_data.append(all_data_plate[4])
                temp_data.append(all_data_plate[5])
                temp_data.append(all_data_plate[6])

        plate_data.append(temp_data)

    return plate_data


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


def plate_layout_to_excel(well_dict, name, folder):
    # for index, plate in enumerate(well_dict):
    #     if index == 0:
    well_col_row, well_type = well_row_col_type(well_dict)
    plate_layout = plate_layout_re_formate(well_dict)

    export_plate_layout(plate_layout, well_col_row, name, folder)


def plate_dilution(config, function, file, dw_amount, add_source_wells, source_well_amount, save_plates, well_layout,
                   plate_layout, pb_source_file):
    output_folder = config["folders"]["main_output_folder"]
    destination_plate_naming_scheme = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]

    if function == "Calculate":
        state = plate_dilution_write_vol_well_amount(config, file, dw_amount, add_source_wells, source_well_amount,
                                                     well_layout)
        return state

    elif function == "Generate":
        plate_layout = plate_layout_re_formate(plate_layout["well_layout"])
        sample_info_dict, replicate_samples_max, replicate_plate_sets, dilution_factor, concentration_counter, \
            control_vol, control_conc = plate_dilution_excel(file, save_plates, dw_amount, well_layout)

        pd = PlateDilution(config, save_plates, pb_source_file, control_vol, control_conc)

        sample_dict, dilution_dict = pd.pd_controller(sample_info_dict, replicate_samples_max, replicate_plate_sets,
                                                      dilution_factor, concentration_counter, plate_layout,
                                                      destination_plate_naming_scheme, dw_amount)


        csv_w = CSVWriter()
        if pb_source_file:
            csv_w.source_plate_dilution(dilution_dict, output_folder)

        csv_w.plate_dilution(sample_dict, output_folder)

        return "CSV file created"


def _sample_peak_dict_creator(peak_table_data):
    sample_peak_dict = {}
    for samples in peak_table_data:
        sample_peak_dict[samples] = {}
        for peaks in peak_table_data[samples]:
            sample_peak_dict[samples][peaks[1]] = {
                "start": peaks[3],
                "end": peaks[4],
                "area": peaks[2],
                "%": peaks[5]
            }
    return sample_peak_dict


def get_peak_information(purity_data, slope_threshold, uv_threshold, rt_solvent_peak, sample_data, wavelength_data,
                         sample=None):
    ms_int = Integration()

    if wavelength_data.endswith("xlsx"):
        wavelength_data = _get_sample_data(wavelength_data)

    peak_information = ms_int.calculate_uv_integrals(purity_data, slope_threshold, uv_threshold, rt_solvent_peak, sample_data, wavelength_data, sample)

    peak_table_data = {}
    for samples in peak_information:
        peak_table_data[samples] = []
        peak_info_dict = peak_information[samples].to_dict("index")
        for row in peak_info_dict:
            peak_info = [
                samples, peak_info_dict[row]["Peak list"], peak_info_dict[row]["Integrals"],
                peak_info_dict[row]["Peak start time"], peak_info_dict[row]["Peak end time"],
                peak_info_dict[row]["purity"]
            ]
            peak_table_data[samples].append(peak_info)

    sample_peak_dict = _sample_peak_dict_creator(peak_table_data)

    return peak_information, peak_table_data, sample_peak_dict


def import_ms_data(folder):

    file_list = get_file_list(folder)
    data = dm_controller(file_list)

    if isinstance(data, str):
        return data

    else:
        purity_data = data[0]
        missing_samples = data[1]

        samples = []
        table_data = []
        for sample in purity_data:
            samples.append(sample)
            temp_data = []
            for data in purity_data[sample]:
                temp_data.append(purity_data[sample][data])
            table_data.append(temp_data)

                # Makes a dict of batches, for adding to the batch database
        all_data = [table_data, samples, purity_data, missing_samples]
    return all_data


def purity_data_to_db(config, purity_data):

    batch_dict = {}
    today = date.today()
    today = today.strftime("%m_%Y")
    folder = f"{config['folders']['main_output_folder']}/{config['folders']['purity_data']}"
    try:
        mkdir(folder)
    except FileExistsError:
        pass

    path = f"{folder}/{today}.txt"
    batch_list = []
    table = "lc_raw"

    row = get_number_of_rows(config, table)

    for samples in purity_data:
        row += 1
        batch_dict[purity_data[samples]["batch"]] = {"batch": purity_data[samples]["batch"],
                                                     "date": purity_data[samples]["date"]}

        for batch in batch_dict:
            # makes sure that only new batches are added to the database
            if batch not in batch_list:
                update_database(batch_dict[batch], "lc_experiment", None, config)
                batch_list.append(batch)

        temp_file_date = {}
        temp_file_date[f"{purity_data[samples]['sample']}_{purity_data[samples]['batch']}"] = {
            "uv": purity_data[samples]["uv"],
            "ms_neg": purity_data[samples]["ms_neg"],
            "ms_pos": purity_data[samples]["ms_pos"]}

        # Writes UV and MS data to a separated file. Name based on sample_batch.
        df_writer(path, temp_file_date)
        temp_data_ditc = {"row_id": row,
                          "sample": purity_data[samples]["sample"],
                          "batch": purity_data[samples]["batch"],
                          "method": purity_data[samples]["method"],
                          "file_name": path,
                          "date": purity_data[samples]["date"]}
        update_database(temp_data_ditc, table, None, config)
    return batch_dict


def get_number_of_rows(config, table):
    dbf = DataBaseFunctions(config)
    return dbf.number_of_rows(table)


def purity_data_compounds_to_db(config, table_data):

    table = "purity"
    dbf = DataBaseFunctions(config)
    row = get_number_of_rows(config, table_data)

    for rows in table_data:

        batch = rows[1]
        temp_table = "lc_experiment"
        search_limiter = {temp_table: {"value": [batch], "operator": "IN", "target_column": "batch", "use": True}}
        row_data = dbf.return_table_data(temp_table, search_limiter)
        date = row_data[batch]["date"]
        row += 1
        temp_data_ditc = {"purity_id": row,
                          "compound_id": rows[0],
                          "batch": batch,
                          "result_max": rows[3],
                          "result_max_ion": rows[4],
                          "result_total": rows[7],
                          "date": date}
        update_database(temp_data_ditc, table, None, config)


def purity_plotting(method, data, canvas, samples, fig_size, ms_mode, rt_start, rt_end, wavelength, bin_numbers,
                    mz_value, canvas_lines):
    if method == "uv_chromatogram":
        try:
            figure_canvas_agg = uv_chromatogram(data, canvas, samples, fig_size, canvas_lines)
        except TypeError:
            return "Missing Values"
    elif method == "ms_chromatogram":
        try:
            figure_canvas_agg = ms_chromatogram(data, canvas, samples, fig_size, ms_mode)
        except TypeError:
            return "Missing Values"
    elif method == "ms_spectrum":
        try:
            figure_canvas_agg = ms_spectrum(data, canvas, samples, fig_size, ms_mode, rt_start)
        except TypeError:
            return "Missing Values"
    elif method == "ms_spectrum_range":
        try:
            figure_canvas_agg = ms_spectrum_range(data, canvas, samples, fig_size, ms_mode, canvas_lines)
        except TypeError:
            return "Missing Values"
    elif method == "heatmap_uv_sample":
        try:
            figure_canvas_agg = heatmap_uv_sample(data, canvas, samples, fig_size)
        except TypeError:
            return "Missing Values"
    elif method == "heatmap_uv_rt":
        try:
            figure_canvas_agg = heatmap_uv_rt(data, canvas, samples, fig_size, rt_start)
        except TypeError:
            return "Missing Values"
    elif method == "heatmap_uv_wavelength":
        try:
            figure_canvas_agg = heatmap_uv_wavelength(data, canvas, samples, fig_size, wavelength)
        except TypeError:
            return "Missing Values"
    elif method == "heatmap_ms_sample_binned":
        try:
            figure_canvas_agg = heatmap_ms_sample_binned(data, canvas, samples, fig_size, ms_mode, bin_numbers)
        except TypeError:
            return "Missing Values"
    elif method == "heatmap_ms_rt_binned":
        try:
            figure_canvas_agg = heatmap_ms_rt_binned(data, canvas, samples, fig_size, ms_mode, bin_numbers, rt_start)
        except TypeError:
            return "Missing Values"
    elif method == "heatmap_ms_mz":
        try:
            figure_canvas_agg = heatmap_ms_mz(data, canvas, samples, fig_size, ms_mode, mz_value)
        except TypeError:
            return "Missing Values"
    else:
        return "No viz selected"
    return figure_canvas_agg


def _purity_mass(purity_data, uv_peak_information, mass_hit):
    """uses area data combined with mass to find purity"""

    for sample in purity_data:

        if mass_hit[sample] == "No Hits":
            pass
        else:
            for peak in mass_hit[sample]:
                for index, row in uv_peak_information[sample].iterrows():
                    if row[0] == peak:
                        mass_hit[sample][peak]["purity"] = row[4]

        purity_data[sample]["peak_hits"] = mass_hit[sample]


def _purity_overview_table_data_creation(purity_data, sample_data):
    purity_overview_table_data = []
    no_hits = False
    purity_peak_list_table_data = {}
    for sample in purity_data:
        purity_peak_list_table_data[sample] = []
        temp_sample_info = []
        temp_purity = []
        temp_ion = []
        for peaks in purity_data[sample]["peak_hits"]:
            purity = purity_data[sample]["peak_hits"][peaks]["purity"]
            ion_info = purity_data[sample]["peak_hits"][peaks]
            temp_purity.append(purity)
            temp_ion.append(ion_info)
            mass = purity_data[sample]["peak_hits"][peaks][list(ion_info)[0]]
            peak_list = [peaks, list(ion_info)[0], mass, purity]
            purity_peak_list_table_data[sample].append(peak_list)
        temp_sample_info.append(sample)
        temp_sample_info.append(purity_data[sample]["batch"])
        temp_sample_info.append(sample_data[sample]["mass"])
        try:
            temp_sample_info.append(max(temp_purity))
        except ValueError:
            temp_sample_info.append("No Hits")
            no_hits = True
        if not no_hits:
            for index_hits, mass_hits in enumerate(temp_purity):
                if mass_hits == max(temp_purity):
                    ion = list(temp_ion[index_hits])[0]
                    temp_sample_info.append(ion)
                    temp_sample_info.append(temp_ion[index_hits][ion])
                    temp_sample_info.append(peaks)
            temp_sample_info.append(sum(temp_purity))
        else:
            temp_sample_info.append("No Hits")
            temp_sample_info.append("No Hits")
            temp_sample_info.append("No Hits")
            temp_sample_info.append("No Hits")
            no_hits = False
        purity_overview_table_data.append(temp_sample_info)

    return purity_overview_table_data, purity_peak_list_table_data


def _get_sample_data(sample_data_file):
    if sample_data_file.endswith(".csv"):
        sample_data = CSVReader.compound_plates(sample_data_file)
    elif sample_data_file.endswith(".xslx"):
        sample_data = pd.read_excel(sample_data_file, index_col=0)
        sample_data = sample_data.to_dict("index")
    else:
        return "Raw data do not match the right formate"
    return sample_data


def add_start_end_time(purity_peak_list_table_data, sample_peak_dict):
    for samples in purity_peak_list_table_data:
        for index, peak_data in enumerate(purity_peak_list_table_data[samples]):
            purity_peak_list_table_data[samples][index].append(sample_peak_dict[samples][peak_data[0]]["start"])
            purity_peak_list_table_data[samples][index].append(sample_peak_dict[samples][peak_data[0]]["end"])


def purity_ops(sample_data, purity_data, peak_information, ms_mode, delta_mass, mz_threshold, peak_amounts,
               mass=None):
    """

    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param sample_data_file:
    :param purity_data:
    :param peak_information:
    :param ms_mode:
    :param delta_mass:
    :param mz_threshold:
    :param peak_amounts:
    :param mass:
    :return:
    """

    mass_hit = mass_search(purity_data, peak_information, ms_mode, sample_data, delta_mass, mz_threshold, peak_amounts,
                           mass)

    _purity_mass(purity_data, peak_information, mass_hit)

    purity_overview_table_data, purity_peak_list_table_data = _purity_overview_table_data_creation(purity_data, sample_data)

    return purity_overview_table_data, purity_peak_list_table_data


def grab_sample_data(sample_data_file, purity_data, config):
    if sample_data_file == "compound_data":
        dbf = DataBaseFunctions(config)
        co = ChemOperators
        table = "compound_main"
        sample_data = {}
        for samples in purity_data:
            Search_limiter = {table:
                                  {"value": [samples],
                                   "operator": "IN",
                                   "target_column": "compound_id",
                                   "use": True}}
            temp_table_data = dbf.return_table_data(table, Search_limiter)
            smiles = temp_table_data[int(samples)]["smiles"]
            mass = co.mw_from_smiles(smiles)
            sample_data[samples] = {"mass": mass}
        db_data = True
    elif sample_data_file:
        sample_data = _get_sample_data(sample_data_file)
        db_data = False

    return sample_data, db_data


def grab_compound_table_data(config, table, sample):
    dbf = DataBaseFunctions(config)
    Search_limiter = {table:
                          {"value": [sample],
                           "operator": "IN",
                           "target_column": "compound_id",
                           "use": True}
                      }

    return dbf.return_table_data(table, Search_limiter)


def name_changer(new_names, raw_data_dict, excel_data_dict, peak_information, sample_peak_dict):

    for names in new_names:
        if names != new_names[names]["raw"]:
            raw_data_dict[names] = raw_data_dict.pop(new_names[names]["raw"])
            peak_information[names] = peak_information.pop(new_names[names]["raw"])
            sample_peak_dict[names] = sample_peak_dict.pop(new_names[names]["raw"])
        elif names != new_names[names]["excel"]:
            excel_data_dict[names] = excel_data_dict.pop(new_names[names]["excel"])


def compound_info_table_data(config, sample):
    all_plate_info_table_data = []
    purity_info_table_data = []
    bio_info_table_data = []
    mp_plate_info_table_data = []
    dp_plate_info_table_data = []
    true_false = {0: False, 1: True}

    # Purity data
    # Batch, Max-purity, Ion, Total-purity
    temp_plate_data = grab_compound_table_data(config, "purity", sample)
    for rows in temp_plate_data:
        purity_info_table_data.append([temp_plate_data[rows]["batch"], temp_plate_data[rows]["result_max"],
                                      temp_plate_data[rows]["result_max_ion"],  temp_plate_data[rows]["result_total"],
                                      temp_plate_data[rows]["date"]])

    # all plate data:
    # mp/dp, plate_id, date, active

    # MP data:
    # MP:
    # plate ID, plate_type, well, vol, freeze_thaw cycles, date, location, active

    temp_plate_data = grab_compound_table_data(config, "compound_mp", sample)
    for rows in temp_plate_data:
        all_plate_info_table_data.append(["MP", temp_plate_data[rows]["mp_barcode"], temp_plate_data[rows]["date"],
                                          true_false[temp_plate_data[rows]["active"]]])

        mp_plate_info_table_data.append([temp_plate_data[rows]["mp_barcode"],
                                         temp_plate_data[rows]["plate_type"], temp_plate_data[rows]["mp_well"],
                                         temp_plate_data[rows]["volume"], temp_plate_data[rows]["freeze_thaw"],
                                         temp_plate_data[rows]["date"], temp_plate_data[rows]["location"],
                                         true_false[temp_plate_data[rows]["active"]]])

    # DP data
    # DP:
    # plate ID, plate_type, well, vol, date, active
    temp_plate_data = grab_compound_table_data(config, "compound_dp", sample)
    for rows in temp_plate_data:
        all_plate_info_table_data.append(["DP", temp_plate_data[rows]["dp_barcode"], temp_plate_data[rows]["date"],
                                          true_false[temp_plate_data[rows]["active"]]])

        dp_plate_info_table_data.append([temp_plate_data[rows]["dp_barcode"],
                                         temp_plate_data[rows]["plate_type"], temp_plate_data[rows]["dp_well"],
                                         temp_plate_data[rows]["volume"], temp_plate_data[rows]["date"],
                                         true_false[temp_plate_data[rows]["active"]]])

    # Bio data:
    # Still missing details about what data to save, and what to show...
    temp_plate_data = grab_compound_table_data(config, "biological", sample)

    for rows in temp_plate_data:
        bio_ex_id = temp_plate_data[rows]["exp_id"]
        temp_bio_exp_data = grab_compound_table_data(config, "bio_experiment", bio_ex_id)
        for row in temp_bio_exp_data:
            bio_info_table_data.append([temp_bio_exp_data[row]["assay_name"], temp_bio_exp_data[row]["responsible"],
                                        temp_bio_exp_data[row]["date"]])

    return all_plate_info_table_data, mp_plate_info_table_data, dp_plate_info_table_data, purity_info_table_data, \
        bio_info_table_data


def database_to_table(config, table, headings):
    dbf = DataBaseFunctions(config)
    row_data = dbf.return_table_data(table, None)

    table_data = []
    for row in row_data:
        temp_table_row = []
        for data in row_data[row]:
            if data in headings:
                temp_table_row.append(row_data[row][data])
        table_data.append(temp_table_row)

    return table_data


def _get_motherplates_with_wells_from_worklist_dict(used_plate_well_dict):
    """
    Makes a new dict, to get MotherPlates
    :param used_plate_well_dict: A dicts of MotherPlates and the wells that have been used
    :type used_plate_well_dict: dict
    :return: motherplate_dict, a dict with MotherPlates as Keys and wells used as a list as value
    :rtype: dict
    """

    motherplate_dict = {}
    for destinations_plates in used_plate_well_dict:

        for wells in used_plate_well_dict[destinations_plates]:
            temp_mp = used_plate_well_dict[destinations_plates][wells]["source_plate"]
            try:
                motherplate_dict[temp_mp]
            except KeyError:
                motherplate_dict[temp_mp] = []

            source_well = used_plate_well_dict[destinations_plates][wells]["source_well"]
            if source_well not in motherplate_dict[temp_mp]:
                motherplate_dict[temp_mp].append(source_well)

    return motherplate_dict


def _get_motherplate_layout(config, mps):
    """
    Generate a dict with the key of each motherplate and the value for all the wells that have compounds.
    :param mps: A list of MotherPlates
    :type: mps: list
    :return: a dict with the key of each motherplate and the value for all the wells that have compounds.
    :rtype: dict
    """

    dbf = DataBaseFunctions(config)
    table = "compound_mp"
    clm_header = "mp_barcode"
    plate_data_dict = {}
    for mp in mps:
        plate_data_dict[mp] = []
        temp_rows = dbf.find_data_single_lookup(table, mp, clm_header)
        for data in temp_rows:
            temp_well = temp_rows[data]["mp_well"]
            plate_data_dict[mp].append(temp_well)

    return plate_data_dict


def _get_free_wells(mps_layout, motherplate_dict):
    """
    Compare MPS with previuse used plates, and finds the wells that have not been used and add them to a dict
    :param mps_layout: A dict with the layout for each motherplate
    :type: mps_layout: dict
    :param motherplate_dict: A dicts of MotherPlates and the wells that have been used
    :type motherplate_dict: dict
    :return: a dict with MotherPlates as keys, and list of free wells as values
    :rtype: dict
    """
    # This depends on the motherplate layout. We use a full 384 layout
    free_wells = {}
    for motherplates in mps_layout:
        free_wells[motherplates] = []
        for wells in mps_layout[motherplates]:

            if motherplate_dict and motherplates in motherplate_dict:
                if wells not in motherplate_dict[motherplates]:
                    free_wells[motherplates].append(wells)
            else:
                free_wells[motherplates].append(wells)
    return free_wells


# def _from_plate_layout_to_destination_dict(plate_layout, destination_name, plate_amount, initial_plate):
#     """
#     Makes a dict out of plate_layout, amount of plates and naming scheme
#     :param plate_amount: Number of plates to produce
#     :type plate_amount: int
#     :param plate_layout: The layout for the plate with values for each well, what state they are in
#     :type plate_layout: dict
#     :param destination_name: the name of the assay, and the name used for the destination plate in the workinglist
#     :type destination_name: str
#     :param initial_plate: the starting number of the plates
#     :type initial_plate: int
#     :return: A dict, with destination plates as the key, and desination wells, as a list:
#     :rtype: dict
#     """
#     destination_dict = {}
#
#     plate_range = plate_amount - initial_plate
#
#     for plate in range(plate_range):
#         destination_plate = f"{destination_name}_{plate + initial_plate + 1}"
#         destination_dict[destination_plate] = []
#         for wells in plate_layout:
#             if plate_layout[wells]["state"] == "sample":
#                 destination_well = plate_layout[wells]["well_id"]
#                 destination_dict[destination_plate].append(destination_well)
#
#     return destination_dict


def generate_worklist(config, plate_amount, mps, plate_layout, used_plate_well_dict, assay_name, initial_plate, volume,
                      worklist_analyse_method, sample_direction, control_layout, control_samples, bonus_compound):
    """
    Generates a worklist based on a plate-layout and data from MotherPlates.

    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param plate_amount: Number of plates to produce
    :type plate_amount: int
    :param mps: A list of MotherPlates
    :type: mps: list
    :param plate_layout: The layout for the plate with values for each well, what state they are in
    :type plate_layout: dict
    :param used_plate_well_dict: A dicts of MotherPlates and the wells that have been used
    :type used_plate_well_dict: dict
    :param assay_name: the name of the assay, and the name used for the destination plate in the workinglist
    :type assay_name: str
    :param initial_plate: the starting number of the plates
    :type initial_plate: int
    :param volume: How much volume to transfere to each well. The same amount of liquid will be transfered to each.
    :type volume: float
    :param worklist_analyse_method: The method use for the sample layout.
    :type worklist_analyse_method: str
    :param sample_direction: The direction for the sample layout. Only relevant if the analyse_method differ from
        "Single Point"
    :type sample_direction: str
    :param control_layout: The file path to the plate layout for the control or other compounds
    :type control_layout: pathlib.WindowsPath or None
    :param control_samples: A dict for samples for positive and negative control
    :type: control_samples: None or dict
    :param bonus_compound: If there are a compound that needs to be added to multiple well_states
    :type bonus_compound: dict
    :return:
    """

    # Gets a dict of motherplates that have been used, and what wells in each motherplate,
    # based on the worklist provided.
    if used_plate_well_dict:
        motherplate_dict = _get_motherplates_with_wells_from_worklist_dict(used_plate_well_dict)
    else:
        motherplate_dict = None

    # Get the layout of the motherplates. incase the motherplate is not a full 384 plate.
    mps_layout = _get_motherplate_layout(config, mps)

    # gets all the free wells, based on the layout of each motherplate, and what wells have been used before,
    # based on the worklist provided
    free_well_dict = _get_free_wells(mps_layout, motherplate_dict)

    if not control_layout:
        control_bonus_source = None
    elif control_layout.suffix == ".xlsx":
        control_bonus_source = well_compound_list(control_layout)
        if type(control_bonus_source) == str:
            return control_bonus_source
    elif control_layout.suffix == ".txt":
        pass            # ToDo add these ?
    elif control_layout.suffix == ".csv":
        pass            # ToDo add these ?
    else:
        return "error: wrong file format for Control Layout"
    # destination_dict = _from_plate_layout_to_destination_dict(plate_layout, assay_name, plate_amount, initial_plate)
    csv_w = CSVWriter()
    file, msg = csv_w.worklist_writer(config, plate_layout, mps, free_well_dict, assay_name, plate_amount,
                                      initial_plate, volume, sample_direction, worklist_analyse_method,
                                      control_bonus_source, control_samples, bonus_compound)

    if type(msg) == str:
        print(msg)
    return file, msg


def _get_plate_names(dbf, table, column_headline):
    """
    Gets a list of the plate names from the database
    :param dbf: The DataBaseFunction
    :type dbf: class
    :param table: The table where the data is storage
    :type table: str
    :param column_headline: The headline of the column that is used for finding data
    :type column_headline: str
    :return:
    """
    temp_rows = dbf.find_column_data(table, column_headline)
    plate_names = [names for names in temp_rows]
    return plate_names


def _get_plate_archive(dbf, table, column_headline, plate_names):
    """
    Gets a dict of the plate layouts from the database, based on a list of plate names
    :param dbf: The DataBaseFunction
    :type dbf: class
    :param plate_names: A list of plate names
    :type plate_names: list
    :param table: The table where the data is storage
    :type table: str
    :param column_headline: The headline of the column that is used for finding data
    :type column_headline: str
    :return: a dict of the plate layouts
    :rtype: dict
    """
    archive_plates_dict = {}

    # loops over all the plates
    for plates in plate_names:
        temp_row_data = dbf.find_data_single_lookup(table, plates, column_headline)

        # grabs data from the return rown
        plate_name = temp_row_data[0][1]
        well_layout = temp_row_data[0][2]
        plate_type = temp_row_data[0][3]

        # Generates the dict
        archive_plates_dict[plate_name] = {"well_layout": eval(well_layout),
                                           "plate_type": plate_type}
    return archive_plates_dict


def get_plate_layout(config):
    """
    Gets a list of plate names and their layout from the database
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :return: plate_names, archive_plates_dict
    :rtype: list, dict
    """

    # Connects to the database and setting up standard values
    dbf = DataBaseFunctions(config)
    table = "plate_layout"
    column_headline = "plate_name"

    # gets a list of the plate names
    plate_names = _get_plate_names(dbf, table, column_headline)

    # Gets a dict over all the plate_layouts
    archive_plates_dict = _get_plate_archive(dbf, table, column_headline, plate_names)

    return plate_names, archive_plates_dict


def delete_records_from_database(config, table, headline, data_value):
    """
    Deletes a record from the database
    :param table: What table are the plates in
    :type table: str
    :param data_value: The value of the thing you are looking for
    :type data_value: str
    :param headline: Headline for the coloumn where  the data is, in the table
    :type headline: str
    :return:
    """
    dbf = DataBaseFunctions(config)
    dbf.delete_records(table, headline, data_value)


def rename_record_in_the_database(config, table, headline, old_value, new_value):
    """
    Deletes a record from the database
    :param table: What table are the plates in
    :type table: str
    :param data_value: The value of the thing you are looking for
    :type data_value: str
    :param headline: Headline for the coloumn where  the data is, in the table
    :type headline: str
    :return:
    """
    dbf = DataBaseFunctions(config)
    dbf.rename_record_value(table, headline, old_value, new_value)

if __name__ == "__main__":

    sg = "hej"
    csv_file = "C:/Users/phch/Desktop/more_data_files/50sets_picklist_updated_worklists_new.txt"
    bio_sample_list = [csv_file]
    sample_dict = bio_compound_info_from_worklist(sg, bio_sample_list)
    print(sample_dict)
