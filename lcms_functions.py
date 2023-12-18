from datetime import date
from os import mkdir
from PySimpleGUI import TreeData
import pandas as pd

from bio_data_functions import well_row_col_type
from chem_operators import mw_from_smiles
from file_type_handler_csv import CSVWriter, CSVReader
from database_controller import FetchData
from database_functions import get_number_of_rows, update_database, grab_compound_table_data
from database_handler import DataBaseFunctions
from file_type_handler_excel import export_plate_layout
from file_handler import get_file_list
from lcms_data_handler import LCMSHandler
from lcms_data_miner import dm_controller
from lcms_ms_search import mass_search
from lcms_uv_integration import Integration
from import_pickle_handler import df_writer
from compound_plate_formatting import plate_layout_to_well_ditc, daughter_plate_generator, plate_layout_re_formate
from lcms_visualization import uv_chromatogram, ms_chromatogram, ms_spectrum, ms_spectrum_range, heatmap_uv_sample, \
    heatmap_uv_rt, heatmap_uv_wavelength, heatmap_ms_sample_binned, heatmap_ms_rt_binned, heatmap_ms_mz


def calculate_sample_amount(plate_amount, samples_per_plate=384):
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


def _compound_list(config, mp_amount, min_mp, samples_per_plate, ignore_active, source_table, fd, search_limiter):
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
        sample_amount = calculate_sample_amount(mp_amount, samples_per_plate)
    else:
        sample_amount = None
    plated_compounds = []
    if not ignore_active:
        # plated_compounds = [compounds for compounds in fd.data_search(config["Tables"]["compound_mp_table"], None)]
        plated_compounds = [[test for row, test in data.items() if row == "compound_id"][0] for _, data in
                            fd.data_search(config["Tables"]["compound_mp_table"], None).items()]

    # Gets a list of compounds, based on search criteria
    items = fd.list_limiter(sample_amount, min_mp, samples_per_plate, source_table, ignore_active, plated_compounds,
                            search_limiter)

    return items


def table_update_tree(mp_amount, min_mp, samples_per_plate, ignore_active, source_table, search_limiter, config):
    """
    Updates the compound table with compounds depending on search criteria

    :param mp_amount: amount of mother plates to find compounds from
    :type mp_amount: int
    :param min_mp: If all the samples needs  to be from as few plates as possible
    :type min_mp: bool
    :param samples_per_plate: amount of sample per plate
    :type samples_per_plate: int
    :param ignore_active: If it should take into account witch compounds are already in MotherPlates.
        True = All compounds, False = only compounds not found in MotherPlates
    :type ignore_active: bool
    :param source_table: what table to look for compounds in. (not sure if this one makes sense...)
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

    temp_all_data = _compound_list(config, mp_amount, min_mp, samples_per_plate, ignore_active, source_table, fd,
                                   search_limiter)
    if not temp_all_data:
        return None, None, None, None
    else:
        for data_index, values in enumerate(temp_all_data):
            all_data[all_data_headlines[data_index]] = values

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
        else:
            print("HEJ!!!!!! something fucked up :D ")

        rows = {}
        temp_dict = fd.data_search(source_table, search_limiter_tree)
        for key, value in temp_dict.items():
            # print(f"{key} - {value}")
            rows[key] = value
        counter = 0
        treedata = TreeData()

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


def lcms_raw_data_handler(folder, uv_one, uv_same_wavelength, wavelength, uv_threshold, rt_solvent_peak, ms_delta, ms_mode,
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


def dp_creator(config, plate_layout, sample_amount, mp_data, transferee_volume, dp_name, output_folder):
    csv_w = CSVWriter()

    dp_layout = plate_layout_to_well_ditc(plate_layout)
    # generate a dict for Daughter_Plates
    dp_dict = daughter_plate_generator(mp_data, sample_amount, dp_name, dp_layout, transferee_volume)

    # generate CSV-file for PlateButler
    csv_w.dp_writer(config, dp_dict, output_folder)

    # Generate list over mp needed
    csv_w.plate_list_writer(mp_data, output_folder)


def plate_layout_to_excel(config, well_dict, name, folder):
    # for index, plate in enumerate(well_dict):
    #     if index == 0:
    well_col_row, well_type = well_row_col_type(well_dict)
    plate_layout = plate_layout_re_formate(config, well_dict)

    export_plate_layout(plate_layout, well_col_row, name, folder)


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

    peak_information = ms_int.calculate_uv_integrals(purity_data, slope_threshold, uv_threshold, rt_solvent_peak,
                                                     sample_data, wavelength_data, sample)

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

    if type(folder) == str:
        file_list = get_file_list(folder)
    else:
        file_list = folder.glob("*/**")

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


def lcms_data_to_db(config, purity_data):

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


def lcms_plotting(method, data, canvas, samples, fig_size, ms_mode, rt_start, rt_end, wavelength, bin_numbers,
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


def _uv_mass_to_purity(purity_data, uv_peak_information, mass_hit):
    """uses area data combined with mass to find purity"""

    for sample in purity_data:
        if "blank" in sample.casefold():
            continue

        try:
            mass_hit[sample]
        except KeyError:
            continue

        if mass_hit[sample] == "No Hits":
            pass
        else:
            for peak in mass_hit[sample]:
                for index, row in uv_peak_information[sample].iterrows():
                    if row[0] == peak:
                        mass_hit[sample][peak]["purity"] = row[4]

        purity_data[sample]["peak_hits"] = mass_hit[sample]


def _lcms_overview_table_data_creation(purity_data, sample_data):
    """

    :param purity_data:
    :param sample_data:
    :return:
    """

    purity_overview_table_data = []
    no_hits = False
    purity_peak_list_table_data = {}
    for sample in purity_data:
        if "blank" in sample.casefold():
            continue

        try:
            purity_data[sample]["peak_hits"]
        except KeyError:
            continue

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
                    sample_data[sample]["max_hit"] = {"ion": list(temp_ion[index_hits])[0],
                                                      "ion_mass": temp_ion[index_hits][ion],
                                                      "peak_id": peaks,
                                                      "purity": temp_purity}
            temp_sample_info.append(sum(temp_purity))
        else:
            sample_data[sample]["max_hit"] = False
            temp_sample_info.append("No Hits")
            temp_sample_info.append("No Hits")
            temp_sample_info.append("No Hits")
            temp_sample_info.append("No Hits")
            no_hits = False

        purity_overview_table_data.append(temp_sample_info)

    return purity_overview_table_data, purity_peak_list_table_data


def _get_sample_data(sample_data_file):
    """
    Gets data from either CSV file or xslx file
    :param sample_data_file:
    :return:
    """

    if sample_data_file.endswith(".csv"):
        sample_data = CSVReader.compound_plates(sample_data_file)
    elif sample_data_file.endswith(".xlsx"):
        sample_data = pd.read_excel(sample_data_file, index_col=0)
        sample_data = sample_data.to_dict("index")
    else:
        return "Raw data do not match the right formate"
    return sample_data


def add_start_end_time(purity_peak_list_table_data, sample_peak_dict):
    """
    Adds the start and end time to the peak data
    :param purity_peak_list_table_data:
    :param sample_peak_dict:
    :return:
    """

    for samples in purity_peak_list_table_data:
        for index, peak_data in enumerate(purity_peak_list_table_data[samples]):
            purity_peak_list_table_data[samples][index].append(sample_peak_dict[samples][peak_data[0]]["start"])
            purity_peak_list_table_data[samples][index].append(sample_peak_dict[samples][peak_data[0]]["end"])


def lcms_ops(sample_data, purity_data, peak_information, ms_mode, delta_mass, mz_threshold, peak_amounts,
             mass=None):
    """
    #todo write this ?
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

    _uv_mass_to_purity(purity_data, peak_information, mass_hit)

    purity_overview_table_data, purity_peak_list_table_data = _lcms_overview_table_data_creation(purity_data, sample_data)

    return purity_overview_table_data, purity_peak_list_table_data


def grab_sample_data(sample_data_file, purity_data, config):
    if sample_data_file == "compound_data":
        dbf = DataBaseFunctions(config)
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
            mass = mw_from_smiles(smiles)
            sample_data[samples] = {"mass": mass}
        db_data = True
    elif sample_data_file:
        sample_data = _get_sample_data(sample_data_file)
        db_data = False

    return sample_data, db_data


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


def lcms_to_compounds(config, table_data):

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


def name_changer(new_names, raw_data_dict, sample_data, peak_information, sample_peak_dict):

    for names in new_names:
        if names != new_names[names]["raw"]:
            raw_data_dict[names] = raw_data_dict.pop(new_names[names]["raw"])
            peak_information[names] = peak_information.pop(new_names[names]["raw"])
            sample_peak_dict[names] = sample_peak_dict.pop(new_names[names]["raw"])
        elif names != new_names[names]["excel"]:
            sample_data[names] = sample_data.pop(new_names[names]["excel"])
