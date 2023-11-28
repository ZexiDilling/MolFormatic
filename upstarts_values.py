import os


from database_handler import DataBaseFunctions
from database_functions import grab_table_data
from info import clm_to_row_96, clm_to_row_384, clm_to_row_1536, row_to_clm_1536, row_to_clm_384, row_to_clm_96

# Makes a dict over all tables in the software. It is used to make tables sortable.
# Any tables in this dict, will be sortable by clicking the top bar.
all_table_data = {"-COMPOUND_INFO_PLATE_TABLE-": None,
                  "-BIO_INFO_OVERVIEW_TABLE-": None,
                  "-BIO_INFO_OVERVIEW_AVG_TABLE-": None,
                  "-BIO_INFO_OVERVIEW_STDEV_TABLE-": None,
                  "-BIO_INFO_OVERVIEW_Z_PRIME_TABLE-": None,
                  "-BIO_INFO_Z_PRIME_LIST_TABLE-": None,
                  "-BIO_INFO_HIT_LIST_LOW_TABLE-": None,
                  "-BIO_INFO_HIT_LIST_MID_TABLE-": None,
                  "-BIO_INFO_HIT_LIST_HIGH_TABLE-": None,
                  "-BIO_INFO_MATRIX_TABLE-": None,
                  "-PURITY_INFO_OVERVIEW_TABLE-": None,
                  "-PURITY_INFO_PURITY_OVERVIEW_TABLE-": None,
                  "-PURITY_INFO_PEAK_TABLE-": None,
                  "-PURITY_INFO_PURITY_PEAK_LIST_TABLE-": None,
                  "-BIO_EXP_PLATE_TABLE-": None,
                  "-BIO_EXP_COMPOUND_TABLE-": None,
                  "-LC_MS_SAMPLE_TABLE-": None,
                  "-PLATE_TABLE_TABLE-": None,
                  "-PURITY_INFO_RAW_DATA_TABLE-": None,
                  "-COMPOUND_INFO_ALL_PLATE_INFO_TABLE-": None,
                  "-COMPOUND_INFO_MP_PLATE_INFO_TABLE-": None,
                  "-COMPOUND_INFO_DP_PLATE_INFO_TABLE-": None,
                  "-COMPOUND_INFO_BIO_INFO_TABLE-": None,
                  "-COMPOUND_INFO_PURITY_INFO_TABLE-": None,
                  "-EXTRA_DATABASE_CUSTOMERS_TABLE-": None,
                  "-EXTRA_DATABASE_VENDORS_TABLE-": None,
                  "-EXTRA_DATABASE_AC_TABLE-": None,
                  "-EXTRA_DATABASE_RESPONSIBLE_TABLE-": None,
                  "-EXTRA_DATABASE_LOCATIONS_TABLE-": None,
                  }


plate_type_count = {"plate_96": 96, "plate_384": 384, "plate_1536": 1536}
row_to_clm_converter = {"plate_96": row_to_clm_96, "plate_384": row_to_clm_384, "plate_1536": row_to_clm_1536}
clm_to_row_converter = {"plate_96": clm_to_row_96, "plate_384": clm_to_row_384, "plate_1536": clm_to_row_1536}
ms_mode_selector = {"Positive": "ms_pos", "Negative": "ms_neg"}

draw_tool_values = {
    "start_point": None,
    "end_point": None,
    "prior_rect": None,
    "x": None,
    "y": None,
    "min_x": None,
    "max_x": None,
    "min_y": None,
    "max_y": None,
    "dragging": False,
    "temp_draw_tool": None,
    "temp_selector": False,
    "well_off_set": 0,
    }

bio_info_set_settings = {
    "bio_info_sub_setting_tab_mapping_calc": False,
    "bio_info_sub_setting_tab_matrix_calc": False,
    "bio_info_sub_setting_tab_list_calc": False,
    "bio_info_sub_setting_tab_z_prime_calc": False,
    "bio_info_sub_setting_tab_plate_overview_calc": False,
    "bio_info_sub_setting_tab_overview_calc": False,
    "bio_info_sub_setting_tab_hit_list_calc": False,
}

# PLATE TABLE CONSTANTS #
plate_table_table_heading_mp = ["Barcode", "Compound", "Well", "Volume", "Date", "Active",
                                "Freeze/Thaw", "Plate Type", "location"]
plate_table_table_heading_dp = ["Barcode", "Compound", "Well", "Volume", "Date", "Active", "Freeze/Thaw",
                                "Plate Type", "location", "Source Plate", "Source Well"]

window_1_search = {
    "ac_use": False,
    "origin_use": False,
    "transferee_volume": None,
    "compound_table_clear": False,
    "current_table_data": None
}

window_1_bio = {
    "bio_export_folder": None
}

window_1_plate_layout = {
    "plate_type": None,
    "plate_active": False,
    "temp_draw_tool": "sample",
    "total_sample_spots": 0,
    "graph_plate": 0,
    "temp_draw_tool_tracker": "-RECT_SAMPLES-",
    "draw_tool_dict": {
        "-RECT_SAMPLES-": "sample",
        "-RECT_BLANK-": "blank",
        "-RECT_MAX-": "max",
        "-RECT_MIN-": "minimum",
        "-RECT_NEG-": "negative",
        "-RECT_POS-": "positive",
        "-RECT_EMPTY-": "empty",
        "-COLOUR-": "paint",
        "-RECT_DOSE-": "dose"
    }
}

window_1_lcms = {
    "purity_data": None,
    "purity_data_added_to_db": False,
    "update_purity_info_peak_table": True,
    "canvas_lines": {
       "uv": None,
       "peak_lines": {}
    },
    "sample_data_file": None,
    "plot_style": None,
    "toolbar": None,
    "purity_info_rt_start": None,
    "purity_info_rt_end": None,
    "temp_purity_info_canvas": None,
    "purity_info_values": False,
    "purity_info_samples": None}

window_1_extra = {
    "customers_data": None,
    "vendors_data": None,
    "ac_data": None,
    "plate_types_data": None,
    "responsible_data": None,
    "location_data": None,
}

# BIO EXP TABLE CONSTANTS:
window_tables = {"all_assays": None, "plate_bio_info": None}


def database_guard(config, cw):
    try:
        if os.path.exists(config["Database"]["database"]):
            db_active = True
        else:
            cw.delete_all_info("Database")
            db_active = False
    except KeyError:
        db_active = False
    return db_active


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
        plate_type = temp_row_data[0][2]

        temp_sub_row_data = dbf.find_data_single_lookup("plate_layout_sub", plate_name, "plate_sub")
        try:
            well_layout = eval(temp_sub_row_data[0][3])
        except IndexError:
            well_layout = None

        try:
            sample_type = temp_sub_row_data[0][4]
        except IndexError:
            sample_type = None

        # Generates the dict
        archive_plates_dict[plate_name] = {"well_layout": well_layout,
                                           "plate_type": plate_type,
                                           "sample_type": sample_type,
                                           "sample": [],
                                           "blank": [],
                                           "max": [],
                                           "minimum": [],
                                           "positive": [],
                                           "negative": [],
                                           "empty": []}

        # Makes list of the different well-types and adds them
        if well_layout is not None:
            temp_well_dict = well_layout
            for counter in temp_well_dict:
                well_id = temp_well_dict[counter]["well_id"]
                if temp_well_dict[counter]["state"] == "sample":
                    archive_plates_dict[plate_name]["sample"].append(well_id)
                elif temp_well_dict[counter]["state"] == "blank":
                    archive_plates_dict[plate_name]["blank"].append(well_id)
                elif temp_well_dict[counter]["state"] == "max":
                    archive_plates_dict[plate_name]["max"].append(well_id)
                elif temp_well_dict[counter]["state"] == "minimum":
                    archive_plates_dict[plate_name]["minimum"].append(well_id)
                elif temp_well_dict[counter]["state"] == "positive":
                    archive_plates_dict[plate_name]["positive"].append(well_id)
                elif temp_well_dict[counter]["state"] == "negative":
                    archive_plates_dict[plate_name]["negative"].append(well_id)
                elif temp_well_dict[counter]["state"] == "empty":
                    archive_plates_dict[plate_name]["empty"].append(well_id)

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


def get_assay_list(config):
    assay_table_data = grab_table_data(config, "assay")
    assay_list = []
    for row in assay_table_data[0]:
        assay_list.append(row[1])

    return assay_list


def gui_data_fetcher(db_active, config):
    if db_active:
        plate_list, archive_plates_dict = get_plate_layout(config)
        assay_list = get_assay_list(config)
    else:
        archive_plates_dict = {}
        plate_list = []
        assay_list = []

    return plate_list, assay_list, archive_plates_dict



