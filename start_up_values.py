import os

from info import clm_to_row_96, clm_to_row_384, clm_to_row_1536, row_to_clm_1536, row_to_clm_384, row_to_clm_96

all_table_data_extra = {
    "-COMPOUND_INFO_INFO_MP_TABLE-": {"name": "Assay Table",
                                      "headings": ["plate", "well", "vol"]},
    "-COMPOUND_INFO_INFO_DP_TABLE-": {"name": "Plate Table",
                                      "headings": ["Plate", "well", "Vol"]},
    "-COMPOUND_INFO_INFO_ASSAY_TABLE-": {"name": "Assay Compound Table",
                                         "headings": ["Assay", "run", "Plate", "Well", "Score", "Approved"]},
    "-COMPOUND_INFO_INFO_HITS_TABLE-": {"name": "",
                                        "headings": ["Assay", "Score", "Conc."]},
    "-COMPOUND_INFO_INFO_PURITY_USED_TABLE-": {"name": "",
                                               "headings": ["Purity", "Replicates", "Date"]},
}

plate_type_count = {"plate_96": 96, "plate_384": 384, "plate_1536": 1536}
row_to_clm_converter = {"plate_96": row_to_clm_96, "plate_384": row_to_clm_384, "plate_1536": row_to_clm_1536}
clm_to_row_converter = {"plate_96": clm_to_row_96, "plate_384": clm_to_row_384, "plate_1536": clm_to_row_1536}
ms_mode_selector = {"Positive": "ms_pos", "Negative": "ms_neg"}
search_reverse = {}

bio_info_clicked = False

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

all_data = None

# BIO EXP TABLE CONSTANTS:
window_tables = {"all_assays": None, "plate_bio_info": None}
bio_info_tables = ["-BIO_EXP_TABLE_ASSAY_LIST_BOX-+-double click-",
                   "-BIO_EXP_ASSAY_RUN_TABLE-+-double click-",
                   "-BIO_EXP_PLATE_TABLE-+-double click-"]

colour_chooser_buttons = ["-BIO_INFO_HIT_MAP_TH_1_COLOUR_TARGET-",
                          "-BIO_INFO_HIT_MAP_TH_2_COLOUR_TARGET-",
                          "-BIO_INFO_HIT_MAP_TH_3_COLOUR_TARGET-",
                          "-BIO_INFO_HIT_MAP_TH_4_COLOUR_TARGET-",
                          "-BIO_INFO_HIT_MAP_TH_5_COLOUR_TARGET-",
                          "-BIO_INFO_HIT_MAP_TH_6_COLOUR_TARGET-",
                          "-BIO_INFO_HIT_MAP_TH_7_COLOUR_TARGET-",
                          "-BIO_INFO_HIT_MAP_TH_8_COLOUR_TARGET-",
                          "-BIO_INFO_HIT_MAP_TH_9_COLOUR_TARGET-",
                          "-BIO_INFO_HEATMAP_LOW_COLOUR_TARGET-",
                          "-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-",
                          "-BIO_INFO_HEATMAP_HIGH_COLOUR_TARGET-"]

assay_updater_list = ["-BIO_INFO_ASSAY_DROPDOWN-",
                      "-BIO_INFO_APPROVED_CHECK-",
                      "-BIO_INFO_ANALYSE_METHOD-"]

compound_info_tables = ["-COMPOUND_INFO_INFO_MP_TABLE-+-double click-",
                        "-COMPOUND_INFO_INFO_DP_TABLE-+-double click-",
                        "-COMPOUND_INFO_INFO_ASSAY_TABLE-+-double click-",
                        "-COMPOUND_INFO_INFO_HITS_TABLE-+-double click-",
                        "-COMPOUND_INFO_INFO_PURITY_USED_TABLE-+-double click-"]


def database_guard(config, cw):
    try:
        if os.path.exists(config["Database"]["database"]):
            if config["Database"]["database"] is not None:
                db_active = True
            else:
                db_active = False
        else:
            cw.delete_all_info("Database")
            db_active = False
    except KeyError:
        db_active = False
    return db_active


def start_up_gui(config, window):
    well_dict_bio_info = {}
    well_dict = {}
    dose_colour_dict = {}
    colour_select = {}
    sub_search_info = {}
    for keys in list(config["plate_colouring"].keys()):
        colour_select[keys] = config["plate_colouring"][keys]

    graph_bio_exp = window["-BIO_INFO_CANVAS-"]
    window_1_plate_layout["graph_plate"] = window["-RECT_BIO_CANVAS-"]

    #   WINDOW 2 - PURITY INFO  #
    lc_graph_showing = [keys for keys in list(config["lc_mapping"].keys())]

    # Makes it possible to double-click on the table

    window["-COMPOUND_INFO_INFO_MP_TABLE-"].bind('<Double-Button-1>', "+-double click-")
    window["-COMPOUND_INFO_INFO_DP_TABLE-"].bind('<Double-Button-1>', "+-double click-")
    window["-COMPOUND_INFO_INFO_ASSAY_TABLE-"].bind('<Double-Button-1>', "+-double click-")
    window["-COMPOUND_INFO_INFO_HITS_TABLE-"].bind('<Double-Button-1>', "+-double click-")
    window["-COMPOUND_INFO_INFO_PURITY_USED_TABLE-"].bind('<Double-Button-1>', "+-double click-")
    window["-MAIN_COMPOUND_TABLE-"].bind('<Double-Button-1>', "+-double click-")

    # Search Menu
    window["-SUB_SEARCH_TABLE-"].bind('<Double-Button-1>', "+-double click-")

    # Plate Tables
    window["-PLATE_TABLE_TABLE-"].bind('<Double-Button-1>', "+-double click-")

    # Bio tables
    window["-BIO_EXP_TABLE_ASSAY_LIST_BOX-"].bind('<Double-Button-1>', "+-double click-")
    window["-BIO_EXP_ASSAY_RUN_TABLE-"].bind('<Double-Button-1>', "+-double click-")
    window["-BIO_EXP_PLATE_TABLE-"].bind('<Double-Button-1>', "+-double click-")
    window["-BIO_EXP_COMPOUND_TABLE-"].bind('<Double-Button-1>', "+-double click-")

    # LC tables
    window["-LC_MS_SAMPLE_TABLE-"].bind('<Double-Button-1>', "+-double click-")

    # Bio Analyser

    # Table stuff
    window.Element("-BIO_INFO_MATRIX_TABLE-").Widget.configure(displaycolumns=[])
    window.Element("-PLATE_TABLE_TABLE-").Widget.configure(displaycolumns=plate_table_table_heading_mp)

    return well_dict_bio_info, well_dict, dose_colour_dict, colour_select, graph_bio_exp, lc_graph_showing, \
           sub_search_info
