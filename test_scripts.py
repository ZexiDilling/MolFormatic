import copy
import os.path
import configparser

from gui_layout import GUILayout
from gui_settings_control import GUISettingsController
from gui_functions import *
from bio_data_functions import org, norm, pora, pora_internal
from json_handler import plate_dict_reader, dict_writer, dict_reader
from plate_formatting import plate_layout_re_formate
from gui_popup import matrix_popup
from gui_help_info_controller import help_info_controller
from config_writer import ConfigWriter
from database_startup import DatabaseSetUp



def main(config):
    """
    The main control modul for the GUI.

    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :return: This is a gui control modul. Returns depending on what is being done in the GUI
    """

    # importing config writer, to write data to the config file
    cw = ConfigWriter(config)
    try:
        if os.path.exists(config["Database"]["database"]):
            db_active = True
        else:
            cw.delete_all_info("Database")
            db_active = False
    except KeyError:
        db_active = False



    # File names, for files with dict over different kind of data.
    plate_file = config["files"]["plate_layouts"]
    bio_files = config["files"]["bio_experiments"]

    try:
        plate_list, archive_plates_dict = plate_dict_reader(plate_file)
    except TypeError:
        plate_list = []
        archive_plates_dict = {}

    gl = GUILayout(config, plate_list)

    window = gl.full_layout()
    window.maximize()

    #   Deaful Values / Simple Settings #
    simple_settings = {
        "plate_colouring": {
            "sample": config["plate_colouring"]["sample"],
            "blank": config["plate_colouring"]["blank"],
            "max": config["plate_colouring"]["max"],
            "minimum": config["plate_colouring"]["minimum"],
            "positive": config["plate_colouring"]["positive"],
            "negative": config["plate_colouring"]["negative"],
            "empty": config["plate_colouring"]["empty"]
        }
    }
    #   WINDOW 1 - SEARCH   #
    ac_use = False
    origin_use = False
    transferee_volume = None
    compound_table_clear = False
    current_table_data = None

    #   WINDOW 1 - BIO  #
    graph_bio = window["-BIO_CANVAS-"]
    bio_export_folder = None
    bio_final_report_setup = {
        "methods": {"original": config["Settings_bio"].getboolean("final_report_methods_original"),
                    "normalised": config["Settings_bio"].getboolean("final_report_methods_normalised"),
                    "pora": config["Settings_bio"].getboolean("final_report_methods_pora")},
        "analyse": {"sample": config["Settings_bio"].getboolean("final_report_analyse_sample"),
                    "minimum": config["Settings_bio"].getboolean("final_report_analyse_minimum"),
                    "max": config["Settings_bio"].getboolean("final_report_analyse_max"),
                    "empty": config["Settings_bio"].getboolean("final_report_analyse_empty"),
                    "negative": config["Settings_bio"].getboolean("final_report_analyse_negative"),
                    "positive": config["Settings_bio"].getboolean("final_report_analyse_positive"),
                    "blank": config["Settings_bio"].getboolean("final_report_analyse_blank")},
        "calc": {"original": {"overview": config["Settings_bio"].getboolean("final_report_calc_original_overview"),
                              "sample": config["Settings_bio"].getboolean("final_report_calc_original_sample"),
                              "minimum": config["Settings_bio"].getboolean("final_report_calc_original_minimum"),
                              "max": config["Settings_bio"].getboolean("final_report_calc_original_max"),
                              "empty": config["Settings_bio"].getboolean("final_report_calc_original_empty"),
                              "negative": config["Settings_bio"].
                                  getboolean("final_report_calc_original_negative"),
                              "positive": config["Settings_bio"].
                                  getboolean("final_report_calc_original_positive"),
                              "blank": config["Settings_bio"].getboolean("final_report_calc_original_blank")},
                 "normalised": {"overview": config["Settings_bio"].getboolean("final_report_calc_normalised_overview"),
                                "sample": config["Settings_bio"].getboolean("final_report_calc_normalised_sample"),
                                "minimum": config["Settings_bio"].getboolean("final_report_calc_normalised_minimum"),
                                "max": config["Settings_bio"].getboolean("final_report_calc_normalised_max"),
                                "empty": config["Settings_bio"].getboolean("final_report_calc_normalised_empty"),
                                "negative": config["Settings_bio"].
                                    getboolean("final_report_calc_normalised_negative"),
                                "positive": config["Settings_bio"].
                                    getboolean("final_report_calc_normalised_positive"),
                                "blank": config["Settings_bio"].getboolean("final_report_calc_normalised_blank")},
                 "pora": {"overview": config["Settings_bio"].getboolean("final_report_calc_pora_overview"),
                          "sample": config["Settings_bio"].getboolean("final_report_calc_pora_sample"),
                          "minimum": config["Settings_bio"].getboolean("final_report_calc_pora_minimum"),
                          "max": config["Settings_bio"].getboolean("final_report_calc_pora_max"),
                          "empty": config["Settings_bio"].getboolean("final_report_calc_pora_empty"),
                          "negative": config["Settings_bio"].getboolean("final_report_calc_pora_negative"),
                          "positive": config["Settings_bio"].getboolean("final_report_calc_pora_positive"),
                          "blank": config["Settings_bio"].getboolean("final_report_calc_pora_blank")},
                 "z_prime": config["Settings_bio"].getboolean("final_report_calc_Z_prime")},
        "pora_threshold": {"low": {"min": config["Settings_bio"].getfloat("final_report_pora_threshold_low_min"),
                                   "max": config["Settings_bio"].getfloat("final_report_pora_threshold_low_max")},
                           "mid": {"min": config["Settings_bio"].getfloat("final_report_pora_threshold_mid_min"),
                                   "max": config["Settings_bio"].getfloat("final_report_pora_threshold_mid_max")},
                           "high": {"min": config["Settings_bio"].getfloat("final_report_pora_threshold_high_min"),
                                    "max": config["Settings_bio"].getfloat("final_report_pora_threshold_high_max")}},
        "data": {"sample": {"matrix": config["Settings_bio"].getboolean("final_report_data_sample_matrix"),
                            "list": config["Settings_bio"].getboolean("final_report_data_sample_list"),
                            "max_min": config["Settings_bio"].getboolean("final_report_data_sample_max_min")},
                 "minimum": {"matrix": config["Settings_bio"].getboolean("final_report_data_minimum_matrix"),
                             "list": config["Settings_bio"].getboolean("final_report_data_minimum_list"),
                             "max_min": config["Settings_bio"].getboolean("final_report_data_minimum_max_min")},
                 "max": {"matrix": config["Settings_bio"].getboolean("final_report_data_max_matrix"),
                         "list": config["Settings_bio"].getboolean("final_report_data_max_list"),
                         "max_min": config["Settings_bio"].getboolean("final_report_data_max_max_min")},
                 "empty": {"matrix": config["Settings_bio"].getboolean("final_report_data_empty_matrix"),
                           "list": config["Settings_bio"].getboolean("final_report_data_empty_list"),
                           "max_min": config["Settings_bio"].getboolean("final_report_data_empty_max_min")},
                 "negative": {"matrix": config["Settings_bio"].getboolean("final_report_data_negative_matrix"),
                              "list": config["Settings_bio"].getboolean("final_report_data_negative_list"),
                              "max_min": config["Settings_bio"].getboolean("final_report_data_negative_max_min")},
                 "positive": {"matrix": config["Settings_bio"].getboolean("final_report_data_positive_matrix"),
                              "list": config["Settings_bio"].getboolean("final_report_data_positive_list"),
                              "max_min": config["Settings_bio"].getboolean("final_report_data_positive_max_min")},
                 "blank": {"matrix": config["Settings_bio"].getboolean("final_report_data_blank_matrix"),
                           "list": config["Settings_bio"].getboolean("final_report_data_blank_list"),
                           "max_min": config["Settings_bio"].getboolean("final_report_data_blank_max_min")},
                 "z_prime": {"matrix": config["Settings_bio"].getboolean("final_report_data_z_prime_matrix"),
                             "list": config["Settings_bio"].getboolean("final_report_data_z_prime_list"),
                             "max_min": config["Settings_bio"].getboolean("final_report_data_z_prime_max_min")}}}

    bio_plate_report_setup = {
        "well_states_report_method": {"original": config["Settings_bio"].
            getboolean("well_states_report_method_original"),
                                      "normalised": config["Settings_bio"].
                                          getboolean("well_states_report_method_normalised"),
                                      "pora": config["Settings_bio"].getboolean("well_states_report_method_pora"),
                                      "pora_internal": config["Settings_bio"].
                                          getboolean("well_states_report_method_pora_internal")},
        "well_states_report": {'sample': config["Settings_bio"].getboolean("plate_report_well_states_report_sample"),
                               'blank': config["Settings_bio"].getboolean("plate_report_well_states_report_blank"),
                               'max': config["Settings_bio"].getboolean("plate_report_well_states_report_max"),
                               'minimum': config["Settings_bio"].getboolean("plate_report_well_states_report_minimum"),
                               'positive': config["Settings_bio"].getboolean("plate_report_well_states_report_positive")
            ,
                               'negative': config["Settings_bio"].getboolean("plate_report_well_states_report_negative")
            ,
                               'empty': config["Settings_bio"].getboolean("plate_report_well_states_report_empty")},
        "calc_dict": {"original": {"use": config["Settings_bio"].getboolean("plate_report_calc_dict_original_use"),
                                   "avg": config["Settings_bio"].getboolean("plate_report_calc_dict_original_avg"),
                                   "stdev": config["Settings_bio"].getboolean("plate_report_calc_dict_original_stdev"),
                                   "state": {"sample": config["Settings_bio"].
                                       getboolean("plate_report_calc_dict_original_state_sample"),
                                             "minimum": config["Settings_bio"].
                                                 getboolean("plate_report_calc_dict_original_state_minimum"),
                                             "max": config["Settings_bio"].
                                                 getboolean("plate_report_calc_dict_original_state_max"),
                                             "empty": config["Settings_bio"].
                                                 getboolean("plate_report_calc_dict_original_state_empty"),
                                             "negative": config["Settings_bio"].
                                                 getboolean("plate_report_calc_dict_original_state_negative"),
                                             "positive": config["Settings_bio"].
                                                 getboolean("plate_report_calc_dict_original_state_positive"),
                                             "blank": config["Settings_bio"].
                                                 getboolean("plate_report_calc_dict_original_state_blank")}},
                      "normalised": {"use": config["Settings_bio"].getboolean("plate_report_calc_dict_normalised_use"),
                                     "avg": config["Settings_bio"].
                                         getboolean("plate_report_calc_dict_normalised_avg"),
                                     "stdev": config["Settings_bio"].
                                         getboolean("plate_report_calc_dict_normalised_stdev"),
                                     "state": {"sample": config["Settings_bio"].
                                         getboolean("plate_report_calc_dict_normalised_"
                                                    "state_sample"),
                                               "minimum": config["Settings_bio"].
                                                   getboolean("plate_report_calc_dict_normalised_"
                                                              "state_minimum"),
                                               "max": config["Settings_bio"].
                                                   getboolean("plate_report_calc_dict_normalised_"
                                                              "state_max"),
                                               "empty": config["Settings_bio"].
                                                   getboolean("plate_report_calc_dict_normalised_"
                                                              "state_empty"),
                                               "negative": config["Settings_bio"].
                                                   getboolean("plate_report_calc_dict_normalised_"
                                                              "state_negative"),
                                               "positive": config["Settings_bio"].
                                                   getboolean("plate_report_calc_dict_normalised_"
                                                              "state_positive"),
                                               "blank": config["Settings_bio"].
                                                   getboolean("plate_report_calc_dict_normalised_"
                                                              "state_blank")}},
                      "pora": {"use": config["Settings_bio"].getboolean("plate_report_calc_dict_pora_use"),
                               "avg": config["Settings_bio"].getboolean("plate_report_calc_dict_pora_avg"),
                               "stdev": config["Settings_bio"].getboolean("plate_report_calc_dict_pora_stdev"),
                               "state": {"sample": config["Settings_bio"].
                                   getboolean("plate_report_calc_dict_pora_state_sample"),
                                         "minimum": config["Settings_bio"].
                                             getboolean("plate_report_calc_dict_pora_state_minimum"),
                                         "max": config["Settings_bio"].getboolean(
                                             "plate_report_calc_dict_pora_state_max"),
                                         "empty": config["Settings_bio"].
                                             getboolean("plate_report_calc_dict_pora_state_empty"),
                                         "negative": config["Settings_bio"].
                                             getboolean("plate_report_calc_dict_pora_state_negative"),
                                         "positive": config["Settings_bio"].
                                             getboolean("plate_report_calc_dict_pora_state_positive"),
                                         "blank": config["Settings_bio"].
                                             getboolean("plate_report_calc_dict_pora_state_blank")}},
                      "pora_internal": {"use": config["Settings_bio"].
                          getboolean("plate_report_calc_dict_pora_internal_use"),
                                        "avg": config["Settings_bio"].
                                            getboolean("plate_report_calc_dict_pora_internal_avg"),
                                        "stdev": config["Settings_bio"].
                                            getboolean("plate_report_calc_dict_pora_internal_stdev"),
                                        "state": {"sample": config["Settings_bio"].
                                            getboolean("plate_report_calc_dict_pora_internal_state_sample"),
                                                  "minimum": config["Settings_bio"].
                                                      getboolean("plate_report_calc_dict_pora_internal_state_minimum"),
                                                  "max": config["Settings_bio"].
                                                      getboolean("plate_report_calc_dict_pora_internal_state_max"),
                                                  "empty": config["Settings_bio"].
                                                      getboolean("plate_report_calc_dict_pora_internal_state_empty"),
                                                  "negative": config["Settings_bio"].
                                                      getboolean("plate_report_calc_dict_pora_internal_state_negative"),
                                                  "positive": config["Settings_bio"].
                                                      getboolean("plate_report_calc_dict_pora_internal_state_positive"),
                                                  "blank": config["Settings_bio"].
                                                      getboolean("plate_report_calc_dict_pora_internal_state_blank")}},
                      "other": {"use": config["Settings_bio"].getboolean("plate_report_calc_dict_other_use"),
                                "calc": {"z_prime": config["Settings_bio"].
                                    getboolean("plate_report_calc_dict_other_calc_z_prime")}}},
        "plate_calc_dict": {
            "original": {"use": config["Settings_bio"].getboolean("plate_report_plate_calc_dict_original_use"),
                         "avg": config["Settings_bio"].getboolean("plate_report_plate_calc_dict_original_avg"),
                         "stdev": config["Settings_bio"].getboolean("plate_report_plate_calc_dict_original_stdev"),
                         "state": {"sample": config["Settings_bio"].
                             getboolean("plate_report_plate_calc_dict_original_state_sample"),
                                   "minimum": config["Settings_bio"].
                                       getboolean("plate_report_plate_calc_dict_original_state_minimum"),
                                   "max": config["Settings_bio"].
                                       getboolean("plate_report_plate_calc_dict_original_state_max"),
                                   "empty": config["Settings_bio"].
                                       getboolean("plate_report_plate_calc_dict_original_state_empty"),
                                   "negative": config["Settings_bio"].
                                       getboolean("plate_report_plate_calc_dict_original_state_negative"),
                                   "positive": config["Settings_bio"].
                                       getboolean("plate_report_plate_calc_dict_original_state_positive"),
                                   "blank": config["Settings_bio"].
                                       getboolean("plate_report_plate_calc_dict_original_state_blank")}},
            "normalised": {"use": config["Settings_bio"].getboolean("plate_report_plate_calc_dict_normalised_use"),
                           "avg": config["Settings_bio"].getboolean("plate_report_plate_calc_dict_normalised_avg"),
                           "stdev": config["Settings_bio"].getboolean("plate_report_plate_calc_dict_normalised_stdev"),
                           "state": {"sample": config["Settings_bio"].
                               getboolean("plate_report_plate_calc_dict_normalised_state_sample"),
                                     "minimum": config["Settings_bio"].
                                         getboolean("plate_report_plate_calc_dict_normalised_state_minimum"),
                                     "max": config["Settings_bio"].
                                         getboolean("plate_report_plate_calc_dict_normalised_state_max"),
                                     "empty": config["Settings_bio"].
                                         getboolean("plate_report_plate_calc_dict_normalised_state_empty"),
                                     "negative": config["Settings_bio"].
                                         getboolean("plate_report_plate_calc_dict_normalised_state_negative"),
                                     "positive": config["Settings_bio"].
                                         getboolean("plate_report_plate_calc_dict_normalised_state_positive"),
                                     "blank": config["Settings_bio"].
                                         getboolean("plate_report_plate_calc_dict_normalised_state_blank")}},
            "pora": {"use": config["Settings_bio"].getboolean("plate_report_plate_calc_dict_pora_use"),
                     "avg": config["Settings_bio"].getboolean("plate_report_plate_calc_dict_pora_avg"),
                     "stdev": config["Settings_bio"].getboolean("plate_report_plate_calc_dict_pora_stdev"),
                     "state": {"sample": config["Settings_bio"].
                         getboolean("plate_report_plate_calc_dict_pora_state_sample"),
                               "minimum": config["Settings_bio"].
                                   getboolean("plate_report_plate_calc_dict_pora_state_minimum"),
                               "max": config["Settings_bio"].
                                   getboolean("plate_report_plate_calc_dict_pora_state_max"),
                               "empty": config["Settings_bio"].
                                   getboolean("plate_report_plate_calc_dict_pora_state_empty"),
                               "negative": config["Settings_bio"].
                                   getboolean("plate_report_plate_calc_dict_pora_state_negative"),
                               "positive": config["Settings_bio"].
                                   getboolean("plate_report_plate_calc_dict_pora_state_positive"),
                               "blank": config["Settings_bio"].
                                   getboolean("plate_report_plate_calc_dict_pora_state_blank")}},
            "pora_internal": {"use": config["Settings_bio"].
                getboolean("plate_report_plate_calc_dict_pora_internal_use"),
                              "avg": config["Settings_bio"].
                                  getboolean("plate_report_plate_calc_dict_pora_internal_avg"),
                              "stdev": config["Settings_bio"].
                                  getboolean("plate_report_plate_calc_dict_pora_internal_stdev"),
                              "state": {"sample": config["Settings_bio"].
                                  getboolean("plate_report_plate_calc_dict_pora_internal_state_sample"),
                                        "minimum": config["Settings_bio"].
                                            getboolean("plate_report_plate_calc_dict_pora_internal_state_minimum"),
                                        "max": config["Settings_bio"].
                                            getboolean("plate_report_plate_calc_dict_pora_internal_state_max"),
                                        "empty": config["Settings_bio"].
                                            getboolean("plate_report_plate_calc_dict_pora_internal_state_empty"),
                                        "negative": config["Settings_bio"].
                                            getboolean("plate_report_plate_calc_dict_pora_internal_state_negative"),
                                        "positive": config["Settings_bio"].
                                            getboolean("plate_report_plate_calc_dict_pora_internal_state_positive"),
                                        "blank": config["Settings_bio"].
                                            getboolean("plate_report_plate_calc_dict_pora_internal_state_blank")}},
        },
        "plate_analysis_dict": {"original": {"use": config["Settings_bio"].
            getboolean("plate_report_plate_analysis_dict_original_use"),
                                             "methode": org,
                                             "state_map": config["Settings_bio"].
                                                 getboolean("plate_report_plate_analysis_dict_original_state_map"),
                                             "heatmap": config["Settings_bio"].
                                                 getboolean("plate_report_plate_analysis_dict_original_heatmap"),
                                             "hit_map": config["Settings_bio"].
                                                 getboolean("plate_report_plate_analysis_dict_original_hit_map"),
                                             "none": config["Settings_bio"].
                                                 getboolean("plate_report_plate_analysis_dict_original_none")},
                                "normalised": {"use": config["Settings_bio"].
                                    getboolean("plate_report_plate_analysis_dict_normalised_use"),
                                               "methode": norm,
                                               "state_map": config["Settings_bio"].
                                                   getboolean("plate_report_plate_analysis_dict_normalised_state_map"),
                                               "heatmap": config["Settings_bio"].
                                                   getboolean("plate_report_plate_analysis_dict_normalised_heatmap"),
                                               "hit_map": config["Settings_bio"].
                                                   getboolean("plate_report_plate_analysis_dict_normalised_hit_map"),
                                               "none": config["Settings_bio"].
                                                   getboolean("plate_report_plate_analysis_dict_normalised_none")},
                                "pora": {"use": config["Settings_bio"].
                                    getboolean("plate_report_plate_analysis_dict_pora_use"),
                                         "methode": pora,
                                         "state_map": config["Settings_bio"].
                                             getboolean("plate_report_plate_analysis_dict_pora_state_map"),
                                         "heatmap": config["Settings_bio"].
                                             getboolean("plate_report_plate_analysis_dict_pora_heatmap"),
                                         "hit_map": config["Settings_bio"].
                                             getboolean("plate_report_plate_analysis_dict_pora_hit_map"),
                                         "none": config["Settings_bio"].
                                             getboolean("plate_report_plate_analysis_dict_pora_none")},
                                "pora_internal": {"use": config["Settings_bio"].
                                    getboolean("plate_report_plate_analysis_dict_pora_internal_use"),
                                                  "methode": pora_internal,
                                                  "state_map": config["Settings_bio"].
                                                      getboolean(
                                                      "plate_report_plate_analysis_dict_pora_internal_state_map")
                                    ,
                                                  "heatmap": config["Settings_bio"].
                                                      getboolean(
                                                      "plate_report_plate_analysis_dict_pora_internal_heatmap"),
                                                  "hit_map": config["Settings_bio"].
                                                      getboolean(
                                                      "plate_report_plate_analysis_dict_pora_internal_hit_map"),
                                                  "none": config["Settings_bio"].
                                                      getboolean("plate_report_plate_analysis_dict_pora_internal_none")}
                                },
        "z_prime_calc": config["Settings_bio"].getboolean("plate_report_z_prime_calc"),
        "heatmap_colours": {'low': config["Settings_bio"]["plate_report_heatmap_colours_low"],
                            'mid': config["Settings_bio"]["plate_report_heatmap_colours_mid"],
                            'high': config["Settings_bio"]["plate_report_heatmap_colours_high"]},
        "pora_threshold": {"low": {"min": config["Settings_bio"].getfloat("plate_report_pora_threshold_low_min"),
                                   "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_low_max")},
                           "mid": {"min": config["Settings_bio"].getfloat("plate_report_pora_threshold_mid_min"),
                                   "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_mid_max")},
                           "high": {"min": config["Settings_bio"].getfloat("plate_report_pora_threshold_high_min"),
                                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_high_max")},
                           "colour": {"low": config["Settings_bio"]["plate_report_pora_threshold_colour_low"],
                                      "mid": config["Settings_bio"]["plate_report_pora_threshold_colour_mid"],
                                      "high": config["Settings_bio"]["plate_report_pora_threshold_colour_high"]}
                           }
    }

    #   WINDOW 1 - PURITY DATA
    ms_settings = {
        "ions": {"positive": {
            "m+3h": bool(config["Positive ion mode"]["m+3h"].split(",")[-1]),
            "m+2h+na": bool(config["Positive ion mode"]["m+2h+na"].split(",")[-1]),
            "m+h+2na": bool(config["Positive ion mode"]["m+h+2na"].split(",")[-1]),
            "m+3na": bool(config["Positive ion mode"]["m+3na"].split(",")[-1]),
            "m+2h": bool(config["Positive ion mode"]["m+2h"].split(",")[-1]),
            "m+h+nh4": bool(config["Positive ion mode"]["m+h+nh4"].split(",")[-1]),
            "m+h+na": bool(config["Positive ion mode"]["m+h+na"].split(",")[-1]),
            "m+h+k": bool(config["Positive ion mode"]["m+h+k"].split(",")[-1]),
            "m+acn+2h": bool(config["Positive ion mode"]["m+acn+2h"].split(",")[-1]),
            "m+2na": bool(config["Positive ion mode"]["m+2na"].split(",")[-1]),
            "m+2acn+2h": bool(config["Positive ion mode"]["m+2acn+2h"].split(",")[-1]),
            "m+3acn+2h": bool(config["Positive ion mode"]["m+3acn+2h"].split(",")[-1]),
            "m+h": bool(config["Positive ion mode"]["m+h"].split(",")[-1]),
            "m+nh4": bool(config["Positive ion mode"]["m+nh4"].split(",")[-1]),
            "m+na": bool(config["Positive ion mode"]["m+na"].split(",")[-1]),
            "m+ch3oh+h": bool(config["Positive ion mode"]["m+ch3oh+h"].split(",")[-1]),
            "m+k": bool(config["Positive ion mode"]["m+k"].split(",")[-1]),
            "m+acn+h": bool(config["Positive ion mode"]["m+acn+h"].split(",")[-1]),
            "m+2na-h": bool(config["Positive ion mode"]["m+2na-h"].split(",")[-1]),
            "m+isoprop+h": bool(config["Positive ion mode"]["m+isoprop+h"].split(",")[-1]),
            "m+acn+na": bool(config["Positive ion mode"]["m+acn+na"].split(",")[-1]),
            "m+2k+h": bool(config["Positive ion mode"]["m+2k+h"].split(",")[-1]),
            "m+dmso+h": bool(config["Positive ion mode"]["m+dmso+h"].split(",")[-1]),
            "m+2acn+h": bool(config["Positive ion mode"]["m+2acn+h"].split(",")[-1]),
            "m+isoprop+na+h": bool(config["Positive ion mode"]["m+isoprop+na+h"].split(",")[-1]),
            "2m+h": bool(config["Positive ion mode"]["2m+h"].split(",")[-1]),
            "2m+nh4": bool(config["Positive ion mode"]["2m+nh4"].split(",")[-1]),
            "2m+na": bool(config["Positive ion mode"]["2m+na"].split(",")[-1]),
            "2m+3h2o+2h": bool(config["Positive ion mode"]["2m+3h2o+2h"].split(",")[-1]),
            "2m+k": bool(config["Positive ion mode"]["2m+k"].split(",")[-1]),
            "2m+acn+h": bool(config["Positive ion mode"]["2m+acn+h"].split(",")[-1]),
            "2m+acn+na": bool(config["Positive ion mode"]["2m+acn+na"].split(",")[-1])
        },
            "negative": {
                "m-3h": bool(config["Negative ion mode"]["m-3h"].split(",")[-1]),
                "m-2h": bool(config["Negative ion mode"]["m-2h"].split(",")[-1]),
                "m-h2o-h": bool(config["Negative ion mode"]["m-h2o-h"].split(",")[-1]),
                "m-h": bool(config["Negative ion mode"]["m-h"].split(",")[-1]),
                "m+na-2h": bool(config["Negative ion mode"]["m+na-2h"].split(",")[-1]),
                "m+cl": bool(config["Negative ion mode"]["m+cl"].split(",")[-1]),
                "m+k-2h": bool(config["Negative ion mode"]["m+k-2h"].split(",")[-1]),
                "m+fa-h": bool(config["Negative ion mode"]["m+fa-h"].split(",")[-1]),
                "m+hac-h": bool(config["Negative ion mode"]["m+hac-h"].split(",")[-1]),
                "m+br": bool(config["Negative ion mode"]["m+br"].split(",")[-1]),
                "m+tfa-h": bool(config["Negative ion mode"]["m+tfa-h"].split(",")[-1]),
                "2m-h": bool(config["Negative ion mode"]["2m-h"].split(",")[-1]),
                "2m+fa-h": bool(config["Negative ion mode"]["2m+fa-h"].split(",")[-1]),
                "2m+hac-h": bool(config["Negative ion mode"]["2m+hac-h"].split(",")[-1]),
                "3m-h": bool(config["Negative ion mode"]["3m-h"].split(",")[-1]),
            }
        }
    }
    ms_mode_selector = {"Positive": "ms_pos", "Negative": "ms_neg"}

    #   WINDOW 1 - PLATE LAYOUT #
    graph_plate = window["-CANVAS-"]
    dragging = False
    temp_selector = False
    plate_active = False

    #   WINDOW 2 - BIO EXP  #
    graph_bio_exp = window["-BIO_INFO_CANVAS-"]

    #   WINDOW 2 - PURITY INFO  #
    temp_purity_info_canvas = None
    purit_info_values = False
    purity_info_samples = None
    lc_graph_showing = [keys for keys in list(config["lc_mapping"].keys())]
    color_select = {}
    for keys in list(config["plate_colouring"].keys()):
        color_select[keys] = config["plate_colouring"][keys]

    start_point = end_point = prior_rect = temp_tool = None
    well_dict = {}

    bio_info_sub_setting_tab_mapping_calc = False
    bio_info_sub_setting_tab_matrix_calc = False
    bio_info_sub_setting_tab_list_calc = False
    bio_info_sub_setting_tab_z_prime_calc = False
    bio_info_sub_setting_tab_plate_overview_calc = False
    bio_info_sub_setting_tab_overview_calc = False
    bio_info_sub_setting_tab_hit_list_calc = False

    # COMPOUND TABLE CONSTANTS #
    all_data = None
    treedata = None

    # BIO EXP TABLE CONSTANTS:
    bio_exp_table_data = None
    temp_well_id_bio_info = None
    plate_bio_info = None
    well_dict_bio_info = {}
    well_dict_bio_info_calc_handler = {}

    # Table stuff
    gsc = GUISettingsController(config, bio_final_report_setup, bio_plate_report_setup, ms_settings, simple_settings)
    window.Element("-BIO_INFO_MATRIX_TABLE-").Widget.configure(displaycolumns=[])
    plate_table_headings = ["Barcode", "Compound", "Well", "Volume", "Date", "Source Barcode", "Source Well"]
    window.Element("-PLATE_TABLE_TABLE-").Widget.configure(displaycolumns=plate_table_headings)

    # all_table_data = {"-COMPOUND_INFO_PLATE_TABLE-": "",
    #                   "-BIO_INFO_OVERVIEW_TABLE-": "",
    #                   "-BIO_INFO_OVERVIEW_AVG_TABLE-": "",
    #                   "-BIO_INFO_OVERVIEW_STDEV_TABLE-": "",
    #                   "-BIO_INFO_OVERVIEW_Z_PRIME_TABLE-": "",
    #                   "-BIO_INFO_Z_PRIME_LIST_TABLE-": "",
    #                   "-BIO_INFO_HIT_LIST_LOW_TABLE-": "",
    #                   "-BIO_INFO_HIT_LIST_MID_TABLE-": "",
    #                   "-BIO_INFO_HIT_LIST_HIGH_TABLE-": "",
    #                   "-BIO_INFO_MATRIX_TABLE-": "",
    #                   "-PURITY_INFO_OVERVIEW_TABLE-": "",
    #                   "-PURITY_INFO_PURITY_OVERVIEW_TABLE-": "",
    #                   "-PURITY_INFO_PEAK_TABLE-": "",
    #                   "-PURITY_INFO_PURITY_PEAK_LIST_TABLE-": "",
    #                   "-BIO_EXP_TABLE-": "",
    #                   "-LC_MS_SAMPLE_TABLE-": "",
    #                   "-PLATE_TABLE_TABLE-": "",
    #                   }
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
                      "-BIO_EXP_TABLE-": None,
                      "-LC_MS_SAMPLE_TABLE-": None,
                      "-PLATE_TABLE_TABLE-": None
                      }


    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == "Exit":
            break


        #   WINDOW 1 - Purity           ###
        if event == "-PURITY_DATA_IMPORT-":

            # #clearing purity info window:
            # window["-PURITY_INFO_RT_START-"].update(value="")
            # window["-PURITY_INFO_RT_END-"].update(value="")
            # window["-PURITY_INFO_WAVELENGTH-"].update(value="")
            # window["-PURITY_INFO_BIN-"].update(value="")
            # window["-PURITY_INFO_MZ-"].update(value="")
            # window["-PURITY_INFO_BATCH_BOX-"].update(values=[])
            # window["-PURITY_INFO_SAMPLE_BOX-"].update(values=[])
            # window["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"].update(values="")
            # window["-PURITY_INFO_PEAK_TABLE-"].update(values="")
            # window["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"].update(values="")
            # all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"] = None
            # all_table_data["-PURITY_INFO_PEAK_TABLE-"] = None
            # all_table_data["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"] = None
            window["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"].update(select_rows=[])
            window["-PURITY_INFO_PEAK_TABLE-"].update(select_rows=[])
            # window["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"].update(select_rows=[])
            #
            # window["-PURITY_DATA_IMPORT-"].update(text="Import Data")
            # purit_info_values = False




        #   WINDOW 2 - DRAWING THINGY!!!    ###
        print(event)



if __name__ == "__main__":
    config = configparser.ConfigParser()
    config.read("config.ini")
    main(config)

    # sg.main_get_debug_data()
