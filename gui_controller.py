import copy
import os.path
import configparser

from rdkit.Chem import MolToSmiles
import multiprocessing as mp
from multiprocessing import Queue

from excel_handler import purity_sample_layout_import, purity_sample_layout_export
from gui_layout import GUILayout
from gui_settings_control import GUISettingsController
from gui_functions import *
from gui_guards import *
from bio_data_functions import org, norm, pora, pora_internal
from plate_formatting import plate_layout_re_formate
from gui_popup import matrix_popup, sample_to_compound_name_controller, ms_raw_name_guard, bio_data_approval_table, \
    assay_generator, assay_run_naming
from gui_help_info_controller import help_info_controller
from config_writer import ConfigWriter
from database_startup import DatabaseSetUp


def main(config, queue_gui, queue_mol):
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

    # The archive could be removed from here, but then there would be a call to the DB everytime a layout is used...
    # Grabs the updated data from the database
    if db_active:
        plate_list, archive_plates_dict = get_plate_layout(config)
    else:
        plate_list = []

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
    if db_active:
        assay_table_data = grab_table_data(config, "assay")
        assay_list = []
        for row in assay_table_data[0]:
            assay_list.append(row[1])
    else:
        assay_list = []
    window["-BIO_ASSAY_NAME-"].update(values=assay_list)
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
        "pora_threshold": {"th_1": {"min": config["Settings_bio"].getfloat("final_report_pora_threshold_th_1_min"),
                                    "max": config["Settings_bio"].getfloat("final_report_pora_threshold_th_1_max"),
                                    "use": config["Settings_bio"].getboolean("final_report_pora_threshold_th_1_use")},
                           "th_2": {"min": config["Settings_bio"].getfloat("final_report_pora_threshold_th_2_min"),
                                    "max": config["Settings_bio"].getfloat("final_report_pora_threshold_th_2_max"),
                                    "use": config["Settings_bio"].getboolean("final_report_pora_threshold_th_2_use")},
                           "th_3": {"min": config["Settings_bio"].getfloat("final_report_pora_threshold_th_3_min"),
                                    "max": config["Settings_bio"].getfloat("final_report_pora_threshold_th_3_max"),
                                    "use": config["Settings_bio"].getboolean("final_report_pora_threshold_th_3_use")},
                           "th_4": {"min": config["Settings_bio"].getfloat("final_report_pora_threshold_th_4_min"),
                                    "max": config["Settings_bio"].getfloat("final_report_pora_threshold_th_4_max"),
                                    "use": config["Settings_bio"].getboolean("final_report_pora_threshold_th_4_use")},
                           "th_5": {"min": config["Settings_bio"].getfloat("final_report_pora_threshold_th_5_min"),
                                    "max": config["Settings_bio"].getfloat("final_report_pora_threshold_th_5_max"),
                                    "use": config["Settings_bio"].getboolean("final_report_pora_threshold_th_5_use")},
                           "th_6": {"min": config["Settings_bio"].getfloat("final_report_pora_threshold_th_6_min"),
                                    "max": config["Settings_bio"].getfloat("final_report_pora_threshold_th_6_max"),
                                    "use": config["Settings_bio"].getboolean("final_report_pora_threshold_th_6_use")},
                           "th_7": {"min": config["Settings_bio"].getfloat("final_report_pora_threshold_th_7_min"),
                                    "max": config["Settings_bio"].getfloat("final_report_pora_threshold_th_7_max"),
                                    "use": config["Settings_bio"].getboolean("final_report_pora_threshold_th_7_use")},
                           "th_8": {"min": config["Settings_bio"].getfloat("final_report_pora_threshold_th_8_min"),
                                    "max": config["Settings_bio"].getfloat("final_report_pora_threshold_th_8_max"),
                                    "use": config["Settings_bio"].getboolean("final_report_pora_threshold_th_8_use")},
                           "th_9": {"min": config["Settings_bio"].getfloat("final_report_pora_threshold_th_9_min"),
                                    "max": config["Settings_bio"].getfloat("final_report_pora_threshold_th_9_max"),
                                    "use": config["Settings_bio"].getboolean("final_report_pora_threshold_th_9_use")},
                           "th_10": {"min": config["Settings_bio"].getfloat("final_report_pora_threshold_th_10_min"),
                                     "max": config["Settings_bio"].getfloat("final_report_pora_threshold_th_10_max"),
                                     "use": config["Settings_bio"].getboolean(
                                         "final_report_pora_threshold_th_10_use")}},
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
                                   "pstdev": config["Settings_bio"].getboolean(
                                       "plate_report_calc_dict_original_pstdev"),
                                   "pvariance": config["Settings_bio"].getboolean(
                                       "plate_report_calc_dict_original_pvariance"),
                                   "variance": config["Settings_bio"].getboolean(
                                       "plate_report_calc_dict_original_variance"),
                                   "st_dev_%": config["Settings_bio"].getboolean(
                                       "plate_report_calc_dict_original_st_dev_%"),
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
                                     "pstdev": config["Settings_bio"].
                                         getboolean("plate_report_calc_dict_normalised_pstdev"),
                                     "pvariance": config["Settings_bio"].
                                         getboolean("plate_report_calc_dict_normalised_pvariance"),
                                     "variance": config["Settings_bio"].
                                         getboolean("plate_report_calc_dict_normalised_variance"),
                                     "st_dev_%": config["Settings_bio"].
                                         getboolean("plate_report_calc_dict_normalised_st_dev_%"),
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
                               "pstdev": config["Settings_bio"].getboolean("plate_report_calc_dict_pora_pstdev"),
                               "pvariance": config["Settings_bio"].getboolean("plate_report_calc_dict_pora_pvariance"),
                               "variance": config["Settings_bio"].getboolean("plate_report_calc_dict_pora_variance"),
                               "st_dev_%": config["Settings_bio"].getboolean("plate_report_calc_dict_pora_st_dev_%"),
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
                         "pstdev": config["Settings_bio"].getboolean("plate_report_plate_calc_dict_original_pstdev"),
                         "pvariance": config["Settings_bio"].getboolean(
                             "plate_report_plate_calc_dict_original_pvariance"),
                         "variance": config["Settings_bio"].getboolean(
                             "plate_report_plate_calc_dict_original_variance"),
                         "st_dev_%": config["Settings_bio"].getboolean(
                             "plate_report_plate_calc_dict_original_st_dev_%"),
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
                           "pstdev": config["Settings_bio"].getboolean(
                               "plate_report_plate_calc_dict_normalised_pstdev"),
                           "pvariance": config["Settings_bio"].getboolean(
                               "plate_report_plate_calc_dict_normalised_pvariance"),
                           "variance": config["Settings_bio"].getboolean(
                               "plate_report_plate_calc_dict_normalised_variance"),
                           "st_dev_%": config["Settings_bio"].getboolean(
                               "plate_report_plate_calc_dict_normalised_st_dev_%"),
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
                     "pstdev": config["Settings_bio"].getboolean("plate_report_plate_calc_dict_pora_pstdev"),
                     "pvariance": config["Settings_bio"].getboolean("plate_report_plate_calc_dict_pora_pvariance"),
                     "variance": config["Settings_bio"].getboolean("plate_report_plate_calc_dict_pora_variance"),
                     "st_dev_%": config["Settings_bio"].getboolean("plate_report_plate_calc_dict_pora_st_dev_%"),
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
                              "pstdev": config["Settings_bio"].
                                  getboolean("plate_report_plate_calc_dict_pora_internal_pstdev"),
                              "pvariance": config["Settings_bio"].
                                  getboolean("plate_report_plate_calc_dict_pora_internal_pvariance"),
                              "variance": config["Settings_bio"].
                                  getboolean("plate_report_plate_calc_dict_pora_internal_variance"),
                              "st_dev_%": config["Settings_bio"].
                                  getboolean("plate_report_plate_calc_dict_pora_internal_st_dev_%"),
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
        "pora_threshold": {"th_1": {"min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_1_min"),
                                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_1_max"),
                                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_1_use")},
                           "th_2": {"min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_2_min"),
                                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_2_max"),
                                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_2_use")},
                           "th_3": {"min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_3_min"),
                                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_3_max"),
                                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_3_use")},
                           "th_4": {"min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_4_min"),
                                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_4_max"),
                                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_4_use")},
                           "th_5": {"min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_5_min"),
                                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_5_max"),
                                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_5_use")},
                           "th_6": {"min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_6_min"),
                                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_6_max"),
                                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_6_use")},
                           "th_7": {"min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_7_min"),
                                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_7_max"),
                                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_7_use")},
                           "th_8": {"min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_8_min"),
                                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_8_max"),
                                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_8_use")},
                           "th_9": {"min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_9_min"),
                                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_9_max"),
                                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_9_use")},
                           "th_10": {"min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_10_min"),
                                     "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_10_max"),
                                     "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_10_use")},
                           "colour": {"th_1": config["Settings_bio"]["plate_report_pora_threshold_colour_th_1"],
                                      "th_2": config["Settings_bio"]["plate_report_pora_threshold_colour_th_2"],
                                      "th_3": config["Settings_bio"]["plate_report_pora_threshold_colour_th_3"],
                                      "th_4": config["Settings_bio"]["plate_report_pora_threshold_colour_th_4"],
                                      "th_5": config["Settings_bio"]["plate_report_pora_threshold_colour_th_5"],
                                      "th_6": config["Settings_bio"]["plate_report_pora_threshold_colour_th_6"],
                                      "th_7": config["Settings_bio"]["plate_report_pora_threshold_colour_th_7"],
                                      "th_8": config["Settings_bio"]["plate_report_pora_threshold_colour_th_8"],
                                      "th_9": config["Settings_bio"]["plate_report_pora_threshold_colour_th_9"],
                                      "th_10": config["Settings_bio"]["plate_report_pora_threshold_colour_th_10"]}
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
    temp_wave_test = True

    #   WINDOW 1 - PLATE LAYOUT #
    graph_plate = window["-RECT_BIO_CANVAS-"]
    dragging = False
    temp_selector = False
    plate_active = False

    #   WINDOW 1 - EXTRA    #
    customers_data = None
    vendors_data = None
    ac_data = None
    plate_types_data = None
    responsible_data = None
    location_data = None

    #   WINDOW 2 - BIO EXP  #
    graph_bio_exp = window["-BIO_INFO_CANVAS-"]

    #   WINDOW 2 - PURITY INFO  #
    purity_data = None
    purity_data_added_to_db = False
    update_purity_info_peak_table = True
    canvas_lines = {
        "uv": None,
        "peak_lines": {}
    }
    sample_data_file = None
    plot_style = None
    toolbar = None
    purity_info_rt_start = None
    purity_info_rt_end = None
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
    compound_info_tables = {"-COMPOUND_INFO_INFO_MP-": "-COMPOUND_INFO_INFO_MP_TABLE-",
                            "-COMPOUND_INFO_INFO_DP-": "-COMPOUND_INFO_INFO_DP_TABLE-",
                            "-COMPOUND_INFO_INFO_ASSAY-": "-COMPOUND_INFO_INFO_ASSAY_TABLE-",
                            "-COMPOUND_INFO_INFO_HITS-": "-COMPOUND_INFO_INFO_HITS_TABLE-",
                            "-COMPOUND_INFO_INFO_TRANSFERS-": "-COMPOUND_INFO_INFO_TRANSFERS_TABLE-",
                            "-COMPOUND_INFO_INFO_PURITY-": "-COMPOUND_INFO_INFO_PURITY_USED_TABLE-"}

    # BIO EXP TABLE CONSTANTS:
    all_assays = None
    bio_exp_table_data = None
    temp_well_id_bio_info = None
    plate_bio_info = None
    well_dict_bio_info = {}
    well_dict_bio_info_calc_handler = {}

    # PLATE TABLE CONSTANTS #
    plate_table_table_heading_mp = ["Barcode", "Compound", "Well", "Volume", "Date", "Active",
                                    "Freeze/Thaw", "Plate Type", "location"]
    plate_table_table_heading_dp = ["Barcode", "Compound", "Well", "Volume", "Date", "Active", "Freeze/Thaw",
                                    "Plate Type", "location", "Source Plate", "Source Well"]

    # Table stuff
    gsc = GUISettingsController(config, bio_final_report_setup, bio_plate_report_setup, ms_settings, simple_settings)
    window.Element("-BIO_INFO_MATRIX_TABLE-").Widget.configure(displaycolumns=[])
    window.Element("-PLATE_TABLE_TABLE-").Widget.configure(displaycolumns=plate_table_table_heading_mp)
    search_reverse = {}

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

        # if event == "-SUB_SEARCH_DRAW_MOL-":
        #     smiles = values["-SUB_SEARCH_SMILES-"]
        #     queue_mol.put(("start_draw_tool", smiles))


        if event == "-START_UP_DB-":
            if db_active:
                dsu = DatabaseSetUp(config, config["Database"]["database"])
            else:
                db = sg.PopupGetText("Choose database name")
                if db is not None:
                    db += ".db"
                    dsu = DatabaseSetUp(config, db)
            try:
                dsu.controller()
            except UnboundLocalError:
                print("No dp selected")
            else:
                config_update(config)
                window.close()
                window = gl.full_layout()

        #   WINDOW MENU         ###
        if event == "Open    Ctrl-O":
            db = sg.PopupGetFile("Choose database", file_types=(("Database", ".db"),))
            setting_dict = {"Database": {"database": db}}
            cw.run(setting_dict, "simple_settings", True)
            config_update(config)
            window.close()
            window = gl.full_layout()

        if event == "Save    Ctrl-S":
            sg.popup_error("Not working - Will never work... as all data is written to the DB... "
                           "unless I want to add a temp DB...")

        if event == "About...":
            sg.popup_error("Not working")

        if event == "Info":
            help_info_controller(config)

        #   WINDOW 1 - SEARCH     ###
        if event == "-SEARCH_AC-":
            ac = values["-SEARCH_AC-"]
            origin_values = []
            if ac:
                for acs in ac:
                    acs = acs.casefold()
                    for values in config[f"database_specific_{acs}"]:
                        origin_values.append(config[f"database_specific_{acs}"][values])
            window["-SEARCH_ORIGIN-"].update(values=origin_values)

        if event == "-SUB_SEARCH_METHOD-":
            if values["-SUB_SEARCH_METHOD-"] == "morgan":
                window["-SUB_SEARCH_MORGAN_OPTIONS-"].update(visible=True)
                window["-SUB_SEARCH_MORGAN_CHIRALITY-"].update(visible=True)
                window["-SUB_SEARCH_MORGAN_FEATURES-"].update(visible=True)
                window["-SUB_SEARCH_BITS_TEXT-"].update(visible=True)
                window["-SUB_SEARCH_MORGAN_BITS-"].update(visible=True)
                window["-SUB_SEARCH_BOUND_TEXT-"].update(visible=True)
                window["-SUB_SEARCH_MORGAN_RANGE-"].update(visible=True)
            else:
                window["-SUB_SEARCH_MORGAN_OPTIONS-"].update(visible=False)
                window["-SUB_SEARCH_MORGAN_CHIRALITY-"].update(visible=False)
                window["-SUB_SEARCH_MORGAN_FEATURES-"].update(visible=False)
                window["-SUB_SEARCH_BITS_TEXT-"].update(visible=False)
                window["-SUB_SEARCH_MORGAN_BITS-"].update(visible=False)
                window["-SUB_SEARCH_BOUND_TEXT-"].update(visible=False)
                window["-SUB_SEARCH_MORGAN_RANGE-"].update(visible=False)

        if event == "-SEARCH_PLATE_PRODUCTION-" and values["-SEARCH_PLATE_PRODUCTION-"] == "Daughter Plates":
            window["-SEARCH_PLATE_LAYOUT-"].update(disabled=False)
            window["-SEARCH_MP_MINIMIZED-"].update(disabled=False)
            window["-SEARCH_IGNORE_PLATED_COMPOUNDS-"].update(disabled=True)

        if event == "-SEARCH_PLATE_PRODUCTION-" and values["-SEARCH_PLATE_PRODUCTION-"] == "Mother Plates":
            window["-SEARCH_PLATE_LAYOUT-"].update(disabled=True)
            window["-SEARCH_MP_MINIMIZED-"].update(disabled=True)
            window["-SEARCH_IGNORE_PLATED_COMPOUNDS-"].update(disabled=False)
            window["-SEARCH_PLATE_LAYOUT_SAMPLE_AMOUNT-"].update(value=384)
            window["-SEARCH_PLATE_LAYOUT-"].update(value="")

        if event == "-SEARCH_PLATE_LAYOUT-":
            temp_counter = []
            for counter in archive_plates_dict[values["-SEARCH_PLATE_LAYOUT-"]]["well_layout"]:
                if archive_plates_dict[values["-SEARCH_PLATE_LAYOUT-"]]["well_layout"][counter]["state"] == "sample":
                    temp_counter.append(archive_plates_dict[values["-SEARCH_PLATE_LAYOUT-"]]["well_layout"][counter]
                                        ["well_id"])
            window["-SEARCH_PLATE_LAYOUT_SAMPLE_AMOUNT-"].update(value=len(temp_counter))

        if event == "-SUB_SEARCH_DRAW_MOL-":
            smiles = values["-SUB_SEARCH_SMILES-"]
            queue_mol.put(("start_draw_tool", smiles))

        #     WINDOW 1 - BIO DATA         ###
        if event == "-BIO_COMBINED_REPORT-" and values["-BIO_COMBINED_REPORT-"] is True:
            window["-BIO_FINAL_REPORT_INCLUDE_HITS-"].update(disabled=False)
            window["-BIO_FINAL_REPORT_INCLUDE_SMILES-"].update(disabled=False)

        if event == "-BIO_COMBINED_REPORT-" and values["-BIO_COMBINED_REPORT-"] is False:
            window["-BIO_FINAL_REPORT_INCLUDE_HITS-"].update(disabled=True)
            window["-BIO_FINAL_REPORT_INCLUDE_SMILES-"].update(disabled=True)
            window["-BIO_FINAL_REPORT_USE_THRESHOLD-"].update(disabled=True)
            window["-BIO_FINAL_REPORT_USE_AMOUNT-"].update(disabled=True)
            window["-BIO_FINAL_REPORT_THRESHOLD-"].update(value="")
            window["-BIO_FINAL_REPORT_THRESHOLD-"].update(disabled=True)
            window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(value="")
            window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(disabled=True)

        if event == "-BIO_FINAL_REPORT_INCLUDE_HITS-" and values["-BIO_FINAL_REPORT_INCLUDE_HITS-"] is True:
            window["-BIO_FINAL_REPORT_USE_THRESHOLD-"].update(disabled=False)
            window["-BIO_FINAL_REPORT_USE_THRESHOLD-"].update(value=False)
            window["-BIO_FINAL_REPORT_USE_AMOUNT-"].update(disabled=False)
            window["-BIO_FINAL_REPORT_USE_AMOUNT-"].update(value=False)

        if event == "-BIO_FINAL_REPORT_INCLUDE_HITS-" and values["-BIO_FINAL_REPORT_INCLUDE_HITS-"] is False:
            window["-BIO_FINAL_REPORT_USE_THRESHOLD-"].update(disabled=True)
            window["-BIO_FINAL_REPORT_USE_AMOUNT-"].update(disabled=True)
            window["-BIO_FINAL_REPORT_THRESHOLD-"].update(value="")
            window["-BIO_FINAL_REPORT_THRESHOLD-"].update(disabled=True)
            window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(value="")
            window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(disabled=True)

        if event == "-BIO_FINAL_REPORT_INCLUDE_SMILES-" and values["-BIO_FINAL_REPORT_INCLUDE_SMILES-"] is True:
            window["-BIO_FINAL_REPORT_INCLUDE_STRUCTURE-"].update(disabled=True)

        if event == "-BIO_FINAL_REPORT_INCLUDE_SMILES-" and values["-BIO_FINAL_REPORT_INCLUDE_SMILES-"] is False:
            window["-BIO_FINAL_REPORT_INCLUDE_STRUCTURE-"].update(disabled=False)

        if event == "-BIO_FINAL_REPORT_USE_THRESHOLD-" and values["-BIO_FINAL_REPORT_USE_THRESHOLD-"] is True:
            window["-BIO_FINAL_REPORT_USE_AMOUNT-"].update(value=False)
            window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(value="")
            window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(disabled=True)
            window["-BIO_FINAL_REPORT_THRESHOLD-"].update(disabled=False)

        if event == "-BIO_FINAL_REPORT_USE_THRESHOLD-" and values["-BIO_FINAL_REPORT_USE_THRESHOLD-"] is False:
            window["-BIO_FINAL_REPORT_THRESHOLD-"].update(value="")
            window["-BIO_FINAL_REPORT_THRESHOLD-"].update(disabled=True)

        if event == "-BIO_FINAL_REPORT_USE_AMOUNT-" and values["-BIO_FINAL_REPORT_USE_AMOUNT-"] is True:
            window["-BIO_FINAL_REPORT_USE_THRESHOLD-"].update(value=False)
            window["-BIO_FINAL_REPORT_THRESHOLD-"].update(value="")
            window["-BIO_FINAL_REPORT_THRESHOLD-"].update(disabled=True)
            window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(disabled=False)

        if event == "-BIO_FINAL_REPORT_USE_AMOUNT-" and values["-BIO_FINAL_REPORT_USE_AMOUNT-"] is False:
            window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(value="")
            window["-BIO_FINAL_REPORT_HIT_AMOUNT-"].update(disabled=True)

        if event == "-BIO_PLATE_LAYOUT-":
            well_dict.clear()
            well_dict = copy.deepcopy(archive_plates_dict[values["-BIO_PLATE_LAYOUT-"]]["well_layout"])
            well_dict = plate_layout_re_formate(well_dict)
            plate_size = archive_plates_dict[values["-BIO_PLATE_LAYOUT-"]]["plate_type"]
            archive = True
            gui_tab = "bio"
            sample_type = values["-BIO_SAMPLE_TYPE-"]
            draw_plate(config, graph_bio, plate_size, well_dict, gui_tab, archive, sample_layout=sample_type)

        if event == "-BIO_SAMPLE_TYPE-":
            if values["-BIO_PLATE_LAYOUT-"]:
                well_dict.clear()
                well_dict = copy.deepcopy(archive_plates_dict[values["-BIO_PLATE_LAYOUT-"]]["well_layout"])
                well_dict = plate_layout_re_formate(well_dict)
                plate_size = archive_plates_dict[values["-BIO_PLATE_LAYOUT-"]]["plate_type"]
                archive = True
                gui_tab = "bio"
                sample_type = values["-BIO_SAMPLE_TYPE-"]
                if sample_type != "Custom":
                    try:
                        draw_plate(config, graph_bio, plate_size, well_dict, gui_tab, archive, sample_layout=sample_type)
                    except KeyError:
                        print(f"There is no place layout for {sample_type}")

        if event == "-BIO_COMPOUND_DATA-" and values["-BIO_COMPOUND_DATA-"]:
            window["-BIO_EXPERIMENT_ADD_TO_DATABASE-"].update(value=True)

        if event == "-BIO_REPORT_ADD_COMPOUND_IDS-" and values["-BIO_REPORT_ADD_COMPOUND_IDS-"]:
            window["-BIO_COMPOUND_DATA-"].update(value=True)

        if event == "-BIO_COMPOUND_DATA-" and not values["-BIO_COMPOUND_DATA-"]:
            window["-BIO_REPORT_ADD_COMPOUND_IDS-"].update(value=False)

        if event == "-BIO_EXPERIMENT_ADD_TO_DATABASE-" and values["-BIO_EXPERIMENT_ADD_TO_DATABASE-"]:
            window["-BIO_COMPOUND_DATA-"].update(value=True)

        if event == "-BIO_COMBINED_REPORT-" and not values["-FINAL_BIO_NAME-"] and \
                values["-BIO_COMBINED_REPORT-"] is True:
            final_report_name = sg.popup_get_text("Final Report Name?")
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

        if event == "-BIO_REPORT_SETTINGS-" or event == "-PURITY_ADVANCED_SETTINGS-":
            reports = gsc.main_settings_controller(bio_final_report_setup, bio_plate_report_setup, ms_settings)
            if reports:
                bio_final_report_setup, bio_plate_report_setup, ms_settings, simple_settings = reports
                set_colours(window, reports)

        if event == "-BIO_ANALYSE_TYPE-":
            sg.popup("This functions does nothing ATM ")

        # Add a new assay to the database
        if event == "-BIO_NEW_ASSAY-":
            assay_check = assay_generator(config, plate_list)

            if assay_check:
                assay_table_data = grab_table_data(config, "assay")
                assay_list = []
                for row in assay_table_data[0]:
                    assay_list.append(row[1])

                window["-BIO_ASSAY_NAME-"].update(values=assay_list)

        if event == "-BIO_ASSAY_NAME-":
            temp_assay_data = grab_table_data(config, "assay", single_row=True, data_value=values["-BIO_ASSAY_NAME-"],
                                              headline="assay_name")
            window["-BIO_PLATE_LAYOUT-"].update(value=temp_assay_data[0][4])

        if event == "-BIO_FINAL_REPORT_CONCENTRATION_MULTIPLE-":
            sg.PopupError("Dose not work atm")
        if event == "-BIO_FINAL_REPORT_CONCENTRATION_MULTIPLE-":
            sg.PopupError("MISSING A WAY GET MULTIPLE CONCENTRATIONS!!! ")

        if event == "-BIO_CALCULATE-":
            bio_breaker = False
            if not values["-BIO_PLATE_LAYOUT-"]:
                sg.popup_error("Please choose a plate layout")
            elif not values["-BIO_IMPORT_FOLDER-"]:
                sg.popup_error("Please choose an import folder")
            elif values["-BIO_COMBINED_REPORT-"] and not values["-BIO_EXPORT_FOLDER-"]:
                sg.popup_error("Please choose an export folder")
            elif values["-BIO_EXPORT_TO_EXCEL-"] and not values["-BIO_EXPORT_FOLDER-"]:
                sg.popup_error("Please choose an export folder")
            elif values["-BIO_COMBINED_REPORT-"] and not values["-FINAL_BIO_NAME-"]:
                sg.popup_error("Please choose an Report name")
            elif values["-BIO_EXPERIMENT_ADD_TO_DATABASE-"] and not values["-BIO_ASSAY_NAME-"]:
                sg.popup_error("Please choose an Assay name")
            elif values["-BIO_EXPERIMENT_ADD_TO_DATABASE-"] and not values["-BIO_RESPONSIBLE-"]:
                sg.popup_error("Please choose an Responsible ")
            elif values["-BIO_EXPERIMENT_ADD_TO_DATABASE-"] and not values["-BIO_FINAL_REPORT_CONCENTRATION_NUMBER-"]:
                sg.popup_error("Please write a concentration for the run")
            elif values["-BIO_FINAL_REPORT_INCLUDE_HITS-"] and not values["-BIO_FINAL_REPORT_HIT_AMOUNT-"] \
                    and values["-BIO_FINAL_REPORT_INCLUDE_HITS-"] and not values["-BIO_FINAL_REPORT_THRESHOLD-"]:
                sg.popup_error("Please either select amount of hits or the threshold for the score, if "
                               "Hits are to be included in the report")
            # Missing setting move moving files after analyse is done.
            # elif not values["-BIO_ANALYSE_TYPE-"]:
            #     sg.popup_error("Please choose an analyse type")
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

                if not bio_export_folder:
                    bio_export_folder = values["-BIO_EXPORT_FOLDER-"]

                if not same_layout:
                    # If there are difference between what layout each plate is using, or if you know some data needs
                    # to be dismissed, you can choose different plate layout for each plate.
                    plate_to_layout = plate_layout_setup(bio_import_folder, values["-BIO_PLATE_LAYOUT-"], plate_list)
                    if plate_to_layout is None:
                        bio_breaker = True
                else:
                    # If all plate uses the same plate layout
                    plate_to_layout = default_plate_layout
                if not bio_breaker:
                    # If this is checked, it will ask for worklist, that can be converted to a sample dict,
                    # that can be used for finding sample info in the database.
                    if values["-BIO_COMPOUND_DATA-"]:
                        bio_sample_list = sg.popup_get_file("Please select worklist files", multiple_files=True)
                        if bio_sample_list is not None:
                            bio_sample_list = bio_sample_list.split(";")
                            bio_sample_dict, all_destination_plates = bio_compound_info_from_worklist(config, sg,
                                                                                                      bio_sample_list)
                        else:
                            bio_sample_dict = None
                            bio_breaker = True
                            all_destination_plates = []
                    else:
                        all_destination_plates = None
                        bio_sample_dict = None
                    if not bio_breaker:
                        try:
                            float(values["-BIO_FINAL_REPORT_CONCENTRATION_NUMBER-"])
                        except ValueError:
                            temp_concentration = float(sg.popup_get_text("Please provide a concentration"
                                                                         "\n numbers only"))
                        else:
                            temp_concentration = float(values["-BIO_FINAL_REPORT_CONCENTRATION_NUMBER-"])

                        analyse_method = values["-BIO_ANALYSE_TYPE-"]     # not used atm...
                        add_compound_ids = values["-BIO_REPORT_ADD_COMPOUND_IDS-"]
                        combined_report_check = values["-BIO_COMBINED_REPORT-"]
                        import_to_database_check = values["-BIO_EXPERIMENT_ADD_TO_DATABASE-"]
                        responsible = values["-BIO_RESPONSIBLE-"]
                        assay_name = values["-BIO_ASSAY_NAME-"]

                        bio_import_report_handler(config, bio_import_folder, plate_to_layout, archive_plates_dict,
                                                  bio_plate_report_setup, analyse_method, bio_sample_dict,
                                                  bio_export_folder, add_compound_ids, export_to_excel,
                                                  all_destination_plates, combined_report_check,
                                                  import_to_database_check, bio_final_report_setup, final_report_name,
                                                  include_hits, threshold, hit_amount, include_smiles,
                                                  include_structure, assay_name, responsible, temp_concentration)

        if event == "-BIO_SEND_TO_INFO-":
            sg.popup_error("Needs to be updated to do something else")
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

        #   WINDOW 1 - Purity           ###
        if event == "-PURITY_DATA_IMPORT-":
            if purit_info_values:
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
                purity_data_added_to_db = False
                purity_data = None

                window["-PURITY_DATA_IMPORT-"].update(text="Import Data")
                purit_info_values = False

            elif not purit_info_values:  # ToDo add threading
                temp_wavelength = values["-PURITY_DATA_UV_WAVE-"]
                temp_wave_test, wavelength_data = guard_purity_data_wavelength(temp_wavelength)

                if not values["-PURITY_DATA_IMPORT_FOLDER-"]:
                    sg.PopupError("Please select a folder")
                # Check for excel file if needed
                elif values["-PURITY_DATA_CALC_PURITY-"] and not values["-PURITY_DATA_COMPOUND_DATA-"] and not values[
                    "-PURITY_DATA_USE_COMPOUNDS-"]:
                    sg.PopupError("Please select an import file")

                elif not temp_wave_test:
                    sg.PopupError(wavelength_data)

                else:

                    folder = values["-PURITY_DATA_IMPORT_FOLDER-"]
                    all_data = import_ms_data(folder)

                    # GAURD# Checking if there are UV data
                    if isinstance(all_data, str):
                        sg.Popup("Missing UV Data")
                    else:
                        _, purity_samples, purity_data, missing_samples = all_data
                        if missing_samples:
                            guard = sg.PopupYesNo(
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
                                new_names = sample_to_compound_name_controller(config, purity_data, fd,
                                                                               purity_sample_layout_import,
                                                                               purity_sample_layout_export, sort_table)

                                if new_names:
                                    for sample in new_names:
                                        if new_names[sample] == "Delete":
                                            purity_data.pop(sample)
                                        else:
                                            purity_data[new_names[sample]] = purity_data.pop(sample)
                                purity_samples = []
                                for samples in purity_data:
                                    purity_samples.append(samples)

                            if values["-PURITY_DATA_ADD_TO_DATABASE-"]:
                                _ = purity_data_to_db(config, purity_data)
                                purity_data_added_to_db = True

                            peak_information, peak_table_data, sample_peak_dict = get_peak_information(
                                purity_data, slope_threshold, uv_threshold, rt_solvent_peak, sample_data,
                                wavelength_data)

                            if values["-PURITY_DATA_CALC_PURITY-"]:
                                if not values["-PURITY_DATA_USE_COMPOUNDS-"]:
                                    sample_data_file = values["-PURITY_DATA_COMPOUND_DATA-"]
                                else:
                                    sample_data_file = "compound_data"

                                ms_mode = ms_mode_selector[values["-PURITY_DATA_MS_MODE-"]]
                                delta_mass = float(values["-PURITY_DATA_MS_DELTA-"])
                                mz_threshold = int(values["-PURITY_DATA_MS_THRESHOLD-"])
                                peak_amounts = int(values["-PURITY_DATA_MS_PEAKS-"])

                                sample_data, db_data = grab_sample_data(sample_data_file, purity_data, config)

                                if not compound_data:
                                    # GUARD Check if file names the same. #ToDO this is duplicated earlier. FIX ! ! !
                                    raw_data_samples = natsort.natsorted([keys for keys in purity_data])
                                    excel_data_samples = natsort.natsorted([keys for keys in sample_data])

                                    if raw_data_samples != excel_data_samples:
                                        # If they are not the same, a popup will show, where you can set names for the data.
                                        new_names = ms_raw_name_guard(raw_data_samples, excel_data_samples, db_data,
                                                                      config)

                                    name_changer(new_names, purity_data, sample_data, peak_information,
                                                 sample_peak_dict)

                                all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"], \
                                purity_peak_list_table_data = purity_ops(sample_data, purity_data, peak_information,
                                                                         ms_mode, delta_mass, mz_threshold,
                                                                         peak_amounts)

                                add_start_end_time(purity_peak_list_table_data, sample_peak_dict)
                                window["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"]. \
                                    update(values=all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"])

                            window["-PURITY_INFO_SAMPLE_BOX-"].update(values=purity_samples)
                            window["-PURITY_DATA_IMPORT-"].update(text="Clear Purity Info")
                            purit_info_values = True
                            # window["-PURITY_INFO_OVERVIEW_TABLE-"].update(values=all_table_data["-PURITY_INFO_OVERVIEW_TABLE-"])

        if event == "-PURITY_INFO_PURITY_OVERVIEW_IMPORT-" and purity_data:  # ToDo Check if this works
            if not "-PURITY_DATA_USE_COMPOUNDS-":
                sg.PopupError("Not compound data. Can't be added to the database. ")
            else:
                if not purity_data_added_to_db:
                    batch_dict = purity_data_to_db(config, purity_data)

                purity_data_compounds_to_db(config, all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"])

        if event == "-PURITY_DATA_REPORT-":
            sg.Popup("Not working atm")  # ToDo Make this work. create a report based on data.

        #     WINDOW 1 - PLATE LAYOUT     ###
        if event == "-PLATE_LAYOUT_COLOUR_CHOSE_TARGET-":
            if values["-PLATE_LAYOUT_COLOUR_CHOSE_TARGET-"] != "None":
                window["-PLATE_LAYOUT_COLOUR_CHOSE-"].update(button_color=values["-PLATE_LAYOUT_COLOUR_CHOSE_TARGET-"])

        if event == "-DRAW-":
            well_dict.clear()
            # sets the size of the well for when it draws the plate
            graph = graph_plate
            plate_type = values["-PLATE-"]
            archive_plates = values["-ARCHIVE-"]
            gui_tab = "plate_layout"
            sample_type = values["-RECT_SAMPLE_TYPE-"]

            if values["-ARCHIVE-"]:
                try:
                    well_dict = copy.deepcopy(archive_plates_dict[values["-ARCHIVE_PLATES-"]]["well_layout"])
                    well_dict = plate_layout_re_formate(well_dict)
                    window["-PLATE-"].update(archive_plates_dict[values["-ARCHIVE_PLATES-"]]["plate_type"])
                    plate_type = archive_plates_dict[values["-ARCHIVE_PLATES-"]]["plate_type"]

                except KeyError:
                    window["-ARCHIVE-"].update(False)
                    values["-ARCHIVE-"] = False

            well_dict, min_x, min_y, max_x, max_y = draw_plate(config, graph, plate_type, well_dict, gui_tab,
                                                               archive_plates, sample_layout=sample_type)
            plate_active = True

        if event == "-EXPORT_LAYOUT-":
            if not well_dict:
                sg.PopupError("Please create a layout to Export")
            name = sg.PopupGetText("Name the file")
            if name:
                folder = sg.PopupGetFolder("Choose save location")
                if folder:
                    plate_layout_to_excel(well_dict, name, folder)

                    sg.Popup("Done")

        if event == "-SAVE_LAYOUT-":
            if not well_dict:
                sg.PopupError("Please create a layout to save")
            elif any("paint" in stuff.values() for stuff in well_dict.values()):
                sg.PopupError("Can't save layout with paint as well states")
            else:
                temp_well_dict = {}
                temp_dict_name = sg.PopupGetText("Name plate layout")

                if temp_dict_name:
                    archive_plates_dict[temp_dict_name] = {}
                    for index, well_counter in enumerate(well_dict):
                        temp_well_dict[index + 1] = copy.deepcopy(well_dict[well_counter])

                    archive_plates_dict[temp_dict_name]["well_layout"] = temp_well_dict
                    archive_plates_dict[temp_dict_name]["plate_type"] = values["-PLATE-"]

                    # saves the layout to the Database
                    # setting up the data for importing the new plate_layout to the database
                    temp_table = "plate_layout"
                    # ToDo add plate model ??
                    temp_plate_layout_data = {
                        "plate_name": temp_dict_name,
                        "plate_type": values["-PLATE-"],
                        "plate_model": "placeholder"
                    }

                    update_database(config, temp_table, temp_plate_layout_data)

                    temp_sub_plate_layout_data = {
                        "plate_sub": temp_dict_name,
                        "plate_main": temp_dict_name,
                        "plate_layout": f"{temp_well_dict}"
                    }
                    update_database(config, "plate_layout_sub", temp_sub_plate_layout_data)

                    # Updates the plate_list and archive_plates_dict with the new plate
                    plate_list, archive_plates_dict = get_plate_layout(config)

                    # Updates the window with new values
                    window["-ARCHIVE_PLATES-"].update(values=sorted(plate_list), value=plate_list[0])
                    window["-BIO_PLATE_LAYOUT-"].update(values=sorted(plate_list), value="")
                    window["-SEARCH_PLATE_LAYOUT-"].update(values=sorted(plate_list), value="")

        if event == "-DELETE_LAYOUT-":
            if not values["-ARCHIVE_PLATES-"]:
                sg.PopupError("Please select a layout to delete")
            else:
                # Set up values for the database, and deletes the record
                table = "plate_layout"
                headline = "plate_name"
                data_value = values["-ARCHIVE_PLATES-"]
                delete_records_from_database(config, table, headline, data_value)

                # Grabs the updated data from the database
                plate_list, archive_plates_dict = get_plate_layout(config)

                # Updates the window with new values
                try:
                    window["-ARCHIVE_PLATES-"].update(values=sorted(plate_list), value=plate_list[0])
                    window["-BIO_PLATE_LAYOUT-"].update(values=sorted(plate_list), value="")
                except IndexError:
                    window["-ARCHIVE_PLATES-"].update(values=[], value="")
                    window["-BIO_PLATE_LAYOUT-"].update(values=[], value="")

        if event == "-RENAME_LAYOUT-":
            if not values["-ARCHIVE_PLATES-"]:
                sg.PopupError("Please select a layout to rename")
            else:
                temp_dict_name = sg.PopupGetText("Name plate layout")
                if temp_dict_name:
                    # Updates the database with new values
                    table = "plate_layout"
                    headline = "plate_name"
                    old_value = values["-ARCHIVE_PLATES-"]
                    new_value = temp_dict_name
                    rename_record_in_the_database(config, table, headline, old_value, new_value)

                    # Grabs the updated data from the database
                    plate_list, archive_plates_dict = get_plate_layout(config)

                    # Removes the old name and data from the plate_list and archive_plates_dict and adds the new one
                    # archive_plates_dict[temp_dict_name] = archive_plates_dict[values["-ARCHIVE_PLATES-"]]
                    # archive_plates_dict.pop(values["-ARCHIVE_PLATES-"])
                    # plate_list.remove(values["-ARCHIVE_PLATES-"])
                    # plate_list.append(temp_dict_name)

                    # Updates the window with new values
                    window["-ARCHIVE_PLATES-"].update(values=sorted(plate_list), value=plate_list[0])
                    window["-BIO_PLATE_LAYOUT-"].update(values=sorted(plate_list), value="")

        # Used both for Plate layout and Bio Info
        # prints coordinate and well under the plate layout
        try:
            if event.endswith("+MOVE") and type(event) != tuple:

                if values["-BIO_INFO_CANVAS-"][0] and values["-BIO_INFO_CANVAS-"][1]:
                    try:
                        temp_well_bio_info = graph_bio_exp.get_figures_at_location(values['-BIO_INFO_CANVAS-'])[0]
                        temp_well_id_bio_info = well_dict_bio_info[temp_well_bio_info]["well_id"]
                    except IndexError:
                        temp_well_id_bio_info = ""
                    window["-INFO_BIO_GRAPH_TARGET-"].update(value=f"{temp_well_id_bio_info}")

                    if temp_well_id_bio_info:
                        temp_plate_name = values["-BIO_INFO_PLATES-"]
                        temp_analyse_method = values["-BIO_INFO_ANALYSE_METHOD-"]

                        temp_plate_wells_bio_info = plate_bio_info[temp_plate_name]["plates"][temp_analyse_method][
                            "wells"]
                        window["-INFO_BIO_WELL_VALUE-"].update(value=temp_plate_wells_bio_info[temp_well_id_bio_info])

                if values["-RECT_BIO_CANVAS-"][0] and values["-RECT_BIO_CANVAS-"][1]:
                    try:
                        temp_well = graph_plate.get_figures_at_location(values['-RECT_BIO_CANVAS-'])[0]
                        temp_well_id = well_dict[temp_well]["well_id"]
                    except (IndexError or KeyError) as error:
                        print(f"Canvas Error: {error}")
                        temp_well_id = ""
                    window["-INFO-"].update(value=f"{values['-RECT_BIO_CANVAS-']} {temp_well_id}")

        except AttributeError:
            pass

        if event == "-RECT_BIO_CANVAS-":
            x, y = values["-RECT_BIO_CANVAS-"]
            if not dragging:
                start_point = (x, y)
                dragging = True
            else:
                end_point = (x, y)
            if prior_rect:
                graph_plate.delete_figure(prior_rect)

            # Choosing which tool to pain the plate with.
            if None not in (start_point, end_point):
                # ToDo clean up the code
                if values["-RECT_SAMPLES-"]:
                    temp_draw_tool = "sample"
                elif values["-RECT_BLANK-"]:
                    temp_draw_tool = "blank"
                elif values["-RECT_MAX-"]:
                    temp_draw_tool = "max"
                if values["-RECT_MIN-"]:
                    temp_draw_tool = "minimum"
                elif values["-RECT_NEG-"]:
                    temp_draw_tool = "negative"
                elif values["-RECT_POS-"]:
                    temp_draw_tool = "positive"
                elif values["-RECT_EMPTY-"]:
                    temp_draw_tool = "empty"
                elif values["-COLOUR-"]:
                    temp_draw_tool = "paint"
                temp_selector = True
                prior_rect = graph_plate.draw_rectangle(start_point, end_point, fill_color="",
                                                        line_color="white")

        # it does not always detect this event:
        try:
            if event.endswith("+UP"):

                if temp_selector and plate_active:

                    # if you drag and let go too fast, the values are set to None. this is to handle that bug
                    if not start_point:
                        start_point = (0, 0)
                    if not end_point:
                        end_point = (0, 0)

                    # get a list of coordination within the selected area
                    temp_x = []
                    temp_y = []

                    if start_point[0] < end_point[0]:
                        for x_cord in range(start_point[0], end_point[0]):
                            temp_x.append(x_cord)
                    if start_point[0] > end_point[0]:
                        for x_cord in range(end_point[0], start_point[0]):
                            temp_x.append(x_cord)

                    if start_point[1] < end_point[1]:
                        for y_cord in range(start_point[1], end_point[1]):
                            temp_y.append(y_cord)
                    if start_point[1] > end_point[1]:
                        for y_cord in range(end_point[1], start_point[1]):
                            temp_y.append(y_cord)

                    # This is to enable clicking on wells to mark them
                    if not temp_x:
                        temp_x = [x]
                    if not temp_y:
                        temp_y = [y]

                    # makes a set, for adding wells, to avoid duplicates
                    graphs_list = set()

                    # goes over the coordinates and if they are within the bounds of the plate
                    # if that is the case, then the figure for that location is added the set
                    for index_x, cords_x in enumerate(temp_x):
                        for index_y, cords_y in enumerate(temp_y):
                            if min_x <= temp_x[index_x] <= max_x and min_y <= temp_y[index_y] <= max_y:
                                graphs_list.add(
                                    graph_plate.get_figures_at_location((temp_x[index_x], temp_y[index_y]))[0])

                    # colours the wells in different colour, depending on if they are samples or blanks
                    for wells in graphs_list:
                        color = color_select[temp_draw_tool]
                        well_state = temp_draw_tool
                        if color == "paint":
                            color = values["-PLATE_LAYOUT_COLOUR_CHOSE_TARGET-"]
                        graph_plate.Widget.itemconfig(wells, fill=color)
                        well_dict[wells]["colour"] = color
                        well_dict[wells]["state"] = well_state

                # deletes the rectangle used for selection
                if prior_rect:
                    graph_plate.delete_figure(prior_rect)

                # reset everything
                start_point, end_point = None, None
                dragging = False
                prior_rect = None
                temp_selector = False
                temp_draw_tool = None
        except AttributeError:
            pass

        #     WINDOW 1 - UPDATE Database      ###
        if event == "-UPDATE_COMPOUND-":
            if not values["-UPDATE_FOLDER-"]:
                sg.popup_error("Please select a folder containing compound data")
            else:
                update_database(config, "compound_main", values["-UPDATE_FOLDER-"])
                config_update(config)
                window.close()
                window = gl.full_layout()

        if event == "-UPDATE_MP-":
            if not values["-UPDATE_FOLDER-"]:
                sg.popup_error("Please select a folder containing MotherPlate data")
            else:
                update_database(config, "compound_mp", values["-UPDATE_FOLDER-"], "pb_mp_output")
                sg.popup("Done")

        if event == "-UPDATE_DP-":
            if not values["-UPDATE_FOLDER-"]:
                sg.popup_error("Please select a folder containing AssayPlate data")
            else:
                update_database(config, "compound_dp", values["-UPDATE_FOLDER-"], "echo_dp_out")

        if event == "-UPDATE_BIO-":
            if not values["-UPDATE_FOLDER-"]:
                sg.popup_error("Please select a folder containing bio data")
            sg.popup_error("this is not working atm")

        if event == "-UPDATE_AUTO-":
            sg.PopupOKCancel("this is not working atm")

        #     WINDOW 1 - Worklist     ###
        if event == "-TAB_GROUP_ONE-" and values["-TAB_GROUP_ONE-"] == "Worklist":
            temp_mp_plates, _ = grab_table_data(config, "mp_plates")
            worklist_mp_plates_list = []
            for rows in temp_mp_plates:
                worklist_mp_plates_list.append(rows[0])

            # sortes the table
            worklist_mp_plates_list = natsorted(worklist_mp_plates_list)

            window["-WORKLIST_MP_LIST-"].update(values=worklist_mp_plates_list)
            # window["-WORKLIST_ASSAY_LIST-"].update(values=worklist_mp_plates_list)    # ToDO add the right data here

            temp_assay_list, _ = grab_table_data(config, "assay")

        if event == "-WORKLIST_CONTROL_LAYOUT-":
            worklist_layout = sg.PopupGetFile("Please select a worklist layout file")
            worklist_data = well_compound_list(worklist_layout)
            for well_state in worklist_data:
                well_state = well_state.casefold()
                if well_state == "positive":
                    window["-WORKLIST_USE_POSITIVE_CONTROL-"].update(value=True)
                    window["-WORKLIST_POSITIVE_CONTROL_ID-"].update(disabled=False,
                                                                    value=worklist_data[well_state]["compound"])
                if well_state == "negative":
                    window["-WORKLIST_USE_NEGATIVE_CONTROL-"].update(value=True)
                    window["-WORKLIST_NEGATIVE_CONTROL_ID-"].update(disabled=False,
                                                                    value=worklist_data[well_state]["compound"])
                if well_state == "bonus":
                    window["-WORKLIST_USE_BONUS_COMPOUND-"].update(value=True)
                    window["-WORKLIST_BONUS_COMPOUND_ID-"].update(disabled=False,
                                                                  value=worklist_data[well_state]["compound"])
            window["-WORKLIST_CONTROL_LAYOUT_TARGET-"].update(value=worklist_layout)

        if event == "-WORKLIST_USE_POSITIVE_CONTROL-":
            window["-WORKLIST_POSITIVE_CONTROL_ID-"].update(disabled=not values["-WORKLIST_USE_POSITIVE_CONTROL-"])

        if event == "-WORKLIST_USE_NEGATIVE_CONTROL-":
            window["-WORKLIST_NEGATIVE_CONTROL_ID-"].update(disabled=not values["-WORKLIST_USE_NEGATIVE_CONTROL-"])

        if event == "-WORKLIST_USE_BONUS_COMPOUND-":
            window["-WORKLIST_BONUS_COMPOUND_ID-"].update(disabled=not values["-WORKLIST_USE_BONUS_COMPOUND-"])

        if event == "-WORKLIST_GENERATE-":
            if not values["-WORKLIST_PLATE_LAYOUT-"]:
                sg.PopupError("Please select a plate layout")

            elif not values["-WORKLIST_MP_LIST-"] and not values["-WORKLIST_USE_ALL_MOTHERPLATES-"]:
                sg.PopupError("Please select witch MotherPlates to use")

            elif not values["-WORKLIST_PLATE_AMOUNT-"]:
                sg.PopupError("Please write how many plates are needed")

            elif not values["-WORKLIST_INITIAL_PLATE-"]:
                sg.PopupError("Please write what the starting number is")

            elif not values["-WORKLIST_VOLUME-"]:
                sg.PopupError("Please write how much volume is needed per plate")

            elif values["-WORKLIST_USE_POSITIVE_CONTROL-"] and not values["-WORKLIST_POSITIVE_CONTROL_ID-"]:
                sg.PopupError("Please write an ID for the Positive control.")

            elif values["-WORKLIST_USE_NEGATIVE_CONTROL-"] and not values["-WORKLIST_NEGATIVE_CONTROL_ID-"]:
                sg.PopupError("Please write an ID for the Negative control.")

            elif values["-WORKLIST_USE_POSITIVE_CONTROL-"] and not values["-WORKLIST_CONTROL_LAYOUT_TARGET-"]:
                sg.PopupError("Please select a Plate Layout for the Positive control.")

            elif values["-WORKLIST_USE_NEGATIVE_CONTROL-"] and not values["-WORKLIST_CONTROL_LAYOUT_TARGET-"]:
                sg.PopupError("Please select a Plate Layout for the Negative control.")

            elif values["-WORKLIST_USE_BONUS_COMPOUND-"] and not values["-WORKLIST_BONUS_COMPOUND_ID-"]:
                sg.PopupError("Please write an ID for the Bonus Compound.")

            elif values["-WORKLIST_USE_BONUS_COMPOUND-"] and not values["-WORKLIST_CONTROL_LAYOUT_TARGET-"]:
                sg.PopupError("Please select a Plate Layout for the Negative control.")
            # ToDo make this one work
            # elif values["-WORKLIST_USE_BONUS_COMPOUND-"] and not values["-WORKLIST_BONUS_MAX-"]\
            #     or values["-WORKLIST_USE_BONUS_COMPOUND-"] and not values["-WORKLIST_BONUS_POSITIVE-"]\
            #     or values["-WORKLIST_USE_BONUS_COMPOUND-"] and not values["-WORKLIST_BONUS_EMPTY-"]\
            #     or values["-WORKLIST_USE_BONUS_COMPOUND-"] and not values["-WORKLIST_BONUS_MIN-"]\
            #     or values["-WORKLIST_USE_BONUS_COMPOUND-"] and not values["-WORKLIST_BONUS_NEGATIVE-"]\
            #     or values["-WORKLIST_USE_BONUS_COMPOUND-"] and not values["-WORKLIST_BONUS_BLANK-"]:
            #         sg.PopupError("Please select minimum one Well states where the Bonus compounds should be added")
            else:
                if values["-WORKLIST_USE_ALL_MOTHERPLATES-"]:
                    # This is a list of all MotherPlates. It is generated when the tab for "Worklist" is clicked
                    mps = worklist_mp_plates_list
                else:
                    mps = values["-WORKLIST_MP_LIST-"]

                if not values["-WORKLIST_ASSAY_LIST-"]:
                    assays = None
                    mp_check = sg.PopupYesNo("Have any of the MotherPlates been used for for this production before?")
                    if mp_check.casefold() == "yes":
                        worklist = sg.popup_get_file("Please select worklist files", multiple_files=True)
                        if worklist:
                            worklist = worklist.split(";")
                        else:
                            worklist = "cancelled"
                    else:
                        worklist = None
                else:
                    worklist = None
                    assays = values["-WORKLIST_ASSAY_LIST-"]
                if worklist != "cancelled":
                    plate_layout = archive_plates_dict[values["-WORKLIST_PLATE_LAYOUT-"]]

                    if worklist:
                        # Get data from the worklist, to see what plate and witch wells have been used before
                        used_plate_well_dict = bio_compound_info_from_worklist(config, sg, worklist)
                    elif assays:
                        print("Find assay data")  # ToDo make this work
                        used_plate_well_dict = None
                    else:
                        used_plate_well_dict = None

                    if values["-WORKLIST_USE_POSITIVE_CONTROL-"] or values["-WORKLIST_USE_NEGATIVE_CONTROL-"] or \
                            values["-WORKLIST_USE_BONUS_COMPOUND-"]:
                        use_control = True
                    else:
                        use_control = False

                    if use_control:
                        control_layout = Path(values["-WORKLIST_CONTROL_LAYOUT_TARGET-"])

                    else:
                        control_layout = None

                    control_samples = {"positive":
                                           {"use": values["-WORKLIST_USE_POSITIVE_CONTROL-"],
                                            "sample": values["-WORKLIST_POSITIVE_CONTROL_ID-"]},
                                       "negative":
                                           {"use": values["-WORKLIST_USE_NEGATIVE_CONTROL-"],
                                            "sample": values["-WORKLIST_NEGATIVE_CONTROL_ID-"]},
                                       "max": {"use": False},
                                       "minimum": {"use": False},
                                       "blank": {"use": False},
                                       "empty": {"use": False},
                                       "sample": {"use": False}
                                       }

                    bonus_compound = {"sample_name": values["-WORKLIST_BONUS_COMPOUND_ID-"].casefold(),
                                      "max": values["-WORKLIST_BONUS_MAX-"],
                                      "minimum": values["-WORKLIST_BONUS_MIN-"],
                                      "positive": values["-WORKLIST_BONUS_POSITIVE-"],
                                      "negative": values["-WORKLIST_BONUS_NEGATIVE-"],
                                      "blank": values["-WORKLIST_BONUS_BLANK-"],
                                      "empty": values["-WORKLIST_BONUS_EMPTY-"],
                                      "sample": values["-WORKLIST_BONUS_SAMPLE-"]}

                    worklist_analyse_method = values["-WORKLIST_SAMPLE_STYLE-"]
                    sample_direction = values["-WORKLIST_DROPDOWN_SAMPLE_DIRECTION-"]

                    assay_name = values["-WORKLIST_ASSAY_NAME-"]
                    plate_amount = int(values["-WORKLIST_PLATE_AMOUNT-"])
                    initial_plate = int(values["-WORKLIST_INITIAL_PLATE-"])
                    volume = float(values["-WORKLIST_VOLUME-"])
                    file, msg = generate_worklist(config, plate_amount, mps, plate_layout, used_plate_well_dict,
                                                  assay_name, initial_plate, volume, worklist_analyse_method,
                                                  sample_direction, control_layout, control_samples,
                                                  bonus_compound)

                    if not file:
                        sg.PopupError("Something crashed up")
                    elif type(msg) == str:
                        sg.Popup(f"{msg} - File still created with fewer plates, saved here {file}")
                    else:
                        sg.Popup(f"Worklist have been created and saved here: {file}")

        #       WINDOW 1 - EXTRA            ###
        if event == "-PD_METHOD_DD-" and values["-PD_METHOD_DD-"] == "Generate":
            window["-PD_SAVE_PLATES-"].update(disabled=False)
            window["-PD_WELL_AMOUNT-"].update(disabled=False)
            window["-PD_ADD_SOURCE_WELLS-"].update(value=False, disabled=True)
            window["-DP_PLATE_LAYOUT-"].update(disabled=False)
            window["-DP_SOURCE_FILE_GENERATE-"].update(disabled=False)
            window["-PD_WELL_LAYOUT-"].update(disabled=False)

        if event == "-PD_METHOD_DD-" and values["-PD_METHOD_DD-"] == "Calculate":
            window["-PD_SAVE_PLATES-"].update(disabled=True)
            window["-PD_WELL_AMOUNT-"].update(disabled=True)
            window["-PD_ADD_SOURCE_WELLS-"].update(disabled=False)
            window["-PD_SOURCE_WELL_AMOUNT-"].update(disabled=True)
            window["-DP_PLATE_LAYOUT-"].update(disabled=True)
            window["-DP_SOURCE_FILE_GENERATE-"].update(disabled=True)
            window["-PD_WELL_LAYOUT-"].update(disabled=True)

        if event == "-PD_ADD_SOURCE_WELLS-":
            window["-PD_SOURCE_WELL_AMOUNT-"].update(disabled=not values["-PD_ADD_SOURCE_WELLS-"])
            window["-PD_SOURCE_WELL_AMOUNT_INPUT-"].update(disabled=not values["-PD_ADD_SOURCE_WELLS-"])
            window["-PD_WELL_LAYOUT-"].update(disabled=not values["-PD_ADD_SOURCE_WELLS-"])

        if event == "-PD_EXECUTE_BUTTON-":
            if not values["-PD_FILE-"]:
                sg.PopupError("Please select an Excel file")
            elif values["-PD_METHOD_DD-"] == "Generate" and not values["-DP_PLATE_LAYOUT-"]:
                sg.PopupError("Please select a Plate-Layout")
            else:
                function = values["-PD_METHOD_DD-"]
                file = values["-PD_FILE-"]
                if values["-PD_WELL_AMOUNT-"] == "Input":
                    dw_amount = int(values["-PD_WELL_AMOUNT_INPUT-"])
                else:
                    dw_amount = values["-PD_WELL_AMOUNT-"]
                save_plates = values["-PD_SAVE_PLATES-"]
                well_layout = values["-PD_WELL_LAYOUT-"]
                try:
                    plate_layout = archive_plates_dict[values["-DP_PLATE_LAYOUT-"]]
                except KeyError:
                    plate_layout = None
                add_source_wells = values["-PD_ADD_SOURCE_WELLS-"]
                if values["-PD_SOURCE_WELL_AMOUNT-"] == "Minimum":
                    source_well_amount = values["-PD_SOURCE_WELL_AMOUNT-"]
                else:
                    source_well_amount = int(values["-PD_SOURCE_WELL_AMOUNT_INPUT-"])
                pb_source_file = values["-DP_SOURCE_FILE_GENERATE-"]

                state = plate_dilution(config, function, file, dw_amount, add_source_wells, source_well_amount,
                                       save_plates, well_layout, plate_layout, pb_source_file)

                sg.Popup(state)

        if event == "-EXTRA_SUB_TABS-":
            ...
        if event == "-EXTRA_SUB_DATABASE_TABS-":
            if values["-EXTRA_SUB_DATABASE_TABS-"] == "Responsible" and not responsible_data:
                table = "responsible"
                headings = config["Extra_tab_database_headings"]["responsible"].split(",")
                if db_active:
                    responsible_data = database_to_table(config, table, headings)
                else:
                    responsible_data = []
                window["-EXTRA_DATABASE_RESPONSIBLE_TABLE-"].update(values=[responsible_data])
                print(f"responsible_data: {responsible_data}")
            elif values["-EXTRA_SUB_DATABASE_TABS-"] == "Customers" and not customers_data:
                table = "customers"
                headings = config["Extra_tab_database_headings"]["customers"].split(",")
                customers_data = database_to_table(config, table, headings)
                window["-EXTRA_DATABASE_CUSTOMERS_TABLE-"].update(values=[customers_data])
                print(f"customers_data: {customers_data}")
            elif values["-EXTRA_SUB_DATABASE_TABS-"] == "Vendors" and not vendors_data:
                table = "vendors"
                headings = config["Extra_tab_database_headings"]["vendors"].split(",")
                vendors_data = database_to_table(config, table, headings)
                window["-EXTRA_DATABASE_VENDORS_TABLE-"].update(values=[vendors_data])
                print(f"vendors_data: {vendors_data}")
            elif values["-EXTRA_SUB_DATABASE_TABS-"] == "AC" and not ac_data:
                table = "origin"
                headings = config["Extra_tab_database_headings"]["origin"].split(",")
                ac_data = database_to_table(config, table, headings)
                window["-EXTRA_DATABASE_AC_TABLE-"].update(values=[ac_data])
                print(f"ac_data/origine: {ac_data}")
            elif values["-EXTRA_SUB_DATABASE_TABS-"] == "Plate Types" and not plate_types_data:
                table = "plate_type"
                headings = config["Extra_tab_database_headings"]["plate_type"].split(",")
                plate_types_data = database_to_table(config, table, headings)
                window["-EXTRA_PLATE_TYPE_LISTBOX-"].update(values=[plate_types_data])

                vendors = database_to_table(config, "vendors", "name")
                window["-EXTRA_PLATE_TYPE_VENDOR-"].update(values=vendors)
                print(f"plate_types_data: {plate_types_data}")

            elif values["-EXTRA_SUB_DATABASE_TABS-"] == "Location" and not location_data:
                table = "locations"
                headings = config["Extra_tab_database_headings"]["locations"].split(",")
                location_data = database_to_table(config, table, headings)
                window["-EXTRA_DATABASE_LOCATIONS_TABLE-"].update(values=[location_data])
                print(f"location_data: {location_data}")

        if event == "-EXTRA_DATABASE_RESPONSIBLE_IMPORT_DB-":
            if not values["-EXTRA_DATABASE_RESPONSIBLE_NAME-"]:
                sg.popup_error("Please fill in name")
            elif not values["-EXTRA_DATABASE_RESPONSIBLE_E_MAIL-"]:
                sg.popup_error("Please fill in E-mail")
            else:
                table = "responsible"
                row_id = get_number_of_rows(config, table) + 1
                name = values["-EXTRA_DATABASE_RESPONSIBLE_NAME-"]
                e_mail = values["-EXTRA_DATABASE_RESPONSIBLE_E_MAIL-"]
                info = values["-EXTRA_DATABASE_RESPONSIBLE_INFO-"]
                data = {"row_id": row_id,
                        "name": name,
                        "e_mail": e_mail,
                        "info": info}
                update_database(config, table, data)

        if event == "-EXTRA_DATABASE_CUSTOMERS_IMPORT_DB-":
            if not values["-EXTRA_DATABASE_CUSTOMERS_NAME-"]:
                sg.popup_error("Please fill in name")
            elif not values["-EXTRA_DATABASE_CUSTOMERS_E_MAIL-"]:
                sg.popup_error("Please fill in E-mail")
            else:
                table = "customers"
                row_id = get_number_of_rows(config, table) + 1
                name = values["-EXTRA_DATABASE_CUSTOMERS_NAME-"]
                e_mail = values["-EXTRA_DATABASE_CUSTOMERS_E_MAIL-"]
                info = values["-EXTRA_DATABASE_CUSTOMERS_INFO-"]
                data = {"row_id": row_id,
                        "name": name,
                        "e_mail": e_mail,
                        "info": info}
                update_database(config, table, data)

        if event == "-EXTRA_DATABASE_VENDORS_IMPORT_DB-":
            if not values["-EXTRA_DATABASE_VENDORS_NAME-"]:
                sg.popup_error("Please fill in name")
            elif not values["-EXTRA_DATABASE_VENDORS_E_MAIL-"]:
                sg.popup_error("Please fill in E-mail")
            else:
                table = "vendors"
                row_id = get_number_of_rows(config, table) + 1
                name = values["-EXTRA_DATABASE_VENDORS_NAME-"]
                e_mail = values["-EXTRA_DATABASE_VENDORS_E_MAIL-"]
                info = values["-EXTRA_DATABASE_VENDORS_INFO-"]
                data = {"row_id": row_id,
                        "name": name,
                        "e_mail": e_mail,
                        "info": info}
                update_database(config, table, data)

        if event == "-EXTRA_DATABASE_AC_IMPORT_DB-":
            if not values["-EXTRA_DATABASE_AC_NAME-"]:
                sg.popup_error("Please fill in name")
            elif not values["-EXTRA_DATABASE_AC_AC-"]:
                sg.popup_error("Please Choose A/C")
            else:
                table = "origin"
                row_id = get_number_of_rows(config, table) + 1
                origin = values["-EXTRA_DATABASE_AC_NAME-"]
                ac = values["-EXTRA_DATABASE_AC_AC-"]
                contact_person = values["-EXTRA_DATABASE_AC_CONTACT_NAME-"]
                e_mail = values["-EXTRA_DATABASE_AC_E_MAIL-"]
                info = values["-EXTRA_DATABASE_AC_INFO-"]
                data = {"ac_id": row_id,
                        "origin": origin,
                        "ac": ac,
                        "contact_person": contact_person,
                        "e_mail": e_mail,
                        "info": info}
                update_database(config, table, data)

        if event == "-EXTRA_DATABASE_PLACE_TYPE_IMPORT_DB-":
            if not values["-EXTRA_PLATE_TYPE_NAME-"]:
                sg.popup_error("Please fill in Plate Type")
            elif not values["-EXTRA_PLATE_TYPE_VENDOR-"]:
                sg.popup_error("Please Choose a Vendor")
            elif not values["-EXTRA_PLATE_TYPE_PRODUCT_NUMBER-"]:
                sg.popup_error("Please fill in product number")
            else:
                table = "plate_type"
                try:
                    data = {"row_id": "",
                            "plate_type": values["-EXTRA_PLATE_TYPE_NAME-"],
                            "vendor": values["-EXTRA_PLATE_TYPE_VENDOR-"][0],
                            "product_number": values["-EXTRA_PLATE_TYPE_PRODUCT_NUMBER-"],
                            "sterile": values["-EXTRA_PLATE_TYPE_STERILE-"],
                            "info": values["-EXTRA_PLATE_TYPE_INFO-"],
                            "size": int(values["-EXTRA_PLATE_TYPE_SIZE-"]),
                            "well_offset_x": float(values["-EXTRA_PLATE_TYPE_WELL_OFFSET_X-"]),
                            "well_offset_y": float(values["-EXTRA_PLATE_TYPE_WELL_OFFSET_Y-"]),
                            "well_spacing_x": float(values["-EXTRA_PLATE_TYPE_WELL_SPACING_X-"]),
                            "well_spacing_y": float(values["-EXTRA_PLATE_TYPE_WELL_SPACING_Y-"]),
                            "plate_height": float(values["-EXTRA_PLATE_TYPE_PLATE_HEIGHT-"]),
                            "plate_height_lid": float(values["EXTRA_PLATE_TYPE_PLATE_HEIGHT_LID-"]),
                            "flang_height": float(values["-EXTRA_PLATE_TYPE_FLANGE_HEIGHT-"]),
                            "well_depth": (float(values["-EXTRA_PLATE_TYPE_PLATE_HEIGHT-"]) -
                                           float(values["-EXTRA_PLATE_TYPE_FLANGE_HEIGHT-"])),
                            "well_width": float(values["-EXTRA_PLATE_TYPE_WELL_WIDTH-"]),
                            "max_volume": float(values["-EXTRA_PLATE_TYPE_MAX_VOL-"]),
                            "working_volume": float(values["-EXTRA_PLATE_TYPE_WORKING_VOL-"]),
                            "dead_volume": float(values["-EXTRA_PLATE_TYPE_WELL_DEAD_VOL-"])}
                except ValueError:
                    data = {"row_id": 0,
                            "plate_type": values["-EXTRA_PLATE_TYPE_NAME-"],
                            "vendor": values["-EXTRA_PLATE_TYPE_VENDOR-"][0],
                            "product_number": values["-EXTRA_PLATE_TYPE_PRODUCT_NUMBER-"],
                            "sterile": values["-EXTRA_PLATE_TYPE_STERILE-"],
                            "info": values["-EXTRA_PLATE_TYPE_INFO-"],
                            "well_offset_x": 0,
                            "well_offset_y": 0,
                            "well_spacing_x": 0,
                            "well_spacing_y": 0,
                            "plate_height": 0,
                            "plate_height_lid": 0,
                            "flang_height": 0,
                            "well_depth": 0,
                            "well_width": 0,
                            "max_volume": 0,
                            "working_volume": 0,
                            "dead_volume": 0}
                update_database(config, table, data)

        if event == "-EXTRA_DATABASE_LOCATION_IMPORT_DB-":
            if not values["-EXTRA_DATABASE_LOCATIONS_ROOM-"]:
                sg.popup_error("Please fill in Room")
            elif not values["-EXTRA_DATABASE_LOCATIONS_BUILDING-"]:
                sg.popup_error("Please fill in Building")
            else:
                table = "locations"
                row_id = get_number_of_rows(config, table) + 1
                room = values["-EXTRA_DATABASE_LOCATIONS_ROOM-"]
                building = values["-EXTRA_DATABASE_LOCATIONS_BUILDING-"]
                spot = values["-EXTRA_DATABASE_LOCATIONS_SPOT-"]
                data = {"loc_id": row_id,
                        "room": room,
                        "building": building,
                        "spot": spot}
                update_database(config, table, data)

        #       WINDOW 1 - SIMULATE         ###
        if event == "-SIM_INPUT_EQ-":

            if values["-SIM_INPUT_EQ-"] == "comPOUND":
                window["-SIM_COMPOUND_FRAME-"].update(visible=True)
                window["-SIM_MP_FRAME-"].update(visible=False)
                window["-SIM_DP_FRAME-"].update(visible=False)

            elif values["-SIM_INPUT_EQ-"] == "MP Production":
                window["-SIM_COMPOUND_FRAME-"].update(visible=False)
                window["-SIM_MP_FRAME-"].update(visible=True)
                window["-SIM_DP_FRAME-"].update(visible=False)

            elif values["-SIM_INPUT_EQ-"] == "DP production":
                window["-SIM_COMPOUND_FRAME-"].update(visible=False)
                window["-SIM_MP_FRAME-"].update(visible=False)
                window["-SIM_DP_FRAME-"].update(visible=True)

        if event == "-SIM_RUN-":
            if values["-SIM_INPUT_EQ-"] == "comPOUND":
                if not values["-SIM_INPUT_COMPOUND_FILE-"]:
                    sg.popup_error("Missing Compound file")
                else:
                    tube_file = values["-SIM_INPUT_COMPOUND_FILE-"]
                    output_folder = values["-SIM_OUTPUT-"]

                    compound_freezer_to_2d_simulate(tube_file, output_folder)

                    sg.Popup("Done")

            elif values["-SIM_INPUT_EQ-"] == "MP Production":
                if not values["-SIM_INPUT_MP_FILE-"]:
                    sg.popup_error("Missing 2D barcode file")
                else:
                    output_folder = values["-SIM_OUTPUT-"]
                    barcodes_2d = values["-SIM_INPUT_MP_FILE-"]
                    mp_name = values["-SIM_MP_NAME-"]
                    trans_vol = values["-SIM_MP_VOL-"]

                    mp_production_2d_to_pb_simulate(output_folder, barcodes_2d, mp_name, trans_vol)

                    sg.Popup("Done")

            elif values["-SIM_INPUT_EQ-"] == "DP production":
                if not values["-SIM_INPUT_DP_FILE-"]:
                    sg.popup_error("Missing PlateButler file")
                else:
                    sg.Popup("not working atm")

            else:
                print(f"SIM_INPUT_EQ: {values['-SIM_INPUT_EQ-']}")

        #     WINDOW TABLES - COMPOUND TABLE      ###
        if event == "-TREE_DB-":
            try:
                temp_id = window.Element("-TREE_DB-").SelectedRows[0]
            except IndexError:
                pass
            # temp_info = window.Element("-TREE_DB-").TreeData.tree_dict[temp_id].values'
            tree_sample = compound_data[temp_id]["compound_id"]
            window["-COMPOUND_INFO_ID-"].update(value=compound_data[temp_id]["compound_id"])
            window["-COMPOUND_INFO_SMILES-"].update(value=compound_data[temp_id]["smiles"])
            window["-COMPOUND_INFO_MP_VOLUME-"].update(value=compound_data[temp_id]["volume"])
            window["-COMPOUND_INFO_PIC-"].update(data=compound_data[temp_id]["png"])
            window["-COMPOUND_INFO_ORIGIN_ID-"].update(value=compound_data[temp_id]["origin_id"])
            window["-COMPOUND_INFO_CONCENTRATION-"].update(value=compound_data[temp_id]["concentration"])

            search_limiter_origin = {"academic_commercial": {"value": values["-SEARCH_AC-"],
                                                             "operator": "IN",
                                                             "target_column": "ac",
                                                             "use": ac_use},
                                     "vendor_center": {"value": values["-SEARCH_ORIGIN-"],
                                                       "operator": "IN",
                                                       "target_column": "origin",
                                                       "use": origin_use}}
            all_data_origin, _ = grab_table_data(config, "origin", search_limiter_origin)
            window["-COMPOUND_INFO_AC-"].update(value=all_data_origin[0][1])
            window["-COMPOUND_INFO_ORIGIN-"].update(value=all_data_origin[0][2])

            compound_id = compound_data[temp_id]["compound_id"]

            # Table updates:
            all_table_data["-COMPOUND_INFO_ALL_PLATE_INFO_TABLE-"], \
            all_table_data["-COMPOUND_INFO_MP_PLATE_INFO_TABLE-"], \
            all_table_data["-COMPOUND_INFO_DP_PLATE_INFO_TABLE-"], \
            all_table_data["-COMPOUND_INFO_BIO_INFO_TABLE-"], \
            all_table_data["-COMPOUND_INFO_PURITY_INFO_TABLE-"] = compound_info_table_data(config, tree_sample)

            window["-COMPOUND_INFO_ALL_PLATE_INFO_TABLE-"].update(
                values=all_table_data["-COMPOUND_INFO_ALL_PLATE_INFO_TABLE-"])
            window["-COMPOUND_INFO_MP_PLATE_INFO_TABLE-"].update(
                values=all_table_data["-COMPOUND_INFO_MP_PLATE_INFO_TABLE-"])
            window["-COMPOUND_INFO_DP_PLATE_INFO_TABLE-"].update(
                values=all_table_data["-COMPOUND_INFO_DP_PLATE_INFO_TABLE-"])
            window["-COMPOUND_INFO_BIO_INFO_TABLE-"].update(
                values=all_table_data["-COMPOUND_INFO_BIO_INFO_TABLE-"])
            window["-COMPOUND_INFO_PURITY_INFO_TABLE-"].update(
                values=all_table_data["-COMPOUND_INFO_PURITY_INFO_TABLE-"])

        if event == "-C_TABLE_REFRESH-":
            if not compound_table_clear:
                if values["-SEARCH_ALL_COMPOUNDS-"]:
                    values["-IGNORE_ACTIVE-"] = True
                    values["-SUB_SEARCH-"] = False
                    values["-SEARCH_PLATE_AMOUNT_MAX-"] = True
                    values["-SEARCH_IGNORE_VOLUME-"] = True
                    values["-SEARCH_AC-"] = None
                    values["-SEARCH_ORIGIN-"] = None

                if values["-SEARCH_PLATE_PRODUCTION-"] == "Daughter Plates":
                    table = "join_main_mp"

                elif values["-SEARCH_PLATE_PRODUCTION-"] == "Mother Plates":
                    table = config["Tables"]["compound_main"]
                current_table_data = values["-SEARCH_PLATE_PRODUCTION-"]

                if values["-SEARCH_PLATE_AMOUNT-"] == "" and not values["-SEARCH_PLATE_AMOUNT_MAX-"]:
                    sg.popup_error("Please fill out plate amount")
                elif not values["-SEARCH_TRANS_VOL-"] and not values["-SEARCH_IGNORE_VOLUME-"]:
                    sg.popup_error("Please specify transferee amount")
                else:
                    if not values["-SEARCH_PLATE_AMOUNT_MAX-"]:
                        mp_amount = int(values["-SEARCH_PLATE_AMOUNT-"])
                    else:
                        mp_amount = None
                    if not values["-SEARCH_IGNORE_VOLUME-"]:
                        vol_converter = {"mL": 1, "uL": 100, "nL": 1000000}

                        transferee_volume = float(values["-SEARCH_TRANS_VOL-"]) / vol_converter[
                            values["-SEARCH_VOL_PARAMETERS-"]]
                    else:
                        transferee_volume = None

                    ignore_active = values["-SEARCH_IGNORE_PLATED_COMPOUNDS-"]
                    sub_search = values["-SUB_SEARCH-"]
                    smiles = values["-SUB_SEARCH_SMILES-"]
                    sub_search_methode = values["-SUB_SEARCH_METHOD-"]
                    threshold = float(values["-SUB_SEARCH_THRESHOLD-"])
                    source_table = table

                    if values["-SEARCH_AC-"]:
                        ac_use = True
                    else:
                        ac_use = False

                    if values["-SEARCH_ORIGIN-"]:
                        origin_use = True
                    else:
                        origin_use = False

                    samples_per_plate = int(values["-SEARCH_PLATE_LAYOUT_SAMPLE_AMOUNT-"])
                    search_limiter = {
                        config["Tables"]["compound_source"]: {"academic_commercial": {"value": values["-SEARCH_AC-"],
                                                                                      "operator": "IN",
                                                                                      "target_column": "ac",
                                                                                      "use": ac_use},
                                                              "vendor_center": {"value": values["-SEARCH_ORIGIN-"],
                                                                                "operator": "IN",
                                                                                "target_column": "origin",
                                                                                "use": origin_use}},
                        config["Tables"]["compound_main"]: {"origin_id": {"value": "",
                                                                          "operator": "IN",
                                                                          "target_column": "ac_id",
                                                                          "use": ac_use},
                                                            "volume": {"value": transferee_volume,
                                                                       "operator": "<",
                                                                       "target_column": "volume",
                                                                       "use": not values["-SEARCH_IGNORE_VOLUME-"]}},
                        "join_tables": {config["Tables"]["compound_main"]: {},
                                        config["Tables"]["compound_mp_table"]: {
                                            "compound_id": {"value": "",
                                                            "operator": "IN",
                                                            "target_column": "compound_id",
                                                            "use": False}},
                                        "shared_data": "compound_id"}
                    }
                    min_mp = values["-SEARCH_MP_MINIMIZED-"]
                    table_data = table_update_tree(mp_amount, min_mp, samples_per_plate, ignore_active, sub_search,
                                                   smiles,
                                                   sub_search_methode, threshold, source_table, search_limiter, config)
                    if table_data:
                        treedata, all_data, compound_data, counter = table_data
                        window['-TREE_DB-'].image_dict.clear()
                        window["-TREE_DB-"].update(treedata)
                        window["-C_TABLE_COUNT-"].update(f"Compounds: {counter}")
                        window["-C_TABLE_REFRESH-"].update(text="Clear Table")
                        compound_table_clear = True
                    else:
                        if values["-SEARCH_ALL_COMPOUNDS-"]:
                            sg.popup_error("All compounds are in MotherPlates")
                        else:
                            sg.popup_error("Database is empty")
                    #
                    # except ValueError:
                    #     sg.Popup("Fill in missing data")

            elif compound_table_clear:
                window['-TREE_DB-'].image_dict.clear()
                treedata = sg.TreeData()
                window['-TREE_DB-'].update(treedata)
                window["-C_TABLE_REFRESH-"].update(text="Refresh")
                window["-C_TABLE_COUNT-"].update(f"Compounds: 0")
                window['-TREE_DB-'].image_dict.clear()
                compound_table_clear = False

        if event == "-C_TABLE_EXPORT-":
            if not all_data:
                sg.popup_error("Missing table data")
            elif not values["-SEARCH_OUTPUT_FOLDER-"]:
                sg.popup_error("missing folder")
            else:
                if current_table_data == "Mother Plates":
                    output_folder = values["-SEARCH_OUTPUT_FOLDER-"]
                    all_compound_data = all_data["compound_list"]

                    compound_export(output_folder, all_compound_data)
                    file_location = f"{values['-SEARCH_OUTPUT_FOLDER-']}/comPOUND"

                elif current_table_data == "Daughter Plates":

                    if not values["-SEARCH_PLATE_AMOUNT-"]:
                        sg.PopupError("Please fill out plate Amount")
                    elif not values["-SEARCH_TRANS_VOL-"]:
                        sg.PopupError("Please fill out Transferee volume")
                    elif not values["-SEARCH_PLATE_LAYOUT-"]:
                        sg.PopupError("Please select a plate Layout for the DP production")
                    dp_name = sg.PopupGetText("Dp name? ")
                    if dp_name:
                        plate_layout = archive_plates_dict[values["-SEARCH_PLATE_LAYOUT-"]]
                        sample_amount = int(values["-SEARCH_PLATE_LAYOUT_SAMPLE_AMOUNT-"])
                        vol_converter = {"mL": 1000000, "uL": 10000, "nL": 1}

                        transferee_volume = values["-SEARCH_TRANS_VOL-"] * vol_converter[
                            values["-SEARCH_VOL_PARAMETERS-"]]
                        output_folder = values["-SEARCH_OUTPUT_FOLDER-"]
                        mp_data = all_data["mp_data"]

                        dp_creator(config, plate_layout, sample_amount, mp_data, transferee_volume, dp_name,
                                   output_folder)

                        file_location = f"{values['-SEARCH_OUTPUT_FOLDER-']}/dp_output/"

                sg.popup(f"Done - files are located {file_location}")

        #   WINDOW TABLE - BIO EXPERIMENT TABLE     ###

        # ToDo Re-design
        # if event == "-BIO_EXP_PLATE_TABLE-":
        #     if bio_exp_table_data:
        #         if not values["-BIO_INFO_ANALYSE_METHOD-"]:
        #             window["-BIO_INFO_ANALYSE_METHOD-"].Update(value="original")
        #         if not values["-BIO_INFO_MAPPING-"]:
        #             window["-BIO_INFO_MAPPING-"].Update(value="state mapping")
        #
        #     gui_tab = "bio_exp"
        #     archive = True
        #
        #     file_name = "bio_experiments.txt"
        #     plate_dict_name = bio_exp_table_data[values["-BIO_EXP_PLATE_TABLE-"][0]][2]
        #     plate_bio_info = ...
        #
        #     bio_info_plate_layout = bio_exp_table_data[values["-BIO_EXP_PLATE_TABLE-"][0]][3]
        #     bio_info_plate_size = archive_plates_dict[bio_info_plate_layout]["plate_type"]
        #     bio_info_state_dict = copy.deepcopy(archive_plates_dict[bio_info_plate_layout]["well_layout"])
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
        #     window["-BIO_INFO_ANALYSE_METHOD-"].update(values=bio_info_analyse_method, value=bio_info_analyse_method[0])
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

        # if event == "-BIO_EXP_TABLE_REFRESH-":
        #     # ToDO make this work better...
        #     start_date = values["-BIO_EXP_TABLE_DATE_START_TARGET-"]
        #     end_date = values["-BIO_EXP_TABLE_DATE_END_TARGET-"]
        #     bio_exp_table_responsible = values["-BIO_EXP_TABLE_RESPONSIBLE-"]
        #
        #     if start_date:
        #         use_start_date = True
        #     else:
        #         use_start_date = False
        #     if end_date:
        #         use_end_date = True
        #     else:
        #         use_end_date = False
        #     if bio_exp_table_responsible:
        #         use_responsible = True
        #     else:
        #         use_responsible = False
        #
        #     search_limiter = {
        #         "start_date": {"value": start_date, "operator": "<", "target_column": "date", "use": use_start_date},
        #         "end_date": {"value": end_date, "operator": ">", "target_column": "date", "use": use_end_date},
        #         # "responsible": {"value": bio_exp_table_responsible, "operator": "=", "target_column": "responsible",
        #         #                 "use": use_responsible}
        #     }
        #
        #     table_name = "assay_runs"
        #
        #     all_table_data["-BIO_EXP_PLATE_TABLE-"], headlines = grab_table_data(config, table_name, search_limiter)
        #     # print(table_data)
        #
        #     window["-BIO_EXP_PLATE_TABLE-"].update(values=all_table_data["-BIO_EXP_PLATE_TABLE-"])

        if event == "-TABLE_TAB_GRP-" and values["-TABLE_TAB_GRP-"] == "Bio Experiment table":

            if all_assays is None:
                all_assays_data, _ = grab_table_data(config, "assay")
                all_assays = []
                for rows in all_assays_data:
                    all_assays.append(rows[1])

                window["-BIO_EXP_TABLE_ASSAY_LIST_BOX-"].update(values=all_assays)

        if event == "-BIO_EXP_TABLE_ASSAY_LIST_BOX-" and values["-BIO_EXP_TABLE_ASSAY_LIST_BOX-"]:
            table_name = "assay_runs"
            selected_assays = values["-BIO_EXP_TABLE_ASSAY_LIST_BOX-"]
            selected_headlines = ["run_name", "batch", "date", "note"]
            search_list_clm = "assay_name"
            bio_exp_assay_runs, _ = grab_table_data(config, table_name, selected_assays,
                                                    specific_rows=selected_headlines, search_list_clm=search_list_clm)

            window["-BIO_EXP_ASSAY_RUN_TABLE-"].update(values=bio_exp_assay_runs)

        if event == "-BIO_EXP_ASSAY_RUN_TABLE-" and values["-BIO_EXP_ASSAY_RUN_TABLE-"]:
            table_name = "biological_plate_data"
            bio_exp_selected_runs = []
            for values in values["-BIO_EXP_ASSAY_RUN_TABLE-"]:
                bio_exp_selected_runs.append(bio_exp_assay_runs[values][0])
            search_list_clm = "assay_run"

            selected_headlines = ["plate_name", "z_prime", "approval", "note", "analysed_method", "assay_run"]

            bio_exp_plate_data, _ = grab_table_data(config, table_name, bio_exp_selected_runs,
                                                    specific_rows=selected_headlines, search_list_clm=search_list_clm)

            window["-BIO_EXP_PLATE_TABLE-"].update(values=bio_exp_plate_data)

        if event == "-BIO_EXP_PLATE_TABLE-" and values["-BIO_EXP_PLATE_TABLE-"]:
            bio_exp_selected_plates = []
            table_name = "biological_compound_data"
            for values in values["-BIO_EXP_PLATE_TABLE-"]:

                bio_exp_selected_plates.append(bio_exp_plate_data[values][0])

            search_list_clm = "assay_plate"
            selected_headlines = ["compound_id", "score", "hit", "concentration", "approved", "assay_well", "note"]

            bio_exp_compound_data, _ = grab_table_data(config, table_name, bio_exp_selected_plates,
                                                       specific_rows=selected_headlines, search_list_clm=search_list_clm)
            #
            window["-BIO_EXP_COMPOUND_TABLE-"].update(values=bio_exp_compound_data)

        #   WINDOW TABLE - LC EXPERIMENT    ###
        if event == "-TABLE_TAB_GRP-" and values["-TABLE_TAB_GRP-"] == "LC Experiment table":
            # print("update listbox with data, if list box is empty")
            lc_exp_data, headlines = grab_table_data(config, "lc_experiment")
            window["-LC_MS_TABLE_BATCH_LIST_BOX-"].update(values=lc_exp_data)

        if event == "-LC_MS_TABLE_DATE_START_TARGET-" or event == "-LC_MS_TABLE_DATE_END_TARGET-":
            start_date = values["-LC_MS_TABLE_DATE_START_TARGET-"]
            end_date = values["-LC_MS_TABLE_DATE_END_TARGET-"]

            if start_date:
                use_start_date = True
            else:
                use_start_date = False
            if end_date:
                use_end_date = True
            else:
                use_end_date = False

            search_limiter = {
                "start_date": {"value": start_date, "operator": "<", "target_column": "date", "use": use_start_date},
                "end_date": {"value": end_date, "operator": ">", "target_column": "date", "use": use_end_date},
            }

            table_name = "lc_experiment"

            table_data, _ = grab_table_data(config, table_name, search_limiter)
            table_data_2, _ = grab_table_data(config, table_name)

            window["-LC_MS_TABLE_BATCH_LIST_BOX-"].update(values=table_data)

        if event == "-LC_MS_TABLE_BATCH_LIST_BOX-":
            batch_date = values["-LC_MS_TABLE_BATCH_LIST_BOX-"]
            batch = []
            for data in batch_date:
                batch.append(data[0])

            if batch:
                search_limiter = {
                    "batch": {"value": batch, "operator": "IN", "target_column": "batch", "use": True},
                }
                table_name = "lc_raw"
                all_table_data["-LC_MS_SAMPLE_TABLE-"], _ = grab_table_data(config, table_name)

                window["-LC_MS_SAMPLE_TABLE-"].update(values=all_table_data["-LC_MS_SAMPLE_TABLE-"])

        #   WINDOW TABLE - PLATE TABLE      ###
        if event == "-TABLE_TAB_GRP-" and values["-TABLE_TAB_GRP-"] == "Plate tables":
            temp_mp_plates, _ = grab_table_data(config, "mp_plates")
            mp_plates_list = []
            for rows in temp_mp_plates:
                mp_plates_list.append(rows[0])

            # sortes the table
            mp_plates_list = natsorted(mp_plates_list)

            window["-PLATE_TABLE_BARCODE_LIST_BOX-"].update(values=mp_plates_list)

        if event == "-PLATE_TABLE_CLEAR-":
            window["-PLATE_TABLE_TABLE-"].update(values=[])
            window["-PLATE_TABLE_START_DATE_TARGET-"].update(value="")
            window["-PLATE_TABLE_END_DATE_TARGET-"].update(value="")
            window["-PLATE_TABLE_CHOOSER-"].update(value="Mother Plates")

            temp_mp_plates, _ = grab_table_data(config, "mp_plates")
            mp_plates_list = []
            for rows in temp_mp_plates:
                mp_plates_list.append(rows[0])

            # sortes the table
            mp_plates_list = natsorted(mp_plates_list)

            window["-PLATE_TABLE_BARCODE_LIST_BOX-"].update(values=mp_plates_list)
            window.Element("-PLATE_TABLE_TABLE-").Widget.configure(displaycolumns=plate_table_table_heading_mp)

        if event == "-PLATE_TABLE_CHOOSER-" \
                or event == "-PLATE_TABLE_START_DATE_TARGET-" and values["-PLATE_TABLE_CHOOSER-"] \
                or event == "-PLATE_TABLE_END_DATE_TARGET-" and values["-PLATE_TABLE_CHOOSER-"]:

            window["-PLATE_TABLE_TABLE-"].update(values=[])

            table_dict = {"Mother Plates": "mp_plates", "Daughter Plates": "dp_plates"}

            if values["-PLATE_TABLE_START_DATE_TARGET-"]:
                use_start_date = True
            else:
                use_start_date = False
            if values["-PLATE_TABLE_END_DATE_TARGET-"]:
                use_end_date = True
            else:
                use_end_date = False

            search_limiter = {
                "start_date": {"value": values["-PLATE_TABLE_START_DATE_TARGET-"], "operator": "<", "target_column":
                    "date", "use": use_start_date},
                "end_date": {"value": values["-PLATE_TABLE_END_DATE_TARGET-"], "operator": ">", "target_column":
                    "date", "use": use_end_date},
            }
            plate_data, _ = grab_table_data(config, table_dict[values["-PLATE_TABLE_CHOOSER-"]], search_limiter)
            if plate_data:
                plates = []
                for plate in plate_data:
                    plates.append(plate[0])

                window["-PLATE_TABLE_BARCODE_LIST_BOX-"].update(values=plates)
            else:
                window["-PLATE_TABLE_BARCODE_LIST_BOX-"].update(values=[[]])

            if values["-PLATE_TABLE_CHOOSER-"] == "Mother Plates":
                window.Element("-PLATE_TABLE_TABLE-").Widget.configure(displaycolumns=plate_table_table_heading_mp)
            else:
                window.Element("-PLATE_TABLE_TABLE-").Widget.configure(displaycolumns=plate_table_table_heading_dp)

        if event == "-PLATE_TABLE_BARCODE_LIST_BOX-":
            table_dict = {"Mother Plates": {"clm": "mp_barcode", "table": "compound_mp"},
                          "Daughter Plates": {"clm": "dp_barcode", "table": "compound_dp"}}
            search_limiter = {"academic_commercial": {"value": values["-PLATE_TABLE_BARCODE_LIST_BOX-"],
                                                      "operator": "IN",
                                                      "target_column": table_dict[values["-PLATE_TABLE_CHOOSER-"]][
                                                          "clm"],
                                                      "use": True}}

            all_table_data["-PLATE_TABLE_TABLE-"], _ = grab_table_data(config,
                                                                       table_dict[values["-PLATE_TABLE_CHOOSER-"]][
                                                                           "table"], search_limiter)

            # print(headings)
            if values["-PLATE_TABLE_BARCODE_LIST_BOX-"]:
                window["-PLATE_TABLE_TABLE-"].update(values=all_table_data["-PLATE_TABLE_TABLE-"])
            else:
                window["-PLATE_TABLE_TABLE-"].update(values=[])

        if event == "-PLATE_TABLE_BUTTON_LIMITER-":
            mp_limiter = values["-PLATE_TABLE_TEXT_LIMITER-"]
            mp_plates_list = []
            for rows in temp_mp_plates:
                if mp_limiter in rows[0].casefold():
                    mp_plates_list.append(rows[0])

            window["-PLATE_TABLE_BARCODE_LIST_BOX-"].update(values=mp_plates_list)

        #   WINDOW 2 - COMPOUND INFO    ###
        if event == "-COMPOUND_INFO_SEARCH_COMPOUND_ID-":
            compound_id = values["-COMPOUND_INFO_ID-"]

            if compound_id:
                sample_row = grab_table_data(config, "compound_main", single_row=True, data_value=compound_id, headline="compound_id")
                if sample_row:
                    # Get Academic/commercial information:
                    ac_id = sample_row[0][6]
                    all_data_ac = grab_table_data(config, "origin", single_row=True, data_value=ac_id, headline="ac_id")

                    # Update info frame:
                    window["-COMPOUND_INFO_AC-"].update(value=all_data_ac[0][1])
                    window["-COMPOUND_INFO_ORIGIN-"].update(value=all_data_ac[0][2])
                    window["-COMPOUND_INFO_ORIGIN_ID-"].update(value=sample_row[0][7])
                    window["-COMPOUND_INFO_CONCENTRATION-"].update(value=sample_row[0][5])
                    window["-COMPOUND_INFO_TUBE_VOLUME-"].update(value=sample_row[0][4])
                    # ToDo find info from DB if that is there, else:
                    window["-COMPOUND_INFO_TUBE_VOLUME-"].update(value="Missing info")

                    # Update Picture frame:
                    window["-COMPOUND_INFO_SMILES-"].update(value=sample_row[0][2])
                    window["-COMPOUND_INFO_PIC-"].update(data=sample_row[0][3])

                    mp_table_data = grab_table_data(config, "compound_mp", single_row=True, data_value=compound_id,
                                                    headline="compound_id")

                    dp_table_data = grab_table_data(config, "compound_dp", single_row=True, data_value=compound_id,
                                                    headline="compound_id")
                    assay_compound_table_data = grab_table_data(config, "biological_compound_data", single_row=True,
                                                       data_value=compound_id, headline="compound_id")
                    print(mp_table_data)
                    print(dp_table_data)
                    print(assay_compound_table_data)
                    assay_plate = []
                    for assays in assay_compound_table_data:
                        assay_plate.append(assays[3])
                    # ToDo Figure out what to show and how to show the data for the compound
                    assay_plate_table_data, _ = grab_table_data(config, "biological_plate_data", assay_plate,
                                                             specific_rows=None,
                                                            search_list_clm="plate_name")

                    print(assay_plate_table_data)

                    assay_runs = []
                    for plates in assay_plate_table_data:
                        assay_runs.append(plates[1])

                    assay_run_table_data, _ = grab_table_data(config, "assay_runs", assay_runs,
                                                       specific_rows=None, search_list_clm="run_name")

                    # print(assay_run_table_data)

                    assays = []
                    for runs in assay_run_table_data:
                        assays.append(runs[1])
                    assay_table_data, _ = grab_table_data(config, "assay", assays, specific_rows=None,
                                                          search_list_clm="assay_name")

                    print(assay_table_data)

        if event in compound_info_tables:
            print(compound_info_tables[event])



        #   WINDOW 2 - BIO INFO         ###
        if event == "-BIO_INFO_STATES-" and values["-BIO_INFO_STATES-"]:
            update_bio_info_values(values, window, plate_bio_info)
        if event == "-BIO_INFO_ANALYSE_METHOD-" and values["-BIO_INFO_STATES-"]:
            update_bio_info_values(values, window, plate_bio_info)
        if event == "-BIO_INFO_PLATES-" and values["-BIO_INFO_STATES-"]:
            update_bio_info_values(values, window, plate_bio_info)

        # Updating Sub setting data
        if event == "-BIO_INFO_HEATMAP_LOW_COLOUR_TARGET-":
            if values["-BIO_INFO_HEATMAP_LOW_COLOUR_TARGET-"] != "None":
                window["-BIO_INFO_HEATMAP_LOW_COLOUR_BOX-"]. \
                    update(background_color=values["-BIO_INFO_HEATMAP_LOW_COLOUR_TARGET-"])
        if event == "-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-":
            if values["-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-"] != "None":
                window["-BIO_INFO_HEATMAP_MID_COLOUR_BOX-"]. \
                    update(background_color=values["-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-"])
        if event == "-BIO_INFO_HEATMAP_HIGH_COLOUR_TARGET-":
            if values["-BIO_INFO_HEATMAP_HIGH_COLOUR_TARGET-"] != "None":
                window["-BIO_INFO_HEATMAP_HIGH_COLOUR_BOX-"]. \
                    update(background_color=values["-BIO_INFO_HEATMAP_HIGH_COLOUR_TARGET-"])

        if event == "-BIO_INFO_HIT_MAP_LOW_COLOUR_TARGET-":
            if values["-BIO_INFO_HIT_MAP_LOW_COLOUR_TARGET-"] != "None":
                window["-BIO_INFO_HIT_MAP_LOW_COLOUR_BOX-"]. \
                    update(background_color=values["-BIO_INFO_HIT_MAP_LOW_COLOUR_TARGET-"])
        if event == "-BIO_INFO_HIT_MAP_MID_COLOUR_TARGET-":
            if values["-BIO_INFO_HIT_MAP_MID_COLOUR_TARGET-"] != "None":
                window["-BIO_INFO_HIT_MAP_MID_COLOUR_BOX-"]. \
                    update(background_color=values["-BIO_INFO_HIT_MAP_MID_COLOUR_TARGET-"])
        if event == "-BIO_INFO_HIT_MAP_HIGH_COLOUR_TARGET-":
            if values["-BIO_INFO_HIT_MAP_HIGH_COLOUR_TARGET-"] != "None":
                window["-BIO_INFO_HIT_MAP_HIGH_COLOUR_BOX-"]. \
                    update(background_color=values["-BIO_INFO_HIT_MAP_HIGH_COLOUR_TARGET-"])
        if event == "-BIO_INFO_BOUNDS_BUTTON-":
            sg.PopupError(
                "This is not working")  # ToDo make this button work. Should get a small popup, to choose all the bins for the bio analysis.

        if event == "-BIO_INFO_SUB_SETTINGS_TABS-" and values["-BIO_INFO_SUB_SETTINGS_TABS-"] == "Plate Overview" \
                and bio_info_sub_setting_tab_plate_overview_calc or event == "-BIO_INFO_PLATE_OVERVIEW_METHOD_LIST-" \
                or event == "-BIO_INFO_PLATE_OVERVIEW_STATE_LIST-" or event == "-BIO_INFO_PLATE_OVERVIEW_PLATE-":
            method = values["-BIO_INFO_PLATE_OVERVIEW_METHOD_LIST-"]
            state = values["-BIO_INFO_PLATE_OVERVIEW_STATE_LIST-"]
            plate = values["-BIO_INFO_PLATE_OVERVIEW_PLATE-"]
            all_table_data["-BIO_INFO_OVERVIEW_TABLE-"] = sub_settings_plate_overview(plate_bio_info, method, plate,
                                                                                      state)
            window["-BIO_INFO_OVERVIEW_TABLE-"].update(values=all_table_data["-BIO_INFO_OVERVIEW_TABLE-"])
            bio_info_sub_setting_tab_plate_overview_calc = False

        if event == "-BIO_INFO_SUB_SETTINGS_TABS-" and values["-BIO_INFO_SUB_SETTINGS_TABS-"] == "Overview" \
                and bio_info_sub_setting_tab_overview_calc or event == "-BIO_INFO_OVERVIEW_METHOD-" \
                or event == "-BIO_INFO_OVERVIEW_STATE-":
            method = values["-BIO_INFO_OVERVIEW_METHOD-"]
            state = values["-BIO_INFO_OVERVIEW_STATE-"]
            sub_settings_overview_table_data = sub_settings_overview(plate_bio_info, method, state)
            all_table_data["-BIO_INFO_OVERVIEW_AVG_TABLE-"], all_table_data["-BIO_INFO_OVERVIEW_STDEV_TABLE-"], \
            all_table_data["-BIO_INFO_OVERVIEW_Z_PRIME_TABLE-"] = sub_settings_overview_table_data
            window["-BIO_INFO_OVERVIEW_AVG_TABLE-"].update(values=all_table_data["-BIO_INFO_OVERVIEW_AVG_TABLE-"])
            window["-BIO_INFO_OVERVIEW_STDEV_TABLE-"].update(values=all_table_data["-BIO_INFO_OVERVIEW_STDEV_TABLE-"])
            window["-BIO_INFO_OVERVIEW_Z_PRIME_TABLE-"].update(
                values=all_table_data["-BIO_INFO_OVERVIEW_Z_PRIME_TABLE-"])
            bio_info_sub_setting_tab_overview_calc = False

        # if event == "-BIO_INFO_SUB_SETTINGS_TABS-" and values["-BIO_INFO_SUB_SETTINGS_TABS-"] == "List" \
        #         and bio_info_sub_setting_tab_list_calc or event == "-BIO_INFO_LIST_METHOD-" \
        #         or event == "-BIO_INFO_LIST_STATE-" or event == "-BIO_INFO_LIST_CALC-":
        #
        #     method = values["-BIO_INFO_LIST_METHOD-"]
        #     state = values["-BIO_INFO_LIST_STATE-"]
        #     calc = values["-BIO_INFO_LIST_CALC-"]
        #     sub_setting_list_table_data = sub_settings_list(plate_bio_info, method, state, calc)
        #     window["-BIO_INFO_LIST_TABLE-"].update(values=sub_setting_list_table_data)
        #     bio_info_sub_setting_tab_list_calc = False

        if event == "-BIO_INFO_SUB_SETTINGS_TABS-" and values["-BIO_INFO_SUB_SETTINGS_TABS-"] == "Z-Prime" \
                and bio_info_sub_setting_tab_z_prime_calc:
            z_prime_data = sub_settings_z_prime(plate_bio_info)
            all_table_data["-BIO_INFO_Z_PRIME_LIST_TABLE-"], z_prime_max_barcode, z_prime_max_value, z_prime_min_barcode \
                , z_prime_min_value = z_prime_data
            window["-BIO_INFO_Z_PRIME_LIST_TABLE-"].update(values=all_table_data["-BIO_INFO_Z_PRIME_LIST_TABLE-"])
            window["-BIO_INFO_Z_PRIME_MAX_BARCODE-"].update(value=z_prime_max_barcode)
            window["-BIO_INFO_Z_PRIME_MAX_VALUE-"].update(value=z_prime_max_value)
            window["-BIO_INFO_Z_PRIME_MIN_BARCODE-"].update(value=z_prime_min_barcode)
            window["-BIO_INFO_Z_PRIME_MIN_VALUE-"].update(value=z_prime_min_value)

            bio_info_sub_setting_tab_z_prime_calc = False

        if event == "-BIO_INFO_SUB_SETTINGS_TABS-" and values["-BIO_INFO_SUB_SETTINGS_TABS-"] == "Hit List" \
                and bio_info_sub_setting_tab_hit_list_calc:
            plate = values["-BIO_INFO_HIT_LIST_PLATES-"]
            method = values["-BIO_INFO_HIT_LIST_METHOD-"]
            state = values["-BIO_INFO_HIT_LIST_STATE-"]

            pora_thresholds = {
                "low": {"min": float(values["-BIO_INFO_PORA_LOW_MIN_HIT_THRESHOLD-"]),
                        "max": float(values["-BIO_INFO_PORA_LOW_MAX_HIT_THRESHOLD-"])},
                "mid": {"min": float(values["-BIO_INFO_PORA_MID_MIN_HIT_THRESHOLD-"]),
                        "max": float(values["-BIO_INFO_PORA_MID_MAX_HIT_THRESHOLD-"])},
                "high": {"min": float(values["-BIO_INFO_PORA_HIGH_MIN_HIT_THRESHOLD-"]),
                         "max": float(values["-BIO_INFO_PORA_HIGH_MAX_HIT_THRESHOLD-"])}
            }

            sub_settings_hit_list_table_data = sub_settings_hit_list(plate_bio_info, plate, method, state,
                                                                     bio_info_state_dict, pora_thresholds)

            all_table_data["-BIO_INFO_HIT_LIST_LOW_TABLE-"], all_table_data["-BIO_INFO_HIT_LIST_MID_TABLE-"], \
            all_table_data["-BIO_INFO_HIT_LIST_HIGH_TABLE-"] = sub_settings_hit_list_table_data
            window["-BIO_INFO_HIT_LIST_LOW_TABLE-"].update(values=all_table_data["-BIO_INFO_HIT_LIST_LOW_TABLE-"])
            window["-BIO_INFO_HIT_LIST_MID_TABLE-"].update(values=all_table_data["-BIO_INFO_HIT_LIST_MID_TABLE-"])
            window["-BIO_INFO_HIT_LIST_HIGH_TABLE-"].update(values=all_table_data["-BIO_INFO_HIT_LIST_HIGH_TABLE-"])

            bio_info_sub_setting_tab_hit_list_calc = False

        if event == "-BIO_INFO_PORA_LOW_MIN_HIT_THRESHOLD-" or event == "-BIO_INFO_PORA_LOW_MAX_HIT_THRESHOLD-" or \
                event == "-BIO_INFO_PORA_MID_MIN_HIT_THRESHOLD-" or event == "-BIO_INFO_PORA_MID_MAX_HIT_THRESHOLD-" or \
                event == "-BIO_INFO_PORA_HIGH_MIN_HIT_THRESHOLD-" or event == "-BIO_INFO_PORA_HIGH_MAX_HIT_THRESHOLD-":

            if not bio_info_sub_setting_tab_hit_list_calc:
                bio_info_sub_setting_tab_hit_list_calc = True

        if event == "-BIO_INFO_MAPPING-" or event == "-BIO_INFO_ANALYSE_METHOD-" or event == "-BIO_INFO_PLATES-" \
                or event == "-BIO_INFO_RE_DRAW-" \
                or event == "-BIO_INFO_SUB_SETTINGS_TABS-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
                or event == "-BIO_INFO_HIT_MAP_LOW_COLOUR_TARGET-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
                or event == "-BIO_INFO_HIT_MAP_MID_COLOUR_TARGET-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
                or event == "-BIO_INFO_HIT_MAP_HIGH_COLOUR_TARGET-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
                or event == "-BIO_INFO_PORA_LOW_MIN_HIT_THRESHOLD-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
                and values["-BIO_INFO_PORA_LOW_MIN_HIT_THRESHOLD-"] \
                or event == "-BIO_INFO_PORA_LOW_MAX_HIT_THRESHOLD-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
                and values["-BIO_INFO_PORA_LOW_MAX_HIT_THRESHOLD-"] \
                or event == "-BIO_INFO_PORA_MID_MIN_HIT_THRESHOLD-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
                and values["-BIO_INFO_PORA_MID_MIN_HIT_THRESHOLD-"] \
                or event == "-BIO_INFO_PORA_MID_MAX_HIT_THRESHOLD-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
                and values["-BIO_INFO_PORA_MID_MAX_HIT_THRESHOLD-"] \
                or event == "-BIO_INFO_PORA_HIGH_MIN_HIT_THRESHOLD-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
                and values["-BIO_INFO_PORA_HIGH_MIN_HIT_THRESHOLD-"] \
                or event == "-BIO_INFO_PORA_HIGH_MAX_HIT_THRESHOLD-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
                and values["-BIO_INFO_PORA_HIGH_MAX_HIT_THRESHOLD-"] \
                or event == "-BIO_INFO_STATE_LIST_BOX-" and values["-BIO_INFO_MAPPING-"] == "Hit Map" \
                or event == "-BIO_INFO_STATE_LIST_BOX-" and values["-BIO_INFO_MAPPING-"] == "Heatmap" \
                or event == "-BIO_INFO_HEATMAP_LOW_COLOUR_TARGET-" and values["-BIO_INFO_MAPPING-"] == "Heatmap" \
                or event == "-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-" and values["-BIO_INFO_MAPPING-"] == "Heatmap" \
                or event == "-BIO_INFO_HEATMAP_HIGH_COLOUR_TARGET-" and values["-BIO_INFO_MAPPING-"] == "Heatmap" \
                or event == "-BIO_INFO_HEAT_PERCENTILE_LOW-" and values["-BIO_INFO_MAPPING-"] == "Heatmap" \
                and values["-BIO_INFO_HEAT_PERCENTILE_LOW-"] \
                or event == "-BIO_INFO_HEAT_PERCENTILE_MID-" and values["-BIO_INFO_MAPPING-"] == "Heatmap" \
                and values["-BIO_INFO_HEAT_PERCENTILE_MID-"] \
                or event == "-BIO_INFO_HEAT_PERCENTILE_HIGH-" and values["-BIO_INFO_MAPPING-"] == "Heatmap" \
                and values["-BIO_INFO_HEAT_PERCENTILE_HIGH-"]:

            if plate_bio_info:
                if values["-BIO_INFO_MAPPING-"] == "State Mapping":
                    gui_tab = "bio_exp"
                    archive = True

                    well_dict_bio_info, bio_info_min_x, bio_info_min_y, bio_info_max_x, bio_info_max_y \
                        = draw_plate(config, graph_bio_exp, bio_info_plate_size, bio_info_state_dict, gui_tab, archive)

                if values["-BIO_INFO_MAPPING-"] == "Heatmap":
                    mapping = {
                        "mapping": values["-BIO_INFO_MAPPING-"],
                        "colours": {"low": [values["-BIO_INFO_HEATMAP_LOW_COLOUR_TARGET-"],
                                            values["-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-"]],
                                    "high": [values["-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-"],
                                             values["-BIO_INFO_HEATMAP_HIGH_COLOUR_TARGET-"]]},
                        "states": values["-BIO_INFO_STATE_LIST_BOX-"],
                        "percentile": {"low": float(values["-BIO_INFO_HEAT_PERCENTILE_LOW-"]),
                                       "mid": float(values["-BIO_INFO_HEAT_PERCENTILE_MID-"]),
                                       "high": float(values["-BIO_INFO_HEAT_PERCENTILE_HIGH-"])}
                    }

                    gui_tab = "bio_exp"
                    plate = values["-BIO_INFO_PLATES-"]
                    analyse_method = values["-BIO_INFO_ANALYSE_METHOD-"]

                    temp_plate_bio_info = plate_bio_info[plate]["plates"][analyse_method]["wells"]
                    well_dict_bio_info, bio_info_min_x, bio_info_min_y, bio_info_max_x, bio_info_max_y \
                        = draw_plate(config, graph_bio_exp, bio_info_plate_size, temp_plate_bio_info, gui_tab,
                                     mapping=mapping, state_dict=bio_info_state_dict)

                if values["-BIO_INFO_MAPPING-"] == "Hit Map":
                    mapping = {
                        "mapping": values["-BIO_INFO_MAPPING-"],
                        "lower_bound_start": float(values["-BIO_INFO_PORA_LOW_MIN_HIT_THRESHOLD-"]),
                        "lower_bound_end": float(values["-BIO_INFO_PORA_LOW_MAX_HIT_THRESHOLD-"]),
                        "middle_bound_start": float(values["-BIO_INFO_PORA_MID_MIN_HIT_THRESHOLD-"]),
                        "middle_bound_end": float(values["-BIO_INFO_PORA_MID_MAX_HIT_THRESHOLD-"]),
                        "higher_bound_start": float(values["-BIO_INFO_PORA_HIGH_MIN_HIT_THRESHOLD-"]),
                        "higher_bound_end": float(values["-BIO_INFO_PORA_HIGH_MAX_HIT_THRESHOLD-"]),
                        "low_colour": values["-BIO_INFO_HIT_MAP_LOW_COLOUR_TARGET-"],
                        "mid_colour": values["-BIO_INFO_HIT_MAP_MID_COLOUR_TARGET-"],
                        "high_colour": values["-BIO_INFO_HIT_MAP_HIGH_COLOUR_TARGET-"],
                        "states": values["-BIO_INFO_STATE_LIST_BOX-"]
                    }

                    plate = values["-BIO_INFO_PLATES-"]
                    analyse_method = values["-BIO_INFO_ANALYSE_METHOD-"]
                    gui_tab = "bio_exp"

                    temp_plate_bio_info = plate_bio_info[plate]["plates"][analyse_method]["wells"]
                    well_dict_bio_info, bio_info_min_x, bio_info_min_y, bio_info_max_x, bio_info_max_y \
                        = draw_plate(config, graph_bio_exp, bio_info_plate_size, temp_plate_bio_info, gui_tab,
                                     mapping=mapping, state_dict=bio_info_state_dict)

        if event == "-BIO_INFO_MATRIX_POPUP-" and plate_bio_info:
            method = values["-BIO_INFO_MATRIX_METHOD-"]
            state = values["-BIO_INFO_MATRIX_STATE-"]
            calc = values["-BIO_INFO_MATRIX_CALC-"]

            matrix_popup(plate_bio_info, calc_values, state_values, method_values, calc, sub_settings_matrix, state,
                         method)

        if event == "-BIO_INFO_Z_PRIME_MATRIX_BUTTON-" and plate_bio_info:
            matrix_popup(plate_bio_info, calc_values, state_values, method_values, "z_prime", sub_settings_matrix)

        if event == "-BIO_INFO_MATRIX_BUTTON-" and plate_bio_info:
            data_dict = plate_bio_info
            state = values["-BIO_INFO_MATRIX_STATE-"]
            if state == "z_prime":
                calc = None
                method = None
            else:
                calc = values["-BIO_INFO_MATRIX_CALC-"]
                method = values["-BIO_INFO_MATRIX_METHOD-"]

            try:
                all_table_data["-BIO_INFO_MATRIX_TABLE-"], display_columns = sub_settings_matrix(data_dict, calc,
                                                                                                 method, state)
            except KeyError:
                sg.popup_error("Please select all information")
            else:
                window.Element("-BIO_INFO_MATRIX_TABLE-").Widget.configure(displaycolumns=display_columns)

                window["-BIO_INFO_MATRIX_TABLE-"].update(values=all_table_data["-BIO_INFO_MATRIX_TABLE-"])
            # window["-BIO_INFO_MATRIX_TABLE-"].update(headings=headings)

        if event == "-BIO_INFO_EXPORT-":
            sg.popup("This functions does nothing ATM ")

        #   WINDOW 2 - PURITY INFO  ###
        if event == "-PURITY_INFO_SAMPLE_SELECTION-":
            if values["-PURITY_INFO_SAMPLE_SELECTION-"]:
                window["-PURITY_INFO_SAMPLE_BOX-"].update(select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE)
            else:
                window["-PURITY_INFO_SAMPLE_BOX-"].update(select_mode=sg.LISTBOX_SELECT_MODE_SINGLE)

        if event == "-PURITY_INFO_RE_CALC-" and values["-PURITY_INFO_SAMPLE_BOX-"]:
            samples = values["-PURITY_INFO_SAMPLE_BOX-"]
            slope_threshold = int(values["-PURITY_INFO_SLOPE_THRESHOLD-"])
            uv_threshold = int(values["-PURITY_INFO_UV_THRESHOLD-"])
            rt_solvent_peak = float(values["-PURITY_INFO_RT_SOLVENT-"])
            wavelength_data = values["-PURITY_INFO_UV_WAVE-"]

            for sample in samples:
                temp_peak_information, temp_peak_table_data, temp_sample_peak_dict = get_peak_information(
                    purity_data, slope_threshold, uv_threshold, rt_solvent_peak, None, wavelength_data, sample)

                peak_information[sample] = temp_peak_information[sample]
                peak_table_data[sample] = temp_peak_table_data[sample]
                sample_peak_dict[sample] = temp_sample_peak_dict[sample]

            if values["-PURITY_INFO_USE_MS_DATA-"]:
                if not sample_data_file:
                    mass = sg.popup_get_text("What is the mass of the compound?")
                else:
                    mass = None

                if values["-PURITY_INFO_MS_MODE_POS-"]:
                    ms_mode = "ms_pos"
                elif values["-PURITY_INFO_MS_MODE_NEG-"]:
                    ms_mode = "ms_neg"

                delta_mass = float(values["-PURITY_INFO_MS_DELTA-"])
                mz_threshold = int(values["-PURITY_INFO_MS_THRESHOLD-"])
                peak_amounts = int(values["-PURITY_INFO_MS_PEAKS-"])

                sample_data, _ = grab_sample_data(sample_data_file, purity_data, config)

                all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"], \
                purity_peak_list_table_data = purity_ops(sample_data, purity_data, peak_information, ms_mode,
                                                         delta_mass, mz_threshold, peak_amounts, mass)

                add_start_end_time(purity_peak_list_table_data, sample_peak_dict)
                window["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"]. \
                    update(values=all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"])

        if event == "-PURITY_INFO_DRAW_STUFF-":
            sg.Popup("DO NOT WORK ATM :D Still figuring out how to do stuff")

        if event == "-PURITY_INFO_CANVAS-":
            print(f"PURITY_INFO_CANVAS: {values}")

        # if event == "-PURITY_INFO_PURITY_OVERVIEW_TABLE-":
        #     print(purity_peak_list_table_data)
        #   WINDOW 2 - DRAWING THINGY!!!    ###
        if event == "-PURITY_INFO_SAMPLE_BOX-" and values["-PURITY_INFO_SAMPLE_BOX-"] \
                or event == "-PURITY_INFO_PURITY_OVERVIEW_TABLE-" and purit_info_values and \
                all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"] \
                or event == "-PURITY_INFO_GRAPH_SHOWING-" and purity_info_samples \
                or event == "-PURITY_INFO_PURITY_PEAK_LIST_TABLE-" and purit_info_values \
                or event == "-PURITY_INFO_DRAW_PEAKS-" and plot_style \
                or event == "-PURITY_INFO_PEAK_TABLE-" and plot_style:

            if event == "-PURITY_INFO_SAMPLE_BOX-":
                purity_info_samples = values["-PURITY_INFO_SAMPLE_BOX-"]
                lc_method = values["-PURITY_INFO_GRAPH_SHOWING-"]

                if lc_method == lc_graph_showing[3]:
                    start = values["-PURITY_INFO_RT_START-"]
                    end = values["-PURITY_INFO_RT_END-"]
                    peak = "None"
                    canvas_lines["peak_lines"][peak] = {"start": start, "end": end}

            elif event == "-PURITY_INFO_GRAPH_SHOWING-":
                lc_method = values["-PURITY_INFO_GRAPH_SHOWING-"]

            elif event == "-PURITY_INFO_PURITY_OVERVIEW_TABLE-" and values["-PURITY_INFO_GRAPH_SHOWING-"]:
                lc_method = lc_graph_showing[0]
                window["-PURITY_INFO_GRAPH_SHOWING-"].update(value=lc_method)
                print(f'All table data - purity info: {all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"]}')
                window["-PURITY_INFO_MZ-"].update(value=all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"][values[
                    "-PURITY_INFO_PURITY_OVERVIEW_TABLE-"][0]][1])
                purity_info_samples = [all_table_data["-PURITY_INFO_PURITY_OVERVIEW_TABLE-"][values[
                    "-PURITY_INFO_PURITY_OVERVIEW_TABLE-"][0]][0]]
                all_table_data["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"] = \
                    purity_peak_list_table_data[purity_info_samples[0]]
                window["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"].update(
                    values=all_table_data["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"])
                window["-PURITY_INFO_PEAK_LIST_SAMPLE_TEXT-"].update(value=purity_info_samples[0])
                window["-PURITY_INFO_PEAK_LIST_SAMPLE-"].update(value=purity_info_samples)

            elif event == "-PURITY_INFO_PURITY_PEAK_LIST_TABLE-" and values["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"]:
                lc_method = lc_graph_showing[3]
                window["-PURITY_INFO_GRAPH_SHOWING-"].update(value=lc_method)
                purity_info_samples = values["-PURITY_INFO_PEAK_LIST_SAMPLE-"]
                purity_info_rt_start = all_table_data["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"][
                    values["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"][0]][4]
                window["-PURITY_INFO_RT_START-"].update(value=purity_info_rt_start)

                purity_info_rt_end = all_table_data["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"][
                    values["-PURITY_INFO_PURITY_PEAK_LIST_TABLE-"][0]][5]
                window["-PURITY_INFO_RT_END-"].update(value=purity_info_rt_end)

                purity_info_samples = purity_info_samples.strip("',)")
                purity_info_samples = purity_info_samples.split("'")[1]
                purity_info_samples = [purity_info_samples]
                canvas_lines["peak_lines"][peak] = {"start": purity_info_rt_start, "end": purity_info_rt_end}
            # Set size of the canvas figure

            elif event == "-PURITY_INFO_DRAW_PEAKS-":
                for data in all_table_data["-PURITY_INFO_PEAK_TABLE-"]:
                    peak = data[1]
                    start = data[3]
                    end = data[4]
                    canvas_lines["peak_lines"][peak] = {"start": start, "end": end}
                    update_purity_info_peak_table = False

            elif event == "-PURITY_INFO_PEAK_TABLE-" and values["-PURITY_INFO_PEAK_TABLE-"]:
                purity_info_samples = [
                    all_table_data["-PURITY_INFO_PEAK_TABLE-"][values["-PURITY_INFO_PEAK_TABLE-"][0]][0]]
                peak = all_table_data["-PURITY_INFO_PEAK_TABLE-"][values["-PURITY_INFO_PEAK_TABLE-"][0]][1]
                start = all_table_data["-PURITY_INFO_PEAK_TABLE-"][values["-PURITY_INFO_PEAK_TABLE-"][0]][3]
                end = all_table_data["-PURITY_INFO_PEAK_TABLE-"][values["-PURITY_INFO_PEAK_TABLE-"][0]][4]
                canvas_lines["peak_lines"][peak] = {"start": start, "end": end}

                if values["-PURITY_INFO_RADIO_PEAKS_UV-"]:
                    lc_method = lc_graph_showing[0]
                elif values["-PURITY_INFO_RADIO_PEAKS_MS_SPECTRA-"]:
                    lc_method = lc_graph_showing[3]

                window["-PURITY_INFO_GRAPH_SHOWING-"].update(value=lc_method)
                window["-PURITY_INFO_RT_START-"].update(value=start)
                window["-PURITY_INFO_RT_END-"].update(value=end)
                update_purity_info_peak_table = False

            if values["-PURITY_INFO_DRAW_THRESHOLD-"]:
                uv_line = float(values["-PURITY_INFO_UV_THRESHOLD-"])
                canvas_lines["uv"] = uv_line

            fig_size = (7, 3)
            print("MISSING GUARD TO PREVENT DIFFERENT LENGTH DATA CRASHING THE PROGRAM!!!")

            purity_info_canvas = window["-PURITY_INFO_CANVAS-"]

            if values["-PURITY_INFO_MS_MODE_POS-"]:
                ms_mode = "ms_pos"
            elif values["-PURITY_INFO_MS_MODE_NEG-"]:
                ms_mode = "ms_neg"

            if not purity_info_rt_start:
                try:
                    purity_info_rt_start = float(values["-PURITY_INFO_RT_START-"])
                except ValueError:
                    purity_info_rt_start = None
            try:
                wavelength = float(values["-PURITY_INFO_WAVELENGTH-"])
            except ValueError:
                wavelength = None
            try:
                bin_numbers = int(values["-PURITY_INFO_BIN-"])
            except ValueError:
                bin_numbers = None
            try:
                mz_value = float(values["-PURITY_INFO_MZ-"])
            except ValueError:
                mz_value = None

            if purity_info_samples:
                plot_style = purity_plotting(lc_method, purity_data, purity_info_canvas, purity_info_samples, fig_size,
                                             ms_mode, purity_info_rt_start, purity_info_rt_end, wavelength, bin_numbers,
                                             mz_value, canvas_lines)

            else:
                fig = plt.gcf()
                plt.close(fig)
                try:
                    temp_purity_info_canvas.get_tk_widget().forget()
                except AttributeError:
                    print("Attribute Error for info canvas")

            if type(plot_style) == str:
                sg.PopupError(plot_style)
                plot_style = None
            else:

                if temp_purity_info_canvas is not None:
                    if not temp_purity_info_canvas == plot_style:
                        fig = plt.gcf()
                        plt.close(fig)
                        temp_purity_info_canvas.get_tk_widget().forget()
                    else:
                        plot_style.get_tk_widget().forget()
                if purity_info_samples:
                    plot_style.draw()

                    if toolbar:
                        toolbar.destroy()

                    toolbar = Toolbar(plot_style, window["-PURITY_INFO_CANVAS_TOOLBAR-"].TKCanvas)
                    toolbar.update()
                    plot_style.get_tk_widget().pack()
                try:
                    window.refresh()
                except AttributeError:
                    print("Canvas - AttributeError on window.refresh")
                temp_purity_info_canvas = plot_style
                temp_peak_table_data = []
                if len(purity_info_samples) > 1:
                    for sample in purity_info_samples:
                        for rows in peak_table_data[sample]:
                            temp_peak_table_data.append(rows)
                else:
                    try:
                        temp_peak_table_data = peak_table_data[purity_info_samples[0]]
                    except IndexError:
                        temp_peak_table_data = ""

                if update_purity_info_peak_table:
                    all_table_data["-PURITY_INFO_PEAK_TABLE-"] = temp_peak_table_data
                    window["-PURITY_INFO_PEAK_TABLE-"].update(values=all_table_data["-PURITY_INFO_PEAK_TABLE-"])
                all_table_data["-PURITY_INFO_RAW_DATA_TABLE-"] = purity_data[purity_info_samples[0]]["peak_table_raw"][
                    1]
                window["-PURITY_INFO_RAW_DATA_TABLE-"].update(values=all_table_data["-PURITY_INFO_RAW_DATA_TABLE-"])

                update_purity_info_peak_table = True
                purity_info_rt_start = None
                purity_info_rt_end = None
                canvas_lines = {
                    "uv": None,
                    "peak_lines": {}
                }

        # Sorting when clicking on Table headers. All tables should be in here execpt compound table, as it is a tree
        # Table...
        if isinstance(event, tuple):
            # TABLE CLICKED Event has value in format ('-TABLE=', '+CLICKED+', (row,col))
            # You can also call Table.get_last_clicked_position to get the cell clicked
            if event[0] in all_table_data and all_table_data[event[0]]:
                clicked_table = event[0]
                try:
                    search_reverse[clicked_table]
                except KeyError:
                    search_reverse[clicked_table] = {}
                if event[2][0] == -1 and event[2][1] != -1:  # Header was clicked and wasn't the "row" column
                    col_num_clicked = event[2][1]
                    try:
                        search_reverse[clicked_table][col_num_clicked]
                    except KeyError:
                        search_reverse[clicked_table][col_num_clicked] = False

                    new_table, search_reverse[clicked_table][col_num_clicked] = \
                        sort_table(all_table_data[clicked_table][0:][:], (col_num_clicked, 0),
                                   search_reverse[clicked_table][col_num_clicked])
                    window[clicked_table].update(new_table)
                    # all_table_data[clicked_table] = [all_table_data[clicked_table][0]] + new_table
                    all_table_data[clicked_table] = new_table


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


