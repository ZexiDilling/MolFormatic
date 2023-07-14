from gui_functions import grab_table_data, bio_compound_info_from_worklist, bio_data, get_plate_layout, draw_plate
import PySimpleGUI as sg

from bio_data_functions import org, norm, pora, pora_internal
from gui_popup import bio_data_approval_table


def testing(config, folder, default_plate_layout, bio_plate_report_setup, work_sheet, color_select):
    plate_list, archive_plates_dict = get_plate_layout(config)
    plate_to_layout = default_plate_layout
    # Sets values for different parametors
    bio_import_folder = folder
    # default_plate_layout = archive_plates_dict[values["-BIO_PLATE_LAYOUT-"]]
    default_plate_layout = default_plate_layout
    include_hits = False

    threshold = False

    hit_amount = int(10)

    include_smiles = True
    final_report_name = "testing"
    export_to_excel = False
    same_layout = True
    include_structure = False

    bio_export_folder = folder

    bio_sample_list = work_sheet
    bio_sample_list = bio_sample_list.split(";")
    bio_sample_dict = bio_compound_info_from_worklist(config, sg, bio_sample_list)

    analyse_method = "single point"
    add_compound_ids = False

    worked, all_plates_data, date, used_plates, plate_to_layout = \
        bio_data(config, bio_import_folder, plate_to_layout, archive_plates_dict,
                 bio_plate_report_setup, analyse_method, bio_sample_dict, bio_export_folder,
                 add_compound_ids, export_to_excel)

    assay_name = "Alpha_so"

    # # Grabs the value for the assay from the database
    assay_data_row = grab_table_data(config, table_name="assay", single_row=True,
                                     data_value=assay_name, headline="assay_name")

    assay_data = {"assay_name": assay_data_row[0][2],
                  "plate_layout": assay_data_row[0][4],
                  "z_prime_threshold": assay_data_row[0][7],
                  "hit_threshold": assay_data_row[0][8],
                  "plate_size": assay_data_row[0][9]}

    # Open a popup window where data can be checked before being added to the database.
    all_plates_data, plate_analyse_methods = \
        bio_data_approval_table(draw_plate, config, all_plates_data, assay_data, same_layout,
                                plate_to_layout, archive_plates_dict, color_select, bio_plate_report_setup)

if __name__ == "__main__":
    import configparser
    config = configparser.ConfigParser()
    config.read("config.ini")
    folder = r"C:\Users\phch\Desktop\test\test_data\data"
    work_sheet = r"C:\Users\phch\OneDrive - Danmarks Tekniske Universitet\Mapper\Python_data\alpha_SO\worklisted_used_with_good_data\Worklist_alpha_so_225_to_234.txt"

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
    default_plate_layout = "Daniells_Alpha_So"



    color_select = {}
    for keys in list(config["plate_colouring"].keys()):
        color_select[keys] = config["plate_colouring"][keys]

    testing(config, folder, default_plate_layout, bio_plate_report_setup, work_sheet, color_select)
