from bio_data_functions import norm, pora, pora_internal, org


def plate_colouring_fetch(config):
    # simple_settings
    plate_colouring = {
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
    return plate_colouring


def ms_settings_fetch(config):
    ms_setting = {
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
    return ms_setting


def bio_plate_report_setup_fetch(config):
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
    return bio_plate_report_setup


def bio_final_report_setup_fetch(config):
    temp_bio_report = {
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

    return temp_bio_report


def color_select_fetch(config):
    color_select = {}
    for keys in list(config["plate_colouring"].keys()):
        color_select[keys] = config["plate_colouring"][keys]
    return color_select

