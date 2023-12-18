import copy

from compound_plate_formatting import plate_layout_re_formate
from draw_basic import draw_plate
from helpter_functions import eval_guard_dict, eval_guard_list, plate_layout_to_state_dict
from start_up_values import draw_tool_values


def colour_chooser_update(window, values, event):
    if values[event] != "None":
        window_name = event.replace("_TARGET-", "-")
        window[window_name].update(button_color=values[event])


def bio_info_grab_data(dbf, plate):
    plate_data = dbf.find_data_single_lookup("biological_plate_data", plate, "plate_name")[0]
    assay_run = plate_data[2]
    raw_data = eval_guard_dict(plate_data[4])
    specific_plate_layout = eval_guard_dict(plate_data[5])
    if not specific_plate_layout:
        return None
    z_prime = plate_data[6]
    responsible = plate_data[7]
    approval = plate_data[8]
    note = plate_data[9]
    plate_layout = plate_data[10]
    skipped_wells = eval_guard_list(plate_data[11])
    if type(skipped_wells) != list:
        return None
    analysed_method = plate_data[12]

    plate_data = {"assay_run": assay_run,
                  "raw_data": raw_data,
                  "specific_plate_layout": specific_plate_layout,
                  "z_prime": z_prime,
                  "responsible": responsible,
                  "approval": approval,
                  "note": note,
                  "plate_layout": plate_layout,
                  "analysed_method": analysed_method,
                  "skipped_wells": skipped_wells}
    return plate_data


def _bio_info_update_window(analyse_method, plate_list, well_dict):

    pass


def bio_info_window_update(dbf, window, values, assay=None):
    if not assay:
        assay = values["-BIO_INFO_ASSAY_DROPDOWN-"]

    if assay:
        runs = ["All"]
        run_row_data = dbf.find_data_single_lookup("assay_runs", assay, "assay_name")

        for rows in run_row_data:
            runs.append(rows[1])

        skipped_runs, plates = bio_info_plate_list_update(dbf, window, values, runs)

        for skipping_counter in reversed(skipped_runs):
            runs.pop(skipping_counter)

        window["-BIO_INFO_RUN_DROPDOWN-"].update(values=runs, value=runs[0])


def bio_info_plate_list_update(dbf, window, values, runs):

    if not runs:
        if values["-BIO_INFO_RUN_DROPDOWN-"].casefold() == "all":
            bio_info_window_update(dbf, window, values)
            return
        else:
            runs = ["All", values["-BIO_INFO_RUN_DROPDOWN-"]]

    only_approved = values["-BIO_INFO_APPROVED_CHECK-"]
    method = values["-BIO_INFO_ANALYSE_METHOD-"]
    plates = ["All"]
    skipped_runs = []
    for counter, temp_run in enumerate(runs):
        if counter > 0:
            if only_approved:
                temp_row_data = dbf.find_data_double_lookup("biological_plate_data",
                                                            temp_run, "True", "assay_run", "approval")
            else:
                temp_row_data = dbf.find_data_single_lookup("biological_plate_data", temp_run, "assay_run")
            try:
                temp_row_data[0]
            except IndexError:
                skipped_runs.append(counter)
            else:
                for rows in temp_row_data:
                    if method == "All":
                        plates.append(rows[3])
                    else:
                        if rows[12].casefold() == method.casefold():
                            plates.append(rows[3])

    window["-BIO_INFO_PLATES_DROPDOWN-"].update(values=plates, value=[plates[0]])

    return skipped_runs, plates


def bio_info_plate_update(dbf, config, window, values, event, well_dict_bio_info, plate=None):
    if not plate:
        plate = values[event]
    if plate.casefold() != "all":
        plate_bio_data = bio_info_grab_data(dbf, plate)
        if plate_bio_data:

            well_dict_bio_info.clear()
            plate_layout_data = dbf.find_data_single_lookup("plate_layout", plate_bio_data["plate_layout"],
                                                            "layout_name")
            try:
                plate_layout_data = plate_layout_data[0]
            except IndexError:
                return well_dict_bio_info
            else:

                plate_layout = eval(plate_layout_data[5])
                well_dict_bio_info = copy.deepcopy(plate_layout)
                well_dict_bio_info = plate_layout_re_formate(config, well_dict_bio_info)
                plate_size = plate_layout_data[2]

                gui_tab = "bio_exp"
                graph_bio = window["-BIO_INFO_CANVAS-"]

            if values["-BIO_INFO_MAPPING-"] == "State Mapping":
                mapping = None
                archive = True
                temp_well_ditch = well_dict_bio_info
                skipped_wells = None
                state_dict = None
            elif values["-BIO_INFO_MAPPING-"] == "Heatmap":
                temp_well_ditch = plate_bio_data["specific_plate_layout"]
                archive = False
                skipped_wells = plate_bio_data["skipped_wells"]
                state_dict = plate_layout_to_state_dict(plate_layout)
                mapping = {
                    "mapping": "Heatmap",
                    "colours": {"low": [values["-BIO_INFO_HEATMAP_LOW_COLOUR_TARGET-"],
                                        values["-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-"]],
                                "high": [values["-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-"],
                                         values["-BIO_INFO_HEATMAP_HIGH_COLOUR_TARGET-"]]},
                    "percentile": {"low": float(values["-BIO_INFO_HEAT_PERCENTILE_LOW-"]),
                                   "mid": float(values["-BIO_INFO_HEAT_PERCENTILE_MID-"]),
                                   "high": float(values["-BIO_INFO_HEAT_PERCENTILE_HIGH-"])}
                }
            elif values["-BIO_INFO_MAPPING-"] == "Hit Map":
                temp_well_ditch = plate_bio_data["specific_plate_layout"]
                archive = False
                skipped_wells = plate_bio_data["skipped_wells"]
                state_dict = plate_layout_to_state_dict(plate_layout)
                th_use = {
                    "TH_1": False,
                    "TH_2": False,
                    "TH_3": False,
                    "TH_4": False,
                    "TH_5": False,
                    "TH_6": False,
                    "TH_7": False,
                    "TH_8": False,
                    "TH_9": False,
                    "TH_10": False}
                for temp_th_value in th_use:
                    if values[f"-BIO_INFO_PORA_{temp_th_value}_MIN_HIT_THRESHOLD-"] != \
                            values[f"-BIO_INFO_PORA_{temp_th_value}_MAX_HIT_THRESHOLD-"]:
                        th_use[temp_th_value] = True
                mapping = {
                    "mapping": "Hit Map",
                    "bins": {"th_1": {
                        "use": th_use["TH_1"],
                        "min": float(values["-BIO_INFO_PORA_TH_1_MIN_HIT_THRESHOLD-"]),
                        "max": float(values["-BIO_INFO_PORA_TH_1_MAX_HIT_THRESHOLD-"]),
                        "colour": values["-BIO_INFO_HIT_MAP_TH_1_COLOUR_TARGET-"]},
                        "th_2": {
                            "use": th_use["TH_2"],
                            "min": float(values["-BIO_INFO_PORA_TH_2_MIN_HIT_THRESHOLD-"]),
                            "max": float(values["-BIO_INFO_PORA_TH_2_MAX_HIT_THRESHOLD-"]),
                            "colour": values["-BIO_INFO_HIT_MAP_TH_2_COLOUR_TARGET-"]},
                        "th_3": {
                            "use": th_use["TH_3"],
                            "min": float(values["-BIO_INFO_PORA_TH_3_MIN_HIT_THRESHOLD-"]),
                            "max": float(values["-BIO_INFO_PORA_TH_3_MAX_HIT_THRESHOLD-"]),
                            "colour": values["-BIO_INFO_HIT_MAP_TH_3_COLOUR_TARGET-"]},
                        "th_4": {
                            "use": th_use["TH_4"],
                            "min": float(values["-BIO_INFO_PORA_TH_4_MIN_HIT_THRESHOLD-"]),
                            "max": float(values["-BIO_INFO_PORA_TH_4_MAX_HIT_THRESHOLD-"]),
                            "colour": values["-BIO_INFO_HIT_MAP_TH_4_COLOUR_TARGET-"]},
                        "th_5": {
                            "use": th_use["TH_5"],
                            "min": float(values["-BIO_INFO_PORA_TH_5_MIN_HIT_THRESHOLD-"]),
                            "max": float(values["-BIO_INFO_PORA_TH_5_MAX_HIT_THRESHOLD-"]),
                            "colour": values["-BIO_INFO_HIT_MAP_TH_5_COLOUR_TARGET-"]},
                        "th_6": {
                            "use": th_use["TH_6"],
                            "min": float(values["-BIO_INFO_PORA_TH_6_MIN_HIT_THRESHOLD-"]),
                            "max": float(values["-BIO_INFO_PORA_TH_6_MAX_HIT_THRESHOLD-"]),
                            "colour": values["-BIO_INFO_HIT_MAP_TH_6_COLOUR_TARGET-"]},
                        "th_7": {
                            "use": th_use["TH_7"],
                            "min": float(values["-BIO_INFO_PORA_TH_7_MIN_HIT_THRESHOLD-"]),
                            "max": float(values["-BIO_INFO_PORA_TH_7_MAX_HIT_THRESHOLD-"]),
                            "colour": values["-BIO_INFO_HIT_MAP_TH_7_COLOUR_TARGET-"]},
                        "th_8": {
                            "use": th_use["TH_8"],
                            "min": float(values["-BIO_INFO_PORA_TH_8_MIN_HIT_THRESHOLD-"]),
                            "max": float(values["-BIO_INFO_PORA_TH_8_MAX_HIT_THRESHOLD-"]),
                            "colour": values["-BIO_INFO_HIT_MAP_TH_8_COLOUR_TARGET-"]},
                        "th_9": {
                            "use": th_use["TH_9"],
                            "min": float(values["-BIO_INFO_PORA_TH_9_MIN_HIT_THRESHOLD-"]),
                            "max": float(values["-BIO_INFO_PORA_TH_9_MAX_HIT_THRESHOLD-"]),
                            "colour": values["-BIO_INFO_HIT_MAP_TH_9_COLOUR_TARGET-"]},
                        "th_10": {
                            "use": th_use["TH_10"],
                            "min": float(values["-BIO_INFO_PORA_TH_10_MIN_HIT_THRESHOLD-"]),
                            "max": float(values["-BIO_INFO_PORA_TH_10_MAX_HIT_THRESHOLD-"]),
                            "colour": values["-BIO_INFO_HIT_MAP_TH_10_COLOUR_TARGET-"]}
                    }
                }
            else:
                return well_dict_bio_info

            well_dict_bio_info, _, _, _, _, off_set = draw_plate(config, graph_bio, plate_size,
                                                                 temp_well_ditch, gui_tab, archive, mapping=mapping,
                                                                 skipped_well=skipped_wells, state_dict=state_dict)
            draw_tool_values["well_off_set"] = off_set
            return well_dict_bio_info
        else:
            return well_dict_bio_info
    else:
        return well_dict_bio_info


def bio_info_canvas_clicked(window, values, event, bio_info_clicked):
    if bio_info_clicked:
        bio_info_clicked = False
        bio_info_well = values["-BIO_INFO_CANVAS-"]
        print(f"Canvas on BIO info was clicked. This is the value: {bio_info_well}")
    return bio_info_clicked

if __name__ == "__main__":
    # pass
    import configparser
    #
    from database_functions import DataBaseFunctions
    # _get_list_of_names_from_database_double_lookup
    #
    config = configparser.ConfigParser()
    config.read("config.ini")
    dbf = DataBaseFunctions(config)
    temp_run = "Alpha_so_run_3"
    temp_row_data = dbf.find_data_double_lookup("biological_plate_data",
                                                temp_run, "True", "assay_run", "approval")
    print(temp_row_data)

    # plate = "alpha_so80"
    # plate_data = bio_info_grab_data(dbf, window=None, values=None, event=None, plate=plate)
    # for data in plate_data:
    #     print(data)
    # print(plate_data["plate_layout"])
    # # values = {"-BIO_INFO_ASSAY_DROPDOWN-": "Alpha_so"}
    # # bio_info_window_update(dbf, None, values)


