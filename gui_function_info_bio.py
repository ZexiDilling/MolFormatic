import copy

from compound_plate_formatting import plate_layout_re_formate
from draw_basic import draw_plate
from helpter_functions import eval_guard_dict, eval_guard_list
from start_up_values import draw_tool_values


def colour_chooser_update(window, values, event):
    if event != "None":
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
                  "z_prime": z_prime,
                  "responsible": responsible,
                  "approval": approval,
                  "note": note,
                  "plate_layout": plate_layout,
                  "analysed_method": analysed_method}
    return plate_data


def _bio_info_update_window(analyse_method, plate_list, well_dict):

    pass


def bio_info_window_update(dbf, window, values):

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


def bio_info_plate_update(dbf, config, window, values, event, well_dict_bio_info):
    plate = values[event]
    if plate.casefold() != "all":
        plate_bio_data = bio_info_grab_data(dbf, plate)
        if plate_bio_data:
            well_dict_bio_info.clear()
            plate_layout_data = dbf.find_data_single_lookup("plate_layout", plate_bio_data["plate_layout"], "layout_name")
            try:
                plate_layout_data = plate_layout_data[0]
            except IndexError:
                return well_dict_bio_info
            else:
                plate_layout = eval(plate_layout_data[5])
                well_dict_bio_info = copy.deepcopy(plate_layout)
                well_dict_bio_info = plate_layout_re_formate(config, well_dict_bio_info)
                plate_size = plate_layout_data[2]
                archive = True
                gui_tab = "bio_exp"
                graph_bio = window["-BIO_INFO_CANVAS-"]

                well_dict_bio_info, _, _, _, _, off_set = draw_plate(config, graph_bio, plate_size, well_dict_bio_info,
                                                               gui_tab, archive)

                draw_tool_values["well_off_set"] = off_set

            return well_dict_bio_info
        else:
            return well_dict_bio_info
    else:
        return well_dict_bio_info


if __name__ == "__main__":
    pass
    # import configparser
    #
    # from database_functions import DataBaseFunctions, _get_list_of_names_from_database, \
    # _get_list_of_names_from_database_double_lookup
    #
    # config = configparser.ConfigParser()
    # config.read("config.ini")

    # dbf = DataBaseFunctions(config)
    # plate = "alpha_so80"
    # plate_data = bio_info_grab_data(dbf, window=None, values=None, event=None, plate=plate)
    # for data in plate_data:
    #     print(data)
    # print(plate_data["plate_layout"])
    # # values = {"-BIO_INFO_ASSAY_DROPDOWN-": "Alpha_so"}
    # # bio_info_window_update(dbf, None, values)


