from math import floor

from bio_heatmap import Heatmap
from database_functions import _get_list_of_names_from_database_double_lookup
from helpter_functions import int_guard
from start_up_values import window_1_plate_layout


def colour_target_update(window, values):
    if values["-PLATE_LAYOUT_COLOUR_CHOSE_TARGET-"] != "None":
        window["-PLATE_LAYOUT_COLOUR_CHOSE-"].update(
            button_color=values["-PLATE_LAYOUT_COLOUR_CHOSE_TARGET-"])


def plate_archive(window, values):
    pass


def plate_list_updater(dbf, window, values, event):
    if event == "-RECT_SAMPLE_TYPE-":
        target_window = "-ARCHIVE_PLATES-"
    elif event == "-BIO_ANALYSE_TYPE-":
        target_window = "-BIO_PLATE_LAYOUT-"

    table = "plate_layout"
    column_headline = "layout_name"
    limiting_value = values[event].casefold()
    limiting_header = "style"

    plate_list = _get_list_of_names_from_database_double_lookup(dbf, table, column_headline, limiting_value,
                                                                limiting_header)
    window[target_window].update(values=plate_list, value=plate_list[0])


def plate_layout_draw_groups(window, values, well_dict):
    if values["-PLATE_LAYOUT_DRAW_GROUPS-"] == "Dose Response":
        for tools in window_1_plate_layout["draw_tool_dict"]:
            if values[tools]:
                temp_draw_tool_tracker = tools

        window["-RECT_DOSE-"].update(value=True)
        total_sample_spots = 0
        for wells in well_dict:
            if well_dict[wells]["state"] == "sample":
                total_sample_spots += 1

        window["-SAMPLE_SPOTS-"].update(value=total_sample_spots)
    elif values["-PLATE_LAYOUT_DRAW_GROUPS-"] == "State":
        window[window_1_plate_layout["temp_draw_tool_tracker"]].update(value=True)


def _update_dose_tool(window, event, values, temp_sample_amount):

    try:
        int(temp_sample_amount)
    except ValueError:
        temp_sample_amount = 0

    sample_group_list = [number for number in range(1, temp_sample_amount + 1)]
    if sample_group_list:
        window["-SAMPLE_CHOOSER_DROPDOWN-"].update(values=sample_group_list, value=sample_group_list[0])

    try:
        int(values["-DOSE_REPLICATES-"])
    except ValueError:
        pass
    else:
        replicate_list = [x for x in range(1, int(values["-DOSE_REPLICATES-"]) + 1)]
        if replicate_list:
            window["-REPLICATE_CHOOSER_DROPDOWN-"].update(values=replicate_list, value=replicate_list[0])

    heatmap = Heatmap()
    colour_list = [values["-DOSE_COLOUR_LOW-"], values["-DOSE_COLOUR_HIGH-"]]

    dose_colour_dict = heatmap.poly_linear_gradient(colour_list, temp_sample_amount)
    for colour_index, colour in enumerate(dose_colour_dict["hex"]):
        if colour_index % 2 == 0:
            dose_colour_dict["hex"][colour_index] = heatmap.get_complementary(colour)

    return dose_colour_dict


def dose_sample_amount(window, event, values, dose_colour_dict):

    temp_replicates = int_guard(window, values["-DOSE_REPLICATES-"], 1)
    temp_sample_amount = int_guard(window, values["-DOSE_SAMPLE_AMOUNT-"], 0)

    if values["-DOSE_REPLICATES-"] is None or values["-DOSE_SAMPLE_AMOUNT-"] is None or \
            window_1_plate_layout["total_sample_spots"] is None or temp_replicates is None or \
            temp_sample_amount is None or window_1_plate_layout["total_sample_spots"] == 0 or \
            temp_sample_amount == 0 or temp_replicates == 0:
        return dose_colour_dict

    try:
        temp_dilutions = int(floor(window_1_plate_layout["total_sample_spots"] / (temp_sample_amount * temp_replicates)))
    except (ZeroDivisionError or TypeError):
        return dose_colour_dict

    else:
        window["-DOSE_DILUTIONS-"].update(value=temp_dilutions)
        empty_sample_spots = window_1_plate_layout["total_sample_spots"] - \
                             (temp_dilutions * temp_sample_amount * temp_replicates)

        window["-DOSE_EMPTY_SAMPLE_SPOTS-"].update(value=empty_sample_spots)

        dose_colour_dict = _update_dose_tool(window, event, values, temp_sample_amount)
    return dose_colour_dict


def dose_dilution_replicates(window, event, values, dose_colour_dict):

    temp_replicates = int_guard(window, values["-DOSE_REPLICATES-"], 1)
    temp_dilutions = int_guard(window, values["-DOSE_DILUTIONS-"], 0)

    if values["-DOSE_REPLICATES-"] is None or values["-DOSE_DILUTIONS-"] is None or window_1_plate_layout["total_sample_spots"] is None \
            or temp_dilutions == 0 or temp_replicates == 0 or window_1_plate_layout["total_sample_spots"] is None or temp_dilutions is None \
            or temp_replicates is None:
        return dose_colour_dict

    try:
        int(floor(window_1_plate_layout["total_sample_spots"] / (temp_dilutions * temp_replicates)))
    except (ZeroDivisionError or TypeError):
        return dose_colour_dict
    else:
        temp_sample_amount = int(floor(window_1_plate_layout["total_sample_spots"] / (temp_dilutions * temp_replicates)))
        window["-DOSE_SAMPLE_AMOUNT-"].update(value=temp_sample_amount)
        empty_sample_spots = window_1_plate_layout["total_sample_spots"] - (temp_dilutions * temp_sample_amount * temp_replicates)
        window["-DOSE_EMPTY_SAMPLE_SPOTS-"].update(value=empty_sample_spots)

        dose_colour_dict = _update_dose_tool(window, event, values, temp_sample_amount)

        return dose_colour_dict


def dose_colouring(window, event, values):
    if "HIGH" in event:
        window["-DOSE_COLOUR_BUTTON_HIGH-"].update(button_color=event)
    else:
        window["-DOSE_COLOUR_BUTTON_LOW-"].update(button_color=event)
    dose_colour_dict = _update_dose_tool(window, event, values, window_1_plate_layout["temp_sample_amount"])
    return dose_colour_dict


