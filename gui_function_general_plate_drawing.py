import copy
from PySimpleGUI import PopupError, PopupGetText, PopupGetFolder, Popup, PopupYesNo

from draw_basic import draw_plate
from gui_popup import popup_three_box_solution
from lcms_functions import plate_layout_to_excel
from database_functions import update_database, _get_list_of_names_from_database, rename_record_in_the_database
from helpter_functions import int_guard
from compound_plate_formatting import plate_layout_re_formate
from start_up_values import window_1_plate_layout, draw_tool_values, clm_to_row_converter, plate_type_count, \
    plate_layout_dropdowns


def _on_up_grab_graph_list(values, temp_x, temp_y):
    # makes a set, for adding wells, to avoid duplicates
    graphs_list = set()

    # goes over the coordinates and if they are within the bounds of the plate
    # if that is the case, then the figure for that location is added the set
    for index_x, cords_x in enumerate(temp_x):
        for index_y, cords_y in enumerate(temp_y):
            if draw_tool_values["min_x"] <= temp_x[index_x] <= draw_tool_values["max_x"] and \
                    draw_tool_values["min_y"] <= temp_y[index_y] <= draw_tool_values["max_y"]:
                graphs_list.add(
                    window_1_plate_layout["graph_plate"].get_figures_at_location((temp_x[index_x], temp_y[index_y]))[0])

    # colours the wells in different colour, depending on if they are samples or blanks
    if values["-DOSE_VERTICAL-"]:
        new_graphs_list = []
        temp_graph_list = []
        temp_converter = clm_to_row_converter[values["-PLATE-"]]

        # Change the wells numbers to what they would be, if they counted differently.
        temp_well_saver = {}
        for wells in graphs_list:
            temp_well = wells % plate_type_count[window_1_plate_layout["plate_type"]]
            converted_well = temp_converter[temp_well + 1]
            temp_graph_list.append(converted_well)
            temp_well_saver[converted_well] = wells

        # sorts the list
        temp_graph_list = sorted(temp_graph_list)

        for wells in temp_graph_list:
            new_graphs_list.append(temp_well_saver[wells])

    else:
        new_graphs_list = sorted(graphs_list)

    return new_graphs_list


def _on_up_x_y(start_point, end_point, x, y):
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

    return temp_x, temp_y


def _on_up_well_handler(values, well_dict, new_graphs_list, temp_sample_amount, temp_dilutions, temp_sample_group,
                        replicate_loop, replicates, dose_colour_dict, concentration_count, colour_select):

    for well_counter, wells in enumerate(new_graphs_list):
        temp_well = wells
        wells = wells - draw_tool_values["well_off_set"]

        try:
            int(values["-PLATE_LAYOUT_GROUP-"].strip("Group "))
        except ValueError:
            group_number = 0
        else:
            group_number = int(values["-PLATE_LAYOUT_GROUP-"].strip("Group "))

        if window_1_plate_layout["temp_draw_tool"] == "dose":
            try:
                well_dict[wells]["state"]
            except KeyError:
                PopupError("Sorry, can't deal with different plate types in a row")
                break
            else:

                if well_dict[wells]["state"] == "sample":
                    # if temp_replicate > replicates:
                    #     replicates = 1
                    #     temp_replicate = 1

                    if well_counter % temp_dilutions == 0 and well_counter > 0:
                        temp_replicate += 1
                        if temp_replicate > replicates:
                            temp_sample_group += 1
                            temp_replicate = 1

                        concentration_count = 1
                        if temp_sample_group > (temp_sample_amount * replicates):

                            colour = colour_select["sample"]
                            window_1_plate_layout["graph_plate"].Widget.itemconfig(wells, fill=colour)

                            well_dict[wells]["colour"] = colour
                            well_dict[wells]["group"] = group_number
                            well_dict[wells]["replicate"] = 0
                            well_dict[wells]["concentration"] = 0
                            continue

                    else:
                        concentration_count += 1

                    try:
                        dose_colour_dict["hex"][(temp_sample_group - 1)]
                    except (IndexError or TypeError):
                        print(f"Index error on dose colour dict for this group: {temp_sample_group}")
                        colour = colour_select["sample"]
                        window_1_plate_layout["graph_plate"].Widget.itemconfig(wells, fill=colour)

                        well_dict[wells]["colour"] = colour
                        well_dict[wells]["group"] = group_number
                        well_dict[wells]["replicate"] = 0
                        well_dict[wells]["concentration"] = 0
                        continue
                    else:
                        colour = dose_colour_dict["hex"][(temp_sample_group - 1)]

                    window_1_plate_layout["graph_plate"].Widget.itemconfig(wells, fill=colour)
                    try:
                        well_dict[wells]["colour"]
                    except IndexError as e:
                        print(f"error - {e} - for:")
                        print(well_dict)
                    else:
                        well_dict[wells]["colour"] = colour
                    well_dict[wells]["group"] = temp_sample_group
                    well_dict[wells]["replicate"] = temp_replicate
                    well_dict[wells]["concentration"] = concentration_count

        else:
            colour = colour_select[draw_tool_values["temp_draw_tool"]]
            well_state = draw_tool_values["temp_draw_tool"]
            if colour == "paint":
                colour = values["-PLATE_LAYOUT_COLOUR_CHOSE_TARGET-"]
            window_1_plate_layout["graph_plate"].Widget.itemconfig(temp_well, fill=colour)
            well_dict[wells]["colour"] = colour
            well_dict[wells]["state"] = well_state
            try:
                well_dict[wells]["group"]
            except KeyError:
                pass
            else:
                well_dict[wells]["group"] = group_number
                well_dict[wells]["replicate"] = 0
                well_dict[wells]["concentration"] = 0

    return well_dict


def on_up(window, values, well_dict, dose_colour_dict, colour_select):

    try:
        draw_tool_values["temp_selector"] and draw_tool_values["plate_active"]
    except KeyError:
        pass
    else:
        if draw_tool_values["temp_selector"] and draw_tool_values["plate_active"]:

            # Dose response values
            if draw_tool_values["temp_draw_tool"] == "dose":
                temp_dilutions = int_guard(window, values["-DOSE_DILUTIONS-"], 0)
                temp_sample_amount = int_guard(window, values["-DOSE_SAMPLE_AMOUNT-"], 0)
                temp_sample_group = int_guard(window, values["-SAMPLE_CHOOSER_DROPDOWN-"], 1)
                replicates = int_guard(window, values["-DOSE_REPLICATES-"], 1)
                concentration_count = 0
                replicate_loop = int_guard(window, values["-REPLICATE_CHOOSER_DROPDOWN-"], 1)
            else:
                temp_dilutions = temp_sample_amount = temp_sample_group = \
                    replicates = replicate_loop = concentration_count = 0

            if draw_tool_values["temp_draw_tool"] == "dose" and temp_sample_amount == 0 and temp_dilutions == 0:
                pass
            else:

                temp_x, temp_y = _on_up_x_y(draw_tool_values["start_point"], draw_tool_values["end_point"],
                                            draw_tool_values["x"], draw_tool_values["y"])

                new_graphs_list = _on_up_grab_graph_list(values, temp_x, temp_y)

                well_dict = _on_up_well_handler(values, well_dict, new_graphs_list, temp_sample_amount, temp_dilutions,
                                                temp_sample_group, replicate_loop, replicates, dose_colour_dict,
                                                concentration_count, colour_select)

                # Reset or sets the sample tool counter
                if temp_sample_group >= temp_sample_amount:
                    temp_sample_group = 1
                else:
                    temp_sample_group += 1
                window["-SAMPLE_CHOOSER_DROPDOWN-"].update(value=temp_sample_group)

                if window_1_plate_layout["temp_draw_tool"] == "dose":
                    window["-RECT_SAMPLE_TYPE-"].update(value="Dose Response")

    # deletes the rectangle used for selection
    if draw_tool_values["prior_rect"]:
        window_1_plate_layout["graph_plate"].delete_figure(draw_tool_values["prior_rect"])

    # reset everything
    draw_tool_values["start_point"], draw_tool_values["end_point"] = None, None
    draw_tool_values["dragging"] = False
    draw_tool_values["prior_rect"] = None
    draw_tool_values["temp_selector"] = False
    draw_tool_values["temp_draw_tool"] = None

    return well_dict


def on_move(window, values, graph_bio_exp, well_dict, well_dict_bio_info):

    try:
        values["-BIO_INFO_CANVAS-"]
    except KeyError:
        pass
    else:
        if values["-BIO_INFO_CANVAS-"][0] and values["-BIO_INFO_CANVAS-"][1]:
            try:
                temp_well_bio_info = graph_bio_exp.get_figures_at_location(values['-BIO_INFO_CANVAS-'])[0]
                temp_well_bio_info = temp_well_bio_info - draw_tool_values["well_off_set"]
            except IndexError:
                temp_well_id_bio_info = ""
            else:
                try:
                    temp_well_id_bio_info = well_dict_bio_info[temp_well_bio_info]["well_id"]
                except KeyError:
                    temp_well_id_bio_info = ""
            window["-INFO_BIO_GRAPH_TARGET-"].update(value=f"{temp_well_id_bio_info}")

    try:
        values["-RECT_BIO_CANVAS-"]
    except KeyError:
        pass
    else:
        if values["-RECT_BIO_CANVAS-"][0] and values["-RECT_BIO_CANVAS-"][1]:
            try:
                window_1_plate_layout["graph_plate"].get_figures_at_location(values['-RECT_BIO_CANVAS-'])[0]
            except (IndexError or KeyError) as error:
                # print(f"Canvas Error: {error}")
                temp_well_id = ""
                temp_well_group = ""
                temp_rep = ""
                temp_conc = ""
            else:
                temp_well = window_1_plate_layout["graph_plate"].get_figures_at_location(values['-RECT_BIO_CANVAS-'])[0]
                look_up_well = temp_well - draw_tool_values["well_off_set"]
                try:
                    well_dict[look_up_well]["well_id"]
                except (KeyError or UnboundLocalError):
                    temp_well_id = ""
                    temp_well_group = ""
                else:
                    temp_well_id = well_dict[look_up_well]["well_id"]
                    temp_well_group = well_dict[look_up_well]["group"]

                try:
                    well_dict[look_up_well]["replicate"]
                except KeyError:
                    temp_rep = ""
                    temp_conc = ""
                else:
                    temp_rep = well_dict[look_up_well]["replicate"]
                    temp_conc = well_dict[look_up_well]["concentration"]

            window["-CANVAS_INFO_WELL-"].update(value=f"Well: {temp_well_id}")
            window["-CANVAS_INFO_GROUP-"].update(value=f"Group: {temp_well_group}")
            window["-CANVAS_INFO_CONC-"].update(value=f"conc: {temp_conc}")
            window["-CANVAS_INFO_REP-"].update(value=f"rep: {temp_rep}")


def draw_layout(dbf, config, window, values, well_dict):
    well_dict.clear()
    # sets the size of the well for when it draws the plate
    graph = window_1_plate_layout["graph_plate"]
    plate_type = values["-PLATE-"]
    archive_plates = values["-ARCHIVE-"]
    gui_tab = "plate_layout"
    sample_type = values["-RECT_SAMPLE_TYPE-"]

    if values["-ARCHIVE-"]:

        try:
            plate_data = dbf.find_data_single_lookup("plate_layout", values["-ARCHIVE_PLATES-"], "layout_name")[0]
        except KeyError:
            window["-ARCHIVE-"].update(False)
            values["-ARCHIVE-"] = False
        else:
            well_dict = eval(plate_data[5])
            well_dict = plate_layout_re_formate(config, well_dict)
            window["-PLATE-"].update(plate_data[2])
            plate_type = plate_data[2]
            sample_type = plate_data[4]

    well_dict, min_x, min_y, max_x, max_y, off_set = draw_plate(config, graph, plate_type, well_dict, gui_tab,
                                                                archive_plates, sample_layout=sample_type)
    plate_active = True

    draw_tool_values["min_x"] = min_x
    draw_tool_values["min_y"] = min_y
    draw_tool_values["max_x"] = max_x
    draw_tool_values["max_y"] = max_y
    draw_tool_values["well_off_set"] = off_set
    draw_tool_values["plate_type"] = plate_type
    draw_tool_values["plate_active"] = plate_active
    draw_tool_values["last_plate"] = plate_type

    return well_dict


def export_layout(config, window, values, well_dict):
    if not well_dict:
        PopupError("Please create a layout to Export")
    name = PopupGetText("Name the file")
    if name:
        folder = PopupGetFolder("Choose save location")
        if folder:
            plate_layout_to_excel(config, well_dict, name, folder)

            Popup("Done")


def save_layout(dbf, config, window, values, well_dict):
    print(well_dict)
    if not well_dict:
        PopupError("Please create a layout to save")
    elif any("paint" in stuff.values() for stuff in well_dict.values()):
        PopupError("Can't save layout with paint as well states")
    else:
        name = "Overwrite or New layout?"
        box_1 = "Overwrite"
        box_2 = "New"
        question=f"Do you wish to overwrite {values['-ARCHIVE_PLATES-']}-layout?"
        if values["-ARCHIVE_PLATES-"]:
            overwrite_check = popup_three_box_solution(
                config, name=name, box_1=box_1, box_2=box_2, question=question)
        else:
            overwrite_check = box_2

        if overwrite_check == box_2:
            sample_type_check = PopupYesNo(f"You are about to save the layout with the style: \n "
                                           f"{values['-RECT_SAMPLE_TYPE-']} \n"
                                           f"Is that fine?")

            if sample_type_check.casefold() == "yes":
                # ToDo For Dose Response, add som guard functions for checking samples and so on
                temp_well_dict = {}
                temp_dict_name = PopupGetText("Name plate layout")
                if temp_dict_name:
                    for index, well_counter in enumerate(well_dict):
                        temp_well_dict[index + 1] = copy.deepcopy(well_dict[well_counter])

                    # saves the layout to the Database
                    # setting up the data for importing the new plate_layout to the database
                    temp_table = "plate_layout"
                    # ToDo add plate model ??
                    temp_plate_layout_data = {
                        "layout_name": temp_dict_name,
                        "plate_size": values["-PLATE-"],
                        "plate_mode": "placeholder",
                        "style": values["-RECT_SAMPLE_TYPE-"].casefold(),
                        "plate_layout": f"{temp_well_dict}"
                    }
                    update_database(config, temp_table, temp_plate_layout_data)
                else:
                    return
        elif overwrite_check == box_1:
            temp_well_dict = {}

            for index, well_counter in enumerate(well_dict):
                temp_well_dict[index + 1] = copy.deepcopy(well_dict[well_counter])
            
            table = "plate_layout"
            headline_for_changing_value = "plate_layout"
            headline_for_indicator_value = "layout_name"
            indicator_value = values['-ARCHIVE_PLATES-']
            new_value = f"{temp_well_dict}"
            rename_record_in_the_database(config, table, headline_for_changing_value, headline_for_indicator_value,
                                          indicator_value, new_value)
        else:
            return

        update_plate_layout_dropdowns(window, dbf)


def delete_layout(dbf, window, values):
    # ToDO FIX THIS so it can delete layouts
    if not values["-ARCHIVE_PLATES-"]:
        PopupError("Please select a layout to delete")
    else:
        # Set up values for the database, and deletes the record
        table = "plate_layout"
        headline = "layout_name"
        data_value = values["-ARCHIVE_PLATES-"]
        dbf.delete_records(table, headline, data_value)

        update_plate_layout_dropdowns(window, dbf)


def rename_layout(dbf, config, window, values):
    if not values["-ARCHIVE_PLATES-"]:
        PopupError("Please select a layout to rename")
    else:
        temp_dict_name = PopupGetText("Name plate layout")
        if temp_dict_name:
            # Updates the database with new values
            table = "plate_layout"
            headline_for_changing_value = "layout_name"
            headline_for_indicator_value = headline_for_changing_value
            indicator_value = values["-ARCHIVE_PLATES-"]
            new_value = temp_dict_name

            rename_record_in_the_database(config, table, headline_for_changing_value, headline_for_indicator_value,
                                          indicator_value, new_value)

            update_plate_layout_dropdowns(window, dbf)


def bio_canvas_group_labeling(event, values, well_dict):
    """
    Changes the group liable of a well by right-clicking on it
    :param event:
    :param values:
    :param well_dict:
    :return:
    """
    group_name = event
    temp_well = window_1_plate_layout["graph_plate"].get_figures_at_location(values['-RECT_BIO_CANVAS-'])[0]
    look_up_well = temp_well - draw_tool_values["well_off_set"]
    try:
        well_dict[look_up_well]["well_id"]
    except (KeyError or UnboundLocalError):
        return well_dict
    else:

        well_dict[look_up_well]["group"] = group_name.strip("Group ")

    return well_dict


def bio_canvas(values, draw_tool=None):

    if draw_tool:
        draw_tool_value = draw_tool
    else:
        draw_tool_value = draw_tool_values

    draw_tool_value["x"], draw_tool_value["y"] = values["-RECT_BIO_CANVAS-"]
    if not draw_tool_value["dragging"]:
        draw_tool_value["start_point"] = (draw_tool_value["x"], draw_tool_value["y"])
        draw_tool_value["dragging"] = True
    else:
        draw_tool_value["end_point"] = (draw_tool_value["x"], draw_tool_value["y"])
    if draw_tool_value["prior_rect"]:
        window_1_plate_layout["graph_plate"].delete_figure(draw_tool_value["prior_rect"])

    # Choosing which tool to pain the plate with.
    if None not in (draw_tool_value["start_point"], draw_tool_value["end_point"]):
        for temp_draw_value in window_1_plate_layout["draw_tool_dict"]:
            if values[temp_draw_value]:
                draw_tool_value["temp_draw_tool"] = window_1_plate_layout["draw_tool_dict"][temp_draw_value]
        draw_tool_value["temp_selector"] = True
        draw_tool_value["prior_rect"] = window_1_plate_layout["graph_plate"].\
            draw_rectangle(draw_tool_value["start_point"], draw_tool_value["end_point"],
                           fill_color="", line_color="white")


def update_plate_layout_dropdowns(window, dbf):
    plate_list = _get_list_of_names_from_database(dbf, "plate_layout", "layout_name")

    for dropdown in plate_layout_dropdowns:
        window[dropdown].update(values=sorted(plate_list), value=plate_list[0])


if __name__ == "__main__":
    pass
    # import configparser
    # from database_controller import DataBaseFunctions
    # config = configparser.ConfigParser()
    # config.read("config.ini")
    # dbf = DataBaseFunctions(config)
    # plate_list = _get_list_of_names_from_database(dbf, "plate_layout", "layout_name")
