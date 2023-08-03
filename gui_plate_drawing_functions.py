import copy
from math import floor

from heatmap import Heatmap
from gui_functions import get_plate_layout, draw_plate, plate_layout_to_excel, update_database, \
    delete_records_from_database, rename_record_in_the_database
from plate_formatting import plate_layout_re_formate


def _on_up_grab_graph_list(values, clm_to_row_converter, plate_type_count, plate_type,
                           temp_x, temp_y, min_x, max_x, min_y, max_y, graph_plate):
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
    if values["-DOSE_VERTICAL-"]:
        new_graphs_list = []
        temp_graph_list = []
        temp_converter = clm_to_row_converter[values["-PLATE-"]]

        # Change the wells numbers to what they would be, if they counted differently.
        temp_well_saver = {}
        print(plate_type)
        for wells in graphs_list:
            temp_well = wells % plate_type_count[plate_type]
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


def _on_up_well_handler(sg, values, well_dict, graph_plate, new_graphs_list, temp_draw_tool, color_select,
                        temp_sample_amount, temp_dilutions, temp_sample_group, temp_replicate, replicates,
                        dose_colour_dict, concentration_count):

    print(f"temp_sample_group - {temp_sample_group}")
    print(f"temp_sample_amount - {temp_sample_amount}")
    print(f"replicates - {replicates}")
    for well_counter, wells in enumerate(new_graphs_list):
        if temp_draw_tool == "dose":
            try:
                well_dict[wells]["state"]
            except KeyError:
                sg.PopupError("Sorry, can't deal with different plate types in a row")
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

                            colour = color_select["sample"]
                            graph_plate.Widget.itemconfig(wells, fill=colour)

                            well_dict[wells]["colour"] = colour
                            well_dict[wells]["group"] = 0
                            well_dict[wells]["replicate"] = 0
                            well_dict[wells]["concentration"] = 0
                            continue

                    else:
                        concentration_count += 1

                    try:
                        dose_colour_dict["hex"][(temp_sample_group - 1)]
                    except (IndexError or TypeError):
                        print(f"Index error on dose colour dict for this group: {temp_sample_group}")
                        colour = color_select["sample"]
                        graph_plate.Widget.itemconfig(wells, fill=colour)

                        well_dict[wells]["colour"] = colour
                        well_dict[wells]["group"] = 0
                        well_dict[wells]["replicate"] = 0
                        well_dict[wells]["concentration"] = 0
                        continue
                    else:
                        colour = dose_colour_dict["hex"][(temp_sample_group - 1)]

                    graph_plate.Widget.itemconfig(wells, fill=colour)
                    try:
                        well_dict[wells]["colour"]
                    except IndexError:
                        print(well_dict)
                    else:
                        well_dict[wells]["colour"] = colour
                    well_dict[wells]["group"] = temp_sample_group
                    well_dict[wells]["replicate"] = temp_replicate
                    well_dict[wells]["concentration"] = concentration_count

        else:
            colour = color_select[temp_draw_tool]
            well_state = temp_draw_tool
            if colour == "paint":
                colour = values["-PLATE_LAYOUT_COLOUR_CHOSE_TARGET-"]
            graph_plate.Widget.itemconfig(wells, fill=colour)
            well_dict[wells]["colour"] = colour
            well_dict[wells]["state"] = well_state

    return well_dict


def on_up(config, sg, window, event, values, temp_selector, plate_active, start_point, end_point, x, y, min_x, max_x,
          min_y, max_y, graph_plate, clm_to_row_converter, plate_type_count, plate_type, color_select, temp_draw_tool,
          well_dict, prior_rect, dose_colour_dict):
    if temp_selector and plate_active:

        # Dose response values
        if temp_draw_tool == "dose":
            temp_dilutions = _int_guard(window, values["-DOSE_DILUTIONS-"], 0)
            temp_sample_amount = _int_guard(window, values["-DOSE_SAMPLE_AMOUNT-"], 0)
            temp_sample_group = _int_guard(window, values["-SAMPLE_CHOOSER_DROPDOWN-"], 1)
            replicates = _int_guard(window, values["-DOSE_REPLICATES-"], 1)
            concentration_count = 0
            replicate_loop = 1
        else:
            temp_dilutions = temp_sample_amount = temp_sample_group = \
                replicates = replicate_loop = concentration_count = 0

        if temp_draw_tool == "dose" and temp_sample_amount == 0 or temp_draw_tool == "dose" and temp_dilutions == 0:
            pass
        else:

            temp_x, temp_y = _on_up_x_y(start_point, end_point, x, y)

            new_graphs_list = _on_up_grab_graph_list(values, clm_to_row_converter, plate_type_count, plate_type,
                                                     temp_x, temp_y, min_x, max_x, min_y, max_y, graph_plate)

            well_dict = _on_up_well_handler(sg, values, well_dict, graph_plate, new_graphs_list, temp_draw_tool,
                                            color_select, temp_sample_amount, temp_dilutions, temp_sample_group,
                                            replicate_loop, replicates, dose_colour_dict, concentration_count)

            # Reset or sets the sample tool counter
            if temp_sample_group >= temp_sample_amount:
                temp_sample_group = 1
            else:
                temp_sample_group += 1
            window["-SAMPLE_CHOOSER_DROPDOWN-"].update(value=temp_sample_group)

            if temp_draw_tool == "dose":
                window["-RECT_SAMPLE_TYPE-"].update(value="Dose Response")

    # deletes the rectangle used for selection
    if prior_rect:
        graph_plate.delete_figure(prior_rect)

    # reset everything
    start_point, end_point = None, None
    dragging = False
    prior_rect = None
    temp_selector = False
    temp_draw_tool = None

    return start_point, end_point, dragging, prior_rect, temp_selector, temp_draw_tool, well_dict


def on_move(config, sg, window, event, values, graph_bio_exp, well_dict_bio_info, plate_bio_info, graph_plate,
            well_dict):

    try:
        values["-BIO_INFO_CANVAS-"]
    except KeyError:
        pass
    else:
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
    try:
        values["-RECT_BIO_CANVAS-"]
    except KeyError:
        pass
    else:
        if values["-RECT_BIO_CANVAS-"][0] and values["-RECT_BIO_CANVAS-"][1]:
            try:
                graph_plate.get_figures_at_location(values['-RECT_BIO_CANVAS-'])[0]
            except (IndexError or KeyError) as error:
                # print(f"Canvas Error: {error}")
                temp_well_id = ""
                temp_well_group = ""
                temp_rep = ""
                temp_conc = ""
            else:
                temp_well = graph_plate.get_figures_at_location(values['-RECT_BIO_CANVAS-'])[0]
                try:
                    well_dict[temp_well]["well_id"]
                except (KeyError or UnboundLocalError):
                    temp_well_id = ""
                    temp_well_group = ""
                else:
                    temp_well_id = well_dict[temp_well]["well_id"]
                    temp_well_group = well_dict[temp_well]["group"]

                try:
                    well_dict[temp_well]["replicate"]
                except KeyError:
                    temp_rep = ""
                    temp_conc = ""
                else:
                    temp_rep = well_dict[temp_well]["replicate"]
                    temp_conc = well_dict[temp_well]["concentration"]

            window["-CANVAS_INFO_WELL-"].update(value=f"Well: {temp_well_id}")
            window["-CANVAS_INFO_GROUP-"].update(value=f"Group: {temp_well_group}")
            window["-CANVAS_INFO_CONC-"].update(value=f"conc: {temp_conc}")
            window["-CANVAS_INFO_REP-"].update(value=f"rep: {temp_rep}")


def draw_layout(config, sg, window, event, values, well_dict, graph_plate, archive_plates_dict):
    well_dict.clear()
    # sets the size of the well for when it draws the plate
    graph = graph_plate
    plate_type = values["-PLATE-"]
    archive_plates = values["-ARCHIVE-"]
    gui_tab = "plate_layout"
    sample_type = values["-RECT_SAMPLE_TYPE-"]

    if values["-ARCHIVE-"]:
        try:
            copy.deepcopy(archive_plates_dict[values["-ARCHIVE_PLATES-"]]["well_layout"])
        except KeyError:
            window["-ARCHIVE-"].update(False)
            values["-ARCHIVE-"] = False
        else:
            well_dict = copy.deepcopy(archive_plates_dict[values["-ARCHIVE_PLATES-"]]["well_layout"])
            well_dict = plate_layout_re_formate(config, well_dict)
            window["-PLATE-"].update(archive_plates_dict[values["-ARCHIVE_PLATES-"]]["plate_type"])
            plate_type = archive_plates_dict[values["-ARCHIVE_PLATES-"]]["plate_type"]
            sample_type = archive_plates_dict[values["-ARCHIVE_PLATES-"]]["sample_type"]

    well_dict, min_x, min_y, max_x, max_y = draw_plate(config, graph, plate_type, well_dict, gui_tab,
                                                       archive_plates, sample_layout=sample_type)
    plate_active = True

    return well_dict, min_x, min_y, max_x, max_y, plate_active, graph, plate_type, archive_plates, gui_tab, sample_type


def export_layout(config, sg, window, event, values, well_dict):
    if not well_dict:
        sg.PopupError("Please create a layout to Export")
    name = sg.PopupGetText("Name the file")
    if name:
        folder = sg.PopupGetFolder("Choose save location")
        if folder:
            plate_layout_to_excel(config, well_dict, name, folder)

            sg.Popup("Done")


def save_layout(config, sg, window, event, values, well_dict, archive_plates_dict):
    if not well_dict:
        sg.PopupError("Please create a layout to save")
    elif any("paint" in stuff.values() for stuff in well_dict.values()):
        sg.PopupError("Can't save layout with paint as well states")
    else:
        sample_type_check = sg.PopupYesNo(f"You are about to save the layout with the style: \n "
                                          f"{values['-RECT_SAMPLE_TYPE-']} \n"
                                          f"Is that fine?")

        if sample_type_check.casefold() == "yes":
            # ToDo For Dose Response, add som guard functions for checking samples and so on
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
                    "plate_layout": f"{temp_well_dict}",
                    "style": values["-RECT_SAMPLE_TYPE-"]
                }
                update_database(config, "plate_layout_sub", temp_sub_plate_layout_data)

                # Updates the plate_list and archive_plates_dict with the new plate
                plate_list, archive_plates_dict = get_plate_layout(config)

                # Updates the window with new values
                window["-ARCHIVE_PLATES-"].update(values=sorted(plate_list), value=plate_list[0])

                # window["-BIO_PLATE_LAYOUT-"].update(values=sorted(plate_list), value="")
                # window["-SEARCH_PLATE_LAYOUT-"].update(values=sorted(plate_list), value="")


def delete_layout(config, sg, window, event, values):
    # ToDO FIX THIS so it can delete layouts
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
        window["-ARCHIVE_PLATES-"].update(values=sorted(plate_list), value=plate_list[0])
        # window["-BIO_PLATE_LAYOUT-"].update(values=sorted(plate_list), value="")


def rename_layout(config, sg, window, event, values):
    # ToDO FIX THIS so it can Rename layouts - It needs to rename main and all sub layouts -
    #  Needs to take into account if layouts have been used by assays or other.
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
            try:
                window["-BIO_PLATE_LAYOUT-"]
            except KeyError:
                pass
            else:
                window["-BIO_PLATE_LAYOUT-"].update(values=sorted(plate_list), value="")


def _update_dose_tool(window, event, values, temp_sample_amount):

    try:
        int(temp_sample_amount)
    except ValueError:
        temp_sample_amount = 0

    sample_group_list = [number for number in range(1, temp_sample_amount + 1)]
    window["-SAMPLE_CHOOSER_DROPDOWN-"].update(values=sample_group_list, value=sample_group_list[0])

    heatmap = Heatmap()
    colour_list = [values["-DOSE_COLOUR_LOW-"], values["-DOSE_COLOUR_HIGH-"]]

    dose_colour_dict = heatmap.poly_linear_gradient(colour_list, temp_sample_amount)
    for colour_index, colour in enumerate(dose_colour_dict["hex"]):
        if colour_index % 2 == 0:
            dose_colour_dict["hex"][colour_index] = heatmap.get_complementary(colour)

    return dose_colour_dict


def _int_guard(window, value, fall_back_value):
    if value:
        try:
            int(value)
        except ValueError:
            new_value = fall_back_value
            window[value].update(value=fall_back_value)
        else:
            new_value = int(value)

        return new_value
    else:
        return None


def dose_sample_amount(window, event, values, total_sample_spots, dose_colour_dict):
    temp_replicates = _int_guard(window, values["-DOSE_REPLICATES-"], 1)
    temp_sample_amount = _int_guard(window, values["-DOSE_SAMPLE_AMOUNT-"], 0)

    try:
        temp_dilutions = int(floor(total_sample_spots / (temp_sample_amount * temp_replicates)))
    except ZeroDivisionError:
        return dose_colour_dict
    else:
        window["-DOSE_DILUTIONS-"].update(value=temp_dilutions)
        empty_sample_spots = total_sample_spots - (temp_dilutions * temp_sample_amount * temp_replicates)
        window["-DOSE_EMPTY_SAMPLE_SPOTS-"].update(value=empty_sample_spots)

        dose_colour_dict = _update_dose_tool(window, event, values, temp_sample_amount)
    return dose_colour_dict


def dose_dilution_replicates(window, event, values, total_sample_spots, dose_colour_dict):
    temp_replicates = _int_guard(window, values["-DOSE_REPLICATES-"], 1)
    temp_dilutions = _int_guard(window, values["-DOSE_DILUTIONS-"], 0)

    try:
        int(floor(total_sample_spots / (temp_dilutions * temp_replicates)))
    except (ZeroDivisionError or TypeError):
        return dose_colour_dict
    else:
        temp_sample_amount = int(floor(total_sample_spots / (temp_dilutions * temp_replicates)))
        window["-DOSE_SAMPLE_AMOUNT-"].update(value=temp_sample_amount)
        empty_sample_spots = total_sample_spots - (temp_dilutions * temp_sample_amount * temp_replicates)
        window["-DOSE_EMPTY_SAMPLE_SPOTS-"].update(value=empty_sample_spots)

        dose_colour_dict = _update_dose_tool(window, event, values, temp_sample_amount)

        return dose_colour_dict


def dose_colouring(window, event, values, temp_sample_amount, button, value):
    window[button].update(button_color=values[value])
    dose_colour_dict = _update_dose_tool(window, event, values, temp_sample_amount)
    return dose_colour_dict

