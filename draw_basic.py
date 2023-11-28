from heatmap import Heatmap


def draw_plate(config, graph, plate_type, well_data_dict, gui_tab, archive_plate=False, sample_layout=None, mapping=None
               , state_dict=None, skipped_well=None):
    """
    Draws different plate type in on a canvas/graph.

    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param graph: The canvas/graph that is setup in sg.window
    :type graph: PySimpleGUI.PySimpleGUI.Graph
    :param plate_type: what platetype that it needs to draw. there are 3 options 96, 384 or 1536.
    :type plate_type: str
    :param well_dict: A dict over wells, this is used for drawing saved plate layouts. The dict hold information for
        what type each well is (sample, blank or paint) and the colour of the well to draw, or the value of the sample
        if the dict is from experimental data.
    :type well_dict: dict
    :param archive_plate: bool to see if it needs to draw a blank plate or a saved plate
    :type archive_plate: bool
    :param gui_tab: what tab the plate is drawn on. differet tabs differe sizes:
    :type gui_tab: str
    :param sample_layout: This is for single point, or multiple samples with same ID.
    :type sample_layout: str
    :param mapping: Information to colour wells in specific colours, depending on what state mapping is used.
        There are 3 states - state Mapping, heatmap and hit mapping
    :type mapping: dict
    :param state_dict: A dict over what state each sample is in
    :type state_dict: dict
    :param skipped_well: A list of wells that should be skipped
    :type skipped_well: list
    :return:
        - well_dict: a dict over the wells, name, state, colour and number.
        - min_x: coordinate boundaries for the plate on the canvas
        - min_y: coordinate boundaries for the plate on the canvas
        - max_x: coordinate boundaries for the plate on the canvas
        - max_y: coordinate boundaries for the plate on the canvas
    :rtype:
        - dict
        - int
        - int
        - int
        - int
    """
    if mapping and mapping["mapping"] == "Heatmap":

        heatmap = Heatmap()

        # heatmap_dict = heatmap.dict_convert(well_data_dict, state_dict, mapping["states"])

        colour_dict, well_percentile, max_values, min_values = heatmap.heatmap_colours(well_data_dict["value"],
                                                                                       mapping["percentile"],
                                                                                       mapping["colours"])

    well_dict = {}
    graph.erase()
    if gui_tab == "bio":
        well_size = 20
        start_x = 5
        start_y = 165
    elif gui_tab == "bio_exp":
        well_size = 20
        start_x = 5
        start_y = 165
    elif gui_tab == "bio_approval_popup":
        well_size = 20
        start_x = 5
        start_y = 165
    else:
        well_size = 40
        start_x = 10
        start_y = 335

    fill_colour = config["plate_colouring"]["empty"]
    line_colour = "black"
    well_state = "empty"
    size = {"plate_96": well_size, "plate_384": well_size / 2, "plate_1536": well_size / 4}
    start_x = start_x + size[plate_type]
    rows = {"plate_96": 12, "plate_384": 24, "plate_1536": 48}
    columns = {"plate_96": 8, "plate_384": 16, "plate_1536": 32}
    well_id_col = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T",
                   "U", "V", "W", "X", "Y", "Z", "AA", "AB", "AC", "AD", "AE", "AF"]
    sample_layout_dict = {
        "single point": 1,
        "duplicate": 2,
        "triplicate": 3
    }
    counter = 0
    sample_counter = 0
    all_wells_draw = {}
    for row in range(rows[plate_type]):
        for column in range(columns[plate_type]):
            bottom_left = (start_x + row * size[plate_type],
                           start_y - column * size[plate_type])
            top_right = (bottom_left[0] - size[plate_type],
                         bottom_left[1] - size[plate_type])
            well_id = f"{well_id_col[column]}{row+1}"

            if skipped_well and well_id in skipped_well:
                fill_colour = "#FFFFFF"
                group = 0
                well_state = well_data_dict[well_id]["state"]
            else:
                if archive_plate:
                    counter += 1
                    well_state = well_data_dict[well_id]["state"]
                    if sample_layout == "single point":
                        fill_colour = config["plate_colouring"][well_state]
                    else:
                        fill_colour = well_data_dict[well_id]["colour"]
                    group = well_data_dict[well_id]["group"]

                elif sample_layout and sample_layout.casefold() != "single point":
                    group = 1
                    if well_state == "sample":
                        sample_counter += 1
                        temp_colour = group % 200
                        if group % 2 == 0:
                            fill_colour = f"#FFFF{format(temp_colour, '02x')}"
                        else:
                            fill_colour = f"#FF{format(temp_colour, '02x')}FF"
                        if sample_counter % sample_layout_dict[sample_layout] == 0:
                            group += 1
                else:
                    group = 0

                if mapping:
                    try:
                        well_data_dict["value"][well_id]
                    except KeyError:
                        map_well = False
                    else:
                        map_well = True
                    if mapping["mapping"] == "Heatmap":
                        if map_well:
                            try:
                                fill_colour = heatmap.get_well_colour(colour_dict, well_percentile, well_id)
                            except ZeroDivisionError:
                                fill_colour = "#FFFFFF"
                        else:
                            fill_colour = "#FFFFFF"
                    elif mapping["mapping"] == "Hit Map":
                        if map_well:
                            for _, params in mapping["bins"].items():

                                if params["use"] and params["min"] <= well_data_dict["value"][well_id] <= params["max"]:
                                    fill_colour = params["colour"]
                                else:
                                    fill_colour = "black"
                        else:
                            fill_colour = "#FFFFFF"
            all_wells_draw[well_id] = (bottom_left, top_right)
            temp_well = graph.DrawRectangle(bottom_left, top_right, line_color=line_colour, fill_color=fill_colour)
            if column == 0 and row == 0:
                off_set = temp_well - 1
            temp_well = temp_well - off_set
            well_dict[temp_well] = {}
            if group == "dose":
                group = well_data_dict[well_id]["group"]
            well_dict[temp_well]["group"] = group
            well_dict[temp_well]["well_id"] = well_id
            well_dict[temp_well]["state"] = well_state
            well_dict[temp_well]["colour"] = fill_colour
    min_x = start_x - size[plate_type]
    min_y = start_y - (columns[plate_type] * size[plate_type])
    max_x = start_x + (rows[plate_type] * size[plate_type])-size[plate_type]
    max_y = start_y

    return well_dict, min_x, min_y, max_x, max_y, off_set