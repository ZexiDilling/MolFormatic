import PySimpleGUI as sg
from copy import deepcopy, copy

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from info import matrix_header, plate_384_row, plate_96_row
from visualization import Toolbar
# from gui_functions import sub_settings_matrix
# from database_controller import FetchData
# from excel_handler import purity_sample_layout_import, purity_sample_layout_export
# from gui_functions import sort_table

# todo FIX THIS WHOLE THING!!!!!!!! - easy fix is to import the module that each function needs to run!!!
from plate_formatting import plate_layout_re_formate


def _matrix_popup_layout(calc, state=None, method=None):
    table_headings = matrix_header
    table_data = []
    col_matrix_tabl = sg.Column([
        [sg.Table(table_data, headings=table_headings, auto_size_columns=False,
                  key="-POPUP_MATRIX_TABLE-", vertical_scroll_only=False)],
        [sg.T("", size=10)],
    ])

    row_matrix = sg.Frame("Matrix", [
        [col_matrix_tabl],
        [sg.T("", key="-POPUP_TABLE_INFO-")],
        [sg.T("Method", size=14), sg.T("State", size=14), sg.T("Calc", size=14)],
        [sg.DropDown([], key="-POPUP_MATRIX_METHOD-", size=14, default_value=calc),
         sg.DropDown([], key="-POPUP_MATRIX_STATE-", size=14, default_value=state),
         sg.DropDown([], key="-POPUP_MATRIX_CALC-", size=14, default_value=method)],
        [sg.B("Generate Matrix", key="-POPUP_MATRIX_GENERATOR_BUTTON-"),
         sg.B("close", key="-CLOSE-")]
    ], size=(1500, 500))

    layout = [[row_matrix]]

    return sg.Window("Matrix", layout, finalize=True, resizable=True)


def matrix_popup(data_dict, calc_values, state_values, method_values, calc, state=None, method=None):
    window = _matrix_popup_layout(calc, state, method)
    window["-POPUP_MATRIX_METHOD-"].update(values=method_values)
    window["-POPUP_MATRIX_STATE-"].update(values=state_values)
    window["-POPUP_MATRIX_CALC-"].update(values=calc_values)

    if calc:
        table_data, display_columns = sub_settings_matrix(data_dict, calc, method, state)
        window["-POPUP_MATRIX_TABLE-"].update(values=table_data)
    else:
        calc = "z_prime"
        table_data, display_columns = sub_settings_matrix(data_dict, calc, method, state)
        window["-POPUP_MATRIX_TABLE-"].update(values=table_data)
    if calc == "z_prime":
        temp_text = f"Table is showing matrix for {calc}"
    else:
        temp_text = f"Table is showing matrix of {calc} for {method}-data, {state}-wells"
    window["-POPUP_TABLE_INFO-"].update(value=temp_text)
    window.Element("-POPUP_MATRIX_TABLE-").Widget.configure(displaycolumns=display_columns)

    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == "-CLOSE-":
            break
        if event == "-POPUP_MATRIX_GENERATOR_BUTTON-":
            if not values["-POPUP_MATRIX_CALC-"]:
                sg.PopupError("Missing Calculation choice ")
            elif values["-POPUP_MATRIX_CALC-"] != "z_prime" and not values["-POPUP_MATRIX_STATE-"]:
                sg.PopupError("Missing State selection")
            elif values["-POPUP_MATRIX_CALC-"] != "z_prime" and not values["-POPUP_MATRIX_METHOD-"]:
                sg.PopupError("Missing Method selection")
            else:

                calc = values["-POPUP_MATRIX_CALC-"]

                if calc == "z_prime":
                    temp_text = f"Table is showing matrix for {calc}"
                    state = None
                    method = None
                else:
                    method = values["-POPUP_MATRIX_METHOD-"]
                    state = values["-POPUP_MATRIX_STATE-"]
                    temp_text = f"Table is showing matrix of {calc} for {method}-data and well state: {state}"

                table_data, display_columns = sub_settings_matrix(data_dict, calc, method, state)
                window["-POPUP_MATRIX_TABLE-"].update(values=table_data)

                window["-POPUP_TABLE_INFO-"].update(value=temp_text)
                window.Element("-POPUP_MATRIX_TABLE-").Widget.configure(displaycolumns=display_columns)


def _sample_checker_layout(table_data, headings, db_data):
    raw_table_col = sg.Frame("Compound Samples", [[
        sg.Column([
            [sg.Table(values=table_data, headings=headings,
                      key="-POP_SAMPLE_CHECKER_TABLE-", enable_events=True, enable_click_events=True)]
        ])
    ]])

    layout = [
        [raw_table_col],
        [sg.Button("Done", key="-POP_SAMPLE_CHECKER_OK-", expand_x=True),
         sg.Button("Check Database", key="-POP_SAMPLE_CHECKER_DB_CHECK-", expand_x=True, visible=db_data),
         sg.Button("Colour", key="-POP_SAMPLE_CHECKER_COLOUR-", expand_x=True, visible=db_data),
         sg.Button("Import", key="-POP_SAMPLE_CHECKER_IMPORT-", expand_x=True, visible=db_data),
         sg.Button("Export", key="-POP_SAMPLE_CHECKER_EXPORT-", expand_x=True, visible=db_data),
         sg.Button("Rename To Raw", key="-POP_SAMPLE_CHECKER_RE_NAME_TO_RAW-", expand_x=True, visible=not db_data),
         sg.Button("Rename To Excel", key="-POP_SAMPLE_CHECKER_RE_NAME_TO_EXCEL-", expand_x=True, visible=not db_data),
         sg.Button("Back To Default", key="-POP_SAMPLE_CHECKER_BACK_TO_DEFAULT-", expand_x=True, visible=not db_data),
         sg.Button("Cancel", key="-WINDOW_TWO_CANCEL-", expand_x=True)]
    ]
    return sg.Window("Samples", layout, finalize=True, resizable=True), table_data


def sample_to_compound_name_controller(config, data_dict, db_data=True):
    fd = FetchData(config)
    samples = [keys for keys in data_dict]
    found_dict = fd.list_to_rows(samples)
    found_samples = []
    not_found_samples = []
    blanks = []
    for sample in found_dict:
        found_samples.append(sample)
    for sample in samples:
        if sample.lower() == "blank":
            blanks.append(samples)
        elif sample not in found_samples:
            not_found_samples.append(sample)

    if not not_found_samples:
        return False

    default_colouring = []
    search_reverse = {}
    export_file = None
    table_headings = ["Raw Name", "compound ID", "In DB"]

    table_data = []
    for samples in found_samples:
        temp_data = [samples, samples, "Found"]
        table_data.append(temp_data)

    for samples in not_found_samples:
        temp_data = [samples, "None", "Not in DB"]
        table_data.append(temp_data)

    for samples in blanks:
        temp_data = [samples, "Blank", "Blank"]
        table_data.append(temp_data)

    window, table_data = _sample_checker_layout(table_data, table_headings, db_data)
    window["-POP_SAMPLE_CHECKER_TABLE-"].bind('<Double-Button-1>', "+-double click-")

    while True:
        event, values = window.read()

        if event == sg.WIN_CLOSED or event == "-WINDOW_TWO_CANCEL-":
            break

        if event == "-POP_SAMPLE_CHECKER_TABLE-+-double click-":
            try:
                table_row = values["-POP_SAMPLE_CHECKER_TABLE-"][0]
            except IndexError:
                pass
            else:
                new_name = sg.PopupGetText("New name for the sample")
                if new_name:
                    table_data[table_row][1] = new_name
                    if fd.list_to_rows([new_name]):
                        table_data[table_row][2] = "Found"
                    else:
                        table_data[table_row][2] = "Not in DB"
                window["-POP_SAMPLE_CHECKER_TABLE-"].update(values=table_data, row_colors=default_colouring)

        if event == "-POP_SAMPLE_CHECKER_COLOUR-":
            table_colouring = []

            for row_index, row_data in enumerate(table_data):
                temp_colour = [row_index]
                default_colouring.append([row_index, None])
                if row_data[2] == "Found":
                    temp_colour.append("green")
                elif row_data[2] == "Blank":
                    temp_colour.append("yellow")
                else:
                    temp_colour.append("red")
                table_colouring.append(temp_colour)

            window['-POP_SAMPLE_CHECKER_TABLE-'].update(row_colors=table_colouring)

        if event == "-POP_SAMPLE_CHECKER_DB_CHECK-":
            for data in table_data:
                if data[1] == "Blank":
                    data[2] = "Blank"
                elif fd.list_to_rows([data[1]]):
                    data[2] = "Found"
                else:
                    data[2] = "Not in DB"

            window["-POP_SAMPLE_CHECKER_TABLE-"].update(values=table_data, row_colors=default_colouring)

        if event == "-POP_SAMPLE_CHECKER_IMPORT-":
            if export_file:
                import_file = export_file
            else:
                import_file = sg.popup_get_file("Please select file")

            table_data = purity_sample_layout_import(import_file, table_headings)
            window["-POP_SAMPLE_CHECKER_TABLE-"].update(values=table_data, row_colors=default_colouring)
            export_file = None

        if event == "-POP_SAMPLE_CHECKER_EXPORT-":
            export_file = purity_sample_layout_export(config, table_data, table_headings)

        if event == "-POP_SAMPLE_CHECKER_OK-":

            if any("Not in DB" in data for data in table_data):
                pop_sample_export = sg.PopupYesNo("No all sample names are in the Database. "
                                                  "Any sample not in DB or blank will be deletede. Continue?")
            else:
                pop_sample_export = "Yes"

            if pop_sample_export == "Yes":
                new_sample_names = {}
                for data in table_data:
                    if data[2] != "Found" and data[1].LOWER() != "blank":
                        new_sample_names[data[0]] = "Delete"
                    else:
                        if data[0] != data[1]:
                            new_sample_names[data[0]] = data[1]
                window.Close()
                return new_sample_names

        if isinstance(event, tuple) and event[0] == "-POP_SAMPLE_CHECKER_TABLE-":
            # TABLE CLICKED Event has value in format ('-TABLE=', '+CLICKED+', (row,col))
            # You can also call Table.get_last_clicked_position to get the cell clicked
            try:
                search_reverse[event[0]]
            except KeyError:
                search_reverse[event[0]] = {}
            if event[2][0] == -1 and event[2][1] != -1:  # Header was clicked and wasn't the "row" column
                col_num_clicked = event[2][1]

                try:
                    search_reverse[event[0]][col_num_clicked]
                except KeyError:
                    search_reverse[event[0]][col_num_clicked] = False
                new_table, search_reverse[event[0]][col_num_clicked] = \
                    sort_table(table_data[0:][:], (col_num_clicked, 0), search_reverse[event[0]][col_num_clicked])

                window["-POP_SAMPLE_CHECKER_TABLE-"].update(new_table)
                # all_table_data[clicked_table] = [all_table_data[clicked_table][0]] + new_table
                table_data = new_table


def ms_raw_name_guard(raw_data_samples, excel_data_samples, db_data, config):
    table_data = []
    new_name = {}
    for raw_sample_index, raw_sample in enumerate(raw_data_samples):
        if raw_sample in excel_data_samples:
            temp_data = [raw_sample, raw_sample, raw_sample]
        else:
            temp_data = [raw_sample, excel_data_samples[raw_sample_index], raw_sample]
        table_data.append(temp_data)

    default_table_data = deepcopy(table_data)

    default_colouring = []
    search_reverse = {}
    export_file = None
    table_headings = ["Raw Name", "Excel name", "Final Name"]

    window, table_data = _sample_checker_layout(table_data, table_headings, db_data)
    window["-POP_SAMPLE_CHECKER_TABLE-"].bind('<Double-Button-1>', "+-double click-")

    while True:
        event, values = window.read()

        if event == sg.WIN_CLOSED or event == "-WINDOW_TWO_CANCEL-" or event == "-POP_SAMPLE_CHECKER_OK-":
            for row_index, row in enumerate(table_data):
                new_name[table_data[row_index][2]] = {"raw": table_data[row_index][0],
                                                      "excel": table_data[row_index][1]}

            return new_name
            window.close()

        if event == "-POP_SAMPLE_CHECKER_RE_NAME_TO_RAW-":
            for row_index, row in enumerate(table_data):
                table_data[row_index][2] = table_data[row_index][0]

            window["-POP_SAMPLE_CHECKER_TABLE-"].update(values=table_data, row_colors=default_colouring)

        if event == "-POP_SAMPLE_CHECKER_RE_NAME_TO_EXCEL-":
            for row_index, row in enumerate(table_data):
                table_data[row_index][2] = table_data[row_index][1]

            window["-POP_SAMPLE_CHECKER_TABLE-"].update(values=table_data, row_colors=default_colouring)

        if event == "-POP_SAMPLE_CHECKER_BACK_TO_DEFAULT-":
            table_data = deepcopy(default_table_data)

            window["-POP_SAMPLE_CHECKER_TABLE-"].update(values=table_data, row_colors=default_colouring)

        if event == "-POP_SAMPLE_CHECKER_TABLE-+-double click-":
            try:
                table_row = values["-POP_SAMPLE_CHECKER_TABLE-"][0]
            except IndexError:
                pass
            else:
                new_name = sg.PopupGetText("New name for the sample")
                if new_name:
                    table_data[table_row][2] = new_name

                window["-POP_SAMPLE_CHECKER_TABLE-"].update(values=table_data, row_colors=default_colouring)

        if isinstance(event, tuple) and event[0] == "-POP_SAMPLE_CHECKER_TABLE-":
            # TABLE CLICKED Event has value in format ('-TABLE=', '+CLICKED+', (row,col))
            # You can also call Table.get_last_clicked_position to get the cell clicked
            try:
                search_reverse[event[0]]
            except KeyError:
                search_reverse[event[0]] = {}
            if event[2][0] == -1 and event[2][1] != -1:  # Header was clicked and wasn't the "row" column
                col_num_clicked = event[2][1]

                try:
                    search_reverse[event[0]][col_num_clicked]
                except KeyError:
                    search_reverse[event[0]][col_num_clicked] = False
                new_table, search_reverse[event[0]][col_num_clicked] = \
                    sort_table(table_data[0:][:], (col_num_clicked, 0), search_reverse[event[0]][col_num_clicked])

                window["-POP_SAMPLE_CHECKER_TABLE-"].update(new_table)
                # all_table_data[clicked_table] = [all_table_data[clicked_table][0]] + new_table
                table_data = new_table


def _new_headlines_layout(table_data, headings):
    raw_table_col = sg.Frame("Please select new headlines", [[
        sg.Column([
            [sg.Table(values=table_data, headings=headings,
                      key="-POP_HEADLINE_TABLE-", enable_events=True, enable_click_events=True,
                      tooltip='double click to select a new headline in the "NEW headline" column ')]
        ])
    ]])

    layout = [
        [raw_table_col],
        [sg.Button("Done", key="-POP_SAMPLE_CHECKER_OK-", expand_x=True),
         sg.Button("Cancel", key="-WINDOW_TWO_CANCEL-", expand_x=True)]
    ]
    return sg.Window("Samples", layout, finalize=True, resizable=True), table_data


def new_headlines_popup(right_headlines, wrong_headlines):
    search_reverse = {}
    table_data = []
    new_headlines = {}

    for headline_index, headline in enumerate(wrong_headlines):
        temp_data = [headline, wrong_headlines[headline_index]]
        table_data.append(temp_data)

    table_headings = ["Headline in file", "New Headline"]

    window, table_data = _new_headlines_layout(table_data, table_headings)
    window["-POP_HEADLINE_TABLE-"].bind('<Double-Button-1>', "+-double click-")

    while True:
        event, values = window.read()

        if event == sg.WIN_CLOSED or event == "-WINDOW_TWO_CANCEL-" or event == "-POP_SAMPLE_CHECKER_OK-":
            for row_index, row in enumerate(table_data):
                temp_wrong_headline = table_data[row_index][0]
                if type(table_data[row_index][1]) == list:
                    temp_right_headline = table_data[row_index][1][0]
                else:
                    temp_right_headline = table_data[row_index][1]

                new_headlines[temp_wrong_headline] = temp_right_headline

            # check if there are duplicates
            flag = False
            hash_val = dict()
            for keys in new_headlines:
                if new_headlines[keys] in hash_val:
                    flag = True
                    break
                else:
                    hash_val[new_headlines[keys]] = 1
            if flag:
                sg.popup_error("There are repeated values, please try again")
            else:
                window.close()
                return new_headlines

        if event == "-POP_HEADLINE_TABLE-+-double click-":
            try:
                table_row = values["-POP_HEADLINE_TABLE-"][0]
            except IndexError:
                pass
            else:
                event, values = sg.Window('Select One',
                                          layout=[[sg.Listbox(right_headlines, key='_LIST_',
                                                              size=(max([len(str(v)) for v in right_headlines]) + 2,
                                                                    len(right_headlines)), select_mode='extended',
                                                              bind_return_key=True), sg.OK()]]).read(close=True)

                temp_new_headlines = values['_LIST_'] if event is not None else None
                if temp_new_headlines:
                    table_data[table_row][1] = temp_new_headlines

                window["-POP_HEADLINE_TABLE-"].update(values=table_data)

        if isinstance(event, tuple) and event[0] == "-POP_HEADLINE_TABLE-":
            # TABLE CLICKED Event has value in format ('-TABLE=', '+CLICKED+', (row,col))
            # You can also call Table.get_last_clicked_position to get the cell clicked
            try:
                search_reverse[event[0]]
            except KeyError:
                search_reverse[event[0]] = {}
            if event[2][0] == -1 and event[2][1] != -1:  # Header was clicked and wasn't the "row" column
                col_num_clicked = event[2][1]

                try:
                    search_reverse[event[0]][col_num_clicked]
                except KeyError:
                    search_reverse[event[0]][col_num_clicked] = False
                new_table, search_reverse[event[0]][col_num_clicked] = \
                    sort_table(table_data[0:][:], (col_num_clicked, 0), search_reverse[event[0]][col_num_clicked])

                window["-POP_HEADLINE_TABLE-"].update(new_table)
                # all_table_data[clicked_table] = [all_table_data[clicked_table][0]] + new_table
                table_data = new_table


def _plate_layout_chooser_layout(table_data, headings):
    raw_table_col = sg.Frame("Please select plate layout for each plate", [[
        sg.Column([
            [sg.Table(values=table_data, headings=headings,
                      key="-POP_HEADLINE_TABLE-", enable_events=True, enable_click_events=True,
                      tooltip='double click to select a new headline in the "NEW headline" column ')]
        ])
    ]])

    layout = [
        [raw_table_col],
        [sg.Button("Done", key="-POP_SAMPLE_CHECKER_OK-", expand_x=True),
         sg.Button("Cancel", key="-WINDOW_TWO_CANCEL-", expand_x=True)]
    ]
    return sg.Window("Samples", layout, finalize=True, resizable=True), table_data


def plate_layout_chooser(files, default_plate_layout, all_plate_layouts):
    table_data = []
    plate_to_layouy = {}
    # Add the options to skip files.
    all_plate_layouts.append("skip")

    for file in files:
        temp_data = [file, default_plate_layout]
        table_data.append(temp_data)

    table_headings = ["Plate", "Layout"]

    window, table_data = _plate_layout_chooser_layout(table_data, table_headings)
    window["-POP_HEADLINE_TABLE-"].bind('<Double-Button-1>', "+-double click-")

    while True:
        event, values = window.read()

        if event == sg.WIN_CLOSED or event == "-WINDOW_TWO_CANCEL-" or event == "-POP_SAMPLE_CHECKER_OK-":
            # Grabs all the data from the table
            for row_index, row in enumerate(table_data):
                if type(row[1]) == list:
                    plate_to_layouy[row[0]] = row[1][0]
                else:
                    plate_to_layouy[row[0]] = row[1]

            window.close()
            return plate_to_layouy

        if event == "-POP_HEADLINE_TABLE-+-double click-":
            try:
                table_row = values["-POP_HEADLINE_TABLE-"][0]
            except IndexError:
                pass
            else:
                event, values = sg.Window('Select One',
                                          layout=[[sg.Listbox(all_plate_layouts, key='_LIST_',
                                                              size=(max([len(str(v)) for v in all_plate_layouts]) + 2,
                                                                    len(all_plate_layouts)), select_mode='extended',
                                                              bind_return_key=True), sg.OK()]]).read(close=True)

                temp_new_plate_layout = values['_LIST_'] if event is not None else None
                if temp_new_plate_layout:
                    table_data[table_row][1] = temp_new_plate_layout

                window["-POP_HEADLINE_TABLE-"].update(values=table_data)


def _bio_data_approval_table_layout(config, plate_table_data, plate_headings, compound_table_data, compound_headings,
                                    plate_analyse_methods, plate_calculations, analyse_methods, draw_options,
                                    well_state_overview):
    """
    The layout for the Bio Approval table that comes up when you import bio data to the database
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param plate_table_data: The data for the plate table
    :type plate_table_data: list
    :param plate_headings: the headings for the plate table
    :type plate_headings: list
    :param compound_table_data: The data for the compound table
    :type compound_table_data: list
    :param compound_headings: the headings for the compound table
    :type compound_headings: list
    :param plate_analyse_methods: The way the plate have been calculated
    :type plate_analyse_methods: list
    :param plate_calculations: The different calculations made for each plate
    :type plate_calculations: list
    :param analyse_methods: Different ways to analyse a plate, like historgrams or gragps...
    :type analyse_methods: list
    :param well_state_overview: The states included in the plate-map. It is used to disabling states that are not used
    :type well_state_overview: dict
    :return: The layout for the popup
    """
    standard_size = 12

    raw_table_col = sg.Frame("Please approve or dimiss plates", [[
        sg.Column([
            [sg.Table(values=plate_table_data, headings=plate_headings,
                      key="-BIO_APPROVAL_PLATE_TABLE-", enable_events=True, enable_click_events=True)],
            [sg.Table(values=compound_table_data, headings=compound_headings,
                      key="-BIO_APPROVAL_COMPOUND_TABLE-", enable_events=True, enable_click_events=True)]
        ])
    ]])

    graph_window = sg.Frame("Plate Layout", [[
        sg.Graph(canvas_size=(500, 300), graph_bottom_left=(0, 0), graph_top_right=(250, 175),
                 background_color='grey', key="-BIO_APPROVAL_TABL_GRAPH-", enable_events=False,
                 drag_submits=False, motion_events=True)],
        sg.vtop([
            sg.TabGroup([[
                sg.Tab("State", [[
                    sg.Column([
                        [sg.T("States", size=10), sg.T("Colour", size=10), sg.T("Include")],
                        [sg.HorizontalSeparator()],
                        [sg.Text("Sample:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["sample"], size=10,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_SAMPLE-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABL_SAMPLE-", enable_events=True,
                                     disabled=not well_state_overview["sample"]
                                     , default=True)],
                        [sg.Text("Blank:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["blank"], size=10,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_BLANK-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABL_BLANK-", enable_events=True,
                                     disabled=not well_state_overview["blank"])],
                        [sg.Text("Maximum:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["max"], size=10,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_NAX-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABL_MAX-", enable_events=True,
                                     disabled=not well_state_overview["max"])],
                        [sg.Text("Minimum:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["minimum"], size=10,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_MINIMUM-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABL_MIN-", enable_events=True,
                                     disabled=not well_state_overview["minimum"])],
                        [sg.Text("Positive Control:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["positive"], size=10,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_POSITIVE-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABL_POSITIVE-", enable_events=True,
                                     disabled=not well_state_overview["pos"])],
                        [sg.Text("Negative Control:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["negative"], size=10,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_NEGATIVE-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABL_NEGATIVE-", enable_events=True,
                                     disabled=not well_state_overview["neg"])],
                        [sg.Text("Empty:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["empty"], size=10,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_EMPTY-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABL_EMPTY-", disabled=not well_state_overview["empty"])]
                    ])
                ]]),
                sg.Tab("HeatMap", [[
                    sg.Column([
                        [sg.T("States", size=10), sg.T("Colour", size=6)],
                        [sg.HorizontalSeparator()],
                        [sg.Text("Min:", size=standard_size),
                         sg.T(background_color=config["Settings_bio"]["plate_report_heatmap_colours_low"], size=10,
                              key="-BIO_PLATE_LAYOUT_COLOUR_HEATMAP_MIN-", relief="groove")],
                        [sg.Text("Mid:", size=standard_size),
                         sg.T(background_color=config["Settings_bio"]["plate_report_heatmap_colours_mid"], size=10,
                              key="-BIO_PLATE_LAYOUT_COLOUR_HEATMAP_MID-", relief="groove")],
                        [sg.Text("Max:", size=standard_size),
                         sg.T(background_color=config["Settings_bio"]["plate_report_heatmap_colours_high"], size=10,
                              key="-BIO_PLATE_LAYOUT_COLOUR_HEATMAP_MAX-", relief="groove")]
                    ])
                ]]),
                sg.Tab("Hit Map", [[
                    sg.Column([
                        [sg.T("Bins", size=5), sg.T("Min", size=5), sg.T("Max", size=5), sg.T("Colour", size=6)],
                        [sg.HorizontalSeparator()],
                        [sg.Text("T1:", size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_1_min"], size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_1_max"], size=5),
                         sg.T(background_color=config["Settings_bio"]["plate_report_pora_threshold_colour_th_1"],
                              size=5, key="-BIO_PLATE_LAYOUT_COLOUR_HEATMAP_MIN-", relief="groove")],
                        [sg.Text("T2:", size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_2_min"], size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_2_max"], size=5),
                         sg.T(background_color=config["Settings_bio"]["plate_report_pora_threshold_colour_th_2"],
                              size=5, key="-BIO_PLATE_LAYOUT_COLOUR_HEATMAP_MIN-", relief="groove")],
                        [sg.Text("T3:", size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_3_min"], size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_3_max"], size=5),
                         sg.T(background_color=config["Settings_bio"]["plate_report_pora_threshold_colour_th_3"],
                              size=5, key="-BIO_PLATE_LAYOUT_COLOUR_HEATMAP_MIN-", relief="groove")],
                        [sg.Text("T4:", size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_4_min"], size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_4_max"], size=5),
                         sg.T(background_color=config["Settings_bio"]["plate_report_pora_threshold_colour_th_4"],
                              size=5, key="-BIO_PLATE_LAYOUT_COLOUR_HEATMAP_MIN-", relief="groove")],
                        [sg.Text("T5:", size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_5_min"], size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_5_max"], size=5),
                         sg.T(background_color=config["Settings_bio"]["plate_report_pora_threshold_colour_th_5"],
                              size=5, key="-BIO_PLATE_LAYOUT_COLOUR_HEATMAP_MIN-", relief="groove")],
                        [sg.Text("T6:", size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_6_min"], size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_6_max"], size=5),
                         sg.T(background_color=config["Settings_bio"]["plate_report_pora_threshold_colour_th_6"],
                              size=5, key="-BIO_PLATE_LAYOUT_COLOUR_HEATMAP_MIN-", relief="groove")],
                        [sg.Text("T7:", size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_7_min"], size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_7_max"], size=5),
                         sg.T(background_color=config["Settings_bio"]["plate_report_pora_threshold_colour_th_7"],
                              size=5, key="-BIO_PLATE_LAYOUT_COLOUR_HEATMAP_MIN-", relief="groove")],
                        [sg.Text("T8:", size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_8_min"], size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_8_max"], size=5),
                         sg.T(background_color=config["Settings_bio"]["plate_report_pora_threshold_colour_th_8"],
                              size=5, key="-BIO_PLATE_LAYOUT_COLOUR_HEATMAP_MIN-", relief="groove")],
                        [sg.Text("T9:", size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_9_min"], size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_9_max"], size=5),
                         sg.T(background_color=config["Settings_bio"]["plate_report_pora_threshold_colour_th_9"],
                              size=5, key="-BIO_PLATE_LAYOUT_COLOUR_HEATMAP_MIN-", relief="groove")],
                        [sg.Text("T01:", size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_10_min"], size=5),
                         sg.Text(config["Settings_bio"]["plate_report_pora_threshold_th_10_max"], size=5),
                         sg.T(background_color=config["Settings_bio"]["plate_report_pora_threshold_colour_th_10"],
                              size=5, key="-BIO_PLATE_LAYOUT_COLOUR_HEATMAP_MIN-", relief="groove")],
                    ])
                ]])
            ]], selected_background_color=config["GUI"]["tab_colour"], key="-BIO_APPROVAL_TABLE_COLOUR_TAB-",
                enable_events=False, expand_x=True, expand_y=True),
            sg.Column([
                [sg.T("", key="-BIO_APPROVAL_WELL_ID-")],
                [sg.Frame("", [
                    [sg.T("Calculations")],
                    [sg.HorizontalSeparator()],
                    [sg.T("avg", size=10), sg.T("Calc", key="-BIO_APPROVAL_TABLE_CALC_AVG-", size=10)],
                    [sg.T("stdev", size=10), sg.T("Calc", key="-BIO_APPROVAL_TABLE_CALC_STDEV-", size=10)],
                    [sg.T("pstdev", size=10), sg.T("Calc", key="-BIO_APPROVAL_TABLE_CALC_PSTDEV-", size=10)],
                    [sg.T("pvariance", size=10), sg.T("Calc", key="-BIO_APPROVAL_TABLE_CALC_PVARIANCE-", size=10)],
                    [sg.T("variance", size=10), sg.T("Calc", key="-BIO_APPROVAL_TABLE_CALC_VARIANCE-", size=10)],
                    [sg.T("st_dev_%", size=10), sg.T("Calc", key="-BIO_APPROVAL_TABLE_CALC_ST_DEV-", size=10)],
                    [sg.T("S/B", size=10), sg.T("Calc", key="-BIO_APPROVAL_TABLE_CALC_SB-", size=10)]
                ])],
                [sg.Frame("Note", [
                    [sg.Multiline(default_text="", key="-BIO_APPROVAL_PLATE_NOTE-", size=(25, 2))]
                ])]
            ])
        ]),
        [sg.T("Draw style", size=14), sg.T("Data", size=14), sg.T("Calculations")],
        [sg.DropDown(values=draw_options, key="-BIO_APPROVAL_TABLE_DROPDOWN_DRAW_OPTIONS-", size=14,
                     enable_events=True, default_value=draw_options[2]),
         sg.DropDown(values=plate_analyse_methods, key="-BIO_APPROVAL_TABLE_DROPDOWN_PLATE_ANALYSE_METHODS-",
                     size=14, enable_events=True, default_value=plate_analyse_methods[2],
                     tooltip="The different ways that he plate have been calculated"),
         sg.DropDown(values=plate_calculations, key="-BIO_APPROVAL_TABLE_DROPDOWN_CALCULATIONS-", size=14,
                     enable_events=True, default_value=plate_calculations[3],
                     tooltip="Shows different calculations for the plate")],
        [sg.Button("Re-calculate", key="-BIO_APPROVAL_TABLE_RE_CALCULATE-", size=14,
                   tooltip="Will re-calculate with the new layout. OBS"),
         sg.Button("Apply", key="-BIO_APPROVAL_TABLE_APPLY-", size=14,
                   tooltip="Will apply the new calculations the the data. "
                           "This will overwrite the old calculated data, but raw data is still saved"),
         sg.Button("Add Note", key="-BIO_APPROVAL_TABLE_ADD_NOTE-", size=14,
                   tooltip="Adds a note to the plate! Any plate not approved needs to have a note."),
         sg.Button("Add Note to all", key="-BIO_APPROVAL_TABLE_ADD_NOTE_ALL-", size=14,
                   tooltip="Adds the same note to all the plates.")
         ]

    ])

    canvas_window = sg.Frame("Data", [[
        sg.Column([
            [sg.Canvas(key="-BIO_APPROVAL_TABL_TOOLBAR-")],
            [sg.Canvas(key="-BIO_APPROVAL_TABL_CANVAS-", background_color="Grey", size=(500, 300))],
            [sg.DropDown(values=analyse_methods, key="-BIO_APPROVAL_TABLE_ANALYSE_METHODS-", enable_events=True,
                         default_value=analyse_methods[0])]
        ])
    ]])

    layout = [sg.vtop([raw_table_col, graph_window, canvas_window]),
              [sg.Button("Done", key="-BIO_DATA_APPROVED-", expand_x=True),
               sg.Button("Cancel", key="-WINDOW_TWO_CANCEL-", expand_x=True)]
              ]
    return sg.Window("Samples", layout, finalize=True, resizable=True), plate_table_data


def _generate_well_dict(config, plate_type, plate_data, analysed_method, show_state_list):
    well_dict = {"value": {}}

    well_data_dict = {}
    # Get a list of wells:
    if plate_type == "plate_384":
        plate_layout = plate_384_row
    elif plate_type == "plate_96":
        plate_layout = plate_96_row
    else:
        print("MISSING LAYOUT FOR PLATE_1536")

    # Get states
    all_states = []
    for states in plate_data[analysed_method]:
        if states != "wells" and states not in all_states:
            all_states.append(states)
    # if analysed_method == "pora":
    #     ...
    # elif analysed_method == "normalised":
    #     ...
    #
    # else:
    # This should always be the state mapping of the assay
    for wells in plate_layout:
        for state in all_states:
            if wells in plate_data[analysed_method][state]:
                temp_state = state
                colour = config["plate_colouring"][temp_state]
                if show_state_list[temp_state]:
                    well_dict["value"][wells] = plate_data[analysed_method]["wells"][wells]
                well_dict[wells] = {"group": 0,
                                    "value": plate_data[analysed_method]["wells"][wells],
                                    "state": temp_state,
                                    "colour": colour}
                continue

    return well_dict, all_states


def _draw_plate_on_graph(draw_plate, config, plate_data, analysed_method, graph, plate_size, graph_placement,
                         show_state_list, draw_option):
    well_dict, all_states = _generate_well_dict(config, plate_size, plate_data, analysed_method, show_state_list)
    if draw_option == "heatmap":
        mapping = {
            "mapping": "Heatmap",
            "colours": {"low": [config["Settings_bio"]["plate_report_heatmap_colours_low"],
                                config["Settings_bio"]["plate_report_heatmap_colours_mid"]],
                        "high": [config["Settings_bio"]["plate_report_heatmap_colours_mid"],
                                 config["Settings_bio"]["plate_report_heatmap_colours_high"]]},
            "states": all_states,
            "percentile": {"low": float(0),
                           "mid": float(50),
                           "high": float(100)}
        }
    elif draw_option == "hit_map":
        mapping = {
            "mapping": "Hit Map",
            "bins": {"th_1": {
                "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_1_use"),
                "min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_1_min"),
                "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_1_max"),
                "colour": config["Settings_bio"]["plate_report_pora_threshold_colour_th_1"]},
                "th_2": {
                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_2_use"),
                    "min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_2_min"),
                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_2_max"),
                    "colour": config["Settings_bio"]["plate_report_pora_threshold_colour_th_2"]},
                "th_3": {
                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_3_use"),
                    "min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_3_min"),
                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_3_max"),
                    "colour": config["Settings_bio"]["plate_report_pora_threshold_colour_th_3"]},
                "th_4": {
                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_4_use"),
                    "min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_4_min"),
                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_4_max"),
                    "colour": config["Settings_bio"]["plate_report_pora_threshold_colour_th_4"]},
                "th_5": {
                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_5_use"),
                    "min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_5_min"),
                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_5_max"),
                    "colour": config["Settings_bio"]["plate_report_pora_threshold_colour_th_5"]},
                "th_6": {
                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_6_use"),
                    "min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_6_min"),
                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_6_max"),
                    "colour": config["Settings_bio"]["plate_report_pora_threshold_colour_th_6"]},
                "th_7": {
                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_7_use"),
                    "min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_7_min"),
                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_7_max"),
                    "colour": config["Settings_bio"]["plate_report_pora_threshold_colour_th_7"]},
                "th_8": {
                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_8_use"),
                    "min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_8_min"),
                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_8_max"),
                    "colour": config["Settings_bio"]["plate_report_pora_threshold_colour_th_8"]},
                "th_9": {
                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_9_use"),
                    "min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_9_min"),
                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_9_max"),
                    "colour": config["Settings_bio"]["plate_report_pora_threshold_colour_th_9"]},
                "th_10": {
                    "use": config["Settings_bio"].getboolean("plate_report_pora_threshold_th_10_use"),
                    "min": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_10_min"),
                    "max": config["Settings_bio"].getfloat("plate_report_pora_threshold_th_10_max"),
                    "colour": config["Settings_bio"]["plate_report_pora_threshold_colour_th_10"]}
            }
        }
    else:
        mapping = None
    draw_plate(config, graph, plate_size, well_dict, graph_placement, archive_plate=True, mapping=mapping)


def __get_figure_for_drawing(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg


def _draw_histogram(window, data, histogram_canvas=None, toolbar=None):
    num_bins = 100
    canvas = window["-BIO_APPROVAL_TABL_CANVAS-"].TKCanvas
    fig, ax = plt.subplots()

    ax.hist(data, num_bins, density=True)

    # Format the plots
    plt.xlabel("Activity")
    plt.ylabel("Amount")
    plt.title("Histogram")
    # plt.legend(loc='upper right')

    plot_style = __get_figure_for_drawing(canvas, fig)

    fig = plt.gcf()
    plt.close(fig)

    try:
        histogram_canvas.get_tk_widget().forget()
    except AttributeError:
        print("Attribute Error for info canvas")

        if histogram_canvas is not None:
            if not histogram_canvas == plot_style:
                fig = plt.gcf()
                plt.close(fig)
                histogram_canvas.get_tk_widget().forget()
            else:
                plot_style.get_tk_widget().forget()

    plot_style.draw()

    if toolbar:
        toolbar.destroy()

    toolbar = Toolbar(plot_style, window["-BIO_APPROVAL_TABL_TOOLBAR-"].TKCanvas)
    toolbar.update()
    plot_style.get_tk_widget().pack()

    try:
        window.refresh()
    except AttributeError:
        print("Canvas - AttributeError on window.refresh")

    histogram_canvas = plot_style
    return histogram_canvas, toolbar


def bio_data_approval_table(draw_plate, config, all_plate_data, z_prime_threshold, plate_size, same_layout=True):
    """
    The controller for the pop-up that comes when you are importing plate-reader data to the database

    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param all_plate_data: The raw and calculated data for all the plates
    :type all_plate_data: dict
    :param z_prime_threshold: The threshold for auto approved plates by the software. Should be set in the assay
    :type z_prime_threshold: float
    :return: The layout for the popup
    """
    # ToDO add the option to click colours and change them in the different mappings!
    table_row = None  # Make sure that it does not try to calculate stuff on empty data
    plate_analyse_methods = []  # What method are used for analysing plates
    plate_calculations = []  # What calculations have been done for each plate
    plate_table_data = []  # The data for all the plates
    compound_table_data = []  # The data for all the wells in the plates
    hist_data = {"all": []}  # The data for drawing the histogram
    well_state_mapping = {}  # The mapping of states to a list of wells
    well_state_overview = {  # A dict over each state if it is included in the analysis or not
        "sample": False,
        "blank": False,
        "max": False,
        "minimum": False,
        "pos": False,
        "neg": False,
        "empty": False,
    }

    # Sets stuff for drawing the histogram
    histogram_canvas = None
    toolbar = None
    last_plate_name = None  # Make sure to only update the histogram if a new plate is selected

    BLANK_BOX = '☐'
    CHECKED_BOX = '☑'

    draw_options = ["state_mapping", "heatmap", "hit_map"]

    approved_data = {"plate_name": str,
                     "z_prime": float,
                     "notes": str,
                     "approved": bool,
                     "plate_layout": str,
                     "new_layout": dict}

    # Not sure if there needs to be other methods, but for now there is one
    analyse_methods = ["Histogram"]
    compound_data_dict = {}
    # Grabs data for the table
    for plate_index, plate in enumerate(all_plate_data):
        compound_data_dict[plate] = {}
        for analysed_methods in all_plate_data[plate]["plates"]:
            # Finds the way the plate have been analysed
            if analysed_methods not in plate_analyse_methods:
                plate_analyse_methods.append(analysed_methods)
            compound_data_dict[plate][analysed_methods] = []
            for status in all_plate_data[plate]["plates"][analysed_methods]:
                if status == "wells":
                    for well in all_plate_data[plate]["plates"][analysed_methods][status]:
                        well_score = round(float(all_plate_data[plate]["plates"][analysed_methods]["wells"][well]), 2)
                        temp_well_data = [plate, well, well_score]
                        compound_data_dict[plate][analysed_methods].append(temp_well_data)
                        compound_table_data.append(temp_well_data)
                else:
                    # Stop af the first plate, as all the plates should have same layout from the start .
                    # If they are different layouts I need to re-think this :D
                    if plate_index == 0 and same_layout:
                        if status not in well_state_mapping:
                            well_state_overview[status] = True
                            if same_layout:
                                temp_plate = "all"
                            else:
                                temp_plate = plate
                            well_state_mapping[temp_plate] = {status: []}
                            for wells in all_plate_data[plate]["plates"][analysed_methods][status]:
                                well_state_mapping[temp_plate][status].append(wells)
        if plate_index == 0:
            for calc_index, calculations in enumerate(all_plate_data[plate]["calculations"]):
                if calc_index == 0:
                    # Finds the different calculations that have been done on the plate
                    if calculations != "other":
                        for temp_counter, analysed_method in enumerate(
                                all_plate_data[plate]["calculations"][calculations]):
                            if analysed_method not in calculations and analysed_method not in plate_calculations:
                                if analysed_method != "other":
                                    plate_calculations.append(analysed_method)
        temp_plate_name = plate
        temp_z_prime = round(float(all_plate_data[plate]["calculations"]["other"]["z_prime"]), 2)
        if temp_z_prime > z_prime_threshold:
            temp_approval = CHECKED_BOX
        else:
            temp_approval = BLANK_BOX
        temp_row_data = [temp_plate_name, temp_z_prime, "", temp_approval]
        plate_table_data.append(temp_row_data)

    # Gets the data for each well that is marked as a sample on the last analysed method, for each and all plates
    for plate in all_plate_data:
        hist_data[plate] = []
        if same_layout:
            temp_plate = "all"
        else:
            temp_plate = plate
        for well in well_state_mapping[temp_plate]["sample"]:
            data = all_plate_data[plate]["plates"][plate_analyse_methods[-1]]["wells"][well]
            hist_data[plate].append(data)
            hist_data["all"].append(data)

    # Headlines for the table
    plate_table_headings = ["Plate", "Z-prime", "Notes", "Approved"]
    compound_table_headings = ["Plate", "Well", "Score"]

    # Gets the layout
    window, plate_table_data = _bio_data_approval_table_layout(config, plate_table_data, plate_table_headings,
                                                               compound_table_data, compound_table_headings,
                                                               plate_analyse_methods, plate_calculations,
                                                               analyse_methods, draw_options, well_state_overview)

    # Makes it possible to double-click on the table
    window["-BIO_APPROVAL_PLATE_TABLE-"].bind('<Double-Button-1>', "+-double click-")

    graph = window["-BIO_APPROVAL_TABL_GRAPH-"]
    plate_archive_draw = True
    graph_placement = "bio_approval_popup"

    # Draw the histogram for all the data:
    histogram_canvas, toolbar = _draw_histogram(window, hist_data["all"], histogram_canvas, toolbar)

    while True:
        event, values = window.read()

        if event == sg.WIN_CLOSED or event == "-WINDOW_TWO_CANCEL-" or event == "-POP_SAMPLE_CHECKER_OK-":
            # Grabs all the data from the table
            exit_check = sg.PopupYesNo(f"You are about to close the window.\n "
                                       f"None of the data will be saved and you have to start the import over.\n\n"
                                       f"Do you wish to continue ?")

            if exit_check.casefold() == "yes":
                window.close()
                return "cancelled"

        # ToDo make it possible to de-select wells on the plate-layout and then re-calculate everything with this layout
        if event == "-BIO_APPROVAL_TABLE_RE_CALCULATE-":
            sg.Popup("Does not work. missing a way to de-select data from the plate mapping or the compound table")

        # Not sure what to use this button for ? Delte it ?
        if event == "-BIO_APPROVAL_TABLE_APPLY-":
            sg.Popup("Does not work - not sure what it suppose to do ? ")

        # Sends the data to the database
        if event == "-BIO_DATA_APPROVED-":
            temp_plate_errors = []
            for rows in plate_table_data:
                if rows[3] == "☐" and rows[2] == "":
                    temp_plate_errors.append(rows[0])

            if temp_plate_errors:
                sg.PopupError(f"The following plates are missing notes:\n\n"
                              f"{temp_plate_errors}\n\n"
                              f"You need to add notes to them, before they can be added to the database")

            else:
                # Grabs all the data from the table
                for row_index, row in enumerate(plate_table_data):
                    temp_plate = row[0]

                    # Add the note field to all plates, to make it more consisten to loop through later on
                    if row[2] == "":
                        all_plate_data[temp_plate]["Note"] = None

                    if row[3] == "☐":
                        approved = True
                    else:
                        approved = False

                    all_plate_data[temp_plate]["Approved"] = approved

                window.close()

                return all_plate_data

        # Add a note to the selected plate
        if event == "-BIO_APPROVAL_TABLE_ADD_NOTE-":
            if table_row is None:
                sg.PopupError("Please select a plate")
            else:
                selected_plate = plate_table_data[table_row][0]
                if values["-BIO_APPROVAL_PLATE_NOTE-"]:
                    temp_note = values["-BIO_APPROVAL_PLATE_NOTE-"]
                else:
                    temp_note = sg.PopupGetText("Note:")

                if temp_note:
                    plate_table_data[table_row][2] = "Note"
                    all_plate_data[selected_plate]["Note"] = temp_note

                    window["-BIO_APPROVAL_PLATE_TABLE-"].update(values=plate_table_data)
                    window["-BIO_APPROVAL_PLATE_NOTE-"].update(value=temp_note)

        # Add notes to all plates
        if event == "-BIO_APPROVAL_TABLE_ADD_NOTE_ALL-":
            if values["-BIO_APPROVAL_PLATE_NOTE-"]:
                temp_note = values["-BIO_APPROVAL_PLATE_NOTE-"]
            else:
                temp_note = sg.PopupGetText("Note:")

            if temp_note:

                guard = sg.PopupYesNo(f"This will add the following to all plates:"
                                      f"\n {temp_note}\n\n "
                                      f"Do you wish to continue?")

                if guard.casefold() == "yes":
                    for row_index, plates in enumerate(all_plate_data):
                        plate_table_data[row_index][2] = "Note"
                        all_plate_data[plates]["Note"] = temp_note

                    window["-BIO_APPROVAL_PLATE_TABLE-"].update(values=plate_table_data)
                    window["-BIO_APPROVAL_PLATE_NOTE-"].update(value=temp_note)

        # Upate the calculations in the window, when you change the calculation in the dropdown menu
        if event == "-BIO_APPROVAL_TABLE_DROPDOWN_CALCULATIONS-" and table_row is not None:
            current_plate = plate_table_data[table_row][0]
            temp_analysed_method = values["-BIO_APPROVAL_TABLE_DROPDOWN_PLATE_ANALYSE_METHODS-"]

            states = values["-BIO_APPROVAL_TABLE_DROPDOWN_CALCULATIONS-"]
            all_calc = all_plate_data[current_plate]["calculations"][temp_analysed_method][states]
            window["-BIO_APPROVAL_TABLE_CALC_AVG-"].update(value=all_calc["avg"])
            window["-BIO_APPROVAL_TABLE_CALC_STDEV-"].update(value=all_calc["stdev"])
            window["-BIO_APPROVAL_TABLE_CALC_PSTDEV-"].update(value=all_calc["pstdev"])
            window["-BIO_APPROVAL_TABLE_CALC_PVARIANCE-"].update(value=all_calc["pvariance"])
            window["-BIO_APPROVAL_TABLE_CALC_VARIANCE-"].update(value=all_calc["variance"])
            window["-BIO_APPROVAL_TABLE_CALC_ST_DEV-"].update(value=all_calc["st_dev_%"])

        # Makes it possible to check and un-check a checkbox.
        if event[0] == "-BIO_APPROVAL_PLATE_TABLE-" and event[2][1] == 3:
            row_data = event[2][0]

            if plate_table_data[row_data][3] == CHECKED_BOX:

                plate_table_data[row_data][3] = BLANK_BOX
            else:
                plate_table_data[row_data][3] = CHECKED_BOX
            # if table_data[row + 1][0]:
            #     table_data[row + 1][0] = False
            # else:
            #
            #     table_data[row + 1][0] = True

            window['-BIO_APPROVAL_PLATE_TABLE-'].update(values=plate_table_data)

        # updates the window, when double-clicking the plate table, or changing what data to look at
        if event == "-BIO_APPROVAL_PLATE_TABLE-+-double click-" or event == "-BIO_APPROVAL_TABL_SAMPLE-" or \
                event == "-BIO_APPROVAL_TABL_BLANK-" or event == "-BIO_APPROVAL_TABL_MAX-" or \
                event == "-BIO_APPROVAL_TABL_MIN-" or event == "-BIO_APPROVAL_TABL_POSITIVE-" or \
                event == "-BIO_APPROVAL_TABL_NEGATIVE-" or event == "-BIO_APPROVAL_TABLE_DROPDOWN_DRAW_OPTIONS-" or \
                event == "-BIO_APPROVAL_TABLE_DROPDOWN_PLATE_ANALYSE_METHODS-":

            try:
                table_row = values["-BIO_APPROVAL_PLATE_TABLE-"][0]
            except IndexError:
                pass
            else:

                # Grab data
                temp_plate_name = plate_table_data[table_row][0]
                temp_analysed_method = values["-BIO_APPROVAL_TABLE_DROPDOWN_PLATE_ANALYSE_METHODS-"]

                # Get calculations:
                states = values["-BIO_APPROVAL_TABLE_DROPDOWN_CALCULATIONS-"]
                all_calc = all_plate_data[temp_plate_name]["calculations"][temp_analysed_method][states]
                window["-BIO_APPROVAL_TABLE_CALC_AVG-"].update(value=all_calc["avg"])
                window["-BIO_APPROVAL_TABLE_CALC_STDEV-"].update(value=all_calc["stdev"])
                window["-BIO_APPROVAL_TABLE_CALC_PSTDEV-"].update(value=all_calc["pstdev"])
                window["-BIO_APPROVAL_TABLE_CALC_PVARIANCE-"].update(value=all_calc["pvariance"])
                window["-BIO_APPROVAL_TABLE_CALC_VARIANCE-"].update(value=all_calc["variance"])
                window["-BIO_APPROVAL_TABLE_CALC_ST_DEV-"].update(value=all_calc["st_dev_%"])
                window["-BIO_APPROVAL_TABLE_CALC_SB-"]. \
                    update(value=all_plate_data[temp_plate_name]["calculations"][temp_analysed_method]["other"]["S/B"])

                # Gets what well states should be included in the mapping.
                show_state_list = {"sample": values["-BIO_APPROVAL_TABL_SAMPLE-"],
                                   "blank": values["-BIO_APPROVAL_TABL_BLANK-"],
                                   "max": values["-BIO_APPROVAL_TABL_MAX-"],
                                   "minimum": values["-BIO_APPROVAL_TABL_MIN-"],
                                   "positive": values["-BIO_APPROVAL_TABL_POSITIVE-"],
                                   "negative": values["-BIO_APPROVAL_TABL_NEGATIVE-"],
                                   "empty": values["-BIO_APPROVAL_TABL_EMPTY-"]}
                draw_options = values["-BIO_APPROVAL_TABLE_DROPDOWN_DRAW_OPTIONS-"]

                if not any(show_state_list.values()):
                    sg.PopupError("Please select states to include in the graph")
                else:
                    # Updates the plate map with selected values
                    _draw_plate_on_graph(draw_plate, config, all_plate_data[temp_plate_name]["plates"],
                                         temp_analysed_method, graph, plate_size, graph_placement, show_state_list,
                                         draw_options)

                    # Updates the compound table, with compounds corresponding to the selected
                    # values for the plate mapping
                    new_plate_data = []
                    for states in show_state_list:
                        if show_state_list[states]:
                            for data in compound_data_dict[temp_plate_name][temp_analysed_method]:
                                if same_layout:
                                    temp_plate = "all"
                                else:
                                    temp_plate = temp_plate_name
                                if data[1] in well_state_mapping[temp_plate][states]:
                                    new_plate_data.append(data)

                    window["-BIO_APPROVAL_COMPOUND_TABLE-"].update(values=new_plate_data)

                # Draw histrogram
                if values["-BIO_APPROVAL_TABLE_ANALYSE_METHODS-"] == "Histogram" and temp_plate_name != last_plate_name:
                    histogram_canvas, toolbar = _draw_histogram(window, hist_data[temp_plate_name], histogram_canvas,
                                                                toolbar)
                    last_plate_name = temp_plate_name


if __name__ == "__main__":
    right_headlines = ["source_plates", "destination_plates", "source_well", "destination_well"]
    wrong_headlines = ["source_plates", "destination_plates", "source_well", "destination_well"]
    new_headlines_popup(right_headlines, wrong_headlines)
