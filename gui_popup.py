import PySimpleGUI as sg
from copy import deepcopy, copy


from info import matrix_header, plate_384_row, plate_96_row

# from gui_functions import sub_settings_matrix
# from database_controller import FetchData
# from excel_handler import purity_sample_layout_import, purity_sample_layout_export
# from gui_functions import sort_table

# todo FIX THIS WHOLE THING!!!!!!!!
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
                new_name[table_data[row_index][2]] = {"raw": table_data[row_index][0], "excel": table_data[row_index][1]}

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


def _bio_data_approval_table_layout(config, table_data, headings, plate_analyse_methods, plate_calculations,
                                    analyse_methods, draw_options):
    """
    The layout for the Bio Approval table that comes up when you import bio data to the database
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param table_data: The data for the table
    :type table_data: list
    :param headings: the headings for the table
    :type headings: list
    :param plate_analyse_methods: The way the plate have been calculated
    :type plate_analyse_methods: list
    :param plate_calculations: The different calculations made for each plate
    :type plate_calculations: list
    :param analyse_methods: Different ways to analyse a plate, like historgrams or gragps...
    :type analyse_methods: list
    :return: The layout for the popup
    """
    standard_size = 12

    raw_table_col = sg.Frame("Please approve or dimiss plates", [[
        sg.Column([
            [sg.Table(values=table_data, headings=headings,
                      key="-POP_HEADLINE_TABLE-", enable_events=True, enable_click_events=True,
                      tooltip='double click to select a new headline in the "NEW headline" column ')]
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
                         sg.Checkbox("", key="-BIO_APPROVAL_TABL_SAMPLE-", default=True)],
                        [sg.Text("Blank:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["blank"], size=10,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_BLANK-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABL_BLANK-")],
                        [sg.Text("Maximum:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["max"], size=10,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_NAX-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABL_MAX-")],
                        [sg.Text("Minimum:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["minimum"], size=10,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_MINIMUM-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABL_MIN-")],
                        [sg.Text("Positive Control:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["positive"], size=10,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_POSITIVE-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABL_POSITIVE-")],
                        [sg.Text("Negative Control:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["negative"], size=10,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_NEGATIVE-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABL_NEGATIVE-")],
                        [sg.Text("Empty:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["empty"], size=10,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_EMPTY-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABL_EMPTY-")]
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
        sg.Frame("", [
            [sg.T("Calculations")],
            [sg.T("avg", size=10), sg.T("calc", key="-BIO_APPROVAL_TABLE_CALC_AVG-", size=10)],
            [sg.T("stdev", size=10), sg.T("calc", key="-BIO_APPROVAL_TABLE_CALC_STDEV-", size=10)],
            [sg.T("pstdev", size=10), sg.T("calc", key="-BIO_APPROVAL_TABLE_CALC_PSTDEV-", size=10)],
            [sg.T("pvariance", size=10), sg.T("calc", key="-BIO_APPROVAL_TABLE_CALC_PVARIANCE-", size=10)],
            [sg.T("variance", size=10), sg.T("calc", key="-BIO_APPROVAL_TABLE_CALC_VARIANCE-", size=10)],
            [sg.T("st_dev_%", size=10), sg.T("calc", key="-BIO_APPROVAL_TABLE_CALC_ST_DEV-", size=10)],
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
                           "This will overwrite the old calculated data, but raw data is still saved")]

    ])

    canvas_window = sg.Frame("Data", [[
        sg.Column([
            [sg.Canvas(key="-BIO_APPROVAL_TABL_GRAPH-", background_color="Grey", size=(500, 300))],
            [sg.DropDown(values=analyse_methods, key="-BIO_APPROVAL_TABLE_ANALYSE_METHODS-", enable_events=True,
                         default_value=analyse_methods[0])]
        ])
    ]])


    layout = [sg.vtop([raw_table_col, graph_window, canvas_window]),
        [sg.Button("Done", key="-BIO_DATA_APPROVED-", expand_x=True),
         sg.Button("Cancel", key="-WINDOW_TWO_CANCEL-", expand_x=True)]
    ]
    return sg.Window("Samples", layout, finalize=True, resizable=True), table_data


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


def bio_data_approval_table(draw_plate, config, all_plate_data, z_prime_threshold, plate_size):
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
    plate_analyse_methods = []
    plate_calculations = []
    table_data = []

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

    # Grabs data for the table
    for plate in all_plate_data:
        for analysed_methods in all_plate_data[plate]["plates"]:
            # Finds the way the plate have been alaysed
            plate_analyse_methods.append(analysed_methods)
        for calculations in all_plate_data[plate]["calculations"]:

            # Finds the different calculations that have been done on the plate
            if calculations != "other":
                for analysed_method in all_plate_data[plate]["calculations"][calculations]:
                    if analysed_method not in calculations:
                        plate_calculations.append(analysed_method)
        temp_plate_name = plate
        temp_z_prime = all_plate_data[plate]["calculations"]["other"]["z_prime"]
        if temp_z_prime > z_prime_threshold:
            temp_approval = CHECKED_BOX
        else:
            temp_approval = BLANK_BOX
        temp_row_data = [temp_plate_name, temp_z_prime, "", temp_approval]
        table_data.append(temp_row_data)

    # Headlines for the table
    table_headings = ["Plate", "Z-prime", "Notes", "Approved"]

    window, table_data = _bio_data_approval_table_layout(config, table_data, table_headings, plate_analyse_methods,
                                                         plate_calculations, analyse_methods, draw_options)
    window["-POP_HEADLINE_TABLE-"].bind('<Double-Button-1>', "+-double click-")

    graph = window["-BIO_APPROVAL_TABL_GRAPH-"]
    plate_archive_draw = True
    graph_placement = "bio_approval_popup"

    while True:
        event, values = window.read()

        if event == sg.WIN_CLOSED or event == "-WINDOW_TWO_CANCEL-" or event == "-POP_SAMPLE_CHECKER_OK-":
            # Grabs all the data from the table
            for row_index, row in enumerate(table_data):
                print(row)
                if type(row[1]) == list:
                    approved_data[row[0]] = row[1][0]
                else:
                    approved_data[row[0]] = row[1]

            window.close()
            return approved_data

        # Makes it possible to check and un-check a checkbox.
        if event[0] == "-POP_HEADLINE_TABLE-" and event[2][1] == 3:
            row_data = event[2][0]

            if table_data[row_data][3] == CHECKED_BOX:

                table_data[row_data][3] = BLANK_BOX
            else:
                table_data[row_data][3] = CHECKED_BOX
            # if table_data[row + 1][0]:
            #     table_data[row + 1][0] = False
            # else:
            #
            #     table_data[row + 1][0] = True

            window['-POP_HEADLINE_TABLE-'].update(values=table_data)

        if event == "-POP_HEADLINE_TABLE-+-double click-":

            try:
                table_row = values["-POP_HEADLINE_TABLE-"][0]
            except IndexError:
                pass
            else:
                # print(table_data[table_row])
                # Grab data
                temp_plate_name = table_data[table_row][0]
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
                # Draw the plate:

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
                    _draw_plate_on_graph(draw_plate, config, all_plate_data[temp_plate_name]["plates"],
                                         temp_analysed_method, graph, plate_size, graph_placement, show_state_list,
                                         draw_options)




if __name__ == "__main__":
    right_headlines = ["source_plates", "destination_plates", "source_well", "destination_well"]
    wrong_headlines = ["source_plates", "destination_plates", "source_well", "destination_well"]
    new_headlines_popup(right_headlines, wrong_headlines)