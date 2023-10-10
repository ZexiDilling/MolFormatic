import PySimpleGUI as sg
import copy
import matplotlib

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from extra_functions import increment_text_string
from xml_handler import convert_echo_to_db
from bio_date_handler import BIOAnalyser
from database_handler import DataBaseFunctions
from info import matrix_header, plate_384_row, plate_96_row
from visualization import Toolbar

matplotlib.use('TkAgg')


def popup_select(the_list, select_multiple=False):
    layout = [[
        sg.Listbox(the_list,key="-LIST_ITEMS-",size=(45,len(the_list)),
                   select_mode="extended" if select_multiple else "single",bind_return_key=True),sg.OK()
    ]]
    window = sg.Window('Select One',layout=layout)
    event, values = window.read()
    window.close()
    del window
    if select_multiple or values['-LIST_ITEMS-'] is None:
        return values['-LIST_ITEMS-']
    else:
        return values['-LIST_ITEMS-'][0]


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


def matrix_popup(data_dict, calc_values, state_values, method_values, calc, sub_settings_matrix, state=None,
                 method=None):
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


def sample_to_compound_name_controller(config, data_dict, fd, purity_sample_layout_import, purity_sample_layout_export,
                                       sort_table, db_data=True):
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

    default_table_data = copy.deepcopy(table_data)

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


def new_headlines_popup(sort_table, right_headlines, wrong_headlines):
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
    # ToDo Make it possible to choose multiple file to change layout for at the same time
    table_data = []
    plate_to_layout = {}
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
                    plate_to_layout[row[0]] = row[1][0]
                else:
                    plate_to_layout[row[0]] = row[1]

            window.close()
            return plate_to_layout

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


def _assay_generator_layout(plate_layout):
    text_size = 12
    input_size = 10
    assay_layout = sg.Frame("New Assay", [[
        sg.Column([
            [sg.T("Assay Name", size=text_size),
             sg.Input("", key="-NEW_ASSAY_NAME-", size=input_size),
             sg.T("Plate Layout", size=text_size),
             sg.DropDown(plate_layout, key="-NEW_ASSAY_PLATE_LAYOUT-", size=input_size)],
            [sg.T("Z-Prime Threshold", size=text_size),
             sg.Input("", key="-NEW_ASSAY_Z_PRIME-", size=input_size,
                      tooltip="The default threshold that plates will compare it self to. Can be changed. "
                              "Can be left empty"), sg.T("Hit Threshold", size=text_size),
             sg.Input("", key="-NEW_ASSAY_HIT-", size=input_size,
                      tooltip="The default Hit Score that samples will compare it self to. Can be changed. "
                              "Can be left empty")],
            [sg.HorizontalSeparator()],
            [sg.T("SOP")],
            [sg.Multiline(key="-NEW_ASSAY_SOP-", size=(50, 25))],

        ])
    ]])

    layout = [
        [assay_layout],
        [sg.Button("Done", key="-NEW_ASSAY_DONE-", expand_x=text_size), sg.Push(),
         sg.Button("Cancel", key="-WINDOW_TWO_CANCEL-", expand_x=text_size)]
    ]

    return sg.Window("New Assay", layout, finalize=True, resizable=True)


def assay_generator(config, plate_list):
    window = _assay_generator_layout(plate_list)

    while True:
        event, values = window.read()

        if event == sg.WIN_CLOSED or event == "-WINDOW_TWO_CANCEL-":

            window.close()
            return False

        elif event == "-NEW_ASSAY_DONE-":
            if not values["-NEW_ASSAY_NAME-"]:
                sg.PopupError("Please provide an assay name")
            elif not values["-NEW_ASSAY_PLATE_LAYOUT-"]:
                sg.PopupError("Please Select a plate-layout")
            else:
                dbf = DataBaseFunctions(config)

                # Gets values
                assay_table = "assay"
                row_id = dbf.number_of_rows(assay_table)
                assay_name = values["-NEW_ASSAY_NAME-"]
                sop = values["-NEW_ASSAY_SOP-"]
                assay_plate_layout = values["-NEW_ASSAY_PLATE_LAYOUT-"]
                z_prime_threshold = values["-NEW_ASSAY_Z_PRIME-"]
                hit_threshold = values["-NEW_ASSAY_HIT-"]
                plate_data = dbf.find_data_single_lookup("plate_layout", assay_plate_layout, "plate_name")
                plate_size = plate_data[0][3]

                assay_data = {"row_id": row_id,
                              "assay_name": assay_name,
                              "SOP": sop,
                              "plate_layout": assay_plate_layout,
                              "plates_run": 0,
                              "compounds_run": 0,
                              "z_prime_threshold": z_prime_threshold,
                              "hit_threshold": hit_threshold,
                              "plate_size": plate_size}

                dbf.add_records_controller(assay_table, assay_data)

                window.close()
                return True


def _assay_run_naming_layout(plate_table_headline, run_name, previous_runs, assay_name, used_plates_table_data,
                             all_batch_numbers, batch_number):
    text_size = 15
    input_size = 15
    previous_runs.append("New")
    top_layout = sg.Column([
        [sg.T("Run Name", size=text_size),
         sg.Input(run_name, key="-ASSAY_RUN_NAME-", size=input_size)],
        [sg.T("Chose old Run", size=text_size),
         sg.DropDown(previous_runs, key="-ASSAY_RUN_PREVIOUS-",
                     size=input_size, enable_events=True, default_value=previous_runs[-1])],
        [sg.B("New batch", key="-ASSAY_RUN_NEW_BATCH-", size=text_size),
         sg.T(batch_number, key="-ASSAY_RUN_CURRENT_BATCH-", size=text_size, visible=False),
         sg.DropDown(values=all_batch_numbers, key="-ASSAY_RUN_ALL_BATCHES-", size=text_size,
                     default_value=batch_number, enable_events=True,
                     tooltip="A way to split up data, "
                             "if there are changes to the data that are specific to some run.")],
        [sg.CalendarButton("Start date", key="-ASSAY_RUN_DATE-", format="%Y-%m-%d", enable_events=True,
                           target="-ASSAY_RUN_DATE_TARGET-", size=(text_size, 1)),
         sg.Input(key="-ASSAY_RUN_DATE_TARGET-", enable_events=True, size=text_size)],
        [sg.B("Worklist", key="-ASSAY_RUN_WORKLIST-", size=text_size),
         sg.T("No Worklist", key="-ASSAY_RUN_WORKLIST_INDICATOR-", relief="sunken", size=text_size),
         sg.Multiline(visible=False, key="-ASSAY_RUN_WORKLIST_DATA-")],
        [sg.B("Echo Data", key="-ASSAY_RUN_ECHO-", size=text_size),
         sg.T("No Echo Data", key="-ASSAY_RUN_ECHO_INDICATOR-", relief="sunken", size=text_size),
         sg.Multiline(visible=False, key="-ASSAY_RUN_ECHO_DATA-")],

        [sg.HorizontalSeparator()]
    ])
    bottom_layer = sg.vtop([
        sg.Column([
            [sg.T("Notes")],
            [sg.Multiline(key="-ASSAY_RUN_NOTES-", size=(18, 15))],
            [sg.Button("Update Run", key="-ASSAY_RUN_UPDATE-", size=text_size)]
        ]),
        sg.Column([
            [sg.T("All Plates")],
            [sg.Table(values=used_plates_table_data, size=(18, 25), headings=plate_table_headline,
                      key="-ASSAY_RUN_USED_PLATES_TABLE-")],
            [sg.B("Selected", key="-ASSAY_RUN_APPLY_SELECTED-", size=10,
                  tooltip="Apply the Run name to the selected plates"),
             sg.B("All", key="-ASSAY_RUN_APPLY_ALL-", size=7,
                  tooltip="Apply the Run name to all plates")]
        ])
    ])

    assay_layout = sg.Frame(f"New Run for {assay_name}", [
        [top_layout],
        bottom_layer
    ])

    layout = [
        [assay_layout],
        [sg.Button("Apply", key="-ASSAY_RUN_DONE-", expand_x=text_size,
                   tooltip="Apply the assay_run name to all plates, and continue with importing the data"),
         sg.Button("Dismiss Runs", key="-ASSAY_RUN_DISMISS-", expand_x=text_size,
                   tooltip="Will ask for assay_run name for witch plates to dismiss,"
                           "if there is more than one run name"
                           "Notes needs to be filled in for the assay_run"),
         sg.Button("Cancel", key="-WINDOW_TWO_CANCEL-", expand_x=text_size)]
    ]

    return sg.Window("Assay Run", layout, finalize=True, resizable=True)


def _previous_runs_data(dbf, assay_name):
    temp_rows = dbf.find_data_single_lookup("assay_runs", assay_name, "assay_name")
    temp_run_name = None
    previous_runs = []
    all_batch_numbers = []
    if temp_rows:
        for data in temp_rows:
            temp_run_name = data[1]
            temp_batch_number = data[3]
            previous_runs.append(temp_run_name)
            all_batch_numbers.append(temp_batch_number)

    if temp_run_name is not None:
        run_name = increment_text_string(temp_run_name)
        batch_number = all_batch_numbers[-1]

    else:
        run_name = f"{assay_name}_run_1"
        batch_number = "Batch_1"
        all_batch_numbers.append(batch_number)

    return previous_runs, run_name, all_batch_numbers, batch_number


def _previous_single_run_data(dbf, assay_name, run_name):
    single_run_data = dbf.find_data_double_lookup("assay_runs", assay_name, run_name, "assay_name", "run_name")
    worklist = single_run_data[0][4]
    echo_data = single_run_data[0][5]
    date = single_run_data[0][6]
    note = single_run_data[0][7]
    return worklist, echo_data, date, note


def _handle_worklist(config, worklist_file, bio_compound_info_from_worklist, worklist_echo_data, run_name):
    worklist_data, _ = bio_compound_info_from_worklist(config, sg, [worklist_file])
    worklist_echo_data["worklist"][run_name] = worklist_data
    worklist_txt_data = str(worklist_data)
    return worklist_txt_data


def _handle_echo_data(echo_files_string, worklist_echo_data, run_name, transfer_dict):
    echo_file = echo_files_string.split(";")
    echo_data = convert_echo_to_db(echo_file)
    if transfer_dict:
        for plates in echo_data:
            transfer_dict[plates] = echo_data[plates]
    else:
        transfer_dict = copy.deepcopy(echo_data)
    worklist_echo_data["echo"][run_name] = echo_data
    echo_data = str(echo_data)

    # echo_data layout:
    # temp = {"plate_name":
    #             {"destination_well":
    #                  {"mp_source_plate": str,
    #                   "mp_source_well": str,
    #                   "reason": str}}}

    return echo_data, transfer_dict


def _update_assay_run_database_data(dbf, current_run_name, batch, worklist, echo_data, date, note):
    source_table = "assay_runs"
    table_data = {
        "batch": batch,
        "worklist": worklist,
        "echo_data": echo_data,
        "date": date,
        "note": note
    }
    table_index_value = current_run_name
    headline = "run_name"

    dbf.update_database_items(source_table, table_data, table_index_value, headline)


def assay_run_naming(config, all_plates_data, analysis_method, used_plates, assay_name, all_destination_plates,
                     dead_run_check, bio_compound_info_from_worklist):

    dbf = DataBaseFunctions(config)
    previous_runs, run_name, all_batch_numbers, batch_number = _previous_runs_data(dbf, assay_name)
    used_plates_table_data = []
    if dead_run_check.casefold() == "yes":
        plate_table_headline = ["Plate", "Run", "Data"]
        plate_listed = all_destination_plates
        for plates in all_destination_plates:
            if plates in used_plates:
                data_check = "Got Data"
            else:
                data_check = "No Data"
            temp_data = [plates, "", data_check]
            used_plates_table_data.append(temp_data)
    else:
        plate_listed = used_plates
        plate_table_headline = ["Plate", "Run"]
        for plates in used_plates:
            temp_data = [plates, ""]
            used_plates_table_data.append(temp_data)

    window = _assay_run_naming_layout(plate_table_headline, run_name, previous_runs, assay_name, used_plates_table_data,
                                      all_batch_numbers, batch_number)
    run_notes = {}
    plates_checked = []
    dismissed_plates = {}
    different_assay_runs = []
    transfer_dict = None
    echo_data = None

    worklist_echo_data = {"worklist": {},
                          "echo": {}}

    while True:
        event, values = window.read()

        if event == sg.WIN_CLOSED or event == "-WINDOW_TWO_CANCEL-":
            window.close()
            return False, None, None

        if event == "-ASSAY_RUN_DONE-":
            can_close = True
            if len(plates_checked) == len(plate_listed):
                for plates in all_plates_data:
                    try:
                        transfer_dict[plates]
                    except KeyError:
                        sg.PopupError(f"{plates} is missing echo data")
                        can_close = False
                    else:
                        all_plates_data[plates]["transferred_wells"] = []
                        for wells in transfer_dict[plates]["transferred_wells"]:
                            all_plates_data[plates]["transferred_wells"].append(wells)
                            if plates in dismissed_plates:
                                dismissed_plates[plates].append(wells)

                if can_close:
                    window.close()
                    return True, transfer_dict, dismissed_plates
            else:
                sg.PopupError("Not all plates have gotten an assay run assigned to them")

        if event == "-ASSAY_RUN_DISMISS-":
            if len(different_assay_runs) == 0:
                sg.PopupError("Please assign a run to the plates that needs to be dismissed")

            elif len(different_assay_runs) == 1:
                if run_notes[values["-ASSAY_RUN_NAME-"]]:
                    dismissed_plates_list = plate_listed

                    for temp_plates in dismissed_plates_list:
                        dismissed_plates[temp_plates] = []

                else:
                    sg.PopupError("Please fill out the Note")

            else:
                used_runs = []
                for rows in used_plates_table_data:
                    temp_assay_runs = rows[1]
                    if temp_assay_runs not in used_runs:
                        used_runs.append(temp_assay_runs)
                dismissed_plates_list = popup_select(used_runs, select_multiple=True)

                for temp_plates in dismissed_plates_list:
                    dismissed_plates[temp_plates] = []

        if event == "-ASSAY_RUN_APPLY_SELECTED-" or event == "-ASSAY_RUN_APPLY_ALL-":
            plate_table_index = values["-ASSAY_RUN_USED_PLATES_TABLE-"]
            run = "assay_run_1"
            temp_plate_data = []
            if event == "-ASSAY_RUN_APPLY_ALL-":
                for plates in used_plates_table_data:
                    temp_plate_data.append(plates[0])
            else:
                for index in plate_table_index:
                    temp_plate_data.append(used_plates_table_data[index][0])
                    used_plates_table_data[index][1] = run

            window["-ASSAY_RUN_USED_PLATES_TABLE-"].update(values=used_plates_table_data)

            temp_event = event
            if not values["-ASSAY_RUN_NAME-"]:
                sg.popup_error("Please select a name for the run")
            elif not values["-ASSAY_RUN_DATE_TARGET-"]:
                sg.popup_error("Please select a date for the run")
            elif not values["-ASSAY_RUN_ECHO_DATA-"]:
                sg.popup_error("Please provide echo data for the transferee")
            elif not values["-ASSAY_RUN_WORKLIST_DATA-"]:
                sg.popup_error("Please provide the worklist used for the run")
            # ToDo have a check to see if the echo data fit with the plate
            else:
                assay_run_name = values["-ASSAY_RUN_NAME-"]
                if assay_run_name not in different_assay_runs:
                    different_assay_runs.append(assay_run_name)
                    run_notes[assay_run_name] = values["-ASSAY_RUN_NOTES-"]
                # Test if the run name is in the database already, if it isn't add the name to the run data
                if assay_run_name not in previous_runs or len(previous_runs) == 1:
                    assay_run_data = {"run_name": assay_run_name,
                                      "assay_name": assay_name,
                                      "batch": values["-ASSAY_RUN_ALL_BATCHES-"],
                                      "worklist": values["-ASSAY_RUN_WORKLIST_DATA-"],
                                      "echo_data": values["-ASSAY_RUN_ECHO_DATA-"],
                                      "date": values["-ASSAY_RUN_DATE_TARGET-"],
                                      "note": values["-ASSAY_RUN_NOTES-"]}

                    dbf.add_records_controller("assay_runs", assay_run_data)
                    # add the plates to assay_plates
                    print("ARE WE HERE ? ")
                    for plates in temp_plate_data:
                        print(plates)
                        assay_plate_data = {"plate_name": plates, "assay_run": assay_run_name}
                        dbf.add_records_controller("assay_plates", assay_plate_data)

                    # Add the new run to the list of old runs,
                    # and updates the dropdown to include it, and the value to it
                    if len(previous_runs) != 1:
                        previous_runs.append(assay_run_name)
                        window["-ASSAY_RUN_PREVIOUS-"].update(values=previous_runs, value=assay_run_name)

                if temp_event == "-ASSAY_RUN_APPLY_SELECTED-":
                    temp_plate_list = []
                    for index in values["-ASSAY_RUN_USED_PLATES_TABLE-"]:
                        temp_plate_list.append(used_plates_table_data[index][0])
                        used_plates_table_data[index][1] = run

                    # Update the table with the assay run next to the plate name
                    window["-ASSAY_RUN_USED_PLATES_TABLE-"].update(values=used_plates_table_data)
                else:
                    temp_plate_list = plate_listed
                    for index, _ in enumerate(used_plates_table_data):
                        used_plates_table_data[index][1] = run

                    # Update the table with the assay run next to the plate name
                    window["-ASSAY_RUN_USED_PLATES_TABLE-"].update(values=used_plates_table_data)

                for plates in temp_plate_list:
                    all_plates_data[plates]["run_name"] = assay_run_name     # Update the dict with the run_name
                    all_plates_data[plates]["analysis_method"] = analysis_method
                    if plates not in plates_checked:        # Make sure that already checked plates are not added
                        plates_checked.append(plates)

                # Reset all the info, and count up one for the name:
                assay_run_name = increment_text_string(assay_run_name)
                window["-ASSAY_RUN_NAME-"].update(value=assay_run_name)
                window["-ASSAY_RUN_DATE_TARGET-"].update(value="Choose date")
                window["-ASSAY_RUN_WORKLIST_INDICATOR-"].update(value="No Worklist")
                window["-ASSAY_RUN_ECHO_INDICATOR-"].update(value="No Echo Data")
                window["-ASSAY_RUN_NOTES-"].update(value="")
                window["-ASSAY_RUN_WORKLIST_DATA-"].update(value="")
                window["-ASSAY_RUN_ECHO_DATA-"].update(value="")

        if event == "-ASSAY_RUN_UPDATE-":
            current_assay_name = values["-ASSAY_RUN_NAME-"]
            if current_assay_name not in previous_runs:
                sg.popup_error(f"{current_assay_name} is not in the Database. "
                               f"Please either select plates or add it to all ")
            else:

                current_run_name = values["-ASSAY_RUN_NAME-"]
                batch = values["-ASSAY_RUN_ALL_BATCHES-"]
                worklist = values["-ASSAY_RUN_WORKLIST_DATA-"]
                echo_data = values["-ASSAY_RUN_ECHO_DATA-"]
                date = values["-ASSAY_RUN_DATE_TARGET-"]
                note = values["-ASSAY_RUN_NOTES-"]

                _update_assay_run_database_data(dbf, current_run_name, batch, worklist, echo_data, date, note)

        if event == "-ASSAY_RUN_WORKLIST-":
            worklist_file = sg.popup_get_file("Please select worklist")
            # ToDo add guard for wrong file-format
            if worklist_file:
                temp_run_name = values["-ASSAY_RUN_NAME-"]
                worklist_txt_data = _handle_worklist(config, worklist_file, bio_compound_info_from_worklist,
                                                     worklist_echo_data, temp_run_name)
                window["-ASSAY_RUN_WORKLIST_DATA-"].update(value=worklist_txt_data)
                window["-ASSAY_RUN_WORKLIST_INDICATOR-"].update(value="Got Worklist")

        if event == "-ASSAY_RUN_ECHO-":
            echo_files_string = sg.popup_get_file("Please Select the files with the Echo Data", multiple_files=True)
            # ToDo add guard for wrong file-format
            if echo_files_string:
                temp_run_name = values["-ASSAY_RUN_NAME-"]
                echo_txt_data, transfer_dict = _handle_echo_data(echo_files_string, worklist_echo_data, temp_run_name,
                                                                 transfer_dict)
                window["-ASSAY_RUN_ECHO_DATA-"].update(value=echo_txt_data)
                window["-ASSAY_RUN_ECHO_INDICATOR-"].update(value="Got ECHO DATA")

        if event == "-ASSAY_RUN_PREVIOUS-":
            run_name = values["-ASSAY_RUN_PREVIOUS-"]
            if run_name.casefold() == "new":
                worklist = ""
                date = ""
                note = ""
            else:
                worklist, echo_data, date, note = _previous_single_run_data(dbf, assay_name, run_name)
                if transfer_dict:
                    for plates in eval(echo_data):
                        transfer_dict[plates] = echo_data[plates]
                else:
                    transfer_dict = copy.deepcopy(eval(echo_data))

            window["-ASSAY_RUN_NAME-"].update(value=run_name)
            window["-ASSAY_RUN_DATE_TARGET-"].update(value=date)
            if worklist:
                window["-ASSAY_RUN_WORKLIST_INDICATOR-"].update(value="Got Worklist")
            window["-ASSAY_RUN_WORKLIST_DATA-"].update(value=worklist)
            if echo_data:
                window["-ASSAY_RUN_ECHO_INDICATOR-"].update(value="Got Echo Data")
            window["-ASSAY_RUN_ECHO_DATA-"].update(value=worklist)
            window["-ASSAY_RUN_NOTES-"].update(value=note)

        if event == "-ASSAY_RUN_NEW_BATCH-":
            new_batch_name = increment_text_string(all_batch_numbers[-1])
            all_batch_numbers.append(new_batch_name)
            window["-ASSAY_RUN_ALL_BATCHES-"].update(values=all_batch_numbers)
            window["-ASSAY_RUN_CURRENT_BATCH-"].update(value=new_batch_name)

        if event == "-ASSAY_RUN_ALL_BATCHES-":
            selected_batch = values["-ASSAY_RUN_ALL_BATCHES-"]
            window["-ASSAY_RUN_CURRENT_BATCH-"].update(value=selected_batch)


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

    state_colour_size = 5
    checkbox_spacer = 35

    graph_window = sg.Frame("Plate Layout", [[
        sg.Graph(canvas_size=(500, 300), graph_bottom_left=(0, 0), graph_top_right=(250, 175),
                 background_color='grey', key="-BIO_APPROVAL_TABLE_GRAPH-", enable_events=True,
                 drag_submits=True, motion_events=True)],
        sg.vtop([
            sg.TabGroup([[
                sg.Tab("State", [[
                    sg.Column([
                        [sg.T("States", size=10), sg.T("Colour", size=10), sg.T("Include", size=6),
                         sg.T("Draw")],
                        [sg.HorizontalSeparator()],

                        [sg.Text("Sample:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["sample"], size=state_colour_size,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_SAMPLE-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABLE_SAMPLE-", enable_events=True,
                                     disabled=not well_state_overview["sample"], default=True,
                                     pad=(checkbox_spacer, 0)),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABLE_SAMPLE_DRAW-", enable_events=True, default=True)],

                        [sg.Text("Blank:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["blank"], size=state_colour_size,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_BLANK-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABLE_BLANK-", enable_events=True,
                                     disabled=not well_state_overview["blank"], pad=(checkbox_spacer, 0)),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABLE_BLANK_DRAW-", enable_events=True)],

                        [sg.Text("Maximum:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["max"], size=state_colour_size,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_NAX-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABLE_MAX-", enable_events=True,
                                     disabled=not well_state_overview["max"], pad=(checkbox_spacer, 0)),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABLE_MAX_DRAW-", enable_events=True)],

                        [sg.Text("Minimum:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["minimum"], size=state_colour_size,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_MINIMUM-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABLE_MIN-", enable_events=True,
                                     disabled=not well_state_overview["minimum"], pad=(checkbox_spacer, 0)),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABLE_MIN_DRAW-", enable_events=True)],

                        [sg.Text("Positive Control:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["positive"], size=state_colour_size,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_POSITIVE-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABLE_POSITIVE-", enable_events=True,
                                     disabled=not well_state_overview["pos"], pad=(checkbox_spacer, 0)),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABLE_POSITIVE_DRAW-", enable_events=True)],

                        [sg.Text("Negative Control:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["negative"], size=state_colour_size,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_NEGATIVE-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABLE_NEGATIVE-", enable_events=True,
                                     disabled=not well_state_overview["neg"], pad=(checkbox_spacer, 0)),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABLE_NEGATIVE_DRAW-", enable_events=True)],

                        [sg.Text("Empty:", size=standard_size),
                         sg.T(background_color=config["plate_colouring"]["empty"], size=state_colour_size,
                              key="-BIO_PLATE_LAYOUT_COLOUR_BOX_EMPTY-", relief="groove"),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABLE_EMPTY-", enable_events=True,
                                     disabled=not well_state_overview["empty"], pad=(checkbox_spacer, 0)),
                         sg.Checkbox("", key="-BIO_APPROVAL_TABLE_EMPTY_DRAW-", enable_events=True)]
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
                [sg.T("Well", key="-BIO_APPROVAL_WELL_ID-", size=8),
                 sg.T("Score:", key="-BIO_APPROVAL_WELL_SCORE-", size=12)],
                [sg.Frame("", [
                    [sg.T("Calculations")],
                    [sg.HorizontalSeparator()],
                    [sg.T("avg", size=10), sg.T("Calc", key="-BIO_APPROVAL_TABLE_CALC_AVG-", size=10)],
                    [sg.T("stdev", size=10), sg.T("Calc", key="-BIO_APPROVAL_TABLE_CALC_STDEV-", size=10)],
                    [sg.T("pstdev", size=10), sg.T("Calc", key="-BIO_APPROVAL_TABLE_CALC_PSTDEV-", size=10)],
                    [sg.T("pvariance", size=10), sg.T("Calc", key="-BIO_APPROVAL_TABLE_CALC_PVARIANCE-", size=10)],
                    [sg.T("variance", size=10), sg.T("Calc", key="-BIO_APPROVAL_TABLE_CALC_VARIANCE-", size=10)],
                    [sg.T("st_dev_%", size=10), sg.T("Calc", key="-BIO_APPROVAL_TABLE_CALC_ST_DEV-", size=10)],
                    [sg.T("S/B", size=10), sg.T("Calc", key="-BIO_APPROVAL_TABLE_CALC_SB-", size=10)],
                    [sg.T("Z-Prime", size=10), sg.T("Calc", key="-BIO_APPROVAL_TABLE_CALC_Z_PRIME-", size=10)]
                ])],
                [sg.Frame("Note", [
                    [sg.Multiline(default_text="", key="-BIO_APPROVAL_PLATE_NOTE-", size=(25, 2))]
                ])]
            ])
        ]),
        [sg.T("Draw style:", size=10), sg.T("Data:", size=10), sg.T("Calculations:", size=10)],
        [sg.DropDown(values=draw_options, key="-BIO_APPROVAL_TABLE_DROPDOWN_DRAW_OPTIONS-", size=10,
                     enable_events=True, default_value=draw_options[2]),
         sg.DropDown(values=plate_analyse_methods, key="-BIO_APPROVAL_TABLE_DROPDOWN_PLATE_ANALYSE_METHODS-",
                     size=10, enable_events=True, default_value=plate_analyse_methods[2],
                     tooltip="The different ways that he plate have been calculated"),
         sg.DropDown(values=plate_calculations, key="-BIO_APPROVAL_TABLE_DROPDOWN_CALCULATIONS-", size=10,
                     enable_events=True, default_value=plate_calculations[3],
                     tooltip="Shows different calculations for the plate")],
        [sg.Button("Re-calculate", key="-BIO_APPROVAL_TABLE_RE_CALCULATE-", size=12,
                   tooltip="Will re-calculate with the new layout. OBS"),
         sg.Button("Apply", key="-BIO_APPROVAL_TABLE_APPLY-", size=12,
                   tooltip="Will apply the new calculations the the data. "
                           "This will overwrite the old calculated data, but raw data is still saved"),
         sg.Button("Reset", key="-BIO_APPROVAL_TABLE_RESET-", size=12,
                   tooltip="This will reset the drawing"),
         sg.Button("Add Note", key="-BIO_APPROVAL_TABLE_ADD_NOTE-", size=12,
                   tooltip="Adds a note to the plate! Any plate not approved needs to have a note."),
         sg.Button("Add Note to all", key="-BIO_APPROVAL_TABLE_ADD_NOTE_ALL-", size=12,
                   tooltip="Adds the same note to all the plates.")
         ]

    ])

    canvas_window = sg.Frame("Data", [[
        sg.Column([
            [sg.Canvas(key="-BIO_APPROVAL_TABLE_TOOLBAR-")],
            [sg.Canvas(key="-BIO_APPROVAL_TABLE_CANVAS-", background_color="Grey", size=(500, 300))],
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
                         show_state_list, draw_option, skipped_wells):

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
    well_dict, min_x, min_y, max_x, max_y = \
        draw_plate(config, graph, plate_size, well_dict, graph_placement, archive_plate=True, mapping=mapping,
                   skipped_well=skipped_wells)

    return well_dict, min_x, min_y, max_x, max_y


def __get_figure_for_drawing(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg


def _get_hist_data(temp_plate_dict, plate_analyse_methods, well_type):
    temp_hist_data = []
    for well in well_type["sample"]:
        well_data = temp_plate_dict["plates"][plate_analyse_methods[-1]]["wells"][well]
        temp_hist_data.append(well_data)

    return temp_hist_data


def _draw_histogram(window, hist_data, histogram_canvas=None, toolbar=None):
    num_bins = 100
    canvas = window["-BIO_APPROVAL_TABLE_CANVAS-"].TKCanvas
    fig, ax = plt.subplots()
    ax.hist(hist_data, num_bins, density=True)

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

    toolbar = Toolbar(plot_style, window["-BIO_APPROVAL_TABLE_TOOLBAR-"].TKCanvas)
    toolbar.update()
    plot_style.get_tk_widget().pack()

    try:
        window.refresh()
    except AttributeError:
        print("Canvas - AttributeError on window.refresh")

    histogram_canvas = plot_style
    return histogram_canvas, toolbar


def _get_well_type(plate_layout, selected_wells):
    well_type = {
        "sample": [],
        "blank": [],
        "max": [],
        "minimum": [],
        "pos": [],
        "neg": [],
        "empty": [],
    }
    for counter in plate_layout:
        temp_well = plate_layout[counter]["well_id"]
        if temp_well in selected_wells:
            temp_state = selected_wells[temp_well]["state"]
        else:
            temp_state = plate_layout[counter]["state"]
        well_type[temp_state].append(temp_well)

    if not well_type["sample"] or not well_type["max"] or not well_type["minimum"]:
        return None
    else:
        return well_type


def _update_plate_calculations(window, values, temp_plate_dict, temp_analysed_method):

    print(temp_plate_dict)
    print(temp_analysed_method)
    print(temp_plate_dict["calculations"][temp_analysed_method])
    states = values["-BIO_APPROVAL_TABLE_DROPDOWN_CALCULATIONS-"]
    all_calc = temp_plate_dict["calculations"][temp_analysed_method][states]
    window["-BIO_APPROVAL_TABLE_CALC_AVG-"].update(value=all_calc["avg"])
    window["-BIO_APPROVAL_TABLE_CALC_STDEV-"].update(value=all_calc["stdev"])
    window["-BIO_APPROVAL_TABLE_CALC_PSTDEV-"].update(value=all_calc["pstdev"])
    window["-BIO_APPROVAL_TABLE_CALC_PVARIANCE-"].update(value=all_calc["pvariance"])
    window["-BIO_APPROVAL_TABLE_CALC_VARIANCE-"].update(value=all_calc["variance"])
    window["-BIO_APPROVAL_TABLE_CALC_ST_DEV-"].update(value=all_calc["st_dev_%"])
    window["-BIO_APPROVAL_TABLE_CALC_SB-"]. \
        update(value=temp_plate_dict["calculations"][temp_analysed_method]["other"]["S/B"])

    window["-BIO_APPROVAL_TABLE_CALC_Z_PRIME-"]. \
        update(value=temp_plate_dict["calculations"]["other"]["z_prime"])

    draw_options = values["-BIO_APPROVAL_TABLE_DROPDOWN_DRAW_OPTIONS-"]

    return draw_options

#
# def __change_plate_layout_naming(dbf, default_plate_layout_name):
#     table = "plate_layout_sub"
#     data_value = default_plate_layout_name
#     headline = "plate_main"
#     print(data_value)
#     rows = dbf.find_data_single_lookup(table, data_value, headline)
#     # Generate name for the new sub layout
#     print(rows)
#     highest_sub_name = rows[-1][1]
#     temp_name_list = highest_sub_name.split("_")
#     temp_counter = temp_name_list[-1]
#
#     try:
#         int(temp_counter)
#     except ValueError:
#         temp_name = highest_sub_name
#         temp_counter = "1"
#     else:
#         temp_counter = int(temp_counter) + 1
#         temp_name = highest_sub_name.removesuffix(str(temp_counter))
#
#     new_sub_name = f"{temp_name}_{temp_counter}"
#
#     return new_sub_name


def __change_plate_layout_dict(well_dict, worklist_echo_data, current_plate_name):
    # Change the state and colour of wells where the echo have not transfereed.
    for counter in well_dict:
        temp_well = well_dict[counter]["well_id"]
        try:
            worklist_echo_data[current_plate_name][temp_well]
        except KeyError:
            pass
        else:
            well_dict[counter]["state"] = "empty"
    return well_dict


# def _change_plata_layout_archive_plate_update(temp_archive_plates_dict, new_sub_name, well_dict, current_layout):
#     # Generates the dict
#     plate_type = temp_archive_plates_dict[current_layout]["plate_type"]
#     temp_archive_plates_dict[new_sub_name] = {"well_layout": well_dict,
#                                               "plate_type": plate_type,
#                                               "sample": [],
#                                               "blank": [],
#                                               "max": [],
#                                               "minimum": [],
#                                               "positive": [],
#                                               "negative": [],
#                                               "empty": []}
#
#     # Makes list of the different well-types and adds them
#     temp_well_dict = well_dict
#     for counter in temp_well_dict:
#         well_id = temp_well_dict[counter]["well_id"]
#         temp_state = temp_well_dict[counter]["state"]
#         temp_archive_plates_dict[new_sub_name][temp_state].append(well_id)
#
#     return temp_archive_plates_dict


def _change_plate_layout(dbf, well_dict, current_plate_name, default_plate_layout_name):
    # Get sub plate outs to add a new
    # new_sub_name = __change_plate_layout_naming(dbf, default_plate_layout_name)
    # get the assay_run for the current plate and fetch the echo_data for that plate to make a plate_layout from.
    rows = dbf.find_data_single_lookup(table="assay_plates", data_value=current_plate_name, headline="plate_name")

    assay_run = rows[0][2]
    rows = dbf.find_data_single_lookup(table="assay_runs", data_value=assay_run, headline="run_name")
    worklist_echo_data = eval(rows[0][4])

    # Make a new plate layout dict
    new_well_dict = __change_plate_layout_dict(well_dict, worklist_echo_data, current_plate_name)

    # return new_sub_name, new_well_dict
    return new_well_dict


def _get_well_dict_data(dbf, plate_layout_name):
    table = "plate_layout_sub"
    data_value = plate_layout_name
    headline = "plate_sub"
    rows = dbf.find_data_single_lookup(table, data_value, headline)
    well_dict = eval(rows[0][3])

    return well_dict


def _update_plate_table(plate_table_data, table_row, z_prime, plate_layout, notes, run_name, approved):
    plate_table_data[table_row][1] = z_prime
    plate_table_data[table_row][2] = plate_layout
    plate_table_data[table_row][3] = notes
    plate_table_data[table_row][4] = run_name
    plate_table_data[table_row][5] = approved

    return plate_table_data


def _update_compound_table(compound_data_dict, temp_archive_plates_dict, current_plate_name, current_layout_name):
    new_compound_table_data = []
    for data in compound_data_dict[current_plate_name]:
        if data[1] in temp_archive_plates_dict[current_layout_name]["sample"] and data[4] == "☑":
            new_compound_table_data.append(data)
        elif data[1] in temp_archive_plates_dict[current_layout_name]["sample"] and data[4] == "☐":
            data[4] = "☑"
            new_compound_table_data.append(data)
        elif data[1] not in temp_archive_plates_dict[current_layout_name]["sample"]:
            data[4] = "☐"
            new_compound_table_data.append(data)
            if not data[3]:
                data[3] = "dismissed"

    return new_compound_table_data


def _plate_layout_controller(dbf, temp_archive_plates_dict, plate_layout_check, current_plate_name, plate_to_layout,
                             all_plates_data):
    temp_plate_layout = plate_to_layout[current_plate_name]

    try:
        temp_archive_plates_dict[temp_plate_layout]
    except KeyError:
        well_dict = _get_well_dict_data(dbf, temp_plate_layout)
    else:
        well_dict = copy.deepcopy(temp_archive_plates_dict[temp_plate_layout]["well_layout"])

        if plate_layout_check.casefold() == "changed":
            well_dict = _change_plate_layout(dbf, well_dict, current_plate_name, temp_plate_layout)

    # Old version where there was made a lot of sub plate layouts
    # if plate_layout_check.casefold() == "changed":
        # new_sub_name, well_dict = _change_plate_layout(dbf, well_dict, current_plate_name, temp_plate_layout)
        # # update temp_archive_plates_dict
        # temp_archive_plates_dict = _change_plata_layout_archive_plate_update(temp_archive_plates_dict,
        #                                                                      new_sub_name, well_dict,
        #                                                                      temp_plate_layout)
        # plate_to_layout[current_plate_name] = new_sub_name

    temp_plate_dict = copy.deepcopy(all_plates_data[current_plate_name])

    return temp_plate_dict, well_dict, temp_archive_plates_dict, temp_plate_layout


def _sort_plate_data(plate_table_data, all_plates_data, sub_layouts):
    # Grabs all the data from the table
    for row_index, row in enumerate(plate_table_data):
        temp_plate = row[0]

        # Add the note field to all plates, to make it more consisten to loop through later on
        if row[3] == "":
            all_plates_data[temp_plate]["note"] = None

        if row[5] == "☑":
            approved = True
        else:
            approved = False

        all_plates_data[temp_plate]["approved"] = approved

    return all_plates_data


def _sort_compound_data(compound_data_dict):
    all_compound_data = {}
    for plates in compound_data_dict:
        all_compound_data[plates] = {}
        for rows in compound_data_dict[plates]:
            well = rows[1]
            note = rows[3]
            approved = rows[4]
            transfereed = rows[5]

            if approved == "☑":
                temp_approved = True
            else:
                temp_approved = False
            all_compound_data[plates][well] = {"note": note,
                                               "approved": temp_approved,
                                               "transfereed": transfereed}
    return all_compound_data


def bio_data_approval_table(draw_plate, config, all_plates_data, assay_data, plate_to_layout, archive_plates_dict,
                            bio_plate_report_setup, transfer_dict, dismissed_plates,
                            all_plates_are_dismissed, dead_plates):
    """
    The controller for the pop-up that comes when you are importing plate-reader data to the database
    :param draw_plate: A function for drawing plates
    :type draw_plate: function
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param all_plates_data: The raw and calculated data for all the plates
    :type all_plates_data: dict
    :param plate_to_layout: a dicts for the plate_layout
    :type plate_to_layout: dict
    :param archive_plates_dict: the dict over the layouys
    :type archive_plates_dict: dict
    :param bio_plate_report_setup: dict over what state wells should be in, to be printed on the report sheet.
    :type bio_plate_report_setup: dict
    :param transfer_dict: data over wells that have been transferred and wells that have been skipped for each plate
    :type transfer_dict: dict
    :param dismissed_plates: A list of plates that have been dismissed
    :type dismissed_plates: list
    :param all_plates_are_dismissed: An idicator to see if all plates are dissed or not. If True, skips approval window
    :type all_plates_are_dismissed: bool
    :param dead_plates: A list of plates with no data
    :type dead_plates: list
    :return: The layout for the popup
    """
    # Initialize Bio analyse tools:
    bioa = BIOAnalyser(config, bio_plate_report_setup)

    # ToDO add the option to click colours and change them in the different mappings!
    # Makes a copy of the dict as it will be changed doing the analyses of the data
    temp_archive_plates_dict = copy.deepcopy(archive_plates_dict)
    table_row = None  # Make sure that it does not try to calculate stuff on empty data
    plate_analyse_methods = []  # What method are used for analysing plates
    plate_calculations = []  # What calculations have been done for each plate
    plate_table_data = []  # The data for all the plates
    compound_table_data = []  # The data for all the wells in the plates
    hist_data = {"all": []}  # The data for drawing the histogram
    well_state_overview = {  # A dict over each state if it is included in the analysis or not
        "sample": False,
        "blank": False,
        "max": False,
        "minimum": False,
        "pos": False,
        "neg": False,
        "empty": False,
    }
    status_checker = []

    # Sets stuff for drawing the histogram and plates
    histogram_canvas = None
    toolbar = None
    current_plate_name = None  # Make sure to only update the histogram if a new plate is selected
    color_select = {}
    for keys in list(config["plate_colouring"].keys()):
        color_select[keys] = config["plate_colouring"][keys]

    # Stuff for the tables
    BLANK_BOX = "☐"
    CHECKED_BOX = "☑"

    # sub_layout ditch
    sub_layouts = {}

    # Get z_prime_threshold:
    z_prime_threshold = assay_data["z_prime_threshold"]
    plate_size = assay_data["plate_size"]

    draw_options = ["state_mapping", "heatmap", "hit_map"]

    # Not sure if there needs to be other methods, but for now there is one
    analyse_methods = ["Histogram"]
    compound_data_dict = {}
    colour_row_index = []
    # Grabs data for the table
    for plate_index, plate in enumerate(all_plates_data):
        temp_plate_name = plate
        compound_data_dict[plate] = []
        sub_layouts[plate] = {"new": False,
                              "name": None,
                              "layout": None}
        run_name = all_plates_data[plate]["run_name"]
        try:
            temp_z_prime = round(float(all_plates_data[plate]["calculations"]["other"]["z_prime"]), 2)
        except KeyError:
            temp_z_prime = 0

        if temp_z_prime > z_prime_threshold:
            temp_approval = CHECKED_BOX
        else:
            temp_approval = BLANK_BOX
        if plate not in dead_plates:
            if transfer_dict[plate]["skipped_wells"]:
                temp_layout = "Changed"
                temp_approval = BLANK_BOX
            else:
                temp_layout = plate_to_layout[temp_plate_name]

        try:
            all_plates_data[plate]["plates"]
        except KeyError:
            for trans_skip_status in transfer_dict[plate]:
                for well in transfer_dict[plate][trans_skip_status]:
                    if trans_skip_status == "skipped_wells":
                        well_score = 0
                        note = transfer_dict[plate]["skipped_wells"][well]["reason"]
                        well_approval = BLANK_BOX
                        transferred = False
                        print(f"no plate {well} - {plate}")
                        try:
                            all_plates_data[plate]["skipped_wells"]
                        except KeyError:
                            all_plates_data[plate]["skipped_wells"] = [well]
                        else:
                            all_plates_data[plate]["skipped_wells"].append(well)
                    else:
                        well_score = 0
                        note = "No transfer data"
                        well_approval = temp_approval
                        transferred = True

                    if plate not in dismissed_plates:
                        temp_well_data = [plate, well, well_score, note, well_approval, transferred]
                    else:
                        if well in dismissed_plates[plate]:
                            transferred = True
                        else:
                            transferred = False
                        temp_well_data = [plate, well, "Na", "Dismissed", BLANK_BOX, transferred]
                    compound_data_dict[plate].append(temp_well_data)
                    compound_table_data.append(temp_well_data)

        else:
            for analysed_methods in all_plates_data[plate]["plates"]:
                # Finds the way the plate have been analysed
                if analysed_methods not in plate_analyse_methods:
                    plate_analyse_methods.append(analysed_methods)

                for status in all_plates_data[plate]["plates"][analysed_methods]:
                    if status == "wells":
                        for well_index, well in enumerate(all_plates_data[plate]["plates"][analysed_methods][status]):
                            try:
                                transfer_dict[plate]["skipped_wells"][well]
                            except KeyError:
                                well_score = round(float(all_plates_data[plate]["plates"][analysed_methods]["wells"]
                                                         [well]), 2)
                                note = ""
                                well_approval = temp_approval
                                transferred = True
                            else:
                                print(f"Plate {well} - {plate}")
                                well_score = 0
                                try:
                                    transfer_dict[plate]["skipped_wells"][well]["reason"]
                                except KeyError:
                                    note = "Missing reason"
                                else:
                                    note = transfer_dict[plate]["skipped_wells"][well]["reason"]
                                well_approval = BLANK_BOX
                                transferred = False
                                try:
                                    all_plates_data[plate]["skipped_wells"]
                                except KeyError:
                                    all_plates_data[plate]["skipped_wells"] = [well]
                                else:
                                    all_plates_data[plate]["skipped_wells"].append(well)
                            if plate not in dismissed_plates:
                                temp_well_data = [plate, well, well_score, note, well_approval, transferred]
                            else:
                                if well in dismissed_plates[plate]:
                                    transferred = True
                                else:
                                    transferred = False
                                temp_well_data = [plate, well, "Na", "Dismissed", BLANK_BOX, transferred]
                            compound_data_dict[plate].append(temp_well_data)
                            compound_table_data.append(temp_well_data)
                    elif status not in status_checker:
                        status_checker.append(status)
                        well_state_overview[status] = True

        if plate_index == 0:
            for calc_index, calculations in enumerate(all_plates_data[plate]["calculations"]):
                if calc_index == 0:
                    # Finds the different calculations that have been done on the plate
                    if calculations != "other":
                        for temp_counter, analysed_method in enumerate(
                                all_plates_data[plate]["calculations"][calculations]):
                            if analysed_method not in calculations and analysed_method not in plate_calculations:
                                if analysed_method != "other":
                                    plate_calculations.append(analysed_method)

        if plate in dismissed_plates:
            temp_row_data = [temp_plate_name, "Na", "Na", "Dismissed", run_name, BLANK_BOX]
            colour_row_index.append(plate_index)
        elif plate in dead_plates:
            temp_row_data = [temp_plate_name, "Dead Plate", "Dead Plate", "Dismissed", run_name, BLANK_BOX]
        else:
            temp_row_data = [temp_plate_name, temp_z_prime, temp_layout, "", run_name, temp_approval]

        plate_table_data.append(temp_row_data)

    # Gets the data for each well that is marked as a sample on the last analysed method, for each and all plates
    for plate in all_plates_data:
        if plate not in dead_plates:
            hist_data[plate] = []
            for well in temp_archive_plates_dict[plate_to_layout[plate]]["sample"]:
                data = all_plates_data[plate]["plates"][plate_analyse_methods[-1]]["wells"][well]
                hist_data[plate].append(data)
                hist_data["all"].append(data)

    # Headlines for the table
    plate_table_headings = ["Plate", "Z-prime", "Layout", "Notes", "run_name", "Approved"]
    compound_table_headings = ["Plate", "Well", "Score", "Notes", "Approved", "transferred"]

    if all_plates_are_dismissed:
        all_plates_data = _sort_plate_data(plate_table_data, all_plates_data, sub_layouts)
        all_compound_data = _sort_compound_data(compound_data_dict)
        return all_plates_data, all_compound_data, plate_analyse_methods, sub_layouts
    else:

        # Gets the layout
        window, plate_table_data = _bio_data_approval_table_layout(config, plate_table_data, plate_table_headings,
                                                                   compound_table_data, compound_table_headings,
                                                                   plate_analyse_methods, plate_calculations,
                                                                   analyse_methods, draw_options, well_state_overview)

        # Graph setup
        graph_plate = window["-BIO_APPROVAL_TABLE_GRAPH-"]
        well_dict = {}
        current_draw_tool = "sample"
        draw_tool_chooser = {"-BIO_APPROVAL_TABLE_SAMPLE_DRAW-": "sample",
                             "-BIO_APPROVAL_TABLE_BLANK_DRAW-": "blank",
                             "-BIO_APPROVAL_TABLE_MAX_DRAW-": "max",
                             "-BIO_APPROVAL_TABLE_MIN_DRAW-": "minimum",
                             "-BIO_APPROVAL_TABLE_POSITIVE_DRAW-": "positive",
                             "-BIO_APPROVAL_TABLE_NEGATIVE_DRAW-": "negative",
                             "-BIO_APPROVAL_TABLE_EMPTY_DRAW-": "empty"}
        dragging = False
        prior_selector = None
        selection_active = False
        start_point, end_point, x, y = None, None, None, None
        temp_well_selected = {}
        temp_plate_dict = {}
        plate_type_draw_counter = {"plate_96": 96,
                                   "plate_384": 384,
                                   "plate_1536": 1536}

        # Makes it possible to double-click on the table
        window["-BIO_APPROVAL_PLATE_TABLE-"].bind('<Double-Button-1>', "+-double click-")

        graph = window["-BIO_APPROVAL_TABLE_GRAPH-"]
        plate_archive_draw = True
        graph_placement = "bio_approval_popup"

        # Draw the histogram for all the data:
        histogram_canvas, toolbar = _draw_histogram(window, hist_data["all"], histogram_canvas, toolbar)

        # Initialize the Database-functions
        dbf = DataBaseFunctions(config)

        while True:
            event, values = window.read()

            if event == sg.WIN_CLOSED or event == "-WINDOW_TWO_CANCEL-" or event == "-POP_SAMPLE_CHECKER_OK-":
                # Grabs all the data from the table

                window.close()
                return "cancelled", plate_analyse_methods

            # ToDo make it possible to de-select wells on the plate-layout and then re-calculate everything with this
            #  layout - Make sure that the layout for the "plate-data" and the plate_layout are both changed
            #  for the specific plates
            # Re-calculates the current plate
            if event == "-BIO_APPROVAL_TABLE_RE_CALCULATE-" and current_plate_name:
                well_type = _get_well_type(well_dict, temp_well_selected)
                if well_type is not None:
                    new_plate_calc, _ = bioa.data_converter(temp_plate_dict, well_type)

                    # Get calculations:
                    _update_plate_calculations(window, values, temp_plate_dict, temp_analysed_method)

                else:
                    print("Missing well marked as something specific...")

                # Draw a new histogram with the new data
                if values["-BIO_APPROVAL_TABLE_ANALYSE_METHODS-"] == "Histogram":
                    temp_hist_data = _get_hist_data(temp_plate_dict, plate_analyse_methods, well_type)
                    histogram_canvas, toolbar = _draw_histogram(window, temp_hist_data, histogram_canvas, toolbar)

                # Re-draw the plate map with the newly selected well-states
                skipped_wells = all_plates_data[current_plate_name]["skipped_wells"]
                _draw_plate_on_graph(draw_plate, config, temp_plate_dict["plates"],
                                     temp_analysed_method, graph, plate_size, graph_placement, show_state_list,
                                     draw_options, skipped_wells)

            if event == "-BIO_APPROVAL_TABLE_APPLY-" and current_plate_name:

                all_plates_data[current_plate_name] = copy.deepcopy(temp_plate_dict)

                temp_z_prime = temp_plate_dict["calculations"]["other"]["z_prime"]
                new_plate_layout = plate_to_layout[current_plate_name]
                plate_layout_check = "changed"

                try:
                    all_plates_data[current_plate_name]["note"]
                except KeyError:
                    temp_notes = ""
                else:
                    temp_notes = "Note"

                temp_run_name = all_plates_data[current_plate_name]["run_name"]
                if temp_z_prime > z_prime_threshold:
                    temp_approval = CHECKED_BOX
                else:
                    temp_approval = BLANK_BOX

                temp_plate_dict.clear()
                temp_plate_dict, well_dict, temp_archive_plates_dict, current_layout_name = \
                    _plate_layout_controller(dbf, temp_archive_plates_dict, plate_layout_check,
                                             current_plate_name, plate_to_layout, all_plates_data)

                new_sub_name = plate_to_layout[current_plate_name]
                sub_layouts[current_plate_name] = {"new": True,
                                                   "name": new_sub_name,
                                                   "layout": temp_archive_plates_dict[new_sub_name]}

                plate_table_data = _update_plate_table(plate_table_data, table_row, temp_z_prime, new_plate_layout,
                                                       temp_notes, temp_run_name, temp_approval)

                new_compound_table_data = _update_compound_table(compound_data_dict, temp_archive_plates_dict,
                                                                 current_plate_name, current_layout_name)

                window["-BIO_APPROVAL_PLATE_TABLE-"].update(values=plate_table_data)
                window["-BIO_APPROVAL_COMPOUND_TABLE-"].update(values=new_compound_table_data)
            # Sends the data to the database
            if event == "-BIO_DATA_APPROVED-":
                temp_plate_errors = []
                for rows in plate_table_data:
                    if rows[5] == "☐" and rows[3] == "":
                        temp_plate_errors.append(rows[0])

                if temp_plate_errors:
                    sg.PopupError(f"The following plates are missing notes:\n\n"
                                  f"{temp_plate_errors}\n\n"
                                  f"You need to add notes to them, before they can be added to the database")

                else:

                    all_plates_data = _sort_plate_data(plate_table_data, all_plates_data, sub_layouts)
                    all_compound_data = _sort_compound_data(compound_data_dict)
                    window.close()

                    return all_plates_data, all_compound_data, plate_analyse_methods, sub_layouts

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
                        all_plates_data[selected_plate]["note"] = temp_note

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
                        for row_index, plates in enumerate(all_plates_data):
                            plate_table_data[row_index][3] = "Note"
                            all_plates_data[plates]["note"] = temp_note

                        window["-BIO_APPROVAL_PLATE_TABLE-"].update(values=plate_table_data)
                        window["-BIO_APPROVAL_PLATE_NOTE-"].update(value=temp_note)

            # Update the calculations in the window, when you change the calculation in the dropdown menu
            if event == "-BIO_APPROVAL_TABLE_DROPDOWN_CALCULATIONS-" and table_row is not None:
                temp_analysed_method = values["-BIO_APPROVAL_TABLE_DROPDOWN_PLATE_ANALYSE_METHODS-"]

                states = values["-BIO_APPROVAL_TABLE_DROPDOWN_CALCULATIONS-"]
                all_calc = all_plates_data[current_plate_name]["calculations"][temp_analysed_method][states]
                window["-BIO_APPROVAL_TABLE_CALC_AVG-"].update(value=all_calc["avg"])
                window["-BIO_APPROVAL_TABLE_CALC_STDEV-"].update(value=all_calc["stdev"])
                window["-BIO_APPROVAL_TABLE_CALC_PSTDEV-"].update(value=all_calc["pstdev"])
                window["-BIO_APPROVAL_TABLE_CALC_PVARIANCE-"].update(value=all_calc["pvariance"])
                window["-BIO_APPROVAL_TABLE_CALC_VARIANCE-"].update(value=all_calc["variance"])
                window["-BIO_APPROVAL_TABLE_CALC_ST_DEV-"].update(value=all_calc["st_dev_%"])

            # Makes it possible to check and un-check a checkbox.
            if event[0] == "-BIO_APPROVAL_PLATE_TABLE-" and event[2][1] == 5:
                row_data = event[2][0]
                if plate_table_data[row_data][3].casefold() == "dismissed":
                    sg.popup_error("This plate have been dismissed, and can't be approved")
                else:
                    if plate_table_data[row_data][5] == CHECKED_BOX:

                        plate_table_data[row_data][5] = BLANK_BOX
                    else:
                        plate_table_data[row_data][5] = CHECKED_BOX

                    window['-BIO_APPROVAL_PLATE_TABLE-'].update(values=plate_table_data)

            # updates the window, when double-clicking the plate table, or changing what data to look at
            if event == "-BIO_APPROVAL_PLATE_TABLE-+-double click-" or event == "-BIO_APPROVAL_TABLE_RESET-":
                try:
                    table_row = values["-BIO_APPROVAL_PLATE_TABLE-"][0]
                except IndexError:
                    pass
                else:
                    # Grab data
                    current_plate_name = plate_table_data[table_row][0]
                    if current_plate_name in dead_plates:
                        sg.PopupError("No Plate Data on selected plate")
                    else:

                        plate_layout_check = plate_table_data[table_row][2]
                        temp_analysed_method = values["-BIO_APPROVAL_TABLE_DROPDOWN_PLATE_ANALYSE_METHODS-"]

                        # Generates a dict over the wells on the graph for drawing
                        well_dict.clear()
                        temp_well_selected.clear()
                        temp_plate_dict.clear()

                        temp_plate_dict, well_dict, temp_archive_plates_dict, current_layout_name = \
                            _plate_layout_controller(dbf, temp_archive_plates_dict, plate_layout_check,
                                                     current_plate_name, plate_to_layout,
                                                     all_plates_data)

                        # Gets what well states should be included in the mapping.
                        show_state_list = {"sample": values["-BIO_APPROVAL_TABLE_SAMPLE-"],
                                           "blank": values["-BIO_APPROVAL_TABLE_BLANK-"],
                                           "max": values["-BIO_APPROVAL_TABLE_MAX-"],
                                           "minimum": values["-BIO_APPROVAL_TABLE_MIN-"],
                                           "positive": values["-BIO_APPROVAL_TABLE_POSITIVE-"],
                                           "negative": values["-BIO_APPROVAL_TABLE_NEGATIVE-"],
                                           "empty": values["-BIO_APPROVAL_TABLE_EMPTY-"]}

                        # Draw histrogram
                        if values["-BIO_APPROVAL_TABLE_ANALYSE_METHODS-"] == "Histogram":
                            histogram_canvas, toolbar = _draw_histogram(window, hist_data[current_plate_name],
                                                                        histogram_canvas, toolbar)

                        # Updates the compound table, with compounds corresponding to the selected
                        # values for the plate mapping
                        new_compound_table = []
                        for states in show_state_list:
                            if show_state_list[states]:
                                for data in compound_data_dict[current_plate_name]:
                                    if data[1] in temp_archive_plates_dict[plate_to_layout[current_plate_name]][states]:
                                        new_compound_table.append(data)
                        window["-BIO_APPROVAL_COMPOUND_TABLE-"].update(values=new_compound_table)

                        # Get calculations:
                        draw_options = _update_plate_calculations(window, values, temp_plate_dict, temp_analysed_method)

                        if not any(show_state_list.values()):
                            sg.PopupError("Please select states to include in the graph")
                        else:
                            # Updates the plate map with selected values
                            skipped_wells = all_plates_data[current_plate_name]["skipped_wells"]
                            _, min_x, min_y, max_x, max_y = \
                                _draw_plate_on_graph(draw_plate, config, temp_plate_dict["plates"],
                                                     temp_analysed_method, graph, plate_size, graph_placement,
                                                     show_state_list, draw_options, skipped_wells)

            if event == "-BIO_APPROVAL_TABLE_SAMPLE-" or \
                    event == "-BIO_APPROVAL_TABLE_BLANK-" or event == "-BIO_APPROVAL_TABLE_MAX-" or \
                    event == "-BIO_APPROVAL_TABLE_MIN-" or event == "-BIO_APPROVAL_TABLE_POSITIVE-" or \
                    event == "-BIO_APPROVAL_TABLE_NEGATIVE-" or event == "-BIO_APPROVAL_TABLE_DROPDOWN_DRAW_OPTIONS-" or \
                    event == "-BIO_APPROVAL_TABLE_DROPDOWN_PLATE_ANALYSE_METHODS-":

                # Grab data
                try:
                    plate_table_data[table_row]
                except KeyError:
                    pass
                else:
                    temp_plate_name = plate_table_data[table_row][0]
                    temp_analysed_method = values["-BIO_APPROVAL_TABLE_DROPDOWN_PLATE_ANALYSE_METHODS-"]

                    # Gets what well states should be included in the mapping.
                    show_state_list = {"sample": values["-BIO_APPROVAL_TABLE_SAMPLE-"],
                                       "blank": values["-BIO_APPROVAL_TABLE_BLANK-"],
                                       "max": values["-BIO_APPROVAL_TABLE_MAX-"],
                                       "minimum": values["-BIO_APPROVAL_TABLE_MIN-"],
                                       "positive": values["-BIO_APPROVAL_TABLE_POSITIVE-"],
                                       "negative": values["-BIO_APPROVAL_TABLE_NEGATIVE-"],
                                       "empty": values["-BIO_APPROVAL_TABLE_EMPTY-"]}

                    # Get calculations:
                    draw_options = _update_plate_calculations(window, values, temp_plate_dict, temp_analysed_method)

                    if not any(show_state_list.values()):
                        sg.PopupError("Please select states to include in the graph")
                    else:
                        # Updates the plate map with selected values
                        skipped_wells = all_plates_data[temp_plate_name]["skipped_wells"]
                        _, min_x, min_y, max_x, max_y = \
                            _draw_plate_on_graph(draw_plate, config, temp_plate_dict["plates"],
                                                 temp_analysed_method, graph, plate_size, graph_placement,
                                                 show_state_list, draw_options, skipped_wells)

            # update drawing tool
            if event in draw_tool_chooser:
                last_draw_tool = current_draw_tool
                current_draw_tool = draw_tool_chooser[event]
                last_event = [_ for _ in draw_tool_chooser if draw_tool_chooser[_] == last_draw_tool]
                window[last_event[0]].update(value=False)
                window[event].update(value=True)

            # Canvas controls
            try:
                event.endswith("+MOVE")
            except AttributeError:
                pass
            else:
                if event.endswith("+MOVE"):
                    if values["-BIO_APPROVAL_TABLE_GRAPH-"][0] and values["-BIO_APPROVAL_TABLE_GRAPH-"][1]:
                        try:
                            temp_well = graph_plate.get_figures_at_location(values['-BIO_APPROVAL_TABLE_GRAPH-'])[0]
                            plate_type = temp_archive_plates_dict[plate_to_layout[current_plate_name]]["plate_type"]
                            temp_well = temp_well % int(plate_type_draw_counter[plate_type])
                            if temp_well == 0:
                                temp_well = str(plate_type_draw_counter[plate_type])
                            else:
                                temp_well = str(temp_well)
                            try:
                                well_dict[temp_well]["well_id"]
                            except KeyError:
                                temp_well = int(temp_well)

                            temp_well_id = (well_dict[temp_well]["well_id"])
                        except (IndexError or KeyError) as error:
                            # print(f"Canvas Error: {error}")
                            temp_well_id = ""

                        graph_coordinates = values['-BIO_APPROVAL_TABLE_GRAPH-']
                        graph_string = f"Well: {temp_well_id}"
                        score_string = ""
                        if values["-BIO_APPROVAL_TABLE_DROPDOWN_DRAW_OPTIONS-"] == "hit_map" and temp_well_id:
                            temp_well_score = \
                                round(float(all_plates_data[temp_plate_name]["plates"][plate_analyse_methods[-1]]
                                            ["wells"][temp_well_id]), 2)
                            score_string += f" Score: {temp_well_score}"
                        elif values["-BIO_APPROVAL_TABLE_DROPDOWN_DRAW_OPTIONS-"] == "heatmap":
                            graph_string += "Placeholder"
                        else:
                            # graph_string += f" Draw-Tool: {current_draw_tool}"
                            score_string = f" Draw-Tool: {current_draw_tool}"

                        window["-BIO_APPROVAL_WELL_ID-"].update(value=graph_string)
                        window["-BIO_APPROVAL_WELL_SCORE-"].update(value=score_string)

            if event == "-BIO_APPROVAL_TABLE_GRAPH-":
                x, y = values["-BIO_APPROVAL_TABLE_GRAPH-"]

                # Sets a start point if mouse button have not been clicked, and sets an end point if it is released
                if not dragging:
                    start_point = (x, y)
                    dragging = True
                else:
                    end_point = (x, y)

                # Delete prior selection if present
                if prior_selector:
                    graph_plate.delete_figure(prior_selector)

                if None not in (start_point, end_point):
                    selection_active = True
                    prior_selector = graph_plate.draw_rectangle(start_point, end_point, fill_color="", line_color="white")

            try:
                event.endswith("+UP")
            except AttributeError:
                pass
            else:
                if event.endswith("+UP"):
                    if selection_active and current_plate_name:

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

                        # makes a set, for adding wells, to avoid duplicates
                        graphs_list = set()

                        # goes over the coordinates and if they are within the bounds of the plate
                        # if that is the case, then the figure for that location is added the set
                        for index_x, cords_x in enumerate(temp_x):
                            for index_y, cords_y in enumerate(temp_y):
                                if min_x <= temp_x[index_x] <= max_x and min_y <= temp_y[index_y] <= max_y:
                                    graphs_list.add(
                                        graph_plate.get_figures_at_location((temp_x[index_x], temp_y[index_y]))[0])

                        # colours the wells in different colour, depending on their status
                        for well in graphs_list:
                            # ToDo Get this selection tool, to only change the colour of the border not the whole cell
                            # ToDO make it work for all plate size not only 384, as it is hard coded to atm.
                            colour = color_select[current_draw_tool]
                            well_state = current_draw_tool
                            # if colour == "paint":
                            #     colour = values["-PLATE_LAYOUT_COLOUR_CHOSE_TARGET-"]
                            graph_plate.Widget.itemconfig(well, fill=colour)
                            plate_type = temp_archive_plates_dict[plate_to_layout[current_plate_name]]["plate_type"]
                            well = well % plate_type_draw_counter[plate_type]
                            if well == 0:
                                well = str(plate_type_draw_counter[plate_type])
                            else:
                                well = str(well)
                            try:
                                (well_dict[well]["well_id"])
                            except KeyError:
                                well = well_dict[int(well)]["well_id"]
                            else:
                                well = (well_dict[well]["well_id"])
                            temp_well_selected[well] = {"state": well_state,
                                                        "colour": colour}
                            # well_dict[wells]["colour"] = colour
                            # well_dict[wells]["state"] = well_state

                    # deletes the rectangle used for selection
                    if prior_selector:
                        graph_plate.delete_figure(prior_selector)

                        # reset everything
                    start_point, end_point = None, None
                    dragging = False
                    prior_selector = None
                    temp_selector = False
                    temp_draw_tool = None


if __name__ == "__main__":
    import configparser
    config = configparser.ConfigParser()
    config.read("config.ini")

    dbf = DataBaseFunctions(config)
    current_plate_name = "alpha_so1"

    rows = dbf.find_data_single_lookup(table="assay_plates", data_value=current_plate_name, headline="plate_name")
    assay_run = rows[0][2]
    rows = dbf.find_data_single_lookup(table="assay_runs", data_value=assay_run, headline="run_name")

    worklist_echo_data = eval(rows[0][4])
