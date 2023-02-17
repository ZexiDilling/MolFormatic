import PySimpleGUI as sg

# from copy import deepcopy
# from info import matrix_header
# from gui_functions import sub_settings_matrix
# from database_controller import FetchData
# from excel_handler import purity_sample_layout_import, purity_sample_layout_export
# from gui_functions import sort_table


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
    raw_table_col = sg.Frame("Headlines", [[
        sg.Column([
            [sg.Table(values=table_data, headings=headings,
                      key="-POP_HEADLINE_TABLE-", enable_events=True, enable_click_events=True)]
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

    for headline_index, headline in enumerate(right_headlines):
        temp_data = [headline, wrong_headlines[headline_index]]
        table_data.append(temp_data)

    table_headings = ["Name in file", "Right Name"]

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
            print(new_headlines)
            for keys in new_headlines:
                if new_headlines[keys] in hash_val:
                    flag = True
                    break
                else:
                    hash_val[new_headlines[keys]] = 1
            if flag:
                sg.popup_error("There are repeated values, please try again")
            else:
                return new_headlines
                window.close()

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

if __name__ == "__main__":
    right_headlines = ["source_plates", "destination_plates", "source_well", "destination_well"]
    wrong_headlines = ["source_plates", "destination_plates", "source_well", "destination_well"]
    new_headlines_popup(right_headlines, wrong_headlines)