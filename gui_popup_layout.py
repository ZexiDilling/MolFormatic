import PySimpleGUI as sg

from info import unit_converter_list_mol, unit_converter_list_liquids, matrix_header


def morgan_popup_layout(config, main_values):
    sg.theme(config["GUI"]["theme"])

    sub_searchs_text = 10
    sub_search_input_fields = 15

    col = sg.Frame("Morgan Values", [[
        sg.Column([
            [sg.Text("Morgan specific options", key="-MORGAN_POPUP_OPTIONS-")],
            [sg.Checkbox(text="chirality", key="-MORGAN_POPUP_CHIRALITY-"),
             sg.Checkbox(text="Features", key="-MORGAN_POPUP_FEATURES-")],
            [sg.Text("n bits", key="-MORGAN_POPUP_BITS_TEXT-", size=sub_searchs_text),
             sg.InputText(key="-MORGAN_POPUP_BITS-", size=sub_search_input_fields)],
            [sg.Text("bound range", key="-MORGAN_POPUP_BOUND_TEXT-", size=sub_searchs_text),
             sg.InputText(key="-MORGAN_POPUP_RANGE-", size=sub_search_input_fields)],
            [sg.B("Apply", key="-MORGAN_POPUP_APPLY-"),
             sg.B("Cancel", key="-CLOSE_POPUP-")]
        ])
    ]])

    layout = [[col]]
    return sg.Window("Morgan Values", layout, finalize=True, resizable=True)


def popup_three_box_solution_layout(config, name, question, box_1, box_2):
    sg.theme(config["GUI"]["theme"])
    col = sg.Column([
            [sg.Text(question)],
            [sg.Button(box_1, key="-TABLE_POPUP_BOX_1-"),
             sg.Button(box_2, key="-TABLE_POPUP_BOX_2-"),
             sg.Button("Cancel", key="-TABLE_POPUP_CANCEL-")]
    ])

    layout = [[col]]
    return sg.Window(name, layout, finalize=True, resizable=True)


def table_popup_layout(config, table_name, table_headings, table_data):
    sg.theme(config["GUI"]["theme"])
    col = sg.Frame(table_name, [[
        sg.Column([
            [sg.Table(table_data, headings=table_headings, auto_size_columns=True, key="-POPUP_TABLE-")],
            [sg.Button("Close", key="-TABLE_POPUP_DONE-")]
        ])
    ]])

    layout = [[col]]
    return sg.Window("Table overview", layout, finalize=True, resizable=True)


def matrix_popup_layout(config, calc, state=None, method=None):
    sg.theme(config["GUI"]["theme"])
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


def sample_checker_layout(config, table_data, headings, db_data):
    sg.theme(config["GUI"]["theme"])
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


def new_headlines_layout(config, table_data, headings):
    sg.theme(config["GUI"]["theme"])
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


def plate_layout_chooser_layout(config, table_data, headings):
    sg.theme(config["GUI"]["theme"])
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


def assay_generator_layout(config, plate_layout):
    sg.theme(config["GUI"]["theme"])
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


def assay_run_naming_layout(config, plate_table_headline, run_name, previous_runs, assay_name, used_plates_table_data,
                            all_batch_numbers, batch_number):
    sg.theme(config["GUI"]["theme"])
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


def dead_run_naming_layout(config, plate_table_headline, run_name, previous_runs, assay_name, plates_table_data,
                           all_batch_numbers, batch_number, worklist):
    sg.theme(config["GUI"]["theme"])
    text_size = 15
    button_size = 15
    input_size = 15
    if worklist:
        temp_worklist = worklist
    else:
        temp_worklist = "No Worklist"

    previous_runs.append("New")
    top_layout = sg.Column([
        [sg.T("Run Name", size=text_size),
         sg.Input(run_name, key="-ASSAY_RUN_NAME-", size=input_size)],
        [sg.T("Chose old Run", size=text_size),
         sg.DropDown(previous_runs, key="-ASSAY_RUN_PREVIOUS-",
                     size=input_size, enable_events=True, default_value=previous_runs[-1])],
        [sg.B("New batch", key="-ASSAY_RUN_NEW_BATCH-", size=button_size),
         sg.T(batch_number, key="-ASSAY_RUN_CURRENT_BATCH-", size=text_size, visible=False),
         sg.DropDown(values=all_batch_numbers, key="-ASSAY_RUN_ALL_BATCHES-", size=text_size,
                     default_value=batch_number, enable_events=True,
                     tooltip="A way to split up data, "
                             "if there are changes to the data that are specific to some run.")],
        [sg.CalendarButton("Start date", key="-ASSAY_RUN_DATE-", format="%Y-%m-%d", enable_events=True,
                           target="-ASSAY_RUN_DATE_TARGET-", size=(button_size, 1)),
         sg.Input(key="-ASSAY_RUN_DATE_TARGET-", enable_events=True, size=text_size)],
        [sg.B("Worklist", key="-ASSAY_RUN_WORKLIST-", size=button_size),
         sg.T(temp_worklist, key="-ASSAY_RUN_WORKLIST_INDICATOR-", relief="sunken", size=text_size),
         sg.Multiline(visible=False, key="-ASSAY_RUN_WORKLIST_DATA-")],
        [sg.B("Echo Data", key="-ASSAY_RUN_ECHO-", size=button_size),
         sg.T("No Echo Data", key="-ASSAY_RUN_ECHO_INDICATOR-", relief="sunken", size=text_size),
         sg.Multiline(visible=False, key="-ASSAY_RUN_ECHO_DATA-")],

        [sg.HorizontalSeparator()]
    ])
    bottom_layer = sg.vtop([
        sg.Column([
            [sg.T("Notes")],
            [sg.Multiline(default_text="Dead Plates - ", key="-ASSAY_RUN_NOTES-", size=(18, 15))],
            [sg.Button("Update Run", key="-ASSAY_RUN_UPDATE-", size=text_size)]
        ]),
        sg.Column([
            [sg.T("All Plates")],
            [sg.Table(values=plates_table_data, size=(18, 25), headings=plate_table_headline,
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
         # sg.Button("Dismiss Runs", key="-ASSAY_RUN_DISMISS-", expand_x=text_size,
         #           tooltip="Will ask for assay_run name for witch plates to dismiss,"
         #                   "if there is more than one run name"
         #                   "Notes needs to be filled in for the assay_run"),
         sg.Button("Cancel", key="-WINDOW_TWO_CANCEL-", expand_x=text_size)]
    ]

    return sg.Window("Assay Run", layout, finalize=True, resizable=True)


def bio_data_approval_table_layout(config, plate_table_data, plate_headings, compound_table_data, compound_headings,
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
    :param draw_options: The different option for drawing the plate
    :type draw_options: list
    :param well_state_overview: The states included in the plate-map. It is used to disabling states that are not used
    :type well_state_overview: dict
    :return: The layout for the popup
    """
    sg.theme(config["GUI"]["theme"])
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
                 sg.T("Score:", key="-BIO_APPROVAL_WELL_SCORE-", size=14)],
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


def bio_dose_response_set_up_layout(config, worklist, assay_name, run_name, dose_response_calc, used_plates_table_data,
                                    plate_table_headline, calc_methodes):
    sg.theme(config["GUI"]["theme"])
    text_size = 15
    button_size = 15
    input_size = 15
    if worklist:
        temp_worklist = worklist
    else:
        temp_worklist = "No Worklist"

    overview_headings = ["Conc.R", "Conc.T", "vol (nL)", "% solvent", "Stock", "D-Fold"]
    stock_headings = ["Conc.", "D-Fold"]

    top_layout = sg.vtop([
        sg.Column([
            [sg.T("Run Name", size=text_size),
             sg.Input(run_name, key="-DOSE_RESPONSE_RUN_NAME-", size=input_size)],
            [sg.T("Curve Shape", size=text_size, tooltip="Z-shape is for active assays and S-Shape is for inactive"),
             sg.Radio("Z - TODO", group_id=1, key="-DOSE_RESPONSE_CURVE_SHAPE_Z-",
                      enable_events=True, default=True),
             sg.Radio("S - TODO", group_id=1, key="-DOSE_RESPONSE_CURVE_SHAPE_S-",
                      enable_events=True),
             # TODO check if this is correct
             sg.Input("Z", key="-DOSE_RESPONSE_CURVE_SHAPE-", visible=False)],
            [sg.T("Calc Method:", size=text_size,
                  tooltip="This is the method use for analysing the curve generated"),
             sg.DropDown(calc_methodes, key="-CALC_METHOD-", size=input_size, enable_events=True)],
            [sg.B("New Method", key="-NEW_CALC_METHOD-", size=button_size),
             sg.T("", key="-DOSE_RESPONSE_CALC_METHOD-")],
            [sg.T("Dilution setup:", size=text_size,
                  tooltip="This is the setup for the dilution of the samples, concentration and so."),
             sg.DropDown(dose_response_calc, key="-DOSE_RESPONSE_CALC-", size=input_size, enable_events=True)],
            [sg.T("", key="-DOSE_RESPONSE_SETUP-", visible=True)],
            [sg.CalendarButton("Start date", key="-DOSE_RESPONSE_DATE-", format="%Y-%m-%d", enable_events=True,
                               target="-DOSE_RESPONSE_DATE_TARGET-", size=(button_size, 1)),
             sg.Input(key="-DOSE_RESPONSE_DATE_TARGET-", enable_events=True, size=text_size)],
            [sg.B("Worklist", key="-DOSE_RESPONSE_WORKLIST-", size=button_size),
             sg.T(temp_worklist, key="-DOSE_RESPONSE_WORKLIST_INDICATOR-", relief="sunken", size=text_size),
             sg.Multiline(visible=False, key="-DOSE_RESPONSE_WORKLIST_DATA-")],
            [sg.B("Echo Data", key="-DOSE_RESPONSE_ECHO-", size=button_size),
             sg.T("No Echo Data", key="-DOSE_RESPONSE_ECHO_INDICATOR-", relief="sunken", size=text_size),
             sg.Multiline(visible=False, key="-DOSE_RESPONSE_ECHO_DATA-")],
            [sg.HorizontalSeparator()]
        ]),
    ])

    bottom_layer = sg.vtop([
        sg.Column([
            [sg.T("Notes")],
            [sg.Multiline(default_text="Dead Plates - ", key="-DOSE_RESPONSE_NOTES-", size=(18, 15))],
            [sg.Button("Update Run", key="-ASSAY_RUN_UPDATE-", size=text_size)],
            [sg.T("Current run-notes")],
            [sg.T("", key="-DOSE_RESPONSE_CURRENT_NOTES-", relief="sunken", size=18),
             sg.InputText("", key="-DOSE_RESPONSE_CURRENT_NOTES_TARGET-", visible=False)],
            [sg.B("Show full text", key="-DOSE_RESPONSE_SHOW_NOTES-", size=button_size)]
        ]),
        sg.Column([
            [sg.T("All Plates")],
            [sg.Table(values=used_plates_table_data, size=(18, 25), headings=plate_table_headline,
                      key="-DOSE_RESPONSE_USED_PLATES_TABLE-")],
            [sg.B("Selected", key="-DOSE_RESPONSE_APPLY_SELECTED-", size=10,
                  tooltip="Apply the Run name to the selected plates"),
             sg.B("All", key="-DOSE_RESPONSE_APPLY_ALL-", size=7,
                  tooltip="Apply the Run name to all plates")]
        ])
    ])

    calc_layout = sg.Column([
        [sg.T("Stock:", size=text_size),
         sg.Input("", key="-CALC_DOSE_STOCK-", size=text_size),
         sg.DropDown(values=unit_converter_list_mol, key="-CALC_DOSE_STOCK_UNIT-",
                     default_value="mM")],
        [sg.T("Stock Dilution:", size=text_size),
         sg.Input("", key="-CALC_DOSE_STOCK_DILUTION-", size=text_size)],
        [sg.T("Max % solvent conc.:", size=text_size),
         sg.Input("", key="-CALC_MAX_SOLVENT_CONCENTRATION-", size=text_size)],
        [sg.T("Max Concentration:", size=text_size),
         sg.Input("", key="-CALC_DOSE_MAX_CONC-", size=text_size),
         sg.DropDown(values=unit_converter_list_mol, key="-CALC_DOSE_MAX_CONC_UNIT-",
                     default_value="mM")],
        [sg.T("Min Concentration:", size=text_size),
         sg.Input("", key="-CALC_DOSE_MIN_CONC-", size=text_size),
         sg.DropDown(values=unit_converter_list_mol, key="-CALC_DOSE_MIN_CONC_UNIT-",
                     default_value="mM")],
        [sg.T("Final Volume:", size=text_size),
         sg.Input("", key="-CALC_DOSE_FINAL_VOL-", size=text_size),
         sg.DropDown(values=unit_converter_list_liquids, key="-CALC_DOSE_FINAL_VOL_UNIT-",
                     default_value="uL")],
        [sg.T("Min trans Volume:", size=text_size),
         sg.Input(2.5, key="-CALC_DOSE_MIN_TRANS_VOL-", size=text_size),
         sg.DropDown(values=unit_converter_list_liquids, key="-CALC_DOSE_MIN_TRANS_VOL_UNIT-",
                     default_value="nL")],
        [sg.T("Dilution factor:", size=text_size),
         sg.Input("", key="-CALC_DOSE_DILUTION_FACTOR-", size=text_size)],
        [sg.T("Dilution Steps:", size=text_size, visible=False),
         sg.Input("", key="-CALC_DOSE_DILUTION_STEPS-", size=text_size, visible=False)],
        [sg.Table(values=[], headings=overview_headings, key="-CALC_TABLE_OVERVIEW-",
                  justification="center", auto_size_columns=False, enable_click_events=True,
                  num_rows=15, visible=True)],
        [sg.Table(values=[], headings=stock_headings, key="-CALC_TABLE_STOCK-",
                  justification="center", auto_size_columns=False, enable_click_events=True,
                  num_rows=5, visible=True)],
        [sg.B("Calculate", key="-DOSE_RESPONSE_CALC_BUTTON-", size=button_size),
         sg.VPush(),
         sg.B("Add to DB", key="-DOSE_RESPONSE_CALC_ADD_TO_DB-", size=button_size,
              tooltip="Adds these calculations to the database"),
         sg.B("Clear", key="-DOSE_RESPONSE_CALC_CLEAR-", size=button_size)]
    ]),

    assay_layout = sg.Frame(f" Dose response for {assay_name}", [
        sg.vtop([
            sg.Column([
                top_layout,
                bottom_layer,
            ]),
            sg.Frame("Calculations", [
                [sg.Column([
                    calc_layout
                ])]
            ])
        ])
    ])

    layout = [
        [assay_layout],
        [sg.Button("Analyse", key="-DOSE_RESPONSE_SETUP_DONE-", expand_x=text_size),
         sg.Button("Cancel", key="-WINDOW_DOSE_CANCELLED-", expand_x=text_size)]
    ]

    return sg.Window("Dose Response", layout, finalize=True, resizable=True)


def export_chooser_popup_layout(config):
    sg.theme(config["GUI"]["theme"])
    layout = [
        [sg.Checkbox("Excel", key="-EXPORT_EXCEL-"),
         sg.Checkbox("CSV", key="-EXPORT_CSV-")],
        [sg.B("OK", key="-EXPORT_OK-"),
         sg.B("Cancel", key="-EXPORT_CANCEL")]
    ]

    return sg.Window("Export", layout, finalize=True, resizable=False)


def bio_dose_response_approval_layout(config, plate_table_data, plate_headings,
                                      compound_overview_headings, compound_dose_headings,
                                      method, dose_units, analyse_methods, calc_method):
    sg.theme(config["GUI"]["theme"])
    raw_table_col = sg.Frame("Please approve or dismiss plates", [[
        sg.Column([
            [sg.Table(values=plate_table_data, headings=plate_headings, auto_size_columns=False,
                      key="-DOSE_APPROVAL_PLATE_TABLE-", enable_events=True, enable_click_events=True,
                      justification="center")],
            [sg.Table(values=[], headings=compound_overview_headings, auto_size_columns=False,
                      key="-DOSE_APPROVAL_COMPOUND_OVERVIEW_TABLE-", enable_events=True, enable_click_events=True,
                      justification="center")],
            [sg.Table(values=[], headings=compound_dose_headings, auto_size_columns=False,
                      key="-DOSE_APPROVAL_DOSE_COMPOUND_TABLE-", enable_events=True, enable_click_events=True,
                      justification="center")],

            [sg.T("TO DO LIST:")],
            [sg.T("GET difference between each point and the intended line")],
            [sg.T("Make a difference avg")],
            [sg.T("Make it possible to remove outliers")],
            [sg.T("Mark samples as Hits")]
        ])
    ]])

    text_size = 15
    input_size = 7
    if type(dose_units) == list:
        dose_units = dose_units[0]

    dose_data_info = sg.Frame("Basic", [[
        sg.Column([
            [sg.T(f"{method} ({dose_units})", key="-DOSE_APPROVAL_METHOD_UNITS-", size=text_size),
             sg.Input("", key="-DOSE_APPROVAL_METHOD_UNIT-", size=input_size)],
            [sg.T("Rsquared", size=text_size),
             sg.Input("", key="-DOSE_APPROVAL_R_SQUARED-", size=input_size)],
            [sg.T("Hillslope", size=text_size),
             sg.Input("", key="-DOSE_APPROVAL_HILLSLOPE-", size=input_size)],
            [sg.T("n_low dose", size=text_size),
             sg.Input("", key="-DOSE_APPROVAL_N_LOW_DOSE_DATA-", size=input_size)],
            [sg.T("Std low dose", size=text_size),
             sg.Input("", key="-DOSE_APPROVAL_STD_LOW_DOSE_DATA-", size=input_size)],
            [sg.T("n high dose", size=text_size),
             sg.Input("", key="-DOSE_APPROVAL_N_HIGH_DOSE_DATA-", size=input_size)],
            [sg.T("Std high dose", size=text_size),
             sg.Input("", key="-DOSE_APPROVAL_STD_HIGH_DOSE_DATA-", size=input_size)],
            [sg.T("Dose conc. stepsize", size=text_size),
             sg.Input("", key="-DOSE_APPROVAL_DOSE_CONC-", size=input_size)],
            [sg.T("Slope at lowdose", size=text_size),
             sg.Input("", key="-DOSE_APPROVAL_SLOP_LOW_DOSE-", size=input_size)],
            [sg.T("Slope at highdose", size=text_size),
             sg.Input("", key="-DOSE_APPROVAL_SLOP_HIGH_DOSE-", size=input_size)],
            [sg.T("Sample", size=text_size),
             sg.Input("", key="-DOSE_APPROVAL_SAMPLE-", size=input_size)],
        ])
    ]])

    canvas_window = sg.Frame("Data", [[
        sg.Column([
            [sg.Canvas(key="-DOSE_APPROVAL_TOOLBAR-")],
            [sg.Canvas(key="-DOSE_APPROVAL_CANVAS-", background_color="Grey", size=(500, 300))],
            [sg.DropDown(values=analyse_methods, key="-DOSE_APPROVAL_ANALYSE_METHODS-", enable_events=True,
                         default_value=analyse_methods[0]),
             sg.DropDown(values=calc_method, key="-DOSE_APPROVAL_CALC_METHODS-", enable_events=True,
                         default_value=calc_method[0]),
             sg.Checkbox("Include Control", key="-DOSE_APPROVAL_INCLUDE_CONTROL-", default=True)]
        ])
    ]])

    layout = [sg.vtop([raw_table_col, dose_data_info, canvas_window]),
              [sg.Button("Done", key="-DOSE_DATA_APPROVED-", expand_x=True),
               sg.B("DOSE", key="-DOSE_BUTTON-"),
               sg.Button("Cancel", key="-WINDOW_TWO_CANCEL-", expand_x=True)]
              ]

    return sg.Window("Dose Response", layout, finalize=True, resizable=True)





