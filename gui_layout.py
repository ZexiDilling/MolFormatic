import PySimpleGUI as sg
from info import matrix_header

#ToDo add tooltips to everything!!!!!!!!! ARG!!!!!!

class GUILayout:
    def __init__(self, config, plate_list):
        self.config = config
        self.standard_size = 20
        self.button_height = 1
        self.tab_colour = config["GUI"]["tab_colour"]
        self.plate_list = plate_list
        self.sample_style = ["Single Point", "Duplicate", "Triplicate", "Custom"]
        self.analyse_style = ["Single", "Dose Response"]

    @staticmethod
    def menu_top():
        """
        :return: The layout for the top menu
        :rtype: list
        """
        menu_top_def = [
            ["&File", ["&Open    Ctrl-O", "&Save    Ctrl-S", "---", '&Properties',  "&Exit", ]],
            ["&Edit", ["Paste", ["Special", "Normal", ], "Undo"], ],
            ["&Help", ["Info", "About..."]],
            ]
        layout = [[sg.Menu(menu_top_def)]]
        return layout

    @staticmethod
    def menu_mouse():
        """

        :return: The mouse menu
        :rtype: list
        """
        menu_mouse_def = [['File', ['Open', 'Save', 'Exit', ]],
                    ['Edit', ['Paste', ['Special', 'Normal', ], 'Undo'], ],
                    ['Help', 'About...'], ]

        layout = [[sg.ButtonMenu("place_holder", menu_mouse_def, key="-RMB-")]]

        return layout, menu_mouse_def

    def setup_1_search(self):
        """

        :return: A layour for the search-module in the top box
        :rtype: list
        """
        ac = ["Academic", "Commercial"]
        # origin = [self.config["database_specific_commercial"][values] for values in self.config["database_specific_commercial"]]
        subs_search_methods = list(self.config["structure_search_methode"].keys())
        plate_production = ["Mother Plates", "Daughter Plates"]

        col_1 = sg.Frame("Search", [[
            sg.Column([
                [sg.T("Out put folder"),
                 sg.FolderBrowse(key="-SEARCH_OUTPUT_FOLDER-", target="-SEARCH_OUTPUT_FOLDER_TARGET-"
                                 , initial_folder="output_files")],
                [sg.Text(key="-SEARCH_OUTPUT_FOLDER_TARGET-", size=50)],
                [sg.Listbox(values=ac, key="-SEARCH_AC-", enable_events=True, size=(10, 5),
                            select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE),
                 sg.Listbox([], key="-SEARCH_ORIGIN-", size=(10, 5), select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE)],
                [sg.Checkbox("All compounds", key="-SEARCH_ALL_COMPOUNDS-", enable_events=True),
                 sg.Checkbox(text="Ignore plated compounds?", key="-SEARCH_IGNORE_PLATED_COMPOUNDS-")],
                [sg.Text("Amount of Plates", size=self.standard_size),
                 sg.InputText(key="-SEARCH_PLATE_AMOUNT-", size=10),
                 sg.Checkbox("Max amount of plates/tube-racks", key="-SEARCH_PLATE_AMOUNT_MAX-"),
                 sg.Checkbox("Minimize amount of MP use", key="-SEARCH_MP_MINIMIZED-")],
                [sg.Text("Transferee vol:", size=self.standard_size),
                 sg.InputText(key="-SEARCH_TRANS_VOL-", size=10),
                 sg.DropDown(["mL", "uL", "nL"], key="-SEARCH_VOL_PARAMETERS-", default_value="mL", size=5),
                 sg.Checkbox("Ignore volume", key="-SEARCH_IGNORE_VOLUME-")],
                [sg.Text("Plate Production"),
                 sg.DropDown(plate_production, key="-SEARCH_PLATE_PRODUCTION-", default_value=plate_production[0],
                             enable_events=True),
                 sg.DropDown(sorted(self.plate_list), key="-SEARCH_PLATE_LAYOUT-", disabled=True, enable_events=True),
                 sg.T("Samples per plate:"), sg.Input("384", key="-SEARCH_PLATE_LAYOUT_SAMPLE_AMOUNT-", size=5,
                                                      readonly=True, disabled_readonly_text_color="#FFFFFF",
                                                      disabled_readonly_background_color="#4D4D4D")]
            ])
        ]])

        col_sub_search = sg.Frame("Structure Search", [[
            sg.Column([
                [sg.Checkbox(text="Structure Search", key="-SUB_SEARCH-")],
                [sg.Text("Smiles", size=self.standard_size), sg.InputText(key="-SUB_SEARCH_SMILES-", size=self.standard_size),
                 sg.Button("Draw molecule", key="-SUB_SEARCH_DRAW_MOL-")],
                [sg.Text("Search Method", size=self.standard_size),
                 sg.DropDown(subs_search_methods, key="-SUB_SEARCH_METHOD-", default_value=subs_search_methods[0],
                             enable_events=True)],
                [sg.Text("Similarity Threshold", size=self.standard_size),
                 sg.InputText(key="-SUB_SEARCH_THRESHOLD-", default_text=0, size=self.standard_size)],
                [sg.HorizontalSeparator()],
                [sg.Text("Morgan specific options", key="-SUB_SEARCH_MORGAN_OPTIONS-", visible=False)],
                [sg.Checkbox(text="chirality", key="-SUB_SEARCH_MORGAN_CHIRALITY-", visible=False),
                 sg.Checkbox(text="Features", key="-SUB_SEARCH_MORGAN_FEATURES-", visible=False)],
                [sg.Text("n bits", key="-SUB_SEARCH_BITS_TEXT-", size=self.standard_size, visible=False),
                 sg.InputText(key="-SUB_SEARCH_MORGAN_BITS-", size=self.standard_size, visible=False)],
                [sg.Text("bound range", key="-SUB_SEARCH_BOUND_TEXT-", size=self.standard_size, visible=False),
                 sg.InputText(key="-SUB_SEARCH_MORGAN_RANGE-", size=self.standard_size, visible=False)],
            ])
        ]])

        layout = [sg.vtop([col_1,  col_sub_search])]

        return layout

    # def setup_1_mp_dp_update(self):
    #     """
    #
    #     :return: A layour for the bio-module in the top box
    #     :rtype: list
    #     """
    #     print("is this working?")
    #     col_bio_analysis = sg.Frame("Analyse setup", [[
    #         sg.Column([
    #             [sg.Text("Plate Layout", size=self.standard_size),
    #              sg.DropDown(sorted(self.plate_list), key="-MP_DP_UPDATE_PLATE_LAYOUT-", enable_events=True,
    #                          size=self.standard_size)],
    #
    #         ])
    #     ]])
    #
    #     col_extra = sg.Frame("Extra Settings", [[
    #         sg.Column([
    #         ])
    #     ]])
    #
    #     col_graph = sg.Frame("Plate Layout", [[
    #         sg.Column([
    #             [sg.Graph(canvas_size=(250, 175), graph_bottom_left=(0, 0), graph_top_right=(250, 175),
    #                       background_color='grey', key="-MP_DP_UPDATE_CANVAS-", enable_events=False, drag_submits=False,
    #                       motion_events=False)],
    #             [sg.Text("Sample:", size=self.standard_size),
    #              sg.T(background_color=self.config["plate_colouring"]["sample"], size=10,
    #                   key="-MP_DP_UPDATE_PLATE_LAYOUT_COLOUR_BOX_SAMPLE-", relief="groove")],
    #             [sg.Text("Blank:", size=self.standard_size),
    #              sg.T(background_color=self.config["plate_colouring"]["blank"], size=10,
    #                   key="-MP_DP_UPDATE_PLATE_LAYOUT_COLOUR_BOX_BLANK-", relief="groove")],
    #             [sg.Text("Maximum:", size=self.standard_size),
    #              sg.T(background_color=self.config["plate_colouring"]["max"], size=10,
    #                   key="-MP_DP_UPDATE_PLATE_LAYOUT_COLOUR_BOX_NAX-", relief="groove")],
    #             [sg.Text("Minimum:", size=self.standard_size),
    #              sg.T(background_color=self.config["plate_colouring"]["minimum"], size=10,
    #                   key="-MP_DP_UPDATE_PLATE_LAYOUT_COLOUR_BOX_MINIMUM-", relief="groove")],
    #             [sg.Text("Positive Control:", size=self.standard_size),
    #              sg.T(background_color=self.config["plate_colouring"]["positive"], size=10,
    #                   key="-MP_DP_UPDATE_PLATE_LAYOUT_COLOUR_BOX_POSITIVE-", relief="groove")],
    #             [sg.Text("Negative Control:", size=self.standard_size),
    #              sg.T(background_color=self.config["plate_colouring"]["negative"], size=10,
    #                   key="-MP_DP_UPDATE_PLATE_LAYOUT_COLOUR_BOX_NEGATIVE-", relief="groove")],
    #             [sg.Text("Empty:", size=self.standard_size),
    #              sg.T(background_color=self.config["plate_colouring"]["empty"], size=10,
    #                   key="-MP_DP_UPDATE_PLATE_LAYOUT_COLOUR_BOX_EMPTY-", relief="groove")],
    #         ])
    #     ]])
    #
    #
    #     layout = [sg.vtop([col_bio_analysis, col_extra, col_graph])]
    #
    #     return layout

    def setup_1_bio(self):
        """

        :return: A layour for the bio-module in the top box
        :rtype: list
        """

        # maybe make it a config file. and update the config file when new method are added ?
        responsible = [keys for keys in list(self.config["Responsible"].keys())]

        # Colours for the heatmap, if added back in
        # colours = [keys for keys in list(self.config["colours to hex"].keys())]

        col_bio_analysis = sg.Frame("Analyse setup", [[
            sg.Column([
                [sg.FolderBrowse(button_text="Import Folder", key="-BIO_IMPORT_FOLDER-", target="-BIO_IMPORT_TARGET-")],
                [sg.Text(key="-BIO_IMPORT_TARGET-", size=self.standard_size*2)],
                [sg.FolderBrowse(button_text="Export Folder", key="-BIO_EXPORT_FOLDER-", target="-BIO_EXPORT_TARGET-")],
                [sg.Text(key="-BIO_EXPORT_TARGET-", size=self.standard_size*2)],
                [sg.Checkbox("Same layout for all plates?", key="-BIO_PLATE_LAYOUT_CHECK-", default=True,
                             tooltip="Will use the chosen platelayout in the dropdown for all plates if True, "
                                     "else there will be a popup where you can choose the layout for each plate. "
                                     "Default will be the layout chosen in the platelayouyt dropdown")],
                [sg.Text("Plate Layout", size=self.standard_size),
                 sg.DropDown(sorted(self.plate_list), key="-BIO_PLATE_LAYOUT-", enable_events=True,
                             size=self.standard_size)],
                [sg.Text("Sample Type (Not working)", size=self.standard_size),
                 sg.DropDown(self.sample_style, key="-BIO_SAMPLE_TYPE-", default_value=self.sample_style[0],
                             size=self.standard_size, enable_events=True,
                             tooltip="This indicates how many times each sample is on the plate. "
                                     "Choosing custome, will let you choose more than 3 times.")],
                [sg.Text("Analyse Style (Not working)", size=self.standard_size),
                 sg.DropDown(self.analyse_style, key="-BIO_ANALYSE_TYPE-", size=self.standard_size, enable_events=True,
                             default_value=self.analyse_style[0],
                             tooltip="What style to use for analysing the samples")],
                [sg.T("Sample Type:"),
                 sg.Radio("Use Layout", group_id="-BIO_INFO_SAMPLE_TYPE_RADIO-",
                          key="-BIO_INFO_SAMPLE_TYPE_RADIO_LAYOUT-"),
                 sg.Radio("Use Sample ID", group_id="-BIO_INFO_SAMPLE_TYPE_RADIO-",
                          key="-BIO_INFO_SAMPLE_TYPE_RADIO_ID-")],
                # [sg.Checkbox("heatmap", key="-BIO_HEATMAP-")],
                # [sg.Text("start Colour:", size=self.standard_size),
                #  sg.DropDown(colours, key="-HEAT_START-", size=self.standard_size,
                #              default_value=colours[0])],
                # [sg.Text("Mid Colour:", size=self.standard_size),
                #  sg.DropDown(colours, key="-HEAT_MID-", size=self.standard_size,
                #              default_value=colours[13])],
                # [sg.Text("End Colour:", size=self.standard_size),
                #  sg.DropDown(colours, key="-HEAT_END-", size=self.standard_size,
                #              default_value=colours[4])],
                # [sg.Checkbox("State colours", key="-BIO_STATE-")],
                [sg.Button("Calculate", key="-BIO_CALCULATE-",
                           tooltip="Will do calculations on all the files in the chosen folder. "
                                   "Data can be exported to Excel and/or imported to the Database"), sg.Push(),
                 sg.Button("Send to Info", key="-BIO_SEND_TO_INFO-")]
            ])
        ]])

        col_extra = sg.Frame("Extra Settings", [[
            sg.Column([
                [sg.Checkbox("Compound Related Assay", key="-BIO_COMPOUND_DATA-", enable_events=True,
                             tooltip="Will ask the user for a worklist when analysing the data, to track what compounds"
                                     "is in each well")
                 , sg.Push(), sg.B("Report Settings", key="-BIO_REPORT_SETTINGS-", size=12)],
                [sg.Checkbox("Export to Excel", key="-BIO_EXPORT_TO_EXCEL-", default=True)],
                [sg.Checkbox("Add To Database", key="-BIO_EXPERIMENT_ADD_TO_DATABASE-", enable_events=True,
                             tooltip="Will add the data to the database")],
                [sg.Checkbox("Add Compound ID", key="-BIO_REPORT_ADD_COMPOUND_IDS-", enable_events=True,
                             tooltip="Will add ID's to each report")],
                [sg.T("Assay:", size=12),
                 sg.DropDown(values=[], key="-BIO_ASSAY_NAME-", size=14,
                             tooltip="If you want to add new data to an assay that have been run before"),
                 sg.Button("New Assay", size=12, key="-BIO_NEW_ASSAY-",
                           tooltip="Will give a pop-up where you can fill in the data needed for an assay")],
                [sg.T("Responsible:", size=12),
                 sg.DropDown(responsible, key="-BIO_RESPONSIBLE-", size=14,
                             tooltip="The main responsible for the data")],
                [sg.HorizontalSeparator()],
                [sg.Checkbox("Combined Report ", key="-BIO_COMBINED_REPORT-", default=False, enable_events=True)],
                [sg.T("Report Name:", size=12),
                 sg.InputText(key="-FINAL_BIO_NAME-", size=14, tooltip="The Name the final report will be saved as")],
                [sg.Checkbox("Include Hits", key="-BIO_FINAL_REPORT_INCLUDE_HITS-",
                             tooltip="Include a list of compounds with a score lower than the threshold sat, "
                                     "or x-amount of the lowest once, depending on the Hit Amount sat",
                             disabled=True, enable_events=True),
                 sg.Checkbox("Include smiles", key="-BIO_FINAL_REPORT_INCLUDE_SMILES-",
                             tooltip="Include the smiles for the Hits.",
                             disabled=True, enable_events=True)],
                [sg.Checkbox("Use Threshold", key="-BIO_FINAL_REPORT_USE_THRESHOLD-",
                             disabled=True, enable_events=True),
                 sg.Checkbox("Use Amount", key="-BIO_FINAL_REPORT_USE_AMOUNT-",
                             disabled=True, enable_events=True)],
                [sg.T("Threshold", size=12,
                      tooltip="Will use threshold, for hits. Can't use threshold and 'Hit Amount' at the same time"),
                 sg.T("Hit Amount", size=12,
                      tooltip="Will use amount of hits. Can't use threshold and 'Hit Amount' at the same time")],
                [sg.InputText(key="-BIO_FINAL_REPORT_THRESHOLD-", size=14, disabled=True,
                              tooltip="The threshold where samples should be included. Any sample with a score "
                                      "lower than the threshold will be included in the report"),
                 sg.InputText(key="-BIO_FINAL_REPORT_HIT_AMOUNT-", size=14, disabled=True,
                              tooltip="The Amount of sample to list. Will be sorted after lowest score")],
                # [sg.Button("Select Worklist", key="-BIO_SELECT_WORKLIST-",
                #            tooltip="Use Control til select multiple files")],
                # [sg.T(key="-BIO_SAMPLE_LIST_TARGET-")]
            ])
        ]])

        col_graph = sg.Frame("Plate Layout", [[
            sg.Column([
                [sg.Graph(canvas_size=(250, 175), graph_bottom_left=(0, 0), graph_top_right=(250, 175),
                          background_color='grey', key="-BIO_CANVAS-", enable_events=False, drag_submits=False,
                          motion_events=False)],
                [sg.Text("Sample:", size=self.standard_size),
                 sg.T(background_color=self.config["plate_colouring"]["sample"], size=10,
                      key="-BIO_PLATE_LAYOUT_COLOUR_BOX_SAMPLE-", relief="groove")],
                [sg.Text("Blank:", size=self.standard_size),
                 sg.T(background_color=self.config["plate_colouring"]["blank"], size=10,
                      key="-BIO_PLATE_LAYOUT_COLOUR_BOX_BLANK-", relief="groove")],
                [sg.Text("Maximum:", size=self.standard_size),
                 sg.T(background_color=self.config["plate_colouring"]["max"], size=10,
                      key="-BIO_PLATE_LAYOUT_COLOUR_BOX_NAX-", relief="groove")],
                [sg.Text("Minimum:", size=self.standard_size),
                 sg.T(background_color=self.config["plate_colouring"]["minimum"], size=10,
                      key="-BIO_PLATE_LAYOUT_COLOUR_BOX_MINIMUM-", relief="groove")],
                [sg.Text("Positive Control:", size=self.standard_size),
                 sg.T(background_color=self.config["plate_colouring"]["positive"], size=10,
                      key="-BIO_PLATE_LAYOUT_COLOUR_BOX_POSITIVE-", relief="groove")],
                [sg.Text("Negative Control:", size=self.standard_size),
                 sg.T(background_color=self.config["plate_colouring"]["negative"], size=10,
                      key="-BIO_PLATE_LAYOUT_COLOUR_BOX_NEGATIVE-", relief="groove")],
                [sg.Text("Empty:", size=self.standard_size),
                 sg.T(background_color=self.config["plate_colouring"]["empty"], size=10,
                      key="-BIO_PLATE_LAYOUT_COLOUR_BOX_EMPTY-", relief="groove")],
            ])
        ]])


        layout = [sg.vtop([col_bio_analysis, col_extra, col_graph])]

        return layout

    def setup_1_purity(self):
        """

        :return: A layour for the purity-module in the top box
        :rtype: list
        """
        ms_mode = ["Positive", "Negative"]
        responsible = [keys for keys in list(self.config["Responsible"].keys())]

        col_ms_buttons = sg.Frame("bah", [[
            sg.Column([
                [sg.FolderBrowse(button_text="Import Folder", key="-PURITY_DATA_IMPORT_FOLDER-",
                                 target="-PURITY_DATA_IMPORT_TARGET-")],
                [sg.Text(key="-PURITY_DATA_IMPORT_TARGET-", size=self.standard_size * 2)],
                [sg.FileBrowse(button_text="Compound data", key="-PURITY_DATA_COMPOUND_DATA-",
                               target="-PURITY_DATA_COMPOUND_DATA_TARGET-")],
                [sg.Text(key="-PURITY_DATA_COMPOUND_DATA_TARGET-", size=self.standard_size * 2)],
                [sg.T("Responsible"), sg.DropDown(responsible, key="-PURITY_DATA_RESPONSIBLE-")],
                [sg.Checkbox("Compound data?", key="-PURITY_DATA_USE_COMPOUNDS-"),
                 sg.Checkbox("Add to Database", key="-PURITY_DATA_ADD_TO_DATABASE-")],
                [sg.Checkbox("Calculate purity?", key="-PURITY_DATA_CALC_PURITY-")],
                [sg.B("Generate Report", key="-PURITY_DATA_REPORT-"), sg.Push(),
                 sg.B("Import to info", key="-PURITY_DATA_IMPORT-")]
            ])
        ]])

        col_ms_settings = sg.Frame("MS setup", [[
            sg.Column([
                [sg.Push(), sg.B("Advanced setting", key="-PURITY_ADVANCED_SETTINGS-")],
                [sg.T("UV wavelength", size=self.standard_size),
                 sg.InputText(key="-PURITY_DATA_UV_WAVE-",
                              default_text=self.config["MS_default"]["uv_wavelength"], size=10)],
                [sg.HorizontalSeparator()],
                [sg.Text("UV threshold", size=self.standard_size),
                 sg.InputText(key="-PURITY_DATA_UV_THRESHOLD-",
                              default_text=int(self.config["MS_default"]["uv_threshold"]), size=10)],
                [sg.T("Slope Threshold", size=self.standard_size),
                 sg.InputText(key="-PURITY_DATA_SLOPE_THRESHOLD-",
                              default_text=int(self.config["MS_default"]["slop_threshold"]), size=10)],
                [sg.Text("Solvent peak retention time", size=self.standard_size),
                 sg.InputText(key="-PURITY_DATA_RT_SOLVENT-",
                              default_text=float(self.config["MS_default"]["rt_solvent"]), size=10)],
                [sg.HorizontalSeparator()],
                [sg.DropDown(ms_mode, key="-PURITY_DATA_MS_MODE-", default_value=ms_mode[0])],
                [sg.Text("Delta MS", size=self.standard_size),
                 sg.InputText(key="-PURITY_DATA_MS_DELTA-",
                              default_text=float(self.config["MS_default"]["ms_delta"]), size=10)],
                [sg.Text("MS threshold", size=self.standard_size),
                 sg.InputText(key="-PURITY_DATA_MS_THRESHOLD-",
                              default_text=int(self.config["MS_default"]["ms_threshold"]), size=10)],
                [sg.Text("MS peak amounts", size=self.standard_size),
                 sg.InputText(key="-PURITY_DATA_MS_PEAKS-",
                              default_text=int(self.config["MS_default"]["ms_peak_amount"]), size=10)]

            ])
        ]])

        layout = [sg.vtop([col_ms_buttons, col_ms_settings])]

        return layout

    def setup_1_plate_layout(self):
        """

        :return: A layour for the plate-module in the top box
        :rtype: list
        """
        color_select = {}
        for keys in list(self.config["plate_colouring"].keys()):
            color_select[keys] = self.config["plate_colouring"][keys]

        plate_type = ["plate_96", "plate_384", "plate_1536"]
        sample_type = self.sample_style


        col_graph = sg.Frame("Plate Layout", [[
            sg.Column([
                [sg.Graph(canvas_size=(500, 350), graph_bottom_left=(0, 0), graph_top_right=(500, 350),
                      background_color='grey', key="-RECT_BIO_CANVAS-", enable_events=True, drag_submits=True,
                      motion_events=True)],
                [sg.DropDown(values=plate_type, default_value=plate_type[1], key="-PLATE-"),
                 sg.B("Draw Plate", key="-DRAW-"),
                 # sg.B("Add sample layout", key="-DRAW_SAMPLE_LAYOUT-"),
                 sg.Text(key="-INFO-")]
            ])
        ]])

        col_options = sg.Frame("Options", [[
            sg.Column([
                [sg.Text("Choose what clicking a figure does", enable_events=True)],
                [sg.Radio(f"Select Sample", 1, key="-RECT_SAMPLES-", size=15, enable_events=True,
                          default=True),
                 sg.T(background_color=self.config["plate_colouring"]["sample"], size=10,
                      key="-PLATE_LAYOUT_COLOUR_BOX_SAMPLE-", relief="groove"),
                 sg.DropDown(sample_type, key="-RECT_SAMPLE_TYPE-", default_value=sample_type[0])],
                [sg.Radio(f"Select Blank", 1, key="-RECT_BLANK-", size=15, enable_events=True),
                 sg.T(background_color=self.config["plate_colouring"]["blank"], size=10,
                      key="-PLATE_LAYOUT_COLOUR_BOX_BLANK-", relief="groove")],
                [sg.Radio(f"Select Max Signal", 1, key="-RECT_MAX-", size=15, enable_events=True),
                 sg.T(background_color=self.config["plate_colouring"]["max"], size=10,
                      key="-PLATE_LAYOUT_COLOUR_BOX_NAX-", relief="groove")],
                [sg.Radio(f"Select Minimum Signal", 1, key="-RECT_MIN-", size=15,
                          enable_events=True),
                 sg.T(background_color=self.config["plate_colouring"]["minimum"], size=10,
                      key="-PLATE_LAYOUT_COLOUR_BOX_MINIMUM-", relief="groove")],
                [sg.Radio(f"Select Positive Control", 1, key="-RECT_POS-", size=15,
                          enable_events=True),
                 sg.T(background_color=self.config["plate_colouring"]["positive"], size=10,
                      key="-PLATE_LAYOUT_COLOUR_BOX_POSITIVE-", relief="groove")],
                [sg.Radio(f"Select Negative Control", 1, key="-RECT_NEG-", size=15,
                          enable_events=True),
                 sg.T(background_color=self.config["plate_colouring"]["negative"], size=10,
                      key="-PLATE_LAYOUT_COLOUR_BOX_NEGATIVE-", relief="groove")],
                [sg.Radio(f"Select Empty", 1, key="-RECT_EMPTY-", size=15, enable_events=True),
                 sg.T(background_color=self.config["plate_colouring"]["empty"], size=10,
                      key="-PLATE_LAYOUT_COLOUR_BOX_EMPTY-", relief="groove")],
                [sg.Radio(f"Colour", 1, key="-COLOUR-", enable_events=True, size=self.standard_size),
                 sg.ColorChooserButton("Colour", key="-PLATE_LAYOUT_COLOUR_CHOSE-",
                                       target="-PLATE_LAYOUT_COLOUR_CHOSE_TARGET-"),
                 sg.Input(key="-PLATE_LAYOUT_COLOUR_CHOSE_TARGET-", visible=False, enable_events=True, disabled=True,
                          default_text="#ffffff")],
                # [sg.Radio('Erase', 1, key='-ERASE-', enable_events=True)],
                # [sg.Radio('Move Stuff', 1, key='-MOVE-', enable_events=True)],
                [sg.Checkbox("Use Archive", default=False, key="-ARCHIVE-", size=self.standard_size),
                 sg.DropDown(sorted(self.plate_list), key="-ARCHIVE_PLATES-", size=self.standard_size)],
                [sg.Button("Delete Layout", key="-DELETE_LAYOUT-"), sg.Button("Rename Layout", key="-RENAME_LAYOUT-"),
                 sg.Button("Export", key="-EXPORT_LAYOUT-")],
                [sg.Push(), sg.Button("Save Layout", key="-SAVE_LAYOUT-")],
            ])
        ]])

        layout = [[col_graph, col_options]]

        return layout

    def set_1_worklist(self):
        text_size_short = 10
        text_size_long = 14
        input_size_short = 5
        input_size_long = 10
        
        dropdown_size = 13
        col_basic_setup = sg.Frame("Setup", [[
            sg.Column([
                [sg.T("Assay Name:", size=text_size_short),
                 sg.InputText(key="-WORKLIST_ASSAY_NAME-", size=input_size_long,
                              tooltip="The name of the assay. "
                                      "Will be used for destination plate names, and folder name")],
                [sg.Text("Plate Amount:", size=text_size_short),
                 sg.InputText(key="-WORKLIST_PLATE_AMOUNT-", size=input_size_short,
                              tooltip="How many Destination plates should there be. "
                                      "This do not take into account the initial plate value")],
                [sg.T("Initial Plate:", size=text_size_short),
                 sg.InputText(key="-WORKLIST_INITIAL_PLATE-", size=input_size_short, default_text="1",
                              tooltip="What number should the Destination Plate start at for this worklist")],
                [sg.T("Volume:", size=text_size_short),
                 sg.InputText(key="-WORKLIST_VOLUME-", size=input_size_short,
                              tooltip="How much volume needs to be added to each well, in nL")],
                [sg.Text("Plate Layout:", size=text_size_short),
                 sg.DropDown(sorted(self.plate_list), key="-WORKLIST_PLATE_LAYOUT-", size=dropdown_size,
                             tooltip="What layout to use. Check under Bio Data or Plate Layout to see what "
                                     "state each well is in")],
                [sg.Button("Generate", key="-WORKLIST_GENERATE-", tooltip="This will generate the worklist")]
            ])
        ]])


        sample_directions = ["Vertical", "Horizontale"]
        col_advance_setup = sg.Frame("Extra settings", [[
            sg.Column([
                [sg.T("Sample Style", size=text_size_long),
                 sg.DropDown(values=self.sample_style, key="-WORKLIST_SAMPLE_STYLE-", size=dropdown_size,
                             default_value=self.sample_style[0],
                             tooltip="This determines how many wells each compounds goes to."
                                     "Using custom will give you options to chose any number.")],
                [sg.T("Analyse Style", size=text_size_long),
                 sg.DropDown(values=self.analyse_style, key="-WORKLIST_ANALYSE_STYLE-", size=dropdown_size,
                             default_value=self.analyse_style[0],
                             tooltip="This is a choose between how the sample are layout.")],
                [sg.Text("Sample Direction", size=text_size_long),
                 sg.DropDown(values=sample_directions, default_value=sample_directions[0], size=dropdown_size,
                             key="-WORKLIST_DROPDOWN_SAMPLE_DIRECTION-")],
                [sg.Checkbox("Use Positive Control?", key="-WORKLIST_USE_POSITIVE_CONTROL-", enable_events=True)],
                [sg.T("Positive Control ID:", size=text_size_long),
                 sg.InputText(key="-WORKLIST_POSITIVE_CONTROL_ID-", size=input_size_long, disabled=True,
                              tooltip="Use the ID for compound. "
                                      "The ID should fit with the naming scheme in the Plate Layout file")],
                [sg.Checkbox("Use Negative Control?", key="-WORKLIST_USE_NEGATIVE_CONTROL-", enable_events=True)],
                [sg.T("Positive Control ID:", size=text_size_long),
                 sg.InputText(key="-WORKLIST_NEGATIVE_CONTROL_ID-", size=input_size_long, disabled=True,
                              tooltip="Use the ID for compound. "
                                      "The ID should fit with the naming scheme in the Plate Layout file")],
                [sg.Checkbox("Use Bonus Compound", key="-WORKLIST_USE_BONUS_COMPOUND-", enable_events=True,
                             tooltip="This will add this compound to all selected well states")],
                [sg.Text("Compound Name:", size=text_size_long),
                 sg.InputText(key="-WORKLIST_BONUS_COMPOUND_ID-", size=input_size_long, disabled=True,
                              tooltip="Name Needs to fit with a name in the Control Layout.")],
                [sg.Checkbox("Max", key="-WORKLIST_BONUS_MAX-"),
                 sg.Checkbox("Positive", key="-WORKLIST_BONUS_POSITIVE-"),
                 sg.Checkbox("Empty", key="-WORKLIST_BONUS_EMPTY-")],
                [sg.Checkbox("Min", key="-WORKLIST_BONUS_MIN-"),
                 sg.Checkbox("Negative", key="-WORKLIST_BONUS_NEGATIVE-"),
                 sg.Checkbox("Blank", key="-WORKLIST_BONUS_BLANK-"),
                 sg.Checkbox("Sample", key="-WORKLIST_BONUS_SAMPLE-")],
                [sg.Button("Control Layout",
                               key="-WORKLIST_CONTROL_LAYOUT-",
                               tooltip="Choose a Plate Layout where the controls are located. "
                                   "It should be either CSV formate or excel."),
                 sg.InputText(key="-WORKLIST_CONTROL_LAYOUT_TARGET-", visible=False)]
            ])
        ]])

        motherplates=[]

        col_mp_settings = sg.Frame("Mother Plates to use", [[
            sg.Column([
                [sg.Listbox(values=motherplates, key="-WORKLIST_MP_LIST-", size=(15, 7),
                            select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE,
                            tooltip="Multiple plates can be selected")],
                [sg.Checkbox("Use all", key="-WORKLIST_USE_ALL_MOTHERPLATES-",
                             tooltip="Will use as many MP as needed, in order. Will compare to previous worklist, "
                                     "and skipp any duplicated compounds.")]
            ])
        ]])

        col_assay_settings = sg.Frame("Previuse assays", [[
            sg.Column([
                [sg.Listbox(values=motherplates, key="-WORKLIST_ASSAY_LIST-", size=(15, 7),
                            select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE,
                            tooltip='Multiple plates can be selected')],
                [sg.DropDown(values=[], key="-WORKLIST_ASSAY_LIST_DROPDOWN-", size=dropdown_size)]
            ])
        ]])

        layout = [sg.vtop([col_basic_setup, col_advance_setup, col_mp_settings, col_assay_settings])]

        return layout

    def setup_1_update(self):
        """

        :return:A layour for the update-modul in the top box
        :rtype: list
        """

        col_buttons = sg.Frame("Update", [[
            sg.Column([
                [sg.T("Import folder:"),
                 sg.FolderBrowse(key="-UPDATE_FOLDER-", target="-FOLDER_PATH-",
                                 size=(self.standard_size, self.button_height))],
                [sg.Text("output_file", key="-FOLDER_PATH-", size=50)],
                [sg.Button("Add compounds", key="-UPDATE_COMPOUND-", size=(self.standard_size, self.button_height)),
                 sg.Button("Auto", key="-UPDATE_AUTO-", size=(self.standard_size, self.button_height))],
                [sg.Button("Add MotherPlates", key="-UPDATE_MP-", size=(self.standard_size, self.button_height)),
                 sg.Button("Add DaughterPlates", key="-UPDATE_DP-", size=(self.standard_size, self.button_height))],
            ])
        ]])

        layout = [sg.vtop([col_buttons])]

        return layout
        #return sg.Window("Update database", layout, size, finalize=True)

    def setup_1_extra_files(self):

        # Plate Dilution
        pd_dd = ["Calculate", "Generate"]
        pd_well_layout = ["Row", "Column"]
        pd_dilution_well_amount = ["Copy", "Minimum", "Input"]
        pd_source_well_amount = ["Minimum", "Input"]
        plate_dilution = sg.Frame("Plate Dilution", [[
            sg.Column([
                [sg.FileBrowse(key="-PD_FILE-", target="-PD_FILE_TARGET-")],
                [sg.T(key="-PD_FILE_TARGET-")],
                [sg.T("Method", size=18),
                 sg.DropDown(pd_dd, key="-PD_METHOD_DD-", default_value=pd_dd[0], size=10,
                             enable_events=True)],
                [sg.Checkbox("Add Source-Wells?", key="-PD_ADD_SOURCE_WELLS-", enable_events=True)],
                [sg.T("Source Well amount:", size=18),
                 sg.DropDown(pd_source_well_amount, key="-PD_SOURCE_WELL_AMOUNT-", disabled=True,
                             default_value=pd_source_well_amount[0], size=10),
                 sg.Input(key="-PD_SOURCE_WELL_AMOUNT_INPUT-", size=10)],
                [sg.T("Sample well layout", size=18),
                 sg.DropDown(pd_well_layout, key="-PD_WELL_LAYOUT-", default_value=pd_well_layout[0], size=10,
                             disabled=True)],
                [sg.T("Dilution well amount:", size=18),
                 sg.DropDown(pd_dilution_well_amount, key="-PD_WELL_AMOUNT-", default_value=pd_dilution_well_amount[0],
                             size=10, disabled=True),
                 sg.Input(key="-PD_WELL_AMOUNT_INPUT-", size=10)],
                [sg.T("Plate Layout:", size=18),
                 sg.DropDown(sorted(self.plate_list), key="-DP_PLATE_LAYOUT-", disabled=True)],
                [sg.Checkbox("Save Plates?", key="-PD_SAVE_PLATES-", disabled=True),
                 sg.Checkbox("Generate PB-file for Source plates?", key="-DP_SOURCE_FILE_GENERATE-", disabled=True)],
                [sg.Button("Execute", key="-PD_EXECUTE_BUTTON-")]
            ])
        ]], expand_y=True, expand_x=True)

        # Echo data file convert
        echo_dropdown = ["only non-zero", "destin+source plate", "plate counter", "all data"]
        echo_data = sg.Frame("Echo Data", [[
            sg.Column([
                [sg.FileBrowse("Choose file", key="-EXTRA_ECHO_DATA_FILE_BROWS-",
                               target="-EXTRA_ECHO_DATA_FILE_BROWS_TARGET-"),
                 sg.FolderBrowse("Choose Folder", key="-EXTRA_ECHO_DATA_FOLDER_BROWS-",
                                 target="-EXTRA_ECHO_DATA_FOLDER_BROWS_TARGET-")],
                [sg.T("File:"), sg.T("",key="-EXTRA_ECHO_DATA_FILE_BROWS_TARGET-")],
                [sg.T("Folder:"), sg.T("", key="-EXTRA_ECHO_DATA_FOLDER_BROWS_TARGET-")],
                [sg.DropDown(echo_dropdown, key="-EXTRA_ECHO_DATA_DROPDOWN-", default_value=echo_dropdown[0])],
                [sg.B("Excel", key="-EXTRA_ECHO_DATA_EXCEL-"), sg.B("CSV", key="-EXTRA_ECHO_DATA_CSV-")]
            ])
        ]], expand_y=True, expand_x=True)

        tab_plate_dilution = sg.Tab("Plate Dilution", [[plate_dilution]])
        tab_echo_data = sg.Tab("Echo Data", [[echo_data]])

        tab_list = [tab_plate_dilution, tab_echo_data]

        tab_groups_files = sg.Tab("File Handler", [[sg.TabGroup([tab_list], selected_background_color=self.tab_colour,
                                  key="-EXTRA_SUB_FILES_TABS-",  enable_events=True, tab_location="lefttop",
                                  expand_x=True, expand_y=True)]])

        return tab_groups_files

    def setup_1_extra_database(self):
        text_size = 15

        #ToDo Add "remove and update" options to all the tabs

        resp_headings = self.config["Extra_tab_database_headings"]["responsible"].split(",")
        responsible = sg.Frame("Responsible", [
            sg.vtop([
                sg.Column([
                    [sg.T("Name", size=text_size), sg.Input(key="-EXTRA_DATABASE_RESPONSIBLE_NAME-", size=10)],
                    [sg.T("E-mail", size=text_size), sg.Input(key="-EXTRA_DATABASE_RESPONSIBLE_E_MAIL-", size=10)],
                    [sg.T("Info", size=text_size), sg.Input(key="-EXTRA_DATABASE_RESPONSIBLE_INFO-", size=10)],
                    [sg.Button("Import to DB", key="-EXTRA_DATABASE_RESPONSIBLE_IMPORT_DB-"),
                     sg.Button("Update selected", key="-EXTRA_DATABASE_RESPONSIBLE_UPDATE_SELECTED-"),
                     sg.Button("Import from file", key="-EXTRA_DATABASE_RESPONSIBLE_IMPORT_FILE-")]
                ]),
                sg.Column([
                   [sg.Table(values=[], headings=resp_headings, key="-EXTRA_DATABASE_RESPONSIBLE_TABLE-")]
                ])
            ])
        ])

        customers_headings = self.config["Extra_tab_database_headings"]["customers"].split(",")
        customers = sg.Frame("Customers", [
            sg.vtop([
                sg.Column([
                    [sg.T("Name", size=text_size), sg.Input(key="-EXTRA_DATABASE_CUSTOMERS_NAME-", size=10)],
                    [sg.T("Contact Person", size=text_size), sg.Input(key="-EXTRA_DATABASE_CUSTOMERS_CONTACT_NAME-", size=10)],
                    [sg.T("E-mail", size=text_size), sg.Input(key="-EXTRA_DATABASE_CUSTOMERS_E_MAIL-", size=10)],
                    [sg.T("Info", size=text_size), sg.Input(key="-EXTRA_DATABASE_CUSTOMERS_INFO-", size=10)],
                    [sg.Button("Import to DB", key="-EXTRA_DATABASE_CUSTOMERS_IMPORT_DB-"),
                     sg.Button("Update selected", key="-EXTRA_DATABASE_CUSTOMERS_UPDATE_SELECTED-"),
                     sg.Button("Import from file", key="-EXTRA_DATABASE_CUSTOMERS_IMPORT_FILE-")]
                ]),
                sg.Column([
                    [sg.Table(values=[], headings=customers_headings, key="-EXTRA_DATABASE_CUSTOMERS_TABLE-")]
                ])
            ])
        ])

        vendor_headings = self.config["Extra_tab_database_headings"]["vendors"].split(",")
        vendors = sg.Frame("Vendors", [
            sg.vtop([
                sg.Column([
                    [sg.T("Name", size=text_size), sg.Input(key="-EXTRA_DATABASE_VENDORS_NAME-", size=10)],
                    [sg.T("Contact Person", size=text_size), sg.Input(key="-EXTRA_DATABASE_VENDORS_CONTACT_NAME-", size=10)],
                    [sg.T("E-mail", size=text_size), sg.Input(key="-EXTRA_DATABASE_VENDORS_E_MAIL-", size=10)],
                    [sg.T("Info", size=text_size), sg.Input(key="-EXTRA_DATABASE_VENDORS_INFO-", size=10)],
                    [sg.Button("Import to DB", key="-EXTRA_DATABASE_VENDORS_IMPORT_DB-"),
                     sg.Button("Update selected", key="-EXTRA_DATABASE_VENDORS_UPDATE_SELECTED-"),
                     sg.Button("Import from file", key="-EXTRA_DATABASE_VENDORS_IMPORT_FILE-")]
                ]),
                sg.Column([
                    [sg.Table(values=[], headings=vendor_headings, key="-EXTRA_DATABASE_VENDORS_TABLE-")]
                ])
            ])
        ])

        ac_headings = self.config["Extra_tab_database_headings"]["origin"].split(",")
        drop_down = ["Academia", "Commercial"]
        ac = sg.Frame("AC", [
            sg.vtop([
                sg.Column([
                    [sg.T("Name", size=text_size), sg.Input(key="-EXTRA_DATABASE_AC_NAME-", size=10)],
                    [sg.T("Contact Person", size=text_size), sg.Input(key="-EXTRA_DATABASE_AC_CONTACT_NAME-", size=10)],
                    [sg.T("E-mail", size=text_size), sg.Input(key="-EXTRA_DATABASE_AC_E_MAIL-", size=10)],
                    [sg.T("Info", size=text_size), sg.Input(key="-EXTRA_DATABASE_AC_INFO-", size=10)],
                    [sg.T("A/C", size=text_size),
                     sg.DropDown(drop_down, key="-EXTRA_DATABASE_AC_AC-", default_value=drop_down[0])],
                    [sg.Button("Import to DB", key="-EXTRA_DATABASE_AC_IMPORT_DB-"),
                     sg.Button("Update selected", key="-EXTRA_DATABASE_AC_UPDATE_SELECTED-"),
                     sg.Button("Import from file", key="-EXTRA_DATABASE_AC_IMPORT_FILE-")]
                ]),
                sg.Column([
                    [sg.Table(values=[], headings=ac_headings, key="-EXTRA_DATABASE_AC_TABLE-")]
                ])
            ])
        ])

        plate_types = sg.Frame("Plate Types", [
            sg.vtop([
                sg.Column([
                    [sg.T("Plate Type", size=text_size), sg.Input(key="-EXTRA_PLATE_TYPE_NAME-", size=10)],
                    [sg.T("Vendor", size=text_size), sg.DropDown([], key="-EXTRA_PLATE_TYPE_VENDOR-", size=10)],
                    [sg.T("Product Number", size=text_size), sg.Input(key="-EXTRA_PLATE_TYPE_PRODUCT_NUMBER-", size=10)],

                    [sg.T("Size", size=text_size), sg.Input(key="-EXTRA_PLATE_TYPE_SIZE-", size=5),
                     sg.Checkbox("Sterile", key="-EXTRA_PLATE_TYPE_STERILE-")],
                    [sg.T("Info", size=text_size), sg.Input(key="-EXTRA_PLATE_TYPE_INFO-", size=10)],
                    [sg.VPush()],
                    [sg.VPush()],
                    [sg.VPush()],
                    [sg.Button("Import to DB", key="-EXTRA_DATABASE_PLACE_TYPE_IMPORT_DB-"),
                     sg.Button("Update selected", key="-EXTRA_DATABASE_PLACE_TYPE_UPDATE_SELECTED-"),
                     sg.Button("Import from file", key="-EXTRA_DATABASE_PLACE_TYPE_IMPORT_FILE-")]
                ]),
                sg.Column([
                    [sg.T("Well Offset X (A)", size=text_size), sg.Input(key="-EXTRA_PLATE_TYPE_WELL_OFFSET_X-", size=10)],
                    [sg.T("Well Offset Y (B)", size=text_size), sg.Input(key="-EXTRA_PLATE_TYPE_WELL_OFFSET_Y-", size=10)],
                    [sg.T("Well Spacing X (C)", size=text_size), sg.Input(key="-EXTRA_PLATE_TYPE_WELL_SPACING_X-", size=10)],
                    [sg.T("Well Spacing Y (D)", size=text_size), sg.Input(key="-EXTRA_PLATE_TYPE_WELL_SPACING_Y-", size=10)],
                    [sg.T("Plate Height (E)", size=text_size), sg.Input(key="-EXTRA_PLATE_TYPE_PLATE_HEIGHT-", size=10)],
                    [sg.T("Plate Height + Lid", size=text_size), sg.Input(key="-EXTRA_PLATE_TYPE_PLATE_HEIGHT_LID-", size=10)],
                    [sg.T("Flange Height (F)", size=text_size), sg.Input(key="-EXTRA_PLATE_TYPE_FLANGE_HEIGHT-", size=10)],
                    [sg.T("Well Width (G)", size=text_size), sg.Input(key="-EXTRA_PLATE_TYPE_WELL_WIDTH-", size=10)],
                    [sg.T("Max Volume", size=text_size), sg.Input(key="-EXTRA_PLATE_TYPE_MAX_VOL-", size=10)],
                    [sg.T("Working Volume", size=text_size), sg.Input(key="-EXTRA_PLATE_TYPE_WORKING_VOL-", size=10)],
                    [sg.T("Dead volume", size=text_size), sg.Input(key="-EXTRA_PLATE_TYPE_WELL_DEAD_VOL-", size=10)]
                ]),
                sg.Column([
                   [sg.Listbox([], key="-EXTRA_PLATE_TYPE_LISTBOX-", size=(10, 10))],

                ]),
                sg.Column([
                    [sg.Canvas(size=(100, 100), key="-EXTRA_PLATE_TYPE_CANVAS-")]
                ])
            ])
        ])

        location_headings = self.config["Extra_tab_database_headings"]["locations"].split(",")
        locations = sg.Frame("Locations", [[
            sg.Column([
                [sg.T("Room", size=text_size), sg.Input(key="-EXTRA_DATABASE_LOCATIONS_ROOM-", size=10)],
                [sg.T("Building", size=text_size), sg.Input(key="-EXTRA_DATABASE_LOCATIONS_BUILDING-", size=10)],
                [sg.T("Extra Info/spot", size=text_size), sg.Input(key="-EXTRA_DATABASE_LOCATIONS_SPOT-", size=10)],
                [sg.Button("Import to DB", key="-EXTRA_DATABASE_LOCATION_IMPORT_DB-"),
                 sg.Button("Update selected", key="-EXTRA_DATABASE_LOCATION_UPDATE_SELECTED-"),
                 sg.Button("Import from file", key="-EXTRA_DATABASE_LOCATION_IMPORT_FILE-")]
            ]),
            sg.Column([
                [sg.Table(values=[], headings=location_headings, key="-EXTRA_DATABASE_LOCATIONS_TABLE-")]
            ])
        ]])

        tab_responsible = sg.Tab("Responsible", [[responsible]],
                                 tooltip="Add, remove or update people that can be responsible for an assay or other")
        tab_customers = sg.Tab("Customers", [[customers]],
                               tooltip="Add, remove or update customers for screens")
        tab_vendors = sg.Tab("Vendors", [[vendors]],
                             tooltip="Add, remove or update vendors that sells other stuff than compounds")
        tab_ac = sg.Tab("Origin", [[ac]],
                        tooltip="Add or update Academic or Commercial sites to the database, for compound origin")
        tab_plate_types = sg.Tab("Plate Types", [[plate_types]],
                                 tooltip="Add, remove or update plate-types")
        tab_location = sg.Tab("Location", [[locations]], tooltip="Add, remove or update storrage locations")


        tab_list = [tab_responsible, tab_customers, tab_vendors, tab_ac, tab_plate_types, tab_location]
        # tab_list = []
        tab_groups_database = sg.Tab("Database Handler", [[sg.TabGroup([tab_list],
                                                                       selected_background_color=self.tab_colour,
                                                                       key="-EXTRA_SUB_DATABASE_TABS-",
                                                                       enable_events=True, tab_location="lefttop",
                                                                       expand_x=True, expand_y=True)]])

        return tab_groups_database

    def setup_1_extra(self):

        tab_list_groups = [[sg.TabGroup([[self.setup_1_extra_files(), self.setup_1_extra_database()]],
                                        key="-EXTRA_SUB_TABS-", selected_background_color=self.tab_colour,
                                        expand_x=True, expand_y=True, enable_events=True)]]

        layout = tab_list_groups

        return layout

    def setup_1_simulator(self):
        """

        :return: A layour for the simulation-module in the top box
        :rtype: list
        """
        import_list = ["comPOUND", "MP Production", "DP production"]

        compound_col = sg.Frame("2D files (Compound -> 2D barcode)", [[
            sg.Column([
                [sg.Text("Input folder containing the compound txt file's"),
                 sg.FileBrowse(key="-SIM_INPUT_COMPOUND_FILE-", target="-SIM_COMPOUND_TARGET-")],
                [sg.T(key="-SIM_COMPOUND_TARGET-")],
            ])
        ]], visible=True, key="-SIM_COMPOUND_FRAME-")

        mp_production_col = sg.Frame("MP Production (2D barcode -> PB files)", [[
            sg.Column([
                [sg.Text("Input folder containing the 2D barcode txt file's (; separated)")],
                [sg.FolderBrowse(key="-SIM_INPUT_MP_FILE-", target="-SIM_MP_TARGET-")],
                [sg.T(key="-SIM_MP_TARGET-")],
                [sg.Text("MotherPlate initials", size=15),
                 sg.InputText(key="-SIM_MP_NAME-", size=10)],
                [sg.Text("Volume", size=15),
                 sg.InputText(key="-SIM_MP_VOL-", size=10)]
            ])
        ]], visible=False, key="-SIM_MP_FRAME-")

        dp_production = sg.Frame("DP production (PB files -> ECHO files)", [[
            sg.Column([
                [sg.Text("Input folder containing the MP CSV file's (; separated)"),
                 sg.FolderBrowse(key="-SIM_INPUT_DP_FILE-", target="-SIM_DP_TARGET-")],
                [sg.T(key="-SIM_DP_TARGET-")],
                [sg.Text("DP initials"),
                 sg.InputText(key="-SIM_DP_NAME-")]
            ])
        ]], visible=False, key="-SIM_DP_FRAME-")

        layout = [[sg.DropDown(import_list, size=self.standard_size, key="-SIM_INPUT_EQ-", enable_events=True,
                               default_value=import_list[0])],
                  sg.vtop([compound_col, mp_production_col, dp_production]),
                  [sg.Button("Simulate", key="-SIM_RUN-"),
                   sg.FolderBrowse("Output folder", key="-SIM_OUTPUT-", target="-SIM_OUTPUT_TArGET-"),
                   sg.T(key="-SIM_OUTPUT_TArGET-"), sg.Push(), sg.B("Start Up the Database", key="-START_UP_DB-")]]

        return layout

    def setup_2_compound(self):
        """

        :return: A layour for the info-module in the right box
        :rtype: list
        """

        col_picture = sg.Frame("picture", [[
            sg.Column([
                [sg.Image(key="-COMPOUND_INFO_PIC-", size=(500, 150))],
                [sg.Text(key="-COMPOUND_INFO_SMILES-", size=50)]
                ])
        ]])

        col_info = sg.Frame("Info", [[
            sg.Column([
                [sg.Input("", key="-COMPOUND_INFO_ID-", size=15), sg.Push(),
                 sg.Button("Search", key="-COMPOUND_INFO_SEARCH_COMPOUND_ID-")],
                [sg.Text("Academic/Commercial", size=self.standard_size),
                 sg.Text(key="-COMPOUND_INFO_AC-", size=self.standard_size)],
                [sg.Text("Origin", size=self.standard_size),
                 sg.Text(key="-COMPOUND_INFO_ORIGIN-", size=self.standard_size)],
                [sg.Text("Origin ID", size=self.standard_size),
                 sg.Text(key="-COMPOUND_INFO_ORIGIN_ID-", size=self.standard_size)],
                [sg.Text("Concentration", size=self.standard_size),
                 sg.Text(key="-COMPOUND_INFO_CONCENTRATION-", size=self.standard_size)],
                [sg.Text("volume left in MP", size=self.standard_size),
                 sg.Text(key="-COMPOUND_INFO_MP_VOLUME-", size=self.standard_size)]
                ]),
        ]])

        plate_headlines = ["name", "type", "well", "volume", "date"]
        row_plate_table = sg.Frame("Plate Table", [[
            sg.Column([
                [sg.Table(values=[], headings=plate_headlines, key="-COMPOUND_INFO_PLATE_TABLE-",
                          auto_size_columns=False, enable_click_events=True)]
            ])
        ]], expand_x=True)
        plate_info_table_headings = ["test1", "test2", "Test3"]
        tab_plates_table = sg.Tab("Plate info", [[
            sg.TabGroup([[
                sg.Tab("All Plates", [[
                    sg.Frame("Plate Table", [[
                        sg.Column([
                            [sg.Table([], headings=plate_info_table_headings,
                                      key="-COMPOUND_INFO_ALL_PLATE_INFO_TABLE-")]
                        ])
                    ]])
                ]]),
                sg.Tab("MP plates", [[
                    sg.Frame("MP Table",
                             [[sg.Column(
                                 [[sg.Table([], headings=plate_info_table_headings,
                                            key="-COMPOUND_INFO_MP_PLATE_INFO_TABLE-")]]
                             )]])
                ]]),
                sg.Tab("DP plates", [[
                    sg.Frame("DP Table",
                             [[sg.Column(
                                 [[sg.Table([], headings=plate_info_table_headings,
                                            key="-COMPOUND_INFO_DP_PLATE_INFO_TABLE-")]]
                             )]])
                ]]),
            ]], selected_background_color=self.tab_colour, key="-COMPOUND_INFO_TABLE_TABS-", enable_events=True,
                                 expand_x=True, tab_location="righttop")
        ]])

        bio_info_table_headings = ["Assay", "Responsible", "Date"]
        tab_bio_exp_table = sg.Tab("Bio Info", [[sg.Frame("Bio Experimental", [[
            sg.Column([
                [sg.Table([], headings=bio_info_table_headings, key="-COMPOUND_INFO_BIO_INFO_TABLE-")]
            ])
        ]])
                                    ]])
        purity_info_table_headings = ["Batch", "Result Max", "Result Ion", "Result Total", "Date"]
        tab_purity_table = sg.Tab("Purity Info", [[sg.Frame("Purity", [[
            sg.Column([
                [sg.Table([], headings=purity_info_table_headings, key="-COMPOUND_INFO_PURITY_INFO_TABLE-")]
            ])
        ]])
                                   ]])

        tab_groups = sg.TabGroup([[tab_plates_table, tab_bio_exp_table, tab_purity_table]],
                                 selected_background_color=self.tab_colour, key="-COMPOUND_INFO_SUB_DATA-", enable_events=True,
                                 expand_x=True)

        layout = [sg.vtop([col_info, col_picture]), [sg.VPush()], [row_plate_table], [tab_groups]]

        return layout

    def setup_2_bio(self):
        # analyse_method = ["original", "normalised", "pora"]
        analyse_method = []
        mapping = ["State Mapping", "Heatmap", "Hit Map"]

        row_settings = sg.Frame("Bio Information", [[
            sg.Column([
                        [sg.Push(),
                         sg.Checkbox("Plate Report", key="-BIO_INFO_EXPORT_PLATE_REPORT-"),
                         sg.Checkbox("Final Report", key="-BIO_INFO_EXPORT_FINAL_REPORT-"),
                         sg.B("Export", key="-BIO_INFO_EXPORT-")],
                        [sg.HorizontalSeparator()],
                        [sg.T("Analyse method", size=14),
                         sg.DropDown(analyse_method, key="-BIO_INFO_ANALYSE_METHOD-", size=14, enable_events=True),
                         sg.T("Plate Mapping", size=14),
                         sg.DropDown(mapping, key="-BIO_INFO_MAPPING-", size=14, enable_events=True)],
                        [sg.T("Plate", size=14),
                         sg.InputCombo([], key="-BIO_INFO_PLATES-", size=14, enable_events=True),
                         sg.T("State", size=14),
                         sg.InputCombo([], key="-BIO_INFO_STATES-", size=14, enable_events=True)]
                    ])
        ]], expand_x=True)

        col_cal_info = sg.Column([
            [sg.T("avg", size=14), sg.T(key="-INFO_BIO_AVG-", size=14)],
            [sg.T("stdev", size=14), sg.T(key="-INFO_BIO_STDEV-", size=14)],
            [sg.T("Z-Prime", size=14), sg.T(key="-INFO_BIO_Z_PRIME-", size=14)]
        ])

        col_well_info = sg.Column([
            [sg.T("Well id", size=14), sg.T(key="-INFO_BIO_GRAPH_TARGET-", size=14)],
            [sg.T("Well value", size=14), sg.T(key="-INFO_BIO_WELL_VALUE-", size=14)],
            [sg.T("Compound name", size=14), sg.T(key="-INFO_BIO_GRAPH_COMPOUND_NAME-", size=14)]
        ])

        row_graph = sg.Frame("Plate Layout", [[
            sg.Column([
                [sg.Graph(canvas_size=(500, 350), graph_bottom_left=(0, 0), graph_top_right=(250, 175),
                          background_color='grey', key="-BIO_INFO_CANVAS-", enable_events=True, drag_submits=True,
                          motion_events=True)],
                [col_well_info, col_cal_info]
            ])
        ]], expand_x=True)

        col_hit_mapping = sg.Column([
            [sg.T("Hit-Map Settings", relief="groove")],
            [sg.HorizontalSeparator()],
            [sg.T("TH-1 bound", size=10),
             sg.InputText(key="-BIO_INFO_PORA_LOW_MIN_HIT_THRESHOLD-", size=8, enable_events=True,
                          default_text=self.config["Settings_bio"].getfloat("final_report_pora_threshold_th_1_min")),
             sg.InputText(key="-BIO_INFO_PORA_LOW_MAX_HIT_THRESHOLD-", size=8, enable_events=True,
                          default_text=self.config["Settings_bio"].getfloat("final_report_pora_threshold_th_1_max"))],
            [sg.T("TH-2 bound", size=10),
             sg.InputText(key="-BIO_INFO_PORA_MID_MIN_HIT_THRESHOLD-", size=8, enable_events=True,
                          default_text=self.config["Settings_bio"].getfloat("final_report_pora_threshold_th_2_min")),
             sg.InputText(key="-BIO_INFO_PORA_MID_MAX_HIT_THRESHOLD-", size=8, enable_events=True,
                          default_text=self.config["Settings_bio"].getfloat("final_report_pora_threshold_th_2_max"))],
            [sg.T("TH-3 bound", size=10),
             sg.InputText(key="-BIO_INFO_PORA_HIGH_MIN_HIT_THRESHOLD-", size=8, enable_events=True,
                          default_text=self.config["Settings_bio"].getfloat("final_report_pora_threshold_th_3_min")),
             sg.InputText(key="-BIO_INFO_PORA_HIGH_MAX_HIT_THRESHOLD-", size=8, enable_events=True,
                          default_text=self.config["Settings_bio"].getfloat("final_report_pora_threshold_th_3_max"))],

            [sg.ColorChooserButton("TH-1 Values Colour", key="-INFO_BIO_HIT_MAP_LOW_COLOUR-", size=(15, None),
                                   target="-BIO_INFO_HIT_MAP_TH_1_COLOUR_TARGET-"),
             sg.T(background_color=self.config["Settings_bio"]["plate_report_pora_threshold_colour_th_1"]
                  , key="-BIO_INFO_HIT_MAP_TH_1_COLOUR_BOX-", size=8, relief="groove"),
             sg.Input(key="-BIO_INFO_HIT_MAP_TH_1_COLOUR_TARGET-", visible=False, enable_events=True, disabled=True,
                      default_text=self.config["Settings_bio"]["plate_report_pora_threshold_colour_th_1"])],

            [sg.ColorChooserButton("TH-2 Values Colour", key="-INFO_BIO_HIT_MAP_MID_COLOUR-", size=(15, None),
                                   target="-BIO_INFO_HIT_MAP_TH_2_COLOUR_TARGET-"),
             sg.T(background_color=self.config["Settings_bio"]["plate_report_pora_threshold_colour_th_2"]
                  , key="-BIO_INFO_HIT_MAP_TH_2_COLOUR_BOX-", size=8, relief="groove"),
             sg.Input(key="-BIO_INFO_HIT_MAP_TH_2_COLOUR_TARGET-", visible=False, enable_events=True, disabled=True,
                      default_text=self.config["Settings_bio"]["plate_report_pora_threshold_colour_th_2"])],

            [sg.ColorChooserButton("TH-3 Values Colour", key="-INFO_BIO_HIT_MAP_HIGH_COLOUR-", size=(15, None),
                                   target="-BIO_INFO_HIT_MAP_TH_3_COLOUR_TARGET-"),
             sg.T(background_color=self.config["Settings_bio"]["plate_report_pora_threshold_colour_th_3"]
                  , key="-BIO_INFO_HIT_MAP_TH_3_COLOUR_BOX-", size=8, relief="groove"),
             sg.Input(key="-BIO_INFO_HIT_MAP_TH_3_COLOUR_TARGET-", visible=False, enable_events=True, disabled=True,
                      default_text=self.config["Settings_bio"]["plate_report_pora_threshold_colour_th_3"])],
            [sg.Button("More bins", key="-BIO_INFO_BOUNDS_BUTTON-"), sg.Listbox(values=["1-bin", "2-bin"])]
            # [sg.T("Low Colour", size=10), sg.DropDown(colours, key="-INFO_BIO_Hit_LOW-", enable_events=True,
            #                                           default_value=self.config["Settings_bio"]
            #                                           ["plate_report_pora_threshold_colour_low"], size=14)],
            # [sg.T("Mid Colour", size=10), sg.DropDown(colours, key="-INFO_BIO_Hit_MID-", enable_events=True,
            #                                           default_value=self.config["Settings_bio"]
            #                                           ["plate_report_pora_threshold_colour_mid"], size=14)],
            # [sg.T("High Colour", size=10), sg.DropDown(colours, key="-INFO_BIO_Hit_HIGH-", enable_events=True,
            #                                            default_value=self.config["Settings_bio"]
            #                                            ["plate_report_pora_threshold_colour_high"], size=14)],
        ], key="-INFO_BIO_ROW_HIT-"
            # , size=(0, 0)
        )

        col_heatmap = sg.Column([
            [sg.T("Heatmap Settings", relief="groove")],
            [sg.HorizontalSeparator()],
            [sg.ColorChooserButton("Low Values Colour", key="-INFO_BIO_HEATMAP_LOW_COLOUR-", size=(15, None),
                                   target="-BIO_INFO_HEATMAP_LOW_COLOUR_TARGET-"),
             sg.T(background_color=self.config["Settings_bio"]["plate_report_heatmap_colours_low"]
                      , key="-BIO_INFO_HEATMAP_LOW_COLOUR_BOX-", size=8, relief="groove"),
             sg.Input(key="-BIO_INFO_HEATMAP_LOW_COLOUR_TARGET-", visible=False, enable_events=True, disabled=True,
                      default_text=self.config["Settings_bio"]["plate_report_heatmap_colours_low"])],
            [sg.ColorChooserButton("Mid Values Colour", key="-INFO_BIO_HEATMAP_MID_COLOUR-", size=(15, None),
                                   target="-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-"),
             sg.T(background_color=self.config["Settings_bio"]["plate_report_heatmap_colours_mid"]
                      , key="-BIO_INFO_HEATMAP_MID_COLOUR_BOX-", size=8, relief="groove"),
             sg.Input(key="-BIO_INFO_HEATMAP_MID_COLOUR_TARGET-", visible=False, enable_events=True, disabled=True,
                      default_text=self.config["Settings_bio"]["plate_report_heatmap_colours_mid"])],
            [sg.ColorChooserButton("High Values Colour", key="-INFO_BIO_HEATMAP_HIGH_COLOUR-", size=(15, None),
                                   target="-BIO_INFO_HEATMAP_HIGH_COLOUR_TARGET-"),
             sg.T(background_color=self.config["Settings_bio"]["plate_report_heatmap_colours_high"]
                      , key="-BIO_INFO_HEATMAP_HIGH_COLOUR_BOX-", size=8, relief="groove"),
             sg.Input(key="-BIO_INFO_HEATMAP_HIGH_COLOUR_TARGET-", visible=False, enable_events=True, disabled=True,
                      default_text=self.config["Settings_bio"]["plate_report_heatmap_colours_high"])],
            [sg.T("Low", size=10), sg.T("Mid", size=8), sg.T("High", size=8)],
            [sg.InputText(0, key="-BIO_INFO_HEAT_PERCENTILE_LOW-", size=10, enable_events=True),
             sg.InputText(50, key="-BIO_INFO_HEAT_PERCENTILE_MID-", size=10, enable_events=True),
             sg.InputText(100, key="-BIO_INFO_HEAT_PERCENTILE_HIGH-", size=10, enable_events=True)]

        ], key="-INFO_BIO_ROW_HEAT-"
            # , size=(0, 0)
        )

        col_list_box = sg.Column([
            [sg.T("States for Mapping")],
            [sg.Listbox([], key="-BIO_INFO_STATE_LIST_BOX-", size=(10, 10), enable_events=True,
                        select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE)]
        ])

        row_options = sg.Frame("Graph Settings", [
            # sg.vtop([col_heatmap, col_hit_mapping])  This does not work for some reason....
            [col_heatmap, sg.VerticalSeparator(), col_hit_mapping, sg.VerticalSeparator(), col_list_box],
            [sg.Push(), sg.Button("Re-Draw", key="-BIO_INFO_RE_DRAW-")]
            ])

        table_overview_headings = ["Method", "States", "calc", "value"]
        row_plate_overview = sg.Frame("Plate Overview", [[
            sg.Column([
                [sg.T("An overview over data per plate.")],
                [sg.Table([], headings=table_overview_headings, key="-BIO_INFO_OVERVIEW_TABLE-",
                          enable_click_events=True),
                 sg.Listbox([], key="-BIO_INFO_PLATE_OVERVIEW_METHOD_LIST-", size=(10, 10),
                            select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, enable_events=True),
                 sg.Listbox([], key="-BIO_INFO_PLATE_OVERVIEW_STATE_LIST-", size=(10, 10),
                            select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, enable_events=True)],
                [sg.DropDown([], key="-BIO_INFO_PLATE_OVERVIEW_PLATE-", size=14, enable_events=True)]
            ])
        ]])

        table_overview_headings_avg = ["Barcode", "avg"]
        table_overview_headings_stdev = ["Barcode", "stdev"]
        table_overview_headings_z_prime = ["Barcode", "Z-Prime"]
        row_overview = sg.Frame("Overview", [[
            sg.Column([
                [sg.T("Overview of plates for the different calculations")],
                [sg.Table([], headings=table_overview_headings_avg, key="-BIO_INFO_OVERVIEW_AVG_TABLE-",
                          size=(10, 10), auto_size_columns=False, enable_click_events=True),
                 sg.Table([], headings=table_overview_headings_stdev, key="-BIO_INFO_OVERVIEW_STDEV_TABLE-",
                          size=(10, 10), auto_size_columns=False, enable_click_events=True),
                 sg.Table([], headings=table_overview_headings_z_prime, key="-BIO_INFO_OVERVIEW_Z_PRIME_TABLE-",
                          size=(10, 10), auto_size_columns=False, enable_click_events=True)],
                [sg.DropDown([], key="-BIO_INFO_OVERVIEW_METHOD-", size=14, enable_events=True),
                 sg.DropDown([], key="-BIO_INFO_OVERVIEW_STATE-", size=14, enable_events=True)]
            ])
        ]])

        table_z_prime_list_headings = ["Barcode", "Z-Prime"]
        row_z_prime = sg.Frame("Z-Prime", [[
            sg.Column([
                [sg.T("Overview of Z-Prime values")],
                [sg.Table([], headings=table_z_prime_list_headings, key="-BIO_INFO_Z_PRIME_LIST_TABLE-",
                          auto_size_columns=False, enable_click_events=True),
                 sg.Column([
                     [sg.T("Max:", size=5), sg.T("", key="-BIO_INFO_Z_PRIME_MAX_BARCODE-", size=10),
                      sg.T("", key="-BIO_INFO_Z_PRIME_MAX_VALUE-", size=10)],
                     [sg.T("Min:", size=5), sg.T("", key="-BIO_INFO_Z_PRIME_MIN_BARCODE-", size=10),
                      sg.T("", key="-BIO_INFO_Z_PRIME_MIN_VALUE-", size=10)]
                 ])],
                [sg.B("Matrix", key="-BIO_INFO_Z_PRIME_MATRIX_BUTTON-")]
            ])
        ]])

        table_hit_list_headings = ["well", "value", "compound"]
        row_hit_list = sg.Frame("Hit List", [[
            sg.Column([
                [sg.T("Tables with wells that are within the hit list boundaries")],
                [sg.Table([], headings=table_hit_list_headings, key="-BIO_INFO_HIT_LIST_LOW_TABLE-"
                          , enable_click_events=True, auto_size_columns=False),
                 sg.Table([], headings=table_hit_list_headings, key="-BIO_INFO_HIT_LIST_MID_TABLE-"
                          , enable_click_events=True, auto_size_columns=False),
                 sg.Table([], headings=table_hit_list_headings, key="-BIO_INFO_HIT_LIST_HIGH_TABLE-"
                          , enable_click_events=True, auto_size_columns=False)],

                [sg.DropDown([], key="-BIO_INFO_HIT_LIST_PLATES-", size=14),
                 sg.DropDown([], key="-BIO_INFO_HIT_LIST_METHOD-", size=14),
                 sg.DropDown([], key="-BIO_INFO_HIT_LIST_STATE-", size=14)]
            ])
        ]])

        table_headings = matrix_header
        table_data = []
        col_matrix_tabl = sg.Column([
            [sg.Table(table_data, headings=table_headings, auto_size_columns=False, enable_click_events=True,
                      key="-BIO_INFO_MATRIX_TABLE-", vertical_scroll_only=False, size=(7, 7))],
            [sg.T("", size=10)],
        ])

        row_matrix = sg.Frame("Matrix", [
            [sg.T("A matrix of calculations for the different plate compared to each other, in %")],
            [col_matrix_tabl],
            [sg.DropDown([], key="-BIO_INFO_MATRIX_METHOD-", size=14),
             sg.DropDown([], key="-BIO_INFO_MATRIX_STATE-", size=14),
             sg.DropDown([], key="-BIO_INFO_MATRIX_CALC-", size=14)],
            [sg.B("Generate Matrix", key="-BIO_INFO_MATRIX_BUTTON-"),
             sg.B("Pop out the Matrix", key="-BIO_INFO_MATRIX_POPUP-")]
        ])

        # This is covered by overview
        # table_list_headings = ["barcode", "calc"]
        # row_list = sg.Frame("List", [
        #     [sg.Column([
        #         [sg.T("Table with a list of calculation values for all the plates.")],
        #         [sg.Table([], headings=table_list_headings, key="-BIO_INFO_LIST_TABLE-", auto_size_columns=False)],
        #         [sg.T("Method", size=14), sg.T("State", size=14), sg.T("Calculation", size=14)],
        #         [sg.DropDown([], key="-BIO_INFO_LIST_METHOD-", size=14, enable_events=True),
        #          sg.DropDown([], key="-BIO_INFO_LIST_STATE-", size=14, enable_events=True),
        #          sg.DropDown([], key="-BIO_INFO_LIST_CALC-", size=14, enable_events=True)],
        #     ])]
        # ])

        tab_mapping = sg.Tab("Mapping", [[row_options]])
        tab_overview = sg.Tab("Overview", [[row_overview]])
        tab_plate_overview = sg.Tab("Plate Overview", [[row_plate_overview]])
        tab_z_prime = sg.Tab("Z-Prime", [[row_z_prime]])
        tab_hit_list = sg.Tab("Hit List", [[row_hit_list]])
        tab_matrix = sg.Tab("Matrix", [[row_matrix]])
        # tab_list = sg.Tab("List", [[row_list]])

        tab_bio_list = [tab_mapping, tab_overview, tab_plate_overview, tab_z_prime, tab_hit_list, tab_matrix]

        tab_groups = [sg.TabGroup([tab_bio_list], selected_background_color=self.tab_colour,
                                  key="-BIO_INFO_SUB_SETTINGS_TABS-", enable_events=True)]

        top_row = [[row_settings], [row_graph]]

        layout = [[sg.Pane([sg.Column(top_row), sg.Column([tab_groups])])]]
        return layout

    def setup_2_purity(self):

        lc_graph_showing = [keys for keys in list(self.config["lc_mapping"].keys())]

        row_settings_col1 = sg.Frame("Purity settings", [[
                sg.Column([
                    [sg.DropDown(lc_graph_showing, key="-PURITY_INFO_GRAPH_SHOWING-", default_value=lc_graph_showing[0],
                                 enable_events=True)],
                    [sg.Radio("MS+", group_id="ms_mode", key="-PURITY_INFO_MS_MODE_POS-", default=True),
                     sg.Radio("MS-", group_id="ms_mode", key="-PURITY_INFO_MS_MODE_NEG-")],
                    [sg.T("Retion Time", size=10)],
                    [sg.T("Start:", size=3), sg.Input(key="-PURITY_INFO_RT_START-", size=5),
                     sg.T("End:", size=3), sg.Input(key="-PURITY_INFO_RT_END-", size=5)],
                    [sg.HorizontalSeparator()],
                    [sg.T("Wavelength:", size=10), sg.Input(key="-PURITY_INFO_WAVELENGTH-", size=10)],
                    [sg.T("bin:", size=10), sg.Input(key="-PURITY_INFO_BIN-", size=10)],
                    [sg.T("Mass:", size=10), sg.Input(key="-PURITY_INFO_MZ-", size=10)],
                    [sg.Push(), sg.B("Draw", key="-PURITY_INFO_DRAW_STUFF-")]
                    ])
        ]])

        row_settings_col2 = sg.Frame("Calculations", [[
            sg.Column([
                [sg.T("UV wavelength", size=self.standard_size),
                 sg.InputText(key="-PURITY_INFO_UV_WAVE-",
                              default_text=self.config["MS_default"]["uv_wavelength"], size=10)],
                [sg.HorizontalSeparator()],
                [sg.Text("UV threshold", size=self.standard_size),
                 sg.InputText(key="-PURITY_INFO_UV_THRESHOLD-",
                              default_text=int(self.config["MS_default"]["uv_threshold"]), size=10)],
                [sg.T("Slope Threshold", size=self.standard_size),
                 sg.InputText(key="-PURITY_INFO_SLOPE_THRESHOLD-",
                              default_text=int(self.config["MS_default"]["slop_threshold"]), size=10)],
                [sg.Text("Solvent peak retention time", size=self.standard_size),
                 sg.InputText(key="-PURITY_INFO_RT_SOLVENT-",
                              default_text=float(self.config["MS_default"]["rt_solvent"]), size=10)],
                [sg.HorizontalSeparator()],
                [sg.Text("Delta MS", size=self.standard_size),
                 sg.InputText(key="-PURITY_INFO_MS_DELTA-",
                              default_text=float(self.config["MS_default"]["ms_delta"]), size=10)],
                [sg.Text("MS threshold", size=self.standard_size),
                 sg.InputText(key="-PURITY_INFO_MS_THRESHOLD-",
                              default_text=int(self.config["MS_default"]["ms_threshold"]), size=10)],
                [sg.Text("MS peak amounts", size=self.standard_size),
                 sg.InputText(key="-PURITY_INFO_MS_PEAKS-",
                              default_text=int(self.config["MS_default"]["ms_peak_amount"]), size=10)],
                [sg.HorizontalSeparator()],
                [sg.Checkbox("Draw Threshold", key="-PURITY_INFO_DRAW_THRESHOLD-")],
                [sg.Checkbox("Use MS-data?", key="-PURITY_INFO_USE_MS_DATA-"),
                 sg.Push(), sg.Button("Re-calculate", key="-PURITY_INFO_RE_CALC-")]
            ])
        ]])

        row_settings_col3 = sg.Frame("Sample Selection", [[
            sg.Column([
                [sg.T("Batch", size=8), sg.T("Samples", size=8)],
                [sg.Listbox("", key="-PURITY_INFO_BATCH_BOX-", enable_events=True, size=(10, 7),
                            select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE),
                 sg.Listbox("", key="-PURITY_INFO_SAMPLE_BOX-", enable_events=True, size=(10, 7),
                            select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE)],
                [sg.Checkbox("Multi Sample selection", key="-PURITY_INFO_SAMPLE_SELECTION-", default=True,
                             enable_events=True)]
            ])
        ]])

        row_canvas = sg.Frame("Canvas", [[
            sg.Column([
                [sg.Canvas(key="-PURITY_INFO_CANVAS_TOOLBAR-")],
                [sg.Canvas(key="-PURITY_INFO_CANVAS-", size=(700, 200))]
            ])
        ]])


        # overview_headings = ["row_id", "sample", "batch"]
        # col_overview_table = sg.Frame("Overview", [[
        #     sg.Column([
        #         [sg.Table("", key="-PURITY_INFO_OVERVIEW_TABLE-", headings=overview_headings, enable_click_events=True,
        #                   enable_events=True, expand_x=True, auto_size_columns=False, col_widths=25)]
        #     ])
        # ]])

        purity_overview_table_headings = ["Sample", "batch", "Mass", "Major Peak %", "Ion", "Ion mass", "Peak ID",
                                          "Pur Sum"]
        table_tab_purity_overview = sg.Frame("Purity Overview", [[
            # sg.vbottom(
            sg.Column([
                [sg.Table("", key="-PURITY_INFO_PURITY_OVERVIEW_TABLE-", headings=purity_overview_table_headings,
                             enable_events=True, expand_x=True, auto_size_columns=False, col_widths=25
                             , enable_click_events=True)]
            ]),
            sg.Column([
                [sg.Button("Import", key="-PURITY_INFO_PURITY_OVERVIEW_IMPORT-")],
                [sg.Button("Report", key="-PURITY_INFO_PURITY_OVERVIEW_REPORT-")]
            ])
        # )
        ]])

        peak_table_headings = ["sample", "peak", "integrals", "start", "end", "%"]
        table_tab_peak = sg.Frame("Peaks", [[
            sg.Column([
                [sg.Table("", key="-PURITY_INFO_PEAK_TABLE-", headings=peak_table_headings, expand_x=True,
                          auto_size_columns=False, col_widths=25, enable_click_events=True, enable_events=True)]
            ]),
            sg.Column([
                [sg.B("Draw All Peaks", key="-PURITY_INFO_DRAW_PEAKS-")],
                [sg.Radio("UV", group_id="PURITY_INFO_RADIO_PEAKS", key="-PURITY_INFO_RADIO_PEAKS_UV-", default=True),
                 sg.Radio("MS Spectra", group_id="PURITY_INFO_RADIO_PEAKS", key="-PURITY_INFO_RADIO_PEAKS_MS_SPECTRA-")],
            ])
        ]])

        purity_peak_list_table_headings = ["Peak", "Ion", "Mass", "Purity", "Start", "End"]
        table_tab_purity_peak_list = sg.Frame("Purity Peak List", [[
            sg.Column([
                [sg.T("Samples:", size=10), sg.T("", key="-PURITY_INFO_PEAK_LIST_SAMPLE_TEXT-"),
                 sg.Input("", key="-PURITY_INFO_PEAK_LIST_SAMPLE-", visible=False, size=0)],
                [sg.Table("", key="-PURITY_INFO_PURITY_PEAK_LIST_TABLE-", headings=purity_peak_list_table_headings,
                          expand_x=True, auto_size_columns=False, col_widths=25, enable_click_events=True,
                          enable_events=True)]
            ])
        ]])

        raw_table_data_headings = ['Peak#', 'R.Time', 'I.Time', 'F.Time', 'Area', 'Height', 'A/H', 'Conc.', 'Mark',
                                   'ID#', 'Name', "k'", 'Plate #', 'Plate Ht.', 'Tailing', 'Resolution', 'Sep.Factor',
                                   'Area Ratio', 'Height Ratio', 'Conc. %', 'Norm Conc.']
        raw_table_data = sg.Frame("Raw Table", [[
            sg.Column([
                [sg.Table("", key="-PURITY_INFO_RAW_DATA_TABLE-", headings=raw_table_data_headings,
                          enable_events=True, expand_x=True, auto_size_columns=False, col_widths=25
                          , enable_click_events=True, vertical_scroll_only=False)]
            ])
        ]])

        tab_purity_overview = sg.Tab("Purity Overview", [[table_tab_purity_overview]])
        tab_peak = sg.Tab("Peaks", [[table_tab_peak]])
        tab_purity_peak_list = sg.Tab("Purity Peak List", [[table_tab_purity_peak_list]])

        tab_raw_table_data = sg.Tab("Raw Data table", [[raw_table_data]])

        tab_group_list = [tab_purity_overview, tab_peak, tab_purity_peak_list, tab_raw_table_data]

        table_tabs = sg.TabGroup([tab_group_list], selected_background_color=self.tab_colour, key="-TAB_GROUP_PURITY_INFO-",
                                 enable_events=True, expand_x=True)

        layout = [sg.vtop([row_settings_col1, row_settings_col2, row_settings_col3]), [row_canvas], [table_tabs]]
        return layout

    # def setup_2_misc(self):
    #     headings = ["Name", "country", "AC"]
    #     tab_customers = sg.Frame("FUCK", [[
    #         sg.Column([
    #             [sg.Table(values=[], headings=headings)]
    #         ])
    #     ]])
    #
    #     tabgroup = [sg.Tab("Custemors", [[tab_customers]])]
    #
    #     tab_groups = sg.TabGroup([tabgroup], selected_background_color=self.tab_colour, key="-COMPOUND_INFO_TABLE_TABS-",
    #                              enable_events=True, expand_x=True, tab_location="righttop")
    #
    #     layout = [[tab_groups]]
    #
    #     return layout

    def setup_table_compound(self):
        """

        :return: A layout for the compound table-module in the table box
        :rtype: list
        """
        headlines = ["compound_id", "smiles", "volume"]
        tables = ["Compound", "Mother Plates", "Assay Plates"]
        treedata = sg.TreeData()

        raw_table_col = sg.Column([
            [sg.Text("Raw data")],
            [sg.Tree(data=treedata, headings=headlines, row_height=90, auto_size_columns=False, num_rows=4,
                     col0_width=30, key="-TREE_DB-", show_expanded=True, expand_x=True,
                     enable_events=True)]
        ])

        layout = [
            [raw_table_col],
            [sg.Button("Export", key="-C_TABLE_EXPORT-", size=self.standard_size),
             sg.Button("Refresh", key="-C_TABLE_REFRESH-", size=self.standard_size),
             # sg.DropDown(values=tables, default_value=tables[0], key="-C_TABLE_FILE_TYPE-"),
             sg.Text(text="Compounds: 0", key="-C_TABLE_COUNT-")]
        ]
        return layout
        #return sg.Window("batch and sample selection", layout, size, finalize=True)

    def setup_table_bio_experiment(self):
        """

        :return: A layout for the compound experiment-module in the table box
        :rtype: list
        """
        responsible = [keys for keys in list(self.config["Responsible"].keys())]

        compound_table_data = []
        compound_table_headings = ["Compound ID", "well", "Raw Data", "Concentration", "Score", "Approved"]
        col_compound_table = sg.Frame("Experimental Data", [[
            sg.Column([
                [sg.Table(values=compound_table_data, headings=compound_table_headings, key="-BIO_EXP_COMPOUND_TABLE-",
                          auto_size_columns=False, col_widths=[10, 10, 10], enable_events=True,
                          enable_click_events=True)],
                [sg.T("Amount of Compounds: "), sg.T("0", key="-BIO_EXP_COMPOUND_COUNTER-")],
                [sg.T("Set Threshold:", size=12), sg.Input(key="-BIO_EXP_SET_THRESHOLD-", size=12)],
                [sg.T("Set Compound Amount:", size=12), sg.Input(key="-BIO_EXP_SET_COMPOUND_AMOUNT-", size=12)],
                [sg.Button("Export Compound List", key="-BIO_EXP_EXPORT_COMPOUNDS-")]
            ])
        ]])

        plate_table_data = []
        plate_table_headings = ["Plate Name", "Z-Prime", "Approved", "Layout", "responsible", "date"]

        col_plate_table = sg.Frame("Plates", [[
            sg.Column([
                [sg.Table(values=plate_table_data, headings=plate_table_headings, key="-BIO_EXP_PLATE_TABLE-",
                          auto_size_columns=False, col_widths=[10, 10, 10],
                          enable_events=True, enable_click_events=True)],
                [sg.T("Amount of Plates: "), sg.T("0", key="-BIO_EXP_PLATE_COUNTER-")],
                [sg.Checkbox("Only approved Plates", key="-BIO_EXP_APPROVED_PLATES_ONLY-", default=True,
                             enable_events=True)]
            ])
        ]])

        col_search = sg.Frame("Search Criteria", [[
            sg.Column([
                [sg.Listbox("", key="-BIO_EXP_TABLE_ASSAY_LIST_BOX-", enable_events=True, size=(18, 10),
                            select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE)],
                [sg.T("Assay Name", size=10), sg.Input(key="-BIO_EXP_ASSAY_NAME_SEARCH-", size=10)],
                [sg.CalendarButton("Start date", key="-BIO_EXP_TABLE_DATE_START-", format="%Y-%m-%d", enable_events=True
                                   , target="-BIO_EXP_TABLE_DATE_START_TARGET-", size=(10, 1)),
                 sg.Input(key="-BIO_EXP_TABLE_DATE_START_TARGET-", size=10, enable_events=True)],
                [sg.CalendarButton("End date", key="-BIO_EXP_TABLE_DATE_END-", format="%Y-%m-%d", enable_events=True,
                                   target="-BIO_EXP_TABLE_DATE_END_TARGET-", size=(10, 1)),
                 sg.Input(key="-BIO_EXP_TABLE_DATE_END_TARGET-", size=10)],
                [sg.T("Responsible", size=10), sg.DropDown(responsible, key="-BIO_EXP_TABLE_RESPONSIBLE-", size=10)],
                [sg.B("Refresh", key="-BIO_EXP_TABLE_REFRESH-")]
            ])
        ]])

        layout = [sg.vtop([col_compound_table, col_plate_table, col_search])]
        return layout

    @staticmethod
    def setup_table_lc_experiment():
        """

        :return: A layout for the compound experiment-module in the table box
        :rtype: list
        """

        table_data = []
        headings = ["exp_id", "assay_name", "raw_data", "plate_layout", "responsible", "date"]

        col_table = sg.Frame("LC Samples", [[
            sg.Column([
                [sg.Table(values=table_data, headings=headings, key="-LC_MS_SAMPLE_TABLE-",
                          auto_size_columns=False, col_widths=[10, 10, 10], size=(20, 20), enable_events=True,
                          enable_click_events=True)],
                # [sg.B("Refresh", key="-LC_MS_SAMPLE_TABLE_REFRESH-")]

        ])]])

        col_search = sg.Frame("Search Criteria", [[
            sg.Column([
                [sg.Listbox("", key="-LC_MS_TABLE_BATCH_LIST_BOX-", enable_events=True, size=(18, 10),
                            select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE)],
                [sg.CalendarButton("Start date", key="-LC_MS_TABLE_DATE_START-", format="%Y-%m-%d", enable_events=True
                                   , target="-LC_MS_TABLE_DATE_START_TARGET-", size=(10, 1)),
                 sg.Input(key="-LC_MS_TABLE_DATE_START_TARGET-", size=10, enable_events=True)],
                [sg.CalendarButton("End date", key="-LC_MS_TABLE_DATE_END-", format="%Y-%m-%d", enable_events=True,
                                   target="-LC_MS_TABLE_DATE_END_TARGET-", size=(10, 1)),
                 sg.Input(key="-LC_MS_TABLE_DATE_END_TARGET-", size=10)]

            ])
        ]])

        layout = [sg.vtop([col_table, col_search])]
        return layout

    @staticmethod
    def setup_table_plate():
        """

        :return: A layout for the plate table-module in the table box
        :rtype: list
        """
        headings = ["row_counter", "Barcode", "Compound", "Well", "Volume", "Date", "Active", "Freeze/Thaw",
                    "Plate Type", "location", "Source Plate", "Source Well"]
        # headings = ["row_counter", "Barcode", "Compound", "Well", "Volume", "Date", ,
        #             "Freeze/Thaw", "Active", "Plate Type", "location"]
        dd = ["Mother Plates", "Daughter Plates"]
        col_table = sg.Frame("Plate Table", [[
            sg.Column([
                [sg.Table([], headings=headings, key="-PLATE_TABLE_TABLE-", size=(20, 20), enable_click_events=True)]
                # [sg.Input("", key="-PLATE_TABLE_BARCODE_SEARCH-", size=10), sg.B("Search", key="-PLATE_TABLE_SEARCH-")]
            ])]])

        col_search = sg.Frame("Seach Perameters", [[
            sg.Column([
                [sg.Listbox("", key="-PLATE_TABLE_BARCODE_LIST_BOX-", enable_events=True, size=(15, 7),
                            select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE)],
                [sg.DropDown(dd, key="-PLATE_TABLE_CHOOSER-", enable_events=True, default_value=dd[0])],
                [sg.Button("Limit", key="-PLATE_TABLE_BUTTON_LIMITER-", size=10),
                 sg.Input(key="-PLATE_TABLE_TEXT_LIMITER-", size=10)],
                [sg.CalendarButton("Start Date", key="-PLATE_TABLE_START_DATE-", format="%Y-%m-%d", size=10,
                                   enable_events=True, target="-PLATE_TABLE_START_DATE_TARGET-"),
                 sg.Input(key="-PLATE_TABLE_START_DATE_TARGET-", size=10, enable_events=True)],
                [sg.CalendarButton("End Date", key="-PLATE_TABLE_END_DATE-", format="%Y-%m-%d", size=10,
                                   enable_events=True, target="-PLATE_TABLE_END_DATE_TARGET-"),
                 sg.Input(key="-PLATE_TABLE_END_DATE_TARGET-", size=10, enable_events=True)],
                [sg.Push(), sg.B("Clear", key="-PLATE_TABLE_CLEAR-", size=10)]
            ])
        ]], expand_x=True)

        layout = [sg.vtop([col_table, col_search])]
        return layout

    # def setup_table_customers(self):
    #     headings = ["Name", "country", "AC"]
    #     row = sg.Frame("Curstomers", [[
    #         sg.Column([
    #             [sg.Table(values=[], headings=headings)]
    #         ])
    #     ]])
    #
    #     layout = [[row]]
    #
    #     return layout

    def layout_tab_group_1(self):
        """

        :return: the layout for the tab groups in the top box
        :rtype: list
        """
        tab_1_search = sg.Tab("Search", self.setup_1_search(),
                              tooltip="Search for compounds in the database")
        tab_1_bio_data = sg.Tab("Bio Data", self.setup_1_bio(),
                                tooltip="Handles analysing of bio-data from platereader")
        tab_1_purity_data = sg.Tab("Purity Data", self.setup_1_purity(),
                                   tooltip="Handles data from LC/MS to get purity of compounds")
        tab_1_plate_layout = sg.Tab("Plate Layout", self.setup_1_plate_layout(),
                                    tooltip="Where you can draw a plate-layout used for analysing data, "
                                            "or generate worklist ect.")
        tab_1_add = sg.Tab("Update", self.setup_1_update(),
                           tooltip="Updates the Database with new compounds, Mother-Plates or Daugther-Plates")
        tab_1_worklist = sg.Tab("Worklist", self.set_1_worklist(),
                                tooltip="A module for generating worklist for en Echo.")
        tab_1_extra = sg.Tab("Extra", self.setup_1_extra(),
                             tooltip="Updates the database with other stuff, and misc")
        tab_1_sim = sg.Tab("Sim", self.setup_1_simulator(),
                           tooltip="simulate output files from different aspect of the HTS for testing purpose")

        tab_group_1_list = [tab_1_search, tab_1_bio_data, tab_1_purity_data, tab_1_plate_layout, tab_1_add,
                            tab_1_worklist, tab_1_extra, tab_1_sim]

        return [[sg.TabGroup([tab_group_1_list], selected_background_color=self.tab_colour, key="-TAB_GROUP_ONE-",
                             enable_events=True)]]

    def layout_tab_group_2(self):
        """

        :return: the layout for the tab groups in the right box
        :rtype: list
        """

        tab_2_info = sg.Tab("Compound Info", self.setup_2_compound())
        tab_2_bio_bio = sg.Tab("bio info", self.setup_2_bio())
        tab_2_purity_purity = sg.Tab("purity info", self.setup_2_purity())
        # tab_2_customers = sg.Tab("Misc", self.setup_2_misc())

        tab_group_2_list = [tab_2_info, tab_2_bio_bio, tab_2_purity_purity]

        return [[sg.TabGroup([tab_group_2_list], selected_background_color=self.tab_colour, key="-TAB_GROUP_TWO-",
                             enable_events=True)]]

    def layout_tab_group_tables(self):
        """

        :return: the layout for the tab groups in the table box
        :rtype: list
        """

        tab_table_compound = sg.Tab("Compound table", self.setup_table_compound())
        tab_bio_experiment_table = sg.Tab("Bio Experiment table", self.setup_table_bio_experiment())
        tab_lc_experiment_table = sg.Tab("LC Experiment table", self.setup_table_lc_experiment())
        tab_plate_table = sg.Tab("Plate tables", self.setup_table_plate())
        # tab_customers_table = sg.Tab("Customers", self.setup_table_customers())

        tab_group_tables = [tab_table_compound, tab_bio_experiment_table, tab_lc_experiment_table, tab_plate_table]
        return [[sg.TabGroup([tab_group_tables], key="-TABLE_TAB_GRP-", enable_events=True,
                             selected_background_color=self.tab_colour)]]

    def full_layout(self):
        """

        :return: The final layout
        :rtype: PySimpleGUI.PySimpleGUI.Window
        """
        sg.theme(self.config["GUI"]["theme"])
        window_size = (int(self.config["GUI"]["size_x"]), int(self.config["GUI"]["size_y"]))
        x_size = window_size[0]/3
        y_size = window_size[1]/3

        tab_group_1 = self.layout_tab_group_1()
        tab_group_2 = self.layout_tab_group_2()
        table_block = self.layout_tab_group_tables()
        menu = self.menu_top()
        mouse_right_click, right_click_options = self.menu_mouse()
        # col_1_1 = [[sg.Frame(layout=tab_group_1, title="X-1", size=(x_size*2, y_size))]]
        # col_1_2 = [[sg.Frame(layout=table_block, title="Table", size=(x_size*2, y_size*2))]]
        # col_2 = [[sg.Frame(layout=tab_group_2, title="X-2", size=(x_size, window_size[1]))]]

        layout_complete = [[
            menu,
            sg.Pane([
                sg.Column(
                    [[sg.Pane([sg.Column(tab_group_1), sg.Column(table_block)], orientation="v")]]
                ),
                sg.Column(tab_group_2)
                ], orientation="h")
        ]]

        return sg.Window("SCore", layout_complete, finalize=True, resizable=True, right_click_menu=right_click_options)
