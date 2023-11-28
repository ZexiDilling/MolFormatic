from pathlib import Path

from natsort import natsorted
from PySimpleGUI import PopupGetFile, PopupError, Popup, PopupYesNo

from csv_handler import CSVWriter
from database_handler import DataBaseFunctions

from excel_handler import well_compound_list
from bio_functions import bio_compound_info_from_worklist
from database_functions import grab_table_data


def worklist_control_layout_update(window, values):
    worklist_layout = PopupGetFile("Please select a worklist layout file")
    worklist_data = well_compound_list(worklist_layout)
    for well_state in worklist_data:
        well_state = well_state.casefold()
        if well_state == "positive":
            window["-WORKLIST_USE_POSITIVE_CONTROL-"].update(value=True)
            window["-WORKLIST_POSITIVE_CONTROL_ID-"].update(disabled=False,
                                                            value=worklist_data[well_state]["compound"])
        if well_state == "negative":
            window["-WORKLIST_USE_NEGATIVE_CONTROL-"].update(value=True)
            window["-WORKLIST_NEGATIVE_CONTROL_ID-"].update(disabled=False,
                                                            value=worklist_data[well_state]["compound"])
        if well_state == "bonus":
            window["-WORKLIST_USE_BONUS_COMPOUND-"].update(value=True)
            window["-WORKLIST_BONUS_COMPOUND_ID-"].update(disabled=False,
                                                          value=worklist_data[well_state]["compound"])

    window["-WORKLIST_CONTROL_LAYOUT_TARGET-"].update(value=worklist_layout)


def worklist_tab_clicked(config, window, values):
    if values["-TAB_GROUP_ONE-"] == "Worklist":
        temp_mp_plates, _ = grab_table_data(config, "mp_plates")
        worklist_mp_plates_list = []
        for rows in temp_mp_plates:
            worklist_mp_plates_list.append(rows[0])

        # sortes the table
        worklist_mp_plates_list = natsorted(worklist_mp_plates_list)

        window["-WORKLIST_MP_LIST-"].update(values=worklist_mp_plates_list)
        # window["-WORKLIST_ASSAY_LIST-"].update(values=worklist_mp_plates_list)    # ToDO add the right data here

        temp_assay_list, _ = grab_table_data(config, "assay")
        return temp_assay_list, worklist_mp_plates_list
    else:
        return None, None


def worklist_generator(config, window, values, worklist_mp_plates_list, archive_plates_dict):
    if not values["-WORKLIST_PLATE_LAYOUT-"]:
        PopupError("Please select a plate layout")

    elif not values["-WORKLIST_MP_LIST-"] and not values["-WORKLIST_USE_ALL_MOTHERPLATES-"]:
        PopupError("Please select witch MotherPlates to use")

    elif not values["-WORKLIST_PLATE_AMOUNT-"]:
        PopupError("Please write how many plates are needed")

    elif not values["-WORKLIST_INITIAL_PLATE-"]:
        PopupError("Please write what the starting number is")

    elif not values["-WORKLIST_VOLUME-"]:
        PopupError("Please write how much volume is needed per plate")

    elif values["-WORKLIST_USE_POSITIVE_CONTROL-"] and not values["-WORKLIST_POSITIVE_CONTROL_ID-"]:
        PopupError("Please write an ID for the Positive control.")

    elif values["-WORKLIST_USE_NEGATIVE_CONTROL-"] and not values["-WORKLIST_NEGATIVE_CONTROL_ID-"]:
        PopupError("Please write an ID for the Negative control.")

    elif values["-WORKLIST_USE_POSITIVE_CONTROL-"] and not values["-WORKLIST_CONTROL_LAYOUT_TARGET-"]:
        PopupError("Please select a Plate Layout for the Positive control.")

    elif values["-WORKLIST_USE_NEGATIVE_CONTROL-"] and not values["-WORKLIST_CONTROL_LAYOUT_TARGET-"]:
        PopupError("Please select a Plate Layout for the Negative control.")

    elif values["-WORKLIST_USE_BONUS_COMPOUND-"] and not values["-WORKLIST_BONUS_COMPOUND_ID-"]:
        PopupError("Please write an ID for the Bonus Compound.")

    elif values["-WORKLIST_USE_BONUS_COMPOUND-"] and not values["-WORKLIST_CONTROL_LAYOUT_TARGET-"]:
        PopupError("Please select a Plate Layout for the Negative control.")
    # ToDo make this one work
    # elif values["-WORKLIST_USE_BONUS_COMPOUND-"] and not values["-WORKLIST_BONUS_MAX-"]\
    #     or values["-WORKLIST_USE_BONUS_COMPOUND-"] and not values["-WORKLIST_BONUS_POSITIVE-"]\
    #     or values["-WORKLIST_USE_BONUS_COMPOUND-"] and not values["-WORKLIST_BONUS_EMPTY-"]\
    #     or values["-WORKLIST_USE_BONUS_COMPOUND-"] and not values["-WORKLIST_BONUS_MIN-"]\
    #     or values["-WORKLIST_USE_BONUS_COMPOUND-"] and not values["-WORKLIST_BONUS_NEGATIVE-"]\
    #     or values["-WORKLIST_USE_BONUS_COMPOUND-"] and not values["-WORKLIST_BONUS_BLANK-"]:
    #         PopupError("Please select minimum one Well states where the Bonus compounds should be added")
    else:
        if values["-WORKLIST_USE_ALL_MOTHERPLATES-"]:
            # This is a list of all MotherPlates. It is generated when the tab for "Worklist" is clicked
            mps = worklist_mp_plates_list
        else:
            mps = values["-WORKLIST_MP_LIST-"]

        if not values["-WORKLIST_ASSAY_LIST-"]:
            assays = None
            mp_check = PopupYesNo("Have any of the MotherPlates been used for for this production before?")
            if mp_check.casefold() == "yes":
                worklist = PopupGetFile("Please select worklist files", multiple_files=True)
                if worklist:
                    worklist = worklist.split(";")
                else:
                    worklist = "cancelled"
            else:
                worklist = None
        else:
            worklist = None
            assays = values["-WORKLIST_ASSAY_LIST-"]
        if worklist != "cancelled":
            plate_layout = archive_plates_dict[values["-WORKLIST_PLATE_LAYOUT-"]]

            if worklist:
                # Get data from the worklist, to see what plate and witch wells have been used before
                used_plate_well_dict = bio_compound_info_from_worklist(config, worklist)
            elif assays:
                print("Find assay data")  # ToDo make this work
                used_plate_well_dict = None
            else:
                used_plate_well_dict = None

            if values["-WORKLIST_USE_POSITIVE_CONTROL-"] or values["-WORKLIST_USE_NEGATIVE_CONTROL-"] or \
                    values["-WORKLIST_USE_BONUS_COMPOUND-"]:
                use_control = True
            else:
                use_control = False

            if use_control:
                control_layout = Path(values["-WORKLIST_CONTROL_LAYOUT_TARGET-"])

            else:
                control_layout = None

            control_samples = {"positive":
                                   {"use": values["-WORKLIST_USE_POSITIVE_CONTROL-"],
                                    "sample": values["-WORKLIST_POSITIVE_CONTROL_ID-"]},
                               "negative":
                                   {"use": values["-WORKLIST_USE_NEGATIVE_CONTROL-"],
                                    "sample": values["-WORKLIST_NEGATIVE_CONTROL_ID-"]},
                               "max": {"use": False},
                               "minimum": {"use": False},
                               "blank": {"use": False},
                               "empty": {"use": False},
                               "sample": {"use": False}
                               }

            bonus_compound = {"sample_name": values["-WORKLIST_BONUS_COMPOUND_ID-"].casefold(),
                              "max": values["-WORKLIST_BONUS_MAX-"],
                              "minimum": values["-WORKLIST_BONUS_MIN-"],
                              "positive": values["-WORKLIST_BONUS_POSITIVE-"],
                              "negative": values["-WORKLIST_BONUS_NEGATIVE-"],
                              "blank": values["-WORKLIST_BONUS_BLANK-"],
                              "empty": values["-WORKLIST_BONUS_EMPTY-"],
                              "sample": values["-WORKLIST_BONUS_SAMPLE-"]}

            worklist_analyse_method = values["-WORKLIST_SAMPLE_STYLE-"]
            sample_direction = values["-WORKLIST_DROPDOWN_SAMPLE_DIRECTION-"]

            assay_name = values["-WORKLIST_ASSAY_NAME-"]
            plate_amount = int(values["-WORKLIST_PLATE_AMOUNT-"])
            initial_plate = int(values["-WORKLIST_INITIAL_PLATE-"])
            volume = float(values["-WORKLIST_VOLUME-"])
            file, msg = generate_worklist(config, plate_amount, mps, plate_layout, used_plate_well_dict,
                                          assay_name, initial_plate, volume, worklist_analyse_method,
                                          sample_direction, control_layout, control_samples,
                                          bonus_compound)

            if not file:
                PopupError("Something crashed up")
            elif type(msg) == str:
                Popup(f"{msg} - File still created with fewer plates, saved here {file}")
            else:
                Popup(f"Worklist have been created and saved here: {file}")


def _get_motherplates_with_wells_from_worklist_dict(used_plate_well_dict):
    """
    Makes a new dict, to get MotherPlates
    :param used_plate_well_dict: A dicts of MotherPlates and the wells that have been used
    :type used_plate_well_dict: dict
    :return: motherplate_dict, a dict with MotherPlates as Keys and wells used as a list as value
    :rtype: dict
    """

    motherplate_dict = {}
    for destinations_plates in used_plate_well_dict:

        for wells in used_plate_well_dict[destinations_plates]:
            temp_mp = used_plate_well_dict[destinations_plates][wells]["source_plate"]
            try:
                motherplate_dict[temp_mp]
            except KeyError:
                motherplate_dict[temp_mp] = []

            source_well = used_plate_well_dict[destinations_plates][wells]["source_well"]
            if source_well not in motherplate_dict[temp_mp]:
                motherplate_dict[temp_mp].append(source_well)

    return motherplate_dict


def _get_motherplate_layout(config, mps):
    """
    Generate a dict with the key of each motherplate and the value for all the wells that have compounds.
    :param mps: A list of MotherPlates
    :type: mps: list
    :return: a dict with the key of each motherplate and the value for all the wells that have compounds.
    :rtype: dict
    """

    dbf = DataBaseFunctions(config)
    table = "compound_mp"
    clm_header = "mp_barcode"
    plate_data_dict = {}
    for mp in mps:
        plate_data_dict[mp] = []
        temp_rows = dbf.find_data_single_lookup(table, mp, clm_header)
        for data in temp_rows:
            temp_well = temp_rows[data]["mp_well"]
            plate_data_dict[mp].append(temp_well)

    return plate_data_dict


def _get_free_wells(mps_layout, motherplate_dict):
    """
    Compare MPS with previuse used plates, and finds the wells that have not been used and add them to a dict
    :param mps_layout: A dict with the layout for each motherplate
    :type: mps_layout: dict
    :param motherplate_dict: A dicts of MotherPlates and the wells that have been used
    :type motherplate_dict: dict
    :return: a dict with MotherPlates as keys, and list of free wells as values
    :rtype: dict
    """
    # This depends on the motherplate layout. We use a full 384 layout
    free_wells = {}
    for motherplates in mps_layout:
        free_wells[motherplates] = []
        for wells in mps_layout[motherplates]:

            if motherplate_dict and motherplates in motherplate_dict:
                if wells not in motherplate_dict[motherplates]:
                    free_wells[motherplates].append(wells)
            else:
                free_wells[motherplates].append(wells)
    return free_wells


def generate_worklist(config, plate_amount, mps, plate_layout, used_plate_well_dict, assay_name, initial_plate, volume,
                      worklist_analyse_method, sample_direction, control_layout, control_samples, bonus_compound):
    """
    Generates a worklist based on a plate-layout and data from MotherPlates.

    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param plate_amount: Number of plates to produce
    :type plate_amount: int
    :param mps: A list of MotherPlates
    :type: mps: list
    :param plate_layout: The layout for the plate with values for each well, what state they are in
    :type plate_layout: dict
    :param used_plate_well_dict: A dicts of MotherPlates and the wells that have been used
    :type used_plate_well_dict: dict
    :param assay_name: the name of the assay, and the name used for the destination plate in the workinglist
    :type assay_name: str
    :param initial_plate: the starting number of the plates
    :type initial_plate: int
    :param volume: How much volume to transfere to each well. The same amount of liquid will be transfered to each.
    :type volume: float
    :param worklist_analyse_method: The method use for the sample layout.
    :type worklist_analyse_method: str
    :param sample_direction: The direction for the sample layout. Only relevant if the analyse_method differ from
        "Single Point"
    :type sample_direction: str
    :param control_layout: The file path to the plate layout for the control or other compounds
    :type control_layout: pathlib.WindowsPath or None
    :param control_samples: A dict for samples for positive and negative control
    :type: control_samples: None or dict
    :param bonus_compound: If there are a compound that needs to be added to multiple well_states
    :type bonus_compound: dict
    :return:
    """

    # Gets a dict of motherplates that have been used, and what wells in each motherplate,
    # based on the worklist provided.
    if used_plate_well_dict:
        motherplate_dict = _get_motherplates_with_wells_from_worklist_dict(used_plate_well_dict)
    else:
        motherplate_dict = None

    # Get the layout of the motherplates. incase the motherplate is not a full 384 plate.
    mps_layout = _get_motherplate_layout(config, mps)

    # gets all the free wells, based on the layout of each motherplate, and what wells have been used before,
    # based on the worklist provided
    free_well_dict = _get_free_wells(mps_layout, motherplate_dict)

    if not control_layout:
        control_bonus_source = None
    elif control_layout.suffix == ".xlsx":
        control_bonus_source = well_compound_list(control_layout)
        if type(control_bonus_source) == str:
            return control_bonus_source
    elif control_layout.suffix == ".txt":
        pass            # ToDo add these ?
    elif control_layout.suffix == ".csv":
        pass            # ToDo add these ?
    else:
        return "error: wrong file format for Control Layout"
    # destination_dict = _from_plate_layout_to_destination_dict(plate_layout, assay_name, plate_amount, initial_plate)
    csv_w = CSVWriter()
    file, msg = csv_w.worklist_writer(config, plate_layout, mps, free_well_dict, assay_name, plate_amount,
                                      initial_plate, volume, sample_direction, worklist_analyse_method,
                                      control_bonus_source, control_samples, bonus_compound)

    if type(msg) == str:
        print(msg)
    return file, msg