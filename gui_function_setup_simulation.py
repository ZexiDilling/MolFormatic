from PySimpleGUI import PopupError, Popup

from csv_handler import CSVWriter, CSVReader, CSVConverter
from file_handler import get_file_list


def simulation_input_update(window, values):
    if values["-SIM_INPUT_EQ-"] == "comPOUND":
        window["-SIM_COMPOUND_FRAME-"].update(visible=True)
        window["-SIM_MP_FRAME-"].update(visible=False)
        window["-SIM_DP_FRAME-"].update(visible=False)

    elif values["-SIM_INPUT_EQ-"] == "MP Production":
        window["-SIM_COMPOUND_FRAME-"].update(visible=False)
        window["-SIM_MP_FRAME-"].update(visible=True)
        window["-SIM_DP_FRAME-"].update(visible=False)

    elif values["-SIM_INPUT_EQ-"] == "DP production":
        window["-SIM_COMPOUND_FRAME-"].update(visible=False)
        window["-SIM_MP_FRAME-"].update(visible=False)
        window["-SIM_DP_FRAME-"].update(visible=True)


def simulation_run(window, values):
    if values["-SIM_INPUT_EQ-"] == "comPOUND":
        if not values["-SIM_INPUT_COMPOUND_FILE-"]:
            PopupError("Missing Compound file")
        else:
            tube_file = values["-SIM_INPUT_COMPOUND_FILE-"]
            output_folder = values["-SIM_OUTPUT-"]

            _compound_freezer_to_2d_simulate(tube_file, output_folder)

            Popup("Done")

    elif values["-SIM_INPUT_EQ-"] == "MP Production":
        if not values["-SIM_INPUT_MP_FILE-"]:
            PopupError("Missing 2D barcode file")
        else:
            output_folder = values["-SIM_OUTPUT-"]
            barcodes_2d = values["-SIM_INPUT_MP_FILE-"]
            mp_name = values["-SIM_MP_NAME-"]
            trans_vol = values["-SIM_MP_VOL-"]

            _mp_production_2d_to_pb_simulate(output_folder, barcodes_2d, mp_name, trans_vol)

            Popup("Done")

    elif values["-SIM_INPUT_EQ-"] == "DP production":
        if not values["-SIM_INPUT_DP_FILE-"]:
            PopupError("Missing PlateButler file")
        else:
            Popup("not working atm")

    else:
        print(f"SIM_INPUT_EQ: {values['-SIM_INPUT_EQ-']}")


def _compound_freezer_to_2d_simulate(tube_file, output_folder):
    """
    A simulation modul, for simulating output data

    :param tube_file: The CSV file with all the tube ID's for the comPOUND freezer
    :type tube_file: str
    :param output_folder: The output folder for the CSV files
    :type output_folder: str
    :return: A lot of CSV files. 1 per 96 compounds
    """
    csv_w = CSVWriter()
    csv_r = CSVReader
    tube_dict = csv_r.tube_list_to_dict(tube_file)
    csv_w.compound_freezer_to_2d_csv_simulate(tube_dict, output_folder)


def _mp_production_2d_to_pb_simulate(folder_output, barcodes_2d, mp_name, trans_vol):
    """
    A simulation modul, for simulating output data

    :param folder_output: Output folder
    :type folder_output: str
    :param barcodes_2d: The folder with the 2D barcodes
    :type barcodes_2d: str
    :param mp_name: The name used for the MotherPlates
    :type mp_name: str
    :param trans_vol: The amount of volume to transferee
    :type trans_vol: float
    :return: CSV file resembling the one produced by the PB
    """
    barcodes_2d_files = get_file_list(barcodes_2d)
    csvc = CSVConverter()
    csvc.mp_in_to_out(folder_output, barcodes_2d_files, mp_name, trans_vol)