import configparser
from database_handler import DataBaseFunctions
from pathlib import Path
import xml.etree.ElementTree as ET
from csv_handler import CSVWriter
from file_handler import get_file_list
from gui_functions import _folder_digger_for_file_names
from gui_popup import plate_layout_chooser


def finding_duplicate_compounds_in_mps(config):
    dbf = DataBaseFunctions(config)


    table = "compound_mp"
    table_string = f"SELECT * FROM {table}"
    row = dbf._row_creator(table_string)
    all_compounds = []
    duplicate_compounds = []
    dup_dup_compounds = {}
    dup_plates = {}
    for data in row:
        compound = row[data]["compound_id"]
        if compound not in all_compounds:
            all_compounds.append(compound)
        else:
            if compound in duplicate_compounds:
                try:
                    dup_dup_compounds[compound]
                except KeyError:
                    dup_dup_compounds[compound] = 1
                else:
                    dup_dup_compounds[compound] = dup_dup_compounds[compound] + 1

            duplicate_compounds.append(compound)
            plate = row[data]["mp_barcode"]
            try:
                dup_plates[plate]
            except KeyError:
                dup_plates[plate] = 1
            else:
                dup_plates[plate] = dup_plates[plate] + 1



    print(len(duplicate_compounds))
    print(duplicate_compounds)
    print(dup_plates)
    print(len(dup_plates))
    print(dup_dup_compounds)

def digger(folder):
    p = Path(folder).glob("**/*")
    files = [x for x in p if x.is_file()]
    used_tubes = []
    for file_path in files:
        with open(file_path) as f:
            temp_lines = f.readlines()

            for line in temp_lines:
                tube_id = line.split(";")[-1]
                tube_id = tube_id.removesuffix("\n")
                used_tubes.append(tube_id)

    return used_tubes


def compare(used_tubes, missing_tubes, folder):


    new_list = []
    bah = []
    with open(missing_tubes) as f:
        lines = f.readlines()
        for tube_id in lines:
            tube_id = tube_id.removesuffix("\n").strip("")
            if tube_id not in used_tubes:

                new_list.append(tube_id)
            else:
                bah.append(tube_id)

    CSVWriter.compound_freezer_writer(folder, new_list)


def plate_layout_check(folder, default_plate_layout, all_plate_layouts):

    # file_list = _folder_digger_for_file_names(folder)
    #
    # plate_layout_chooser(file_list, default_plate_layout, all_plate_layouts)
    file_list = get_file_list(folder)
    for files in file_list:
        print(files)


if __name__ == "__main__":
    # config = configparser.ConfigParser()
    # config.read("config.ini")
    folder = r"C:\Users\phch\Desktop\test\test_data\export"
    default_plate_layout = "Place_Holder"
    all_plate_layouts = ["Place_holder", "something_else"]
    plate_layout_check(folder, default_plate_layout, all_plate_layouts)
