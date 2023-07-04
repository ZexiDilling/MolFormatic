import configparser
from openpyxl.styles import PatternFill, Font
from openpyxl import Workbook, load_workbook
from pathlib import Path
import csv

from gui_functions import _get_motherplate_layout


def file_reader(files):

    missing_tubes = []
    file = Path(r"C:\Users\phch\Desktop\test\comPOUND\comPOUND_'2023-06-22'.txt")
    # file = Path(r"D:\MP\TubeRackBarcodes\3000307768.txt")
    with open(file) as temp_csv:
        rows = csv.reader(temp_csv, delimiter=';')
        for row, row_data in enumerate(rows):
            # print(row_data[1])
            missing_tubes.append(row_data[0])
    print(len(missing_tubes))
    tube_location = {}
    files = Path(r"D:\MP\TubeRackBarcodes")
    # files = Path(r"C:\Users\phch\OneDrive - Danmarks Tekniske Universitet\Mapper\Python_data\MP Production\All_plates")
    # files = list(files.iterdir())
    # for file in files:
    #     with open(file) as temp_csv:
    #         rows = csv.reader(temp_csv, delimiter=';')
    #
    #         for row, row_data in enumerate(rows):
    #             if row_data[1] in missing_tubes:
    #
    #                 try:
    #                     tube_location[file]
    #                 except KeyError:
    #                     tube_location[file] = [row_data[1]]
    #                 else:
    #                     tube_location[file].append(row_data[1])
    #
    # print(tube_location)
    # for locations in tube_location:
    #     print(locations)
    #     print(len(tube_location[locations]))



def excel_reader(file, folder):


    files = list(folder.iterdir())
    # files = [folder]
    all_data = {}
    all_data_2 = {}
    for counter, temp_file in enumerate(files):
        if temp_file.suffix == ".csv":
            with open(temp_file) as temp_csv:
                rows = csv.reader(temp_csv, delimiter=';')
                for row, row_data in enumerate(rows):
                    if row != 0:
                        chemist_id = row_data[4]
                        formel = row_data[1]
                        mol_weight = row_data[2]
                        smiles = row_data[3]
                        comp_id = row_data[0]
                        try:
                            qc_plate = row_data[13]
                        except IndexError:
                            qc_plate = "None"

                        all_data[chemist_id] = {
                            "formel": formel,
                            "mol_weight": mol_weight,
                            "smiles": smiles,
                            "comp_id": comp_id,
                            "qc_plate": qc_plate
                        }
                        all_data_2[row_data[5]] = {
                            "formel": row_data[2],
                            "mol_weight": row_data[3],
                            "smiles": row_data[4],
                            "comp_id": row_data[1],
                        }



    wb_final = load_workbook(filename=file)
    all_sheets = wb_final.sheetnames

    lables = []
    non_duplicates = []
    duplicates = []
    found_id = []
    not_found = []

    # for sheets in all_sheets:
    counter = 0
    ws = wb_final["DK"]
    for row, data in enumerate(ws):

        for col, cells in enumerate(data):
            if row == 0:
                if cells.value == "compound-id":
                    compound_id_col = col
            else:
                if col == compound_id_col:
                    temp_comp_id = cells.value

                    if temp_comp_id not in lables:
                        lables.append(temp_comp_id)
                        try:
                            all_data[temp_comp_id]
                        except KeyError:
                            # print(f"The following comp Id could not be found: {temp_comp_id}")
                            try:

                                all_data_2[temp_comp_id]
                            except KeyError:
                                not_found.append(temp_comp_id)
                                # print(sheets)
                            else:
                                found_id.append(temp_comp_id)
                        else:
                            found_id.append(temp_comp_id)

                    else:
                        duplicates.append(temp_comp_id)

    for temp_id in lables:
        if temp_id not in duplicates:
            non_duplicates.append(temp_id)

    try:
        ws = wb_final["All_data"]
    except KeyError:
        ws = wb_final.create_sheet("All_data")

    row_counter = 1
    col_counter = 1

    for counter, compounds in enumerate(found_id):
        if counter == 0:

            ws.cell(column=col_counter, row=row_counter, value="Compound-id")
            ws.cell(column=col_counter + 1, row=row_counter, value="formel")
            ws.cell(column=col_counter + 2, row=row_counter, value="mol_weight")
            ws.cell(column=col_counter + 3, row=row_counter, value="smiles")
        row_counter += 1


        try:
            all_data[compounds]
        except KeyError:
            formel = all_data_2[compounds]["formel"]
            mol_weight = all_data_2[compounds]["mol_weight"]
            smiles = all_data_2[compounds]["smiles"]
            comp_id = all_data_2[compounds]["comp_id"]
        else:
            formel = all_data[compounds]["formel"]
            mol_weight = all_data[compounds]["mol_weight"]
            smiles = all_data[compounds]["smiles"]
            comp_id = all_data[compounds]["comp_id"]

        ws.cell(column=col_counter, row=row_counter, value=compounds)
        ws.cell(column=col_counter + 1, row=row_counter, value=formel)
        ws.cell(column=col_counter + 2, row=row_counter, value=mol_weight)
        ws.cell(column=col_counter + 3, row=row_counter, value=smiles)
    wb_final.save(file)


    #
    # print(not_found)
    # # print(found_id)
    # # print(f"Found: {len(found_id)}, Not Found: {len(not_found)}")
    # print(non_duplicates)
    # print(len(non_duplicates))
    #
    temp_list = []
    for x in not_found:
        if x not in non_duplicates:
            temp_list.append(x)

    print(temp_list)
    print(len(temp_list))

    #
    #
    #
    # print(all_data.keys())

def testing(config, mps):
    _get_motherplate_layout(config, mps)



if __name__ == "__main__":
    config = configparser.ConfigParser()
    config.read("config.ini")
    mps = ["MP2022-004", "MP2022-001"]
    testing(config, mps)



