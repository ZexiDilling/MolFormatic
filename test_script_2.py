import csv
import os
import re
from pathlib import Path

from extra_functions import unit_converter
import numpy as np

from info import unit_converter_dict, unit_converter_list


def testing(well_dict, vol_needed_pure):

    groups = set()
    for well_counter, wells in enumerate(well_dict):
        if well_counter == 250:
            print(well_dict[wells])
        group = well_dict[wells]["group"]
        groups.add(group)
    print(groups)


def dose_response_worklist_writer(config, plate_layout, source_layout, dilution_layout, initial_plate,
                                  assay_name, fill_up=None):
    plate_amount = source_layout["plate_amount"]
    well_state_well_counter = 0
    if fill_up:
        dmso_well_counter = 0
        fill_up_vol = unit_converter(fill_up, old_unit_out=False, new_unit_out="n", as_list=True)
        dead_vol = 2.5

    output_folder = Path(f'{config["output_folders"]["output"]}/worklist')
    try:
        os.mkdir(output_folder)
    except OSError:
        print("directory exist")

    path = output_folder / assay_name
    try:
        os.mkdir(path)
    except OSError:
        print("directory exist")

    headlines = [headlines for headlines in config["worklist_headlines_v1"]]
    temp_file_name = f"Worklist_{assay_name}_{initial_plate}_dose_response"
    file = path / f"{temp_file_name}.csv"
    file_name_counter = 1
    while file.exists():
        file_name = f"{temp_file_name}_{file_name_counter}"
        file = path / f"{file_name}.csv"
        file_name_counter += 1
    file.touch()
    with open(file, "w", newline="\n") as csv_file:

        csv_writer = csv.writer(csv_file, delimiter=";")
        csv_writer.writerow(headlines)

        for plate in range(plate_amount):
            destination_plate = f"{assay_name}_{plate + initial_plate}"
            fill_up_dict = {}
            for wells in plate_layout:
                destination_well = plate_layout[wells]["well_id"]
                well_state = plate_layout[wells]["state"]

                # Check if the well is suppose to have samples in it, and if it does, add sample from MotherPlates

                if well_state == "sample":
                    current_sample = plate_layout[wells]["group"] + plate
                    current_concentration = plate_layout[wells]["concentration"]
                    current_stock = dilution_layout[current_concentration]["stock"]
                    transferee_vol = dilution_layout[current_concentration]["vol"]
                    source_well = source_layout["samples"][current_sample]["stocks"][current_stock]
                    source_plate = source_layout["samples"][current_sample]["plate"]

                    if fill_up:
                        dmso_trans_vol = fill_up_vol - transferee_vol
                        while dmso_trans_vol < (source_plate["dmso"][dmso_well_counter]["vol"] - dead_vol):
                            dmso_well_counter += 1

                        fill_up_dict[wells] = {"destination_plate": destination_plate,
                                               "destination_well": destination_well,
                                               "source_plate": source_plate["dmso"]["plate"],
                                               "soruce_well": source_plate["dmso"]["wells"][dmso_well_counter],
                                               "transferee_vol": dmso_trans_vol}

                elif well_state != "empty" or well_state != "blank":
                    transferee_vol = source_layout["states"][well_state]["vol"]

                    while transferee_vol < (source_plate["states"][well_state]["vol"] - dead_vol):
                        well_state_well_counter += 1

                    source_well = source_layout["states"][well_state]["well"][well_state_well_counter]
                    source_plate = source_layout["states"][well_state]["plate"]

                else:
                    pass
                csv_writer.writerow([destination_plate, destination_well, transferee_vol,
                                     source_well, source_plate])
            if fill_up:
                for counter in fill_up_dict:
                    destination_plate = fill_up_dict[counter]["destination_plate"]
                    destination_well = fill_up_dict[counter]["destination_well"]
                    transferee_vol = fill_up_dict[counter]["source_plate"]
                    source_well = fill_up_dict[counter]["soruce_well"]
                    source_plate = fill_up_dict[counter]["transferee_vol"]
                    csv_writer.writerow([destination_plate, destination_well, transferee_vol,
                                         source_well, source_plate])

def specific_round(value, base):
    new_base, _, _, unit = unit_converter(base, old_unit_out=False, new_unit_out=False, as_list=True)

    try:
        unit_converter_dict[unit]
    except KeyError:
        zeroes = 4
    else:

        temp_unit = str(unit_converter_dict[unit])
        if "e" in temp_unit:
            zeroes = int(temp_unit[-1])
        else:
            zeroes = int(temp_unit.count("0"))

        zeroes += 2
        new_base = float(new_base)

    return round(new_base * round(float(value)/new_base), zeroes)


def number_to_unit_converter(value, unit, rounding=False):
    new_value, new_unit, unit_type, _ = unit_converter(value, new_unit_out=unit, as_list=True)
    if rounding:
        temp_new_value = new_value
        while temp_new_value < 1:
            unit_count = unit_converter_list.index(new_unit)
            try:
                unit_converter_list[unit_count + 1]
            except IndexError:
                break
            else:
                temp_new_unit = unit_converter_list[unit_count + 1]
                temp_new_value = float(unit_converter(value, new_unit_out=temp_new_unit, as_list=True)[0])
                new_unit = temp_new_unit
        new_value = round(temp_new_value,2)

    return f"{new_value}{new_unit}{unit_type}"


def calculate_dilution_series(stock, max_concentration, min_concentration, dilutions_steps,
                              dilutions_factor, echo_min, final_vol, stock_dilution):

    unit_used = False

    stock_concentration_value = float(unit_converter(stock, new_unit_out=unit_used,
                                                     old_unit_out=False,  as_list=True)[0])
    max_concentration_value = float(unit_converter(max_concentration, new_unit_out=unit_used,
                                                   old_unit_out=False,  as_list=True)[0])
    min_concentration_value = float(unit_converter(min_concentration, new_unit_out=unit_used,
                                                   old_unit_out=False,  as_list=True)[0])
    echo_min_volume_value = float(unit_converter(echo_min, new_unit_out=unit_used,
                                                 old_unit_out=False,  as_list=True)[0])
    final_volume_value = float(unit_converter(final_vol, new_unit_out=unit_used,
                                              old_unit_out=False,  as_list=True)[0])

    all_concentration = []
    temp_concentration = max_concentration_value
    while temp_concentration > min_concentration_value:
        all_concentration.append(temp_concentration)
        temp_concentration = temp_concentration/dilutions_factor

    conc_check = []
    diff_check = []
    dmso_conc = []
    vol_needed_string = {}
    vol_needed_pure = {}
    all_stocks = [number_to_unit_converter(stock_concentration_value, "um")]
    temp_stock_concentration = stock_concentration_value
    # Add the minimum concentration at the end to include both endpoints
    for counter, count in enumerate(all_concentration):
        vol_needed_string[counter] = {"vol": "", "stock": "", "new_conc": ""}
        vol_needed_pure[counter] = {"vol": 0, "stock": 0, "new_conc": 0}
        conc_check.append(count)

        temp_vol = (count * final_volume_value) / temp_stock_concentration
        if temp_vol / echo_min_volume_value < 1:
            temp_stock_concentration = temp_stock_concentration/stock_dilution

            all_stocks.append(number_to_unit_converter(temp_stock_concentration, "uM"))
            temp_vol = (count * final_volume_value) / temp_stock_concentration
        new_temp_vol = specific_round(temp_vol, "2.5nL")

        temp_dmso_conc = round((new_temp_vol / final_volume_value) * 100, 4)

        while temp_dmso_conc > 1:
            new_temp_vol -= echo_min_volume_value
            temp_dmso_conc = round((new_temp_vol / final_volume_value) * 100, 4)

        new_temp_conc = (new_temp_vol * temp_stock_concentration) / final_volume_value
        dmso_conc.append(temp_dmso_conc)
        vol_needed_pure[counter] = {"vol": new_temp_vol,
                                    "stock": temp_stock_concentration,
                                    "new_conc": new_temp_conc}

        string_new_temp_vol = number_to_unit_converter(new_temp_vol, "nL", rounding=True)
        string_temp_stock_concentration = number_to_unit_converter(temp_stock_concentration, "mM", rounding=True)
        string_new_temp_conc = number_to_unit_converter(new_temp_conc, "uM", rounding=True)
        vol_needed_string[counter] = {"vol": string_new_temp_vol,
                                      "stock": string_temp_stock_concentration,
                                      "new_conc": string_new_temp_conc}

        if counter > 0:
            diff_1 = float(unit_converter(vol_needed_pure[counter]["new_conc"],
                                          old_unit_out=False, new_unit_out=unit_used, as_list=True)[0])
            diff_2 = float(unit_converter(vol_needed_pure[counter-1]["new_conc"],
                                          old_unit_out=False, new_unit_out=unit_used, as_list=True)[0])
            try:
                diff = round((diff_1 / diff_2) * 100, 1)
            except ZeroDivisionError:
                diff = "Error"
            diff_check.append(diff)
    #
    # print(f"conc_check - {conc_check}")
    # print(f"diff_check - {diff_check}")
    # print(f"vol_needed - {vol_needed_string}")
    # print(f"all_stocks - {all_stocks}")
    # print(f"dmso_conc - {dmso_conc}")

    return vol_needed_pure


if __name__ == "__main__":
    well_dict = {1: {'group': 0, 'well_id': 'A1', 'state': 'empty', 'colour': '#1e0bc8'}, 2: {'group': 0, 'well_id': 'B1', 'state': 'empty', 'colour': '#1e0bc8'}, 3: {'group': 0, 'well_id': 'C1', 'state': 'empty', 'colour': '#1e0bc8'}, 4: {'group': 0, 'well_id': 'D1', 'state': 'empty', 'colour': '#1e0bc8'}, 5: {'group': 0, 'well_id': 'E1', 'state': 'empty', 'colour': '#1e0bc8'}, 6: {'group': 0, 'well_id': 'F1', 'state': 'empty', 'colour': '#1e0bc8'}, 7: {'group': 0, 'well_id': 'G1', 'state': 'empty', 'colour': '#1e0bc8'}, 8: {'group': 0, 'well_id': 'H1', 'state': 'empty', 'colour': '#1e0bc8'}, 9: {'group': 0, 'well_id': 'I1', 'state': 'empty', 'colour': '#1e0bc8'}, 10: {'group': 0, 'well_id': 'J1', 'state': 'empty', 'colour': '#1e0bc8'}, 11: {'group': 0, 'well_id': 'K1', 'state': 'empty', 'colour': '#1e0bc8'}, 12: {'group': 0, 'well_id': 'L1', 'state': 'empty', 'colour': '#1e0bc8'}, 13: {'group': 0, 'well_id': 'M1', 'state': 'empty', 'colour': '#1e0bc8'}, 14: {'group': 0, 'well_id': 'N1', 'state': 'empty', 'colour': '#1e0bc8'}, 15: {'group': 0, 'well_id': 'O1', 'state': 'empty', 'colour': '#1e0bc8'}, 16: {'group': 0, 'well_id': 'P1', 'state': 'empty', 'colour': '#1e0bc8'}, 17: {'group': 0, 'well_id': 'A2', 'state': 'empty', 'colour': '#1e0bc8'}, 18: {'group': 0, 'well_id': 'B2', 'state': 'minimum', 'colour': '#ff8000'}, 19: {'group': 0, 'well_id': 'C2', 'state': 'minimum', 'colour': '#ff8000'}, 20: {'group': 0, 'well_id': 'D2', 'state': 'minimum', 'colour': '#ff8000'}, 21: {'group': 0, 'well_id': 'E2', 'state': 'minimum', 'colour': '#ff8000'}, 22: {'group': 0, 'well_id': 'F2', 'state': 'minimum', 'colour': '#ff8000'}, 23: {'group': 0, 'well_id': 'G2', 'state': 'minimum', 'colour': '#ff8000'}, 24: {'group': 0, 'well_id': 'H2', 'state': 'minimum', 'colour': '#ff8000'}, 25: {'group': 0, 'well_id': 'I2', 'state': 'minimum', 'colour': '#ff8000'}, 26: {'group': 0, 'well_id': 'J2', 'state': 'minimum', 'colour': '#ff8000'}, 27: {'group': 0, 'well_id': 'K2', 'state': 'minimum', 'colour': '#ff8000'}, 28: {'group': 0, 'well_id': 'L2', 'state': 'minimum', 'colour': '#ff8000'}, 29: {'group': 0, 'well_id': 'M2', 'state': 'minimum', 'colour': '#ff8000'}, 30: {'group': 0, 'well_id': 'N2', 'state': 'minimum', 'colour': '#ff8000'}, 31: {'group': 0, 'well_id': 'O2', 'state': 'minimum', 'colour': '#ff8000'}, 32: {'group': 0, 'well_id': 'P2', 'state': 'empty', 'colour': '#1e0bc8'}, 33: {'group': 0, 'well_id': 'A3', 'state': 'empty', 'colour': '#1e0bc8'}, 34: {'group': 0, 'well_id': 'B3', 'state': 'max', 'colour': '#790dc1'}, 35: {'group': 0, 'well_id': 'C3', 'state': 'max', 'colour': '#790dc1'}, 36: {'group': 0, 'well_id': 'D3', 'state': 'max', 'colour': '#790dc1'}, 37: {'group': 0, 'well_id': 'E3', 'state': 'max', 'colour': '#790dc1'}, 38: {'group': 0, 'well_id': 'F3', 'state': 'max', 'colour': '#790dc1'}, 39: {'group': 0, 'well_id': 'G3', 'state': 'max', 'colour': '#790dc1'}, 40: {'group': 0, 'well_id': 'H3', 'state': 'max', 'colour': '#790dc1'}, 41: {'group': 0, 'well_id': 'I3', 'state': 'max', 'colour': '#790dc1'}, 42: {'group': 0, 'well_id': 'J3', 'state': 'max', 'colour': '#790dc1'}, 43: {'group': 0, 'well_id': 'K3', 'state': 'max', 'colour': '#790dc1'}, 44: {'group': 0, 'well_id': 'L3', 'state': 'max', 'colour': '#790dc1'}, 45: {'group': 0, 'well_id': 'M3', 'state': 'max', 'colour': '#790dc1'}, 46: {'group': 0, 'well_id': 'N3', 'state': 'max', 'colour': '#790dc1'}, 47: {'group': 0, 'well_id': 'O3', 'state': 'max', 'colour': '#790dc1'}, 48: {'group': 0, 'well_id': 'P3', 'state': 'empty', 'colour': '#1e0bc8'}, 49: {'group': 0, 'well_id': 'A4', 'state': 'empty', 'colour': '#1e0bc8'}, 50: {'group': 1, 'well_id': 'B4', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 1, 'concentration': 1}, 51: {'group': 1, 'well_id': 'C4', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 2, 'concentration': 10}, 52: {'group': 2, 'well_id': 'D4', 'state': 'sample', 'colour': '#d54b28', 'replicate': 1, 'concentration': 8}, 53: {'group': 2, 'well_id': 'E4', 'state': 'sample', 'colour': '#d54b28', 'replicate': 3, 'concentration': 6}, 54: {'group': 3, 'well_id': 'F4', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 2, 'concentration': 4}, 55: {'group': 4, 'well_id': 'G4', 'state': 'sample', 'colour': '#a54b24', 'replicate': 1, 'concentration': 2}, 56: {'group': 4, 'well_id': 'H4', 'state': 'sample', 'colour': '#a54b24', 'replicate': 2, 'concentration': 11}, 57: {'group': 5, 'well_id': 'I4', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 1, 'concentration': 9}, 58: {'group': 5, 'well_id': 'J4', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 3, 'concentration': 7}, 59: {'group': 6, 'well_id': 'K4', 'state': 'sample', 'colour': '#754b1f', 'replicate': 2, 'concentration': 5}, 60: {'group': 7, 'well_id': 'L4', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 1, 'concentration': 3}, 61: {'group': 7, 'well_id': 'M4', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 3, 'concentration': 1}, 62: {'group': 8, 'well_id': 'N4', 'state': 'sample', 'colour': '#454b1b', 'replicate': 1, 'concentration': 10}, 63: {'group': 8, 'well_id': 'O4', 'state': 'sample', 'colour': '#454b1b', 'replicate': 3, 'concentration': 8}, 64: {'group': 0, 'well_id': 'P4', 'state': 'empty', 'colour': '#1e0bc8'}, 65: {'group': 0, 'well_id': 'A5', 'state': 'empty', 'colour': '#1e0bc8'}, 66: {'group': 1, 'well_id': 'B5', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 1, 'concentration': 2}, 67: {'group': 1, 'well_id': 'C5', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 2, 'concentration': 11}, 68: {'group': 2, 'well_id': 'D5', 'state': 'sample', 'colour': '#d54b28', 'replicate': 1, 'concentration': 9}, 69: {'group': 2, 'well_id': 'E5', 'state': 'sample', 'colour': '#d54b28', 'replicate': 3, 'concentration': 7}, 70: {'group': 3, 'well_id': 'F5', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 2, 'concentration': 5}, 71: {'group': 4, 'well_id': 'G5', 'state': 'sample', 'colour': '#a54b24', 'replicate': 1, 'concentration': 3}, 72: {'group': 4, 'well_id': 'H5', 'state': 'sample', 'colour': '#a54b24', 'replicate': 3, 'concentration': 1}, 73: {'group': 5, 'well_id': 'I5', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 1, 'concentration': 10}, 74: {'group': 5, 'well_id': 'J5', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 3, 'concentration': 8}, 75: {'group': 6, 'well_id': 'K5', 'state': 'sample', 'colour': '#754b1f', 'replicate': 2, 'concentration': 6}, 76: {'group': 7, 'well_id': 'L5', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 1, 'concentration': 4}, 77: {'group': 7, 'well_id': 'M5', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 3, 'concentration': 2}, 78: {'group': 8, 'well_id': 'N5', 'state': 'sample', 'colour': '#454b1b', 'replicate': 1, 'concentration': 11}, 79: {'group': 8, 'well_id': 'O5', 'state': 'sample', 'colour': '#454b1b', 'replicate': 3, 'concentration': 9}, 80: {'group': 0, 'well_id': 'P5', 'state': 'empty', 'colour': '#1e0bc8'}, 81: {'group': 0, 'well_id': 'A6', 'state': 'empty', 'colour': '#1e0bc8'}, 82: {'group': 1, 'well_id': 'B6', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 1, 'concentration': 3}, 83: {'group': 1, 'well_id': 'C6', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 3, 'concentration': 1}, 84: {'group': 2, 'well_id': 'D6', 'state': 'sample', 'colour': '#d54b28', 'replicate': 1, 'concentration': 10}, 85: {'group': 2, 'well_id': 'E6', 'state': 'sample', 'colour': '#d54b28', 'replicate': 3, 'concentration': 8}, 86: {'group': 3, 'well_id': 'F6', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 2, 'concentration': 6}, 87: {'group': 4, 'well_id': 'G6', 'state': 'sample', 'colour': '#a54b24', 'replicate': 1, 'concentration': 4}, 88: {'group': 4, 'well_id': 'H6', 'state': 'sample', 'colour': '#a54b24', 'replicate': 3, 'concentration': 2}, 89: {'group': 5, 'well_id': 'I6', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 1, 'concentration': 11}, 90: {'group': 5, 'well_id': 'J6', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 3, 'concentration': 9}, 91: {'group': 6, 'well_id': 'K6', 'state': 'sample', 'colour': '#754b1f', 'replicate': 2, 'concentration': 7}, 92: {'group': 7, 'well_id': 'L6', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 1, 'concentration': 5}, 93: {'group': 7, 'well_id': 'M6', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 3, 'concentration': 3}, 94: {'group': 8, 'well_id': 'N6', 'state': 'sample', 'colour': '#454b1b', 'replicate': 2, 'concentration': 1}, 95: {'group': 8, 'well_id': 'O6', 'state': 'sample', 'colour': '#454b1b', 'replicate': 3, 'concentration': 10}, 96: {'group': 0, 'well_id': 'P6', 'state': 'empty', 'colour': '#1e0bc8'}, 97: {'group': 0, 'well_id': 'A7', 'state': 'empty', 'colour': '#1e0bc8'}, 98: {'group': 1, 'well_id': 'B7', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 1, 'concentration': 4}, 99: {'group': 1, 'well_id': 'C7', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 3, 'concentration': 2}, 100: {'group': 2, 'well_id': 'D7', 'state': 'sample', 'colour': '#d54b28', 'replicate': 1, 'concentration': 11}, 101: {'group': 2, 'well_id': 'E7', 'state': 'sample', 'colour': '#d54b28', 'replicate': 3, 'concentration': 9}, 102: {'group': 3, 'well_id': 'F7', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 2, 'concentration': 7}, 103: {'group': 4, 'well_id': 'G7', 'state': 'sample', 'colour': '#a54b24', 'replicate': 1, 'concentration': 5}, 104: {'group': 4, 'well_id': 'H7', 'state': 'sample', 'colour': '#a54b24', 'replicate': 3, 'concentration': 3}, 105: {'group': 5, 'well_id': 'I7', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 2, 'concentration': 1}, 106: {'group': 5, 'well_id': 'J7', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 3, 'concentration': 10}, 107: {'group': 6, 'well_id': 'K7', 'state': 'sample', 'colour': '#754b1f', 'replicate': 2, 'concentration': 8}, 108: {'group': 7, 'well_id': 'L7', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 1, 'concentration': 6}, 109: {'group': 7, 'well_id': 'M7', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 3, 'concentration': 4}, 110: {'group': 8, 'well_id': 'N7', 'state': 'sample', 'colour': '#454b1b', 'replicate': 2, 'concentration': 2}, 111: {'group': 8, 'well_id': 'O7', 'state': 'sample', 'colour': '#454b1b', 'replicate': 3, 'concentration': 11}, 112: {'group': 0, 'well_id': 'P7', 'state': 'empty', 'colour': '#1e0bc8'}, 113: {'group': 0, 'well_id': 'A8', 'state': 'empty', 'colour': '#1e0bc8'}, 114: {'group': 1, 'well_id': 'B8', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 1, 'concentration': 5}, 115: {'group': 1, 'well_id': 'C8', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 3, 'concentration': 3}, 116: {'group': 2, 'well_id': 'D8', 'state': 'sample', 'colour': '#d54b28', 'replicate': 2, 'concentration': 1}, 117: {'group': 2, 'well_id': 'E8', 'state': 'sample', 'colour': '#d54b28', 'replicate': 3, 'concentration': 10}, 118: {'group': 3, 'well_id': 'F8', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 2, 'concentration': 8}, 119: {'group': 4, 'well_id': 'G8', 'state': 'sample', 'colour': '#a54b24', 'replicate': 1, 'concentration': 6}, 120: {'group': 4, 'well_id': 'H8', 'state': 'sample', 'colour': '#a54b24', 'replicate': 3, 'concentration': 4}, 121: {'group': 5, 'well_id': 'I8', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 2, 'concentration': 2}, 122: {'group': 5, 'well_id': 'J8', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 3, 'concentration': 11}, 123: {'group': 6, 'well_id': 'K8', 'state': 'sample', 'colour': '#754b1f', 'replicate': 2, 'concentration': 9}, 124: {'group': 7, 'well_id': 'L8', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 1, 'concentration': 7}, 125: {'group': 7, 'well_id': 'M8', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 3, 'concentration': 5}, 126: {'group': 8, 'well_id': 'N8', 'state': 'sample', 'colour': '#454b1b', 'replicate': 2, 'concentration': 3}, 127: {'group': 0, 'well_id': 'O8', 'state': 'sample', 'colour': '#ff00ff', 'replicate': 0, 'concentration': 0}, 128: {'group': 0, 'well_id': 'P8', 'state': 'empty', 'colour': '#1e0bc8'}, 129: {'group': 0, 'well_id': 'A9', 'state': 'empty', 'colour': '#1e0bc8'}, 130: {'group': 1, 'well_id': 'B9', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 1, 'concentration': 6}, 131: {'group': 1, 'well_id': 'C9', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 3, 'concentration': 4}, 132: {'group': 2, 'well_id': 'D9', 'state': 'sample', 'colour': '#d54b28', 'replicate': 2, 'concentration': 2}, 133: {'group': 2, 'well_id': 'E9', 'state': 'sample', 'colour': '#d54b28', 'replicate': 3, 'concentration': 11}, 134: {'group': 3, 'well_id': 'F9', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 2, 'concentration': 9}, 135: {'group': 4, 'well_id': 'G9', 'state': 'sample', 'colour': '#a54b24', 'replicate': 1, 'concentration': 7}, 136: {'group': 4, 'well_id': 'H9', 'state': 'sample', 'colour': '#a54b24', 'replicate': 3, 'concentration': 5}, 137: {'group': 5, 'well_id': 'I9', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 2, 'concentration': 3}, 138: {'group': 6, 'well_id': 'J9', 'state': 'sample', 'colour': '#754b1f', 'replicate': 1, 'concentration': 1}, 139: {'group': 6, 'well_id': 'K9', 'state': 'sample', 'colour': '#754b1f', 'replicate': 2, 'concentration': 10}, 140: {'group': 7, 'well_id': 'L9', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 1, 'concentration': 8}, 141: {'group': 7, 'well_id': 'M9', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 3, 'concentration': 6}, 142: {'group': 8, 'well_id': 'N9', 'state': 'sample', 'colour': '#454b1b', 'replicate': 2, 'concentration': 4}, 143: {'group': 0, 'well_id': 'O9', 'state': 'sample', 'colour': '#ff00ff', 'replicate': 0, 'concentration': 0}, 144: {'group': 0, 'well_id': 'P9', 'state': 'empty', 'colour': '#1e0bc8'}, 145: {'group': 0, 'well_id': 'A10', 'state': 'empty', 'colour': '#1e0bc8'}, 146: {'group': 1, 'well_id': 'B10', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 1, 'concentration': 7}, 147: {'group': 1, 'well_id': 'C10', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 3, 'concentration': 5}, 148: {'group': 2, 'well_id': 'D10', 'state': 'sample', 'colour': '#d54b28', 'replicate': 2, 'concentration': 3}, 149: {'group': 3, 'well_id': 'E10', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 1, 'concentration': 1}, 150: {'group': 3, 'well_id': 'F10', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 2, 'concentration': 10}, 151: {'group': 4, 'well_id': 'G10', 'state': 'sample', 'colour': '#a54b24', 'replicate': 1, 'concentration': 8}, 152: {'group': 4, 'well_id': 'H10', 'state': 'sample', 'colour': '#a54b24', 'replicate': 3, 'concentration': 6}, 153: {'group': 5, 'well_id': 'I10', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 2, 'concentration': 4}, 154: {'group': 6, 'well_id': 'J10', 'state': 'sample', 'colour': '#754b1f', 'replicate': 1, 'concentration': 2}, 155: {'group': 6, 'well_id': 'K10', 'state': 'sample', 'colour': '#754b1f', 'replicate': 2, 'concentration': 11}, 156: {'group': 7, 'well_id': 'L10', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 1, 'concentration': 9}, 157: {'group': 7, 'well_id': 'M10', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 3, 'concentration': 7}, 158: {'group': 8, 'well_id': 'N10', 'state': 'sample', 'colour': '#454b1b', 'replicate': 2, 'concentration': 5}, 159: {'group': 0, 'well_id': 'O10', 'state': 'sample', 'colour': '#ff00ff', 'replicate': 0, 'concentration': 0}, 160: {'group': 0, 'well_id': 'P10', 'state': 'empty', 'colour': '#1e0bc8'}, 161: {'group': 0, 'well_id': 'A11', 'state': 'empty', 'colour': '#1e0bc8'}, 162: {'group': 1, 'well_id': 'B11', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 1, 'concentration': 8}, 163: {'group': 1, 'well_id': 'C11', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 3, 'concentration': 6}, 164: {'group': 2, 'well_id': 'D11', 'state': 'sample', 'colour': '#d54b28', 'replicate': 2, 'concentration': 4}, 165: {'group': 3, 'well_id': 'E11', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 1, 'concentration': 2}, 166: {'group': 3, 'well_id': 'F11', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 2, 'concentration': 11}, 167: {'group': 4, 'well_id': 'G11', 'state': 'sample', 'colour': '#a54b24', 'replicate': 1, 'concentration': 9}, 168: {'group': 4, 'well_id': 'H11', 'state': 'sample', 'colour': '#a54b24', 'replicate': 3, 'concentration': 7}, 169: {'group': 5, 'well_id': 'I11', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 2, 'concentration': 5}, 170: {'group': 6, 'well_id': 'J11', 'state': 'sample', 'colour': '#754b1f', 'replicate': 1, 'concentration': 3}, 171: {'group': 6, 'well_id': 'K11', 'state': 'sample', 'colour': '#754b1f', 'replicate': 3, 'concentration': 1}, 172: {'group': 7, 'well_id': 'L11', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 1, 'concentration': 10}, 173: {'group': 7, 'well_id': 'M11', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 3, 'concentration': 8}, 174: {'group': 8, 'well_id': 'N11', 'state': 'sample', 'colour': '#454b1b', 'replicate': 2, 'concentration': 6}, 175: {'group': 0, 'well_id': 'O11', 'state': 'sample', 'colour': '#ff00ff', 'replicate': 0, 'concentration': 0}, 176: {'group': 0, 'well_id': 'P11', 'state': 'empty', 'colour': '#1e0bc8'}, 177: {'group': 0, 'well_id': 'A12', 'state': 'empty', 'colour': '#1e0bc8'}, 178: {'group': 1, 'well_id': 'B12', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 1, 'concentration': 9}, 179: {'group': 1, 'well_id': 'C12', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 3, 'concentration': 7}, 180: {'group': 2, 'well_id': 'D12', 'state': 'sample', 'colour': '#d54b28', 'replicate': 2, 'concentration': 5}, 181: {'group': 3, 'well_id': 'E12', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 1, 'concentration': 3}, 182: {'group': 3, 'well_id': 'F12', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 3, 'concentration': 1}, 183: {'group': 4, 'well_id': 'G12', 'state': 'sample', 'colour': '#a54b24', 'replicate': 1, 'concentration': 10}, 184: {'group': 4, 'well_id': 'H12', 'state': 'sample', 'colour': '#a54b24', 'replicate': 3, 'concentration': 8}, 185: {'group': 5, 'well_id': 'I12', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 2, 'concentration': 6}, 186: {'group': 6, 'well_id': 'J12', 'state': 'sample', 'colour': '#754b1f', 'replicate': 1, 'concentration': 4}, 187: {'group': 6, 'well_id': 'K12', 'state': 'sample', 'colour': '#754b1f', 'replicate': 3, 'concentration': 2}, 188: {'group': 7, 'well_id': 'L12', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 1, 'concentration': 11}, 189: {'group': 7, 'well_id': 'M12', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 3, 'concentration': 9}, 190: {'group': 8, 'well_id': 'N12', 'state': 'sample', 'colour': '#454b1b', 'replicate': 2, 'concentration': 7}, 191: {'group': 0, 'well_id': 'O12', 'state': 'sample', 'colour': '#ff00ff', 'replicate': 0, 'concentration': 0}, 192: {'group': 0, 'well_id': 'P12', 'state': 'empty', 'colour': '#1e0bc8'}, 193: {'group': 0, 'well_id': 'A13', 'state': 'empty', 'colour': '#1e0bc8'}, 194: {'group': 1, 'well_id': 'B13', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 1, 'concentration': 10}, 195: {'group': 1, 'well_id': 'C13', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 3, 'concentration': 8}, 196: {'group': 2, 'well_id': 'D13', 'state': 'sample', 'colour': '#d54b28', 'replicate': 2, 'concentration': 6}, 197: {'group': 3, 'well_id': 'E13', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 1, 'concentration': 4}, 198: {'group': 3, 'well_id': 'F13', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 3, 'concentration': 2}, 199: {'group': 4, 'well_id': 'G13', 'state': 'sample', 'colour': '#a54b24', 'replicate': 1, 'concentration': 11}, 200: {'group': 4, 'well_id': 'H13', 'state': 'sample', 'colour': '#a54b24', 'replicate': 3, 'concentration': 9}, 201: {'group': 5, 'well_id': 'I13', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 2, 'concentration': 7}, 202: {'group': 6, 'well_id': 'J13', 'state': 'sample', 'colour': '#754b1f', 'replicate': 1, 'concentration': 5}, 203: {'group': 6, 'well_id': 'K13', 'state': 'sample', 'colour': '#754b1f', 'replicate': 3, 'concentration': 3}, 204: {'group': 7, 'well_id': 'L13', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 2, 'concentration': 1}, 205: {'group': 7, 'well_id': 'M13', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 3, 'concentration': 10}, 206: {'group': 8, 'well_id': 'N13', 'state': 'sample', 'colour': '#454b1b', 'replicate': 2, 'concentration': 8}, 207: {'group': 0, 'well_id': 'O13', 'state': 'sample', 'colour': '#ff00ff', 'replicate': 0, 'concentration': 0}, 208: {'group': 0, 'well_id': 'P13', 'state': 'empty', 'colour': '#1e0bc8'}, 209: {'group': 0, 'well_id': 'A14', 'state': 'empty', 'colour': '#1e0bc8'}, 210: {'group': 1, 'well_id': 'B14', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 1, 'concentration': 11}, 211: {'group': 1, 'well_id': 'C14', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 3, 'concentration': 9}, 212: {'group': 2, 'well_id': 'D14', 'state': 'sample', 'colour': '#d54b28', 'replicate': 2, 'concentration': 7}, 213: {'group': 3, 'well_id': 'E14', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 1, 'concentration': 5}, 214: {'group': 3, 'well_id': 'F14', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 3, 'concentration': 3}, 215: {'group': 4, 'well_id': 'G14', 'state': 'sample', 'colour': '#a54b24', 'replicate': 2, 'concentration': 1}, 216: {'group': 4, 'well_id': 'H14', 'state': 'sample', 'colour': '#a54b24', 'replicate': 3, 'concentration': 10}, 217: {'group': 5, 'well_id': 'I14', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 2, 'concentration': 8}, 218: {'group': 6, 'well_id': 'J14', 'state': 'sample', 'colour': '#754b1f', 'replicate': 1, 'concentration': 6}, 219: {'group': 6, 'well_id': 'K14', 'state': 'sample', 'colour': '#754b1f', 'replicate': 3, 'concentration': 4}, 220: {'group': 7, 'well_id': 'L14', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 2, 'concentration': 2}, 221: {'group': 7, 'well_id': 'M14', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 3, 'concentration': 11}, 222: {'group': 8, 'well_id': 'N14', 'state': 'sample', 'colour': '#454b1b', 'replicate': 2, 'concentration': 9}, 223: {'group': 0, 'well_id': 'O14', 'state': 'sample', 'colour': '#ff00ff', 'replicate': 0, 'concentration': 0}, 224: {'group': 0, 'well_id': 'P14', 'state': 'empty', 'colour': '#1e0bc8'}, 225: {'group': 0, 'well_id': 'A15', 'state': 'empty', 'colour': '#1e0bc8'}, 226: {'group': 1, 'well_id': 'B15', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 2, 'concentration': 1}, 227: {'group': 1, 'well_id': 'C15', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 3, 'concentration': 10}, 228: {'group': 2, 'well_id': 'D15', 'state': 'sample', 'colour': '#d54b28', 'replicate': 2, 'concentration': 8}, 229: {'group': 3, 'well_id': 'E15', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 1, 'concentration': 6}, 230: {'group': 3, 'well_id': 'F15', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 3, 'concentration': 4}, 231: {'group': 4, 'well_id': 'G15', 'state': 'sample', 'colour': '#a54b24', 'replicate': 2, 'concentration': 2}, 232: {'group': 4, 'well_id': 'H15', 'state': 'sample', 'colour': '#a54b24', 'replicate': 3, 'concentration': 11}, 233: {'group': 5, 'well_id': 'I15', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 2, 'concentration': 9}, 234: {'group': 6, 'well_id': 'J15', 'state': 'sample', 'colour': '#754b1f', 'replicate': 1, 'concentration': 7}, 235: {'group': 6, 'well_id': 'K15', 'state': 'sample', 'colour': '#754b1f', 'replicate': 3, 'concentration': 5}, 236: {'group': 7, 'well_id': 'L15', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 2, 'concentration': 3}, 237: {'group': 8, 'well_id': 'M15', 'state': 'sample', 'colour': '#454b1b', 'replicate': 1, 'concentration': 1}, 238: {'group': 8, 'well_id': 'N15', 'state': 'sample', 'colour': '#454b1b', 'replicate': 2, 'concentration': 10}, 239: {'group': 0, 'well_id': 'O15', 'state': 'sample', 'colour': '#ff00ff', 'replicate': 0, 'concentration': 0}, 240: {'group': 0, 'well_id': 'P15', 'state': 'empty', 'colour': '#1e0bc8'}, 241: {'group': 0, 'well_id': 'A16', 'state': 'empty', 'colour': '#1e0bc8'}, 242: {'group': 1, 'well_id': 'B16', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 2, 'concentration': 2}, 243: {'group': 1, 'well_id': 'C16', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 3, 'concentration': 11}, 244: {'group': 2, 'well_id': 'D16', 'state': 'sample', 'colour': '#d54b28', 'replicate': 2, 'concentration': 9}, 245: {'group': 3, 'well_id': 'E16', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 1, 'concentration': 7}, 246: {'group': 3, 'well_id': 'F16', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 3, 'concentration': 5}, 247: {'group': 4, 'well_id': 'G16', 'state': 'sample', 'colour': '#a54b24', 'replicate': 2, 'concentration': 3}, 248: {'group': 5, 'well_id': 'H16', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 1, 'concentration': 1}, 249: {'group': 5, 'well_id': 'I16', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 2, 'concentration': 10}, 250: {'group': 6, 'well_id': 'J16', 'state': 'sample', 'colour': '#754b1f', 'replicate': 1, 'concentration': 8}, 251: {'group': 6, 'well_id': 'K16', 'state': 'sample', 'colour': '#754b1f', 'replicate': 3, 'concentration': 6}, 252: {'group': 7, 'well_id': 'L16', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 2, 'concentration': 4}, 253: {'group': 8, 'well_id': 'M16', 'state': 'sample', 'colour': '#454b1b', 'replicate': 1, 'concentration': 2}, 254: {'group': 8, 'well_id': 'N16', 'state': 'sample', 'colour': '#454b1b', 'replicate': 2, 'concentration': 11}, 255: {'group': 0, 'well_id': 'O16', 'state': 'sample', 'colour': '#ff00ff', 'replicate': 0, 'concentration': 0}, 256: {'group': 0, 'well_id': 'P16', 'state': 'empty', 'colour': '#1e0bc8'}, 257: {'group': 0, 'well_id': 'A17', 'state': 'empty', 'colour': '#1e0bc8'}, 258: {'group': 1, 'well_id': 'B17', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 2, 'concentration': 3}, 259: {'group': 2, 'well_id': 'C17', 'state': 'sample', 'colour': '#d54b28', 'replicate': 1, 'concentration': 1}, 260: {'group': 2, 'well_id': 'D17', 'state': 'sample', 'colour': '#d54b28', 'replicate': 2, 'concentration': 10}, 261: {'group': 3, 'well_id': 'E17', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 1, 'concentration': 8}, 262: {'group': 3, 'well_id': 'F17', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 3, 'concentration': 6}, 263: {'group': 4, 'well_id': 'G17', 'state': 'sample', 'colour': '#a54b24', 'replicate': 2, 'concentration': 4}, 264: {'group': 5, 'well_id': 'H17', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 1, 'concentration': 2}, 265: {'group': 5, 'well_id': 'I17', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 2, 'concentration': 11}, 266: {'group': 6, 'well_id': 'J17', 'state': 'sample', 'colour': '#754b1f', 'replicate': 1, 'concentration': 9}, 267: {'group': 6, 'well_id': 'K17', 'state': 'sample', 'colour': '#754b1f', 'replicate': 3, 'concentration': 7}, 268: {'group': 7, 'well_id': 'L17', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 2, 'concentration': 5}, 269: {'group': 8, 'well_id': 'M17', 'state': 'sample', 'colour': '#454b1b', 'replicate': 1, 'concentration': 3}, 270: {'group': 8, 'well_id': 'N17', 'state': 'sample', 'colour': '#454b1b', 'replicate': 3, 'concentration': 1}, 271: {'group': 0, 'well_id': 'O17', 'state': 'sample', 'colour': '#ff00ff', 'replicate': 0, 'concentration': 0}, 272: {'group': 0, 'well_id': 'P17', 'state': 'empty', 'colour': '#1e0bc8'}, 273: {'group': 0, 'well_id': 'A18', 'state': 'empty', 'colour': '#1e0bc8'}, 274: {'group': 1, 'well_id': 'B18', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 2, 'concentration': 4}, 275: {'group': 2, 'well_id': 'C18', 'state': 'sample', 'colour': '#d54b28', 'replicate': 1, 'concentration': 2}, 276: {'group': 2, 'well_id': 'D18', 'state': 'sample', 'colour': '#d54b28', 'replicate': 2, 'concentration': 11}, 277: {'group': 3, 'well_id': 'E18', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 1, 'concentration': 9}, 278: {'group': 3, 'well_id': 'F18', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 3, 'concentration': 7}, 279: {'group': 4, 'well_id': 'G18', 'state': 'sample', 'colour': '#a54b24', 'replicate': 2, 'concentration': 5}, 280: {'group': 5, 'well_id': 'H18', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 1, 'concentration': 3}, 281: {'group': 5, 'well_id': 'I18', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 3, 'concentration': 1}, 282: {'group': 6, 'well_id': 'J18', 'state': 'sample', 'colour': '#754b1f', 'replicate': 1, 'concentration': 10}, 283: {'group': 6, 'well_id': 'K18', 'state': 'sample', 'colour': '#754b1f', 'replicate': 3, 'concentration': 8}, 284: {'group': 7, 'well_id': 'L18', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 2, 'concentration': 6}, 285: {'group': 8, 'well_id': 'M18', 'state': 'sample', 'colour': '#454b1b', 'replicate': 1, 'concentration': 4}, 286: {'group': 8, 'well_id': 'N18', 'state': 'sample', 'colour': '#454b1b', 'replicate': 3, 'concentration': 2}, 287: {'group': 0, 'well_id': 'O18', 'state': 'sample', 'colour': '#ff00ff', 'replicate': 0, 'concentration': 0}, 288: {'group': 0, 'well_id': 'P18', 'state': 'empty', 'colour': '#1e0bc8'}, 289: {'group': 0, 'well_id': 'A19', 'state': 'empty', 'colour': '#1e0bc8'}, 290: {'group': 1, 'well_id': 'B19', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 2, 'concentration': 5}, 291: {'group': 2, 'well_id': 'C19', 'state': 'sample', 'colour': '#d54b28', 'replicate': 1, 'concentration': 3}, 292: {'group': 2, 'well_id': 'D19', 'state': 'sample', 'colour': '#d54b28', 'replicate': 3, 'concentration': 1}, 293: {'group': 3, 'well_id': 'E19', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 1, 'concentration': 10}, 294: {'group': 3, 'well_id': 'F19', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 3, 'concentration': 8}, 295: {'group': 4, 'well_id': 'G19', 'state': 'sample', 'colour': '#a54b24', 'replicate': 2, 'concentration': 6}, 296: {'group': 5, 'well_id': 'H19', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 1, 'concentration': 4}, 297: {'group': 5, 'well_id': 'I19', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 3, 'concentration': 2}, 298: {'group': 6, 'well_id': 'J19', 'state': 'sample', 'colour': '#754b1f', 'replicate': 1, 'concentration': 11}, 299: {'group': 6, 'well_id': 'K19', 'state': 'sample', 'colour': '#754b1f', 'replicate': 3, 'concentration': 9}, 300: {'group': 7, 'well_id': 'L19', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 2, 'concentration': 7}, 301: {'group': 8, 'well_id': 'M19', 'state': 'sample', 'colour': '#454b1b', 'replicate': 1, 'concentration': 5}, 302: {'group': 8, 'well_id': 'N19', 'state': 'sample', 'colour': '#454b1b', 'replicate': 3, 'concentration': 3}, 303: {'group': 0, 'well_id': 'O19', 'state': 'sample', 'colour': '#ff00ff', 'replicate': 0, 'concentration': 0}, 304: {'group': 0, 'well_id': 'P19', 'state': 'empty', 'colour': '#1e0bc8'}, 305: {'group': 0, 'well_id': 'A20', 'state': 'empty', 'colour': '#1e0bc8'}, 306: {'group': 1, 'well_id': 'B20', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 2, 'concentration': 6}, 307: {'group': 2, 'well_id': 'C20', 'state': 'sample', 'colour': '#d54b28', 'replicate': 1, 'concentration': 4}, 308: {'group': 2, 'well_id': 'D20', 'state': 'sample', 'colour': '#d54b28', 'replicate': 3, 'concentration': 2}, 309: {'group': 3, 'well_id': 'E20', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 1, 'concentration': 11}, 310: {'group': 3, 'well_id': 'F20', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 3, 'concentration': 9}, 311: {'group': 4, 'well_id': 'G20', 'state': 'sample', 'colour': '#a54b24', 'replicate': 2, 'concentration': 7}, 312: {'group': 5, 'well_id': 'H20', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 1, 'concentration': 5}, 313: {'group': 5, 'well_id': 'I20', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 3, 'concentration': 3}, 314: {'group': 6, 'well_id': 'J20', 'state': 'sample', 'colour': '#754b1f', 'replicate': 2, 'concentration': 1}, 315: {'group': 6, 'well_id': 'K20', 'state': 'sample', 'colour': '#754b1f', 'replicate': 3, 'concentration': 10}, 316: {'group': 7, 'well_id': 'L20', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 2, 'concentration': 8}, 317: {'group': 8, 'well_id': 'M20', 'state': 'sample', 'colour': '#454b1b', 'replicate': 1, 'concentration': 6}, 318: {'group': 8, 'well_id': 'N20', 'state': 'sample', 'colour': '#454b1b', 'replicate': 3, 'concentration': 4}, 319: {'group': 0, 'well_id': 'O20', 'state': 'sample', 'colour': '#ff00ff', 'replicate': 0, 'concentration': 0}, 320: {'group': 0, 'well_id': 'P20', 'state': 'empty', 'colour': '#1e0bc8'}, 321: {'group': 0, 'well_id': 'A21', 'state': 'empty', 'colour': '#1e0bc8'}, 322: {'group': 1, 'well_id': 'B21', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 2, 'concentration': 7}, 323: {'group': 2, 'well_id': 'C21', 'state': 'sample', 'colour': '#d54b28', 'replicate': 1, 'concentration': 5}, 324: {'group': 2, 'well_id': 'D21', 'state': 'sample', 'colour': '#d54b28', 'replicate': 3, 'concentration': 3}, 325: {'group': 3, 'well_id': 'E21', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 2, 'concentration': 1}, 326: {'group': 3, 'well_id': 'F21', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 3, 'concentration': 10}, 327: {'group': 4, 'well_id': 'G21', 'state': 'sample', 'colour': '#a54b24', 'replicate': 2, 'concentration': 8}, 328: {'group': 5, 'well_id': 'H21', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 1, 'concentration': 6}, 329: {'group': 5, 'well_id': 'I21', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 3, 'concentration': 4}, 330: {'group': 6, 'well_id': 'J21', 'state': 'sample', 'colour': '#754b1f', 'replicate': 2, 'concentration': 2}, 331: {'group': 6, 'well_id': 'K21', 'state': 'sample', 'colour': '#754b1f', 'replicate': 3, 'concentration': 11}, 332: {'group': 7, 'well_id': 'L21', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 2, 'concentration': 9}, 333: {'group': 8, 'well_id': 'M21', 'state': 'sample', 'colour': '#454b1b', 'replicate': 1, 'concentration': 7}, 334: {'group': 8, 'well_id': 'N21', 'state': 'sample', 'colour': '#454b1b', 'replicate': 3, 'concentration': 5}, 335: {'group': 0, 'well_id': 'O21', 'state': 'sample', 'colour': '#ff00ff', 'replicate': 0, 'concentration': 0}, 336: {'group': 0, 'well_id': 'P21', 'state': 'empty', 'colour': '#1e0bc8'}, 337: {'group': 0, 'well_id': 'A22', 'state': 'empty', 'colour': '#1e0bc8'}, 338: {'group': 1, 'well_id': 'B22', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 2, 'concentration': 8}, 339: {'group': 2, 'well_id': 'C22', 'state': 'sample', 'colour': '#d54b28', 'replicate': 1, 'concentration': 6}, 340: {'group': 2, 'well_id': 'D22', 'state': 'sample', 'colour': '#d54b28', 'replicate': 3, 'concentration': 4}, 341: {'group': 3, 'well_id': 'E22', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 2, 'concentration': 2}, 342: {'group': 3, 'well_id': 'F22', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 3, 'concentration': 11}, 343: {'group': 4, 'well_id': 'G22', 'state': 'sample', 'colour': '#a54b24', 'replicate': 2, 'concentration': 9}, 344: {'group': 5, 'well_id': 'H22', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 1, 'concentration': 7}, 345: {'group': 5, 'well_id': 'I22', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 3, 'concentration': 5}, 346: {'group': 6, 'well_id': 'J22', 'state': 'sample', 'colour': '#754b1f', 'replicate': 2, 'concentration': 3}, 347: {'group': 7, 'well_id': 'K22', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 1, 'concentration': 1}, 348: {'group': 7, 'well_id': 'L22', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 2, 'concentration': 10}, 349: {'group': 8, 'well_id': 'M22', 'state': 'sample', 'colour': '#454b1b', 'replicate': 1, 'concentration': 8}, 350: {'group': 8, 'well_id': 'N22', 'state': 'sample', 'colour': '#454b1b', 'replicate': 3, 'concentration': 6}, 351: {'group': 0, 'well_id': 'O22', 'state': 'sample', 'colour': '#ff00ff', 'replicate': 0, 'concentration': 0}, 352: {'group': 0, 'well_id': 'P22', 'state': 'empty', 'colour': '#1e0bc8'}, 353: {'group': 0, 'well_id': 'A23', 'state': 'empty', 'colour': '#1e0bc8'}, 354: {'group': 1, 'well_id': 'B23', 'state': 'sample', 'colour': '#11B4D4', 'replicate': 2, 'concentration': 9}, 355: {'group': 2, 'well_id': 'C23', 'state': 'sample', 'colour': '#d54b28', 'replicate': 1, 'concentration': 7}, 356: {'group': 2, 'well_id': 'D23', 'state': 'sample', 'colour': '#d54b28', 'replicate': 3, 'concentration': 5}, 357: {'group': 3, 'well_id': 'E23', 'state': 'sample', 'colour': '#42B4D9', 'replicate': 2, 'concentration': 3}, 358: {'group': 4, 'well_id': 'F23', 'state': 'sample', 'colour': '#a54b24', 'replicate': 1, 'concentration': 1}, 359: {'group': 4, 'well_id': 'G23', 'state': 'sample', 'colour': '#a54b24', 'replicate': 2, 'concentration': 10}, 360: {'group': 5, 'well_id': 'H23', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 1, 'concentration': 8}, 361: {'group': 5, 'well_id': 'I23', 'state': 'sample', 'colour': '#72B4DE', 'replicate': 3, 'concentration': 6}, 362: {'group': 6, 'well_id': 'J23', 'state': 'sample', 'colour': '#754b1f', 'replicate': 2, 'concentration': 4}, 363: {'group': 7, 'well_id': 'K23', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 1, 'concentration': 2}, 364: {'group': 7, 'well_id': 'L23', 'state': 'sample', 'colour': '#A2B4E2', 'replicate': 2, 'concentration': 11}, 365: {'group': 8, 'well_id': 'M23', 'state': 'sample', 'colour': '#454b1b', 'replicate': 1, 'concentration': 9}, 366: {'group': 8, 'well_id': 'N23', 'state': 'sample', 'colour': '#454b1b', 'replicate': 3, 'concentration': 7}, 367: {'group': 0, 'well_id': 'O23', 'state': 'sample', 'colour': '#ff00ff', 'replicate': 0, 'concentration': 0}, 368: {'group': 0, 'well_id': 'P23', 'state': 'empty', 'colour': '#1e0bc8'}, 369: {'group': 0, 'well_id': 'A24', 'state': 'empty', 'colour': '#1e0bc8'}, 370: {'group': 0, 'well_id': 'B24', 'state': 'empty', 'colour': '#1e0bc8'}, 371: {'group': 0, 'well_id': 'C24', 'state': 'empty', 'colour': '#1e0bc8'}, 372: {'group': 0, 'well_id': 'D24', 'state': 'empty', 'colour': '#1e0bc8'}, 373: {'group': 0, 'well_id': 'E24', 'state': 'empty', 'colour': '#1e0bc8'}, 374: {'group': 0, 'well_id': 'F24', 'state': 'empty', 'colour': '#1e0bc8'}, 375: {'group': 0, 'well_id': 'G24', 'state': 'empty', 'colour': '#1e0bc8'}, 376: {'group': 0, 'well_id': 'H24', 'state': 'empty', 'colour': '#1e0bc8'}, 377: {'group': 0, 'well_id': 'I24', 'state': 'empty', 'colour': '#1e0bc8'}, 378: {'group': 0, 'well_id': 'J24', 'state': 'empty', 'colour': '#1e0bc8'}, 379: {'group': 0, 'well_id': 'K24', 'state': 'empty', 'colour': '#1e0bc8'}, 380: {'group': 0, 'well_id': 'L24', 'state': 'empty', 'colour': '#1e0bc8'}, 381: {'group': 0, 'well_id': 'M24', 'state': 'empty', 'colour': '#1e0bc8'}, 382: {'group': 0, 'well_id': 'N24', 'state': 'empty', 'colour': '#1e0bc8'}, 383: {'group': 0, 'well_id': 'O24', 'state': 'empty', 'colour': '#1e0bc8'}, 384: {'group': 0, 'well_id': 'P24', 'state': 'empty', 'colour': '#1e0bc8'}}
    stock = "10mM"
    max_concentration = "90uM"
    min_concentration = "0.0015uM"
    echo_min = "2.5nL"
    final_vol = "12uL"
    # Test the function
    dilution_steps = 11
    dilution_factor = 3
    stock_dilution = 100
    vol_needed_pure = calculate_dilution_series(stock, max_concentration, min_concentration, dilution_steps,
                                                dilution_factor, echo_min, final_vol, stock_dilution)
    #ToDo Make a popup to choose source_layout for dose_reponse_worklist_writer and make it possible to import from excel & CSV instead

    source_layout = {"plate_amount": int,
                     "states": {state: {"plate": str,
                                        "wells": [wells]}},
                     "samples": {int: {"stocks": {int: well},
                                       "plate": str}}
                     "dmso": {"plate": str,
                              "wells": [wells]}
                     }

    dose_response_worklist_writer(config, well_dict, source_layout, vol_needed_pure, initial_plate,
                                  assay_name, fill_up=None)
    testing(well_dict, vol_needed_pure)

