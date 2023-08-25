from pathlib import Path
import xml.etree.ElementTree as ET
import PySimpleGUI as sg

from gui_functions import bio_compound_info_from_worklist


def getting_echo_data(folder, bio_sample_dict):


    # random data
    skipped_reason = {"empty_wells": {},
                      "over_under_working_vol": {},
                      "precipitated": {},
                      "other": {},
                      "skipped":
                          {"empty_wells": 0,
                           "over_under_working_vol": 0,
                           "precipitated": 0,
                           "other": 0,
                           "total_skipped": 0}
                      }
    skipped_wells = {}
    working_list = {}
    all_data = []
    skip_well_counter = 0

    # counting Transferees
    trans_plate_counter = []
    counter_plates = []

    # counting transferees all data
    all_trans_counter = {}

    # Data for completed plates:
    zero_error_trans_plate = []
    error_trans_plate = []

    all_files = Path(folder).glob("**/*")
    file_list = [file for file in all_files if file.is_file()]
    file_amount = len(file_list)
    for file_index, files in enumerate(file_list):
        if files.name.startswith("Transfer"):
            # path = self.file_names(self.main_folder)
            doc = ET.parse(files)
            root = doc.getroot()

            # for counting plates and transferees
            for plates in root.iter("plate"):
                barcode = plates.get("barcode")
                source_destination = plates.get("type")

                if source_destination == "destination":
                    counter_plates.append(barcode)
                    temp_d_barcode = barcode
                if source_destination == "source":
                    temp_s_barcode = barcode

            try:
                all_trans_counter[temp_d_barcode].append(temp_s_barcode)
            except KeyError:
                all_trans_counter[temp_d_barcode] = [temp_s_barcode]

            if temp_d_barcode in bio_sample_dict and temp_s_barcode.casefold() != "source_2" and temp_s_barcode != "Sourch_2":
            # if temp_s_barcode.casefold() != "source_2" and temp_s_barcode != "Sourch_2":

                # find amount of well that is skipped
                for wells in root.iter("skippedwells"):
                    wells_skipped = wells.get("total")
                    if int(wells_skipped) != 0:
                        all_data.append(wells_skipped)
                        skip_well_counter += int(wells_skipped)

                        # finds barcode for source and destination
                        for dates in root.iter("transfer"):
                            date = dates.get("date")
                            all_data.append(date)

                        # finds barcode for source and destination
                        for plates in root.iter("plate"):
                            barcode = plates.get("barcode")
                            source_destination = plates.get("type")
                            all_data.append(source_destination + ", " + barcode)

                            if source_destination == "source":
                                temp_barcode = barcode
                                try:
                                    skipped_wells[barcode]
                                except KeyError:
                                    skipped_wells[barcode] = {}
                            if source_destination == "destination":
                                temp_dest_barcode = barcode
                                error_trans_plate.append(temp_dest_barcode)
                                try:
                                    working_list[barcode]
                                except KeyError:
                                    working_list[barcode] = {}

                                try:
                                    working_list[barcode][temp_barcode]
                                except KeyError:
                                    working_list[barcode][temp_barcode] = []

                        # finds destination and source wells data
                        for z in range(int(wells_skipped)):
                            temp_trans = []
                            destination_well = wells[z].get("dn")
                            source_well = wells[z].get("n")
                            reason = wells[z].get("reason")
                            vt = wells[z].get("vt")
                            all_data.append("SW: " + source_well + " DW: " + destination_well + " vol: " + vt)
                            all_data.append(" reason: " + reason)
                            temp_trans.append(source_well)
                            temp_trans.append(float(vt))
                            temp_trans.append(destination_well)
                            # Gets only the error code from reason
                            reason = reason.split(":")[0]
                            try:
                                skipped_wells[temp_barcode][source_well]["counter"] += 1
                                skipped_wells[temp_barcode][source_well]["vol"] += float(vt)
                                skipped_wells[temp_barcode][source_well]["reason"].append(reason)
                            except KeyError:
                                skipped_wells[temp_barcode][source_well] = {"counter": 1, "vol": float(vt), "reason": [reason]}

                            working_list[temp_dest_barcode][temp_barcode].append(temp_trans)

                            if reason == "MM0202007":
                                skipped_reason["skipped"]["empty_wells"] += 1
                                try:
                                    skipped_reason["empty_wells"][temp_s_barcode]
                                except KeyError:
                                    skipped_reason["empty_wells"][temp_s_barcode] = [source_well]
                                else:
                                    skipped_reason["empty_wells"][temp_s_barcode].append(source_well)
                            elif reason == "MM0203001":
                                skipped_reason["skipped"]["precipitated"] += 1
                                try:
                                    skipped_reason["precipitated"][temp_s_barcode]
                                except KeyError:
                                    skipped_reason["precipitated"][temp_s_barcode] = [source_well]
                                else:
                                    skipped_reason["precipitated"][temp_s_barcode].append(source_well)
                            elif reason == "MM0202006":
                                skipped_reason["skipped"]["over_under_working_vol"] += 1
                                try:
                                    skipped_reason["over_under_working_vol"][temp_s_barcode]
                                except KeyError:
                                    skipped_reason["over_under_working_vol"][temp_s_barcode] = [source_well]
                                else:
                                    skipped_reason["over_under_working_vol"][temp_s_barcode].append(source_well)
                            else:
                                skipped_reason["skipped"]["other"] += 1
                                try:
                                    skipped_reason["other"][temp_s_barcode]
                                except KeyError:
                                    skipped_reason["other"][temp_s_barcode] = [source_well]
                                else:
                                    skipped_reason["other"][temp_s_barcode].append(source_well)
                            skipped_reason["skipped"]["total_skipped"] += 1

                    else:
                        # finds barcode for destination
                        for plates in root.iter("plate"):
                            barcode = plates.get("barcode")
                            source_destination = plates.get("type")

                            if source_destination == "destination":
                                zero_error_trans_plate.append(barcode)
        print(f"{file_index}/{file_amount} done")
        # counting plates
        # counts the number of repeated barcodes and makes a list with the barcode and amount of
        # instance with the same name

    for plates in counter_plates:
        number = counter_plates.count(plates)
        trans_plate_counter.append(str(plates) + "," + str(number))

    # remove duplicates from the list
    trans_plate_counter = list(dict.fromkeys(trans_plate_counter))
    # print(skip_well_counter)
    plates_to_remove = []
    if zero_error_trans_plate:
        for plates in zero_error_trans_plate:
            if plates in error_trans_plate:
                plates_to_remove.append(plates)

        for plates in plates_to_remove:
            zero_error_trans_plate.remove(plates)

    return all_data, skipped_wells, skip_well_counter, working_list, trans_plate_counter, all_trans_counter, \
           zero_error_trans_plate, temp_d_barcode, skipped_reason


def get_plate_information(config, sg, folder_plates):
    all_files = Path(folder_plates).glob("**/*")
    file_list = [file for file in all_files if file.is_file()]
    return bio_compound_info_from_worklist(config, sg, file_list)


if __name__ == "__main__":
    import configparser
    config = configparser.ConfigParser()
    config.read("config.ini")
    folder_echo = r"C:\Users\phch\Desktop\test\echo_trams_data"
    folder_plates = r"C:\Users\phch\OneDrive - Danmarks Tekniske Universitet\Mapper\Python_data\alpha_SO\worklisted_used_with_good_data"
    bio_sample_dict = get_plate_information(config, sg, folder_plates)

    all_data, skipped_wells, skip_well_counter, working_list, trans_plate_counter, all_trans_counter, \
    zero_error_trans_plate, temp_d_barcode, skipped_reason = getting_echo_data(folder_echo, bio_sample_dict)

    print(skipped_reason)

    skipped_wells = {'empty_wells':
                         {'MP2022-081': ['P16'],
                          'MP2022-095': ['P2', 'P8', 'P10'],
                          'MP2022-111': ['I8'],
                          'MP2022-125': ['H17'],
                          'MP2022-033': ['J2']},
                     'over_under_working_vol':
                         {'MP2022-003': ['P19'],
                          'MP2022-030': ['C20'],
                          'MP2022-048': ['B1', 'D1', 'I2', 'N1', 'P11', 'P24'],
                          'MP2022-058': ['D1', 'I2', 'P1'],
                          'MP2022-021': ['O5', 'O15', 'O4', 'N7', 'N9', 'L13', 'L17', 'L18', 'N22'],
                          'MP2022-038': ['P1', 'N8', 'P10'], 'MP2022-069': ['P3', 'A5'],
                          'MP2022-068': ['M23'],
                          'MP2022-118': ['H12'],
                          'MP2022-132': ['M1', 'M11', 'O11', 'M2', 'O6', 'M8', 'O20', 'I22', 'J1', 'L1', 'N7', 'H11',
                                         'P17', 'P21', 'N23', 'J2', 'L2', 'L12', 'P12', 'P16', 'P20', 'J24', 'N24',
                                         'P24']},
                     'precipitated':
                         {'MP2022-096': ['M9'],
                          'MP2022-007': ['A21'],
                          'MP2022-020': ['A20'],
                          'MP2022-001': ['H11'],
                          'MP2022-013': ['L21'],
                          'MP2022-048': ['P4']},
                     'other': {},
                     'skipped': {'empty_wells': 7,
                                 'over_under_working_vol': 51,
                                 'precipitated': 6,
                                 'other': 0,
                                 'total_skipped': 64}}