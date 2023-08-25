import csv
from datetime import date
from math import ceil
import os
from pathlib import Path

from extra_functions import unit_converter
from info import plate_96_row, plate_96_column, compound_well_layout
from plate_formatting import mother_plate_generator as mpg
from gui_popup import new_headlines_popup


class CSVWriter:
    def __str__(self):
        """
        Writes CSV files
        Not sure this needs to be a class...

        :return: CSV files, in different formate.
        """

    @staticmethod
    def compound_freezer_writer(path, compound_list):
        """
        Writes the tube file for the comPOUND freezer

        :param path: main output folder
        :type path: str
        :param compound_list: list of compounds
        :type compound_list: list
        :return: CSV for the comPOUND freezer, to fetch tubes
        """
        try:
            os.mkdir(f"{path}/comPOUND")
        except OSError:
            print("directory exist")

        path = f"{path}/comPOUND"

        file_name = f"comPOUND_'{date.today()}'.txt"
        complete_file_name = os.path.join(path, file_name)
        with open(complete_file_name, "w+", newline="\n") as f:
            csv_writer = csv.writer(f)
            for compound in compound_list:
                csv_writer.writerow([compound])

    @staticmethod
    def mp_workflow_handler(path, compound_list, plate_layout, tube_rack_list, samples_per_plate=96):
        """
        Generate a dict of tubes and their placement in a plate

        :param path: Main output folder
        :type path: str
        :param compound_list: List of compounds
        :type compound_list: list
        :param plate_layout: Layout for the plate to get well information. is pulled from INFO!
        :type plate_layout: dict
        :param tube_rack_list: Barcode for the tube-racks
        :type tube_rack_list: list
        :param samples_per_plate: How many samples there are per plate. Should always be 96.
        :type samples_per_plate: int
        :return: A dict of tubes with what rack they are in, and what well/spot in that rack.
            and A CSV file for PlateButler to run the MotherPlate protocol.
        :rtype: dict
        """
        tube_dict = {}
        try:
            os.mkdir(f"{path}/mp_files_{date.today()}")
        except OSError:
            print("directory exist")

        path = f"{path}/mp_files_{date.today()}"

        counter_compounds = 0
        if not tube_rack_list:
            plate_amount = len(compound_list)/samples_per_plate
            plate_amount = ceil(plate_amount)
            tube_rack_list = []
            for i in range(plate_amount):
                tube_rack_list.append(i+1)
        for plate in tube_rack_list:
            tube_dict[plate] = {}
            file_name = f"{path}/{plate}.txt"
            with open(file_name, "w", newline="\n") as f:
                csv_writer = csv.writer(f, delimiter=";")
                csv_writer.writerow(["RowCol", "tubeBarcode"])
                for index, _ in enumerate(compound_list):
                    if counter_compounds <= len(compound_list):
                        if index < samples_per_plate:
                            tube_dict[plate][plate_layout[index]] = compound_list[counter_compounds]
                            csv_writer.writerow([plate_layout[index], compound_list[counter_compounds]])
                            counter_compounds += 1

        return tube_dict

    @staticmethod
    def mp_to_pb_writer(path, mp_plates):
        """
        Writes the output file from the PlateButler to double check if things looks fine.
        CAN maybe BE DELEDET!

        :param path: Main Output folder
        :type path: str
        :param mp_plates: Dict with information about what tubes goes into witch well. from 4 x 96 to 384.
        :type mp_plates: dict
        :return: CSV file, to compare with PlateButler output file
        """

        try:
            os.mkdir(f"{path}/pb_output")
        except OSError:
            print("directory exist")

        path = f"{path}/pb_output"
        file_name = f"{path}/pb_mp_generated_output_{date.today()}.csv"
        with open(file_name, "w", newline="\n") as f:
            csv_writer = csv.writer(f, delimiter=";")
            for mp_barcode in mp_plates:
                for values in mp_plates[mp_barcode]:
                    csv_writer.writerow([mp_barcode, *values])

    def compound_freezer_handler(self, path, compound_list, mp_name, tube_rack_list=None, plate_layout=None):
        """
        Generate CSV files for the comPOUND freezer
        And PlateButler for producing MotherPlates ... Maybe... Put a Gate in, and make it optional.

        :param path: Main Output folder
        :type path: str
        :param compound_list: List of compounds
        :type compound_list: list
        :param mp_name: Main name for MotherPlates.
        :type mp_name: str
        :param tube_rack_list: List of barcodes for racks
        :type tube_rack_list: list
        :param plate_layout: Layout for the plate to get well information. is pulled from INFO!
        :type plate_layout: dict
        :return: 3 CSV files. One for the comPOUND freezer to fetch tubes, one for PlateButler to run the protocol,
            1 for comparing output files from PlateButler with Theoretical output
        """
        if not plate_layout:
            plate_layout = plate_96_row
        self.compound_freezer_writer(path, compound_list)

        # tube_dict = self.mp_workflow_handler(path, compound_list, plate_layout, tube_rack_list)
        # _, pb_mp_output = mpg(tube_dict, mp_name)
        # self.mp_to_pb_writer(path, pb_mp_output)

    @staticmethod
    def dp_writer(config, dp_dict, path):
        """
        Writes CSV file for PlateButler protocol for producing DaughterPlates

        :param dp_dict: Dict over compounds. Source and Destination info and volume to transferee
        :type dp_dict: dict
        :param path: Main output folder
        :type path: str
        :return: CSV file for PlateButler to run Assay production Protocol.
        """
        try:
            os.mkdir(f"{path}/dp_output")
        except OSError:
            print("directory exist")

        path = f"{path}/dp_output"
        file_name = f"{path}/dp_{date.today()}.csv"

        with open(file_name, "w", newline="\n") as csv_file:
            csv_writer = csv.writer(csv_file, delimiter="\t")
            headlines = [headlines for headlines in config["worklist_headlines"]]
            csv_writer.writerow(headlines)

            for destination_plate in dp_dict:

                for destination_well, vol, source_well, source_sample, plate in dp_dict[destination_plate]:
                    csv_writer.writerow([destination_plate, destination_well, vol, source_well, plate, source_sample])

    @staticmethod
    def compound_freezer_to_2d_csv_simulate(compound_tube_dict, path):
        """
        Generate a series of CSV files, that should be the same as the one created by the 2D plate scanner

        :param compound_tube_dict: A dict over all the tubes
        :type compound_tube_dict: dict
        :param path: The path for the output file
        :type path: str
        :return: A series of excel files. one file per 96 compounds
        """
        for rack in compound_tube_dict:

            file_name = f"{path}/2D_barcodes_sim_{date.today()}_{rack}.csv"

            with open(file_name, "w", newline="\n") as csv_file:
                csv_writer = csv.writer(csv_file, delimiter=";")
                csv_writer.writerow(["RowCol", "tubeBarcode"])
                for line in compound_tube_dict[rack]:
                    csv_writer.writerow([line, compound_tube_dict[rack][line]])

    @staticmethod
    def plate_list_writer(plate_list, output_folder):
        path = f"{output_folder}/dp_output"
        file_name = f"{path}/Plate_list_{date.today()}.csv"
        with open(file_name, "w", newline="\n") as csv_file:
            csv_writer = csv.writer(csv_file)

            csv_writer.writerow(plate_list)

    @staticmethod
    def source_plate_dilution(sample_dict, output_folder):
        try:
            os.mkdir(f"{output_folder}/plate_dilution_output")
        except OSError:
            print("directory exist")

        path = f"{output_folder}/plate_dilution_output"
        file_name = f"{path}/source_plate_dilution_{date.today()}.csv"

        with open(file_name, "w", newline="\n") as csv_file:
            csv_writer = csv.writer(csv_file, delimiter=";")
            csv_writer.writerow(
                ["concentration", "destination_plates", "destination_well", "compound", "Volume", "source_well",
                 "source_plates", "sour_plate_type"])

            for samples in sample_dict:
                for dilution in sample_dict[samples]:
                    for counter in sample_dict[samples][dilution]:
                        conc = dilution
                        destination_barcode = sample_dict[samples][dilution][counter]["destination_well"]
                        destination_well = sample_dict[samples][dilution][counter]["destination_plate"]
                        compound_id = samples
                        volume = sample_dict[samples][dilution][counter]["sample_vol"]
                        source_well = sample_dict[samples][dilution][counter]["source_well"]
                        source_barcode = sample_dict[samples][dilution][counter]["source_plate"]
                        sour_plate_type = sample_dict[samples][dilution][counter]["plate_type"]
                        csv_writer.writerow([conc, destination_barcode, destination_well, compound_id, volume,
                                             source_well, source_barcode, sour_plate_type])

    @staticmethod
    def plate_dilution(sample_dict, output_folder):

        path = f"{output_folder}/plate_dilution_output"
        try:
            os.mkdir(path)
        except OSError:
            print("directory exist")

        file_name = f"{path}/plate_dilution_{date.today()}.csv"

        with open(file_name, "w", newline="\n") as csv_file:
            csv_writer = csv.writer(csv_file, delimiter="\t")
            csv_writer.writerow(["concentration", "destination_plates", "destination_well", "compound", "Volume",
                                 "source_well", "source_plates", "sour_plate_type"])

            for plate_replica in sample_dict:
                if plate_replica != "other_state":
                    for samples in sample_dict[plate_replica]:
                        for sample_replica in sample_dict[plate_replica][samples]:
                            for dilution in sample_dict[plate_replica][samples][sample_replica]:
                                for wells in sample_dict[plate_replica][samples][sample_replica][dilution]["destination_well"]:
                                    conc = dilution
                                    destination_barcode = sample_dict[plate_replica][samples][sample_replica][dilution][
                                        "destination_plate"]
                                    destination_well = wells
                                    compound_id = samples
                                    volume = sample_dict[plate_replica][samples][sample_replica][dilution]["sample_vol"]
                                    source_well = sample_dict[plate_replica][samples][sample_replica][dilution]["source_well"]
                                    source_barcode = sample_dict[plate_replica][samples][sample_replica][dilution][
                                        "source_plate"]
                                    sour_plate_type = sample_dict[plate_replica][samples][sample_replica][dilution][
                                        "plate_type"]

                                    csv_writer.writerow([conc, destination_barcode, destination_well, compound_id, volume,
                                                         source_well, source_barcode, sour_plate_type])
                else:
                    for states in sample_dict["other_state"]:
                        for plates in sample_dict["other_state"][states]:
                            for wells in sample_dict["other_state"][states][plates]:
                                conc = sample_dict["other_state"][states][plates][wells]["concentration"]
                                destination_barcode = sample_dict["other_state"][states][plates][wells]["destination_barcode"]
                                destination_well = sample_dict["other_state"][states][plates][wells]["destination_well"]
                                compound_id = sample_dict["other_state"][states][plates][wells]["compound_id"]
                                volume = sample_dict["other_state"][states][plates][wells]["volume"]
                                source_well = sample_dict["other_state"][states][plates][wells]["source_well"]
                                source_barcode = sample_dict["other_state"][states][plates][wells]["source_barcode"]
                                sour_plate_type = sample_dict["other_state"][states][plates][wells]["sour_plate_type"]

                                csv_writer.writerow([conc, destination_barcode, destination_well, compound_id, volume,
                                                     source_well, source_barcode, sour_plate_type])

    @staticmethod
    def dose_response_worklist_writer(config, plate_layout, source_layout, plate_amount, dilution_layout, initial_plate,
                                      assay_name, volume, sample_per_plate, fill_up=None):
        """
        Writes the worklist for a dose-response set-up.
        The worklist is for an Echo-liquid handler.
        :param config: The config handler, with all the default information in the config file.
        :type config: configparser.ConfigParser
        :param plate_layout: The layout for a single plate, where each sample goes and the samples dilutions
        :type plate_layout: dict
        :param source_layout: The layout for the source_plates. Where the stock solution of the samples are, and where
            all the controls, and filler-solvent is located
        :type source_layout: dict
        :param plate_amount: The amount of plates to be produced
        :type plate_amount: int
        :param dilution_layout: The guide to tell witch stock and how much should go in each well
        :type dilution_layout: dict
        :param initial_plate: The first plate - used for naming
        :type initial_plate: int
        :param assay_name: The name of the assay / experiment
        :type assay_name: str
        :param volume: How much volume to add for the controls, in string formate including unit type
        :type volume: str
        :param fill_up: If all the wells should have the same amount of liquid in them, can be None.
            It is the amount of liquid to fill up to, in string formate including unit type.
        :type fill_up: str
        :return:
        """

        missing_vol = {"filler": 0,
                       }

        if fill_up:
            fill_up_vol, _, _, _ = unit_converter(fill_up, old_unit_out=False, new_unit_out=False, as_list=True)
            print(fill_up_vol)

        echo_min_transfer = "2.5nL"
        echo_min_transfer, _, _, _ = unit_converter(echo_min_transfer, old_unit_out=False, new_unit_out=False,
                                                    as_list=True)
        volume, _, _, _ = unit_converter(volume, old_unit_out=False, new_unit_out=False, as_list=True)
        dead_vol = "2.5uL"
        dead_vol, _, _, _ = unit_converter(dead_vol, old_unit_out=False, new_unit_out=False, as_list=True)
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
        headlines.append("compound_id")
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
                        current_sample = plate_layout[wells]["group"] + (sample_per_plate * plate)
                        print(current_sample)
                        current_concentration = plate_layout[wells]["concentration"] - 1
                        try:
                            current_stock = dilution_layout[current_concentration]["stock"]
                        except KeyError:
                            print(f"key-error for dilution layout - {plate_layout[wells]}")
                            break
                        transferee_vol = dilution_layout[current_concentration]["vol"]
                        try:
                            source_well = source_layout[well_state][current_sample]["conc"][current_stock]["well"]
                        except KeyError:
                            continue
                        source_plate = source_layout[well_state][current_sample]["conc"][current_stock]["plate"]
                        temp_compound_id = source_layout[well_state][current_sample]["compound_id"]

                        temp_trans_vol, _, _, _ = unit_converter(transferee_vol, old_unit_out=False, new_unit_out="n",
                                                                 as_list=True)
                        temp_fill_vol, _, _, _ = unit_converter(fill_up_vol, old_unit_out=False, new_unit_out="n",
                                                                as_list=True)

                        test_val = temp_fill_vol - temp_trans_vol

                        if fill_up and test_val > 1:

                            dmso_trans_vol = fill_up_vol - transferee_vol
                            well_counter = source_layout["filler"]["well_counter"]
                            filler_wells = len(source_layout["filler"]["source_wells"])

                            if well_counter >= filler_wells:
                                missing_vol["filler"] += dmso_trans_vol
                            else:
                                while dmso_trans_vol > (
                                        source_layout["filler"]["source_wells"][well_counter]["vol"] - dead_vol):
                                    well_counter += 1
                                    source_layout["filler"]["well_counter"] += 1
                                    if well_counter >= filler_wells:
                                        missing_vol["filler"] += dmso_trans_vol
                                        break

                            try:
                                source_layout["filler"]["source_wells"][well_counter]["well_id"]
                            except KeyError:
                                temp_dmso_source_well = "Missing"
                                temp_dmso_source_plate = "Missing"
                            else:
                                temp_dmso_source_well = source_layout["filler"]["source_wells"][well_counter]["well_id"]
                                temp_dmso_source_plate = source_layout["filler"]["source_wells"][well_counter]["plate"]
                                source_layout["filler"]["source_wells"][well_counter]["vol"] -= dmso_trans_vol

                            fill_up_dict[wells] = {"destination_plate": destination_plate,
                                                   "destination_well": destination_well,
                                                   "source_plate": temp_dmso_source_plate,
                                                   "soruce_well": temp_dmso_source_well,
                                                   "transferee_vol": dmso_trans_vol}

                    elif source_layout[well_state]["use"]:

                        transferee_vol = volume
                        well_counter = source_layout[well_state]["well_counter"]
                        temp_amount_well_state = len(source_layout[well_state]["source_wells"])
                        if well_counter >= temp_amount_well_state:
                            missing_vol[well_state] += transferee_vol
                        else:
                            while transferee_vol > (
                                    source_layout[well_state]["source_wells"][well_counter]["vol"] - dead_vol):
                                well_counter += 1
                                if well_counter >= temp_amount_well_state:
                                    try:
                                        missing_vol[well_state]
                                    except KeyError:
                                        missing_vol[well_state] = transferee_vol
                                    else:
                                        missing_vol[well_state] += transferee_vol

                        source_layout[well_state]["well_counter"] = well_counter

                        try:
                            source_layout[well_state]["source_wells"][well_counter]["well_id"]
                        except KeyError:
                            source_well = "Missing_Vol"
                            source_plate = "Missing_Vol"
                        else:
                            source_layout[well_state]["source_wells"][well_counter]["vol"] -= transferee_vol
                            source_well = source_layout[well_state]["source_wells"][well_counter]["well_id"]
                            source_plate = source_layout[well_state]["source_wells"][well_counter]["plate"]
                        temp_compound_id = well_state
                    else:
                        continue

                    transferee_vol, _, _, _ = unit_converter(transferee_vol, old_unit_out=False, new_unit_out="n",
                                                             as_list=True)
                    csv_writer.writerow([destination_plate, destination_well, transferee_vol,
                                         source_well, source_plate, temp_compound_id])
                    destination_well = transferee_vol = source_well = source_plate = None
                if fill_up:
                    for counter in fill_up_dict:
                        destination_plate = fill_up_dict[counter]["destination_plate"]
                        destination_well = fill_up_dict[counter]["destination_well"]
                        transferee_vol = fill_up_dict[counter]["transferee_vol"]
                        transferee_vol, _, _, _ = unit_converter(transferee_vol, old_unit_out=False, new_unit_out="n",
                                                                 as_list=True)
                        source_well = fill_up_dict[counter]["soruce_well"]
                        source_plate = fill_up_dict[counter]["source_plate"]
                        temp_compound_id = "dmso-filler"
                        csv_writer.writerow([destination_plate, destination_well, transferee_vol,
                                             source_well, source_plate, temp_compound_id])
        for states in missing_vol:
            missing_vol[states], _, _, _ = unit_converter(missing_vol[states], old_unit_out=False, new_unit_out="u",
                                                                 as_list=True)
        print(f"missing_vol - {missing_vol}uL")
        print(file)

    @staticmethod
    def worklist_writer(config, plate_layout, mps, free_well_dict, assay_name, plate_amount, initial_plate, volume,
                        sample_direction, worklist_analyse_method, control_bonus_source, control_samples,
                        bonus_compound):
        """
        Writes a worklist to use with plateButler
        :param config: The config handler, with all the default information in the config file.
        :type config: configparser.ConfigParser
        :param plate_layout: The layout for the plate with values for each well, what state they are in
        :type plate_layout: dict
        :param mps: A list of MotherPlates
        :type: mps: list
        :param free_well_dict: A dict over the wells with compounds on motherplates that can be used
        :type free_well_dict: dict
        :param assay_name: the name of the assay, and the name used for the destination plate in the workinglist
        :type assay_name: str
        :param plate_amount: Number of plates to produce
        :type plate_amount: int
        :param initial_plate: the starting number of the plates
        :type initial_plate: int
        :param volume: How much volume to transfere to each well. The same amount of liquid will be transfered to each.
        :type volume: int
        :param sample_direction: The direction for the sample layout. Only relevant if the analyse_method differet from
            "Single Point"
        :type sample_direction: str
        :param worklist_analyse_method: The method use for the sample layout.
        :type worklist_analyse_method: str
        :param control_bonus_source: The layout for where controls and bonus compounds are placed and the name for the
        souceplate
        :type control_bonus_source: dict or None
        :param control_samples: A dict for samples for positive and negative control
        :type: control_samples: None or dict
        :param bonus_compound: If there are a compound that needs to be added to multiple well_states
        :type bonus_compound: None or dict
        :return:
        """

        if control_bonus_source:
            source_plate_bonus = list(control_bonus_source.keys())[0]

        mp_plate_counter = 0
        mp_well_counter = 0

        output_folder = Path(f'{config["output_folders"]["output"]}/worklist')
        try:
            os.mkdir(output_folder)
        except OSError:
            print("directory exist")

        path = output_folder/assay_name
        try:
            os.mkdir(path)
        except OSError:
            print("directory exist")

        headlines = [headlines for headlines in config["worklist_headlines_v1"]]
        temp_file_name = f"Worklist_{assay_name}_{initial_plate}_to_{plate_amount + initial_plate - 1}"
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
                for wells in plate_layout["well_layout"]:
                    destination_well = plate_layout["well_layout"][wells]["well_id"]
                    well_state = plate_layout["well_layout"][wells]["state"]

                    # Check if the well is suppose to have samples in it, and if it does, add sample from MotherPlates
                    if well_state == "sample":
                        # Check if there are any wells left, that have not been used, else it skips to the next plate
                        source_plate = mps[mp_plate_counter]
                        while len(free_well_dict[source_plate]) == 0:
                            print(source_plate)
                            mp_plate_counter += 1
                            mp_well_counter = 0
                            if mp_plate_counter == len(mps):
                                msg = "Not Enough MotherPlates"
                                return file, msg

                            source_plate = mps[mp_plate_counter]

                        source_well = free_well_dict[source_plate][mp_well_counter]

                        # Writes the data to a CSV file
                        csv_writer.writerow([destination_plate, destination_well, volume,
                                             source_well, source_plate])
                        mp_well_counter += 1
                        if mp_well_counter >= len(free_well_dict[source_plate]):
                            mp_plate_counter += 1
                            mp_well_counter = 0
                            if mp_plate_counter == len(mps):
                                msg = "Not Enough MotherPlates"
                                return file, msg

                    else:
                        if control_samples[well_state]["use"] or bonus_compound[well_state]:
                            try:
                                control_compound = control_bonus_source[well_state]["compound"]
                            except KeyError:
                                control_compound = control_bonus_source["bonus"]["compound"]
                                well_state = "bonus"

                            temp_plate_type = control_bonus_source[well_state]["plate_type"]
                            temp_barcode = control_bonus_source[well_state]["barcode"]
                            dead_vol = \
                                config["plate_types_values"][temp_plate_type].split(",")[0]
                            dead_vol = float(dead_vol)
                            got_compound = True
                            for well in control_bonus_source[well_state]["well_vol"]:
                                if control_bonus_source[well_state]["well_vol"][well] >= dead_vol + volume:
                                    source_well = well
                                    source_plate = temp_barcode
                                    control_bonus_source[well_state]["well_vol"][well] -= volume
                                    got_compound = True
                                    break
                                else:
                                    source_well = f"Missing volume for {well_state}"
                                    source_plate = temp_barcode
                                    got_compound = False
                            if not got_compound:
                                print(f"Missing compound for {well_state}")
                            csv_writer.writerow([destination_plate, destination_well, volume,
                                                 source_well, source_plate])

        return file, None


class CSVReader:
    def __str__(self):
        """
        Reads CSV files

        :return: Dict's with information.
        """
    # @staticmethod
    # def _pb_mp_output_files(file):
    #     """
    #     Reads the output file from PlateButler. When running MotherPlates Production
    #
    #     :param file: PlateButler output file
    #     :type file: str
    #     :return: 2 dict. 1 for the data and 1 for the plates
    #     :rtype: dict, dict
    #     """
    #     headline = ["DestinationBarcode", "DestinationWell", "compound_id", "Volume", "Date"]
    #     dict_data = {}
    #     destination_plates = {}
    #     with open(file) as f:
    #         for row_index, line in enumerate(f):
    #             values = line.split(";")
    #             dict_data[f"Transferee_{row_index}"] = {}
    #             for clm_index, value in enumerate(values):
    #                 dict_data[f"Transferee_{row_index}"][headline[clm_index]] = value.strip("\n")
    #                 if headline[clm_index] == "DestinationBarcode":
    #                     destination_plates[value] = {}
    #                     destination_plates[value]["DestinationBarcode"] = value
    #                     destination_plates[value]["Date"] = date.today()
    #             dict_data[f"Transferee_{row_index}"]["Date"] = date.today()
    #             dict_data[f"Transferee_{row_index}"]["Row_Counter"] = row_index
    #
    #
    #     return dict_data, destination_plates

    @staticmethod
    def _pb_mp_output_files(file):
        """
        Reads the output file from PlateButler. When running MotherPlates Production

        :param file: PlateButler output file
        :type file: str
        :return: 2 dict. 1 for the data and 1 for the plates
        :rtype: dict, dict
        """

        headline = ["DestinationBarcode", "DestinationWell", "compound_id", "Volume", "Date"]
        dict_data = {}
        destination_plates = {}
        with open(file) as f:
            for row_index, line in enumerate(f):
                values = line.split(";")
                dict_data[f"Transferee_{row_index}"] = {}
                for clm_index, value in enumerate(values):
                    dict_data[f"Transferee_{row_index}"][headline[clm_index]] = value.strip("\n")
                    if headline[clm_index] == "DestinationBarcode":
                        destination_plates[value] = {}
                        destination_plates[value]["DestinationBarcode"] = value
                        destination_plates[value]["Date"] = date.today()
                dict_data[f"Transferee_{row_index}"]["Date"] = date.today()
                dict_data[f"Transferee_{row_index}"]["Row_Counter"] = row_index

        return dict_data, destination_plates

    @staticmethod
    def pb_tube_files_ind(file, tube_dict):
        """
        Reads Tube files from 2D-scanner

        :param file: Tube files
        :type file: str
        :param tube_dict: Dict for the tubes
        :type tube_dict: dict
        :return: The tube dict, updated with information
        """
        plate = file.split("/")
        plate = plate[-1]
        plate = plate.removesuffix(".txt")
        tube_dict[plate] = {}
        with open(file) as f:
            for line in f:
                values = line.split(";")
                if values != ['RowCol', 'tubeBarcode\n']:
                    tube_dict[plate][values[0]] = values[1].strip("\n")

    @staticmethod
    def _tab_file_reader(file):
        """
        Reads CSV file with tab format that PlateButler prefer to read.

        :param file: CSV file that needs to be read
        :type file: str
        :return: 2 dict. 1 for data, 1 for plates
        :rtype: dict, dict
        """
        counter = -1
        headline = []
        dict_data = {}
        destination_plates = {}
        with open(file) as f:

            for line in f:
                counter += 1
                values = line.split("\t")
                if counter > 0:
                    dict_data[f"Transferee_{counter}"] = {}
                for index, value in enumerate(values):
                    if value != "\n":
                        if counter == 0:

                            headline.append(value.strip("\n"))

                        else:
                            dict_data[f"Transferee_{counter}"][headline[index]] = value.strip("\n")
                            dict_data[f"Transferee_{counter}"]["Date"] = date.today()
                            dict_data[f"Transferee_{counter}"]["Row_Counter"] = counter

                            if headline[index] == "DestinationBarcode":
                                destination_plates[value] = {}
                                destination_plates[value]["DestinationBarcode"] = value.strip("\n")
                                destination_plates[value]["date"] = date.today()
                                #destination_plates[clm_info]["location"] = "Freezer-1"

        return dict_data, destination_plates

    def csv_r_controller(self, csv_file, file_type):
        """
        Handles the CSV reader to get files where they needs to go.
        Could be deleted.¨

        :param csv_file: The CSV file that needs  to be read
        :type csv_file: str
        :param file_type: What kind of CSV file it is
        :type file_type: str
        :return: 2 dict. 1 for data, 1 for plates
        :rtype: dict, dict
        """
        if file_type == "tab":
            dict_data, destination_plates = self._tab_file_reader(csv_file)
        elif file_type == "pb_mp_output":
            dict_data, destination_plates = self._pb_mp_output_files(csv_file)

        return dict_data, destination_plates

    @staticmethod
    def tube_list_to_dict(csv_file):
        """
        Convert a list of tubes from a file, to a dict

        :param csv_file: A CSV file, containing all the tubes
        :type csv_file: str
        :return: a dict with all the tube id's
        :rtype: dict
        """
        temp_dict = {}
        tube_dict = {}
        counter = 1
        well_counter = 0
        temp_dict[counter] = {}
        tube_dict[counter] = {}
        with open(csv_file) as f:
            for index, line in enumerate(f):
                temp_dict[counter][compound_well_layout[well_counter]] = line.strip("\n")
                if [compound_well_layout[well_counter]] == ["A12"]:
                    key_order = plate_96_column
                    tube_dict[counter] = {value: temp_dict[counter][value] for value in key_order}

                    counter += 1
                    temp_dict[counter] = {}
                    well_counter = -1
                well_counter += 1

        return tube_dict

    @staticmethod
    def compound_plates(csv_file):
        temp_dict = {}
        counter = 0
        headlines = []
        with open(csv_file) as file:
            for index, line in enumerate(file):
                if index == 1:
                    temp_headlines = line.split(";")

                    for headline_index, headline in enumerate(temp_headlines):
                        headline = headline.casefold().replace(" ", "_")
                        if "weight" in headline:
                            headline = "mass"
                        if headline == "id":
                            id_number = headline_index
                        headlines.append(headline)

                if index > 1:
                    counter += 1
                    data = line.split(";")
                    sample = data[id_number]
                    temp_dict[sample] = {}
                    for data_index, data in enumerate(data):
                        temp_dict[sample][headlines[data_index]] = data

        return temp_dict

    @staticmethod
    def grab_headlines(csv_file):
        splitter = [";", ",", "."]
        split_indicator = 0
        with open(csv_file) as file:
            for index, lines in enumerate(file):
                lines = lines.removesuffix("\n")
                if index == 0:
                    # Check if the file is a CSV file

                    headlines = lines.split(splitter[split_indicator])
                    return headlines

    @staticmethod
    def echo_worklist_to_dict(config, config_headline, csv_file, right_headlines, new_headline, sample_dict,
                              all_destination_plates):
        splitter = [";", ",", "."]
        split_indicator = 0
        with open(csv_file) as file:
            for index, lines in enumerate(file):
                lines = lines.removesuffix("\n")
                if index == 0:

                    # Check if the file is a CSV file
                    headlines = lines.split(splitter[split_indicator])
                    while len(headlines) < 2:
                        split_indicator += 1
                        headlines = lines.split(splitter[split_indicator])
                        if split_indicator > len(splitter):
                            return "Not a CSV file", headlines, sample_dict

                    # Check if the headlines for the CSV file is correct, if not sends it back to be corrected
                    if not new_headline:
                        for right_headline in right_headlines:
                            if right_headline not in headlines:
                                return "Wrong headlines", headlines, sample_dict

                    for headline_index, headline in enumerate(headlines):
                        if new_headline:
                            headline = new_headline[headline]

                        if headline == config[config_headline]["source_plates"]:
                            source_plates_index = headline_index
                        elif headline == config[config_headline]["destination_plates"]:
                            destination_plates_index = headline_index
                        elif headline == config[config_headline]["source_well"]:
                            source_well_index = headline_index
                        elif headline == config[config_headline]["destination_well"]:
                            destination_well_index = headline_index
                else:
                    line = lines.split(splitter[split_indicator])
                    for data_index, data in enumerate(line):
                        if data_index == source_plates_index:
                            temp_source_plate = data
                        elif data_index == destination_plates_index:
                            temp_destination_plate = data
                        elif data_index == source_well_index:
                            temp_source_well = data
                        elif data_index == destination_well_index:
                            temp_destination_well = data

                    if temp_destination_plate not in all_destination_plates:
                        all_destination_plates.append(temp_destination_plate)

                    try:
                        sample_dict[temp_destination_plate]
                    except KeyError:
                        sample_dict[temp_destination_plate] = {}
                        sample_dict[temp_destination_plate][temp_destination_well] = {"source_plate": temp_source_plate,
                                                                                      "source_well": temp_source_well,
                                                                                      "smiles": "",
                                                                                      "compound_id": ""}
                    else:
                        sample_dict[temp_destination_plate][temp_destination_well] = {"source_plate": temp_source_plate,
                                                                                      "source_well": temp_source_well,
                                                                                      "smiles": "",
                                                                                      "compound_id": ""}

        return "done", headlines, sample_dict, all_destination_plates


class CSVConverter:

    def __str__(self):
        """Convert CSV files to other formate"""

    @staticmethod
    def mp_in_to_out(path, tube_files, mp_name, trans_vol):
        """
        Convert Tube files from MotherPlate Files to a CSV file for PlateButler, incase some files are missing.

        :param path: Main Output folder
        :type path: str
        :param tube_files: file for Tubes
        :type tube_files: list
        :param mp_name: MotherPlate main name.
        :type mp_name: str
        :return: A CSV file for PlateButler
        """
        tube_dict = {}
        for files in tube_files:
            CSVReader.pb_tube_files_ind(files, tube_dict)

        _, pb_mp_output = mpg(tube_dict, mp_name, trans_vol)
        CSVWriter.mp_to_pb_writer(path, pb_mp_output)


if __name__ == "__main__":
    ...

