import re
from csv_handler import CSVWriter
from math import ceil


class PlateDilution:
    def __init__(self, config, save_plates, pb_source_file, control_vol=25, control_conc=1000):
        self.pb_source_file = pb_source_file
        self.csv_w = CSVWriter()
        plate_types_values = [items for items in config["plate_types_values"].keys()]
        self.dead_volume = {}
        self.well_vol = {}

        for plates in plate_types_values:
            plate = plates.upper()
            self.dead_volume[plate] = float(config["plate_types_values"][plates].split(",")[1])
            self.well_vol[plate] = float(config["plate_types_values"][plates].split(",")[2])

        self.plate_type_dict = {}
        plate_types = [items for items in config["plate_types"].keys()]
        for plates in plate_types:
            plates = plates.upper()
            self.plate_type_dict[plates] = config["plate_types"][plates]

        self.control_vol = control_vol
        self.control_conc = control_conc
        # for sample_dict :
        # min vol that echo can do
        self.min_vol = 2.5
        # Chosen coz the value works...
        self.target_conc = 0.1
        # Possible states for a well to have, besides "sample"
        self.state_list = ["blank", "max", "minimum", "positive", "negative"]

        # If plates needs to have one dilution per plate, or save plates:
        self.save_plates = save_plates

    @staticmethod
    def _dilution_plates(concentration_counter, dilution_factor):
        """
        Calculates the plates that needs to be used for different dilutions
        :return:
        """
        plate_dict = {}
        source_plate_name_list = []
        plate_name = "LDV"
        df = []
        for x in range(1, concentration_counter+1):
            df.append(dilution_factor ** (x-1))
            if x % 2 == 0:
                plate_dict[f"{plate_name}_D{int((x / 2) - 1)}"] = df
                df = []
                source_plate_name_list.append(f"{plate_name}_D{int((x / 2) - 1)}")

        return plate_dict, source_plate_name_list

    @staticmethod
    def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
        return [int(text) if text.isdigit() else text.lower()
                for text in _nsre.split(s)]

    def _get_sample_wells(self, plate_layout, direction="row"):
        destination_wells = []
        other_well_types = {
            "blank": [],
            "max": [],
            "minimum": [],
            "positive": [],
            "negative": []
        }

        for counter in plate_layout:
            well = plate_layout["well_id"]
            if plate_layout[counter]["state"] == "sample":
                destination_wells.append(well)
            elif plate_layout[counter]["state"] != "empty":
                other_well_types[plate_layout[counter]["state"]].append(well)
        if direction != "row":
            return destination_wells, other_well_types
        else:
            for states in other_well_types:
                other_well_types[states] = sorted(other_well_types[states], key=self.natural_sort_key)
            return sorted(destination_wells, key=self.natural_sort_key), other_well_types

    @staticmethod
    def _vol_calculation(conc_init, dilution, min_vol, plate_dilution_factor, target_conc):
        final_conc = dilution/target_conc
        return ((min_vol*final_conc)/conc_init)/plate_dilution_factor

    def _sample_dict(self, sample_info_dict, replicate_samples_max, replicate_plate_sets, concentration_counter, dilution_factor,
                     source_plate_name_list, destination_wells, destination_plate_naming_scheme, dw_amount):
        dilution_dict = {}
        sample_dict = {}
        vol_needed = {}
        source_plate_counter = 0
        source_well_counter = {}
        mixer_check = {}
        sample_plate_amount = len(destination_wells)
        plate_list = []
        well_off_set_use = False
        volume_count = {}

        for plate_set in range(1, replicate_plate_sets+1):
            sample_dict[plate_set] = {}
            mixer_check[plate_set] = {}
            temp_well_list = []
            destination_plate_counter = 0
            sample_counter = 0
            for rep in range(1, replicate_samples_max+1):
                mixer_check[plate_set][rep] = {}
                for samples in sample_info_dict:
                    try:
                        volume_count[samples]
                    except KeyError:
                        volume_count[samples] = {}
                    if sample_info_dict[samples]["state"] == "sample":
                        if rep == 1:
                            sample_dict[plate_set][samples] = {}

                            # Keeps track of volume in source wells
                            if plate_set == 1:
                                vol_needed[samples] = {}
                                source_well_counter[samples] = {}
                                dilution_dict[samples] = {}
                        if sample_info_dict[samples]["replicates"] >= rep:
                            sample_dict[plate_set][samples][rep] = {}

                            if sample_info_dict[samples]["mix"]:
                                mix_len = len(sample_info_dict[samples]["combination"])
                            else:
                                mix_len = 1

                            for mix_index in range(mix_len):

                                for x in range(1, concentration_counter+1):
                                    dilution = dilution_factor ** (x - 1)

                                    try:
                                        volume_count[samples][dilution] + 10
                                    except KeyError:
                                        volume_count[samples][dilution] = 0

                                    try:
                                        mixer_check[plate_set][rep][dilution]
                                    except KeyError:
                                        mixer_check[plate_set][rep][dilution] = {}
                                        mixer_check[plate_set][rep][dilution]["tracker"] = []

                                    if sample_info_dict[samples]["combination"]:
                                        temp_mix_name = "".join(sorted([samples, sample_info_dict[samples]["combination"][
                                            mix_index]]))
                                    else:
                                        temp_mix_name = samples

                                    if temp_mix_name not in mixer_check[plate_set][rep][dilution]["tracker"]:
                                        mixer_check[plate_set][rep][dilution]["tracker"].append(temp_mix_name)
                                        mixer_check[plate_set][rep][dilution][temp_mix_name] = ""

                                    # Set up sample dict layout
                                    try:
                                        sample_dict[plate_set][samples][rep][dilution]
                                    except KeyError:
                                        sample_dict[plate_set][samples][rep][dilution] = {"sample_vol": float,
                                                                                          "source_plate": "",
                                                                                          "source_well": "",
                                                                                          "destination_plate": "",
                                                                                          "destination_well": [],
                                                                                          "plate_type": ""}

                                    # Get plate Type
                                    temp_plate_type = self.plate_type_dict[sample_info_dict[samples]["solvent"]]

                                    sample_dict[plate_set][samples][rep][dilution]["plate_type"] = temp_plate_type

                                    # Get source plate
                                    temp_source_plate = source_plate_name_list[source_plate_counter]
                                    sample_dict[plate_set][samples][rep][dilution]["source_plate"] = temp_source_plate

                                    try:
                                        vol_needed[samples][temp_source_plate]
                                    except KeyError:
                                        vol_needed[samples][temp_source_plate] = 0

                                    try:
                                        source_well_counter[samples][temp_source_plate]
                                    except KeyError:
                                        source_well_counter[samples][temp_source_plate] = 0

                                    # Get volume per transferee
                                    conc_init = sample_info_dict[samples]["conc"]

                                    # Sets the transferee volume, depending on what plate is used... not it does not!!!
                                    if self.save_plates and well_off_set_use:
                                        plate_dilution_factor = 100 ** (source_plate_counter+2)
                                    else:
                                        plate_dilution_factor = 100 ** source_plate_counter

                                    if plate_dilution_factor == 0:
                                        plate_dilution_factor = 1
                                    temp_vol = self._vol_calculation(conc_init, dilution, self.min_vol,
                                                                     plate_dilution_factor, self.target_conc)
                                    sample_dict[plate_set][samples][rep][dilution]["sample_vol"] = temp_vol
                                    volume_count[samples][dilution] += temp_vol

                                    # Check if there is volume enough in the source well.
                                    # Else It changes to a new well in the samples well list
                                    if vol_needed[samples][source_plate_name_list[source_plate_counter]] + temp_vol \
                                            >= self.well_vol[temp_plate_type] + self.dead_volume[temp_plate_type]:
                                        vol_needed[samples][source_plate_name_list[source_plate_counter]] = 0
                                        source_well_counter[samples][temp_source_plate] += 1
                                    else:
                                        vol_needed[samples][source_plate_name_list[source_plate_counter]] += temp_vol

                                    # Get source Well
                                    if well_off_set_use:
                                        well_set = 2
                                    else:
                                        well_set = 1
                                    # print(sample_info_dict[samples][f"well_{well_set}"])
                                    # print(temp_source_plate)
                                    # print(samples)
                                    # print(source_well_counter[samples][temp_source_plate])
                                    temp_source_well = sample_info_dict[samples][f"well_{well_set}"][source_well_counter[samples][temp_source_plate]]
                                    sample_dict[plate_set][samples][rep][dilution]["source_well"] = temp_source_well

                                    # Destination well and plate conditions:
                                    if concentration_counter - x >= sample_plate_amount - sample_counter:
                                        sample_counter = 0
                                        destination_plate_counter += 1

                                    # Get destination plate
                                    destination_plate = \
                                        f"{plate_set}_plate_{destination_plate_naming_scheme[destination_plate_counter]}"
                                    sample_dict[plate_set][samples][rep][dilution]["destination_plate"] = destination_plate

                                    if destination_plate not in plate_list:
                                        plate_list.append(destination_plate)

                                    # Get destination well
                                    if not mixer_check[plate_set][rep][dilution][temp_mix_name]:
                                        mixer_check[plate_set][rep][dilution][temp_mix_name] = \
                                            destination_wells[sample_counter]
                                        sample_counter += 1

                                    sample_dict[plate_set][samples][rep][dilution]["destination_well"].append(mixer_check[plate_set][rep][dilution][temp_mix_name])

                                    if x % 2 == 0:
                                        if self.pb_source_file and x < concentration_counter \
                                                and plate_set == 1 and rep == 1 and mix_index == 1:

                                            # Check how many dilution wells needs to be made
                                            if dw_amount == "Copy":
                                                temp_source_well_counter = len(sample_info_dict[samples]["well_1"])
                                            elif dw_amount == "Minimum":
                                                temp_source_well_counter = ceil(
                                                    sample_info_dict[samples]["volume"] / self.well_vol[
                                                        temp_plate_type])
                                            else:
                                                temp_source_well_counter = dw_amount

                                            dilution_dict[samples][dilution] = {}
                                            for i in range(temp_source_well_counter):
                                                dilution_dict[samples][dilution][i] = {"sample_vol": float,
                                                                                       "source_plate": "",
                                                                                       "source_well": "",
                                                                                       "destination_plate": "",
                                                                                       "destination_well": "",
                                                                                       "plate_type": ""}

                                                dilution_dict[samples][dilution][i]["sample_vol"] = sample_info_dict[samples]["volume"] * 0.1 * 1000
                                                dilution_dict[samples][dilution][i]["source_plate"] = temp_source_plate
                                                if i >= len(sample_info_dict[samples]["well_1"]):
                                                    i = i- len(sample_info_dict[samples]["well_1"])

                                                dilution_dict[samples][dilution][i]["source_well"] = sample_info_dict[samples]["well_1"][i]
                                                dilution_dict[samples][dilution][i]["plate_type"] = temp_plate_type

                                        source_plate_counter += 1
                                        # If saves plates, there will be two sets of dilutions per source plate.
                                        # with a factor 10.000 between them.
                                        if self.save_plates and x % 4 == 0:
                                            source_plate_counter -= 2
                                            well_off_set_use = True
                                        else:
                                            well_off_set_use = False
                                        if self.pb_source_file and x < concentration_counter \
                                                and plate_set == 1 and rep == 1 and mix_index == 1:
                                            if well_off_set_use:
                                                well_set = 2
                                            else:
                                                well_set = 1
                                            for i in range(temp_source_well_counter):

                                                dilution_dict[samples][dilution][i]["destination_plate"] = source_plate_name_list[
                                                    source_plate_counter]

                                                dilution_dict[samples][dilution][i]["destination_well"] = \
                                                sample_info_dict[samples][f"well_{well_set}"][i]

                                source_plate_counter = 0

        return sample_dict, plate_list, dilution_dict

    def _other_well_states(self, sample_dict, plate_list, other_well_types, sample_info_dict, source_plate_name_list):
        sample_dict["other_state"] = {}
        for samples in sample_info_dict:
            states = sample_info_dict[samples]["state"]
            source_well_counter = 0
            well_volume_tally = 0
            if states in self.state_list:
                sample_dict["other_state"][states] = {}
                for plates in plate_list:
                    sample_dict["other_state"][states][plates] = {}
                    for wells in other_well_types[states]:
                        temp_plate_type = self.plate_type_dict[sample_info_dict[samples]["solvent"]]
                        sample_dict["other_state"][states][plates][wells] = {
                            "concentration": self.control_conc,
                            "destination_barcode": plates,
                            "destination_well": wells,
                            "compound_id": samples,
                            "volume": self.control_vol,
                            "source_well": sample_info_dict[samples]["well_1"][source_well_counter],
                            "source_barcode": source_plate_name_list[0],
                            "sour_plate_type": temp_plate_type}

                        well_volume_tally += self.control_vol
                        if well_volume_tally > self.dead_volume[temp_plate_type] + self.well_vol[temp_plate_type]:
                            source_well_counter += 1
                            well_volume_tally = 0
        return sample_dict

    def pd_controller(self, sample_info_dict, replicate_samples_max, replicate_plate_sets, dilution_factor,
                      concentration_counter, plate_layout, destination_plate_naming_scheme, dw_amount):
        """

        :param sample_info_dict: All the information for each sample
        :type sample_info_dict: dict
        :param replicate_samples_max: The max number of replication of a single sample
        :type replicate_samples_max: int
        :param replicate_plate_sets: How many sets should be printed
        :type replicate_plate_sets: int
        :param dilution_factor: How much dilutions is there between each dilution step
        :type dilution_factor: int
        :param concentration_counter: How many concentrations needs to be between the two.
        :type concentration_counter: float
        :param plate_layout: Layout for the final plate
        :type plate_layout: dict
        :param destination_plate_naming_scheme: The naming scheme  for the plates.
        :type destination_plate_naming_scheme: list
        :return:dict"""
        plate_dict, source_plate_name_list = self._dilution_plates(concentration_counter, dilution_factor)
        destination_wells, other_well_types = self._get_sample_wells(plate_layout)

        # plate_layout_Check = ... Check if the layout fits with samples + controls to leave zero empty, for full dilution
        # sets. (Dilution set = all dilutions of a sample/mix of samples)
        # if len(destination_wells) % concentration_counter == 0:
        sample_dict, plate_list, dilution_dict = self._sample_dict(sample_info_dict, replicate_samples_max, replicate_plate_sets,
                                                    concentration_counter, dilution_factor, source_plate_name_list,
                                                    destination_wells, destination_plate_naming_scheme, dw_amount)

        sample_dict = self._other_well_states(sample_dict, plate_list, other_well_types, sample_info_dict,
                                              source_plate_name_list)

        return sample_dict, dilution_dict


        # else:
        #     print("Sample layout do not fit with amount of dilutions")


if __name__ == "__main__":
    ...
    # import configparser
    # config = configparser.ConfigParser()
    # config.read("config.ini")
    # sample_amount = 3
    # # sample_list = ["x1", "x2", "x3"]
    # # well_stuff = {"well": [], "mix": bool, "conc": float, "combination": list, "replicates": int, "solvent": str}
    # sample_info_dict = {"X1": {"state": "sample",
    #                            "well_1": ["A1", "A2", "A3"],
    #                            "well_2": ["D1", "D2", "D3"],
    #                            "mix": True,
    #                            "conc": 5,
    #                            "combination": ["X3"],
    #                            "replicates": 3,
    #                            "solvent": "DMSO"},
    #                "X2": {"state": "sample",
    #                       "well_1": ["B1", "B2", "B3"],
    #                       "well_2": ["E1", "E2", "E3"],
    #                       "mix": False,
    #                       "conc": 10,
    #                       "combination": [],
    #                       "replicates": 2,
    #                       "solvent": "DMSO"},
    #                "X3": {"state": "sample",
    #                       "well_1": ["C1", "C2", "C3"],
    #                       "well_2": ["F1", "F2", "F3"],
    #                       "mix": True,
    #                       "conc": 10,
    #                       "combination": ["X1"],
    #                       "replicates": 3,
    #                       "solvent": "DMSO"}}
    # replicate_samples_max = 3
    # replicate_plate_sets = 9
    # concentration_initial = 10
    # # concentration_low = 1
    # # concentration_high = 100000
    # dilution_factor = 10
    # final_conc = 10
    # concentration_counter = 6
    # # Plate layout may only have full sets of sample dilution. add controls enough to fill it up.
    # plate_layouts = {"Standard_384":
    #                      {"well_layout": {"1": {"well_id": "A1", "state": "empty", "colour": "blue"}, "2": {"well_id": "B1", "state": "empty", "colour": "blue"}, "3": {"well_id": "C1", "state": "empty", "colour": "blue"}, "4": {"well_id": "D1", "state": "empty", "colour": "blue"}, "5": {"well_id": "E1", "state": "empty", "colour": "blue"}, "6": {"well_id": "F1", "state": "empty", "colour": "blue"}, "7": {"well_id": "G1", "state": "empty", "colour": "blue"}, "8": {"well_id": "H1", "state": "empty", "colour": "blue"}, "9": {"well_id": "I1", "state": "empty", "colour": "blue"}, "10": {"well_id": "J1", "state": "empty", "colour": "blue"}, "11": {"well_id": "K1", "state": "empty", "colour": "blue"}, "12": {"well_id": "L1", "state": "empty", "colour": "blue"}, "13": {"well_id": "M1", "state": "empty", "colour": "blue"}, "14": {"well_id": "N1", "state": "empty", "colour": "blue"}, "15": {"well_id": "O1", "state": "empty", "colour": "blue"}, "16": {"well_id": "P1", "state": "empty", "colour": "blue"}, "17": {"well_id": "A2", "state": "empty", "colour": "blue"}, "18": {"well_id": "B2", "state": "max", "colour": "purple"}, "19": {"well_id": "C2", "state": "max", "colour": "purple"}, "20": {"well_id": "D2", "state": "max", "colour": "purple"}, "21": {"well_id": "E2", "state": "max", "colour": "purple"}, "22": {"well_id": "F2", "state": "max", "colour": "purple"}, "23": {"well_id": "G2", "state": "max", "colour": "purple"}, "24": {"well_id": "H2", "state": "max", "colour": "purple"}, "25": {"well_id": "I2", "state": "max", "colour": "purple"}, "26": {"well_id": "J2", "state": "max", "colour": "purple"}, "27": {"well_id": "K2", "state": "max", "colour": "purple"}, "28": {"well_id": "L2", "state": "max", "colour": "purple"}, "29": {"well_id": "M2", "state": "max", "colour": "purple"}, "30": {"well_id": "N2", "state": "max", "colour": "purple"}, "31": {"well_id": "O2", "state": "empty", "colour": "blue"}, "32": {"well_id": "P2", "state": "empty", "colour": "blue"}, "33": {"well_id": "A3", "state": "empty", "colour": "blue"}, "34": {"well_id": "B3", "state": "sample", "colour": "orange"}, "35": {"well_id": "C3", "state": "sample", "colour": "orange"}, "36": {"well_id": "D3", "state": "sample", "colour": "orange"}, "37": {"well_id": "E3", "state": "sample", "colour": "orange"}, "38": {"well_id": "F3", "state": "sample", "colour": "orange"}, "39": {"well_id": "G3", "state": "sample", "colour": "orange"}, "40": {"well_id": "H3", "state": "sample", "colour": "orange"}, "41": {"well_id": "I3", "state": "sample", "colour": "orange"}, "42": {"well_id": "J3", "state": "sample", "colour": "orange"}, "43": {"well_id": "K3", "state": "sample", "colour": "orange"}, "44": {"well_id": "L3", "state": "sample", "colour": "orange"}, "45": {"well_id": "M3", "state": "sample", "colour": "orange"}, "46": {"well_id": "N3", "state": "sample", "colour": "orange"}, "47": {"well_id": "O3", "state": "empty", "colour": "blue"}, "48": {"well_id": "P3", "state": "empty", "colour": "blue"}, "49": {"well_id": "A4", "state": "empty", "colour": "blue"}, "50": {"well_id": "B4", "state": "sample", "colour": "orange"}, "51": {"well_id": "C4", "state": "sample", "colour": "orange"}, "52": {"well_id": "D4", "state": "sample", "colour": "orange"}, "53": {"well_id": "E4", "state": "sample", "colour": "orange"}, "54": {"well_id": "F4", "state": "sample", "colour": "orange"}, "55": {"well_id": "G4", "state": "sample", "colour": "orange"}, "56": {"well_id": "H4", "state": "sample", "colour": "orange"}, "57": {"well_id": "I4", "state": "sample", "colour": "orange"}, "58": {"well_id": "J4", "state": "sample", "colour": "orange"}, "59": {"well_id": "K4", "state": "sample", "colour": "orange"}, "60": {"well_id": "L4", "state": "sample", "colour": "orange"}, "61": {"well_id": "M4", "state": "sample", "colour": "orange"}, "62": {"well_id": "N4", "state": "sample", "colour": "orange"}, "63": {"well_id": "O4", "state": "empty", "colour": "blue"}, "64": {"well_id": "P4", "state": "empty", "colour": "blue"}, "65": {"well_id": "A5", "state": "empty", "colour": "blue"}, "66": {"well_id": "B5", "state": "sample", "colour": "orange"}, "67": {"well_id": "C5", "state": "sample", "colour": "orange"}, "68": {"well_id": "D5", "state": "sample", "colour": "orange"}, "69": {"well_id": "E5", "state": "sample", "colour": "orange"}, "70": {"well_id": "F5", "state": "sample", "colour": "orange"}, "71": {"well_id": "G5", "state": "sample", "colour": "orange"}, "72": {"well_id": "H5", "state": "sample", "colour": "orange"}, "73": {"well_id": "I5", "state": "sample", "colour": "orange"}, "74": {"well_id": "J5", "state": "sample", "colour": "orange"}, "75": {"well_id": "K5", "state": "sample", "colour": "orange"}, "76": {"well_id": "L5", "state": "sample", "colour": "orange"}, "77": {"well_id": "M5", "state": "sample", "colour": "orange"}, "78": {"well_id": "N5", "state": "sample", "colour": "orange"}, "79": {"well_id": "O5", "state": "empty", "colour": "blue"}, "80": {"well_id": "P5", "state": "empty", "colour": "blue"}, "81": {"well_id": "A6", "state": "empty", "colour": "blue"}, "82": {"well_id": "B6", "state": "sample", "colour": "orange"}, "83": {"well_id": "C6", "state": "sample", "colour": "orange"}, "84": {"well_id": "D6", "state": "sample", "colour": "orange"}, "85": {"well_id": "E6", "state": "sample", "colour": "orange"}, "86": {"well_id": "F6", "state": "sample", "colour": "orange"}, "87": {"well_id": "G6", "state": "sample", "colour": "orange"}, "88": {"well_id": "H6", "state": "sample", "colour": "orange"}, "89": {"well_id": "I6", "state": "sample", "colour": "orange"}, "90": {"well_id": "J6", "state": "sample", "colour": "orange"}, "91": {"well_id": "K6", "state": "sample", "colour": "orange"}, "92": {"well_id": "L6", "state": "sample", "colour": "orange"}, "93": {"well_id": "M6", "state": "sample", "colour": "orange"}, "94": {"well_id": "N6", "state": "sample", "colour": "orange"}, "95": {"well_id": "O6", "state": "empty", "colour": "blue"}, "96": {"well_id": "P6", "state": "empty", "colour": "blue"}, "97": {"well_id": "A7", "state": "empty", "colour": "blue"}, "98": {"well_id": "B7", "state": "sample", "colour": "orange"}, "99": {"well_id": "C7", "state": "sample", "colour": "orange"}, "100": {"well_id": "D7", "state": "sample", "colour": "orange"}, "101": {"well_id": "E7", "state": "sample", "colour": "orange"}, "102": {"well_id": "F7", "state": "sample", "colour": "orange"}, "103": {"well_id": "G7", "state": "sample", "colour": "orange"}, "104": {"well_id": "H7", "state": "sample", "colour": "orange"}, "105": {"well_id": "I7", "state": "sample", "colour": "orange"}, "106": {"well_id": "J7", "state": "sample", "colour": "orange"}, "107": {"well_id": "K7", "state": "sample", "colour": "orange"}, "108": {"well_id": "L7", "state": "sample", "colour": "orange"}, "109": {"well_id": "M7", "state": "sample", "colour": "orange"}, "110": {"well_id": "N7", "state": "sample", "colour": "orange"}, "111": {"well_id": "O7", "state": "empty", "colour": "blue"}, "112": {"well_id": "P7", "state": "empty", "colour": "blue"}, "113": {"well_id": "A8", "state": "empty", "colour": "blue"}, "114": {"well_id": "B8", "state": "sample", "colour": "orange"}, "115": {"well_id": "C8", "state": "sample", "colour": "orange"}, "116": {"well_id": "D8", "state": "sample", "colour": "orange"}, "117": {"well_id": "E8", "state": "sample", "colour": "orange"}, "118": {"well_id": "F8", "state": "sample", "colour": "orange"}, "119": {"well_id": "G8", "state": "sample", "colour": "orange"}, "120": {"well_id": "H8", "state": "sample", "colour": "orange"}, "121": {"well_id": "I8", "state": "sample", "colour": "orange"}, "122": {"well_id": "J8", "state": "sample", "colour": "orange"}, "123": {"well_id": "K8", "state": "sample", "colour": "orange"}, "124": {"well_id": "L8", "state": "sample", "colour": "orange"}, "125": {"well_id": "M8", "state": "sample", "colour": "orange"}, "126": {"well_id": "N8", "state": "sample", "colour": "orange"}, "127": {"well_id": "O8", "state": "empty", "colour": "blue"}, "128": {"well_id": "P8", "state": "empty", "colour": "blue"}, "129": {"well_id": "A9", "state": "empty", "colour": "blue"}, "130": {"well_id": "B9", "state": "sample", "colour": "orange"}, "131": {"well_id": "C9", "state": "sample", "colour": "orange"}, "132": {"well_id": "D9", "state": "sample", "colour": "orange"}, "133": {"well_id": "E9", "state": "sample", "colour": "orange"}, "134": {"well_id": "F9", "state": "sample", "colour": "orange"}, "135": {"well_id": "G9", "state": "sample", "colour": "orange"}, "136": {"well_id": "H9", "state": "sample", "colour": "orange"}, "137": {"well_id": "I9", "state": "sample", "colour": "orange"}, "138": {"well_id": "J9", "state": "sample", "colour": "orange"}, "139": {"well_id": "K9", "state": "sample", "colour": "orange"}, "140": {"well_id": "L9", "state": "sample", "colour": "orange"}, "141": {"well_id": "M9", "state": "sample", "colour": "orange"}, "142": {"well_id": "N9", "state": "sample", "colour": "orange"}, "143": {"well_id": "O9", "state": "empty", "colour": "blue"}, "144": {"well_id": "P9", "state": "empty", "colour": "blue"}, "145": {"well_id": "A10", "state": "empty", "colour": "blue"}, "146": {"well_id": "B10", "state": "sample", "colour": "orange"}, "147": {"well_id": "C10", "state": "sample", "colour": "orange"}, "148": {"well_id": "D10", "state": "sample", "colour": "orange"}, "149": {"well_id": "E10", "state": "sample", "colour": "orange"}, "150": {"well_id": "F10", "state": "sample", "colour": "orange"}, "151": {"well_id": "G10", "state": "sample", "colour": "orange"}, "152": {"well_id": "H10", "state": "sample", "colour": "orange"}, "153": {"well_id": "I10", "state": "sample", "colour": "orange"}, "154": {"well_id": "J10", "state": "sample", "colour": "orange"}, "155": {"well_id": "K10", "state": "sample", "colour": "orange"}, "156": {"well_id": "L10", "state": "sample", "colour": "orange"}, "157": {"well_id": "M10", "state": "sample", "colour": "orange"}, "158": {"well_id": "N10", "state": "sample", "colour": "orange"}, "159": {"well_id": "O10", "state": "empty", "colour": "blue"}, "160": {"well_id": "P10", "state": "empty", "colour": "blue"}, "161": {"well_id": "A11", "state": "empty", "colour": "blue"}, "162": {"well_id": "B11", "state": "sample", "colour": "orange"}, "163": {"well_id": "C11", "state": "sample", "colour": "orange"}, "164": {"well_id": "D11", "state": "sample", "colour": "orange"}, "165": {"well_id": "E11", "state": "sample", "colour": "orange"}, "166": {"well_id": "F11", "state": "sample", "colour": "orange"}, "167": {"well_id": "G11", "state": "sample", "colour": "orange"}, "168": {"well_id": "H11", "state": "sample", "colour": "orange"}, "169": {"well_id": "I11", "state": "sample", "colour": "orange"}, "170": {"well_id": "J11", "state": "sample", "colour": "orange"}, "171": {"well_id": "K11", "state": "sample", "colour": "orange"}, "172": {"well_id": "L11", "state": "sample", "colour": "orange"}, "173": {"well_id": "M11", "state": "sample", "colour": "orange"}, "174": {"well_id": "N11", "state": "sample", "colour": "orange"}, "175": {"well_id": "O11", "state": "empty", "colour": "blue"}, "176": {"well_id": "P11", "state": "empty", "colour": "blue"}, "177": {"well_id": "A12", "state": "empty", "colour": "blue"}, "178": {"well_id": "B12", "state": "sample", "colour": "orange"}, "179": {"well_id": "C12", "state": "sample", "colour": "orange"}, "180": {"well_id": "D12", "state": "sample", "colour": "orange"}, "181": {"well_id": "E12", "state": "sample", "colour": "orange"}, "182": {"well_id": "F12", "state": "sample", "colour": "orange"}, "183": {"well_id": "G12", "state": "sample", "colour": "orange"}, "184": {"well_id": "H12", "state": "sample", "colour": "orange"}, "185": {"well_id": "I12", "state": "sample", "colour": "orange"}, "186": {"well_id": "J12", "state": "sample", "colour": "orange"}, "187": {"well_id": "K12", "state": "sample", "colour": "orange"}, "188": {"well_id": "L12", "state": "sample", "colour": "orange"}, "189": {"well_id": "M12", "state": "sample", "colour": "orange"}, "190": {"well_id": "N12", "state": "sample", "colour": "orange"}, "191": {"well_id": "O12", "state": "empty", "colour": "blue"}, "192": {"well_id": "P12", "state": "empty", "colour": "blue"}, "193": {"well_id": "A13", "state": "empty", "colour": "blue"}, "194": {"well_id": "B13", "state": "sample", "colour": "orange"}, "195": {"well_id": "C13", "state": "sample", "colour": "orange"}, "196": {"well_id": "D13", "state": "sample", "colour": "orange"}, "197": {"well_id": "E13", "state": "sample", "colour": "orange"}, "198": {"well_id": "F13", "state": "sample", "colour": "orange"}, "199": {"well_id": "G13", "state": "sample", "colour": "orange"}, "200": {"well_id": "H13", "state": "sample", "colour": "orange"}, "201": {"well_id": "I13", "state": "sample", "colour": "orange"}, "202": {"well_id": "J13", "state": "sample", "colour": "orange"}, "203": {"well_id": "K13", "state": "sample", "colour": "orange"}, "204": {"well_id": "L13", "state": "sample", "colour": "orange"}, "205": {"well_id": "M13", "state": "sample", "colour": "orange"}, "206": {"well_id": "N13", "state": "sample", "colour": "orange"}, "207": {"well_id": "O13", "state": "empty", "colour": "blue"}, "208": {"well_id": "P13", "state": "empty", "colour": "blue"}, "209": {"well_id": "A14", "state": "empty", "colour": "blue"}, "210": {"well_id": "B14", "state": "sample", "colour": "orange"}, "211": {"well_id": "C14", "state": "sample", "colour": "orange"}, "212": {"well_id": "D14", "state": "sample", "colour": "orange"}, "213": {"well_id": "E14", "state": "sample", "colour": "orange"}, "214": {"well_id": "F14", "state": "sample", "colour": "orange"}, "215": {"well_id": "G14", "state": "sample", "colour": "orange"}, "216": {"well_id": "H14", "state": "sample", "colour": "orange"}, "217": {"well_id": "I14", "state": "sample", "colour": "orange"}, "218": {"well_id": "J14", "state": "sample", "colour": "orange"}, "219": {"well_id": "K14", "state": "sample", "colour": "orange"}, "220": {"well_id": "L14", "state": "sample", "colour": "orange"}, "221": {"well_id": "M14", "state": "sample", "colour": "orange"}, "222": {"well_id": "N14", "state": "sample", "colour": "orange"}, "223": {"well_id": "O14", "state": "empty", "colour": "blue"}, "224": {"well_id": "P14", "state": "empty", "colour": "blue"}, "225": {"well_id": "A15", "state": "empty", "colour": "blue"}, "226": {"well_id": "B15", "state": "sample", "colour": "orange"}, "227": {"well_id": "C15", "state": "sample", "colour": "orange"}, "228": {"well_id": "D15", "state": "sample", "colour": "orange"}, "229": {"well_id": "E15", "state": "sample", "colour": "orange"}, "230": {"well_id": "F15", "state": "sample", "colour": "orange"}, "231": {"well_id": "G15", "state": "sample", "colour": "orange"}, "232": {"well_id": "H15", "state": "sample", "colour": "orange"}, "233": {"well_id": "I15", "state": "sample", "colour": "orange"}, "234": {"well_id": "J15", "state": "sample", "colour": "orange"}, "235": {"well_id": "K15", "state": "sample", "colour": "orange"}, "236": {"well_id": "L15", "state": "sample", "colour": "orange"}, "237": {"well_id": "M15", "state": "sample", "colour": "orange"}, "238": {"well_id": "N15", "state": "sample", "colour": "orange"}, "239": {"well_id": "O15", "state": "empty", "colour": "blue"}, "240": {"well_id": "P15", "state": "empty", "colour": "blue"}, "241": {"well_id": "A16", "state": "empty", "colour": "blue"}, "242": {"well_id": "B16", "state": "sample", "colour": "orange"}, "243": {"well_id": "C16", "state": "sample", "colour": "orange"}, "244": {"well_id": "D16", "state": "sample", "colour": "orange"}, "245": {"well_id": "E16", "state": "sample", "colour": "orange"}, "246": {"well_id": "F16", "state": "sample", "colour": "orange"}, "247": {"well_id": "G16", "state": "sample", "colour": "orange"}, "248": {"well_id": "H16", "state": "sample", "colour": "orange"}, "249": {"well_id": "I16", "state": "sample", "colour": "orange"}, "250": {"well_id": "J16", "state": "sample", "colour": "orange"}, "251": {"well_id": "K16", "state": "sample", "colour": "orange"}, "252": {"well_id": "L16", "state": "sample", "colour": "orange"}, "253": {"well_id": "M16", "state": "sample", "colour": "orange"}, "254": {"well_id": "N16", "state": "sample", "colour": "orange"}, "255": {"well_id": "O16", "state": "empty", "colour": "blue"}, "256": {"well_id": "P16", "state": "empty", "colour": "blue"}, "257": {"well_id": "A17", "state": "empty", "colour": "blue"}, "258": {"well_id": "B17", "state": "sample", "colour": "orange"}, "259": {"well_id": "C17", "state": "sample", "colour": "orange"}, "260": {"well_id": "D17", "state": "sample", "colour": "orange"}, "261": {"well_id": "E17", "state": "sample", "colour": "orange"}, "262": {"well_id": "F17", "state": "sample", "colour": "orange"}, "263": {"well_id": "G17", "state": "sample", "colour": "orange"}, "264": {"well_id": "H17", "state": "sample", "colour": "orange"}, "265": {"well_id": "I17", "state": "sample", "colour": "orange"}, "266": {"well_id": "J17", "state": "sample", "colour": "orange"}, "267": {"well_id": "K17", "state": "sample", "colour": "orange"}, "268": {"well_id": "L17", "state": "sample", "colour": "orange"}, "269": {"well_id": "M17", "state": "sample", "colour": "orange"}, "270": {"well_id": "N17", "state": "sample", "colour": "orange"}, "271": {"well_id": "O17", "state": "empty", "colour": "blue"}, "272": {"well_id": "P17", "state": "empty", "colour": "blue"}, "273": {"well_id": "A18", "state": "empty", "colour": "blue"}, "274": {"well_id": "B18", "state": "sample", "colour": "orange"}, "275": {"well_id": "C18", "state": "sample", "colour": "orange"}, "276": {"well_id": "D18", "state": "sample", "colour": "orange"}, "277": {"well_id": "E18", "state": "sample", "colour": "orange"}, "278": {"well_id": "F18", "state": "sample", "colour": "orange"}, "279": {"well_id": "G18", "state": "sample", "colour": "orange"}, "280": {"well_id": "H18", "state": "sample", "colour": "orange"}, "281": {"well_id": "I18", "state": "sample", "colour": "orange"}, "282": {"well_id": "J18", "state": "sample", "colour": "orange"}, "283": {"well_id": "K18", "state": "sample", "colour": "orange"}, "284": {"well_id": "L18", "state": "sample", "colour": "orange"}, "285": {"well_id": "M18", "state": "sample", "colour": "orange"}, "286": {"well_id": "N18", "state": "sample", "colour": "orange"}, "287": {"well_id": "O18", "state": "empty", "colour": "blue"}, "288": {"well_id": "P18", "state": "empty", "colour": "blue"}, "289": {"well_id": "A19", "state": "empty", "colour": "blue"}, "290": {"well_id": "B19", "state": "sample", "colour": "orange"}, "291": {"well_id": "C19", "state": "sample", "colour": "orange"}, "292": {"well_id": "D19", "state": "sample", "colour": "orange"}, "293": {"well_id": "E19", "state": "sample", "colour": "orange"}, "294": {"well_id": "F19", "state": "sample", "colour": "orange"}, "295": {"well_id": "G19", "state": "sample", "colour": "orange"}, "296": {"well_id": "H19", "state": "sample", "colour": "orange"}, "297": {"well_id": "I19", "state": "sample", "colour": "orange"}, "298": {"well_id": "J19", "state": "sample", "colour": "orange"}, "299": {"well_id": "K19", "state": "sample", "colour": "orange"}, "300": {"well_id": "L19", "state": "sample", "colour": "orange"}, "301": {"well_id": "M19", "state": "sample", "colour": "orange"}, "302": {"well_id": "N19", "state": "sample", "colour": "orange"}, "303": {"well_id": "O19", "state": "empty", "colour": "blue"}, "304": {"well_id": "P19", "state": "empty", "colour": "blue"}, "305": {"well_id": "A20", "state": "empty", "colour": "blue"}, "306": {"well_id": "B20", "state": "sample", "colour": "orange"}, "307": {"well_id": "C20", "state": "sample", "colour": "orange"}, "308": {"well_id": "D20", "state": "sample", "colour": "orange"}, "309": {"well_id": "E20", "state": "sample", "colour": "orange"}, "310": {"well_id": "F20", "state": "sample", "colour": "orange"}, "311": {"well_id": "G20", "state": "sample", "colour": "orange"}, "312": {"well_id": "H20", "state": "sample", "colour": "orange"}, "313": {"well_id": "I20", "state": "sample", "colour": "orange"}, "314": {"well_id": "J20", "state": "sample", "colour": "orange"}, "315": {"well_id": "K20", "state": "sample", "colour": "orange"}, "316": {"well_id": "L20", "state": "sample", "colour": "orange"}, "317": {"well_id": "M20", "state": "sample", "colour": "orange"}, "318": {"well_id": "N20", "state": "sample", "colour": "orange"}, "319": {"well_id": "O20", "state": "empty", "colour": "blue"}, "320": {"well_id": "P20", "state": "empty", "colour": "blue"}, "321": {"well_id": "A21", "state": "empty", "colour": "blue"}, "322": {"well_id": "B21", "state": "sample", "colour": "orange"}, "323": {"well_id": "C21", "state": "sample", "colour": "orange"}, "324": {"well_id": "D21", "state": "sample", "colour": "orange"}, "325": {"well_id": "E21", "state": "sample", "colour": "orange"}, "326": {"well_id": "F21", "state": "sample", "colour": "orange"}, "327": {"well_id": "G21", "state": "sample", "colour": "orange"}, "328": {"well_id": "H21", "state": "sample", "colour": "orange"}, "329": {"well_id": "I21", "state": "sample", "colour": "orange"}, "330": {"well_id": "J21", "state": "sample", "colour": "orange"}, "331": {"well_id": "K21", "state": "sample", "colour": "orange"}, "332": {"well_id": "L21", "state": "sample", "colour": "orange"}, "333": {"well_id": "M21", "state": "sample", "colour": "orange"}, "334": {"well_id": "N21", "state": "sample", "colour": "orange"}, "335": {"well_id": "O21", "state": "empty", "colour": "blue"}, "336": {"well_id": "P21", "state": "empty", "colour": "blue"}, "337": {"well_id": "A22", "state": "empty", "colour": "blue"}, "338": {"well_id": "B22", "state": "sample", "colour": "orange"}, "339": {"well_id": "C22", "state": "sample", "colour": "orange"}, "340": {"well_id": "D22", "state": "sample", "colour": "orange"}, "341": {"well_id": "E22", "state": "sample", "colour": "orange"}, "342": {"well_id": "F22", "state": "sample", "colour": "orange"}, "343": {"well_id": "G22", "state": "sample", "colour": "orange"}, "344": {"well_id": "H22", "state": "sample", "colour": "orange"}, "345": {"well_id": "I22", "state": "sample", "colour": "orange"}, "346": {"well_id": "J22", "state": "sample", "colour": "orange"}, "347": {"well_id": "K22", "state": "sample", "colour": "orange"}, "348": {"well_id": "L22", "state": "sample", "colour": "orange"}, "349": {"well_id": "M22", "state": "sample", "colour": "orange"}, "350": {"well_id": "N22", "state": "sample", "colour": "orange"}, "351": {"well_id": "O22", "state": "empty", "colour": "blue"}, "352": {"well_id": "P22", "state": "empty", "colour": "blue"}, "353": {"well_id": "A23", "state": "empty", "colour": "blue"}, "354": {"well_id": "B23", "state": "minimum", "colour": "yellow"}, "355": {"well_id": "C23", "state": "minimum", "colour": "yellow"}, "356": {"well_id": "D23", "state": "minimum", "colour": "yellow"}, "357": {"well_id": "E23", "state": "minimum", "colour": "yellow"}, "358": {"well_id": "F23", "state": "minimum", "colour": "yellow"}, "359": {"well_id": "G23", "state": "minimum", "colour": "yellow"}, "360": {"well_id": "H23", "state": "minimum", "colour": "yellow"}, "361": {"well_id": "I23", "state": "minimum", "colour": "yellow"}, "362": {"well_id": "J23", "state": "minimum", "colour": "yellow"}, "363": {"well_id": "K23", "state": "minimum", "colour": "yellow"}, "364": {"well_id": "L23", "state": "minimum", "colour": "yellow"}, "365": {"well_id": "M23", "state": "minimum", "colour": "yellow"}, "366": {"well_id": "N23", "state": "minimum", "colour": "yellow"}, "367": {"well_id": "O23", "state": "empty", "colour": "blue"}, "368": {"well_id": "P23", "state": "empty", "colour": "blue"}, "369": {"well_id": "A24", "state": "empty", "colour": "blue"}, "370": {"well_id": "B24", "state": "empty", "colour": "blue"}, "371": {"well_id": "C24", "state": "empty", "colour": "blue"}, "372": {"well_id": "D24", "state": "empty", "colour": "blue"}, "373": {"well_id": "E24", "state": "empty", "colour": "blue"}, "374": {"well_id": "F24", "state": "empty", "colour": "blue"}, "375": {"well_id": "G24", "state": "empty", "colour": "blue"}, "376": {"well_id": "H24", "state": "empty", "colour": "blue"}, "377": {"well_id": "I24", "state": "empty", "colour": "blue"}, "378": {"well_id": "J24", "state": "empty", "colour": "blue"}, "379": {"well_id": "K24", "state": "empty", "colour": "blue"}, "380": {"well_id": "L24", "state": "empty", "colour": "blue"}, "381": {"well_id": "M24", "state": "empty", "colour": "blue"}, "382": {"well_id": "N24", "state": "empty", "colour": "blue"}, "383": {"well_id": "O24", "state": "empty", "colour": "blue"}, "384": {"well_id": "P24", "state": "empty", "colour": "blue"}}, "plate_type": "plate_384"}}
    # sample_layout = str
    # sample_layout_distribution = "random"
    # destination_plate_naming_scheme = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
    # output_folder = "output_data"
    # save_plates = True
    #
    # from plate_formatting import plate_layout_re_formate
    # plate_layout = plate_layout_re_formate(plate_layouts["Standard_384"]["well_layout"])
    # # print(plate_layout)
    #
    # pd = PlateDilution(config, save_plates)
    # pd.controller(sample_info_dict, replicate_samples_max, replicate_plate_sets, dilution_factor,
    #               concentration_counter, plate_layout, destination_plate_naming_scheme, output_folder)
    #
