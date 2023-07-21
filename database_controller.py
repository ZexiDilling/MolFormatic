from chem_operators import ChemOperators
from math import ceil

from database_handler import DataBaseFunctions
from sdf_handler import SDFReader
from csv_handler import CSVReader
from random import sample
from operator import countOf
from file_handler import get_file_list, file_list_distributor, move_files
from xml_handler import xml_controller
from lc_data_handler import LCMSHandler


class FetchData:
    """
    :param config: The config data, from the config file
    :type config: configparser.ConfigParser
    """
    def __init__(self, config):

        self.config = config
        self.dbf = DataBaseFunctions(config)
        self.dbf.create_connection()
        self.co = ChemOperators()
        self.database = config["Database"]["database"]

    def __str__(self):
        """
        Gets data from the database

        :return: rows of Data
        :rtype: dict
        """

    def get_number_of_rows(self, table):
        return self.dbf.number_of_rows(table)

    @staticmethod
    def _list_generator(sample_amount, compound_list_temp):
        """
        Generate a random list of compounds based on amount of samples, and total list of possible compounds

        :param sample_amount: Amount of compounds needed
        :type sample_amount: int
        :param compound_list_temp: List of all compounds that fits criteria
        :type compound_list_temp: list
        :return: A list of random compounds
        :rtype: list
        """
        compound_list_final = set()
        if compound_list_temp:
            while len(compound_list_final) < sample_amount:
                if len(compound_list_temp) <= sample_amount:
                    sample_amount = len(compound_list_temp)
                compound_list_final = sample(compound_list_temp, sample_amount)

            return compound_list_final
        else:
            return None

    @staticmethod
    def _rows_to_list(rows):
        """
        Makes a list of all compounds

        :param rows: All rows from the database
        :type rows: dict
        :return: A list of all compounds
        :rtype: list
        """

        return [rows[data]["compound_id"] for data in rows]

    @staticmethod
    def _liquid_warning(rows, min_volume):
        """
        Create a list of compounds that are running low or out

        :param rows: Rows of data from the Database
        :type rows: dict
        :param min_volume: Minimum volume needed per transferee. Will also determent the warning threshold
        :type min_volume: int
        :return: A list of compounds that are under the warning threshold with amount of compound left
        :rtype: dict
        """
        dead_volume = 1
        warning_level = min_volume * 3 + dead_volume
        liquid_warning = {}
        for compound_id in rows:
            if rows[compound_id]["volume"] < warning_level:
                liquid_warning[compound_id] = rows[compound_id]["volume"]

        return liquid_warning

    @staticmethod
    def _row_plate_data(rows, compound_id):
        """
        Generate dicts and list of different kind of data, based on data from the Database

        :param rows: Data from the Database
        :type rows: dict
        :param compound_id: List of random compounds.
        :type compound_id: list
        :return: Dict of data, Dict of plates, List of plates and list of plates with amount of each plate.
        :rtype: dict, dict, list, list
        """
        row_data = {}
        mp_data = {}
        plate_set = set()
        plate_count = []
        for row_counter in rows:
            temp_compound_id = rows[row_counter]["compound_id"]
            if temp_compound_id in compound_id:
                row_data[row_counter] = [rows[row_counter]["mp_barcode"], rows[row_counter]["mp_well"], rows[row_counter]["compound_id"]]
                mp_data.setdefault(rows[row_counter]["mp_barcode"], []).append([rows[row_counter]["mp_well"], rows[row_counter]["compound_id"]])
                plate_set.add(rows[row_counter]["mp_barcode"])
                plate_count.append(rows[row_counter]["mp_barcode"])

        return row_data, mp_data, sorted(plate_set), sorted(plate_count)

    @staticmethod
    def _plate_counter(plates):
        """
        Makes a dict of plates and how many of the same plate there is

        :param plates: List of all the plates
        :type plates: list
        :return: Dicts with plates and amount of repeats
        :rtype: dict
        """
        return dict((i, countOf(plates, i)) for i in set(plates))

    def _plate_mapping(self, table, plates, barcode_name):
        """
        Makes a dict of plates, with information, like where they are located

        :param table: The table the plates are in
        :type table: str
        :param plates: list of all the plates without repeats
        :type plates: list
        :param barcode_name: Headlines for the plates in the table
        :type barcode_name: str
        :return: A list of plates where they are located
        :rtype: list
        """
        plate_list = []
        for plate in plates:
            plate_temp = self.dbf.find_data_single_lookup(table, plate, barcode_name)
            plate_list.append(plate_temp)

        return plate_list

    def sub_structure_search(self, search_limiter, smiles, threshold,
                             methode="sub_structure_general", table="compound_main"):
        """
        Structure search controller. Controls witch function to use and what table to look into

        :param smiles: The smiles codes that needs to be compared
        :type smiles: str
        :param threshold: The threshold for how similar the compound needs to be to the main smiles code
        :type threshold: int
        :param search_limiter: A dict with different values when searching for compounds
        :type search_limiter: dict
        :param methode: What search method to use
        :type methode: str
        :param table: What table to find the data
        :type table: str
        :return: Rows of data from the Database with compound that fits the sub structure search criteria.
        :rtype: dict
        """
        if table == "compound_main":
            rows = self.dbf.return_table_data(table, search_limiter)
        elif table == "join_main_mp":
            rows = self.dbf.join_table_controller(search_limiter)

        rows = self.co.structure_search(methode, threshold, rows, smiles)
        return rows

    def data_search(self, table, search_limiter):
        """
        Gets a list of all compounds in the database

        :param search_limiter: A dict over values to search for in the db
        :type search_limiter: dict or int
        :param table: What table to look into.
        :type table: str
        :return: rows of data from the Database
        :rtype: dict
        """

        if table == "join_main_mp":
            rows = self.dbf.join_table_controller(search_limiter)
        else:
            rows = self.dbf.return_table_data(table, search_limiter)

        return rows

    def list_to_rows(self, compound_list, table="compound_main"):
        """
        Gets a list of compounds and finds the corresponding rows in the main table.

        :param compound_list: A list of the compounds
        :type compound_list: list
        :param table: The table to look for the data, should always be the main table
        :type table: str
        :return: the rows for the compounds
        :rtype: dict
        """
        rows = {}
        for compound_id in compound_list:
            temp_dict = self.dbf.find_data_single_lookup(table, compound_id, "compound_id")
            for key, value in temp_dict.items():
                rows[key] = value

        return rows

    def list_limiter(self, sample_amount, min_mp, samples_per_plate, table, sub_search, sub_search_methode, smiles,
                     threshold, ignore_active, plated_compounds, search_limiter):
        """
        Limits the list of compounds based on different criteria.

        :param sample_amount: how many samples is needed
        :type sample_amount: int
        :param min_mp: If all the samples needs  to be from as few plates as possible
        :type min_mp: bool
        :param samples_per_plate: amount of samples per plate
        :type samples_per_plate: int
        :param table: What table to find the data in
        :type table: str
        :param sub_search: True/False if it should be used or not
        :type sub_search: bool
        :param sub_search_methode: What structure method should be used if sub_search is True
        :type sub_search_methode: str or None
        :param smiles: The smiles code to be used for the structure search.
        :type smiles: str or None
        :param threshold: How similar the compound needs to be to the smiles code
        :type threshold: int or  None
        :param ignore_active: If the list needs to take into account compounds already in MotherPlates
        :type ignore_active: bool
        :param plated_compounds: list of compounds in MotherPlates
        :type plated_compounds: list or None
        :param search_limiter: A dict over values to search for in the db
        :type search_limiter: dict
        :return: all_the_things is a list containing the following:
            - limited_compound_list: A list of compounds
            - warnings: A dict of warnings for compounds that are close to being empty
            - row_data: A dict of all the data from the tables
            - mp_data: The MotherPlate data in dict formate
            - mp_mapping: A dict with MotherPlate information
            - plate_count: The amount of plates
        :rtype: list:
            - list
            - dict
            - dict
            - dict
            - dict
            - int
        """
        compound_search = search_limiter[self.config["Tables"]["compound_main"]]

        if compound_search["origin_id"]["use"]:
            origin_tables = self.config["Tables"]["compound_source"]
            rows = self.data_search(origin_tables, search_limiter[origin_tables])
            origin_id = []
            for row in rows:
                origin_id.append(rows[row]["ac_id"])
            compound_search["origin_id"]["value"] = origin_id

            if not compound_search["origin_id"]["value"]:
                return None

        if sub_search:
            rows = self.sub_structure_search(compound_search, smiles, threshold, sub_search_methode, table)
        else:
            if table == "join_main_mp":
                temp_table = "compound_mp"
            else:
                temp_table = table
            rows = self.data_search(temp_table, compound_search)
        if compound_search["volume"]["use"]:
            warnings = self._liquid_warning(rows, compound_search["volume"]["value"])
        else:
            warnings = None

        full_compound_list = self._rows_to_list(rows)

        if not ignore_active and plated_compounds != []:
            try:
                full_compound_list = [compound for compound in full_compound_list if compound not in plated_compounds]

            except TypeError:
                return None

        if sample_amount:

            if min_mp:
                mp_plate_amount = sample_amount / samples_per_plate
                mp_plate_amount = ceil(mp_plate_amount)
                temp_table = "mp_plates"
                row_data = self.data_search(temp_table, None)
                all_mps = [mp for mp in row_data]
                list_of_mps = self._list_generator(mp_plate_amount, all_mps)
                temp_compound_list = []
                for mp in list_of_mps:
                    search_limiter = {
                        "compound_mp": {"value": mp, "operator": "=", "target_column": "mp_barcode", "use": True}}
                    rows = self.return_table_data("compound_mp", search_limiter)
                    for row in rows:
                        temp_compound_list.append(rows[row]["compound_id"])

                limited_compound_list = self._list_generator(sample_amount, temp_compound_list)
            else:
                limited_compound_list = self._list_generator(sample_amount, full_compound_list)

        else:
            limited_compound_list = full_compound_list
        if not limited_compound_list:
            return None
        if table == "compound_main":
            all_the_things = limited_compound_list, warnings
        elif table == "join_main_mp":
            row_data, mp_data, plate_set, plate_count = self._row_plate_data(rows, limited_compound_list)
            mp_mapping = self._plate_mapping("mp_plates", plate_set, "mp_barcode")
            plate_count = self._plate_counter(plate_count)
            all_the_things = limited_compound_list, warnings, row_data, mp_data, mp_mapping, plate_count

        return all_the_things

    def get_tubes(self):
        ...

    def get_mother_plates(self):
        ...

    def info_df_generator(self, compounds):
        ...

    # def find_compound_id(self, table_name, item_id_1, item_id_2, item_name_1, item_name_2):
    #     """
    #     Find the compound ID based on plate for MotherPlates! ! !
    #     :param table: Table to look in
    #     :param plate_barcode: Barcode for the Plate
    #     :param well: Well id for the compound
    #     :return: Compound ID
    #     """
    #     rows = self.dbf.find_data(table_name, item_id_1, item_id_2, item_name_1, item_name_2)
    #     return rows[0][1]


class AddData:
    def __init__(self, config):
        """
        :param config: The config handler, with all the default information in the config file.
        :type config: configparser.ConfigParser
        """
        self.database = config["Database"]["database"]
        self.dbf = DataBaseFunctions(config)
        fd = FetchData(config)
        self.sdf_r = SDFReader(config, fd, self.dbf)
        self.csv_r = CSVReader()
        self.co = ChemOperators()

    def __str__(self):
        """
        Main controller for adding data to the database

        :return: Data added to the database
        """

    @staticmethod
    def _key_name_chang(temp_data):
        """
        Change the name for the keys from the files to fit with the Database naming scheme.

        :param temp_data: A dict with the keys that needs to have a name change
        :type temp_data: dict
        :return: The dict with some of its keys having changed name
        :rtype: dict
        """
        # Create a new dictionary "temp_data" that contains the same key-value pairs as the input dictionary,
        # but with the keys "amount" replaced by "volume"
        temp_data = {key.replace("amount", "volume"): value for key, value in temp_data.items()}

        # Further update "temp_data" by replacing keys "barcode" with "compound_id"
        temp_data = {key.replace("barcode", "compound_id"): value for key, value in temp_data.items()}

        # Return the updated dictionary
        return temp_data

    @staticmethod
    def _popper(table, data_dict):
        """
        Remove data from a dict

        :param table: "Pops" different values depending on what table is being worked on.
        :type table: str
        :param data_dict: Dict of data
        :type data_dict: dict
        :return: A dict of data with less key-values
        :rtype: dict
        """
        # Create a reference to the input dictionary "data_dict" in a new variable "temp_dict"
        temp_dict = data_dict

        # Create a dictionary of ban_lists based on the input table name
        ban_lists = {"compound_mp": ["SourceWell", "SourceBarcode"],
                     "purity": ["Peak_Info", "mass", "time_date", "wavelength"]}

        # Check if the input table name is a key in the ban_lists dictionary
        if table in ban_lists:
            # Get the ban_list for the input table name
            ban_list = ban_lists[table]
            # Check if the first key in the ban list is in the data_dict
            if ban_list[0] in data_dict:
                # Loop through the ban list and remove each key from the data_dict
                for clm in ban_list:
                    data_dict.pop(clm)

        # Return the updated temp_dict
        return temp_dict

    @staticmethod
    def _re_ordering_dict(dict_wrong_order, destination_table):
        """
        Re-arrange the order of the dict to work with the tables in the database

        :param dict_wrong_order: The dicts that needs to be re-arrange
        :type dict_wrong_order: dict
        :param destination_table: What tables the data needs to fit
        :type destination_table: str
        :return: the Dict_wrong_order in the right order to fit the Database
        :rtype:dict
        """
        # Create a dictionary of key orders based on the input destination table name
        key_orders = {
            "compound_main": ["compound_id", "smiles", "png", "volume", "concentration", "ac_id", "origin_id"],
            "compound_mp": ["Row_Counter", "DestinationBarcode", "compound_id", "DestinationWell", "Volume", "Date"],
            "compound_dp": ["Row_Counter", "DestinationBarcode", "compound_id", "DestinationWell", "Volume", "Date",
                            "SourceBarcode", "SourceWell"],
            "purity": ["row_counter", "compound_id", "experiment", "result_max", "result_max_ion", "result_total"]}

        # Check if the input destination table name is a key in the key_orders dictionary
        if destination_table in key_orders:
            # Get the key order for the input destination table name
            key_order = key_orders[destination_table]
            # Return a new dictionary created by looping through the key order and getting the corresponding value
            # from the input dictionary "dict_wrong_order"
            return {value: dict_wrong_order[value] for value in key_order}

    @staticmethod
    def _exp_dict_creator(compound_data, exp_count, exp_type, responsible):
        """
        Generates experiments enteries.
        Have not been tested!!!

        :param compound_data: The data for the compounds that needs to be added to the table
        :type compound_data: dict
        :param exp_count: The experiment counter, to make sure that each experiment have its own row
        :type exp_count: int
        :param exp_type: What type of experiment the data is from.
        :type exp_type: str
        :param responsible: The person that ran the experiment
        :type responsible: str
        :return:
            - experiment_dict: A dict of the experimental data
            - exp_count: The count for the experimental data
        :rtype:
            - dict
            - int
        """
        # Get the first time_date from the compound data dictionary
        time_date = next(iter(compound_data.values()))["time_date"][0]

        # Increment the experiment count
        exp_count += 1

        # Create a dictionary for the experiment details
        experiment_dict = {"exp_id": exp_count, "type": exp_type, "responsible": responsible, "date": time_date}

        # Return the experiment dictionary and experiment count
        return experiment_dict, exp_count

    @staticmethod
    def _purity_unpacker(compound_data, exp_count):
        """
        Unpack purity data and formate them in a way that can be added to a table

        :param compound_data: The compound purity data
        :type compound_data: dict
        :param exp_count: The experiment counter
        :type exp_count: int
        :return: Purity data in a dict formate
        :rtype: dict
        """
        for compound in compound_data:
            temp_purity = []
            temp_ions = []
            for peaks in compound_data[compound]["Peak_Info"]:
                temp_purity_value = compound_data[compound]["Peak_Info"][peaks]["purity"]
                temp_purity.append(temp_purity_value)
                temp_ion = list(compound_data[compound]["Peak_Info"][peaks].keys())[0]
                temp_ions.append(temp_ion)

            purity_max = max(temp_purity)
            temp_total = sum(temp_purity)

            for index, data in enumerate(temp_purity):
                if data == purity_max:
                    ion_max = temp_ions[index]

            compound_data[compound]["result_max"] = purity_max
            compound_data[compound]["result_max_ion"] = ion_max
            compound_data[compound]["result_total"] = temp_total
            compound_data[compound]["experiment"] = exp_count

        return compound_data

    def compound_main(self, file_list, table="compound_main"):
        """
        Adds data to the main table.

        :param file_list: list of sdf_file that contains compound and compound information
        :type file_list: list
        :param table: table to add the data. Should always be "compound_main"
        :type table: str
        :return: Compounds added to the main table.
        """
        for sdf_file in file_list:
            data = self.sdf_r.run(sdf_file)

            for compounds in data:
                temp_compounds_data = self._key_name_chang(data[compounds])
                temp_compounds_data["png"] = self.co.png_string(temp_compounds_data["smiles"])
                temp_compounds_data = self._re_ordering_dict(temp_compounds_data, table)
                temp_compounds_data["active"] = True
                self.dbf.add_records_controller(table, temp_compounds_data)

    def mother_plate(self, file_list, file_type, plate_table_name="mp_plates", destination_table="compound_mp",
                     source_table="compound_main", clm_id="compound_id"):
        """
        Reads a CSV file, Adds the data to "compound_mp" Table. Update the volume of "compound_main" table
        Adds the plates to "mp_plates" with locations.

        :param file_list: list of files that contains csv_files with the data that needs to be added.
        :type file_list: list
        :param file_type: What kind of CSV file the data is in.
        :type file_type: str
        :param plate_table_name: Name of tables  for the plates, should always be "mp_plates"
        :type plate_table_name: str
        :param destination_table: Where the compounds are going, should always be "compound_mp"
        :type destination_table: str
        :param source_table: Where the compounds are coming from, should always be "compound_main"
        :type source_table: str
        :param clm_id: The name of the compound-ID coloumn, should always be "compound_id"
        :type clm_id: str
        :return: Adds MotherPlates to the database and updated the main database with new volumes
        """
        test = 0
        for csv_file in file_list:
            data_dict, plates_dict = self.csv_r.csv_r_controller(csv_file, file_type)
            for plate in plates_dict:
                self.dbf.add_records_controller(plate_table_name, plates_dict[plate], test)
            for transferee in data_dict:
                new_dict = self._popper(destination_table, data_dict[transferee])
                new_dict = self._re_ordering_dict(new_dict, destination_table)
                new_dict["active"] = True
                new_dict["freeze_thaw"] = 0

                #Check if line is allready in the database:
                table = destination_table
                barcode_name = "mp_barcode"
                barcode = new_dict["DestinationBarcode"]
                id_name = "compound_id"
                id_number = new_dict["compound_id"]

                test_db = self.dbf.find_data_double_lookup(table, barcode, id_number, barcode_name, id_name)
                if not test_db:

                    self.dbf.add_records_controller(destination_table, new_dict, test)
                    # self.dbf.update_vol(source_table, data_dict[transferee]["Volume"],
                    #                     data_dict[transferee]["compound_id"], clm_id)
                    test += 1

    def daughter_plate(self, file_list, plate_table_name="dp_plates",
                       destination_table="compound_dp", source_table="compound_mp", clm_id="rowid"):
        """
        Locates all files in a folder with data from the Echo.
        Adds the data to "compound_dp" Table.
        Update the volume of "compound_mp" table
        Adds the plates to "dp_plates" with locations.

        :param file_list: with all the Echo files.
        :type file_list: list
        :param plate_table_name: Name of table for plates, should always be "dp_plates".
        :type plate_table_name: str
        :param destination_table: Where the compounds are going, should always be "compound_dp".
        :type destination_table: str
        :param source_table: Where the compounds are coming from, should always be "compound_mp".
        :type source_table: str
        :param clm_id: Header for clm where row_id is located for updating the right value.
        :type clm_id: str
        :return: Added data to the database and updates volumes.
        """

        data_dict, plates_dict = xml_controller(file_list)

        for plate in plates_dict:
            self.dbf.add_records_controller(plate_table_name, plates_dict[plate])

        # Find compound_ID based on source_plate_barcode and source_well from source_table. and add it to the data_dict
        for transferee in data_dict:
            row_data = self.dbf.find_data_double_lookup(source_table, data_dict[transferee]["SourceBarcode"],
                                                        data_dict[transferee]["SourceWell"], "mp_barcode", "mp_well")
            data_dict[transferee]["compoundID"] = row_data[0][1]

            # re-arrange the order to make it fit to the database setup
            new_dict = self._re_ordering_dict(data_dict[transferee], destination_table)
            new_dict["active"] = True

            # Adding data to the database
            self.dbf.add_records_controller(destination_table, new_dict)

            # Finding row_id for making sure that it is the right data that gets updated...
            # This will fail with multiple copies of the same compounds in the same MotherPlate... ARG!!!!!
            row = self.dbf.find_data_double_lookup(source_table, data_dict[transferee]["SourceBarcode"],
                                                   data_dict[transferee]["compoundID"], "mp_barcode", "compound_id")
            row_id = row[0][0]

            # Update the volume based on row_id.
            self.dbf.update_vol(source_table, data_dict[transferee]["Volume"],
                                row_id, clm_id)

    def purity_data(self, file_list, responsible="PHCH", exp_type="purity", experiment_tabel="experiment",
                    destination_table="purity"):
        """
        Add data to the experiment table, to show when the experiment is done ??
        Add data to the purity table, to have a table over purity, with a reference to the experiment table.
        Update source table to make used plate inactive ?

        :param file_list: contains a list of all the compound_data files
        :type file_list: list
        :param responsible: The responsible person / the person who have produced the data
        :type responsible: str
        :param exp_type: What experimental type the data is from. Should always be purity
        :type exp_type: str
        :param experiment_tabel: The table for experiments
        :type experiment_tabel: str
        :param destination_table: The table where the data goes, should always be purity
        :type destination_table: str
        :return: Adds purity data to the database
        """
        for compound_data in file_list:
            # Get numbers of rows from the experiment database.
            exp_count = self.dbf.number_of_rows(experiment_tabel)
            # Get a dict for the experiment, to make it fit into the table layout for the experimental table.
            experiment_dict, exp_count = self._exp_dict_creator(compound_data, exp_count, exp_type, responsible)
            # add data to the experimental table.
            self.dbf.add_records_controller(experiment_tabel, experiment_dict)
            # unpack the purity data to single values, and add the last columns for the data to fit with the table
            compound_data = self._purity_unpacker(compound_data, exp_count)
            # loop over each compound, removed unnecessary keys, re-order the dict to fit with the table
            # add data to the table
            for compound in compound_data:
                new_dict = self._popper(destination_table, compound_data[compound])
                new_dict = self._re_ordering_dict(new_dict, destination_table)
                self.dbf.add_records_controller(destination_table, new_dict)
        move_files(file_list)

    def auto_controller(self, folder, file_type):
        """
        HAVE NOT BEEN TESTED!!! IS NOT WORKING!!!!!!
        Looks through a folder, and adds all the files to the database, depending on the file type.

        :param folder:
        :type folder: str
        :param file_type:
        :type file_type: str
        :return: Updates the database with data and moves files out of pending folder! ?
        """
        compound_list, mp_list, dp_list, purity_list, bio_list, full_file_list = file_list_distributor(folder)
        lc_h = LCMSHandler()
        if compound_list:
            self.compound_main(compound_list)
        if mp_list:
            self.mother_plate(mp_list, file_type)
        if dp_list:
            self.daughter_plate(dp_list)
        if purity_list:
            # using default values!!
            compound_info = lc_h.lc_controller(purity_list, False, 254, "all", 2500, 2.5, 0.25, "pos", 50000000)
            self.purity_data(compound_info)
        if bio_list:
            pass
        move_files(full_file_list)

    def add_controller(self, table, data, file_type=None):
        """
        Main access point to adding data to the Database.

        :param table: What table the data needs to be added to
        :type table: str
        :param data: The folder with all the data files / For purity data this is compound information!!
        :type data: str or dict
        :param file_type: What kind of file it is
        :type file_type: str
        :return: Update the database with new values
        """
        table_list = ["compound_main", "compound_mp", "compound_dp", "location_table"]

        if table == "auto":
            self.auto_controller(data, file_type)
        elif table == "purity_data":
            self.purity_data(data)
        elif table not in table_list:
            self.dbf.add_records_controller(table, data)

        else:
            file_list = get_file_list(data)

            if table == "compound_main":
                self.compound_main(file_list)
            elif table == "compound_mp":
                self.mother_plate(file_list, file_type)
            elif table == "compound_dp":
                self.daughter_plate(file_list)
            elif table == "location_table":
                pass
            move_files(file_list)






