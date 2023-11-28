from database_controller import AddData, FetchData
from database_handler import DataBaseFunctions


def get_number_of_rows(config, table):
    dbf = DataBaseFunctions(config)
    return dbf.number_of_rows(table)


def update_database(config, table, data, file_type=None):
    """
    Update the database with data. Both for adding data to the different tables, but updating tables with new values

    :param data: Output file from the plate_butler system, tube files for the comPOUND freezer, or other data.
    :type data: str or dict
    :param file_type: If the file is a CSV file, it needs a specific file_type to know how to handle the file. as the
        CSV files are different depending on where they are coming from.
    :type file_type: str
    :param table: What table needs to be updated/where the data is going
    :type table: str
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :return: Updated database with values
    """

    ad_db = AddData(config)
    ad_db.add_controller(table, data, file_type)


def grab_compound_table_data(config, table, sample):
    dbf = DataBaseFunctions(config)
    Search_limiter = {table:
                          {"value": [sample],
                           "operator": "IN",
                           "target_column": "compound_id",
                           "use": True}
                      }

    return dbf.return_table_data(table, Search_limiter)


def database_to_table(config, table, headings):
    dbf = DataBaseFunctions(config)
    row_data = dbf.return_table_data(table, None)

    table_data = []
    for row in row_data:
        temp_table_row = []
        for data in row_data[row]:
            if data in headings:
                temp_table_row.append(row_data[row][data])
        table_data.append(temp_table_row)

    return table_data


def delete_records_from_database(config, table, headline, data_value):
    """
    Deletes a record from the database
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param table: What table are the plates in
    :type table: str
    :param data_value: The value of the thing you are looking for
    :type data_value: str
    :param headline: Headline for the column where  the data is, in the table
    :type headline: str
    :return:
    """
    dbf = DataBaseFunctions(config)
    dbf.delete_records(table, headline, data_value)


def rename_record_in_the_database(config, table, headline, old_value, new_value):
    """
    Deletes a record from the database
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param table: What table are the plates in
    :type table: str
    :param headline: Headline for the coloumn where  the data is, in the table
    :type headline: str
    :param old_value: The value that needs to be changed
    :type old_value: any
    :param new_value: What the value should be changed to
    :type new_value: any
    :return:
    """
    dbf = DataBaseFunctions(config)
    dbf.rename_record_value(table, headline, old_value, new_value)


def grab_table_data(config, table_name, search_limiter=None, single_row=None, data_value=None, headline=None,
                    search_list_clm=None, specific_rows=None):
    """
    Grabs data from tables
    :param config: The config handler, with all the default information in the config file.
    :type config: configparser.ConfigParser
    :param table_name: What table are the plates in
    :type table_name: str
    :param search_limiter: A dict to limit the search
    :type search_limiter: dict or list or None
    :param single_row: bool statement if it is a sling row look up
    :type single_row: bool
    :param data_value: The value of the thing you are looking for
    :type data_value: any or None
    :param headline: Headline for the column where  the data is, in the table
    :type headline: str or None
    :return:
    """

    if single_row:
        dbf = DataBaseFunctions(config)
        row_data = dbf.find_data_single_lookup(table_name, data_value, headline)
        return row_data
    else:
        fd = FetchData(config)
        rows = fd.data_search(table_name, search_limiter, search_list_clm, specific_rows)
        all_table_data = []
        headlines = []
        if rows:
            for row in rows:
                temp_data = []
                try:
                    rows[row]
                except TypeError:
                    print(f"table_name: {table_name}")
                    print(f"search_limiter: {search_limiter}")
                    print(f"search_list_clm: {search_list_clm}")
                    print(f"specific_rows: {specific_rows}")
                else:
                    for data in rows[row]:
                        headlines.append(data)
                        temp_data.append(rows[row][data])

                all_table_data.append(temp_data)

            return all_table_data, headlines
        else:
            return None, None