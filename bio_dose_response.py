from extra_functions import unit_converter, number_to_unit_converter, specific_round


def calculate_dilution_series(stock, max_concentration, min_concentration, dilutions_steps,
                              dilutions_factor, min_transfer_volume, final_vol, stock_dilution, max_solvent_concentration,
                              table_data=False):
    """
    Calculates the stock's needed for diluting a dose response series based on the dilution factor, and what
    concentration the different steps will end up having.
    The dilution is depending on using an Echo-liquid handler for the transfer, so the dilutions steps needs to be
    dividable by 2.5 nL.

    :param stock: The concentration of the stock in string formate including unit type
    :type stock: str
    :param max_concentration: The max concentration for the dilution series, in string formate including unit type
    :type max_concentration: str
    :param min_concentration: The min concentration for the dilution series, in string formate including unit type
    :type min_concentration: str
    :param dilutions_steps: How many steps to get for the dilution - This is not implemented yet as it would then be
        instead of either min_concentration or instead of dilution_factor.
    :type dilutions_steps: int
    :param dilutions_factor: How much each dilution should differe. a dilution_factor of 3 will result in
        3 fold dilution for each step
    :type dilutions_factor: int
    :param min_transfer_volume: The minimum volume the echo can transfer, in string formate including unit type.
        This could be the minimum vol of any item used for making the dilutions.
    :type min_transfer_volume: str
    :param final_vol: The final volume for each well. used to calculate concentration,
        in string formate including unit type
    :type final_vol: str
    :param stock_dilution: how much each stock needs to be diluted, when the transfer goes below "echo_min"
    :type stock_dilution: int
    :param max_solvent_concentration: How high the maximum solvent conctration can be in %
    :type max_solvent_concentration: float
    :param table_data: If the data coming out is needed for a table. will change to return elements
    :type table_data: bool
    :return:
    """

    unit_used = False
    # Converts all units to having no units. for easier calculations.
    stock_concentration_value = float(unit_converter(stock, new_unit_out=unit_used,
                                                     old_unit_out=False,  as_list=True)[0])
    max_concentration_value = float(unit_converter(max_concentration, new_unit_out=unit_used,
                                                   old_unit_out=False,  as_list=True)[0])
    min_concentration_value = float(unit_converter(min_concentration, new_unit_out=unit_used,
                                                   old_unit_out=False,  as_list=True)[0])
    min_transfer_volume_value = float(unit_converter(min_transfer_volume, new_unit_out=unit_used,
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
    all_stocks = [number_to_unit_converter(stock_concentration_value, "mM", rounding=True)]
    temp_stock_concentration = stock_concentration_value
    # Add the minimum concentration at the end to include both endpoints
    for counter, count in enumerate(all_concentration):
        vol_needed_string[counter] = {"real_conc": "", "theoretical_conc": "", "vol": "", "%_dmso": "", "stock": "",
                                      "d_fold": ""}
        vol_needed_pure[counter] = {"vol": 0, "stock": 0, "new_conc": 0}
        conc_check.append(count)

        temp_vol = (count * final_volume_value) / temp_stock_concentration
        if temp_vol / min_transfer_volume_value < 1:
            temp_stock_concentration = temp_stock_concentration/stock_dilution

            all_stocks.append(number_to_unit_converter(temp_stock_concentration, "mM", rounding=True))
            temp_vol = (count * final_volume_value) / temp_stock_concentration
        new_temp_vol = specific_round(temp_vol, "2.5nL")

        temp_dmso_conc = round((new_temp_vol / final_volume_value) * 100, 4)
        new_temp_conc = (new_temp_vol * temp_stock_concentration) / final_volume_value
        vol_needed_string[counter]["theoretical_conc"] = number_to_unit_converter(new_temp_conc, "mM", rounding=True)
        while temp_dmso_conc > max_solvent_concentration:
            new_temp_vol -= min_transfer_volume_value
            temp_dmso_conc = round((new_temp_vol / final_volume_value) * 100, 4)

        new_temp_conc = (new_temp_vol * temp_stock_concentration) / final_volume_value
        dmso_conc.append(temp_dmso_conc)
        vol_needed_pure[counter] = {"vol": new_temp_vol,
                                    "stock": temp_stock_concentration,
                                    "new_conc": new_temp_conc}

        string_new_temp_vol = number_to_unit_converter(new_temp_vol, "nL", rounding=True)
        string_temp_stock_concentration = number_to_unit_converter(temp_stock_concentration, "mM", rounding=True)
        string_new_temp_conc = number_to_unit_converter(new_temp_conc, "uM", rounding=True)
        vol_needed_string[counter]["vol"] = string_new_temp_vol
        vol_needed_string[counter]["stock"] = string_temp_stock_concentration
        vol_needed_string[counter]["real_conc"] = string_new_temp_conc
        vol_needed_string[counter]["%_dmso"] = temp_dmso_conc

        if counter > 0:
            diff_1 = float(unit_converter(vol_needed_pure[counter]["new_conc"],
                                          old_unit_out=False, new_unit_out=unit_used, as_list=True)[0])
            diff_2 = float(unit_converter(vol_needed_pure[counter-1]["new_conc"],
                                          old_unit_out=False, new_unit_out=unit_used, as_list=True)[0])
        else:
            diff_1 = float(unit_converter(vol_needed_pure[counter]["new_conc"],
                                          old_unit_out=False, new_unit_out=unit_used, as_list=True)[0])
            diff_2 = float(unit_converter(all_stocks[0],
                                          old_unit_out=False, new_unit_out=unit_used, as_list=True)[0])


        try:
            diff = round((diff_1 / diff_2) * 100, 1)
        except ZeroDivisionError:
            diff = "Error"
        diff_check.append(diff)
        vol_needed_string[counter]["d_fold"] = diff

    # These are created, but never used... can be deleted or used...
    print(f"conc_check - {conc_check}")
    print(f"diff_check - {diff_check}")
    print(f"dmso_conc - {dmso_conc}")
    print(f"vol_needed_string - {vol_needed_string}")
    print(f"vol_needed_pure - {vol_needed_pure}")
    print(f"all_stocks - {all_stocks}")
    if table_data:
        overview_table_data = []
        stock_table_data = []
        for dilution_steps in vol_needed_string:
            temp_data = []
            for info in vol_needed_string[dilution_steps]:
                temp_data.append(vol_needed_string[dilution_steps][info])
            overview_table_data.append(temp_data)

        for stocks in all_stocks:
            temp_data = [stocks, stock_dilution]
            stock_table_data.append(temp_data)

        return vol_needed_pure, overview_table_data, stock_table_data
    else:
        return vol_needed_pure
