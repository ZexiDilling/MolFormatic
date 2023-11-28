from PySimpleGUI import PopupError

from bio_dose_response import calculate_dilution_series


def calculate_dose(window, values):
    if not values["-CALCULATIONS_INFO_DOSE_STOCK-"]:
        PopupError("Please Fill in all the informations")
    elif not values["-CALCULATIONS_INFO_DOSE_STOCK_DILUTION-"]:
        PopupError("Please Fill in all the informations")
    elif not values["-CALCULATIONS_INFO_DOSE_MAX_CONC-"]:
        PopupError("Please Fill in all the informations")
    elif not values["-CALCULATIONS_INFO_DOSE_MIN_CONC-"]:
        PopupError("Please Fill in all the informations")
    elif not values["-CALCULATIONS_INFO_DOSE_FINAL_VOL-"]:
        PopupError("Please Fill in all the informations")
    elif not values["-CALCULATIONS_INFO_DOSE_DILUTION_FACTOR-"]:
        PopupError("Please Fill in all the informations")
    elif not values["-CALCULATIONS_INFO_MAX_SOLVENT_CONCENTRATION-"]:
        PopupError("Please Fill in all the informations")
    else:

        stock = f'{values["-CALCULATIONS_INFO_DOSE_STOCK-"]}{values["-CALCULATIONS_INFO_DOSE_STOCK_UNIT-"]}'
        max_concentration = f'{values["-CALCULATIONS_INFO_DOSE_MAX_CONC-"]}{values["-CALCULATIONS_INFO_DOSE_MAX_CONC_UNIT-"]}'
        min_concentration = f'{values["-CALCULATIONS_INFO_DOSE_MIN_CONC-"]}{values["-CALCULATIONS_INFO_DOSE_MIN_CONC_UNIT-"]}'
        dilutions_steps = values["-CALCULATIONS_INFO_DOSE_DILUTION_STEPS-"]
        dilutions_factor = int(values["-CALCULATIONS_INFO_DOSE_DILUTION_FACTOR-"])
        echo_min = f'{values["-CALCULATIONS_INFO_DOSE_MIN_TRANS_VOL-"]}{values["-CALCULATIONS_INFO_DOSE_MIN_TRANS_VOL_UNIT-"]}'
        final_vol = f'{values["-CALCULATIONS_INFO_DOSE_FINAL_VOL-"]}{values["-CALCULATIONS_INFO_DOSE_FINAL_VOL_UNIT-"]}'
        stock_dilution = int(values["-CALCULATIONS_INFO_DOSE_STOCK_DILUTION-"])
        max_solvent_concentration = float(values["-CALCULATIONS_INFO_MAX_SOLVENT_CONCENTRATION-"])

        vol_needed_pure, overview_table_data, stock_table_data = \
            calculate_dilution_series(stock, max_concentration, min_concentration, dilutions_steps,
                                      dilutions_factor, echo_min, final_vol, stock_dilution,
                                      max_solvent_concentration, table_data=True)

        window["-CALCULATIONS_TABLE_OVERVIEW-"].update(values=overview_table_data)
        window["-CALCULATIONS_TABLE_STOCK-"].update(values=stock_table_data)