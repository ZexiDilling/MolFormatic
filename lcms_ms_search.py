import numpy as np
import configparser


def _aduct_search(ms_mode):
    """
    Finds aduct from the config file, and add them to a list, that is used when going over the masses found.

    :param ms_mode: what mode the data is in. or is being looked at.
    :type ms_mode: str
    :return: A list of all aduct for either positive or negative ion-mode
    """
    config = configparser.ConfigParser()
    ms_aduct = []
    ions = []

    if ms_mode == 'ms_pos':
        config.read("config.ini")
        for sections in config.sections():
            if sections == "Positive ion mode":
                for ion in config[sections]:
                    ions.append(ion)
                    temp = config[sections][ion].split(",")
                    try:
                        ms_aduct.append([float(temp[0]), float(temp[1])])
                    except ValueError:
                        temp_2 = temp[1].split("/")
                        temp_2 = int(temp_2[0]) / int(temp_2[1])
                        ms_aduct.append([float(temp[0]), float(temp_2)])

    if ms_mode == 'ms_neg':
        config.read("config.ini")
        for sections in config.sections():
            if sections == "Negative ion mode":
                for ion in config[sections]:
                    ions.append(ion)
                    temp = config[sections][ion].split(",")
                    try:
                        ms_aduct.append([float(temp[0]), float(temp[1])])
                    except ValueError:
                        temp_2 = temp[1].split("/")
                        temp_2 = int(temp_2[0]) / int(temp_2[1])
                        ms_aduct.append([float(temp[0]), float(temp_2)])

    return ms_aduct, ions


def mass_search(data, peak_information, ms_mode, sample_data, delta_mass, mz_threshold, peak_amounts, compound_mass):
    """
    M/z search function by top n m/z-values in a peak.
    A m/z value is significant if it is between the top n most intense peaks.

    :param data: All the data for all the compounds
    :type data: dict
    :param peak_information: Information from the uv date.
    :type peak_information: dict
    :param ms_mode: What mode the data is in. positive or negative
    :type ms_mode: str
    :param sample_data: Data for the sample, from a file or compound database.
    :type sample_data: dict or None
    :param compound_mass: Mass for the sample/compound
    :type compound_mass: float or None
    :param delta_mass: +/- range for the mass search field
    :type delta_mass: float
    :param mz_threshold: The minimum amount of signal before the data is being recognised as a peak.
    :type mz_threshold: float
    :param peak_amounts: Amount of peaks to include in the search for the mass
    :type peak_amounts: int

    :return: A dict of peaks for ms-data
    :rtype: dict
    """


    peak_hit = {}
    for sample in data:

        if "blank" in sample.casefold():
            continue

        try:
            sample_data[sample]["mass"]
        except KeyError:
            continue

        try:
            mass = sample_data[sample]["mass"]
        except ValueError:
            mass = compound_mass
        # Guard for mass written with ',' instead of '.'
        try:
            mass = float(mass)
        except ValueError:
            try:
                mass = float(mass.replace(",", "."))
            except:
                print(sample)
                print(sample_data[sample])

        ms_mz = data[sample][ms_mode].columns
        ms_rt = data[sample][ms_mode].index
        peak_hit[sample] = {}
        peak_table = peak_information[sample]
        peak_table_t = np.transpose(peak_table)
        ms_tensor = np.array(data[sample][ms_mode])

        # for index, _ in enumerate(ms_tensor):
        #     ms_data_sample = ms_tensor[index]

        # Loop over all peaks in sample s
        for i in peak_table_t:
            # Find UV RT's (converted to seconds)
            uv_peak_start_rt = peak_table_t[i][2]*60
            uv_peak_end_rt = peak_table_t[i][3]*60
            if type(uv_peak_end_rt) == str:
                continue
            # MS peaks appears later than the corresponding UV peaks in the spectrum
            # Add 12 seconds (0.2 min) to UV_peak_end_RT to correct for the above
            uv_peak_end_rt = uv_peak_end_rt+12

            ms_peak_start_rt = ms_rt[(np.fabs(ms_rt - uv_peak_start_rt)).argmin(axis=0)]
            ms_peak_end_rt = ms_rt[(np.fabs(ms_rt - uv_peak_end_rt)).argmin(axis=0)]

            ms_peak_start_rt_index = [i for i, x in enumerate(ms_rt == ms_peak_start_rt) if x][0]
            ms_peak_end_rt_index = [i for i, x in enumerate(ms_rt == ms_peak_end_rt) if x][0]

            x_temp = ms_tensor[ms_peak_start_rt_index:ms_peak_end_rt_index+1]

            x_sum = np.sum(x_temp, axis=0)

            idx_thress = []
            for index, _ in enumerate(x_sum):
                if x_sum[index] > mz_threshold:
                    idx_thress.append(index)

            idx_amount = np.argpartition(x_sum, -peak_amounts)[-peak_amounts:]

            idx = []
            for values in idx_amount:
                if values in idx_thress:
                    idx.append(values)

            mz_max_values = np.array(ms_mz[idx])
            mz_aduct, ions = _aduct_search(ms_mode)

            for index, _ in enumerate(idx):
                for counter, aduct in enumerate(mz_aduct):
                    temp_mass = (mass * aduct[1]) + aduct[0]
                    if temp_mass-delta_mass < mz_max_values[index] < temp_mass+delta_mass:
                        peak_hit[sample][peak_table_t[i][0]] = {}
                        peak_hit[sample][peak_table_t[i][0]][ions[counter]] = mz_max_values[index]

    return peak_hit


if __name__ == "__main__":

    test = _aduct_search("pos")
    print(test)

