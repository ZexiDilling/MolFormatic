import numpy as np
import scipy
import pandas as pd
import scipy.stats


class Integration:
    def __str__(self):
        """
        A modul for integrating UV data in tensor form. to find peaks area and height, and based on that, find the
        purity of a sample/compound

        :return: A Datafrom with uv-peak information.
        """
    @staticmethod
    def _integration_uv_wavelength(sample, wavelength, slope_threshold, uv_threshold, data, rt_solvent_peak):
        """
        Calculate the integrals of the UV peaks using the trapezoid integration rule.
        Sum over one specific wavelength and calculate the UV peak integrals based on this.
        slope calculation notes:
            Find out when a peak start/stops based on the slope.
            Positive slope value indicate the begging of a peak and negative slope value indicate the end of a peak.
            Uses scipy.stats.linregress to calculate slope.
            Linregress is utilised to compute the linear least-squares regression for two given one-dimensional
            arrays of the same length.
            The slope thresholds can be changed, see '#OBS: Can change slope threshold here'.

        :param sample: sample/compound ID
        :type sample: str
        :param wavelength: what wavelength the integration is working at. can be all (limited to the data)
        :type wavelength: float
        :param slope_threshold: The minimum increase in uv intensity before a slope is consideret a peak
        :type slope_threshold: int
        :param uv_threshold: minimum threshold before data is being looked at. everything below is ignored
        :type uv_threshold: int
        :param data: All data for all the samples
        :type data: dict
        :param rt_solvent_peak: retentions time for the solvent peak, to avoid data from that peak.
        :type rt_solvent_peak: float
        :return: df_integrate_overview is a pandas dataframe with peak list, integrals, peak start time and peak end
            time.
        :rtype: pandas.core.frame.DataFrame
        """

        uv_tensor = np.array(data[sample]["uv"])
        points = 5
        # Save both slope values and the corresponding retentions times
        slope = []
        st_slope = []
        start_slope = []
        end_slope = []

        # Define UV retention times in seconds and wavelengths
        uv_rt = data[sample]["uv"].index*60
        if wavelength == "all":
            uv_tensor = np.sum(uv_tensor, axis=1)
            wave_indx = None

        else:
            wavelength = int(wavelength)
            uv_wavelengths = data[sample]["uv"].columns
            # Find index of selected wavelength
            wave_number = uv_wavelengths[(np.fabs(uv_wavelengths-wavelength)).argmin(axis=0)]
            wave_indx = [i for i, x in enumerate(uv_wavelengths == wave_number) if x][0]


        # Calculate the slope based on 5 data points
        for rt in np.arange(0, len(uv_tensor[points:, wave_indx])):
            x = uv_rt.tolist()[(0+rt):(points+rt)]
            if wave_indx:
                y = uv_tensor[(0 + rt):(points + rt), wave_indx]
            else:
                y = uv_tensor[(0 + rt):(points + rt)]

            if y[0] > uv_threshold:
                temp_slope = scipy.stats.linregress(x, y)[0]
                slope.append(temp_slope)
                # st_slope.append(uv_rt[5+rt])
                index_value = int(round(points/2))
                start_slope.append(uv_rt[index_value + rt])
                end_slope.append(uv_rt[index_value + rt])



        # Find the start and end times based on the slope
        rt_start_peak = []
        rt_end_peak = []

        x = -1  # x makes sure the RT's are order correctly
        for s in np.arange(1, len(slope)-1):
            if slope[s] > slope_threshold > slope[s-1] and x < 0:

                # rt_start_peak.append(st_slope[s])
                rt_start_peak.append(start_slope[s])
                x = 1
            if slope[s] < -slope_threshold < slope[s+1] and x > 0:
                # rt_end_peak.append(st_slope[s])
                rt_end_peak.append(end_slope[s])
                x = -1

        # For the last peak it is possible only to detect rt_start_peak
        # If this happens, delete the last observation in rt_start_peak
        if len(rt_start_peak) > len(rt_end_peak):
            del(rt_start_peak[-1])

        rt_start_peak_new = []
        rt_end_peak_new = []

        # Deleting the solvent peak
        for x in np.arange(0, len(rt_start_peak)):
            if rt_start_peak[x]/60 > rt_solvent_peak:
                rt_start_peak_new.append(rt_start_peak[x])
                rt_end_peak_new.append(rt_end_peak[x])

        # Calculate the integrals numerically with trapezoid integration rule
        integrals = []
        for s in np.arange(0, len(rt_start_peak_new)):

            # Find index of start and end times
            indx_start = [i for i, x in enumerate(uv_rt == rt_start_peak_new[s]) if x][0]
            indx_end = [i for i, x in enumerate(uv_rt == rt_end_peak_new[s]) if x][0]

            # Calculate n, a, b and h
            n = len(uv_rt.tolist()[indx_start:indx_end+1])

            a = uv_rt.tolist()[indx_start]
            b = uv_rt.tolist()[indx_end]
            h = (b - a) / n

            # Define x and y values based on the peak time interval

            if wave_indx:
                y = uv_tensor[indx_start:indx_end+1, wave_indx]
            else:
                y = uv_tensor[indx_start:indx_end + 1]

            # Calculate the integrals using the trapezoid integration rule
            int_list = []
            for i in range(n):
                if i == 0 or i == n - 1:
                    int_list.append(h / 2 * y[i])
                else:
                    int_list.append(h * y[i])
            integral = sum(int_list)
            integrals.append(round(integral))

        # Save peak information and collect it in a dataframe
        peak_list = []
        for rt in np.arange(0, len(rt_start_peak_new)):
            rt_start_peak_new[rt] = round(rt_start_peak_new[rt]/60, 3)
            rt_end_peak_new[rt] = round(rt_end_peak_new[rt]/60, 3)
            peak_list.append(f"Peak {rt + 1}")

        integrals_total = sum(integrals)
        purity = []
        for x in integrals:
            pur = round((x/integrals_total) * 100, 3)
            purity.append(pur)

        # Save in a dataframe
        df = {"Peak list": peak_list, "Integrals": integrals, "Peak start time": rt_start_peak_new,
              "Peak end time": rt_end_peak_new, "purity": purity}

        df_integrate_overview = pd.DataFrame(df)
        if df_integrate_overview.empty:
            df = {"Peak list": ["No peaks detected"], "Integrals": ["No peaks detected"], "Peak start time":
                  ["No peaks detected"], "Peak end time": ["No peaks detected"], "purity": ["No peaks detected"]}
            df_integrate_overview = pd.DataFrame(df)

        return df_integrate_overview

    def calculate_uv_integrals(self, data, slope_threshold, uv_threshold, solvent_peak, sample_data, wavelength_data,
                               sample):
        """
        Calculate integrals for all samples in a plate and combine them to a list of pandas dataframes.
        One dataframe for each sample.
        A dataframe contains the following information, peak list, integrals, peak start time and peak and time.

        :param data: All data for all the samples
        :type data: dict
        :param slope_threshold: The minimum increase in uv intensity before a slope is consideret a peak
        :type slope_threshold: int
        :param uv_threshold: minimum threshold before data is being looked at. everything below is ignored
        :type uv_threshold: int
        :param solvent_peak: retentions time for the solvent peak, to avoid data from that peak.
        :type solvent_peak: float
        :param sample_data: A dict of data for each sample
        :type sample_data: dict or None
        :param wavelength_data: a str or int for all wavelength that are not in sample_data
        :type wavelength_data: str or int
        :param sample: A sample
        :type sample: str or None
        :return: df_combined is a list of pandas dataframes.
        :rtype: dict
        """
        peak_information = {}

        # For now only works for all samples. Needs to make a new dict over samples and batches, to run this over,
        # to make it more flexible. It needs to be added to all methods using the module.
        sample_list = []
        if sample:
            sample_list.append(sample)
        else:
            for samples in data:
                sample_list.append(samples)

        for samples in sample_list:
            peak_information[samples] = {}
            try:
                wavelength = sample_data[samples]["wavelength"]
            except TypeError or KeyError:
                wavelength = wavelength_data
            temp_data_frame = self._integration_uv_wavelength(samples, wavelength, slope_threshold, uv_threshold,
                                                              data, solvent_peak)

            peak_information[samples] = temp_data_frame

        return peak_information
