import numpy as np
import pandas as pd
from statistics import mean
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
# import matplotlib.transforms as mtransforms

ms_title = {"ms_neg": "MS-", "ms_pos": "MS+"}


def uv_chromatogram(data, canvas, samples, fig_size, canvas_lines):
    """uv chromatogram of different sample(s) in the same plate.
    Summed over all wavelength and plots intensity as a function of retention time.
    Samples is a vector with the sample number(s) of intrest.
    P is the plate number.
    uv_tensor is the tensor of intrest, e.g. P1_uv_tensor. When calling the function it is only necessary to define plate number (P).
    Data is the main dictionary containing all data."""
    # Define uv retention times
    uv_tensor = []
    for index_samples, sample in enumerate(samples):
        if index_samples == 0:

            uv_rt = data[sample]["uv"].index
        uv_tensor.append(data[sample]["uv"])

    uv_tensor = np.array(uv_tensor)

    # Sum over all wavelengths
    uv_tensor_sum_wavelengths = np.sum(uv_tensor, axis=2)
    # Plot figure

    fig, ax = plt.subplots(figsize=fig_size)
    for index, sample in enumerate(samples):
        plt.plot(uv_tensor_sum_wavelengths[index, :], label=sample)

    if canvas_lines["uv"]:
        plt.hlines(float(canvas_lines["uv"]), 0, len(uv_rt), color="orange", linestyle="--")

    if canvas_lines["peak_lines"]:

        if data[sample]["scan_speed"][1] == "msec":
            scan_speed = 60000 / int(data[sample]["scan_speed"][0])
        else:
            print(f"No Setting for {data[sample]['scan_speed'][1]}")

        for counter, peaks in enumerate(canvas_lines["peak_lines"]):
            if counter % 2 == 0:
                colour = "red"
            else:
                colour = "green"
            ax.axvspan(canvas_lines["peak_lines"][peaks]["start"] * scan_speed,
                       canvas_lines["peak_lines"][peaks]["end"] * scan_speed,
                       alpha=0.2, color=colour)

    plt.legend(loc='best', ncol=2, shadow=True, fancybox=True)
    ax.set_xticks(list(range(0, len(uv_rt), 75)))
    ax.set_xticklabels([str(x) for x in uv_rt[::75]])
    ax.set_xlabel(r'Retention time [min]')
    ax.set_ylabel(r'Intensity [a.u.]')
    ax.set_xticks
    plt.subplots_adjust(bottom=0.25)
    figure_canvas_agg = FigureCanvasTkAgg(fig, canvas.TKCanvas)

    return figure_canvas_agg


def ms_chromatogram(data, canvas, samples, fig_size, ms_mode):
    """MS(+ or -) chromatogram for one or several sample(s).
    Plots intensity as a function of retention time.
    Samples is a vector with the sample number(s) of intrest.
    P is the plate number.
    MS_mode can only take 'positive' or 'negative' as input.
    MS_tensor is the tensor of intrest, e.g. P4_MS_pos_tensor.
    data is the main dictionary containing all data."""
    ms_tensor = []
    for index_samples, sample in enumerate(samples):
        if index_samples == 0:

            ms_rt = data[sample][ms_mode].index
        ms_tensor.append(data[sample][ms_mode])

    ms_tensor = np.array(ms_tensor)
    # Sum over all mz-values

    ms_tensor_sum_mz = np.sum(ms_tensor, axis=2)
    # Plot figure
    fig, ax = plt.subplots(figsize=fig_size)

    for value, sample in enumerate([value - 1 for value, sample in enumerate(samples)]):
        plt.plot(ms_tensor_sum_mz[value, :], label=samples[value])
    plt.legend(loc='best', ncol=2, shadow=True, fancybox=True)
    ax.set_xticks(list(range(0, len(ms_rt), 5)))
    ax.set_xticklabels([str(x) for x in ms_rt[::5]])
    ax.set_xlabel(r'Retention time [s]')
    ax.set_ylabel(r'Intensity [a.u.]')
    plt.xticks(rotation=90)
    plt.subplots_adjust(bottom=0.25)
    plt.title(f"{ms_title[ms_mode]} Chromatogram")

    figure_canvas_agg = FigureCanvasTkAgg(fig, canvas.TKCanvas)

    return figure_canvas_agg


def ms_spectrum(data, canvas, samples, fig_size, ms_mode, rt_mz):
    """MS(+ or -) spectrum of one or several sample(s) at one retention time.
    Plots intensity as a function of m/z-values.
    samples is a vector with the sample number(s) of intrest.
    P is the plate number.
    RT_mz is the retention time of interest.
    MS_mode can only take 'positive' or 'negative' as input.
    MS_tensor is the tensor of intrest, e.g. P4_MS_pos_tensor.
    data is the main dictionary containing all data."""
    #Convert to seconds
    rt_mz = rt_mz*60

    ms_tensor = []
    for index_samples, sample in enumerate(samples):
        if index_samples == 0:
            ms_mz = data[sample][ms_mode].columns
            ms_rt = data[sample][ms_mode].index
        ms_tensor.append(data[sample][ms_mode])

    ms_tensor = np.array(ms_tensor)

    # Find index for selected retention time
    rt_closest_value = ms_rt[(np.fabs(ms_rt-rt_mz)).argmin(axis=0)]
    rt_indx = [i for i, x in enumerate(ms_rt == rt_closest_value) if x][0]
    # Plot figure
    fig, ax = plt.subplots(figsize=fig_size)

    # for sample in [sample - 1 for sample in samples]:
    for value, sample in enumerate([value - 1 for value, sample in enumerate(samples)]):
        print(value)
        print(rt_indx)
        plt.plot(ms_tensor[value, rt_indx, :], label=samples[value])
    plt.legend(loc='best', ncol=2, shadow=True, fancybox=True)
    ax.set_xticks(list(range(0, len(ms_mz), 500)))
    ax.set_xticklabels([str(x) for x in ms_mz[::500]])
    ax.set_xlabel(r'm/z')
    ax.set_ylabel(r'Intensity [a.u.]')
    plt.xticks(rotation=90)
    plt.subplots_adjust(bottom=0.25)
    plt.title(f"{ms_title[ms_mode]} plot at {round(ms_rt[rt_indx]/60, 1)} minutes")

    figure_canvas_agg = FigureCanvasTkAgg(fig, canvas.TKCanvas)

    return figure_canvas_agg


def ms_spectrum_range(data, canvas, samples, fig_size, ms_mode, canvas_lines):
    if not canvas_lines["peak_lines"]:
        return "Missing data"

    for peaks in canvas_lines["peak_lines"]:
        rt_ms_start = canvas_lines["peak_lines"][peaks]["start"] * 60
        rt_ms_end = canvas_lines["peak_lines"][peaks]["end"] * 60

    ms_tensor = []
    for index, sample in enumerate(samples):
        if index == 0:
            ms_mz = data[sample][ms_mode].columns
            ms_rt = data[sample][ms_mode].index
        ms_tensor.append(data[sample][ms_mode])

    list_rt_diff = []
    temp_rt = None
    for rt in ms_rt:
        if not temp_rt:
            temp_rt = rt
        else:
            list_rt_diff.append(rt-temp_rt)
            temp_rt = rt
    avg_diff_rt = int(round(mean(list_rt_diff)))
    ms_tensor = np.array(ms_tensor)
    rt_close_values = []
    x_counter = []
    rt_ms_start = int(round(rt_ms_start))
    rt_ms_end = int(round(rt_ms_end))
    for index_rts, rts in enumerate(range(rt_ms_start, rt_ms_end, avg_diff_rt)):
        rt_close_values.append(ms_rt[(np.fabs(ms_rt-rts)).argmin(axis=0)])
        x_counter.append(index_rts)
    rt_index_list = []
    for rts_close in rt_close_values:
        rt_index_list.append([i for i, x in enumerate(ms_rt == rts_close) if x][0])

    fig, ax = plt.subplots(figsize=fig_size)
    for index_r, rts in enumerate(rt_index_list):
        if index_r == 0:
            plt.plot(ms_tensor[0, rts, :])
    ax.set_xticks(list(range(0, len(ms_mz), 100)))
    ax.set_xticklabels([str(x) for x in ms_mz[::100]])

    plt.xticks(fontsize=8)
    plt.xticks(rotation="vertical")
    ax.set_xlabel(r'm/z')
    ax.set_ylabel(r'Intensity [a.u.]')
    plt.subplots_adjust(bottom=0.25)
    figure_canvas_agg = FigureCanvasTkAgg(fig, canvas.TKCanvas)

    return figure_canvas_agg


def heatmap_uv_sample(data, canvas, samples, fig_size):
    """uv heatmap for one specific sample.
    sample is the sample number of interest.
    P is the plate number.
    uv_tensor is the tensor of intrest, e.g. P1_uv_tensor.
    data is the main dictionary containing all data."""

    # Define uv wavelengtjs and retention times

    uv_tensor = []
    for index_samples, sample in enumerate(samples):
        if index_samples == 0:
            uv_wavelengths = np.round(data[sample]["uv"].columns, 1)
            uv_rt = data[sample]["uv"].index
        uv_tensor.append(data[sample]["uv"])

    uv_tensor = np.array(uv_tensor)

    # Plot figure
    fig, ax = plt.subplots(figsize=fig_size)
    for value, sample in enumerate(samples):
        sns.heatmap(uv_tensor[value-1, :, :], ax=ax, center=0)
    ax.set_xticks(list(range(0, len(uv_wavelengths), 25)))
    ax.set_xticklabels([str(x) for x in uv_wavelengths[::25]])
    ax.set_yticks(list(range(0, len(uv_rt), 50)))
    ax.set_yticklabels([str(x) for x in uv_rt[::50]])
    ax.set_xlabel(r'Wavelength [nm]')
    ax.set_ylabel(r'Retention time [min]')
    plt.title(f"UV/Vis heat map for {samples}")
    plt.subplots_adjust(bottom=0.25)
    fig.tight_layout()

    figure_canvas_agg = FigureCanvasTkAgg(fig, canvas.TKCanvas)

    return figure_canvas_agg


def heatmap_uv_rt(data, canvas, samples, fig_size, rt_value):
    """uv heatmap for one specific retention time.
    RT_value is the rention time of interest.
    P is the plate number.
    samples is a vector with the sample number(s) of intrest.
    uv_tensor is the tensor of intrest, e.g. P1_uv_tensor.
    data is the main dictionary containing all data."""
    # Extract uv data for the given samples
    uv_tensor = []
    # Define uv wavelengths and retention times
    for index_samples, sample in enumerate(samples):
        if index_samples == 0:
            uv_wavelengths = np.round(data[sample]["uv"].columns, 1)
            uv_rt = data[sample]["uv"].index
        uv_tensor.append(data[sample]["uv"])

    uv_tensor = np.array(uv_tensor)
    uv_tensor_samples = []
    for value, sample in enumerate(samples):
        uv_tensor_samples.append(uv_tensor[value-1])
    uv_tensor_samples = np.array(uv_tensor_samples)

    # Find index for selected retention time
    rt_closest_value = uv_rt[(np.fabs(uv_rt-rt_value)).argmin(axis=0)]
    rt_indx = [i for i, x in enumerate(uv_rt == rt_closest_value) if x][0]

    # Plot figure
    fig, ax = plt.subplots(figsize=fig_size)
    sns.heatmap(uv_tensor_samples[:, rt_indx, :], ax=ax, center=0)
    ax.set_xticks(list(range(0, len(uv_wavelengths), 25)))
    ax.set_xticklabels([str(x) for x in uv_wavelengths[::25]])
    ax.set_yticklabels(samples)
    ax.set_xlabel(r'Wavelength [nm]')
    ax.set_ylabel(r'Sample number')
    plt.title("UV/Vis heat map at {} min".format(uv_rt[rt_indx]))
    plt.subplots_adjust(bottom=0.25)
    fig.tight_layout()

    figure_canvas_agg = FigureCanvasTkAgg(fig, canvas.TKCanvas)

    return figure_canvas_agg


def heatmap_uv_wavelength(data, canvas, samples, fig_size, wavelength):
    """uv heatmap for one specific wavelength.
    RT_value is the rention time of interest.
    P is the plate number.
    samples is a vector with the sample number(s) of intrest.
    uv_tensor is the tensor of intrest, e.g. P1_uv_tensor.
    data is the main dictionary containing all data."""
    # Extract uv data for the given samples
    uv_tensor = []
    # Define uv wavelengths and retention times
    for index_samples, sample in enumerate(samples):
        if index_samples == 0:
            uv_wavelengths = data[sample]["uv"].columns
            uv_rt = data[sample]["uv"].index
        uv_tensor.append(data[sample]["uv"])

    uv_tensor = np.array(uv_tensor)
    uv_tensor_samples = []
    for value, sample in enumerate(samples):
        uv_tensor_samples.append(uv_tensor[value-1])
    uv_tensor_samples = np.array(uv_tensor_samples)
    # Find index for selected wavelength
    wave_number = uv_wavelengths[(np.fabs(uv_wavelengths-wavelength)).argmin(axis=0)]
    wave_indx = [i for i, x in enumerate(uv_wavelengths == wave_number) if x][0]
    # Plot figure
    fig, ax = plt.subplots(figsize=fig_size)
    sns.heatmap(uv_tensor_samples[:, :, wave_indx], ax=ax, center=0)
    # sns.heatmap(uv_tensor_samples[:, :, wave_indx], ax=ax, center=0)
    sns.color_palette("vlag", as_cmap=True)
    ax.set_xticks(list(range(0, len(uv_rt), 50)))
    ax.set_xticklabels([str(x) for x in uv_rt[::50]])
    ax.set_yticklabels(samples)
    ax.set_xlabel(r"Retention time [min]")
    ax.set_ylabel(r"Sample number")
    plt.title(f"UV/Vis heat map at wavelength {uv_wavelengths[wave_indx]} nm")
    plt.subplots_adjust(bottom=0.25)
    fig.tight_layout()

    figure_canvas_agg = FigureCanvasTkAgg(fig, canvas.TKCanvas)

    return figure_canvas_agg


def heatmap_ms_sample_binned(data, canvas, samples, fig_size, ms_mode, bin_numbers):
    """MS heatmap for one specific sample with binned m/z-values.
    sample is the sample number of interest.
    P is the plate number.
    MS_mode is either 'pos' or 'neg'.
    uv_tensor is the tensor of intrest, e.g. P1_uv_tensor.
    data is the main dictionary containing all data."""
    ms_tensor = []
    for index_samples, sample in enumerate(samples):
        if index_samples == 0:
            ms_mz = data[sample][ms_mode].columns
            ms_rt = data[sample][ms_mode].index
        ms_tensor.append(data[sample][ms_mode])

    ms_tensor = np.array(ms_tensor)

    for value, sample in enumerate(samples):
        temp_value = value
    ms_rt = np.round(ms_rt / 60, 1)
    # Bin m/z-values accoridng to the number of bins defined
    ms = []
    ms_columns = []
    bin_size = round(len(ms_mz)/bin_numbers)
    for i in range(0, bin_numbers):
        if i == 0:
            ms.append(np.sum(ms_tensor[temp_value-1, :, 0:bin_size+bin_size], axis=1))
            ms_columns.append(ms_mz[round(bin_size/2)])
        else:
            ms.append(np.sum(ms_tensor[temp_value-1, :, i*bin_size:i*bin_size+bin_size], axis=1))
            ms_columns.append(ms_mz[(bin_size*i)+round(bin_size/2)])
    ms = pd.DataFrame(ms)
    ms.index = ms_columns
    # Plot figure
    fig, ax = plt.subplots(figsize=fig_size)
    sns.heatmap(ms, ax=ax, center=0)
    ax.set_xticks(list(range(0, len(ms_rt), 5)))
    ax.set_xticklabels([str(x) for x in ms_rt[::5]])
    ax.set_ylabel(r'm/z-values (binned)')
    ax.set_xlabel(r'Retention time [min]')
    plt.title(f"{ms_title[ms_mode]} heat map for sample: {samples}")
    plt.subplots_adjust(bottom=0.25)
    fig.tight_layout()

    figure_canvas_agg = FigureCanvasTkAgg(fig, canvas.TKCanvas)

    return figure_canvas_agg


def heatmap_ms_rt_binned(data, canvas, samples, fig_size, ms_mode, bin_numbers, rt_value):
    """MS heatmap for one specific retention time with binned m/z-values.
    RT_value is the rention time of interest.
    P is the plate number.
    samples is a vector with the sample number(s) of intrest.
    MS_mode is either 'pos' or 'neg'.
    uv_tensor is the tensor of intrest, e.g. P1_uv_tensor.
    data is the main dictionary containing all data."""
    ms_tensor = []
    for index_samples, sample in enumerate(samples):
        if index_samples == 0:
            ms_mz = data[sample][ms_mode].columns
            ms_rt = data[sample][ms_mode].index
        ms_tensor.append(data[sample][ms_mode])

    ms_tensor = np.array(ms_tensor)
    # Extract MS data for the given samples
    # ms_tensor_samples = []
    # for value, samples in enumerate(samples):
    #     ms_tensor_samples.append(ms_tensor[value - 1])
    # ms_tensor_samples = np.array(ms_tensor_samples)
    # ms_rt = np.round(ms_rt/60, 1)
    # Find index for selected retention time
    rt_closest_value = ms_rt[(np.fabs(ms_rt-rt_value)).argmin(axis=0)]
    rt_indx = [i for i, x in enumerate(ms_rt == rt_closest_value) if x][0]
    # Bin m/z-values accoridng to the number of bins defined
    ms = []
    ms_columns = []
    bin_size = round(len(ms_mz)/bin_numbers)
    for i in range(0, bin_numbers):
        if i == 0:
            # ms.append(np.sum(ms_tensor_samples[:, rt_indx, 0:bin_size+bin_size], axis=1))
            ms.append(np.sum(ms_tensor[:, rt_indx, 0:bin_size + bin_size], axis=1))
            ms_columns.append(ms_mz[round(bin_size/2)])
        else:
            # ms.append(np.sum(ms_tensor_samples[:, rt_indx, i*bin_size:i*bin_size+bin_size], axis=1))
            ms.append(np.sum(ms_tensor[:, rt_indx, i * bin_size:i * bin_size + bin_size], axis=1))
            ms_columns.append(ms_mz[(bin_size*i)+round(bin_size/2)])
    ms = pd.DataFrame(ms)
    ms.index = ms_columns

    # Plot figure
    fig, ax = plt.subplots(figsize=fig_size)
    sns.heatmap(ms, ax=ax, center=0)
    ax.set_xticklabels(samples)
    ax.set_xlabel(r'Sample number')
    ax.set_ylabel(r'm/z-values (binned)')
    plt.title(f"{ms_title[ms_mode]} heatmap for data at {ms_rt[rt_indx]} minutes")
    plt.subplots_adjust(bottom=0.25)
    fig.tight_layout()

    figure_canvas_agg = FigureCanvasTkAgg(fig, canvas.TKCanvas)

    return figure_canvas_agg


def heatmap_ms_mz(data, canvas, samples, fig_size, ms_mode, mz_value):
    """MS heatmap for one specific m/z-value.
    mz_value is the m/z-value of intrest.
    P is the plate number.
    samples is a vector with the sample number(s) of intrest.
    MS_mode is either 'pos' or 'neg'.
    uv_tensor is the tensor of intrest, e.g. P1_uv_tensor.
    data is the main dictionary containing all data."""
    ms_tensor = []
    for index_samples, sample in enumerate(samples):
        if index_samples == 0:
            ms_mz = data[sample][ms_mode].columns
            ms_rt = data[sample][ms_mode].index
        ms_tensor.append(data[sample][ms_mode])

    ms_tensor = np.array(ms_tensor)
    # Extract MS data for the given samples
    # ms_tensor_samples = []
    # for value, samples in enumerate(samples):
    #     ms_tensor_samples.append(ms_tensor[value-1])
    # ms_pos_tensor_samples = np.array(ms_tensor_samples)

    ms_rt = np.round(ms_rt/60, 1)
    # Find index for selected m/z-value
    mz_closest_value = ms_mz[(np.fabs(ms_mz-mz_value)).argmin(axis=0)]
    mz_indx = [i for i, x in enumerate(ms_mz == mz_closest_value) if x][0]
    # Plot figure
    fig, ax = plt.subplots(figsize=fig_size)
    # sns.heatmap(ms_pos_tensor_samples[:, :, mz_indx], ax=ax, center=0)
    sns.heatmap(ms_tensor[:, :, mz_indx], ax=ax, center=0)
    ax.set_xticks(list(range(0, len(ms_rt), 5)))
    ax.set_xticklabels([str(x) for x in ms_rt[::5]])
    ax.set_ylabel(r'Sample number')
    ax.set_xlabel(r'Retention time [min]')
    plt.title(f"{ms_title[ms_mode]} heatmap for data at {ms_mz[mz_indx]}")
    plt.subplots_adjust(bottom=0.25)
    fig.tight_layout()

    figure_canvas_agg = FigureCanvasTkAgg(fig, canvas.TKCanvas)

    return figure_canvas_agg


class Toolbar(NavigationToolbar2Tk):
    # only display the buttons we need
    toolitems = [t for t in NavigationToolbar2Tk.toolitems if
                 t[0] in ('Home', 'Pan', 'Zoom')]
                # t[0] in ('Home', 'Pan', 'Zoom','Save')]

    def __init__(self, *args, **kwargs):
        super(Toolbar, self).__init__(*args, **kwargs)

