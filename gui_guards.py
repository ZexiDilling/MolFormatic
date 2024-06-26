from pathlib import Path


def guard_purity_data_wavelength(wavelength):
    wavelength_strings = ["all"]
    # Check if value is int or str
    try:
        wavelength_data = int(wavelength)
    except ValueError:
        wavelength_data = wavelength

    # Check if value are within equipment range
    if type(wavelength_data) == int:
        if not 190 < wavelength_data < 800:
            return False, "Only values between 190 and 800"
    elif type(wavelength_data) == str:
        # if string, check if the string that can be used
        if wavelength_data not in wavelength_strings:
            return False, f"Only following values are usable strings: {wavelength_strings} "

    return True, wavelength_data


def file_type_guard(file, types):
    if type(file) == str:
        for file_types in types:
            if file.endswith(file_types):
                return True
    elif isinstance(file, Path):

        for file_types in types:
            if file.suffix == file_types:
                return True
    else:
        return False

