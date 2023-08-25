import openpyxl
from openpyxl import Workbook, load_workbook
from openpyxl.chart import (
    ScatterChart,
    Reference,
    Series,
)


def _translate_temp_name_to_id(all_dose_data, plate_group_to_compound_id):
    temp_samples_dict = {}

    for plates in all_dose_data:
        temp_samples_dict[plates] = []
        for samples in all_dose_data[plates]:
            temp_samples_dict[plates].append(samples)

    for plate in temp_samples_dict:
        for samples in temp_samples_dict[plate]:
            new_name = plate_group_to_compound_id[samples]
            all_dose_data[plate][new_name] = all_dose_data[plate].pop(samples)


def _write_dose_readings(ws_readings, ws_reading_row, ws_reading_clm, temp_data, add_dose, add_state_data, replicates):

    if replicates:
        temp_sample_counter = 1
        temp_rep_counter = 0
        data_placement = {"dose": {"clm": [], "row": []},
                          "sheet": ws_readings}
    else:
        data_placement = {"dose": {"clm": [], "row": []},
                          "samples": {"clm": [], "row": []},
                          "sheet": ws_readings}

    if add_dose:
        temp_dose_clm = ws_reading_clm
        temp_dose_row = ws_reading_clm
        ws_reading_clm += 1

    if add_state_data:
        for states in temp_data["state_data"]:
            ws_readings.cell(column=ws_reading_clm, row=ws_reading_row, value=states)
            for well_value in temp_data["state_data"][states]:
                ws_reading_row += 1
                ws_readings.cell(column=ws_reading_clm, row=ws_reading_row, value=well_value)
            ws_reading_clm += 1
            ws_reading_row = 1

    for sample_index, samples in enumerate(temp_data):
        if samples != "state_data":
            temp_sample_name = samples.split("_")
            temp_sample_name = f"{temp_sample_name[0]}_{temp_sample_name[1]}_{temp_sample_name[2]}_{temp_sample_name[3]}"
            ws_readings.cell(column=ws_reading_clm, row=ws_reading_row, value=samples)
            ws_reading_row += 1
            for reading_data in temp_data[samples]["reading"]["raw"]:

                ws_readings.cell(column=ws_reading_clm, row=ws_reading_row, value=reading_data)

                if replicates:

                    data_placement_name = f"samples_{temp_sample_name}"
                    try:
                        data_placement[data_placement_name]
                    except KeyError:
                        data_placement[data_placement_name] = {"clm": [], "row": []}

                else:
                    data_placement_name = "samples"

                data_placement[data_placement_name]["row"].append(ws_reading_row)
                data_placement[data_placement_name]["clm"].append(ws_reading_clm)

                ws_reading_row += 1

            if add_dose:
                temp_name = f"Dose({temp_data[samples]['dose']['unit']})"
                ws_readings.cell(column=temp_dose_clm, row=temp_dose_row, value=temp_name)
                temp_dose_row += 1

                for reading_data in temp_data[samples]["dose"]["raw"]:
                    ws_readings.cell(column=temp_dose_clm, row=temp_dose_row, value=reading_data)

                    data_placement["dose"]["clm"].append(1)
                    data_placement["dose"]["row"].append(temp_dose_row)
                    temp_dose_row += 1
                add_dose = False
            if replicates:
                temp_rep_counter += 1
                if temp_rep_counter == replicates:
                    temp_sample_counter += 1
                    temp_rep_counter = 0
            ws_reading_clm += 1
            ws_reading_row = 1
    return data_placement
    # reading - raw
    # reading - normalized
    # reading - min
    # reading - max
    # reading - fitted_normalized
    # reading - fitted
    # dose - raw
    # dose - normalized
    # dose - min
    # dose - max
    # dose - fitted_normalized
    # dose - fitted


def _write_dose_data(ws_data, ws_data_row, ws_data_clm, temp_data):

    # write_headlines
    check_list = ["EC50", "rsquared", "hillslope", "n_lowdose_datapoints", "std_resp_lowdose_datapoints",
                  "n_highdose_datapoints", "std_resp_highdose_datapoints", "doseconc_stepsize_at_EC50", "saxe_lowdose",
                  "saxe_highdose"]
    for sample_index, samples in enumerate(temp_data):
        if samples != "state_data" and sample_index == 0:
            for data in temp_data[samples]:
                if data in check_list:
                    ws_data.cell(column=ws_data_clm, row=ws_data_row + 1, value=data)
                    ws_data_row += 1

    ws_data_row = 1
    ws_data_clm += 1

    for sample_index, samples in enumerate(temp_data):

        if samples != "state_data":
            ws_data.cell(column=ws_data_clm, row=ws_data_row, value=samples)
            ws_data_row += 1
            for data in temp_data[samples]:
                if data in check_list:
                    temp_value = temp_data[samples][data]["value"]
                    ws_data.cell(column=ws_data_clm, row=ws_data_row, value=temp_value)
                    ws_data_row += 1
        ws_data_row = 1
        ws_data_clm += 1


def _draw_curves(ws_diagram, ws_diagram_row, ws_diagram_clm, temp_data, data_placement):

    x_values_min_clm = min(data_placement["dose"]["clm"])
    x_values_min_row = min(data_placement["dose"]["row"])
    x_values_max_row = max(data_placement["dose"]["row"])
    chart_data_sheet = data_placement["sheet"]
    x_values = Reference(chart_data_sheet, min_col=x_values_min_clm, min_row=x_values_min_row, max_row=x_values_max_row)
    chart_list = ["B2", "L2", "V2", "B17", "L17", "V17", "B32", "L32", "V32"]
    chart_counter = 0

    for samples in data_placement:
        if samples != "dose" and samples != "sheet":
            chart = ScatterChart()
            chart.title = f"{samples}"
            chart.style = 15
            chart.x_axis.title = 'Dose'
            chart.y_axis.title = 'Signal'
            y_value_min_clm = min(data_placement[samples]["clm"])
            y_value_max_clm = max(data_placement[samples]["clm"])
            Y_values_min_row = min(data_placement[samples]["row"])
            Y_values_max_row = max(data_placement[samples]["row"])

            for clm in range(y_value_min_clm, y_value_max_clm + 1):
                values = Reference(chart_data_sheet, min_col=clm, min_row=Y_values_min_row - 1, max_row=Y_values_max_row)
                series = Series(values, x_values, title_from_data=True)
                series.marker = openpyxl.chart.marker.Marker('x')
                series.graphicalProperties.line.noFill = True
                chart.series.append(series)
            ws_diagram.add_chart(chart, chart_list[chart_counter])
            chart_counter += 1


def dose_excel_controller(all_dose_data, plate_group_to_compound_id, include_id, save_location):


    if include_id:
        _translate_temp_name_to_id(all_dose_data, plate_group_to_compound_id)

    for plates in all_dose_data:
        wb = Workbook()
        ws_readings = wb.create_sheet("Readings")
        # Added to make sure that dose is only added once
        add_dose = True
        add_state_data = True
        replicates = 3

        ws_data = wb.create_sheet("Data")
        ws_diagram = wb.create_sheet("Diagram")

        ws_reading_row = ws_reading_clm = ws_data_row = ws_data_clm = ws_diagram_row = ws_diagram_clm = 1

        temp_data = all_dose_data[plates]
        data_placement = _write_dose_readings(ws_readings, ws_reading_row, ws_reading_clm, temp_data, add_dose,
                                              add_state_data, replicates)
        _write_dose_data(ws_data, ws_data_row, ws_data_clm, temp_data)
        _draw_curves(ws_diagram, ws_diagram_row, ws_diagram_clm, temp_data, data_placement)

        name = save_location/f"{plates}_dose_response.xlsx"
        wb.save(name)


if __name__ == "__main__":
    pass