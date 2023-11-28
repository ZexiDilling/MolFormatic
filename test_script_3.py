# import pandas as pd
# import numpy as np
# import scipy.optimize as opt
# import matplotlib.pyplot as plt
# import seaborn as sns


def ll4(x,b,c,d,e):
    '''This function is basically a copy of the LL.4 function from the R drc package with
     - b: hill slope
     - c: min response
     - d: max response
     - e: EC50'''

    return (c + (d - c) / (1 + np.exp(b * (np.log(x) - np.log(e)))))

def pDose(x):
    '''This is just a helper function, to compute easily log transformed concentrations used in drug discovery'''
    return(-np.log10(1e-6*x))


def main():

    params = [{'compound':'A', 'b':1, 'c':0, 'd':100, 'e':0.4,'startDose':10, 'nDose':8, 'dilution':3},
              {'compound':'B', 'b':0.7, 'c':0, 'd':86, 'e':1.3,'startDose':30, 'nDose':8, 'dilution':3},
              {'compound':'C', 'b':2, 'c':24, 'd':152, 'e':0.02,'startDose':3, 'nDose':8, 'dilution':3}]
    paramsCompound = [item['compound'] for item in params]

    drData=[]

    for curve in params:
        # generate base curve
        curData = pd.DataFrame(data={'compound':curve['compound'],
                                     'dose':curve['startDose']/np.power(curve['dilution'],range(curve['nDose']))})
        curData['logDose'] = pDose(curData.dose)
        curData['response'] = curData.dose.apply(lambda x: ll4(x, *[curve[i] for i in ['b', 'c', 'd', 'e']]))
        # generate replicates
        repData = []
        for i in range(5):
            rep = curData
            rep.response += 0.25*np.random.normal(len(rep.response))
            repData.append(rep.copy())
        repData = pd.concat(repData)
        drData.append(repData)
    # assemble data
    drData = pd.concat(drData)
    drData.head()


    compoundData = drData.groupby(['compound'])
    fitData = []




    for name, group in compoundData:
        fitCoefs, covMatrix = opt.curve_fit(ll4, group.dose, group.response)
        resids = group.response - group.dose.apply(lambda x: ll4(x, *fitCoefs))
        curFit = dict(zip(['b', 'c', 'd', 'e'], fitCoefs))
        curFit['compound'] = name
        curFit['residuals'] = sum(resids ** 2)
        fitData.append(curFit)
    fitCompound = [item['compound'] for item in fitData]


    fitTable = pd.DataFrame(fitData).set_index('compound')
    paramTable = pd.DataFrame(params).set_index('compound')[['b','c','d','e']]
    paramTable.columns = ['ref_'+i for i in paramTable.columns]
    fitTable.join(paramTable)
    refDose = np.linspace(min(drData.logDose)*0.9,max(drData.logDose)*1.1,256)
    refDose = (10**-refDose)*1e6
    sns.lmplot('logDose','response',data=drData,hue='compound',fit_reg=False)
    for fit in fitData:
        plt.plot([pDose(i) for i in refDose],[ll4(i,*[fit[i] for i in ['b','c','d','e']]) for i in refDose])


def testing():
    if event == "-ASSAY_RUN_APPLY_SELECTED-" or event == "-ASSAY_RUN_APPLY_ALL-":
        plate_table_index = values["-ASSAY_RUN_USED_PLATES_TABLE-"]
        run = values["-ASSAY_RUN_NAME-"]
        temp_plate_data = []
        if event == "-ASSAY_RUN_APPLY_ALL-":
            for plates in plates_table_data:
                temp_plate_data.append(plates[0])
        else:
            for index in plate_table_index:
                temp_plate_data.append(plates_table_data[index][0])
                plates_table_data[index][1] = run

        window["-ASSAY_RUN_USED_PLATES_TABLE-"].update(values=plates_table_data)

        temp_event = event
        if not values["-ASSAY_RUN_NAME-"]:
            sg.popup_error("Please select a name for the run")
        elif not values["-ASSAY_RUN_DATE_TARGET-"]:
            sg.popup_error("Please select a date for the run")
        elif not values["-ASSAY_RUN_WORKLIST_DATA-"]:
            sg.popup_error("Please provide the worklist used for the run")
        # ToDo have a check to see if the echo data fit with the plate
        else:
            assay_run_name = values["-ASSAY_RUN_NAME-"]
            latest_assay_run_name = assay_run_name
            if assay_run_name not in different_assay_runs:
                different_assay_runs.append(assay_run_name)
            # Test if the run name is in the database already, if it isn't add the name to the run data
            if assay_run_name not in previous_runs or len(previous_runs) == 1:
                assay_run_data = {"run_name": assay_run_name,
                                  "assay_name": assay_name,
                                  "batch": values["-ASSAY_RUN_ALL_BATCHES-"],
                                  "worklist": values["-ASSAY_RUN_WORKLIST_DATA-"],
                                  "echo_data": values["-ASSAY_RUN_ECHO_DATA-"],
                                  "date": values["-ASSAY_RUN_DATE_TARGET-"],
                                  "note": values["-ASSAY_RUN_NOTES-"]}

                dbf.add_records_controller("assay_runs", assay_run_data)
                # add the plates to assay_plates
                for plates in temp_plate_data:
                    assay_plate_data = {"plate_name": plates, "assay_run": assay_run_name}
                    dbf.add_records_controller("assay_plates", assay_plate_data)

                # Add the new run to the list of old runs,
                # and updates the dropdown to include it, and the value to it
                if len(previous_runs) != 1:
                    previous_runs.append(assay_run_name)
                    window["-ASSAY_RUN_PREVIOUS-"].update(values=previous_runs, value=assay_run_name)

            if temp_event == "-ASSAY_RUN_APPLY_SELECTED-":
                temp_plate_list = []
                for index in values["-ASSAY_RUN_USED_PLATES_TABLE-"]:
                    temp_plate_list.append(plates_table_data[index][0])
                    plates_table_data[index][1] = run

                # Update the table with the assay run next to the plate name
                window["-ASSAY_RUN_USED_PLATES_TABLE-"].update(values=plates_table_data)
            else:
                temp_plate_list = plates_table_data
                for index, _ in enumerate(plates_table_data):
                    plates_table_data[index][1] = run

                # Update the table with the assay run next to the plate name
                window["-ASSAY_RUN_USED_PLATES_TABLE-"].update(values=plates_table_data)

            for rows in temp_plate_list:
                plate = rows[0]
                all_plates_data[plate]["run_name"] = assay_run_name  # Update the dict with the run_name
                if plate not in plates_checked:  # Make sure that already checked plates are not added
                    plates_checked.append(plate)

            # Reset all the info, and count up one for the name:
            assay_run_name = increment_text_string(assay_run_name)
            window["-ASSAY_RUN_NAME-"].update(value=assay_run_name)
            window["-ASSAY_RUN_DATE_TARGET-"].update(value="Choose date")
            window["-ASSAY_RUN_WORKLIST_INDICATOR-"].update(value="No Worklist")
            window["-ASSAY_RUN_ECHO_INDICATOR-"].update(value="No Echo Data")
            window["-ASSAY_RUN_NOTES-"].update(value="")
            window["-ASSAY_RUN_WORKLIST_DATA-"].update(value="")
            window["-ASSAY_RUN_ECHO_DATA-"].update(value="")


if __name__ == "__main__":
    # main()
    min = 1
    max_amouint = 100
    current = 1.4

    if min < current < max_amouint:
        print("HEY")


