

#%%
from NINAconverter import NINA


def se_filter(filter_ : str):
    se_filter = {
                "$id": "67",
                "$type": "NINA.Core.Model.Equipment.FilterInfo, NINA.Core",
                "_name": filter_,
                "_focusOffset": 0,
                "_position": 3,
                "_autoFocusExposureTime": -1.0,
                "_autoFocusFilter": False,
                "FlatWizardFilterSettings": {
                "$id": "1",
                "$type": "NINA.Core.Model.Equipment.FlatWizardFilterSettings, NINA.Core",
                "FlatWizardMode": 0,
                "HistogramMeanTarget": 0.5,
                "HistogramTolerance": 0.1,
                "MaxFlatExposureTime": 30.0,
                "MinFlatExposureTime": 0.01,
                "StepSize": 0.1,
                "MaxAbsoluteFlatDeviceBrightness": 1,
                "MinAbsoluteFlatDeviceBrightness": 0,
                "FlatDeviceAbsoluteStepSize": 1
                },
                "_autoFocusBinning": {
                "$id": "69",
                "$type": "NINA.Core.Model.Equipment.BinningMode, NINA.Core",
                "X": 1,
                "Y": 1
                },
                "_autoFocusGain": -1,
                "_autoFocusOffset": -1
                }
    return se_filter
#%%
count = 1
B = a.se_maker('./','g',count,count,1,1,0,1,1,'BIAS')
#%%

B['Items']['$values'][0]['Filter']
#%%

def NINAcalibration(exptime : list,
                    filter : list,
                    bias : bool,
                    dark : bool,
                    flat : bool,
                    count : int = 9,
                    return_action : bool = True):
    nina = NINA()
    Cd= nina.Condition()
    Act= nina.Action()
    Tr= nina.Trigger()
    maker= nina.Scriptmaker()
    cont_maker= nina.Container()
    

    

#%%