#%%
"""
- History
(23.07.10) Written by Hongjae Moon (V1.1)
(23.07.10) Modified (only NINA class) by Hyeonho Choi for connecting with the ObsScheduler
"""
import json
from astropy.io import ascii
#%%
class NINA_Container: #NINA script에서 가장 기본이 되는 뼈대
    
    def __init__(self):
        pass
    
    def base(self, name:str): #아무것도 없는 상태, 기본 토대
        value={
  "$id": "",
  "$type": "NINA.Sequencer.Container.SequenceRootContainer, NINA.Sequencer",
  "Strategy": {
    "$type": "NINA.Sequencer.Container.ExecutionStrategy.SequentialStrategy, NINA.Sequencer"
  },
  "Name": name, #이 name은 NINA advanced sequence 등에서 가장 위에 표시되는 이름
  "Conditions": {
    "$id": "",
    "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Conditions.ISequenceCondition, NINA.Sequencer]], System",
    "$values": []
  },
  "IsExpanded": True,
  "Items": {
    "$id": "",
    "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.SequenceItem.ISequenceItem, NINA.Sequencer]], System",
    "$values": [
      {
        "$id": "",
        "$type": "NINA.Sequencer.Container.StartAreaContainer, NINA.Sequencer",
        "Strategy": {
          "$type": "NINA.Sequencer.Container.ExecutionStrategy.SequentialStrategy, NINA.Sequencer"
        },
        "Name": "Start",
        "Conditions": {
          "$id": "",
          "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Conditions.ISequenceCondition, NINA.Sequencer]], System",
          "$values": []
        },
        "IsExpanded": True,
        "Items": {
          "$id": "",
          "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.SequenceItem.ISequenceItem, NINA.Sequencer]], System",
          "$values": []
        },
        "Triggers": {
          "$id": "",
          "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Trigger.ISequenceTrigger, NINA.Sequencer]], System",
          "$values": []
        },
        "Parent": {
          "$ref": ""
        },
        "ErrorBehavior": 0,
        "Attempts": 1
      },
      {
        "$id": "",
        "$type": "NINA.Sequencer.Container.TargetAreaContainer, NINA.Sequencer",
        "Strategy": {
          "$type": "NINA.Sequencer.Container.ExecutionStrategy.SequentialStrategy, NINA.Sequencer"
        },
        "Name": "Targets",
        "Conditions": {
          "$id": "",
          "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Conditions.ISequenceCondition, NINA.Sequencer]], System",
          "$values": []
        },
        "IsExpanded": True,
        "Items": {
          "$id": "",
          "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.SequenceItem.ISequenceItem, NINA.Sequencer]], System",
          "$values": []
        },
        "Triggers": {
          "$id": "",
          "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Trigger.ISequenceTrigger, NINA.Sequencer]], System",
          "$values": []
        },
        "Parent": {
          "$ref": ""
        },
        "ErrorBehavior": 0,
        "Attempts": 1
      },
      {
        "$id": "",
        "$type": "NINA.Sequencer.Container.EndAreaContainer, NINA.Sequencer",
        "Strategy": {
          "$type": "NINA.Sequencer.Container.ExecutionStrategy.SequentialStrategy, NINA.Sequencer"
        },
        "Name": "End",
        "Conditions": {
          "$id": "",
          "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Conditions.ISequenceCondition, NINA.Sequencer]], System",
          "$values": []
        },
        "IsExpanded": True,
        "Items": {
          "$id": "",
          "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.SequenceItem.ISequenceItem, NINA.Sequencer]], System",
          "$values": []
        },
        "Triggers": {
          "$id": "",
          "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Trigger.ISequenceTrigger, NINA.Sequencer]], System",
          "$values": []
        },
        "Parent": {
          "$ref": ""
        },
        "ErrorBehavior": 0,
        "Attempts": 1
      }
    ]
  },
  "Triggers": {
    "$id": "",
    "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Trigger.ISequenceTrigger, NINA.Sequencer]], System",
    "$values": []
  },
  "Parent": None,
  "ErrorBehavior": 0,
  "Attempts": 1
}
        return value
    

    def DSO(self, target:str, rotation:float, RAH:float, RAM:float, RAS:float,DecD:float, DecM:float, DecS:float, NegativeDec:bool): #target에 대한 정보가 들어가는 container
            
        value={
              "$id": "",
              "$type": "NINA.Sequencer.Container.DeepSkyObjectContainer, NINA.Sequencer",
              "Target": {
                "$id": "",
                "$type": "NINA.Astrometry.InputTarget, NINA.Astrometry",
                "Expanded": True,
                "TargetName": target,
                "Rotation": rotation,
                "InputCoordinates": {
                  "$id": "",
                  "$type": "NINA.Astrometry.InputCoordinates, NINA.Astrometry",
                  "RAHours": RAH,
                  "RAMinutes": RAM,
                  "RASeconds": RAS,
                  "NegativeDec": NegativeDec,
                  "DecDegrees": DecD,
                  "DecMinutes": DecM,
                  "DecSeconds": DecS
                }
              },
              "ExposureInfoListExpanded": True,
              "ExposureInfoList": {
                "$id": "",
                "$type": "NINA.Core.Utility.AsyncObservableCollection`1[[NINA.Sequencer.Utility.ExposureInfo, NINA.Sequencer]], NINA.Core",
                "$values": []
              },
              "Strategy": {
                "$type": "NINA.Sequencer.Container.ExecutionStrategy.SequentialStrategy, NINA.Sequencer"
              },
              "Name": target, #target이 container의 이름이 된다
              "Conditions": {
                "$id": "",
                "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Conditions.ISequenceCondition, NINA.Sequencer]], System",
                "$values": []
              },
              "IsExpanded": False,
              "Items": {
                "$id": "",
                "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.SequenceItem.ISequenceItem, NINA.Sequencer]], System",
                "$values": []
              },
              "Triggers": {
                "$id": "",
                "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Trigger.ISequenceTrigger, NINA.Sequencer]], System",
                "$values": []
              },
              "Parent": {
                "$ref": ""
              },
              "ErrorBehavior": 0,
              "Attempts": 1
            }

        return value
    
    def parallelContainer(self,cont_name:str): #trigger, condition 없이 action만 포함
            
        value= {
              "$id": "",
              "$type": "NINA.Sequencer.Container.ParallelContainer, NINA.Sequencer",
              "Strategy": {
                "$type": "NINA.Sequencer.Container.ExecutionStrategy.ParallelStrategy, NINA.Sequencer"
              },
              "Name": cont_name, #NINA 창에서 container를 명시하는 이름
              "Conditions": {
                "$id": "",
                "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Conditions.ISequenceCondition, NINA.Sequencer]], System",
                "$values": []
              },
              "IsExpanded": True,
              "Items": {
                "$id": "",
                "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.SequenceItem.ISequenceItem, NINA.Sequencer]], System",
                "$values": []
              },
              "Triggers": {
                "$id": "",
                "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Trigger.ISequenceTrigger, NINA.Sequencer]], System",
                "$values": []
              },
              "Parent": {
                "$ref": ""
              },
              "ErrorBehavior": 0,
              "Attempts": 1
            }

        return value
    
    def sequentialContainer(self,cont_name:str): #trigger, condition, action 3개가 모두 포함되는 기본적인 container
            
        value={
              "$id": "",
              "$type": "NINA.Sequencer.Container.SequentialContainer, NINA.Sequencer",
              "Strategy": {
                "$type": "NINA.Sequencer.Container.ExecutionStrategy.SequentialStrategy, NINA.Sequencer"
              },
              "Name": cont_name, #NINA 창에서 container를 명시하는 이름
              "Conditions": {
                "$id": "",
                "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Conditions.ISequenceCondition, NINA.Sequencer]], System",
                "$values": []
              },
              "IsExpanded": True,
              "Items": {
                "$id": "",
                "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.SequenceItem.ISequenceItem, NINA.Sequencer]], System",
                "$values": []
              },
              "Triggers": {
                "$id": "",
                "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Trigger.ISequenceTrigger, NINA.Sequencer]], System",
                "$values": []
              },
              "Parent": {
                "$ref": ""
              },
              "ErrorBehavior": 0,
              "Attempts": 1
            }

        return value
#%%x
class NINA_Action:
    
    """
    Actions for NINA script
    
    Methods
    =======
    1. smart_exp 
    2. camera_cooling
    3. camera_warming
    4. slew_AA
    5. slew_RD
    6. switch_filt
    7. autofocus
    8. move_focus
    9. fan_cont
    10. delt_cont_second
    11. delt_cont_prime
    12. startPW3
    13. stopPW3
    14. startPW4
    15. stopPW4
    16. wait_safe
    17. park
    18. unpark
    19. home
    20. set_track_side
    21. set_track_stop
    22. solve_sync
    23. ext_script
    24. wait_time
    25. wait_time_dup
    26. wait_timespan
    27. wait_sunaltitude
    28. wait_altitude_ab
    29. wait_altitude_bl
    
    """
    
    def __init__(self):
        pass
    
    def smart_exp(self,iterat:int,dither:int,exp_se:float,gain_se:int,offset_se:int,x_se:int,y_se:int,IT_se:str): #contant plugin을 설치하면 나오는 smart exposure
        value={
                          "$id": "",
                          "$type": "ConstantsPlugin.Constants.SmartExposure, ConstantsPlugin",
                          "ErrorBehavior": 0,
                          "Attempts": 1,
                          "IsValidSmartIterationCount": "#FFADFF2F",
                          "IterationsExpr": str(iterat), #iterat은 몇 장의 이미지를 찍을 지 정함
                          "IterationCountString": str(iterat),
                          "IterationCount": iterat,
                          "IsValidDitherCount": "#FFADFF2F",
                          "DitherExpr": str(dither), #dither는 몇 장의 노출마다 dithering을 진행하는지 정함
                          "DitherCountString": str(dither),
                          "DitherCount": dither,
                          "Strategy": {
                            "$type": "NINA.Sequencer.Container.ExecutionStrategy.SequentialStrategy, NINA.Sequencer"
                          },
                          "Name": "Smart Exposure",
                          "Conditions": {
                            "$id": "",
                            "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Conditions.ISequenceCondition, NINA.Sequencer]], System",
                            "$values": [
                              {
                                "$id": "",
                                "$type": "NINA.Sequencer.Conditions.LoopCondition, NINA.Sequencer",
                                "CompletedIterations": 0,
                                "Iterations": iterat,
                                "Parent": {
                                  "$ref": ""
                                }
                              }
                            ]
                          },
                          "IsExpanded": False,
                          "Items": {
                            "$id": "",
                            "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.SequenceItem.ISequenceItem, NINA.Sequencer]], System",
                            "$values": [
                              {
                                "$id": "",
                                "$type": "NINA.Sequencer.SequenceItem.FilterWheel.SwitchFilter, NINA.Sequencer",
                                "Filter": {},
                                "Parent": {
                                  "$ref": ""
                                },
                                "ErrorBehavior": 0,
                                "Attempts": 1
                              },
                              {
                                "$id": "",
                                "$type": "ConstantsPlugin.Constants.TakeExposure, ConstantsPlugin",
                                "ExposureTimeExpr": str(exp_se),
                                "ExposureTimeString": str(exp_se),
                                "ExposureTime": exp_se,
                                "Gain": gain_se,
                                "Offset": offset_se,
                                "Binning": {
                                  "$id": "",
                                  "$type": "NINA.Core.Model.Equipment.BinningMode, NINA.Core",
                                  "X": x_se,
                                  "Y":y_se
                                },
                                "ImageType": IT_se,
                                "ExposureCount": 0,
                                "IsValidExposureTime": "#FFADFF2F",
                                "Parent": {
                                  "$ref": ""
                                },
                                "ErrorBehavior": 0,
                                "Attempts": 1
                              }
                            ]
                          },
                          "Triggers": {
                            "$id": "",
                            "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Trigger.ISequenceTrigger, NINA.Sequencer]], System",
                            "$values": [
                              {
                                "$id": "",
                                "$type": "NINA.Sequencer.Trigger.Guider.DitherAfterExposures, NINA.Sequencer",
                                "AfterExposures": dither,
                                "Parent": {
                                  "$ref": ""
                                },
                                "TriggerRunner": {
                                  "$id": "",
                                  "$type": "NINA.Sequencer.Container.SequentialContainer, NINA.Sequencer",
                                  "Strategy": {
                                    "$type": "NINA.Sequencer.Container.ExecutionStrategy.SequentialStrategy, NINA.Sequencer"
                                  },
                                  "Name": None,
                                  "Conditions": {
                                    "$id": "",
                                    "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Conditions.ISequenceCondition, NINA.Sequencer]], System",
                                    "$values": []
                                  },
                                  "IsExpanded": True,
                                  "Items": {
                                    "$id": "",
                                    "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.SequenceItem.ISequenceItem, NINA.Sequencer]], System",
                                    "$values": [
                                      {
                                        "$id": "",
                                        "$type": "NINA.Sequencer.SequenceItem.Guider.Dither, NINA.Sequencer",
                                        "Parent": {
                                          "$ref": ""
                                        },
                                        "ErrorBehavior": 0,
                                        "Attempts": 1
                                      }
                                    ]
                                  },
                                  "Triggers": {
                                    "$id": "",
                                    "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Trigger.ISequenceTrigger, NINA.Sequencer]], System",
                                    "$values": []
                                  },
                                  "Parent": None,
                                  "ErrorBehavior": 0,
                                  "Attempts": 1
                                }
                              }
                            ]
                          },
                          "Parent": {
                            "$ref": ""
                          }
                        }

        return value
    
    def camera_cooling(self,temp_cam:float, dur_cam:float): #관측 시작 전 카메라 냉각
        value={
                    "$id": "",
                    "$type": "NINA.Sequencer.SequenceItem.Camera.CoolCamera, NINA.Sequencer",
                    "Temperature": temp_cam, #목표 온도
                    "Duration": dur_cam, #냉각 지속 시간
                    "Parent": {
                      "$ref": ""
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }

        return value
    
    def camera_warming(self,dur_cam_w:float): #관측 완료 후 카메라 가열
        value={
                    "$id": "",
                    "$type": "NINA.Sequencer.SequenceItem.Camera.WarmCamera, NINA.Sequencer",
                    "Duration": dur_cam_w, #가열 지속 시간
                    "Parent": {
                      "$ref": ""
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def slew_AA(self,AzD_t:int, AzM_t:int, AzS_t:float, AltD_t:int, AltM_t:int, AltS_t:float): #망원경 선회, 고도-방위각 기준
        value={"$id": "",
                    "$type": "NINA.Sequencer.SequenceItem.Telescope.SlewScopeToAltAz, NINA.Sequencer",
                    "Coordinates": {
                      "$id": "",
                      "$type": "NINA.Astrometry.InputTopocentricCoordinates, NINA.Astrometry",
                      "AzDegrees": AzD_t,
                      "AzMinutes": AzM_t,
                      "AzSeconds": AzS_t,
                      "AltDegrees": AltD_t,
                      "AltMinutes": AltM_t,
                      "AltSeconds": AltS_t
                    },
                    "Parent": {
                      "$ref": ""
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                    }
        
        return value
    
    def slew_RD(self,RAH_t:int,RAM_t:int,RAS_t:float,NegativeDec_t:bool, DecD_t:int, DecM_t:int,DecS_t:float): #망원경 선회, 적경-적위 기준
        value={
                    "$id": "",
                    "$type": "NINA.Sequencer.SequenceItem.Telescope.SlewScopeToRaDec, NINA.Sequencer",
                    "Inherited": False,
                    "Coordinates": {
                      "$id": "",
                      "$type": "NINA.Astrometry.InputCoordinates, NINA.Astrometry",
                      "RAHours": RAH_t,
                      "RAMinutes": RAM_t,
                      "RASeconds": RAS_t,
                      "NegativeDec": False,
                      "DecDegrees": DecD_t,
                      "DecMinutes": DecM_t,
                      "DecSeconds": DecS_t
                    },
                    "Parent": {
                      "$ref": ""
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def switch_filt(self, filt:str): #filter 바꿈, smart exposure의 경우 이미 filter 정보가 있기에 개별적으로 필터를 바꾸는 경우 사용
        value={
              "$id": "",
              "$type": "NINA.Sequencer.SequenceItem.FilterWheel.SwitchFilter, NINA.Sequencer",
              "Filter": {
                "$id": "",
                "$type": "NINA.Core.Model.Equipment.FilterInfo, NINA.Core",
                "_name": filt,
                "_focusOffset": 0,
                "_position": 1,
                "_autoFocusExposureTime": 4.0,
                "_autoFocusFilter": False,
                "FlatWizardFilterSettings": {
                  "$id": "",
                  "$type": "NINA.Core.Model.Equipment.FlatWizardFilterSettings, NINA.Core",
                  "FlatWizardMode": 0,
                  "HistogramMeanTarget": 0.3,
                  "HistogramTolerance": 0.2,
                  "MaxFlatExposureTime": 10.0,
                  "MinFlatExposureTime": 0.1,
                  "StepSize": 0.1,
                  "MaxAbsoluteFlatDeviceBrightness": 1,
                  "MinAbsoluteFlatDeviceBrightness": 0,
                  "FlatDeviceAbsoluteStepSize": 1
                },
                "_autoFocusBinning": {
                  "$id": "",
                  "$type": "NINA.Core.Model.Equipment.BinningMode, NINA.Core",
                  "X": 1,
                  "Y": 1
                },
                "_autoFocusGain": -1,
                "_autoFocusOffset": -1
              },
              "Parent": {
                "$ref": ""
              },
              "ErrorBehavior": 0,
              "Attempts": 1
            }
 
        return value
    
    def autofocus(self): #auto focus 진행
        value={
                    "$id": "",
                    "$type": "NINA.Sequencer.SequenceItem.Autofocus.RunAutofocus, NINA.Sequencer",
                    "Parent": {
                      "$ref": ""
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def move_focus_relative(self, pos_f:int): #focuser 상대적으로 위치 변경
        value=                  {
                    "$id": "",
                    "$type": "NINA.Sequencer.SequenceItem.Focuser.MoveFocuserRelative, NINA.Sequencer",
                    "RelativePosition": pos_f,
                    "Parent": {
                      "$ref": ""
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }

        return value
      
    def move_focus(self, pos_f:int): #focuser 위치 변경
        value={
                    "$id": "",
                    "$type": "NINA.Sequencer.SequenceItem.Focuser.MoveFocuserAbsolute, NINA.Sequencer",
                    "Position": pos_f, 
                    "Parent": {
                      "$ref": ""
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }

        return value
    
    def fan_cont(self, fanstate:bool): #냉각팬 조절
        value={
                    "$id": "",
                    "$type": "DaleGhent.NINA.PlaneWaveTools.Fans.FanControl, PlaneWave Tools",
                    "FanState": fanstate, #on, off 상태 True or False로 조정
                    "Parent": {
                      "$ref": ""
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
 
        return value
    
    def delt_cont_second(self, Delt_HM:int): #secondary mirror heater 조절
        value={
                    "$id": "",
                    "$type": "DaleGhent.NINA.PlaneWaveTools.DeltaT.DeltaTControl, PlaneWave Tools",
                    "DeltaTHeater": 1,
                    "DeltaTHeaterMode": Delt_HM,#HM=0 is off, 1 is on
                    "Parent": {
                      "$ref": ""
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
      
        return value
    
    def delt_cont_prime(self, Delt_HM:int): #primary backplate heater 조절
        value={
                    "$id": "",
                    "$type": "DaleGhent.NINA.PlaneWaveTools.DeltaT.DeltaTControl, PlaneWave Tools",
                    "DeltaTHeater": 0,
                    "DeltaTHeaterMode": Delt_HM,#HM=0 is off, 1 is on
                    "Parent": {
                      "$ref": ""
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def startPW3(self): #plane wave intergace 3 시작
        value={
                    "$id": "",
                    "$type": "DaleGhent.NINA.PlaneWaveTools.StartStopPwi3.StartPwi3, PlaneWave Tools",
                    "Parent": {
                      "$ref": ""
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def stopPW3(self): #plane wave interface 3 중지
        value={
                    "$id": '',
                    "$type": "DaleGhent.NINA.PlaneWaveTools.StartStopPwi3.StopPwi3, PlaneWave Tools",
                    "Parent": {
                      "$ref": ''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def startPW4(self): #plane wave interface 4 시작
        value={
                    "$id": '',
                    "$type": "DaleGhent.NINA.PlaneWaveTools.StartStopPwi4.StartPwi4, PlaneWave Tools",
                    "Parent": {
                      "$ref":''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
       
        return value
    
    def stopPW4(self): #plane wave interface 4 중지
        value={
                    "$id": '',
                    "$type": "DaleGhent.NINA.PlaneWaveTools.StartStopPwi4.StopPwi4, PlaneWave Tools",
                    "Parent": {
                      "$ref": ''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def wait_safe(self): #wait until safe
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.SequenceItem.SafetyMonitor.WaitUntilSafe, NINA.Sequencer",
                    "Parent": {
                      "$ref": ''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def park(self): #망원경 parking
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.SequenceItem.Telescope.ParkScope, NINA.Sequencer",
                    "Parent": {
                      "$ref": ''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def unpark(self): #망원경 unparking
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.SequenceItem.Telescope.UnparkScope, NINA.Sequencer",
                    "Parent": {
                      "$ref": ''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def home(self): #망원경을 home position으로 움직임, 망원경이 parking된 상태면 움직이지 않음
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.SequenceItem.Telescope.FindHome, NINA.Sequencer",
                    "Parent": {
                      "$ref": ''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def set_track_side(self): #mount의 tracking 모드 결정, sidereal
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.SequenceItem.Telescope.SetTracking, NINA.Sequencer",
                    "TrackingMode": 0,
                    "Parent": {
                      "$ref": ''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def set_track_stop(self): #mount의 tracking 모드 결정, stopped
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.SequenceItem.Telescope.SetTracking, NINA.Sequencer",
                    "TrackingMode": 5,
                    "Parent": {
                      "$ref": ''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def solve_sync(self): #해당 망원경 위치 기준으로 take and solve 이미지, 해당 망원경 위치 기준이지 target 기준이 아니다
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.SequenceItem.Platesolving.SolveAndSync, NINA.Sequencer",
                    "Parent": {
                      "$ref": ''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def ext_script(self,script:str): #원하는 script 추가
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.SequenceItem.Utility.ExternalScript, NINA.Sequencer",
                    "Script": script,
                    "Parent": {
                      "$ref": ''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def wait_time(self,wH:int, wM:int,wMoff:int,wS:int): #특정 시간까지 대기
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.SequenceItem.Utility.WaitForTime, NINA.Sequencer",
                    "Hours": wH,
                    "Minutes": wM,
                    "MinutesOffset": wMoff,
                    "Seconds": wS,
                    "SelectedProvider": {
                      "$id": '',
                      "$type": "NINA.Sequencer.Utility.DateTimeProvider.TimeProvider, NINA.Sequencer"
                    },
                    "Parent": {
                      "$ref": ''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def wait_time_dup(self,wH:int, wM:int,wMoff:int,wS:int): #wati for time가 2번 이상 사용된 경우 'SelectedProvider'를 ref로 받아옴 
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.SequenceItem.Utility.WaitForTime, NINA.Sequencer",
                    "Hours": wH,
                    "Minutes": wM,
                    "MinutesOffset": wMoff,
                    "Seconds": wS,
                    "SelectedProvider": {
                      "$ref": ''
                    },
                    "Parent": {
                      "$ref": ''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def wait_timespan(self,wT:int): #몇 초 대기
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.SequenceItem.Utility.WaitForTimeSpan, NINA.Sequencer",
                    "Time": wT,
                    "Parent": {
                      "$ref": ''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
       
        return value
    
    def wait_sunaltitude(self, wRAH_s:int,wRAM_s:int,wRAS_s:float,wNegativeDec_s:bool, wDecD_s:int,wDecM_s:int,wDecS_s:float,woff_s:float,wcomp_s:int): #태양의 특정 고도 기준으로 <,>가 충족될 때까지 대기
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.SequenceItem.Utility.WaitForSunAltitude, NINA.Sequencer",
                    "Data": {
                      "$id": '',
                      "$type": "NINA.Sequencer.SequenceItem.Utility.WaitLoopData, NINA.Sequencer",
                      "Coordinates": {
                        "$id": '',
                        "$type": "NINA.Astrometry.InputCoordinates, NINA.Astrometry",
                        "RAHours": wRAH_s,
                        "RAMinutes": wRAM_s,
                        "RASeconds": wRAS_s,
                        "NegativeDec": wNegativeDec_s,
                        "DecDegrees": wDecD_s,
                        "DecMinutes": wDecM_s,
                        "DecSeconds": wDecS_s
                      },
                      "Offset": woff_s,
                      "Comparator": wcomp_s # 1이면 <, 3이면 >
                    },
                    "Parent": {
                      "$ref": ''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
    
    def wait_altitude_ab(self,Dso_Parent:bool,wRAH:int,wRAM:int,wRAS:float,wNegativeDec:bool, wDecD:int,wDecM:int,wDecS:float,woff:float): #특정 고도를 넘기까지 대기
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.SequenceItem.Utility.WaitForAltitude, NINA.Sequencer",
                    "AboveOrBelow": '>',
                    "HasDsoParent": Dso_Parent,
                    "Data": {
                      "$id": '',
                      "$type": "NINA.Sequencer.SequenceItem.Utility.WaitLoopData, NINA.Sequencer",
                      "Coordinates": {
                        "$id": '',
                        "$type": "NINA.Astrometry.InputCoordinates, NINA.Astrometry",
                        "RAHours": wRAH,
                        "RAMinutes": wRAM,
                        "RASeconds": wRAS,
                        "NegativeDec": wNegativeDec,
                        "DecDegrees": wDecD,
                        "DecMinutes": wDecM,
                        "DecSeconds": wDecS
                      },
                      "Offset": woff,
                      "Comparator": 3
                    },
                    "Parent": {
                      "$ref": ''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        return value
    
    def wait_altitude_bl(self,Dso_Parent:bool,wRAH:int,wRAM:int,wRAS:float,wNegativeDec:bool, wDecD:int,wDecM:int,wDecS:float,woff:float): #특정 고도 아래로 떨어질 때까지 대기
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.SequenceItem.Utility.WaitForAltitude, NINA.Sequencer",
                    "AboveOrBelow": '<',
                    "HasDsoParent": Dso_Parent,
                    "Data": {
                      "$id": '',
                      "$type": "NINA.Sequencer.SequenceItem.Utility.WaitLoopData, NINA.Sequencer",
                      "Coordinates": {
                        "$id": '',
                        "$type": "NINA.Astrometry.InputCoordinates, NINA.Astrometry",
                        "RAHours": wRAH,
                        "RAMinutes": wRAM,
                        "RASeconds": wRAS,
                        "NegativeDec": wNegativeDec,
                        "DecDegrees": wDecD,
                        "DecMinutes": wDecM,
                        "DecSeconds": wDecS
                      },
                      "Offset": woff,
                      "Comparator": 2
                    },
                    "Parent": {
                      "$ref": ''
                    },
                    "ErrorBehavior": 0,
                    "Attempts": 1
                  }
        
        return value
#%%
class NINA_Condition:

    """
    Conditions for NINA script
    
    Methods
    =======
    1. sun_alt 
    2. safe
    3. unsafe
    4. timespan
    """
    
    def __init__(self):
        pass
    
    def sun_alt(self, RAH_S:float, RAM_S:float,RAS_S:float,NegativeDec_S:bool, DecD_S:float,DecM_S:float,DecS_S:float, offset:float, Comp:int): #태양의 고도가 조건을 만족하면 loop 벗어남
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.Conditions.SunAltitudeCondition, NINA.Sequencer",
                    "Data": {
                      "$id": '',
                      "$type": "NINA.Sequencer.SequenceItem.Utility.WaitLoopData, NINA.Sequencer",
                      "Coordinates": {
                        "$id": '',
                        "$type": "NINA.Astrometry.InputCoordinates, NINA.Astrometry",
                        "RAHours": RAH_S,
                        "RAMinutes": RAM_S,
                        "RASeconds": RAS_S,
                        "NegativeDec": NegativeDec_S,
                        "DecDegrees": DecD_S,
                        "DecMinutes": DecM_S,
                        "DecSeconds": DecS_S
                      },
                      "Offset": offset,
                      "Comparator": Comp
                    },
                    "Parent": {
                      "$ref": ''
                    }
                  }

        return value
      
    def altitude(self, Dso_Parent:bool,wRAH:int,wRAM:int,wRAS:float,wNegativeDec:bool, wDecD:int,wDecM:int,wDecS:float,woff:float):
      value = {
                          "$id": "",
                          "$type": "NINA.Sequencer.Conditions.AboveHorizonCondition, NINA.Sequencer",
                          "HasDsoParent": Dso_Parent,
                          "Data": {
                            "$id": "",
                            "$type": "NINA.Sequencer.SequenceItem.Utility.WaitLoopData, NINA.Sequencer",
                            "Coordinates": {
                              "$id": "",
                              "$type": "NINA.Astrometry.InputCoordinates, NINA.Astrometry",
                              "RAHours": wRAH,
                              "RAMinutes": wRAM,
                              "RASeconds": wRAS,
                              "NegativeDec": wNegativeDec,
                              "DecDegrees": wDecD,
                              "DecMinutes": wDecM,
                              "DecSeconds": wDecS
                            },
                            "Offset": woff,
                            "Comparator": 3
                          },
                          "Parent": {
                            "$ref": "14"
                          }
                        }
      return value
    
    def safe(self): #safety monitor가 safe로 표시되는 동안 계속 loop
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.Conditions.SafetyMonitorCondition, NINA.Sequencer",
                    "Parent": {
                      "$ref": ''
                    }
                  }
        
        return value
    
    def unsafe(self): #safety monitor가 unsafe로 표시되는 동안 계속 loop
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.Conditions.LoopWhileUnsafe, NINA.Sequencer",
                    "Parent": {
                      "$ref": ''
                    }
                  }
       
        return value
    
    def timespan(self,tsH:int,tsM:int,tsS:float): #정해진 시간 동안 loop, 만약 정해진 시간이 끝났는데 loop와 같이 묶인 행동들이 끝나지 않았다면 생략됨
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.Conditions.TimeSpanCondition, NINA.Sequencer",
                    "Hours": tsH,
                    "Minutes": tsM,
                    "Seconds": tsS,
                    "Parent": {
                      "$ref": ''
                    }
                  }
        
        return value
#%%
class NINA_Trigger:
    """
    Trigger for the NINA script
    
    Methods
    =======
    1. AF_filterchan
    2. AF_exp
    3. AF_time
    4. meridianflip
    """
    def __init__(self):
        pass
    
    def AF_filterchan(self):
      value = {
                    "$id": "",
                    "$type": "NINA.Sequencer.Trigger.Autofocus.AutofocusAfterFilterChange, NINA.Sequencer",
                    "Parent": {
                      "$ref": ""
                    },
                    "TriggerRunner": {
                      "$id": "",
                      "$type": "NINA.Sequencer.Container.SequentialContainer, NINA.Sequencer",
                      "Strategy": {
                        "$type": "NINA.Sequencer.Container.ExecutionStrategy.SequentialStrategy, NINA.Sequencer"
                      },
                      "Name": None,
                      "Conditions": {
                        "$id": "",
                        "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Conditions.ISequenceCondition, NINA.Sequencer]], System",
                        "$values": []
                      },
                      "IsExpanded": True,
                      "Items": {
                        "$id": "",
                        "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.SequenceItem.ISequenceItem, NINA.Sequencer]], System",
                        "$values": [
                          {
                            "$id": "",
                            "$type": "NINA.Sequencer.SequenceItem.Autofocus.RunAutofocus, NINA.Sequencer",
                            "Parent": {
                              "$ref": ""
                            },
                            "ErrorBehavior": 0,
                            "Attempts": 1
                          }
                        ]
                      },
                      "Triggers": {
                        "$id": "",
                        "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Trigger.ISequenceTrigger, NINA.Sequencer]], System",
                        "$values": []
                      },
                      "Parent": None,
                      "ErrorBehavior": 0,
                      "Attempts": 1
                    }
                  }
      return value 
    
    #@property
    def AF_exp(self,af_exp:int): #몇 장의 노출 후 autofocus를 할 것인가
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.Trigger.Autofocus.AutofocusAfterExposures, NINA.Sequencer",
                    "AfterExposures": af_exp,
                    "Parent": {
                      "$ref": ''
                    },
                    "TriggerRunner": {
                      "$id": '',
                      "$type": "NINA.Sequencer.Container.SequentialContainer, NINA.Sequencer",
                      "Strategy": {
                        "$type": "NINA.Sequencer.Container.ExecutionStrategy.SequentialStrategy, NINA.Sequencer"
                      },
                      "Name": None,
                      "Conditions": {
                        "$id": '',
                        "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Conditions.ISequenceCondition, NINA.Sequencer]], System",
                        "$values": []
                      },
                      "IsExpanded": True,
                      "Items": {
                        "$id": '',
                        "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.SequenceItem.ISequenceItem, NINA.Sequencer]], System",
                        "$values": [
                          {
                            "$id": '',
                            "$type": "NINA.Sequencer.SequenceItem.Autofocus.RunAutofocus, NINA.Sequencer",
                            "Parent": {
                              "$ref": ''
                            },
                            "ErrorBehavior": 0,
                            "Attempts": 1
                          }
                        ]
                      },
                      "Triggers": {
                        "$id": '',
                        "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Trigger.ISequenceTrigger, NINA.Sequencer]], System",
                        "$values": []
                      },
                      "Parent": None,
                      "ErrorBehavior": 0,
                      "Attempts": 1
                    }
                  }
        
        return value
    
    def AF_time(self,amt:float): #몇 분 뒤에 autofocus를 할 것인가
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.Trigger.Autofocus.AutofocusAfterTimeTrigger, NINA.Sequencer",
                    "Amount": amt,
                    "Parent": {
                      "$ref": ''
                    },
                    "TriggerRunner": {
                      "$id": '',
                      "$type": "NINA.Sequencer.Container.SequentialContainer, NINA.Sequencer",
                      "Strategy": {
                        "$type": "NINA.Sequencer.Container.ExecutionStrategy.SequentialStrategy, NINA.Sequencer"
                      },
                      "Name": None,
                      "Conditions": {
                        "$id": '',
                        "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Conditions.ISequenceCondition, NINA.Sequencer]], System",
                        "$values": []
                      },
                      "IsExpanded": True,
                      "Items": {
                        "$id": '',
                        "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.SequenceItem.ISequenceItem, NINA.Sequencer]], System",
                        "$values": [
                          {
                            "$id": '',
                            "$type": "NINA.Sequencer.SequenceItem.Autofocus.RunAutofocus, NINA.Sequencer",
                            "Parent": {
                              "$ref": ''
                            },
                            "ErrorBehavior": 0,
                            "Attempts": 1
                          }
                        ]
                      },
                      "Triggers": {
                        "$id": '',
                        "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Trigger.ISequenceTrigger, NINA.Sequencer]], System",
                        "$values": []
                      },
                      "Parent": None,
                      "ErrorBehavior": 0,
                      "Attempts": 1
                    }
                  }
        
        return value
    
    def meridianflip(self): #조건 충족시에 meridian flip 진행
        value={
                    "$id": '',
                    "$type": "NINA.Sequencer.Trigger.MeridianFlip.MeridianFlipTrigger, NINA.Sequencer",
                    "Parent": {
                      "$ref": ''
                    },
                    "TriggerRunner": {
                      "$id": '',
                      "$type": "NINA.Sequencer.Container.SequentialContainer, NINA.Sequencer",
                      "Strategy": {
                        "$type": "NINA.Sequencer.Container.ExecutionStrategy.SequentialStrategy, NINA.Sequencer"
                      },
                      "Name": None,
                      "Conditions": {
                        "$id": '',
                        "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Conditions.ISequenceCondition, NINA.Sequencer]], System",
                        "$values": []
                      },
                      "IsExpanded": True,
                      "Items": {
                        "$id": '',
                        "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.SequenceItem.ISequenceItem, NINA.Sequencer]], System",
                        "$values": []
                      },
                      "Triggers": {
                        "$id": '',
                        "$type": "System.Collections.ObjectModel.ObservableCollection`1[[NINA.Sequencer.Trigger.ISequenceTrigger, NINA.Sequencer]], System",
                        "$values": []
                      },
                      "Parent": None,
                      "ErrorBehavior": 0,
                      "Attempts": 1
                    }
                  }
        
        return value
#%%
class NINA_ScriptMaker(NINA_Container, NINA_Action, NINA_Condition, NINA_Trigger):
    """
    
    Methods
    =======
    1. write
    2. edit
    3. se_maker
    4. add_container
    5. add_condition
    6. add_trigger
    7. add_action
    8. update_dict_ids
    9. provider_update

    """
    def __init__(self, filename : str = 'NINAscript.json'):
        super().__init__() # __init__(여기에 name, 등을 넣으면 이전에 정의한 값을 상속)
        self.filename = filename
        
        cont_maker = NINA_Container()
        self.base=cont_maker.base('test')
    
    def __repr__(self):
        #pass
        return f'NINAScriptMaker[filename = {self.filename}]'
    
    def write(self,path:str,name:str,js:dict): #최종적으로 NINA에서 돌릴 수 있는 script 작성
        Name=path+name
        with open(Name, 'w') as outfile:
            json.dump(js,outfile, indent = 2)
            
    def edit(self,filename:str,): #NINA에서 돌리기 이전에 편집을 위해 불러오기
        with open(filename, 'r') as f:
            json_data = json.load(f)
    
    def se_maker(self,path:str,filt:str,iterat:int,dither:int,exp:float,gain:int,offset:int,binx:int,biny:int,imagetype:str): #constant plugin의 smart exposure를 사용하기 위한 최종 단계, action의 smart exposure에 filter 정보 추가

        Ac = NINA_Action()
        #Con=Config()
        se=Ac.smart_exp(iterat,dither,exp,gain,offset,binx,biny,imagetype)
        
        with open(path+'/filter.json','r') as f:
            json_data=json.load(f)
            
        for f in range(len(json_data['Items']['$values'][1]['Items']['$values'])):
            if filt==json_data['Items']['$values'][1]['Items']['$values'][f]['Items']['$values'][0]['Filter']['_name']:
                fil=json_data['Items']['$values'][1]['Items']['$values'][f]['Items']['$values'][0]['Filter']
        Filter=fil
        #Fitler=Conf.find_filter(path,filt)
            
        se['Items']['$values'][0]['Filter']=Filter
        #print(se)
        return se
    
    def add_container(self,container:dict,pos:str,*made:dict): #원하는 container 종류를, 원하는 위치에, 기존 dictionary가 있다면 사용하고 아니면 새로 만들기
        base=self.base
        
        if len(made)==0:
            frame=base
            if frame['Items']['$values'][0]['Name']==pos:
                frame['Items']['$values'][0]['Items']['$values'].append(container)
            #cont=frame
            elif frame['Items']['$values'][1]['Name']==pos:
                frame['Items']['$values'][1]['Items']['$values'].append(container)
            #cont=frame
            elif frame['Items']['$values'][2]['Name']==pos:
                frame['Items']['$values'][2]['Items']['$values'].append(container)
            #cont=frame
        
        else:
            frame=made[0]
            if frame['Items']['$values'][0]['Name']==pos:
                frame['Items']['$values'][0]['Items']['$values'].append(container)
            #cont=frame
            elif frame['Items']['$values'][1]['Name']==pos:
                frame['Items']['$values'][1]['Items']['$values'].append(container)
            #cont=frame
            elif frame['Items']['$values'][2]['Name']==pos:
                frame['Items']['$values'][2]['Items']['$values'].append(container)
            #cont=frame
                
            else:
                for i in range(len(frame['Items']['$values'][0]['Items']['$values'])):
                    if frame['Items']['$values'][0]['Items']['$values'][i]['Name']==pos:
                        frame['Items']['$values'][0]['Items']['$values'][i]['Items']['$values'].append(container)
                        #cont=frame
                    else:
                        pass

                for j in range(len(frame['Items']['$values'][1]['Items']['$values'])):
                    if frame['Items']['$values'][1]['Items']['$values'][j]['Name']==pos:
                        frame['Items']['$values'][1]['Items']['$values'][j]['Items']['$values'].append(container)
                        #cont=frame
                    else:
                        pass

                for k in range(len(frame['Items']['$values'][2]['Items']['$values'])):
                    if frame['Items']['$values'][2]['Items']['$values'][k]['Name']==pos:
                        frame['Items']['$values'][2]['Items']['$values'][k]['Items']['$values'].append(container)
                        #cont=frame
                    else:
                        pass
       
        cont=frame
        #print(cont)
        return cont
    
    def add_condition(self,condition:dict,pos:str,*made:dict): #원하는 condition 종류를, 원하는 위치에, 기존 dictionary가 있다면 사용하고 아니면 새로 만들기
        frame=made[0]
        for i in range(len(frame['Items']['$values'][0]['Items']['$values'])):
            if frame['Items']['$values'][0]['Items']['$values'][i]['Name']==pos:
                frame['Items']['$values'][0]['Items']['$values'][i]['Conditions']['$values'].append(condition)
                condi=frame
            else:
                for ii in range(len(frame['Items']['$values'][0]['Items']['$values'][i]['Items']['$values'])):
                    if frame['Items']['$values'][0]['Items']['$values'][i]['Items']['$values'][ii]['Name']==pos:
                        frame['Items']['$values'][0]['Items']['$values'][i]['Items']['$values'][ii]['Conditions']['$values'].append(condition)
                        condi=frame
                    else:
                        pass

        for j in range(len(frame['Items']['$values'][1]['Items']['$values'])):
            if frame['Items']['$values'][1]['Items']['$values'][j]['Name']==pos:
                frame['Items']['$values'][1]['Items']['$values'][j]['Conditions']['$values'].append(condition)
                condi=frame
            else:
                for jj in range(len(frame['Items']['$values'][1]['Items']['$values'][j]['Items']['$values'])):
                    if frame['Items']['$values'][1]['Items']['$values'][j]['Items']['$values'][jj]['Name']==pos:
                        frame['Items']['$values'][1]['Items']['$values'][j]['Items']['$values'][jj]['Conditions']['$values'].append(condition)
                        condi=frame
                    else:
                        pass

        for k in range(len(frame['Items']['$values'][2]['Items']['$values'])):
            if frame['Items']['$values'][2]['Items']['$values'][k]['Name']==pos:
                frame['Items']['$values'][2]['Items']['$values'][k]['Conditions']['$values'].append(condition)
                condi=frame
            else:
                for kk in range(len(frame['Items']['$values'][2]['Items']['$values'][k]['Items']['$values'])):
                    if frame['Items']['$values'][2]['Items']['$values'][k]['Items']['$values'][kk]['Name']==pos:
                        frame['Items']['$values'][2]['Items']['$values'][k]['Items']['$values'][kk]['Conditions']['$values'].append(condition)
                        condi=frame
                    else:
                        pass
        #print(condi)  
        return condi
    
    def add_trigger(self,trigger:dict,pos:str,*made:dict): #원하는 trigger 종류를, 원하는 위치에, 기존 dictionary가 있다면 사용하고 아니면 새로 만들기
        base=self.base
        
        if len(made)==0:
            frame=base
            frame['Triggers']['$values'].append(trigger)
            trig=frame
            
        else:
            frame=made[0]
            if pos=='Triggers':
                frame['Triggers']['$values'].append(trigger)
                trig=frame
            else:
                for i in range(len(frame['Items']['$values'][0]['Items']['$values'])):
                    if frame['Items']['$values'][0]['Items']['$values'][i]['Name']==pos:
                        frame['Items']['$values'][0]['Items']['$values'][i]['Triggers']['$values'].append(trigger)
                        trig=frame
                    else:
                        for ii in range(len(frame['Items']['$values'][0]['Items']['$values'][i]['Items']['$values'])):
                            if frame['Items']['$values'][0]['Items']['$values'][i]['Items']['$values'][ii]['Name']==pos:
                                frame['Items']['$values'][0]['Items']['$values'][i]['Items']['$values'][ii]['Triggers']['$values'].append(trigger)
                                trig=frame
                            else:
                                pass

                for j in range(len(frame['Items']['$values'][1]['Items']['$values'])):
                    if frame['Items']['$values'][1]['Items']['$values'][j]['Name']==pos:
                        frame['Items']['$values'][1]['Items']['$values'][j]['Triggers']['$values'].append(trigger)
                        trig=frame
                    else:
                        for jj in range(len(frame['Items']['$values'][1]['Items']['$values'][j]['Items']['$values'])):
                          if 'Name' in frame['Items']['$values'][1]['Items']['$values'][j]['Items']['$values'][jj]:
                            if frame['Items']['$values'][1]['Items']['$values'][j]['Items']['$values'][jj]['Name']==pos:
                                frame['Items']['$values'][1]['Items']['$values'][j]['Items']['$values'][jj]['Triggers']['$values'].append(trigger)
                                trig=frame
                            else:
                                pass

                for k in range(len(frame['Items']['$values'][2]['Items']['$values'])):
                    if frame['Items']['$values'][2]['Items']['$values'][k]['Name']==pos:
                        frame['Items']['$values'][2]['Items']['$values'][k]['Triggers']['$values'].append(trigger)
                        trig=frame
                    else:
                        for kk in range(len(frame['Items']['$values'][2]['Items']['$values'][k]['Items']['$values'])):
                            if frame['Items']['$values'][2]['Items']['$values'][k]['Items']['$values'][kk]['Name']==pos:
                                frame['Items']['$values'][2]['Items']['$values'][k]['Items']['$values'][kk]['Triggers']['$values'].append(trigger)
                                trig=frame
                            else:
                                pass
        #print(trig)
        return trig
    
    def add_action(self,action:dict,pos:str,*made:dict): #원하는 action 종류를, 원하는 위치에, 기존 dictionary가 있다면 사용하고 아니면 새로 만들기
        base=self.base
        
        if len(made)==0:
            frame=base
            if frame['Items']['$values'][0]['Name']==pos:
                frame['Items']['$values'][0]['Items']['$values'].append(action)
                act=frame
            elif frame['Items']['$values'][1]['Name']==pos:
                frame['Items']['$values'][1]['Items']['$values'].append(action)
                act=frame
            elif frame['Items']['$values'][2]['Name']==pos:
                frame['Items']['$values'][2]['Items']['$values'].append(action)
                act=frame
            
        else:
            frame=made[0]
            if frame['Items']['$values'][0]['Name']==pos:
                frame['Items']['$values'][0]['Items']['$values'].append(action)
                act=frame
            elif frame['Items']['$values'][1]['Name']==pos:
                frame['Items']['$values'][1]['Items']['$values'].append(action)
                act=frame
            elif frame['Items']['$values'][2]['Name']==pos:
                frame['Items']['$values'][2]['Items']['$values'].append(action)
                act=frame
                
            else:
                for i in range(len(frame['Items']['$values'][0]['Items']['$values'])):
                    if frame['Items']['$values'][0]['Items']['$values'][i]['Name']==pos:
                        frame['Items']['$values'][0]['Items']['$values'][i]['Items']['$values'].append(action)
                        act=frame
                    else:
                        for ii in range(len(frame['Items']['$values'][0]['Items']['$values'][i]['Items']['$values'])):
                            if frame['Items']['$values'][0]['Items']['$values'][i]['Items']['$values'][ii]['Name']==pos:
                                frame['Items']['$values'][0]['Items']['$values'][i]['Items']['$values'][ii]['Items']['$values'].append(action)
                                act=frame
                            else:
                                pass

                for j in range(len(frame['Items']['$values'][1]['Items']['$values'])):
                    if frame['Items']['$values'][1]['Items']['$values'][j]['Name']==pos:
                        frame['Items']['$values'][1]['Items']['$values'][j]['Items']['$values'].append(action)
                        act=frame
                    else:
                        for jj in range(len(frame['Items']['$values'][1]['Items']['$values'][j]['Items']['$values'])):
                            if 'Name' in frame['Items']['$values'][1]['Items']['$values'][j]['Items']['$values'][jj]:
                                if frame['Items']['$values'][1]['Items']['$values'][j]['Items']['$values'][jj]['Name']==pos:
                                    frame['Items']['$values'][1]['Items']['$values'][j]['Items']['$values'][jj]['Items']['$values'].append(action)
                                    act=frame
                            else:
                                pass

                for k in range(len(frame['Items']['$values'][2]['Items']['$values'])):
                    if frame['Items']['$values'][2]['Items']['$values'][k]['Name']==pos:
                        frame['Items']['$values'][2]['Items']['$values'][k]['Items']['$values'].append(action)
                        act=frame
                    else:
                        for kk in range(len(frame['Items']['$values'][2]['Items']['$values'][k]['Items']['$values'])):
                            if frame['Items']['$values'][2]['Items']['$values'][k]['Items']['$values'][kk]['Name']==pos:
                                frame['Items']['$values'][2]['Items']['$values'][k]['Items']['$values'][kk]['Items']['$values'].append(action)
                                act=frame
                            else:
                                pass
        #print(act)
        return act
    
    def update_dict_ids(self, dictionary, id_idx=1): #합쳐진 sequence에 id, ref를 넣어서 NINA를 위한 json 파일 형태 완성
        if isinstance(dictionary, dict):
            if "$id" in dictionary:
                dictionary["$id"] = str(id_idx)
                ref=id_idx
                for key in dictionary.keys():
                    X=dictionary[key]
                    if isinstance(X,dict):
                        for keys in X.keys():
                            Y=X[keys]
                            if isinstance(Y,list):
                                for i in range(len(Y)):
                                    Y[i]['Parent']['$ref']=str(ref)
                            else:
                                pass
                    else:
                        pass
                id_idx += 1
            for value in dictionary.values():
                _,id_idx = self.update_dict_ids(value, id_idx)
        elif isinstance(dictionary, list):
            for item in dictionary:
                _,id_idx = self.update_dict_ids(item, id_idx)
        return dictionary, id_idx
    
    def provider_update(self, dictionary):#일단 wait for time만 고려, 나중에 고칠거면 '$type' 부분 손보기
        if isinstance(dictionary, dict):
            if "$type" in dictionary:
                if dictionary["$type"] == 'NINA.Sequencer.SequenceItem.Utility.WaitForTime, NINA.Sequencer':
                        if '$id' in dictionary['SelectedProvider'].keys():
                            global ref_wt
                            ref_wt=dictionary['SelectedProvider']['$id']
                        else:
                            dictionary['SelectedProvider']['$ref']=ref_wt
                else:
                    pass
            for value in dictionary.values():
                _= self.provider_update(value)
        elif isinstance(dictionary, list):
            #for leng in range(len(dictionary)):
                #if "$type" in dictionary[leng]:
                    #if dictionary[leng]["$type"] == 'NINA.Sequencer.SequenceItem.Utility.WaitForTime, NINA.Sequencer':
                            #if '$id' in dictionary[leng]['SelectedProvider'].keys():
                                #ref_wt=dictionary[leng]['SelectedProvider']['$id']
                                #print(ref_wt)
                            #else:
                                #dictionary[leng]['SelectedProvider']['$ref']=ref_wt
                    #else:
                        #pass
            for item in dictionary:
                _= self.provider_update(item)
        return dictionary
    """
    def remove_duplicate_targetname(self, dictionary):
      if isinstance(dictionary, dict):
        for key, value in dictionary.items():
          if key in ['Name', 'TargetName']:
            if isinstance(value, str):
              split_value = value.split("dup_")[0]
              dictionary[key] = split_value
            elif isinstance(value, (dict, list)):
              self.remove_duplicate_targetname(value)            
      elif isinstance(dictionary, list):
        for item in dictionary:
          if isinstance(item, (dict, list)):
            self.remove_duplicate_targetname(item)
      return dictionary
    """
    def remove_duplicate_targetname(self, dictionary):
      if isinstance(dictionary, dict):
        for key in ['Name', 'TargetName']:
          if key in dictionary:
            val = dictionary[key]
            if isinstance(val, str):
              val_list = val.split('dup_')
              dictionary[key] = val_list[0]
          
        for value in dictionary.values():
          _ = self.remove_duplicate_targetname(value)            
      
      elif isinstance(dictionary, list):
        for item in dictionary:
          _ = self.remove_duplicate_targetname(item)
      return dictionary


      

#%5
class NINA:
  def __init__(self) -> None:
    self.Container = NINA_Container
    self.Action = NINA_Action
    self.Condition = NINA_Condition
    self.Trigger = NINA_Trigger
    self.Scriptmaker = NINA_ScriptMaker

"""
def remove_duplicate_targetname(dictionary):
  if isinstance(dictionary, dict):
    for key, value in dictionary.items():
      if key in ['Name', 'TargetName']:
        if isinstance(value, str):
          split_value = value.split("dup_")[0]
          dictionary[key] = split_value
        elif isinstance(value, (dict, list)):
          remove_duplicate_targetname(value)            
  elif isinstance(dictionary, list):
    for item in dictionary:
      if isinstance(item, (dict, list)):
        remove_duplicate_targetname(item)
  return dictionary
  """