@echo off
rem Run elevated to restore the registration of the WSCs
rem Part of the ACP automatic self-heal feature
rem Parameter must be the ACP install folder, and ACP
rem can't pass it quoted so we get the entire param string
rem and quote it here for the 'CD'
CD "%*"
regsvr32 /s AcquireSupport.wsc
regsvr32 /s Plan.wsc 
regsvr32 /s SequenceToken.wsc
regsvr32 /s WeatherComponents\AAGCloudWatcher.wsc
regsvr32 /s WeatherComponents\AuroraCloudSensor.wsc
regsvr32 /s WeatherComponents\BoltwoodFile.wsc
regsvr32 /s WeatherComponents\BoltwoodServer.wsc
regsvr32 /s WeatherComponents\CloudSensor-II.wsc
regsvr32 /s WeatherComponents\SROWeather.wsc

