VERSION 1.0 CLASS
BEGIN
  MultiUse = -1  'True
  Persistable = 0  'NotPersistable
  DataBindingBehavior = 0  'vbNone
  DataSourceBehavior  = 0  'vbNone
  MTSTransactionMode  = 0  'NotAnMTSObject
END
Attribute VB_Name = "Weather"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = True
Attribute VB_PredeclaredId = False
Attribute VB_Exposed = True
Option Explicit

Private Const ALERT_TITLE As String = "WeatherSim.Weather"
Private Const MSG_START As String = "Property Get "
Private Const MSG_FIN As String = " is not implemented"

Private m_bConnected As Boolean

Private Sub Class_Initialize()

    m_bConnected = False
    If Not g_bRunExecutable Then g_iReferenceCount = g_iReferenceCount + 1
    
End Sub

Private Sub Class_Terminate()

    If Not g_bRunExecutable Then
        If g_iReferenceCount > 1 Then
            g_iReferenceCount = g_iReferenceCount - 1
        Else
            Unload frmMain
        End If
    End If

End Sub

Private Sub CheckConnected()
    If Not m_bConnected Then _
        Err.Raise vbObjectError + &H400, "WeatherSim.Weather", "The weather server is not connected"
End Sub

' =================
' PUBLIC PROPERTIES
' =================

Public Property Get AmbientTemperature() As Single

    CheckConnected
    If frmMain.chkAmbientTemperature.Value = 0 Then _
        Err.Raise vbObjectError + &H400, ALERT_TITLE, _
                MSG_START & "AmbientTemperature" & MSG_FIN
    AmbientTemperature = CSng(frmMain.txtAmbientTemperature.Text)
    
End Property

Public Property Get BarometricPressure() As Single
    
    CheckConnected
    If frmMain.chkBarometricPressure.Value = 0 Then _
        Err.Raise vbObjectError + &H400, ALERT_TITLE, _
                MSG_START & "BarometricPressure" & MSG_FIN
    BarometricPressure = CSng(frmMain.txtBarometricPressure.Text)

End Property

Public Property Get Clouds() As Single

    CheckConnected
    If frmMain.chkClouds.Value = 0 Then _
        Err.Raise vbObjectError + &H400, ALERT_TITLE, _
                MSG_START & "Clouds" & MSG_FIN
    Clouds = CSng(frmMain.txtClouds.Text)
    
End Property

Public Property Get Connected() As Boolean

    Connected = m_bConnected
    
End Property

Public Property Let Connected(ByVal newVal As Boolean)

    m_bConnected = newVal
    
End Property

Public Property Get DewPoint() As Single

    CheckConnected
    If frmMain.chkDewPoint.Value = 0 Then _
        Err.Raise vbObjectError + &H400, ALERT_TITLE, _
                MSG_START & "DewPoint" & MSG_FIN
    DewPoint = CSng(frmMain.txtDewPoint.Text)
    
End Property

Public Property Get InsideTemperature() As Single

    CheckConnected
    If frmMain.chkInsideTemperature.Value = 0 Then _
        Err.Raise vbObjectError + &H400, ALERT_TITLE, _
                MSG_START & "InsideTemperature" & MSG_FIN
    InsideTemperature = CSng(frmMain.txtInsideTemperature.Text)
    
End Property

Public Property Get Name() As String

    Name = "Weather Server Simulator"
    
End Property

Public Property Get Precipitation() As Boolean

    CheckConnected
    If frmMain.chkPrecipitation.Value = 0 Then _
        Err.Raise vbObjectError + &H400, ALERT_TITLE, _
                MSG_START & "Precipitation" & MSG_FIN
    Precipitation = (frmMain.chkHavePrecipitation.Value = 1)
    
End Property

Public Property Get RelativeHumidity() As Single

    CheckConnected
    If frmMain.chkRelativeHumidity.Value = 0 Then _
        Err.Raise vbObjectError + &H400, ALERT_TITLE, _
                MSG_START & "RelativeHumidity" & MSG_FIN
    RelativeHumidity = CSng(frmMain.txtRelativeHumidity.Text)
    
End Property

Public Property Get Safe() As Boolean

    CheckConnected
    Safe = (frmMain.chkSafe.Value = 1)
    
End Property

Public Property Get SkyTemperature() As Single

    CheckConnected
    If frmMain.chkSkyTemperature.Value = 0 Then _
        Err.Raise vbObjectError + &H400, ALERT_TITLE, _
                MSG_START & "SkyTemperature" & MSG_FIN
    SkyTemperature = CSng(frmMain.txtSkyTemperature.Text)
    
End Property

Public Property Get WindDirection() As Single

    CheckConnected
    If frmMain.chkWindDirection.Value = 0 Then _
        Err.Raise vbObjectError + &H400, ALERT_TITLE, _
                MSG_START & "WindDirection" & MSG_FIN
    WindDirection = CSng(frmMain.txtWindDirection.Text)
    
End Property

Public Property Get WindVelocity() As Single

    CheckConnected
    If frmMain.chkWindVelocity.Value = 0 Then _
        Err.Raise vbObjectError + &H400, ALERT_TITLE, _
                MSG_START & "WindVelocity" & MSG_FIN
    WindVelocity = CSng(frmMain.txtWindVelocity.Text)
    
End Property

' ==============
' PUBLIC METHODS
' ==============

Public Sub SetupDialog()

    MsgBox "There is nothing to set up.", (vbOKOnly + vbInformation), ALERT_TITLE
    
End Sub
