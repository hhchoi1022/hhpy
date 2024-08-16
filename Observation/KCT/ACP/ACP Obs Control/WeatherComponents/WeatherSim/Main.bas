Attribute VB_Name = "ModMain"
Option Explicit

Public g_iReferenceCount As Integer
Public g_bRunExecutable As Boolean

Sub Main()

    g_iReferenceCount = 0
    g_bRunExecutable = (App.StartMode = vbSModeStandalone)
    
    Load frmMain
    If g_bRunExecutable Then
        frmMain.WindowState = vbNormal
    Else
        frmMain.WindowState = vbMinimized
    End If
    frmMain.Show
    
End Sub
