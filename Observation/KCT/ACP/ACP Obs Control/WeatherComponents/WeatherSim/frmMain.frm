VERSION 5.00
Begin VB.Form frmMain 
   BorderStyle     =   4  'Fixed ToolWindow
   Caption         =   "ACP Weather Simulator"
   ClientHeight    =   4890
   ClientLeft      =   45
   ClientTop       =   285
   ClientWidth     =   2940
   Icon            =   "frmMain.frx":0000
   LinkTopic       =   "Form1"
   LockControls    =   -1  'True
   MaxButton       =   0   'False
   ScaleHeight     =   4890
   ScaleWidth      =   2940
   StartUpPosition =   2  'CenterScreen
   Begin VB.TextBox txtSkyTemperature 
      Height          =   315
      Left            =   2220
      TabIndex        =   13
      Text            =   "-40"
      ToolTipText     =   "Deg. C"
      Top             =   2835
      Width           =   450
   End
   Begin VB.CheckBox chkSkyTemperature 
      Caption         =   "Sky Temperature:"
      Height          =   210
      Left            =   225
      TabIndex        =   14
      Top             =   2880
      Value           =   1  'Checked
      Width           =   1980
   End
   Begin VB.CheckBox chkHavePrecipitation 
      Height          =   300
      Left            =   2355
      TabIndex        =   11
      Top             =   2490
      Width           =   285
   End
   Begin VB.CheckBox chkWindVelocity 
      Caption         =   "WindVelocity:"
      Height          =   210
      Left            =   225
      TabIndex        =   20
      ToolTipText     =   "Uncheck to make unavailable"
      Top             =   4035
      Value           =   1  'Checked
      Width           =   1980
   End
   Begin VB.TextBox txtWindVelocity 
      Height          =   315
      Left            =   2220
      TabIndex        =   19
      Text            =   "8"
      ToolTipText     =   "Knots"
      Top             =   3990
      Width           =   450
   End
   Begin VB.CheckBox chkWindDirection 
      Caption         =   "WindDirection:"
      Height          =   210
      Left            =   225
      TabIndex        =   18
      ToolTipText     =   "Uncheck to make unavailable"
      Top             =   3630
      Value           =   1  'Checked
      Width           =   1980
   End
   Begin VB.TextBox txtWindDirection 
      Height          =   315
      Left            =   2220
      TabIndex        =   17
      Text            =   "170"
      ToolTipText     =   "Degrees FROM"
      Top             =   3600
      Width           =   450
   End
   Begin VB.CheckBox chkRelativeHumidity 
      Caption         =   "RelativeHumidity (frac):"
      Height          =   210
      Left            =   225
      TabIndex        =   16
      ToolTipText     =   "Uncheck to make unavailable"
      Top             =   3255
      Value           =   1  'Checked
      Width           =   1980
   End
   Begin VB.TextBox txtRelativeHumidity 
      Height          =   315
      Left            =   2220
      TabIndex        =   15
      Text            =   "0.4"
      Top             =   3210
      Width           =   450
   End
   Begin VB.CheckBox chkPrecipitation 
      Caption         =   "Precipitation detected"
      Height          =   210
      Left            =   240
      TabIndex        =   12
      ToolTipText     =   "Uncheck to make unavailable"
      Top             =   2520
      Value           =   1  'Checked
      Width           =   1980
   End
   Begin VB.CheckBox chkInsideTemperature 
      Caption         =   "InsideTemperature:"
      Height          =   210
      Left            =   240
      TabIndex        =   10
      ToolTipText     =   "Uncheck to make unavailable"
      Top             =   2130
      Value           =   1  'Checked
      Width           =   1980
   End
   Begin VB.TextBox txtInsideTemperature 
      Height          =   315
      Left            =   2235
      TabIndex        =   9
      Text            =   "25"
      ToolTipText     =   "Deg. C"
      Top             =   2085
      Width           =   450
   End
   Begin VB.CheckBox chkDewPoint 
      Caption         =   "DewPoint:"
      Height          =   210
      Left            =   240
      TabIndex        =   8
      ToolTipText     =   "Uncheck to make unavailable"
      Top             =   1740
      Value           =   1  'Checked
      Width           =   1980
   End
   Begin VB.TextBox txtDewPoint 
      Height          =   315
      Left            =   2235
      TabIndex        =   7
      Text            =   "13"
      ToolTipText     =   "Deg. C"
      Top             =   1695
      Width           =   450
   End
   Begin VB.CheckBox chkClouds 
      Caption         =   "Clouds (frac 0.0-1.0):"
      Height          =   210
      Left            =   240
      TabIndex        =   6
      ToolTipText     =   "Uncheck to make unavailable"
      Top             =   1365
      Value           =   1  'Checked
      Width           =   1980
   End
   Begin VB.TextBox txtClouds 
      Height          =   315
      Left            =   2235
      TabIndex        =   5
      Text            =   "0.1"
      Top             =   1320
      Width           =   450
   End
   Begin VB.CheckBox chkBarometricPressure 
      Caption         =   "BarometricPressure:"
      Height          =   210
      Left            =   240
      TabIndex        =   4
      ToolTipText     =   "Uncheck to make unavailable"
      Top             =   990
      Value           =   1  'Checked
      Width           =   1980
   End
   Begin VB.TextBox txtBarometricPressure 
      Height          =   315
      Left            =   2235
      TabIndex        =   3
      Text            =   "1011"
      ToolTipText     =   "Millibars"
      Top             =   945
      Width           =   450
   End
   Begin VB.TextBox txtAmbientTemperature 
      Height          =   315
      Left            =   2235
      TabIndex        =   1
      Text            =   "18"
      ToolTipText     =   "Deg. C"
      Top             =   570
      Width           =   450
   End
   Begin VB.CheckBox chkSafe 
      Caption         =   "Weather is safe"
      BeginProperty Font 
         Name            =   "MS Sans Serif"
         Size            =   8.25
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Height          =   270
      Left            =   240
      TabIndex        =   0
      Top             =   240
      Value           =   1  'Checked
      Width           =   1920
   End
   Begin VB.CheckBox chkAmbientTemperature 
      Caption         =   "AmbientTemperature:"
      Height          =   210
      Left            =   240
      TabIndex        =   2
      ToolTipText     =   "Uncheck to make unavailable"
      Top             =   615
      Value           =   1  'Checked
      Width           =   1980
   End
   Begin VB.Label label1 
      Caption         =   "ProgID:   WeatherSim.Weather"
      Height          =   315
      Left            =   225
      TabIndex        =   21
      Top             =   4455
      Width           =   2325
   End
End
Attribute VB_Name = "frmMain"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
Option Explicit

