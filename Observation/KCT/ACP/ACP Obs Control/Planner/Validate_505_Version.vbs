'
' Alert if SN version is less than 5.0.5
'
Sub Main()
    Dim FSO, SNPath, SNVer, VerBits, VerOK, IsPlus, Flavor

    SNPath = Session.Property("CustomActionData")           ' [SNPATH] or [SNPPATH] (ends in \)
    IsPlus = (InStr(SNPath, "Plus") > 0)
    SNPath = SNPath & "starrynight.exe"                     ' Full pathname of EXE
    '
    ' WARNING! SNPATH WILL BE NON_BLANK BUT WRONG IS THIS VERSION OF SN IS NOT INSTALLED
    ' THE INSTALLER CONDITION SHOULD PREVENT THIS SCRIPT FROM BEING RUN, BUT...
    '
    Set FSO = CreateObject("Scripting.FileSystemObject")
    On Error Resume Next                                    ' Catch bad SNPath
    SNVer = FSO.GetFileVersion(SNPath)
    If Err.Number <> 0 Then 
        Set FSO = Nothing
        Exit Sub                                            ' GET OUT NOW
    End If
    Set FSO = Nothing
    VerOK = False                                           ' Assume failure
    If SNVer <> "" Then                                     ' Safety check
        VerBits = Split(SNVer, ".")                         ' Get the version bits
        If UBound(VerBits) >= 2 Then                        ' Another safety check
            If VerBits(0) >= 6 Then
                VerOK = True                                ' Version 6.x.x.x or later
            ElseIf VerBits(1) >= 1 Then
                VerOK = True                                ' Version 5.1.x.x or later
            ElseIf VerBits(2) >= 5 Then
                VerOK = True                                ' Version 5.0.5.x or later
            End If
        End If
    End If

    If Not VerOK Then
        If IsPlus Then
            Flavor = "Pro Plus"
        Else
            Flavor = "Pro"
        End If
        MsgBox "Your Starry Night " & Flavor & " must be updated to at least version 5.0.5." & VbCrLf & _
            "Use the Start menu, Starry Night " & Flavor & ", Check for Starry Night " & Flavor & " updates.", _
                (vbOKOnly + vbInformation), "ACP Planner Setup"
    End If

End Sub

Call Main()