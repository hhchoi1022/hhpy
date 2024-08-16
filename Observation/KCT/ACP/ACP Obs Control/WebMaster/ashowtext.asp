<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN" 
    "http://www.w3.org/TR/REC-html40/loose.dtd"> 
<%@ LANGUAGE="VBSCRIPT" %>
<html>
<head>
    <script language="VBScript" type="text/vbscript" runat="server">
    '
    ' History
    '
    ' 26-Sep-2006   rbd 5.0		From old showtext.asp - should really recode in JS!!
    ' 20-Oct-2006	rbd	5.0		Nicer reporting when no file or file not found
    ' 05-Dec-2006	rbd 5.0 	Additional security, prevent direct access via URI by others
    ' 02-Sep-2014   rbd 7.2     UTF-8 for Cyrillic and special characters
    '
    </script>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
    <meta http-equiv="Content-Style-Type" content="text/css">
    <meta http-equiv="expires" content="Fri, 1 Jan 1990 00:00:00 GMT">
    <meta http-equiv="Pragma" content="no-cache">
    <meta http-equiv="cache-control" content="no-cache, no-store, must-revalidate">
    <style type="text/css">
        body { 
            background-color: #000030;
        } 

        h3  { font-family: Arial,Helvetica,Sans-Serif;
              margin-bottom: 0em;
              color: #e0e0e0; }
        pre { font-family: Lucida Console, Courier, Monospace;
              font-size: 9pt;
              margin-top: 0em; 
              color: #ffff80; }
    </style>
    <% tgt = Trim(Server.URLDecode(Request.ServerVariables("query_string"))) %>
    <title><%= tgt %></title>
</head>
<body>
<h3><%= tgt %></h3>
<hr align="left" width="600" size="1" color="Yellow" noshade>
<pre>
<script language="VBScript" type="text/vbscript" runat="server">
' Catch both missing querystring (otw use) and paths to others' data
deny = (tgt = "")	' Sorry-ass VBS, no short-circuit booleans
If Not deny Then
	sp = Request.ServerVariables("script_name")
	bits = Split(sp, "/")
	deny =  LCase(bits(2)) <> LCase(User.Username)
End If
If deny Then
	Response.Clear
	Response.Status = "401 Unauthorized"
	Response.Write "<h2>You are prohibited from accessing this area</h2>"
	Response.End
End If
fn = Server.MapPath(sp)     										' Our physical file path
fn = Left(fn, (InStr(fn, "ashowtext.asp") - 1))                  	' Our path with trailing /
fn = fn & tgt                                                   	' Path to text file for display
Set fso = CreateObject("Scripting.FileSystemObject")
On Error Resume Next
Set strm = fso.OpenTextFile(fn, 1)
If Err.Number <> 0 Then
	On Error GoTo 0
	Response.Write "There is no such file"
Else
	On Error GoTo 0
	Response.Write Server.HTMLEncode(strm.ReadAll)
	strm.Close
End If
Set fso = Nothing
</script></pre>
<hr align="left" width="600" size="1" color="Yellow" noshade>
</body>
</html>
