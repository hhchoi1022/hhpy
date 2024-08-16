<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN" 
    "http://www.w3.org/TR/REC-html40/loose.dtd"> 
<%@ LANGUAGE="JSCRIPT" %>
<html>
<head>
    <script language="JScript" type="text/jscript" runat="server">
	    //
	    // This runs in an <iframe> since AJAX uploading ain't practical
	    //
	    // history
	    //
	    // 26-Sep-2006  rbd 5.0	Initial version from old upload.asp
	    // 27-Sep-2006	rbd	5.0	New iframe embedding
	    // 28-Sep-2006	rbd	5.0	Self-sizing parent iframe
	    // 24-Oct-2006	rbd 5.0 Major rework, subfolder support
	    // 05-Nov-2006	rbd	5.0 - Always allow admins upload permissions
		// 05-Dec-2006	rbd 5.0 Additional security, prevent direct access via URI by others
		// 10-Dec-2006	rbd	5.0 - Oops, too harsh.
		// 22-Nov-2010  rbd 6.0 - Add class to file upload element for Chrome.
	    //
		var rootPath = Request.ServerVariables("script_name");				// Path part of URI to this ASP page
		var rootBits = rootPath.split("/");									// bit[0] is leading "/"
		var rb1 = rootBits[1].toLowerCase();								// 1st path component
		if(rb1 != "files" && rb1 != "uploads") {							// Allow all to access these
		    if(rootBits.length < 3 && rootBits[2].toLowerCase() != User.Username.toLowerCase()) {		// Deny URI-accesss to others!
				Response.Clear();
				Response.Status = "401 Unauthorized";
				Response.Write("<h2>You are prohibited from accessing this area</h2>");
				Response.End();
			}
		}
		var keyName = rootBits[1] + "_subdir";								// Session subdir key name (e.g., images_subdir)
		var tgtPath = rootPath.substr(0, rootPath.lastIndexOf("/"));		// Strip file name & trailing /
		if(Session(keyName))												// Add target subdir
			tgtPath += "/" + Session(keyName);
    </script>
	<script src="/ac/iframecss.js" runat="client"></script>
	<!-- parent's stylesheets get inserted here -->
	<style>
		/* Now remove borders, margins, and padding 
		 * from parent's tiddler and viewer  classes
		 */
		.tiddler { border: none; margin: 0; padding: 0;}
		.viewer { border: none; margin: 0; padding: 0;}
	</style>
</head>
<body>
<div class="tiddler" id="tdiv"><div class="viewer" id="vdiv">
<% 
	if(!User.CanUpload && !User.IsAdministrator) {
		Response.Write("<p>No permission to upload</p>");
	} else {
		if(Request.ServerVariables("request_method").toLowerCase() != "post") 
		{
%>
Select a file to be uploaded to the <%= tgtPath %> directory:
<form style="margin-bottom: 4px;margin-top:2px;" name="auploadfrm" method="POST" 
					enctype="multipart/form-data" 
					action="<%= Request.ServerVariables("script_name") %>">
<input type="file" class="upload" name="File1" id="IDFile1" size="50"><br>
<input type="hidden" name="DestPath" value="<%= tgtPath %>">
</form>
<a class="button" style="margin-left:0px;" href="javascript: document.auploadfrm.submit()">Upload File</a>
<% 
		} else {
		    //
		    // POST: Move the file from it's temporary location to its final
		    // location.
		    //
		    if(Request.Form("FileName")) {									// Did upload the file
		        var fso = new ActiveXObject("Scripting.FileSystemObject");
		        var dest = Server.MapPath(Request.Form("DestPath")) + "\\" +
		                        fso.GetFileName(Request.Form("FilePath"));
		        if(dest.toLowerCase() != Request.Form("FilePath").toLowerCase()) 
		        {
		        	try {
			            fso.CopyFile(Request.Form("FilePath"), dest, true);	// Replace
			            fso.DeleteFile(Request.Form("FilePath"));
				        Response.Write("\nFile uploaded successfully. " + 
				        			"Click Refresh to update the file list. <a href='" + 
				        			Request.ServerVariables("script_name") + 
				        			"' target='_self'>Upload another.</a>");
			        } catch(ex) {
			        	Response.Write("\nFailed: " + (ex.message ? ex.message : ex) + 
			        			" <a href='" + Request.ServerVariables("script_name") + 
			        			"' target='_self'>Try again.</a>");
			        }
		        }
		        else
			        Response.Write("\nFile uploaded successfully. " + 
			        			"Click Refresh to update the file list. <a href='" + 
			        			Request.ServerVariables("script_name") + 
			        			"' target='_self'>Upload another.</a>");
		    }
		    else
		    	Response.Write("\nNo file. <a href='" + 
		    				Request.ServerVariables("script_name") + 
		    				"' target='_self'>Try again.</a>");
		}
	}
%>
</div></div>
<script language="Javascript" type="text/javascript" runat="client">
	var tdiv = document.getElementById("tdiv");
	document.body.style.backgroundColor = getStyle(tdiv, 'background-color');
	parent.document.getElementById("aupload").height = tdiv.offsetHeight + 12;
</script>
</body></html>
