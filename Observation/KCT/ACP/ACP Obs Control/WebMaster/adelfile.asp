<%@LANGUAGE="JSCRIPT"%><script language="JScript" type="text/jscript" runat="server">
	//
	// Delete a file. Form data must have del=filename
	// History
	//
	// 26-Sep-2006	rbd 5.0	For new web content, initial edit
	// 24-Oct-2006	rbd	5.0	Major overhaul, subdir support
	// 05-Dec-2006	rbd 5.0 Additional security, prevent direct access via URI by others
	// 10-Dec-2006	rbd	5.0 - Oops, too harsh.
	//
	var tgtPath = Request.ServerVariables("script_name");				// Path part of URI to this ASP page
	var tgtBits = tgtPath.split("/");									// bit[0] is leading "/"
	var tb1 = tgtBits[1].toLowerCase();								// 1st path component
	if(tb1 != "files" && tb1 != "uploads") {							// Allow all to access these
	    if(tgtBits.length < 3 && tgtBits[2].toLowerCase() != User.Username.toLowerCase()) {		// Deny URI-accesss to others!
			Response.Clear();
			Response.Status = "401 Unauthorized";
			Response.Write("<h2>You are prohibited from accessing this area</h2>");
			Response.End();
		}
	}
	var keyName = tgtBits[1] + "_subdir";								// Session subdir key name (e.g., images_subdir)
	tgtPath = tgtPath.substr(0, tgtPath.lastIndexOf("/"));				// Strip file name & trailing /
	if(Session(keyName))												// Add target subdir
		tgtPath += "/" + Session(keyName);
		
	Response.Clear();						// Get rid of leading " " (??)
	Response.ContentType = "text/plain";
	
	var reqFile = Request.Form("del");
	if(!reqFile) {
		Response.Write("No 'del' element in form data");
		Response.End();
	}
	
	var fn = Server.MapPath(tgtPath) + "\\" + reqFile;
	try {
		var fso = new ActiveXObject("Scripting.FileSystemObject");
		if(fso.FolderExists(fn))
			fso.DeleteFolder(fn);
		else
			fso.DeleteFile(fn);
		Response.Write("ok");
	} catch(ex) {
		Response.Write(ex.message ? ex.message : ex);
	}
</script>