<%@LANGUAGE="JSCRIPT"%><script language="JScript" type="text/jscript" runat="server">
	//
	// Ceate a folder. Form data must have dir=subdirname
	// History
	//
	// 24-Oct-2006	rbd	5.0 Initial edit
	// 05-Dec-2006	rbd 5.0 Additional security, prevent direct access via URI by others
	// 10-Dec-2006	rbd	5.0 - Oops, too harsh.
	// 12-Apr-2008	rbd 5.1HF5 - Copy support ASP files inito new folder!
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
	
	var reqDir = Request.Form("dir");
	if(!reqDir) {
		Response.Write("No 'dir' element in form data");
		Response.End();
	}
	var path = Server.MapPath(tgtPath) + "\\" + reqDir;
	try {
		var fso = new ActiveXObject("Scripting.FileSystemObject");
		fso.CreateFolder(path);
		// Populate with copies of ASP files in parent
		var e = new Enumerator(fso.GetFolder(fso.GetParentFolderName(path)).Files);
		for(; !e.atEnd(); e.moveNext()) {
			if(e.item().Name.search(/\.asp$/i) != -1) {
				fso.CopyFile(e.item().Path, path + "\\" + e.item().Name);
			}
		}
		Response.Write("ok");
	} catch(ex) {
		Response.Write(ex.message ? ex.message : ex);
	}

</script>