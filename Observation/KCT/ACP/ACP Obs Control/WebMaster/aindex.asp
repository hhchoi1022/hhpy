<%@LANGUAGE="JSCRIPT"%><script language="JScript" type="text/jscript" runat="server">
//
// History
//
// 25-Sep-2006  rbd 5.0	For new web content
// 27-Sep-2006	rbd	5.0	Overhaul, new iframe embedding
// 28-Sep-2006	rbd	5.0	Self-sizing upload iframe
// 29-Sep-2005	rbd	5.0 Custom CSS for the table, much nicer.
// 23-Oct-2006	rbd	5.0 Major rework, subfolder support. This is quite complex!
// 05-Dec-2006	rbd 5.0 Additional security, prevent direct access via URI by others
// 10-Dec-2006	rbd	5.0 - Oops, too harsh.
// 14-Sep-2007	rbd	5.1 - Hide folders whose names start with ~
// 12-Apr-2008	rbd	5.1 (HF5) - Suppress Wikiwording folder names
// 16-Jul-2018  rbd 8.2 - GEM:1601 Expunge FTP from system.
//
function makeViewLink(f)
{
    switch(fso.GetExtensionName(f).toLowerCase())
    {
        case "vbs":
        case "js":
        case "pl":
        case "pls":
        case "py":
        case "pys":
        case "txt":
        case "plan":
        case "log":
        	return "<html><a class='button' href='" + tgtPath + "/ashowtext.asp?" + 
        			f + "' title='View this file in a separate window' target='_blank'>View</a></html>";
        default:
        	return "";
    }
}

function showSize(s)
{
	if(s > 1073741824)
		return Util.FormatVar((s / 1073741824.0), "0.00") + "Gb";
    if(s > 1048576)
        return Util.FormatVar((s / 1048576.0), "0.00") + "Mb";
    else if(s > 1024)
        return Util.FormatVar((s / 1024.0), "0.00") + "Kb";
    else
        return s + 'b';												// Needed to sense a size column in SortableGrid
}

// -----
// SETUP
// -----

var fso = new ActiveXObject("Scripting.FileSystemObject");
//
// We keep the subfolder path from the root to the target folder
// in a Session variable whose name contains the root name.
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
if(!Session(keyName)) Session(keyName) = "";						// Make sure this exists for aux forms

// ---------------
// 'POST' RESPONSE
// ---------------
//
// All this does is change the current (sub) folder for all of the 
// file management ASP pages. 
if(Request.ServerVariables("request_method").toLowerCase() == "post")
{
	Session(keyName) = Request.Form("explore").Item;
	Response.End();
}
// END 'POST'

// -------------
// SETUP (CONT.)
// -------------

rootPath = rootPath.substr(0, rootPath.lastIndexOf("/"));			// Strip file name & trailing /
var physRootPath = Server.MapPath(rootPath);						// Physical path where ASP pages are

var tgtPath = rootPath;												// URL path to target folder
if(Session(keyName))
	tgtPath += "/" + Session(keyName);
var physTgtPath = Server.MapPath(tgtPath);							// Physical path to target folder

var canDelete = fso.FileExists(physRootPath + "\\adelfile.asp");	// Permission flags for this user/folder
var canUpload = (User.IsAdministrator || User.CanUpload) && 
		fso.FileExists(physRootPath + "\\aupload.asp");

// --------------
// 'GET' RESPONSE
// --------------

Response.Clear();													// Get rid of leading " " (??)
Response.ContentType = "text/plain";

Response.Write("<<slider '' 'My Documents Help' Help 'Info on file management'>>\n\n");	// Typical help button
if(canUpload) {														// Upload controls <iframe>, see aupload.asp
	//
	// This must be in an iframe because my AJAX library doesn't know 
	// how to post multipart/form-data
	//
	Response.Write('<html><table width="96%" style="border:none;margin:1px;background:transparent"><iframe id="aupload" src="' + 
		rootPath + '/aupload.asp" width="100%" height="0" frameborder="0" scrolling="no" allowtransparency style="border:none">' +
		'</iframe></table></html>');
}
if(tgtPath != rootPath)												// Parent Folder button
{
	Response.Write("<html><form style='margin:4px 0px 0px 0px' id='aindex_up'><input type='hidden' value='" + 
				Session(keyName).substr(0, Session(keyName).lastIndexOf("/")) +
				"' name='explore'></form></html><<PostForm 'aindex_up' '" + rootPath + 
				"/aindex.asp' 'Parent Folder' 'Up one folder level' '' false>> ");
}
if(canUpload)
	Response.Write("<html><a class='button' href='javascript:;' onClick='DC3.Lib.postMkdir(\"" + 
			rootPath + "/amkdir.asp\", this);'>Create Folder</a></html>");
Response.Write(" {{{" + tgtPath + "}}}");

// ------------
// File listing
// ------------

Response.Write("{{fileList{\n|!Name|!Type |!Date Modified |!Size |! |! |h\n");	// Write common table header (w/custom CSS)
var dirObj = fso.GetFolder(physTgtPath);								// FSO Folder object for current dir
if(dirObj.Files.Count > 0)
{
	//
	// Subfolder rows
	//
	var fColl = new Enumerator(dirObj.SubFolders);
	for(; !fColl.atEnd(); fColl.moveNext())							// Put table lines into array for pre-sort
	{
		var f = fColl.item();
		if(f.Name.substr(0, 1) == "~") continue;					// Hide ~ folders
		var dName = f.Name;
		if(dName.search(/(?:[A-Z][a-z]+){2,}/) != -1)				// WikiWord matcher!
			dName = "~" + dName;									// Suppress WikiWording
		Response.Write("|[" + dName + "]|Folder|" + 
				new Date(f.DateLastModified).toUTCString().replace(/^\w+, /, "") + "|--");
		// This is cheezy - make an Explore button
		// Don't need form name as form is only one in table cell (place)
		Response.Write("|<html><form><input type='hidden' value='" + 
						(Session(keyName) ? Session(keyName) + "/" + f.Name : f.Name) + 
						"' name='explore'></form></html><<PostForm '' '" + rootPath + 
						"/aindex.asp' 'Explore' 'Explore this folder' '' false>>");
		// If find adelfile.asp in folder, add delete column and buttons.
		if(canDelete) {
			// Don't need form name as form is only one in table cell (place)
			Response.Write("|<html><form><input type='hidden' value='" + f.Name + 
						"' name='del'></form></html><<PostForm '' '" + rootPath + 
						"/adelfile.asp' 'Delete' 'Delete this folder' 'Are you sure?' true>>");
		}
		Response.Write("|\n");
	}
	//
	// File rows (filtered)
	//
	fColl = new Enumerator(dirObj.Files);
	for(; !fColl.atEnd(); fColl.moveNext())
	{
		f = fColl.item();
		if(f.Name.substr(0, 1) == "~" || fso.GetExtensionName(f.Name).toLowerCase() == "asp")
			continue;												// Hide some files
		Response.Write("|[[" + f.Name + "|" + Server.UrlEncode(tgtPath + "/" + f.Name) + "]]|" + f.Type + "|" + 
				new Date(f.DateLastModified).toUTCString().replace(/^\w+, /, "") + "|" +
				showSize(f.Size) + "|" + makeViewLink(f.Name));
		// If find adelfile.asp in folder, add delete column and buttons.
		if(canDelete) {
			// Don't need form name as form is only one in table cell (place)
			Response.Write("|<html><form><input type='hidden' value='" + f.Name + 
						"' name='del'></form></html>" +	"<<PostForm '' '" + rootPath + 
						"/adelfile.asp' 'Delete' 'Delete this file' 'Are you sure?' true>>");
		}
		Response.Write("|\n");
	}
	Response.Write("}}}\n");										// Closing syntax for customm CSS
}
</script>
