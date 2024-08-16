//
// Save the existing author.html and index.asp files in a safe way,
// in anticipation of them being replaced by an update installation.
//
// The installer passes %WEBROOT% as the first argument, so this script
// will honor a non-standard web root.
//
// 07-Nov-2006	rbd		Initial edit
//
var fso = new ActiveXObject("Scripting.FileSystemObject");
var webPath = WScript.Arguments(0);				// In case changed DocRoot!
var upgPath = webPath + "\\UpgradeSavedContent";
if(fso.FileExists(webPath + "\\author.html"))	// Only for V5+ upgrades
{
	if(!fso.FolderExists(upgPath))
		fso.CreateFolder(upgPath);
	var i = 0;
	while(fso.FileExists(upgPath + "\\author-bk-" + i + ".html")) i += 1;
	fso.CopyFile(webPath + "\\author.html", upgPath + "\\author-bk-" + i + ".html");
	fso.CopyFile(webPath + "\\index.asp", upgPath + "\\index-bk-" + i + ".asp");
	fso = null;
}
