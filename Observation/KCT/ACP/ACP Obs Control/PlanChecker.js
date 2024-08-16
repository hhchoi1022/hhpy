//tabs=4
//------------------------------------------------------------------------------
//
// Script:		 PlanChecker.js
// Author:		 Robert B. Denny <rdenny@dc3.com>
// Version: 	 8.1.1	  *** CHANGE VERSION BELOW! ***
// Requires:	 ACP 8.1 or later
//				 Windows Script 5.6 or later
//
// Environment:  This script is written to run on the desktop via Windows Script.
//
// Description:  This script reads and checks an ACP observing plan.
//
// Usage:		 1. Put a shortcut to this script on your Desktop.
//				 2. Drag and drop an ACP plan file to the icon for this script,
//					or double-click the icon and browse for the plan file.
//				 4. Notepad will appear with the report or a list of errors.
//
// Revision History:
//
// Date 	 Who	Description
// --------- ---	--------------------------------------------------
// 21-Oct-06 rbd	5.0.1 - Total overhaul, simplification. Now uses the 
//					Plan.wsc component, shared with AcquireImages.js, so the
//					same checking engine is used in both cases.
// 20-Nov-07 rbd    5.1.1 - (5.1HF3) #dush/dawnflats now pseudo-targets, other
//                  minor output cleanups.
// 27-Aug-12 rbd    7.0.1 (7.0) - GEM:783 Fix printout of sun-down angle form
//                  of #waituntil.
// 28-Aug-12 rbd    7.0.1 (7.0) - GEM:898 Add MPCCOMET elements lookup as opposed
//                  to 1-Line comet elements.
// 23-Oct-14 rbd    7.3.1 (7.3) - GEM:1230 Use new ASCOM utilities object.
// 30-Mar-16 rbd    8.1.1 (8.1) - GEM:1534 Include Count in cal sets with only
//                  one filter (bias typically!).
//----------------------------------------------------------------------------
//
var VERSION = "8.1.1";

//
// Enhance String with a trim() method
//
String.prototype.trim = function()
{
	return this.replace(/(^\s*)|(\s*$)/g, "");
};

// ----------
// dumpPlan() - Dump the plan to notepad
// ----------
//
function dumpPlan(fName, oStrm, Pl)
{
	var TT = new Array("Deep Sky Object","Tab-delimited J2000","MPC Asteroid Elements","MPC Comet Elements",
						"MPCORB Elements Lookup","MPCCOMET Elements Lookup", "NEOCP Ephemerides","Major Planet",
						"BIAS", "CHILL", "DARK", "CLOSE", "OPEN", "MANUAL", "DUSKFLATS", "DAWNFLATS", "NOWEATHER");

    oStrm.WriteLine("");
	oStrm.WriteLine("Plan: " + fName);
	oStrm.WriteLine("  Start with set number " + Pl.StartSetNum);
	if(Pl.Sets > 1) oStrm.WriteLine("  Repeat entire plan " + Pl.Sets + " times (" + Pl.Sets + " sets)");
	if(Pl.MinSetTime > 0.0) oStrm.WriteLine("  Minimum repeat set time " + Pl.MinSetTime + " hours");
	if(Pl.AFInterval > 0) oStrm.WriteLine("  Periodic AF every " + Pl.AFInterval + " min.");
	if(Pl.QuitAt !== null) oStrm.WriteLine("  Quit time " + Pl.QuitAt.toUTCString());
	if(Pl.Shutdown) oStrm.WriteLine("  Shut obs down at end");
	if(Pl.ChainPlan !== "") oStrm.WriteLine("  Chain to plan " + FSO.GetFileName(Pl.ChainPlan) + " at end");
	if(Pl.ChainScript !== "") oStrm.WriteLine("  Chain to script " + FSO.GetFileName(Pl.ChainScript) + " at end");

	for(var i in Pl.Targets)
	{
		var buf;
		var tgt = Pl.Targets[i];									// This Target
		oStrm.WriteLine("");
		oStrm.WriteLine("Target " + (parseInt(i) + 1) + ": " + tgt.Name +  // parseInt to remove leading 0 ???
					" (" + TT[tgt.Type] + ")");
		for(var j in tgt.Comments)									// Echo any comment lines
			oStrm.WriteLine("  " + tgt.Comments[j]);
		if(!tgt.NonImage)											// Skip for non-image pseudo-targets
		{
			if(tgt.Directory !== "")
				oStrm.WriteLine("  Images go into folder " + tgt.Directory);
			else if(tgt.CalFrame)
				oStrm.WriteLine("  Images go into normal cal frame folder");
			else
				oStrm.WriteLine("  Images go into normal images folder");
			if(tgt.Repeat > 1)
				oStrm.WriteLine("  Repeat this target " + tgt.Repeat + " times");
		}
		//
		// Order of waits is dependent on order in Acquireimages!
		//
		for(j in tgt.WaitUntil) {
		    var t = tgt.WaitUntil[j];
		    if(t < 0)
			    oStrm.WriteLine("  Wait in set " + (parseInt(j) + 1) + " until sun angle " + tgt.WaitUntil[j] + " deg.");
		    else
			    oStrm.WriteLine("  Wait in set " + (parseInt(j) + 1) + " until " + tgt.WaitUntil[j].toUTCString());
		}
		if(tgt.WaitFor > 0) 
			oStrm.WriteLine("  Wait for " + tgt.WaitFor + " sec. before starting");
		if(!tgt.NonImage)											// Skip all of this for non-immge pseudo-targets
		{
			if(!tgt.CalFrame)										// Skip these for cal frames
			{
				if(tgt.WaitInLimits > 0) 
					oStrm.WriteLine("  Wait until in obs limits or for a max of " + tgt.WaitInLimits + " min.");
				if(tgt.WaitZenDist !== null)
					oStrm.WriteLine("  Wait until zenith distance is " + tgt.WaitZenDist.ZenDist + 
									" or for max of " + tgt.WaitZenDist.TimeLimit + " min.");
				if(tgt.WaitAirMass !== null)
					oStrm.WriteLine("  Wait until air mass is " + tgt.WaitAirMass.AirMass + 
									" or for max of " + tgt.WaitAirMass.TimeLimit + " min.");
				if(tgt.AutoFocus) 
					oStrm.WriteLine("  Force auto-focus before imaging");
				if(tgt.Pointing) 
					oStrm.WriteLine("  Force pointing update before imaging");
				if(tgt.NoPointing) 
					oStrm.WriteLine("  Suppress pointing update before imaging");
				if(tgt.NoSolve) 
					oStrm.WriteLine("  Suppress plate solution and WCS addition on data image(s)");
				if(tgt.OrbTrack) 
					oStrm.WriteLine("  Enable orbital tracking (follow solar system body)");
				oStrm.WriteLine("  Dither	= " + tgt.Dither);
				oStrm.WriteLine("  PA		= " + tgt.PA);
			}
			if(tgt.SubFrame < 1.0) 
				oStrm.WriteLine("  Sub-frame to " + parseInt((tgt.SubFrame * 100)) + "%");
			else
				oStrm.WriteLine("  Full frame");
			if(tgt.ImageSets.length > 1)
			{
				buf = tgt.CalFrame ? "cal" : "filter";
				oStrm.WriteLine("  There are " + tgt.ImageSets.length + " " + buf + " sets:");
				for(j in tgt.ImageSets)
				{
					oStrm.WriteLine("	 Set #" + (parseInt(j) + 1) + ":");
					oStrm.WriteLine("	   Count	= " + tgt.ImageSets[j].Count);
					if(!tgt.CalFrame && tgt.ImageSets[j].FilterName !== "")
						oStrm.WriteLine("	   Filter	= " + tgt.ImageSets[j].FilterName);
					oStrm.WriteLine("	   Binning	= " + tgt.ImageSets[j].Binning);
					oStrm.WriteLine("	   Interval	= " + tgt.ImageSets[j].Interval);
					if(tgt.Stack && tgt.ImageSets[j].Count > 1) 
					{
						if(tgt.Align)
							oStrm.WriteLine("	   Stack these with auto-align of images");
						else
							oStrm.WriteLine("	   Stack without aligning images");
					}
				}
			}
			else
			{
			    oStrm.WriteLine("  Count		= " + tgt.ImageSets[0].Count);
				if(!tgt.CalFrame && tgt.ImageSets[0].FilterName !== "") 
					oStrm.WriteLine("  Filter  = " + tgt.ImageSets[0].FilterName);
				oStrm.WriteLine("  Binning	= " + tgt.ImageSets[0].Binning);
				oStrm.WriteLine("  Interval	= " + tgt.ImageSets[0].Interval);
			}
		}
		buf = tgt.TargetLine;
		if(buf.length > 60)
			buf = buf.substr(0, 60) + "...";
		oStrm.WriteLine("  Target	= \"" + buf + "\"");			// Finally echo the target line
	}
}

// =========
// FAKE UTIL
// =========

var Util = {
	Console: {
		PrintLine: function(msg) {
			oStrm.WriteLine(msg);
		}
	},
    
	HMS_Hours: function(str) {
		return UTL.HMSToHours(str);
	},
    
	DMS_Degrees: function(str) {
		return UTL.DMSToDegrees(str);
	}
};
    
// ===========
// MAIN DRIVER
// ===========

var ARGS = WScript.Arguments;
var planPath;
var sh = new ActiveXObject("WScript.Shell");
if(ARGS.Count() === 0) {											// File not dropped on icon, double-clicked
	try {
		var dlg = new ActiveXObject("MSComDlg.CommonDialog");		// May get "not licensed"
	} catch(e) {
		sh.PopUp("Drag and drop a plan file onto the icon.");
		WScript.Quit();
	}
    
	var ans = sh.Popup("You can drag and drop plans onto the icon. Do you want to browse for the plan now?",
						5, "Ease of use advisory", 36); 			// (4 = Question icon) + (32 = Yes-No buttons)
	if(ans != 6) WScript.Quit();									// Timed out or answered No
    
	dlg.DialogTitle = "Select the ACP Plan to check";				// Put up file browser
	dlg.Filter = "Plan Files (*.txt)|*.txt|All files (*.*)|*.*";
	dlg.InitDir = sh.SpecialFolders("MyDocuments") + "\\ACP Astronomy\\Plans";	// Start in user's default plans folder
	dlg.MaxFileSize = 255;
	dlg.ShowOpen(); 												// (display it now)
	if(dlg.FileName === "") WScript.Quit(); 						// File browser cancelled
	planPath = dlg.FileName;										// Use selected file
} else {
	planPath = ARGS(0); 											// Use drag-dropped file
}

var UTL;
try {
	UTL = new ActiveXObject("ASCOM.Utilities.Util");
} catch(ex) {
	WScript.Echo("The ASCOM Platform 6.3 or later must be installed to use Plan Checker.");
	WScript.Quit();
}
var FSO = new ActiveXObject("Scripting.FileSystemObject");
var oPath = FSO.GetSpecialFolder(2).Path + "\\PlanChecker.txt"; 	// User's temp folder
var oStrm = FSO.CreateTextFile(oPath, true);						// Create temp text file for report
oStrm.WriteLine("ACP Plan Checker " + VERSION + " at " + Date());
var Pl = new ActiveXObject("ACP.Plan");
Pl.Init(Util, null, "", false); 									// Init for offline use
if(Pl.Read(planPath))
	dumpPlan(FSO.GetFileName(planPath), oStrm, Pl);
oStrm.Close();
sh.Run("Notepad \"" + oPath + "\"");
