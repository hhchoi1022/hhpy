//------------------------------------------------------------------------------
//
// ShutdownObs.js - ACP Scheduler Observatory Shutdown Script
//
// Sample shutdown script for Scheduler 3.6. This script runs at the end of an
// observing night, not in response to a weather alert (which is handled by ACP).
// 
// Bob Denny        19-Jan-2009 Initial edit
// Bob Denny        22-Jan-2009 Scheduler now passes dawn cal frame parameter
// Bob Denny        23-Jan-2009 No more cal frames here! Greatly simplified.
// Bob Denny        26-May-2010 GEM:414 Make cooler warmup more robust per 
//                  AcquireSupport
// Bob Denny        29-Feb-2012 GEM: 809 Remove weather disconnect code. This is 
//                  now (as of 3.4) managed by Scheduler.
// Bob Denny        26-Dec-2012 (7.0) Disconnect rotator if needed, remove old
//                  RotControl.exe code.
// Bob Denny        01-Feb-2014 GEM:1093 Extend sanity timer in RunProgram for 
//                  long-running progs like AACmd for FlipFlat. Add power control
//                  sampels and IP*Works functions.
// Bob Denny        04-Apr-2014 Fix log file name (was startup)
// Bob Denny        03-Feb-2015 GEM:1173 Remove old power control, replace with 
//                  code using the ASCOM Digital Loggers switch driver. GEM:1156
//                  code for SynAccess NetBooter NP8 and friends.
// Bob Denny        18-Feb-2015 Typo in object ID and Connected=true for Digital
//                  Logger ASCOM driver.
// Bob Denny        27-Mar-2017 GEM:1533 Support multiple NetBooter power strips.
// Bob Denny        28-Dec-2018 Comm Center #11991 Add notation not to kill APCC
//------------------------------------------------------------------------------

var FSO;                                        // FileSystemObject used in several functions
var SWID = "ASCOM.DigitalLoggers.Switch";        // ID of power switch driver (e.g. DigitalLoggers)
var SWT = null;                                 // [sentinel] ASCOM Switch driver 

function main()
{
    var SEQTOK, buf;
    var doDawnCalFrames = false;
    
    FSO = new ActiveXObject("Scripting.FileSystemObject");

    Console.LogFile = Prefs.LocalUser.DefaultLogDir + "\\Scheduler\\ShutdownObs-" +
                Util.FormatVar(new Date().getVarDate(), "yyyymmdd") + "@" +
                Util.FormatVar(new Date().getVarDate(), "HhNnSs") + ".log";
    Console.Logging = true;
    
    //
    // Similar to shutdown logic in AcquireSupport.wsc
    //
    
    if(Telescope.Connected) {
        Console.PrintLine("Parking scope"); 
        Telescope.Park();                                               // Park the scope if possible, and/or close/park dome
        Console.PrintLine("OK");
        if(Telescope.Connected) {
            Console.PrintLine("Disconnecting Scope");
            Telescope.Connected = false;
        }
    }
    
    if(Util.RotatorConnected) {
        Console.PrintLine("Disconnecting rotator.");
        Util.RotatorConnected = false;
    }
     
    if(Util.CameraConnected) {
        Console.PrintLine("Shutting down imager. Please wait...");
        var z = Camera.TemperatureSetpoint;                             // Remember this, as MaxIm remembers
        Console.PrintLine("  (raising temperature to +5.0C... 20 min max)");
        var tPrev = -273.15;
        Console.PrintLine("  (cooler is now at " + Util.FormatVar(Camera.Temperature, "0.0") + "C");
        Camera.TemperatureSetpoint = 6.0;                               // Raise temperature to +6C
        for(var i = 0; i < 20; i++)                                     // Take 20 minutes max...
        {   
            var tNow = Camera.Temperature;
            if(tNow >= 3.0) break;                                      // Warmed, time to shut down cooler
            if((tNow - tPrev) < 3.0) break;                             // Warming rate < 0.05deg/sec, can shut down
            tPrev = tNow;
            Util.WaitForMilliseconds(60000);                            // Wait another minute
            Console.PrintLine("  (cooler is now at " +
                Util.FormatVar(Camera.Temperature, "0.0") + "C)");
        }
        Camera.CoolerOn = false;
        Camera.TemperatureSetpoint = z;                                 // Restore original setpoint for future
        Util.WaitForMilliseconds(2000);                                 // Give MaxIm a chance to shutdown cooler
        Camera.LinkEnabled = false;
        Util.CameraConnected = false;                                   // Disconnect it from ACP
        Console.PrintLine("OK, imager shutdown complete.");
    }
    
    Util.WaitForMilliseconds(1000);
    
    // ---------------------------------------------------------------------
    // NOTE: DO NOT KILL ASTRO-PHYSICS COMMAND CENTER (APCC) HERE. DOING SO 
    // WILL LEAVE ORPHAN VIRTUAL COM PORTS OPEN. INSTEAD, SET APCC FOR 
    // AUTO-SHUTDOWN AND IT WILL EXIT WHEN THE TELESCOPE IS DISCONNECTED
    // ABOVE.
    // ---------------------------------------------------------------------
    Console.PrintLine("Shutting down programs");
    if(!StopProcess("FocusMax.exe"))
        Console.PrintLine("**Failed to stop FocusMax");
    Util.WaitForMilliseconds(1000);
    if(!StopProcess("MaxIm_DL.exe"))
        Console.PrintLine("**Failed to stop MaxIm");

    // ======================================================
    // HERE IS WHERE YOU ADD CODE TO TURN POWER OFF AS NEEDED
    // ======================================================
    //

    // Sample for Digital Loggers Ethernet Power Controller using the
    // ASCOM.DigltalLoggers.Switch driver, available on the ASCOM Initiative
    // website at http://ascom-standards.org/Downloads/SwitchDrivers.htm
    // See the function definition in the Utility Functions section below
    // Set the SWID variable above for the ID of the specific switch driver.
    //
    // TurnOffPower(0, "Mount and USB hub");
    
    // Sample for SynAccess NetBooter NP8 and work-alikes. See the 
    // See the function definition in the Utility functions below.
    // Ref: http://synaccess-net.com/downloadDoc/NPStartup-B.pdf
    //
    // NBPowerOff("192.168.1.21", 0, "Mount and USB hub");
    //
    
    Console.PrintLine("Observatory shutdown complete");
    Console.Logging = false;
}

//
// Stop a process by executable NAME. WMI magic. Look in TaskManager
// to see what the executable NAMEs are.
//
var WMI = null;                                                         // Avoid creating lots of these
function StopProcess(exeName)
{
    var x = null;
    Console.PrintLine("Stopping " + exeName);
    try {
        // Magic WMI incantations!
        if(!WMI) WMI = new ActiveXObject("WbemScripting.SWbemLocator").ConnectServer(".","root/CIMV2");
        x = new Enumerator(WMI.ExecQuery("Select * From Win32_Process Where Name=\"" + exeName + "\"")).item();
        if(!x) {                                                        // May be 'undefined' as well as null
            Console.PrintLine("(" + exeName + " not running)");
            return true;                                                // This is a success, it's stopped
        }
        x.Terminate();
        Console.PrintLine("OK");
        return true;
    } catch(ex) {
        Console.PrintLine("WMI: " + (ex.message ? ex.message : ex));
        return false;
    }
}

// -----------------
// UTILITY FUNCTIONS
// -----------------

//
// Shell out to an executable. Probably handy for power control.
//  exePath: full path/name of executable to run
//  cmdLine: command line to give to that executable
//
// Waits till executable exits. Throws error if can't start executable,
// if the executable takes longer than 30 seconds to complete, or if 
// the executable exits with error status. The executable is run 
// minimized without getting focus ("background").
//
function RunProgram(exePath, cmdLine)
{
    var exeName = FSO.GetFileName(exePath);                             // Just the file name of the executable
    Console.PrintLine("Running " + exeName + "...");                    // Announce our intentions
    var tid = Util.ShellExec(exePath, cmdLine, 6);                      // Execute command (minimized, no focus, throws on errors)
    var i;
    for(i = 0; i < 30; i++) {                                           // Wait up to 10 sec
        if(!Util.IsTaskActive(tid)) break;
        Util.WaitForMilliseconds(1000);                                 // Wait 1 sec here
    }
    if(i >= 30)                                                         // Wait failed?
        throw "** " + exeName + " failed to exit in 30 sec.";
    if(Util.GetTaskExitStatus(tid) !== 0)                               // Exited with failure status?
        throw "** " + exeName + " exited with error status.";
    Console.PrintLine("...OK");
}

//
// Turn on a switch using the ASCOM switch driver
// Errors are fatal to the script, right?
//
function TurnOffPower(switchNumber, switchName)
{
    if (SWT === null) {
        SWT = new ActiveXObject(SWID);
        if (!SWT.Connected) SWT.Connected = true;
    }
    if(SWT.GetSwitch(switchNumber)) {
        SWT.SetSwitch(switchNumber, false);
        Console.PrintLine("...powered off " + switchName);
    } else {
        Console.PrintLine("..." + switchName + " already powered down");        
    }
}

//
// Turn on a switch for the SynAccess NetBooter NP8
// Ref: http://synaccess-net.com/downloadDoc/NPStartup-B.pdf
// NOTE: Assumes the username/password of all switches are the same
//
function NBPowerOff(IP, switchNumber, switchName)
{
    // --- user config ---
    var USER = "admin";
    var PASS = "admin";
    // -------------------
    var cmdline = "$A3 " + switchNumber + " 0";
    var http = new ActiveXObject("MSXML2.XMLHTTP");
    http.open("GET", "http://" + IP + "/cmd.cgi?" + encodeURI(cmdline), false, USER, PASS);
    http.send();
    if(http.status == 200) {
        Console.PrintLine("...powered off " + switchName);
        Console.PrintLine("   " + http.responseText);
    } else {
        throw "** Power off failed for " + switchName + " " + http.status + " " + http.statusText;
    }
}
