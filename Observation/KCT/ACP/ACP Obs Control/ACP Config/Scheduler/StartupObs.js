//------------------------------------------------------------------------------
//
// StartupObs.js - ACP Scheduler Observatory Shutdown Script
//
// WARNING! THIS SCRIPT MUST RUN HARMLESSLY AND SUCCESSFULLY EVEN IF ALL OF
// THE POWER IS ALREADY ON, ALL PROGRAMS ARE ALREADY RUNNING, AND ALL DEVICES
// ARE ALREADY CONNECTED. IN OTHER WORDS, TEST THIS WITH A COMPLETELY ONLINE
// OBSERVATORY!!!
//
// WARNING! DO NOT OPEN THE DOME IN THIS SCRIPT! THE SCHEDULER OPENS IT IF NEEDED.
//
// Bob Denny        19-Jan-2009 Initial edit
// Bob Denny        23-Apr-2009 DO not open dome here
// Bob Denny        13-Jan-2011 Recognize rotator like shutdown, fix typos, MaxIm 5
// Bob Denny        21-Oct-2011 GEM:685 Do not touch tracking on startup, unsafe
//                  if low hanging roof (causes exception).
// Bob Denny        29-Feb-2012 GEM:809 Explicitly start Clarity as executable.
// Bob Denny        26-Dec-2012 (for ACP 7.0) No more RotControl.exe, connect
//                  rotator via API if configured.
// Bob Denny        30-Jan-2014 GEM:1093 Add RunProgram (for power control),
//                  IP*Works functions, and sample code. Uniquify log name.
// Bob Denny        02-Sep-2014 GEM:1198 Delay between turning on cooler and 
//                  setting the temperature (MaxIm 6.05).
// Bob Denny        03-Feb-2015 GEM:1173 Remove old power control, replace with 
//                  code using the ASCOM Digital Loggers switch driver. GEM:1156
//                  code for SynAccess NetBooter NP8 and friends.
// Bob Denny        18-Feb-2015 Typo in object ID and Connected=true for Digital
//                  Logger ASCOM driver.
// Bob Denny        19-May-2015 Throw errors on failures here so Scheduler goes
//                  to Operator Intervention.
// Bob Denny        27-Mar-2017 GEM:1533 Support multiple NetBooter power strips.
//------------------------------------------------------------------------------

var FSO;                                        // FileSystemObject used in several functions
var SWID = "ASCOM.DigitalLoggers.Switch";        // ID of power switch driver (e.g. DigitalLoggers)
var SWT = null;                                 // [sentinel] ASCOM Switch driver 

function main()
{
    var i;
    
    FSO = new ActiveXObject("Scripting.FileSystemObject");

    Console.LogFile = Prefs.LocalUser.DefaultLogDir + "\\Scheduler\\StartupObs-" +
                Util.FormatVar(new Date().getVarDate(), "yyyymmdd") + "@" +
                Util.FormatVar(new Date().getVarDate(), "HhNnSs") + ".log";
    Console.Logging = true;

    // =====================================================
    // HERE IS WHERE YOU ADD CODE TO TURN POWER ON AS NEEDED
    // =====================================================

    // Sample for Digital Loggers Ethernet Power Controller using the
    // ASCOM.DigltalLoggers.Switch driver, available on the ASCOM Initiative
    // website at http://ascom-standards.org/Downloads/SwitchDrivers.htm
    // See the function definition in the Utility Functions section below
    // Set the SWID variable above for the ID of the specific switch driver.
    //
    // TurnOnPower(0, "Mount and USB hub");
    
    // Sample for SynAccess NetBooter NP8 and work-alikes. See the 
    // See the function definition in the Utility functions below.
    // Ref: http://synaccess-net.com/downloadDoc/NPStartup-B.pdf
    //
    // NBPowerOn("192.168.1.21", 0, "Mount and USB hub");
    
    // ==================================================================
    // CHANGE PATHS TO PROGRAMS AS NEEDED AND ADD OTHER PROGRAMS YOU NEED
    // ==================================================================
    
    //Console.PrintLine("Starting support programs as needed");
    //if(!StartProgram("C:\\Program Files\\Diffraction Limited\\MaxIm DL V5\\MaxIm_DL.exe", 1)) { 
    //    Console.PrintLine("**Failed to start MaxIm");
    //   throw new Error(0x80040001, "==== Startup failed, cannot continue");
    //}
    //if(!StartProgram("C:\\Program Files\\FocusMax\\FocusMax.exe", 1)) {
    //    Console.PrintLine("**Failed to start FocusMax");
    //    throw new Error(0x80040001, "==== Startup failed, cannot continue");
    //}
    // TSX connection (Gu Lim)
    Console.PrintLine("Starting TheSkyX...");
    if(!StartProgram("C:\\Program Files (x86)\\Software Bisque\\TheSkyX Professional\\TheSkyX.exe", 1)){ 
            Console.PrintLine("**Failed to start TheSkyX");
            throw new Error(0x80040001, "==== Startup failed, cannot continue");
            Util.WaitForMilliseconds(3000);
        }


    // ==========================================
    // UNCOMMENT IF YOU USE BOLTWOOD CLOUD SENSOR
    // ==========================================
    
//     if(!StartProgram("C:\\Program Files\\Boltwood Systems\\Clarity II\\ClarityII.exe", 1)) {
//         Console.PrintLine("**Failed to start Clarity");
//         throw new Error(0x80040001, "==== Startup failed, cannot continue");
//     }
//     var C;
//     try {
//         if(!Util.Weather.Available) {
//             Console.PrintLine("Connect to Clarity and wait for data");
//             C = new ActiveXObject("ClarityII.CloudSensorII");           // Start and wait for Clarity data
//             for(i = 0; i < 20; i++)                                     // Wait for up to 5  minutes
//             {
//                 if(C.DataReady()) {
//                     Util.WeatherConnected = true;                       // Connect ACP to Clarity
//                     break;
//                 }
//                 Util.WaitForMilliseconds(15000);
//             }
//             if(i >= 20)
//                 Console.PrintLine("**Clarity didn't connect for 5 minutes");
//         } else {
//             Console.PrintLine("Weather already connected");
//         }
//     } catch(ex) {
//         Console.PrintLine("Clarity trouble:");
//         Console.PrintLine("  " + ex.message);
//     } finally {
//         C = null;                                                       // Assure our Clarity object released
//     }
// 
//     Console.PrintLine("OK");
    
    try {
        if(!Telescope.Connected) {
            Util.WaitForMilliseconds(5000);
            Console.PrintLine("Connect ACP to the telescope, will auto-home if needed");
            Telescope.Connected = true;                                 // (Requires driver "home on connect")
            Console.PrintLine("OK");
            Console.PrintLine("Parking scope...");// TSX Park (Gu Lim)
            Telescope.Park(); // TSX Park (Gu Lim)
            Console.PrintLine("Parked")
        } else {
            Console.PrintLine("Telescope already connected");
        }
    } catch(ex) {
        Console.PrintLine("**Failed to connect to the Telescope:");
        Console.PrintLine("  " + ex.description);
        throw new Error(0x80040001, "==== Startup failed, cannot continue");
    }

    if(Rotator.Configured)                                              // If there is a rotator, get it going
    {
        try {
            if(!Util.RotatorConnected){
                Console.PrintLine("Connecting to rotator.");
                Util.RotatorConnected = true;
            }
        } catch(ex) {
            Console.PrintLine("**Failed to connect to rotator:");
            Console.PrintLine("  " + ex.message);
            throw new Error(0x80040001, "==== Startup failed, cannot continue");
        }
    }

    try {                                                               //  Connect the camera
        if(!Util.CameraConnected) {
            Util.WaitForMilliseconds(5000);
            Console.PrintLine("Connect ACP to the camera");
            Util.CameraConnected = true;
            Console.PrintLine("OK");
        } else {
            Console.PrintLine("Camera already connected");
        }
    } catch(ex) {
        Console.PrintLine("**Failed to connect to MaxIm and camera:");
        Console.PrintLine("  " + ex.description);
        throw new Error(0x80040001, "==== Startup failed, cannot continue");
    }
    
    // ===========================================
    // IF FOR SCHEDULER DO NOT OPEN THE DOME HERE!
    // ===========================================
        
    Console.PrintLine("Chill cooler to -10 by Gu Lim");
    Camera.CoolerOn = true;
    Util.WaitForMilliseconds(5000);                                     // For quirk in MaxIm 6.05
    Camera.TemperatureSetpoint = -10;                                   // Chill the cooler

    // ==========================
    // ADD ANYTHING ELSE YOU NEED
    // ==========================
    
    Console.PrintLine("Startup complete");
    Console.Logging = false;
}

//
// Start program by executable path, only if not already running.
// Waits 15 sec after starting for prog to initialize.
//
function StartProgram(exePath, windowState)
{
    try {
        var f = FSO.GetFile(exePath);
        if(IsProcessRunning(f.Name))                                    // Proc name is full file name
            return true;                                                // Already running
        if(IsProcessRunning(f.ShortName))                               // COM servers can have proc name of file short name
            return true;
        Console.PrintLine("Starting " + f.Name);
        Util.ShellExec(exePath, "", windowState);
        Util.WaitForMilliseconds(15000);                                // Wait for prog to initialize
        Console.PrintLine("OK");
        return true;
    } catch(ex) {
        Console.PrintLine("Exec: " + ex.message);
        return false;
    }
}

//
// Test if program is running (by exe name)
//
var WMI = null;                                                         // Avoid creating lots of these
function IsProcessRunning(exeName)
{
    try {
        // Magic WMI incantations!
        if(!WMI) WMI = new ActiveXObject("WbemScripting.SWbemLocator").ConnectServer(".","root/CIMV2");
    	var x = new Enumerator(WMI.ExecQuery("Select * From Win32_Process Where Name=\"" + exeName + "\"")).item().Name;
    	return true;
    } catch(e) {
    	return false;
    }
}

// -----------------
// UTILITY FUNCTIONS
// -----------------

//
// Shell out to an executable.
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
function TurnOnPower(switchNumber, switchName)
{
    if (SWT === null) {
        SWT = new ActiveXObject(SWID);
        if (!SWT.Connected) SWT.Connected = true;
    }
    if(!SWT.GetSwitch(switchNumber)) {
        SWT.SetSwitch(switchNumber, true);
        Console.PrintLine("...powered on " + switchName);
    } else {
        Console.PrintLine("..." + switchName + " already powered up");        
    }
}

//
// Turn on a switch for the SynAccess NetBooter NP8
// Ref: http://synaccess-net.com/downloadDoc/NPStartup-B.pdf
// NOTE: Assumes the username/password of all switches are the same
//
function NBPowerOn(IP, switchNumber, switchName)
{
    // --- user config ---
    var USER = "admin";
    var PASS = "admin";
    // -------------------
    var cmdline = "$A3 " + switchNumber + " 1";
    var http = new ActiveXObject("MSXML2.XMLHTTP");
    http.open("GET", "http://" + IP + "/cmd.cgi?" + encodeURI(cmdline), false, USER, PASS);
    http.send();
    if(http.status == 200) {
        Console.PrintLine("...powered on " + switchName);
        Console.PrintLine("   " + http.responseText);
    } else {
        throw "** Power on failed for " + switchName + " " + http.status + " " + http.statusText;
    }
}
