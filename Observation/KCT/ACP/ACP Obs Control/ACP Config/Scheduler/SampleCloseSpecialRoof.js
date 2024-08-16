//
// Template script for handling special-needs buildings like roll off
// roofs which hang low and can damage the telescope/mount. This script
// comes into play from Scheduler during Dawn Cal Frames, and oinly if 
// ACP's dome control option "Safe to slew anywhere with dome closed"
// is OFF (in other words, UNSAFE). Then it is this script's responsibility
// to get the roof or whatever closed safely.
//
// Bob Denny	27-Jul-2014
//
function main()
{
    Console.PrintLine("Park the scope then close the special roof...");
    // Depends on Dome setting "Close and park/home AFTER scope is parked
    // in ACP Preferences, Dome control, being ON (checked).
    Telescope.Park();                    // Park the scope, triggering the roof closure afterward
    Console.PrintLine("... completed.");
}