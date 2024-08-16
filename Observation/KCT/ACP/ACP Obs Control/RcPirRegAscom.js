//
// Register RCOS Rotator for ASCOM
//
var P = new ActiveXObject("DriverHelper.Profile");
P.DeviceType = "Rotator";
P.Register("RCOS_AE.Rotator", "RCOS PIR Rotator");
