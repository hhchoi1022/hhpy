            ACP Observatory Control Software V8.2.1
                        Getting Started

          ============================================
          PLEASE DON'T START ACP WITHOUT FIRST OPENING
          THE ONLINE USER'S GUIDE FROM THE  START MENU
          AND FOLLOWING THOSE DIRECTIONS. THANKS!!!
          ============================================

Thanks for your interest in ACP Observatory Control Software.
Requirements and installation info is below in this document.
Once installed, the "getting started" info in the User's Guide
provides a roadmap to success. 

You have 60 days from the time you install ACP to decide if you 
want to purchase a license. During this evaluation period, all 
features, including internet support, are active.

PLEASE PLEASE PLEASE, DON'T START ACP YET! Take the time to
really go through the setup process. I get "ACP is broken" and
"help me" calls resulting from peoples' unwillingness to even
BEGIN to work through the given setup process. Like any new tool
you get, you need to invest a bit of time learning how to
properly use it. 

Here's how to get started:

Requirements:
------------

* PC running Windows 7-10, either 32-bit or 64-bit.

* MaxIm DL version 5.25 or 6.19 or later (contact Diffraction 
  Limited at http://www.cyanogen.com/ or (613) 225-2732 for 
  upgrade info, if needed).

You may test using simulators for the telescope, camera, and
(optionally) guider and a simulation mode in the scripts that
generates synthetic starfields from the PinPoint star catalog.
These synthetic starfields will contain the correct stars for the
coordinates at which the image is taken, and can be plate-solved
by PinPoint. The images also have mechanical and random pointing
errors so you can test ACP's built-in pointing corrector.

To use the script simulation, open the User's Guide from the
Start menu. DO NOT START ACP NOW! Look in the User's Guide for
the page "Using the simulators for Test and Practice" and review
that information. DO NOT START ACP NOW though!! (are you sick of
this yet?) If you decide to use simulation, you don't need a
telescope and dark skies. But you can simulate with your real
telescope and/or camera and still test during daylight or bad
weather.

Otherwise:

* ASCOM-compatible go-to telescope mount. Mechanical problems such
  as excessive backlash, stiction or poor polar alignment will
  greatly reduce your liklihood of success with ACP!

* CCD Camera that can be controlled by MaxIm DL

* CCD Autoguider - You must be able to start autoguiding without 
  fiddling with the camera or guider manually.

* Filter Wheel (optional). Often part of the camera package.

* ASCOM-compatible focuser (optional). If you have one, you can
  use ACP's autofocus support. This requires either the FocusMax 
  program (free or paid version) or the PlaneWave PWI program.
  
* ASCOM-compatible dome or roof controller (optional). Typically 
  this is the most troublesome part of an observatory if present.

* ASCOM-compatible rotator (optional)

* Compatible weather station equipment and software (see ACP Help
  for details)

Installation:
------------

If you're reading this you have already installed the basic 
required software (latest ASCOM Platform, MaxIm DL 5/6).

You'll also need to download and install the Guide Star Catalog 
1.1. It is available at http://gsc.dc3.com/. Follow the 
instructions in the read-me there. Then use the Catalog Checker
(Start menu under ACP) to check its integrity.

Note that later you may want to download the USNO A2.0 for plate 
solving if you have a narrow field (15 arc min or less). While 
newer catalogs such as UCAC4 are recommended for science, the 
A2.0 is still the bet for general plate solving for operational
needs. Go to the DC-3 Dreams Communication Center at

  http://forums.dc3.com/

Open the Company and Product Information section, then open
Downloading and Using PinPoint Reference Catalogs

http://forums.dc3.com/showthread.php?4694

This will tell you how to get the various catalogs and install 
them.

Setup Process:
-------------

DO NOT START ACP YET! Open the ACP User's Guide from the Start
menu. Follow the Getting Started instructions in that document.
Really. Don't try to out-smart it :-)

(22-Jan-2019)

