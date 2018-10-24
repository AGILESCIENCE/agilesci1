AGILE Science Tools
===================

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/8de598011bda45679591a15473a504a0)](https://app.codacy.com/app/bulgarelli/agilesci1?utm_source=github.com&utm_medium=referral&utm_content=AGILESCIENCE/agilesci1&utm_campaign=Badge_Grade_Dashboard)

AGILE is an ASI (Italian Space Agency) Small Scientific Mission dedicated to high-energy astrophysics which was successfully launched on April 23, 2007. The AGILE instrument is composed of three main detectors: a Tungsten-Silicon Tracker designed to detect and image photons in the 30 MeV-50 GeV energy band, an X-ray imager called Super-AGILE operating in the 18-60 keV energy band, and a Mini-Calorimeter that detects gamma-rays and charged particle energy deposits between 300 keV and 100 MeV. The instrument is surrounded by an anti-coincidence (AC) system.
The main purpose of the Silicon Tracker is to provide a compact imager for gamma-ray pho- tons of energy above 30 MeV. The Tracker converts the gamma-rays in heavy-Z material layers (245 μm of Tungsten, 0.07 X0), into an electron/positron pair of Minimum Ionizing Particles (MIP, corresponding to a most probable value of 110 keV ) in the detector, and records the electron/positron tracks through a combination of Silicon microstrip detectors and associated readout.
A GRID event is a collection of all the electron/positron interactions in the microstrip detector (called clusters) with additional information from energy deposits in the MCAL bars and the configuration of the AC plastic scintillators if present. A complete representation of the event topology allows the reconstruction of the incoming direction and energy of the gamma-ray.

The AGILE Science Tools are written in C++ and use the ROOT (CERN) and CFITSIO (NASA) libraries. All the tools can be run from the command line.
The Science Tools use the Parameter Interface Library (PIL), developed by the ISDC, for parameter input, allowing a variety of input methods:
1. the user can start the tool, which then asks for the values of some parameters interactively on the console;
2. the user can input all or some of the parameter values on the command line, in the same order as the provided parameter file;
3. the user can input all or some of the parameter values on the command line, in any order, writing parameter names explicitly as pname=value.
In all task was provided Some
methods the user can accept the defaults (stored in each task parameter file, if the already run defaults are the last values used) for parameter values not explicitly
by adding ’mode=h’ on the command line.
task input parameters are text files of “index” type, reporting in each row an input followed by the minimum and maximum times contained in the file and the file type.
file name
The parameter, when specified on the command line, must be preceeded by a @.

Details can be found here: https://agile.ssdc.asi.it/public/AGILE_SW_5.0_SourceCode/AGILE-IFC-OP-009_Build-21.pdf
