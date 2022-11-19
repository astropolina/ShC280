# ShC280
Sequence for the processing:
1. Reduction: ShC280_reduction.py (using ccdproc package).
2. Running Astrometry.net and SExtractor: ShC280_astrometry_sextractor.py. 
The script uses:
default.param (a text file with SExtractor output parameters).
3. Relative photometry: ShC280_relative_photometry.py. 
The script uses:
Coords_Ephemeris.dat (a text file with equatorial coordinates, orbital periods and initial epoch of the objects), 
Dates_Objects.txt (a text file with dates of observations and observed objects, including number of frames in Green and Red bands),
files in Checks folder (text files with the stars inside field of the objects, which can be used as reference stars).

Examples of the results are presented in relPhotResults folder. 
