# LWFA_codes
Sheet model codes for Laser Wakefield Acceleration of electrons. Includes codes for estimating the normalized quiver velocity for wavebreaking condition. Also estimates the maximum longitudinal electric field attainable before wavebreaking within an under-dense plasma.

NOTE!! All units are cgs. Displayed magnitude of electric-field is in statV per cm.

<electric_field.py> outputs a) electron density and b) longitudinal electric field profiles. Float <num_den> is the plsama-density, float < lamblas > is the laser wavelength in cm, float < v > is used as a factor for a_{0}, to either scale it up or down. Upon execution pyplot window will open, plotting long. E-field in red and e-density in blue vs. the location of a point in system in the lab-frame.

< histogram.py > outputs a) e-density and b) oscillator position vs. lab-frame position profiles. The matplotlib window has 2 panels. Top displays the oscillator postion vs. lab-frame position and the bottom panel displays e-density vs. lab-frame postion profile.

< histogram_rk4.py > does the same job as < histogram.py > only it employs rk-4 integration technique to estimate the maximum net-displacement after the laser passes an oscillator.
