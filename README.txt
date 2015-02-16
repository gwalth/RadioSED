Displays galaxy IR templates with ALMA and VLA bands for a given redshift.
CO emission lines also are plotted, with green representing that they are in a
band and red representing it is not.  The bold dashed lines show the ALMA
sensitivities (as of 2013).  The dash-dotted line represents a 10 mJy
source, the dashed line represents a 0.1 mJy (which is roughly detectable in
2 hours).

Galaxy templates come from Rieke 2009 average templates:
http://mingus.as.arizona.edu/~bjw/ir_templates/

Dependencies:
   numpy, matplotlib

Usage:
   #           redshift
   RadioSED.py 1.2
