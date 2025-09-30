# This is a script for VisIt that can plot a shaping output from Axom's
# shaping driver for the heroic_roses shaping output.
#
# visit -cli -nowin -s visit_heroic_roses.py
#
import sys

# Make a plot and set the colors.
OpenDatabase("shaping.root")
AddPlot("FilledBoundary", "shaping_mesh_material")
fb = GetPlotOptions()
fb.SetMultiColor(0, (0,     0,   0, 255)) # black
fb.SetMultiColor(1, (40,   90, 160, 255)) # blue
fb.SetMultiColor(2, (7,   134,  57, 255)) # brightgreen
fb.SetMultiColor(3, (20,  115,  95, 255)) # darkgreen
fb.SetMultiColor(4, (101,  90,  78, 255)) # greenishbrown
fb.SetMultiColor(5, (238, 117, 136, 255)) # pink
fb.SetMultiColor(6, (140, 100, 126, 255)) # purplish
fb.SetMultiColor(7, (253,   4,  66, 255)) # red
fb.SetMultiColor(8, (189, 147,  25, 255)) # yellow
fb.SetMultiColor(9, (100, 100, 100, 255)) # free
fb.legendFlag = 0
SetPlotOptions(fb)
DrawPlots()

# Turn off some annotations.
annot = GetAnnotationAttributes()
annot.databaseInfoFlag = 0
SetAnnotationAttributes(annot)

# Save an image
SaveWindow()
sys.exit(0)
