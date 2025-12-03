# This is a script for VisIt that can plot a shaping output for the ice_cream example.
#
# visit -cli -nowin -s visit_plot_ice_cream.py
#
import sys

OpenDatabase("shaping.root")

# Plot 1 - The shaped materials
AddPlot("FilledBoundary", "shaping_mesh_material")
AddOperator("MultiresControl", 0)
fb = GetPlotOptions()
fb.SetMultiColor(0, (255, 254, 246, 255))  # air
fb.SetMultiColor(1, (170, 103, 89, 255))   # batter
fb.SetMultiColor(2, (255, 153, 204, 255))  # icecream
fb.SetMultiColor(3, (255, 0, 0, 255))      # sprinkles
fb.SetMultiColor(4, (255, 255, 255, 255))  # free
SetPlotOptions(fb)

# Draw the plots
DrawPlots()

# Set the view
v = GetView2D()
v.windowCoords = (0, 400, 0, 500)
v.viewportCoords = (0.25, 0.98, 0.05, 0.98)
SetView2D(v)

# Set the annotations.
annot = GetAnnotationAttributes()
annot.databaseInfoFlag = 0
annot.backgroundColor = (255,255,255,255)
annot.foregroundColor = (0,0,0,255)
annot.backgroundMode = annot.Solid
annot.userInfoFlag = 0
SetAnnotationAttributes(annot)

# Save an image
SaveWindow()

# Plot 2 - The sprinkle mesh
OpenDatabase("sprinkles.mesh")
AddPlot("Mesh", "main")
m2 = MeshAttributes()
m2.legendFlag = 0
m2.lineWidth = 1
m2.meshColor = (0, 0, 0, 255)
SetPlotOptions(m2)
AddOperator("MultiresControl", 0)
mrc2 = MultiresControlAttributes()
mrc2.resolution = 8
SetOperatorOptions(mrc2)
DrawPlots()

# Save an image
SaveWindow()

# Plot 3 - The shaping mesh domains.
DeleteAllPlots()
OpenDatabase("shaping.root")
AddPlot("Subset", "domains(shaping_mesh)")
DrawPlots()

# Save an image
SaveWindow()
sys.exit(0)
