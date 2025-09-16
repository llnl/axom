# This is a script for VisIt that can plot a shaping output for the circle example.
#
# visit -cli -nowin -s visit_plot_circle.py
#
import sys

OpenDatabase("shaping.root")

# Plot 1 - The shaped materials
AddPlot("FilledBoundary", "shaping_mesh_material")
AddOperator("MultiresControl", 0)
fb = GetPlotOptions()
fb.SetMultiColor(0, (244, 204, 153, 255))  # void
fb.SetMultiColor(1, (128, 128, 128, 255))  # steel
fb.SetMultiColor(2, (0, 0, 255, 255))      # free
SetPlotOptions(fb)
mrc = MultiresControlAttributes()
mrc.resolution = 10
SetOperatorOptions(mrc)

# Plot 2 - The shaping mesh
AddPlot("Mesh", "shaping_mesh")
m = MeshAttributes()
m.legendFlag = 0
m.lineWidth = 0
m.meshColor = (0, 0, 0, 255)
m.opacity = 0.258824
SetPlotOptions(m)

# Plot 3 - The outer circle
OpenDatabase("unit_circle.mesh")
AddPlot("Mesh", "main")
m2 = MeshAttributes()
m2.legendFlag = 0
m2.lineWidth = 6
m2.meshColor = (255, 102, 0, 255)
SetPlotOptions(m2)
AddOperator("Transform", 0)
t = TransformAttributes()
t.doScale = 1
t.scaleOrigin = (0, 0, 0)
t.scaleX = 5
t.scaleY = 5
SetOperatorOptions(t)
AddOperator("MultiresControl", 0)
mrc2 = MultiresControlAttributes()
mrc2.resolution = 15
SetOperatorOptions(mrc2)

# Plot 4 - The inner circle
AddPlot("Mesh", "main")
m3 = MeshAttributes()
m3.legendFlag = 0
m3.lineWidth = 6
m3.meshColor = (255, 0, 255, 255)
SetPlotOptions(m3)
AddOperator("Transform", 0)
t = TransformAttributes()
t.doScale = 1
t.scaleOrigin = (0, 0, 0)
t.scaleX = 2.5
t.scaleY = 2.5
SetOperatorOptions(t)
AddOperator("MultiresControl", 0)
mrc2 = MultiresControlAttributes()
mrc2.resolution = 15
SetOperatorOptions(mrc2)

# Draw the plots
DrawPlots()

# Turn off some annotations.
annot = GetAnnotationAttributes()
annot.databaseInfoFlag = 0
annot.backgroundColor = (255,255,255,255)
annot.foregroundColor = (0,0,0,255)
annot.backgroundMode = annot.Solid
annot.userInfoFlag = 0
SetAnnotationAttributes(annot)

# Save an image
SaveWindow()
sys.exit(0)
