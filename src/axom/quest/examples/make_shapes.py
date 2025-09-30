import math

def lerp_point(p0, p1, t):
    """
    Linearly interpolate 2 points.
    """
    return ((1. - t)*p0[0] + t*p1[0], (1. - t)*p0[1] + t*p1[1])

class mfem:
   """
   This class represents 1D contours that can be written to MFEM format.
   """
   def __init__(self):
      self.elemCount = 0
      self.vertCount = 0
      self.elemIds = []
      self.points = []
      self.points_ho = []

   def add_edge2(self, elem, p0, p1):
      p2 = lerp_point(p0, p1, 1./3.)
      p3 = lerp_point(p0, p1, 2./3.)
      self.add_edge4(elem, (p0, p2, p3, p1))

   def add_edge4(self, elem, p):
      self.elemIds.append(elem)
      self.points.append(p[0])
      self.points.append(p[3])
      self.points_ho.append(p[2])
      self.points_ho.append(p[1])
      self.elemCount = self.elemCount + 1
      self.vertCount = self.vertCount + 2

   def write(self, f):
      f.write("MFEM NURBS mesh v1.0\ndimension\n1\n\n")
      f.write("elements\n")
      f.write(f"{self.elemCount}\n")
      for i in range(self.elemCount):
          f.write(f"{self.elemIds[i]} 1 {2*i} {2*i+1}\n")
      f.write("\nboundary\n0\n\n")
      f.write(f"edges\n{self.elemCount}\n")
      for i in range(self.elemCount):
          f.write(f"{i} 0 1\n")
      f.write(f"\nvertices\n{self.vertCount}\n\n")
      f.write(f"knotvectors\n{self.elemCount}\n")
      for i in range(self.elemCount):
          f.write("3 4 0 0 0 0 1 1 1 1\n")
      f.write("\nweights\n")
      for i in range(self.vertCount):
          f.write("1 1\n")
      f.write("\nFiniteElementSpace\n")
      f.write("FiniteElementCollection: NURBS\n")
      f.write("VDim: 2\n")
      f.write("Ordering: 1\n")
      i = 0
      while i < self.vertCount:
          f.write(f"{self.points[i][0]} {self.points[i][1]} {self.points[(i+1)][0]} {self.points[(i+1)][1]}\n")
          i = i + 2
      i = 0
      while i < self.vertCount:
          f.write(f"{self.points_ho[i][0]} {self.points_ho[i][1]} {self.points_ho[(i+1)][0]} {self.points_ho[(i+1)][1]}\n")
          i = i + 2


def add_polygon(data, elem, pts):
    """
    Add a polygonal shape to the MFEM dataset.
    """
    npts = len(pts)
    if npts < 3:
       raise "Not enough points"
    else:
       for i in range(npts):
           nexti = (i + 1) % npts
           data.add_edge2(elem, pts[i], pts[nexti])

def add_circle(data, elem, origin, radius, reverse = False):
    """
    Add a circle to the MFEM dataset using multiple contours. If reverse is
    true then the contours will be added in the opposite orientation, indicating
    that we want to remove the region from the shape.
    """
    def reverse_points(pts):
        return (pts[3], pts[2], pts[1], pts[0])

    k = 0.5522847498
    p0 = (radius * 1., radius * 0.)
    p1 = (radius * 1., radius * k)
    p2 = (radius * k,  radius * 1.)
    p3 = (radius * 0., radius * 1.)

    segments = ( # segment 0
                ((origin[0] + p0[0], origin[1] + p0[1]),
                 (origin[0] + p1[0], origin[1] + p1[1]),
                 (origin[0] + p2[0], origin[1] + p2[1]),
                 (origin[0] + p3[0], origin[1] + p3[1])),
                 # segment 1
                ((origin[0] - p3[0], origin[1] + p3[1]),
                 (origin[0] - p2[0], origin[1] + p2[1]),
                 (origin[0] - p1[0], origin[1] + p1[1]),
                 (origin[0] - p0[0], origin[1] + p0[1])),
                 # segment 2
                ((origin[0] - p0[0], origin[1] - p0[1]),
                 (origin[0] - p1[0], origin[1] - p1[1]),
                 (origin[0] - p2[0], origin[1] - p2[1]),
                 (origin[0] - p3[0], origin[1] - p3[1])),
                 # segment 3
                ((origin[0] + p3[0], origin[1] - p3[1]),
                 (origin[0] + p2[0], origin[1] - p2[1]),
                 (origin[0] + p1[0], origin[1] - p1[1]),
                 (origin[0] + p0[0], origin[1] - p0[1])))

    if reverse:
        data.add_edge4(elem, reverse_points(segments[3]))
        data.add_edge4(elem, reverse_points(segments[2]))
        data.add_edge4(elem, reverse_points(segments[1]))
        data.add_edge4(elem, reverse_points(segments[0]))
    else:
        data.add_edge4(elem, segments[0])
        data.add_edge4(elem, segments[1])
        data.add_edge4(elem, segments[2])
        data.add_edge4(elem, segments[3])

def gear():
    """
    Generate an MFEM shaping dataset that contains multiple contours that are
    used for adding and subtracting from the shape.
    """
    m = mfem()
    add_circle(m, 1, (0., 0.), 1.)
    add_circle(m, 1, (0., 0.), 0.5, reverse=True)
    # Make smaller cutouts
    n = 8
    for i in range(n):
        a = i * (2 * math.pi) / n
        r = 0.75
        add_circle(m, 1, (r*math.cos(a), r*math.sin(a)), 0.1, reverse=True)
    # Make some teeth
    n = 45
    for i in range(n):
        a = i * (2 * math.pi) / n
        a1 = a - (2 * math.pi) / 90
        a2 = a
        a3 = a + (2 * math.pi) / 90
        r1 = 0.99
        r2 = 1.1
        p0 = (r1*math.cos(a3), r1*math.sin(a3))
        p1 = (r1*math.cos(a1), r1*math.sin(a1))
        p2 = (r2*math.cos(a2), r2*math.sin(a2))
        add_polygon(m, 1, (p0, p1, p2))

    with open("gear.mesh", "wt") as f:
        m.write(f)
        f.close()
    with open("gear.yaml", "wt") as s:
        s.write("dimensions: 2\n")
        s.write("shapes:\n")
        s.write(f" - name: shape1\n")
        s.write(f"   material: gear\n")
        s.write("   geometry:\n")
        s.write("     format: mfem\n")
        s.write(f"     path: gear.mesh\n")
        s.close()

def open_polygons():
    """
    Generate an MFEM shaping dataset that contains open polygons.
    """
    def single_polygon(data, elem, extents, t):
       c0 = (extents[0], extents[2])
       c1 = (extents[1], extents[2])
       c2 = (extents[1], extents[3])
       c3 = (extents[0], extents[3])
       
       p0 = lerp_point(c0, c1, 0.5)
       p1 = c1
       p2 = c2
       p3 = lerp_point(c2, c3, 0.5)
       data.add_edge4(elem, (p0, p1, p2, p3))

       m = lerp_point(c0, c3, 0.5)
       m0 = lerp_point(m, c3, t)
       q0 = lerp_point(p3, m0, 0.)
       q1 = lerp_point(p3, m0, 1./3.)
       q2 = lerp_point(p3, m0, 2./3.)
       q3 = lerp_point(p3, m0, 1.)
       data.add_edge4(elem, (q0, q1, q2, q3))

       m1 = lerp_point(m, c0, t)
       r0 = lerp_point(m1, p0, 0.)
       r1 = lerp_point(m1, p0, 1./3.)
       r2 = lerp_point(m1, p0, 2./3.)
       r3 = lerp_point(m1, p0, 1.)
       data.add_edge4(elem, (r0, r1, r2, r3))

    s = open("open_polygons.yaml", "wt")
    s.write("dimensions: 2\n")
    s.write("shapes:\n")

    ext = [-2, 2, -2, 2]
    n = 4
    for j in range(n):
       for i in range(n):
           idx = j * n + i
           t = float(idx) / (n * n - 1)
           dx = (ext[1] - ext[0]) / n
           dy = (ext[3] - ext[2]) / n
           margin = dx * 0.05
           e = (ext[0] + i * dx + margin,
                ext[0] + (i+1)*dx - margin,
                ext[2] + j * dy + margin,
                ext[2] + (j+1)*dy - margin)

           m = mfem()
           elem = 1 # All shapes have 1 element
           single_polygon(m, elem, e, t)
           mat = "open_polygons%02d" % idx
           filename = mat + ".mesh"
           with open(filename, "wt") as f:
               m.write(f)
               f.close()

           s.write(f" - name: {mat}\n")
           s.write(f"   material: poly\n")
           s.write("   geometry:\n")
           s.write("     format: mfem\n")
           s.write(f"     path: {filename}\n")
    s.close()

def main():
    gear()
    open_polygons()

if __name__ == "__main__":
    main()
