from collections import defaultdict

def extract_triangles(arr):
    """Extract valid triangles from the input array."""
    triangles = []
    for i in range(0, len(arr), 3):
        tri = arr[i:i+3]
        if -1 not in tri and len(tri) == 3:
            triangles.append(tuple(tri))
    return triangles

#def edges(triangle):
#    """Return the set of edges (as sorted tuples) from a triangle."""
#    a, b, c = triangle
#    return {(min(a, b), max(a, b)), (min(b, c), max(b, c)), (min(c, a), max(c, a))}
def edges(polygon):
    return sorted(polygon)


def merge_polygons(poly1, poly2):
    """Merge two polygons that share an edge into a larger polygon."""
    n1, n2 = len(poly1), len(poly2)
    # Find the shared edge
    for i in range(n1):
        a1, a2 = poly1[i], poly1[(i+1)%n1]
        for j in range(n2):
            b1, b2 = poly2[j], poly2[(j+1)%n2]
            if {a1, a2} == {b1, b2}:
                # Merge poly2 into poly1 at the shared edge
                # Remove the shared edge from both
                # poly1: ... a1 a2 ...
                # poly2: ... b1 b2 ...
                # We want: ... a2 ... a1 ... (excluding the shared edge)
                # Find the order to append poly2
                # poly1: a1 a2 ... (rest)
                # poly2: b2 ... (rest) b1
                # Remove the shared edge from both
                idx1 = i
                idx2 = j
                # Build new polygon
                new_poly = []
                # Add poly1 from a2 (next after shared edge) to a1 (before shared edge)
                k = (idx1 + 1) % n1
                while k != idx1:
                    new_poly.append(poly1[k])
                    k = (k + 1) % n1
                # Add poly2 from b2 (next after shared edge) to b1 (before shared edge), in reverse
                k = (idx2 + 1) % n2
                temp = []
                while k != idx2:
                    temp.append(poly2[k])
                    k = (k + 1) % n2
                new_poly += temp
                return tuple(new_poly)
    return None

def combine_triangles(arr):
    triangles = extract_triangles(arr)
    polygons = [tri for tri in triangles]
    merged = True
    while merged:
        merged = False
        n = len(polygons)
        for i in range(n):
            for j in range(i+1, n):
                poly1, poly2 = polygons[i], polygons[j]
                # If they share an edge
                if set(edges(poly1)) & set(edges(poly2)):
                    new_poly = merge_polygons(poly1, poly2)
                    if new_poly:
                        # Replace poly1 and poly2 with new_poly
                        polygons = [p for k, p in enumerate(polygons) if k not in (i, j)] + [new_poly]
                        merged = True
                        break
            if merged:
                break
    return tuple(polygons)

def read_array(filename, name):
   lines = open(filename, "rt").readlines()
   reading = False
   polys = []
   for line in lines:
      if not reading:
         if line.find(name) != -1:
             reading = True
      else:
         if line.find(";") != -1:
             reading = False
             break
         else:
             s = line.find("{") + 1
             e = line.find("}")
             numbers = eval("[" + line[s:e] + "]")
             polys.append(numbers)
   return polys
    
def write_new_table(filename, name, polys):
    edgeNames = ("EA", "EB", "EC", "ED",
                 "EE", "EF", "EG", "EH",
                 "EI", "EJ", "EK", "EL")
    polyNames = ("", "", "",
                 "ST_TRI,  ",
                 "ST_QUAD, ",
                 "ST_POLY5,",
                 "ST_POLY6,",
                 "ST_POLY7,",
                 "ST_POLY8,")

    f = open(filename, "wt")
    f.write("#include \"CutCases.h\"\n\n")

    f.write("namespace axom {\n")
    f.write("namespace bump {\n")
    f.write("namespace cutting {\n")
    f.write("namespace tables {\n\n")
    f.write(f"int numCutCases{name} = {len(polys)};\n\n")   
    f.write("// clang-format off\n")
    f.write(f"unsigned char cutShapes{name}[] = ")
    f.write("{\n")
    offset = 0
    offsets = []
    sizes = []
    i = 0
    for p in polys:
        result = combine_triangles(p)
        f.write(f"  // Case {i}\n")
        offsets.append(offset)
        sizes.append(len(result))
        for shape in result:
           shapeName = polyNames[len(shape)]
           offset = offset + 2 + len(shape)
           edges = [edgeNames[x] for x in shape]
           shape_str = str(edges)[1:-1]
           shape_str = shape_str.replace("'", "")
           s = f"  {shapeName} COLOR0, {shape_str},\n"
           f.write(s)
        i = i + 1
    f.write("};\n")
    f.write("// clang-format on\n\n")

    f.write(f"unsigned char numCutShapes{name}[] = ")
    f.write("{\n")
    f.write(str(sizes)[1:-1])
    f.write("\n};\n\n")
    
    f.write(f"unsigned char startCutShapes{name}[] = ")
    f.write("{\n")
    f.write(str(offsets)[1:-1])
    f.write("\n};\n\n")

    f.write(f"const size_t cutShapes{name}Size = sizeof(cutShapes{name}) / sizeof(unsigned char);\n\n")
    f.write("} // namespace tables\n")
    f.write("} // namespace cutting\n")
    f.write("} // namespace bump\n")
    f.write("} // namespace axom\n")
    f.close()

# Example usage:
arr = [1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
result = combine_triangles(arr)
print(result)

polys = read_array("vtkTriangulationTable.C", "tetTriangulationTable")
write_new_table("CutCasesTet.cpp", "Tet", polys)

polys = read_array("vtkTriangulationTable.C", "pyramidTriangulationTable")
write_new_table("CutCasesPyramid.cpp", "Pyramid", polys)

polys = read_array("vtkTriangulationTable.C", "wedgeTriangulationTable")
write_new_table("CutCasesWedge.cpp", "Wedge", polys)

polys = read_array("vtkTriangulationTable.C", "hexTriangulationTable")
write_new_table("CutCasesHex.cpp", "Hex", polys)
