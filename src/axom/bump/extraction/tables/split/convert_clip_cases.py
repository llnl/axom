import sys

# NOTE: The clipping tables in this directory are adapted from VisIt's
#       polygon clipping cases. The tables contain a mixture of triangle
#       and quad shapes. This script reads in the clipping tables and
#       combines shapes to make up to 8-sided polygons, enabling Axom's
#       clipping algorithms to produce fewer zones and fewer artificially
#       introduced mesh lines.

This script converts the polygonal clipping tables in this directory 
def read_array_cc(filename, name, convert):
   lines = open(filename, "rt").readlines()
   reading = False
   arr = []
   for line in lines:
      if not reading:
         if line.find(name) != -1:
             reading = True
      else:
         if line.find(";") != -1:
             reading = False
             break
         elif line.find("//") == 1:
             # Skip the line.
             continue
         else:
             s = 0
             e = line.find("//")
             if e == -1:
                 e = len(line) - 1
             if convert:
                 numbers = eval("[" + line[s:e] + "]")
                 arr = arr + list(numbers)
             else:
                 tokens1 = [x.replace(' ', '') for x in line[s:e].split(",")]
                 tokens = [x for x in tokens1 if x != '']
                 arr = arr + list(tokens)
   return arr
 
def find_shared_edge(poly1, poly2):
    edges1 = set((poly1[i], poly1[(i+1)%len(poly1)]) for i in range(len(poly1)))
    edges2 = set((poly2[i], poly2[(i+1)%len(poly2)]) for i in range(len(poly2)))
    for e1 in edges1:
        if e1 in edges2 or (e1[1], e1[0]) in edges2:
            return e1
    return None

def merge_polygons(poly1, poly2, shared_edge):
    # Find indexes of shared edge in both polygons
    i1 = poly1.index(shared_edge[0])
    j1 = (i1 + 1) % len(poly1)
    i2 = poly2.index(shared_edge[0])
    j2 = (i2 + 1) % len(poly2)
    # Check direction
    if poly1[j1] == shared_edge[1]:
        poly1_order = poly1
    else:
        poly1_order = poly1[::-1]
    if poly2[j2] == shared_edge[1]:
        poly2_order = poly2
    else:
        poly2_order = poly2[::-1]
    # Remove shared edge from poly2
    idx = poly2_order.index(shared_edge[0])
    merged = poly1_order + poly2_order[idx+2:] + poly2_order[:idx]
    # Remove duplicate vertices at the join
    result = []
    for v in merged:
        if not result or v != result[-1]:
            result.append(v)
    return tuple(result)

def combine_polygons(polygons, max_edges=8):
    polygons = [tuple(p) for p in polygons]
    changed = True
    while changed:
        changed = False
        new_polygons = []
        skip = set()
        for i in range(len(polygons)):
            if i in skip:
                continue
            merged = False
            for j in range(len(polygons)):
                if i == j or j in skip:
                    continue
                shared_edge = find_shared_edge(polygons[i], polygons[j])
                if shared_edge:
                    new_poly = merge_polygons(polygons[i], polygons[j], shared_edge)
                    if len(new_poly) <= max_edges:
                        # We should sort the points so merging works better
                        new_poly = tuple(sort_points(new_poly, 8))
                        #print("Created ", new_poly, "from", polygons[i], polygons[j])
                        new_polygons.append(new_poly)
                        skip.update([i, j])
                        merged = True
                        changed = True
                        break
                    #else:
                    #    print("Could not combine polygons due to length!", polygons[i], polygons[j])
            if not merged:
                new_polygons.append(polygons[i])
        polygons = new_polygons
    return polygons

def example():
    # Example usage:
    polygons = [
    ('A', 'B', 'C', 'D'),
    ('C', 'D', 'E', 'F'),
    ('E', 'F', 'G', 'H'),
    ('X', 'Y', 'Z')
    ]

    result = combine_polygons(polygons, max_edges=8)
    print(result)

def sorted_point_order(numPoints):
    order = {}
    for i in range(numPoints):
       p = f"P{i}"
       e = f"E{chr(ord('A')+i)}"
       order[p] = 2*i
       order[e] = 2*i + 1
    return order

def sort_points(pts, maxPoints):
    order = sorted_point_order(maxPoints)
    sp = {}
    for p in pts:
       sp[order[p]] = p
    sk = sorted(sp.keys())
    out = []
    for k in sk:
       out.append(sp[k])
    return out

def advanceTable():
   return {"ST_LIN" : 2 + 2,
              "ST_TRI": 2 + 3,
              "ST_QUA": 2 + 4,
              "ST_POLY5": 2 + 5,
              "ST_POLY6": 2 + 6,
              "ST_POLY7": 2 + 7,
              "ST_POLY8": 2 + 8,
              # fictional sizes to help for line conversion
              "ST_POLY9": 2 + 9,
              "ST_POLY10": 2 + 10,
              "ST_POLY11": 2 + 11,
              "ST_POLY12": 2 + 12,
             }

def shapeSize(name):
   return advanceTable()[name]

def convert_clip_cases(filename, name, npts, max_edges = 8):
   advance = advanceTable()
   nptsToShapeName = {}
   for k in advance:
       nptsToShapeName[advance[k] - 2] = k

   sizes = read_array_cc(filename, f"numClipShapes{name}", True)
   #print(sizes)

   offsets = read_array_cc(filename, f"startClipShapes{name}", True)
   #print(offsets)

   shapes = read_array_cc(filename, f"clipShapes{name}", False)
   #print(shapes)

   numCases = len(sizes)

   outSizes = [0]*numCases
   outOffsets = [0]*numCases
   outShapes = []
   outOffset = 0

   for ci in range(numCases):
       offset = offsets[ci]
       #print(f"// Case {ci}")
       color0 = []
       color1 = []
       for si in range(sizes[ci]):
           shapeTag = shapes[offset]
           colorTag = shapes[offset + 1]

           shapeSize = advance[shapeTag]
           shapeTokens = shapes[offset:offset+shapeSize]
           if colorTag == "COLOR0":
               color0.append(shapeTokens[2:])
           else:
               color1.append(shapeTokens[2:])
           #print("  ", shapeTokens)
           offset = offset + shapeSize

       outOffsets[ci] = outOffset
       if len(color0) > 0:
          c0 = combine_polygons(color0, max_edges=max_edges)
          #print("COLOR0", c0)
          for c in c0:
             sc = sort_points(c, npts)
             outShapes.append(nptsToShapeName[len(sc)])
             outShapes.append("COLOR0")
             for p in sc:
                outShapes.append(p)
             outOffset = outOffset + 2 + len(sc)
          outSizes[ci] = outSizes[ci] + len(c0)

       if len(color1) > 0:
          c1 = combine_polygons(color1, max_edges=max_edges)
          #print("COLOR1", c1)
          for c in c1:
             sc = sort_points(c, npts)
             outShapes.append(nptsToShapeName[len(sc)])
             outShapes.append("COLOR1")
             for p in sc:
                outShapes.append(p)
             outOffset = outOffset + 2 + len(sc)
          outSizes[ci] = outSizes[ci] + len(c1)

   return (outSizes, outOffsets, outShapes)

def write_new_tables(filename, name, tableNames, sizes, offsets, shapes):

    Clip = tableNames[0]
    clip = tableNames[1]
    clipping = tableNames[2]

    f = open(filename, "wt")
    f.write("// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and\n")
    f.write("// other Axom Project Developers. See the top-level LICENSE file for details.\n")
    f.write("//\n")
    f.write("// SPDX-License-Identifier: (BSD-3-Clause)\n")

    f.write(f"#include \"{Clip}Cases.hpp\"\n\n")

    f.write("namespace axom {\n")
    f.write("namespace bump {\n")
    f.write("namespace extraction {\n")
    f.write("namespace tables {\n")
    f.write(f"namespace {clipping}")
    f.write(" {\n\n")
    f.write(f"int num{Clip}Cases{name} = {len(sizes)};")
    f.write("\n\n")   

    f.write(f"int num{Clip}Shapes{name}[] = ")
    f.write("{\n")
    f.write("  " + str(sizes)[1:-1])
    f.write("\n};\n\n")
    
    f.write(f"int start{Clip}Shapes{name}[] = ")
    f.write("{\n")
    f.write("  " + str(offsets)[1:-1])
    f.write("\n};\n\n")

    f.write("// clang-format off\n")
    f.write(f"unsigned char {clip}Shapes{name}[] = ")
    f.write("{\n")

    numCases = len(sizes)
    for ci in range(numCases):
        offset = offsets[ci]
        f.write(f"  // Case #{ci}\n")
        for si in range(sizes[ci]):
            shapeTag = shapes[offset]
            ss = shapeSize(shapeTag)
            toks = shapes[offset:offset+ss]
            #print(toks)
            toks_str = str(toks)[1:-1].replace("'", "")
            if ci < numCases - 1:
               toks_str = toks_str + ","
            else:
               if si < sizes[ci] - 1:
                  toks_str = toks_str + ","
            f.write("  " + toks_str + "\n")
            offset = offset + ss
    f.write("};\n")
    f.write("// clang-format on\n\n")

    f.write(f"const size_t {clip}Shapes{name}Size = sizeof({clip}Shapes{name}) / sizeof(unsigned char);\n\n")
    f.write("}")
    f.write(f" // namespace {clipping}\n")
    f.write("} // namespace tables\n")
    f.write("} // namespace extraction\n")
    f.write("} // namespace bump\n")
    f.write("} // namespace axom\n")
    f.close()

def make_cut_cases(sizes, offsets, shapes):
    cutSizes = []
    cutOffsets = []
    cutShapes = []

    newOffset = 0
    numCases = len(sizes)
    for ci in range(numCases):
        offset = offsets[ci]
        edges = []
        for si in range(sizes[ci]):
            shapeTag = shapes[offset]
            ss = shapeSize(shapeTag)
            toks = shapes[offset:offset+ss]

            for i in range(len(toks)):
               e0 = toks[i]
               e1 = toks[(i+1)%len(toks)]
               if e0[0] == 'E' and e1[0] == 'E':
                  name = e0 + e1
                  if name not in edges:
                    edges.append(name)

            offset = offset + ss

        for e in edges:
            e0 = e[:2]
            e1 = e[2:]
            cutShapes.append("ST_LIN")
            cutShapes.append("COLOR0")
            cutShapes.append(e0)
            cutShapes.append(e1)
        cutSizes.append(len(edges))
        cutOffsets.append(newOffset)
        newOffset = newOffset + len(edges) * 4 

    return cutSizes, cutOffsets, cutShapes

def make_polygonal_clip_tables():
   tableNames = ["Clip", "clip", "clipping"]

   outSizes, outOffsets, outShapes = convert_clip_cases("ClipCasesQua.cpp", "Qua", 4, max_points=4)
   write_new_tables("../clipping/ClipCasesQua.cpp", "Qua", tableNames, outSizes, outOffsets, outShapes)

   outSizes, outOffsets, outShapes = convert_clip_cases("ClipCasesPoly5.cpp", "Poly5", 5, max_points=5)
   write_new_tables("../clipping/ClipCasesPoly5.cpp", "Poly5", tableNames, outSizes, outOffsets, outShapes)

   outSizes, outOffsets, outShapes = convert_clip_cases("ClipCasesPoly6.cpp", "Poly6", 6, max_points=6)
   write_new_tables("../clipping/ClipCasesPoly6.cpp", "Poly6", tableNames, outSizes, outOffsets, outShapes)

   outSizes, outOffsets, outShapes = convert_clip_cases("ClipCasesPoly7.cpp", "Poly7", 7, max_points=7)
   write_new_tables("../clipping/ClipCasesPoly7.cpp", "Poly7", tableNames, outSizes, outOffsets, outShapes)

   outSizes, outOffsets, outShapes = convert_clip_cases("ClipCasesPoly8.cpp", "Poly8", 8, max_points=8)
   write_new_tables("../clipping/ClipCasesPoly8.cpp", "Poly8", tableNames, outSizes, outOffsets, outShapes)

def make_polygonal_cut_tables():
   tableNames = ["Cut", "cut", "cutting"]

   outSizes, outOffsets, outShapes = convert_clip_cases("ClipCasesPoly5.cpp", "Poly5", 5, max_edges=20)
   outSizes, outOffsets, outShapes = make_cut_cases(outSizes, outOffsets, outShapes)
   write_new_tables("../cutting/CutCasesPoly5.cpp", "Poly5", tableNames, outSizes, outOffsets, outShapes)

   outSizes, outOffsets, outShapes = convert_clip_cases("ClipCasesPoly6.cpp", "Poly6", 6, max_edges=20)
   outSizes, outOffsets, outShapes = make_cut_cases(outSizes, outOffsets, outShapes)
   write_new_tables("../cutting/CutCasesPoly6.cpp", "Poly6", tableNames, outSizes, outOffsets, outShapes)

   outSizes, outOffsets, outShapes = convert_clip_cases("ClipCasesPoly7.cpp", "Poly7", 7, max_edges=20)
   outSizes, outOffsets, outShapes = make_cut_cases(outSizes, outOffsets, outShapes)
   write_new_tables("../cutting/CutCasesPoly7.cpp", "Poly7", tableNames, outSizes, outOffsets, outShapes)

   outSizes, outOffsets, outShapes = convert_clip_cases("ClipCasesPoly8.cpp", "Poly8", 8, max_edges=20)
   outSizes, outOffsets, outShapes = make_cut_cases(outSizes, outOffsets, outShapes)
   write_new_tables("../cutting/CutCasesPoly8.cpp", "Poly8", tableNames, outSizes, outOffsets, outShapes)

def main():
   if "-clip" in sys.argv:
      make_polygonal_clip_tables()
   else:
      make_polygonal_cut_tables()

main()
 
