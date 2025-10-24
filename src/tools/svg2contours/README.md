## Running the svg2contours Python script

The `svg2contours` script converts [SVG](https://developer.mozilla.org/en-US/docs/Web/SVG) images to [MFEM NURBS meshes](https://mfem.org/mesh-format-v1.0/#nurbs-meshes) using the [svgpathtools](https://github.com/mathandy/svgpathtools) Python library. 

The latter can be used with axom's `quest_winding_number` example application to 
sample the winding number field over the generated curves.

Full SVG support requires a (slightly) patched copy of [svgpathtools@1.7.1](https://github.com/mathandy/svgpathtools/releases/tag/v1.7.1), as described in this document.

### Create a virtual environment

```shell
> python3 -m venv venv

# linux
> source venv/bin/activate
# windows bash
> source ./venv/Scripts/activate

> pip3 install -r requirements.txt
```

### Apply patch to svgpathtools to fix how angles are computed from eigenvalues

[svgpathtools@1.7.1](https://github.com/mathandy/svgpathtools/releases/tag/v1.7.1) has a bug in applying the correct rotation angle to a path. The eigenvalues can sometimes be complex numbers. 
This can be resolved by applying the following patch:
```shell
> patch  -p1 venv/lib/python3.9/site-packages/svgpathtools/path.py -i svgpathtools-1.7.1-eigenvec-fix.patch --verbose 
```

#### Developer's note:
These patches were generated from a git diff after modifying a file in svgpathtools:
```shell
 > git diff  bb54e6b15701c1839a6d26bb47893d26a7bde694 93097cc5ec1800267fc04dc6750b837701cde52c > svgpathtools-1.7.1-eigenvec-fix.patch
```

### Run the script on an input SVG mesh

To convert an SVG to an mfem NURBS mesh, run the following command:
```shell
> cd <axom_root>/<build_dir>
> ../src/tools/svg2contours/svg2contours.py -i ../data/contours/svg/shapes.svg 

SVG dimensions: width='210mm' height='297mm' viewBox='0 0 210 297'
Wrote 'drawing.mesh' with 54 vertices and NURBS 27 elements
```
> :information_source: This assumes your Axom clone has the `data` submodule located at `<axom_root>/data`

### Run the quest winding number example
Now that we have an MFEM NURBS mesh, we can run our winding number application

```shell
> cd <axom_root>/<build_dir>/
> ./examples/quest_winding_number_ex      \
    -i ./drawing.mesh  \
    query_mesh --min 0 0 --max 250 250 --res 500 500 
```
