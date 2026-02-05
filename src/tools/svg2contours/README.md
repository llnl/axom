## Running the svg2contours Python script

The `svg2contours` script converts [SVG](https://developer.mozilla.org/en-US/docs/Web/SVG) images to [MFEM NURBS meshes](https://mfem.org/mesh-format-v1.0/#nurbs-meshes) using the [svgpathtools](https://github.com/mathandy/svgpathtools) Python library. 

The latter can be used with Axom's `quest_winding_number` example application to 
sample the winding number field over the generated curves.


### Create a virtual environment using uv

Our tool requires Python 3.9+ (see `pyproject.toml`) and [svgpathtools@1.7.2](https://github.com/mathandy/svgpathtools/releases/tag/v1.7.2) or higher.
Here are instructions for setting up an environment with uv:

```shell
# Optionally add uv if not already present; also ensure that it's on your PATH
# > python3 -m pip install uv

> cd <axom_root>/src/tools/svg2contours

# Optionally, create the initial environment -- this will generate the pyproject.toml
# > uv init
# > uv add svgpathtools

# Create venv w/ deps(using pyproject.toml)
> uv venv
> uv sync
```

### Run the script on an input SVG mesh

To convert an SVG to an mfem NURBS mesh, run the following command:
```shell
> cd <axom_root>/<build_dir>
> uv run --project ../src/tools/svg2contours ../src/tools/svg2contours/svg2contours.py \
    -i ../data/contours/svg/shapes.svg

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
