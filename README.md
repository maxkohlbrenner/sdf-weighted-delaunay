# sdf-weighted-delaunay
Code for:

Kohlbrenner, M. and Alexa, M. (2025), Isosurface Extraction for Signed Distance Functions using Power Diagrams. Computer Graphics Forum e70037. https://doi.org/10.1111/cgf.70037

# Building

```
git clone --recursive git@github.com:maxkohlbrenner/sdf-weighted-delaunay.git
cd sdf-weighted-delaunay/
mkdir build
cmake -B build -DCGAL_DIR=./cgal 
cmake --build build/ --parallel
```

# Usage

```
./build/sdf-weighted-delaunay path_to_obj_file N max_refinement delaunay outpath
```

The arguments are:

path_to_obj_file : Surface mesh in .obj file format <br>
N : number of regular samples per dimension (optional, defaults to 10, a 10x10x10 regular grid) <br>
max_refinement: number of refinement steps  (optional, defaults to 1000) <br>
delaunay: delaunay flag, whether to use delaunay instead of regular triangulation / refinement (optional, defaults to 0, false) <br>
outpath: folder in which an output contour mesh is stored (optional, a viewer is shown instead if not provided and not output is written) <br>
