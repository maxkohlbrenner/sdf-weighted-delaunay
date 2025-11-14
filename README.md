# Code: Isosurface Extraction for Signed Distance Functions using Power Diagrams
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

The cpp file is used as follows:
```
./build/sdf-weighted-delaunay path_to_obj_file (N) (max_refinement) (delaunay) (outpath)
```

The arguments are:
| Argument | Meaning | Default | 
| -------- | ------- | -------- | 
| path_to_obj_file | path to surface mesh in .obj file format                    | required | 
| N | number of regular samples per dimension                                    | optional, defaults to 10 (10x10x10 regular grid) |
|max_refinement| number of refinement steps                                      | optional, defaults to 1000 | 
|delaunay| whether to use delaunay instead of regular triangulation / refinement | optional, defaults to 0 (false) | 
|outpath| folder in which an output contour mesh is stored                       | optional, a viewer is shown instead if not provided and not output is written | 

## Examples:
We show several examples in `Example.ipynb`

#### Short form: 

Assume a surface mesh lies in `data/Armadillo.obj`. 

For marching tets based on the regular triangulation (contouring only, not incremental) on a regular grid of `N=30` samples per dimension call:
```
./build/sdf-weighted-delaunay ./data/Armadillo.obj 30 0 0 ./
```

If you additionally want to perform 1000 incremental steps, call:

```
./build/sdf-weighted-delaunay ./data/Armadillo.obj 30 1000 0 ./
```

In case you are interested in Delaunay refinement (which we used as a baseline in the paper), call

```
./build/sdf-weighted-delaunay ./data/Armadillo.obj 30 1000 1 ./
```
