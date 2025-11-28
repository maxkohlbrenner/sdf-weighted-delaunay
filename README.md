# Code: Isosurface Extraction for Signed Distance Functions using Power Diagrams
Kohlbrenner, M. and Alexa, M. (2025), Isosurface Extraction for Signed Distance Functions using Power Diagrams. Computer Graphics Forum e70037. https://doi.org/10.1111/cgf.70037

We would like to thank Pierre Alliez for letting us use his implementation of Delauany refinement.

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

# Citation
If you use this code in your research, pleace cite
```
@article{https://doi.org/10.1111/cgf.70037,
author = {Kohlbrenner, M. and Alexa, M.},
title = {Isosurface Extraction for Signed Distance Functions using Power Diagrams},
journal = {Computer Graphics Forum},
volume = {44},
number = {2},
pages = {e70037},
keywords = {Signed Distance Function, Contouring, Power Diagram, Weighted Delaunay Triangulation, CCS Concepts, • Computing methodologies → Shape analysis, • Mathematics of computing → Mesh generation},
doi = {https://doi.org/10.1111/cgf.70037},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.70037},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1111/cgf.70037},
abstract = {Abstract Contouring an implicit function typically considers function values in the vicinity of the desired level set, only. In a recent string of works, Sellán at al. have demonstrated that signed distance values contain useful information also if they are further away from the surface. This can be exploited to increase the resolution and amount of detail in surface reconstruction from signed distance values. We argue that the right tool for this analysis is a regular triangulation of the distance samples, with the weights chosen based on the distance values. The resulting triangulation is better suited for reconstructing the surface than a standard Delaunay triangulation of the samples. Moreover, the dual power diagram encodes the envelope enclosing the surface, consisting of spherical caps. We discuss how this information can be exploited for reconstructing the surface. In particular, the approach based on regular triangulations lends itself well to refining the sample set. Refining the sample set based on the power diagram outperforms other reconstruction methods relative to the sample count.},
year = {2025}
}
```
