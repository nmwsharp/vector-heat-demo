# vector-heat-demo
C++ demo of the [Vector Heat Method (Sharp, Soliman, and Crane 2019)](https://www.cs.cmu.edu/~kmcrane/Projects/VectorHeatMethod/index.html) & [Affine Heat Method (Soliman & Sharp 2025)](https://www.yousufsoliman.com/projects/the-affine-heat-method.html).

The main implementation of these algorithms is in [geometry-central](http://geometry-central.net)
- [C++ library implementation](https://geometry-central.net/surface/algorithms/vector_heat_method/)
- [Python bindings](https://github.com/nmwsharp/potpourri3d)

this repository just a simple application, wrapping the C++ library in a UI to demonstrate the functionality.

See also the geometry-central documentation for [the Vector Heat Method](http://geometry-central.net/surface/algorithms/vector_heat_method/) and [centers on surfaces](http://geometry-central.net/surface/algorithms/surface_centers/).

![demo git](https://github.com/nmwsharp/vector-heat-demo/blob/master/vector_heat_demo.gif)

### Building and running

```
git clone --recurse-submodules https://github.com/nmwsharp/vector-heat-demo.git
cd vector-heat-demo
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
./bin/vector_heat /path/to/your/mesh.obj
```

This will open a UI window showing your mesh. The command window in the upper right can be used to run the algorithms, while the window on the left adjusts visualization setting. Note that both transport sources and averaging sites can be selected in these windows; some arbtrirary vertices are selected initially.

## Citation

If these algorithms contribute to academic work, please cite the following paper(s):

```bib
@article{sharp2019vector,
  title={The Vector Heat Method},
  author={Sharp, Nicholas and Soliman, Yousuf and Crane, Keenan},
  journal={ACM Transactions on Graphics (TOG)},
  volume={38},
  number={3},
  pages={24},
  year={2019},
  publisher={ACM}
}

@article{soliman2025affine,
  title={The Affine Heat Method},
  author={Yousuf Soliman, Nicholas Sharp,
  booktitle={Computer Graphics Forum},
  volume={44},
  number={5},
  year={2025}
}
```