# vector-heat-demo
C++ demo of the Vector Heat Method (Sharp, Soliman, and Crane. 2019.)

The Vector Heat Method is implemented in [geometry-central](http://geometry-central.net), this is just a simple application to demonstrate the functionality.

See also the geometry-central documentation for [the Vector Heat Method](http://geometry-central.net/surface/algorithms/vector_heat_method/) and [centers on surfaces](http://geometry-central.net/surface/algorithms/surface_centers/).

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
