```
        ,/
      ,'/
    ,' /        █████╗  ██████╗██████╗  ██████╗  (a)utomated
  ,'  /_____,  ██╔══██╗██╔════╝██╔══██╗██╔════╝  (c)ontamination
.'____    ,'   ███████║██║     ██║  ██║██║       (d)etection and
     /  ,'     ██╔══██║██║     ██║  ██║██║       (c)onfidence estimation for single-cell genome data
    / ,'       ██║  ██║╚██████╗██████╔╝╚██████╗
   /,'         ╚═╝  ╚═╝ ╚═════╝╚═════╝  ╚═════╝
  /'
```

## Synopsis

Acdc is a tool to test single-cell genome data contamination. By using sophisticated dimensionality
reduction and clustering methods, it uses tetramer profiles to differentiate between different
species in a given sample. It automatically detects the number of clusters/species and provides
confidence information.

## Requirements

* C++11 compatible compiler
* [CMake](https://cmake.org/) >= 3.0
* [Boost C++ Libraries](http://www.boost.org) >= 1.55.0 (program_options, system, filesystem)

Other needed libraries are included.

### Optional dependencies

[Kraken](https://github.com/DerrickWood/kraken) >= 0.10.5-beta can be used to improve results by using a database.
[RNAmmer](http://www.cbs.dtu.dk/services/RNAmmer/) = 1.2 is used to highlight and extract 16S regions.

## Installation

Clone the repository, run CMake and build the library

```
# git clone https://github.com/mlux86/acdc.git
# mkdir build && cd build
# cmake ..
# make -j4
```

To install acdc:

```
make install
```

The installation path (prefix path) can be changed by adding the option `-DCMAKE_INSTALL_PREFIX=<prefix-path>` to the cmake command above.

## Running

To show the help, including various parameters of acdc, run

```
# acdc -h
```

To run acdc to check a single fasta file for contamination:

```
# acdc -i <path-to-fasta>
```

To run acdc in batch mode, supply a list of fasta files:

```
# acdc -I <path-to-file-with-list-of-fastas>
```

### Including Kraken results

By default, Kraken is not included in the computation. It is necessary, to specify a database for Kraken to access, i.e.:

```
# acdc -K <path-to-kraken-db> -i <path-to-fasta>
```

Make sure that the `kraken` and `kraken-translate` executables can be found in the $PATH environment variable.

### Highlighting and extraction of 16S regions

Make sure that the `rnammer` executable can be found in the $PATH environment variable and RNAmmer is working properly.

In the visualization, 16S genes are highlighted by a large star shape. A click on it will download the corresponding 16S sequence.

## Viewing results

While acdc is running, it generates results in an output directory defaulting to `./results/` (can be overridden using the `-o` parameter).
It contains the file `index.html` that can be viewed in any modern browser, supporting HTML5 and CSS3.

## Docker Container

In the `docker` folder, acdc provides a Dockerfile script to build a docker container. Building depends on the rnammer sources which can be obtained via the [RNAmmer homepage](http://www.cbs.dtu.dk/services/RNAmmer/). Just put the obtained `rnammer-1.2.src.tar.Z` file into the `docker` folder and from inside the directory, run 

```
# docker build -t acdc .
```

The built acdc docker image shares result files using the `/acdc` and `/usr/share/nginx/html` volumes. The latter can be used by the official [nginx](https://nginx.org/) web server [docker container](https://hub.docker.com/_/nginx/). For example as follows:

```
# docker run --name acdc \
             -v /path/to/assemblies:/assemblies \
             -v /path/to/kraken_db:/krakendb \
             acdc -i test.fasta

# docker run --volumes-from acdc -p 80:80 --rm nginx
```

where `/path/to/assemblies` contains the file `test.fasta` and `/path/to/kraken_db` contains the Kraken database. The second command runs an nginx container which serves the newly created acdc result files on port 80.

## Used Libraries

* [Eigen](http://eigen.tuxfamily.org/), a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
* [SeqAn](http://www.seqan.de/), an open source C++ library of efficient algorithms and data structures for the analysis of sequences with the focus on biological data.
* [t-SNE](https://lvdmaaten.github.io/tsne/), t-Distributed Stochastic Neighbor Embedding (t-SNE) is a technique for dimensionality reduction that is particularly well suited for the visualization of high-dimensional datasets.
* [diptest: Hartigan's Dip Test Statistic for Unimodality](https://cran.r-project.org/web/packages/diptest/), by Martin Maechler (originally from Fortran and S-plus by Dario Ringach, NYU.edu), modified for use in a C++ application.
* [fastcluster](http://danifold.net/fastcluster.html) fastcluster: Fast hierarchical clustering routines for R and Python, extracted relevant code and added pure C++ interface
* [JsonCpp](https://github.com/open-source-parsers/jsoncpp), a C++ library for interacting with JSON.
* [nanoflann](https://github.com/jlblancoc/nanoflann), a C++ header-only library for Nearest Neighbor (NN) search wih KD-trees.
* [Catch](https://github.com/philsquared/Catch), a modern, C++-native, header-only, framework for unit-tests, TDD and BDD.
* [jQuery](https://jquery.com/), a fast, small, and feature-rich JavaScript library.
* [D3.js - Data-Driven Documents](https://d3js.org/), a JavaScript library for visualizing data with HTML, SVG, and CSS.

## License

The MIT License (MIT)

Copyright (c) [2016] [Markus Lux]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
