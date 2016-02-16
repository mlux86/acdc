        ,/
      ,'/
    ,' /        █████╗  ██████╗██████╗  ██████╗  (a)utomated
  ,'  /_____,  ██╔══██╗██╔════╝██╔══██╗██╔════╝  (c)ontamination
.'____    ,'   ███████║██║     ██║  ██║██║       (d)etection and
     /  ,'     ██╔══██║██║     ██║  ██║██║       (c)onfidence estimation for NGS data
    / ,'       ██║  ██║╚██████╗██████╔╝╚██████╗
   /,'         ╚═╝  ╚═╝ ╚═════╝╚═════╝  ╚═════╝
  /'

## Synopsis

Acdc is a tool to test next-generation-sequencing (NGS) data from single-cell sequencing for contamination. By using sophisticated dimensionality reduction and clustering methods, it uses tetramer profiles to differentiate between different species in a given sample. It automatically detects the number of clusters/species and provides confidence information. 

## Requirements

* C++11 compatible compiler
* [CMake](https://cmake.org/) >= 3.0
* [Boost C++ Libraries](http://www.boost.org) >= 1.58.0 (program_options, system, filesystem)
* [Kraken](https://github.com/DerrickWood/kraken) >= 0.10.5-beta

Other needed libraries are included.

## Installation

Clone the repository, run Cmake and build the library

```
# git clone https://github.com/mlux86/acdc.git
# mkdir build && cd build
# cmake ..
# make -j4
```

Installation of the binaries is not supported, yet (but will be in the future).

## Running

To show the help, including various parameters of acdc, run

```
# bin/acdc -h
```

To run acdc to check a single fasta file for contamination:

```
# bin/acdc -i <path-to-fasta>
```

To run acdc in batch mode, supply a list of fasta files:

```
# bin/acdc -I <path-to-file-with-list-of-fastas>
```

## Viewing results

While acdc is running, it generates results in an output directory defaulting to `./results/` (can be overridden using the `-o` parameter).
It contains the file `index.html` that can be viewed in any modern browser, supporting HTML5 and CSS3.

## Used Libraries

* [SeqAn](http://www.seqan.de/), an open source C++ library of efficient algorithms and data structures for the analysis of sequences with the focus on biological data.
* [t-SNE](https://lvdmaaten.github.io/tsne/), t-Distributed Stochastic Neighbor Embedding (t-SNE) is a technique for dimensionality reduction that is particularly well suited for the visualization of high-dimensional datasets.
* [diptest: Hartigan's Dip Test Statistic for Unimodality](https://cran.r-project.org/web/packages/diptest/), by Martin Maechler (originally from Fortran and S-plus by Dario Ringach, NYU.edu), modified for use in a C++ application.
* [JsonCpp](https://github.com/open-source-parsers/jsoncpp), A C++ library for interacting with JSON.
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