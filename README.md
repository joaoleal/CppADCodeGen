# CppADCodeGen [![Build Status](https://travis-ci.org/joaoleal/CppADCodeGen.svg?branch=master)](https://travis-ci.org/joaoleal/CppADCodeGen) [![DOI](https://zenodo.org/badge/20828/joaoleal/CppADCodeGen.svg)](https://zenodo.org/badge/latestdoi/20828/joaoleal/CppADCodeGen)

CppADCodeGen performs **hybrid Automatic Differentiation** (AD), that is, uses 
operator-overloading and produces source-code. Such source-code can be 
statically compiled at runtime using an existing compiler and linked dynamically 
or, alternatively, go through a JIT compilation using Clang/LLVM.

In addition to C source generation, CppADCodeGen can also produce 
 [Latex](http://www.latex-project.org/),
 html+[MathML](http://www.w3.org/Math/), and 
 [dot](https://en.wikipedia.org/wiki/DOT_%28graph_description_language%29)
 source-code files for your algorithm.
Latex sources can be used to create PDF files for documentation purposes, 
html+MathML can be used to display your algorithm in a web browser, and
dot files can be used to create images with a graph of your model 
(see [graphviz](http://graphviz.org/)).

CppADCodeGen can also be used to perform differentiation index reduction of 
Differential Algebraic Equations (DAE) through Pantelides, Soares-Secchi, and Dummy 
Derivatives methods.

CppADCodeGen is built on top of the [CppAD](http://www.coin-or.org/CppAD) 
library, which is a header only C++ AD library using operator overloading.

## License ##

CppADCodeGen is available with both the **EPL** and **GPL** licenses 
(suitable for both open-source and closed-source commercial projects).
See epl-v10.txt and gpl3.txt for a copy of the licenses.

## Requirements ##

CppADCodeGen is a C++11 header only library, therefore there aren't many dependencies:

 - **CppAD** (2017),
 - A **C++11** compiler (such as GCC and Clang),
 - Clang/LLVM (only for JIT compilation), and
 - Eigen 3 (only for DAE differentiation index reduction).

Runtime compilation and dynamic linking:
 - Linux (it might be very easy to support other OSes but it is not implemented yet)

## Installing ##

### General installation ###

Just copy the contents of the folder include to anywhere you would like to 
include from.

### Debian/Ubuntu ###

A debian installation package can be created at the root of the project.
Typically you can create the installer by just typing:

    dpkg-buildpackage

It will create a debian package outside the project's folder.

## Using CppADCodeGen ##

See the [wiki](https://github.com/joaoleal/CppADCodeGen/wiki).

The folder example includes some simple use cases.

---

## Repository Content

|Directories |  Description                                                    |
|------------|-----------------------------------------------------------------|
|bin         | Helper shell and sed scripts used for CppAD development.        |
|bug         | Directory containing demonstration of known bugs (may be empty).|
|doc         | Holds files for generation of developer documentation.          |
|debian      | Debian package creation files (Linux).                          |
|include     | The CppADCodeGen header files.                                  |
|example     | CppADCodegen example files are here.                            |
|pkgconfig   | pkg-config support files.                                       |
|test        | Contains tests for CppADCodeGen.                                |
|speed       | Contains some benchmarks for CppADCodeGen.                      |


| Files         |  Description                                                 |
|---------------|--------------------------------------------------------------|
|AUTHORS        | Statement of authorship and copyright.                       |
|CMakeLists.txt | CppADCodeGen CMake input file.                               |
|COPYING        | Statement of user license to use software.                   |
|epl-v10.txt    | A copy of the Eclipse Public License version 1.              |
|gpl3.txt       | A copy of the GNU General Public License version 3.          |
|INSTALL        | Points to this file.                                         |
