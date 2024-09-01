# CppADCodeGen

CppADCodeGen is a C++ library that extends [CppAD](http://www.coin-or.org/CppAD) 
to allow the generation of C/C++ source code for computing the derivatives of
mathematical models using Algorithmic Differentiation (AD).
Since CppAD uses operator-overloading and CppADCodeGen produces source-code,
the result is **hybrid Automatic Differentiation**.

The evaluation of differential information can be orders of magnitude faster
to compute using a compiled model than using a regular operator-overloading
strategy.

## Key Features ##

 - Code Generation:
   - C/C++: It can be used to compute derivatives of functions/models (see [wiki](https://github.com/joaoleal/CppADCodeGen/wiki/DirectSourceGeneration#c-language)) and generate libraries (see [wiki](https://github.com/joaoleal/CppADCodeGen/wiki/LibGeneration)).
   - [Latex](http://www.latex-project.org/): Latex sources can be used to create PDF files for documentation purposes (see [wiki](https://github.com/joaoleal/CppADCodeGen/wiki/DirectSourceGeneration#latex)).
   - html+[MathML](http://www.w3.org/Math/): tml+MathML can be used to display your algorithm in a web browser (see [wiki](https://github.com/joaoleal/CppADCodeGen/wiki/DirectSourceGeneration#htmlmathmljavascript)).
   - [dot](https://en.wikipedia.org/wiki/DOT_%28graph_description_language%29) source-code: dot files can be used to create images with a graph of your model
     (see [graphviz](http://graphviz.org/)).
 - Support for Multiple Types: Allows the creation of CppAD models using double or float floating-point types.
 - Efficient Sparsity Handling: Supports generation of code optimized for sparse Jacobians and Hessians.
 - Differentiation index reduction of Differential Algebraic Equations (DAE) through:
   - Pantelides,
   - Soares-Secchi, and
   - Dummy Derivatives methods.
 - **Statically compile** generated source-code at runtime can be using an existing compiler and linked
   dynamically or, alternatively, 
 - Generated source-code can be **JIT** compilation using Clang/LLVM.

## License ##

CppADCodeGen is available with both the **EPL** and **GPL** licenses
(suitable for both open-source and closed-source commercial projects).
See epl-v10.txt and gpl3.txt for a copy of the licenses.

## Requirements ##

CppADCodeGen is a C++14 header only library, therefore there aren't many dependencies:

 - Required:
   - [**CppAD**](https://github.com/coin-or/CppAD) (2024),
   - A **C++14** compiler (such as GCC and Clang),
 - Optional:
   - Clang/LLVM (only required for JIT compilation; supported versions <= v10.0), and
   - [Eigen 3](https://gitlab.com/libeigen/eigen) (required when DAE differentiation index reduction is used).

Runtime compilation and dynamic linking:
 - Linux 

It might be very easy to support other OSes, but it is not implemented yet.

## Installing ##

### General installation ###

Get the sources from GitHub:
```sh
git clone https://github.com/joaoleal/CppADCodeGen.git CppADCodeGen
```
Create a new folder to build the project:
```sh
mkdir cppadcg-build
```
Build the project (no compilation of C/C++ occurs, just generation of header files):
```sh
cd cppadcg-build
cmake ../CppADCodeGen
```
Either install the project in your system:
```sh
make install
```
or to some other folder:
```sh
make DESTDIR=someotherfolder install
```

### Debian/Ubuntu ###

A debian installation package can be created at the root of the project.
Typically, you can create the installer by just typing:
```sh
dpkg-buildpackage
```
It will create a debian package outside the project's folder.

## Using CppADCodeGen ##

See the [wiki](https://github.com/joaoleal/CppADCodeGen/wiki).

The folder example includes some simple use cases.

## Testing ##

Get the sources from GitHub:
```sh
git clone https://github.com/joaoleal/CppADCodeGen.git CppADCodeGen
```
Create a new folder for the tests:
```sh
cd make-build-debug
cmake ../CppADCodeGen
```
Testing requires [google-test](https://github.com/google/googletest) (version 1.14.0) which will be downloaded from GitHub.

Then compile the tests:
```sh
make build_tests
```

Run the complete set of tests:
```sh
make test
```
If [valgrind](https://valgrind.org/) is available in your system, CppADCodeGen will also perform memory checks which can
lead to a very lengthy test execution.
It is possible to disable memory validations by turning off the CMake option `USE_VALGRIND`.
For instance, by calling the following command before running the tests:
```sh
cmake -DUSE_VALGRIND=OFF ../CppADCodeGen 
```
---

## Repository Content

|Directories |  Description                                                    |
|------------|-----------------------------------------------------------------|
|bin         | Helper shell and sed scripts used for CppAD development.        |
|bug         | Directory containing demonstration of known bugs (may be empty).|
|debian      | Debian package creation files (Linux).                          |
|doc         | Holds files for generation of developer documentation.          |
|example     | CppADCodegen example files are here.                            |
|include     | The CppADCodeGen header files.                                  |
|pkgconfig   | pkg-config support files.                                       |
|python      | Pretty printers for GDB (debugging).                            |
|speed       | Contains some benchmarks for CppADCodeGen.                      |
|test        | Contains tests for CppADCodeGen.                                |


| Files         |  Description                                                 |
|---------------|--------------------------------------------------------------|
|AUTHORS        | Statement of authorship and copyright.                       |
|CMakeLists.txt | CppADCodeGen CMake input file.                               |
|COPYING        | Statement of user license to use software.                   |
|epl-v10.txt    | A copy of the Eclipse Public License version 1.              |
|gpl3.txt       | A copy of the GNU General Public License version 3.          |
|INSTALL        | Points to this file.                                         |
