GNU Octave -- a high-level language for numerical computations
==============================================================

Copyright (C) 1996-2022 The Octave Project Developers

See the file COPYRIGHT.md in the top-level directory of this
distribution or <https://octave.org/copyright/>.

Code in The No-KLU Directory
----------------------------
### Contains : 
- **Custom Implementations of NVector** 
  - `nvector_octave.h`
  - `nvector_octave.cpp`
- **Custom Implementations of SUNMatrix** 
  - `sunmatrix_oct.h`
  - `sunmatrix_oct.cpp`
- **Examples to test custom implementations on**
  - For NVector :
    - `kry_example.c`
- **Tests to test individual functions of custom implementations**
  - For NVector
    Given in SUNDIALS to test each function's expected output
    - `test_nvector.h`
    - `test_nvector.cpp`
    Main Code to call tests for functions implemented by us
    - `test_nvector_octave`
- **Building**
  - `Makefile`

### Status :
While the NVector implementation can be used in the example shown, it is currently being worked on function-by-function to make them wrappers to ColumnVector methods rather than elementwise operations on ColumnVector entries, to take advantage of vectorization in octave.

### How to build :
Since the makefile is made without using Autotools, 
the paths for required libraries are stored in
personalizable variables that must be changed for each build system.
The variables to look out for are :
```
# include+lib paths to be defined by the user temporarily
KLU_INC_PATH = < path to include files for Suitesparse/KLU >
SUN_INC_PATH = < path to include files for Sundials >
OCT_INC_PATH = < path to include files for Octave >

KLU_LIB_PATH = < path to library files for Suitesparse/KLU >
SUN_LIB_PATH = < path to library files for Sundials >
OCT_LIB_PATH = < path to library files for Octave >
```

The Makefile contains the following rules :
  - `app`
  Is currently linked to build the `kry_example.c` 
  with all the reuired dependencies, creates the executable `app`
  - `test`
  Builds the `test_nvector_app.cpp`, that calls tests for
  implemented functions. **Note :** The commented calls are
  to functions not implemented, either because they are not needed
  or are still a work in progress, creates the executable `test`
  - `clean`
  Removes all object and executable files created in the no-klu/src 
  directory for a fresh build.

After navigating to the correct directory i.e
`/octave/no-klu/src`, the expected output for each make rule should be :
```
/octave/no-klu/src$ make app
Makefile:31: warning: overriding recipe for target 'app'
Makefile:25: warning: ignoring old recipe for target 'app'
Makefile:34: warning: overriding recipe for target 'test'
Makefile:28: warning: ignoring old recipe for target 'test'
g++ -g  -I/usr/include/suitesparse -I/usr/local/include -I/usr/include/octave-8.0.0/octave -L/usr/lib/x86_64-linux-gnu -L/usr/local/lib -L/usr/lib/octave/8.0.0 -lsundials_generic -lsundials_nvecserial -lsundials_ida -lsundials_sunlinsolklu -loctave -c nvector_octave.cpp 
make: Circular app <- app dependency dropped.
g++ -g  kry_example.c nvector_octave.o -I/usr/include/suitesparse -I/usr/local/include -I/usr/include/octave-8.0.0/octave -L/usr/lib/x86_64-linux-gnu -L/usr/local/lib -L/usr/lib/octave/8.0.0 -lsundials_generic -lsundials_nvecserial -lsundials_ida -lsundials_sunlinsolklu -loctave -o app
```
```
octave/no-klu/src$ make test
Makefile:31: warning: overriding recipe for target 'app'
Makefile:25: warning: ignoring old recipe for target 'app'
Makefile:34: warning: overriding recipe for target 'test'
Makefile:28: warning: ignoring old recipe for target 'test'
g++ -g  -I/usr/include/suitesparse -I/usr/local/include -I/usr/include/octave-8.0.0/octave -L/usr/lib/x86_64-linux-gnu -L/usr/local/lib -L/usr/lib/octave/8.0.0 -lsundials_generic -lsundials_nvecserial -lsundials_ida -lsundials_sunlinsolklu -loctave -c test_nvector_octave.cpp 
g++ -g  -I/usr/include/suitesparse -I/usr/local/include -I/usr/include/octave-8.0.0/octave -L/usr/lib/x86_64-linux-gnu -L/usr/local/lib -L/usr/lib/octave/8.0.0 -lsundials_generic -lsundials_nvecserial -lsundials_ida -lsundials_sunlinsolklu -loctave -c test_nvector.cpp 
make: Circular test <- test dependency dropped.
g++ -g  test_nvector_octave.o test_nvector.o nvector_octave.o -I/usr/include/suitesparse -I/usr/local/include -I/usr/include/octave-8.0.0/octave -L/usr/lib/x86_64-linux-gnu -L/usr/local/lib -L/usr/lib/octave/8.0.0 -lsundials_generic -lsundials_nvecserial -lsundials_ida -lsundials_sunlinsolklu -loctave -o test
```
```
octave/no-klu/src$ make clean
Makefile:31: warning: overriding recipe for target 'app'
Makefile:25: warning: ignoring old recipe for target 'app'
Makefile:34: warning: overriding recipe for target 'test'
Makefile:28: warning: ignoring old recipe for target 'test'
rm -f *.o *~ app  test 
```
### Reproducing the code :
Once succesfully built, the executables app and test can be called on the command line as follows:

`./app` : Runs the Kry example for the 2D Heat Equation and shows the solver's output
`./test <vector length> <print timing>` : Creates a new vector of the specified length, and prints the status of each test run for the vector implementation, along with the timing if `<print timing>` is inputted as **1**, or without if it is inputted as **0**

Overview
--------

GNU Octave is a high-level interpreted language, primarily intended
for numerical computations.  It provides capabilities for the
numerical solution of linear and nonlinear problems, and for
performing other numerical experiments.  It also provides extensive
graphics capabilities for data visualization and manipulation.  GNU
Octave is normally used through its interactive interface (CLI and
GUI), but it can also be used to write non-interactive programs.
The GNU Octave language is quite similar to Matlab so that most
programs are easily portable.

GNU Octave is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GNU Octave is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Octave; see the file COPYING.  If not, see
<https://www.gnu.org/licenses/>.

Availability
------------

The latest released version of Octave is always available from
<https://ftp.gnu.org/gnu/octave/> and many mirror sites around the
world.  You may also find links to binary distributions at
<https://www.octave.org/download.html>.  The current development
sources may be found under the Source Code tab on
[Savannah](https://savannah.gnu.org/projects/octave/).

Installation
------------

Octave requires approximately 475 MB of disk storage to unpack and
compile from source (significantly more, 3.8 GB, if you compile with
debugging symbols).  Once installed, Octave requires approximately
75 MB of disk space (again, considerably more, 415 MB, if you don't
build shared libraries or the binaries and libraries include
debugging symbols).

To compile Octave, you will need a recent version of:

- [GNU Make](https://www.gnu.org/software/make/)
- [GNU G++](https://gcc.gnu.org/) or another C++11 compiler
- [GNU Fortran](https://gcc.gnu.org/fortran/), another Fortran 77
  compiler, or [f2c](http://www.netlib.org/f2c/)

Octave's Makefiles use features of GNU Make that are not present in
other versions of make.  If you use `f2c`, you will need a script
like `fort77` that works like a normal Fortran compiler by combining
`f2c` with your C compiler in a single script.

See the notes in the file `INSTALL.OCTAVE` of the Octave source
distribution for more detailed installation instructions.

Bugs and Patches
----------------

The file `BUGS` (or `doc/interpreter/bugs.txi`) explains the recommended
procedure for reporting bugs on the [bug tracker](https://bugs.octave.org)
or contributing patches; online documentation is also available
[here](https://www.gnu.org/software/octave/bugs.html).

Documentation
-------------

* [Octave's manual](https://www.octave.org/doc/interpreter/) is a
  comprehensive user guide covering introductive and more advanced
  topics.
* [Octave's wiki](https://wiki.octave.org) is a user community page,
  covering various topics and answering
  [FAQ](https://wiki.octave.org/FAQ).
* [Octave's Doxygen](https://www.octave.org/doxygen/) documentation
  explains the C++ class libraries.

Partially, the up-to-dateness of the documentation is lagging a bit
behind the development of the software.  If you notice omissions or
inconsistencies, please report them at our bug tracker.  Specific
suggestions for ways to improve Octave and its documentation are
always welcome.  Reports with patches are even more welcome.

Additional Information
----------------------

Up to date information about Octave is available on the WWW at
<https://www.octave.org>, or ask for help via email
<help@octave.org>.
