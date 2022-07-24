octave/no-klu/src$ make app
Makefile:31: warning: overriding recipe for target 'app'
Makefile:25: warning: ignoring old recipe for target 'app'
Makefile:34: warning: overriding recipe for target 'test'
Makefile:28: warning: ignoring old recipe for target 'test'
g++ -g  -I/usr/include/suitesparse -I/usr/local/include -I/usr/include/octave-8.0.0/octave -L/usr/lib/x86_64-linux-gnu -L/usr/local/lib -L/usr/lib/octave/8.0.0 -lsundials_generic -lsundials_nvecserial -lsundials_ida -lsundials_sunlinsolklu -loctave -c nvector_octave.cpp 
nvector_octave.cpp:694:19: error: conflicting declaration of C function ‘octave_value_list N_VDotProd_Octave(N_Vector, N_Vector)’
  694 | octave_value_list N_VDotProd_Octave(N_Vector x, N_Vector y)
      |                   ^~~~~~~~~~~~~~~~~
In file included from nvector_octave.cpp:6:
nvector_octave.h:75:24: note: previous declaration ‘double N_VDotProd_Octave(N_Vector, N_Vector)’
   75 | SUNDIALS_EXPORT double N_VDotProd_Octave(N_Vector x, N_Vector y);
      |                        ^~~~~~~~~~~~~~~~~
nvector_octave.cpp: In function ‘octave_value_list N_VDotProd_Octave(N_Vector, N_Vector)’:
nvector_octave.cpp:714:44: error: no matching function for call to ‘octave_value_list::octave_value_list(<brace-enclosed initializer list>)’
  714 |   return(Fdot(octave_value_list{(*xv),(*yv)},1));
      |                                            ^
In file included from /usr/include/octave-8.0.0/octave/builtin-defun-decls.h:8,
                 from nvector_octave.cpp:8:
/usr/include/octave-8.0.0/octave/ovl.h:59:3: note: candidate: ‘template<template<class ...> class OV_Container> octave_value_list::octave_value_list(const OV_Container<octave_value>&)’
   59 |   octave_value_list (const OV_Container<octave_value>& args)
      |   ^~~~~~~~~~~~~~~~~
/usr/include/octave-8.0.0/octave/ovl.h:59:3: note:   template argument deduction/substitution failed:
nvector_octave.cpp:714:44: note:   mismatched types ‘octave_value’ and ‘double’
  714 |   return(Fdot(octave_value_list{(*xv),(*yv)},1));
      |                                            ^
nvector_octave.cpp:714:44: note:   mismatched types ‘octave_value’ and ‘double’
nvector_octave.cpp:714:44: note:   ‘ColumnVector’ is not derived from ‘const OV_Container<octave_value>’
In file included from /usr/include/octave-8.0.0/octave/builtin-defun-decls.h:8,
                 from nvector_octave.cpp:8:
/usr/include/octave-8.0.0/octave/ovl.h:82:3: note: candidate: ‘octave_value_list::octave_value_list(const std::__cxx11::list<octave_value_list>&)’
   82 |   octave_value_list (const std::list<octave_value_list>&);
      |   ^~~~~~~~~~~~~~~~~
/usr/include/octave-8.0.0/octave/ovl.h:82:3: note:   candidate expects 1 argument, 2 provided
/usr/include/octave-8.0.0/octave/ovl.h:81:3: note: candidate: ‘octave_value_list::octave_value_list(const std::__cxx11::list<octave_value>&)’
   81 |   octave_value_list (const std::list<octave_value>&);
      |   ^~~~~~~~~~~~~~~~~
/usr/include/octave-8.0.0/octave/ovl.h:81:3: note:   candidate expects 1 argument, 2 provided
/usr/include/octave-8.0.0/octave/ovl.h:78:3: note: candidate: ‘octave_value_list::octave_value_list(octave_value_list&&)’
   78 |   octave_value_list (octave_value_list&& obj) = default;
      |   ^~~~~~~~~~~~~~~~~
/usr/include/octave-8.0.0/octave/ovl.h:78:3: note:   candidate expects 1 argument, 2 provided
/usr/include/octave-8.0.0/octave/ovl.h:76:3: note: candidate: ‘octave_value_list::octave_value_list(const octave_value_list&)’
   76 |   octave_value_list (const octave_value_list& obj) = default;
      |   ^~~~~~~~~~~~~~~~~
/usr/include/octave-8.0.0/octave/ovl.h:76:3: note:   candidate expects 1 argument, 2 provided
/usr/include/octave-8.0.0/octave/ovl.h:69:3: note: candidate: ‘octave_value_list::octave_value_list(const Cell&)’
   69 |   octave_value_list (const Cell& c)
      |   ^~~~~~~~~~~~~~~~~
/usr/include/octave-8.0.0/octave/ovl.h:69:3: note:   candidate expects 1 argument, 2 provided
/usr/include/octave-8.0.0/octave/ovl.h:62:3: note: candidate: ‘octave_value_list::octave_value_list(const Array<octave_value>&)’
   62 |   octave_value_list (const Array<octave_value>& a)
      |   ^~~~~~~~~~~~~~~~~
/usr/include/octave-8.0.0/octave/ovl.h:62:3: note:   candidate expects 1 argument, 2 provided
/usr/include/octave-8.0.0/octave/ovl.h:55:3: note: candidate: ‘octave_value_list::octave_value_list(const octave_value&)’
   55 |   octave_value_list (const octave_value& tc)
      |   ^~~~~~~~~~~~~~~~~
/usr/include/octave-8.0.0/octave/ovl.h:55:3: note:   candidate expects 1 argument, 2 provided
/usr/include/octave-8.0.0/octave/ovl.h:52:3: note: candidate: ‘octave_value_list::octave_value_list(octave_idx_type, const octave_value&)’
   52 |   octave_value_list (octave_idx_type n, const octave_value& val)
      |   ^~~~~~~~~~~~~~~~~
/usr/include/octave-8.0.0/octave/ovl.h:52:38: note:   no known conversion for argument 1 from ‘ColumnVector’ to ‘octave_idx_type’ {aka ‘long int’}
   52 |   octave_value_list (octave_idx_type n, const octave_value& val)
      |                      ~~~~~~~~~~~~~~~~^
/usr/include/octave-8.0.0/octave/ovl.h:49:12: note: candidate: ‘octave_value_list::octave_value_list(octave_idx_type)’
   49 |   explicit octave_value_list (octave_idx_type n)
      |            ^~~~~~~~~~~~~~~~~
/usr/include/octave-8.0.0/octave/ovl.h:49:12: note:   candidate expects 1 argument, 2 provided
/usr/include/octave-8.0.0/octave/ovl.h:47:3: note: candidate: ‘octave_value_list::octave_value_list()’
   47 |   octave_value_list (void) = default;
      |   ^~~~~~~~~~~~~~~~~
/usr/include/octave-8.0.0/octave/ovl.h:47:3: note:   candidate expects 0 arguments, 2 provided
make: *** [Makefile:37: nvector_octave.o] Error 1
