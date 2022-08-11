CXX = g++
CXXFLAGS = -g
SOURCES1 = kry_example.c nvector_octave.cpp
SOURCES2 = test_nvector_octave.cpp test_nvector.cpp nvector_octave.cpp
SOURCES3 = klu_example.c nvector_octave.cpp octmatrix_sparse.cpp octlinsol_gen.cc
SOURCES4 = test_sunmatrix_octave.cpp test_sunmatrix.cpp octmatrix_sparse.cpp nvector_octave.cpp
SOURCES5 = test_sunlinsol_oct.cpp test_sunlinsol.cpp octlinsol_gen.cc octmatrix_sparse.cpp nvector_octave.cpp
OBJECTS1 = $(SOURCES1:.cpp=.o)
OBJECTS2 = $(SOURCES2:.cpp=.o)
OBJECTS3 = $(SOURCES3:.cpp=.o)
OBJECTS4 = $(SOURCES4:.cpp=.o)
OBJECTS5 = $(SOURCES5:.cpp=.o)
EXECUTABLE1 = app 
EXECUTABLE2 = test 
EXECUTABLE3 = app2 
EXECUTABLE4 = test2 
EXECUTABLE5 = test3 

# include+lib paths to be defined by the user temporarily
KLU_INC_PATH = /usr/include/suitesparse
SUN_INC_PATH = /usr/local/include
OCT_INC_PATH = /usr/include/octave-8.0.0/octave

PT_LIB_PATH = /usr/lib/i386-linux-gnu
KLU_LIB_PATH = /usr/lib/x86_64-linux-gnu
SUN_LIB_PATH = /usr/local/lib
OCT_LIB_PATH = /usr/lib/octave/8.0.0
LDFLAGS = -I$(KLU_INC_PATH) -I$(SUN_INC_PATH) -I$(OCT_INC_PATH) -L$(KLU_LIB_PATH) -L$(SUN_LIB_PATH) -L$(OCT_LIB_PATH) -Wl,-Bsymbolic 
# -I/usr/local/include -I/usr/include/octave-8.0.0/octave -L/usr/local/lib -L/usr/lib/octave/8.0.0 -lsundials_generic -lsundials_nvecserial -lsundials_ida -lsundials_sunlinsolspgmr -loctave

LIBS = -lsundials_generic -lsundials_nvecserial -lsundials_ida -lklu -lsundials_sunlinsolklu -loctave -loctinterp -lm -lpthread

app: $(EXECUTABLE1)
	@echo Make of example has finished successfully.

test: $(EXECUTABLE2)
	@echo Make of test finished successfully.

app2 : $(EXECUTABLE3)
	@echo Make of klu example finished successfully.

test2 : $(EXECUTABLE4)
	@echo Make of matrix test finished successfully.

test3 : $(EXECUTABLE5)
	@echo Make of matrix test finished successfully.

$(EXECUTABLE1) : $(OBJECTS1)
	$(CXX) $(CXXFLAGS) $(OBJECTS1) $(LDFLAGS) $(LIBS) -o $@

$(EXECUTABLE2) : $(OBJECTS2)
	$(CXX) $(CXXFLAGS) $(OBJECTS2) $(LDFLAGS) $(LIBS) -o $@

$(EXECUTABLE3) : $(OBJECTS3)
	$(CXX) $(CXXFLAGS) $(OBJECTS3) $(LDFLAGS) $(LIBS) -o $@

$(EXECUTABLE4) : $(OBJECTS4)
	$(CXX) $(CXXFLAGS) $(OBJECTS4) $(LDFLAGS) $(LIBS) -o $@

$(EXECUTABLE5) : $(OBJECTS5)
	$(CXX) $(CXXFLAGS) $(OBJECTS5) $(LDFLAGS) $(LIBS) -o $@

%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LIBS) -c $*.cpp 

clean:
	$(RM) *.o *~ $(EXECUTABLE1) $(EXECUTABLE2) $(EXECUTABLE3) $(EXECUTABLE4)

# requires the export LD_LIBRARY_PATH=/usr/lib/octave/8.0.0 variable to be set on my device
# to access shared libraries as octave is built from source