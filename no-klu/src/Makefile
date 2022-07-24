CXX = g++
CXXFLAGS = -g 
SOURCES1 = kry_example.c nvector_octave.cpp
SOURCES2 = test_nvector_octave.cpp test_nvector.cpp nvector_octave.cpp
DEPS = nvector_octave.h test_nvector.h
OBJECTS1 = $(SOURCES1:.cpp=.o)
OBJECTS2 = $(SOURCES2:.cpp=.o)
EXECUTABLE1 = app 
EXECUTABLE2 = test 

# include+lib paths to be defined by the user temporarily
KLU_INC_PATH = /usr/include/suitesparse
SUN_INC_PATH = /usr/local/include
OCT_INC_PATH = /usr/include/octave-8.0.0/octave

KLU_LIB_PATH = /usr/lib/x86_64-linux-gnu
SUN_LIB_PATH = /usr/local/lib
OCT_LIB_PATH = /usr/lib/octave/8.0.0
LDFLAGS =-I$(KLU_INC_PATH) -I$(SUN_INC_PATH) -I$(OCT_INC_PATH) -L$(KLU_LIB_PATH) -L$(SUN_LIB_PATH) -L$(OCT_LIB_PATH)
# -I/usr/local/include -I/usr/include/octave-8.0.0/octave -L/usr/local/lib -L/usr/lib/octave/8.0.0 -lsundials_generic -lsundials_nvecserial -lsundials_ida -lsundials_sunlinsolspgmr -loctave

LIBS = -lsundials_generic -lsundials_nvecserial -lsundials_ida -lsundials_sunlinsolklu -loctave

app: $(EXECUTABLE1)
	@echo Make of example has finished successfully.

test: $(EXECUTABLE2)
	@echo Make of test finished successfully.

$(EXECUTABLE1) : $(OBJECTS1)
	$(CXX) $(CXXFLAGS) $(OBJECTS1) $(LDFLAGS) $(LIBS) -o $@

$(EXECUTABLE2) : $(OBJECTS2)
	$(CXX) $(CXXFLAGS) $(OBJECTS2) $(LDFLAGS) $(LIBS) -o $@

%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LIBS) -c $*.cpp 

clean:
	$(RM) *.o *~ $(EXECUTABLE1) $(EXECUTABLE2)

# requires the export LD_LIBRARY_PATH=/usr/lib/octave/8.0.0 variable to be set on my device
# to access shared libraries as octave is built from source