CXX = g++
CXXFLAGS = -g 
SOURCES = example.cpp nvector_col.cpp
DEPS = nvector_col.h
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = app

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

all: $(EXECUTABLE)
	@echo Make has finished successfully.

$(EXECUTABLE) : $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS) $(LIBS) -o $@

%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LIBS) -c $*.cpp 

clean:
	$(RM) *.o *~ $(EXECUTABLE)

# requires the export LD_LIBRARY_PATH=/usr/lib/octave/8.0.0 variable to be set on my device
# to access shared libraries as octave is built from source