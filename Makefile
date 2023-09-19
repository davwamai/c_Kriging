# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall

# Include paths for Eigen and NLopt
INCLUDES = -I./includes/eigen -I/usr/local/include

# Library paths for NLopt
LDFLAGS = -L/usr/local/lib
LDLIBS = -lnlopt -lm

# Source and object files
SRCS = OrdinaryKriging.cpp UniversalKriging.cpp Variograms.cpp main.cpp
OBJS = $(SRCS:.cpp=.o)

# Target executable
TARGET = kriging

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LDLIBS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(TARGET)
