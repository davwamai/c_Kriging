# compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall

# include paths for Eigen and NLopt
INCLUDES = -I./includes/eigen -I/usr/local/include #make sure these point to the correct directories

# library paths for NLopt
LDFLAGS = -L/usr/local/lib
LDLIBS = -lnlopt -lm

# source and object files
SRCS = OrdinaryKriging.cpp UniversalKriging.cpp Variograms.cpp main.cpp
OBJS = $(SRCS:.cpp=.o)

# target executable
TARGET = kriging

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LDLIBS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(TARGET)
