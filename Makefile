# compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall

# include paths for Eigen and NLopt
# INCLUDES = -I./includes/eigen -I/usr/local/include #make sure these point to the correct directories
# INCLUDES = -I./includes/eigen -I/mnt/c/Users/dwamai/nloptinstall/include #lab pc
INCLUDES = -I./includes/eigen -I/opt/homebrew/Cellar/nlopt/2.7.1/include #mac

# library paths for NLopt
# LDFLAGS = -L/usr/local/lib
# LDFLAGS =  -L/mnt/c/Users/dwamai/nloptinstall/lib #lab pc
LDFLAGS =  -L/opt/homebrew/Cellar/nlopt/2.7.1/lib #mac
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
