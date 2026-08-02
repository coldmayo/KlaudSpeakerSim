TARGET = speaker_sim

CXX = g++
CXXFLAGS = -O2 -std=c++20 -Wall

PYTHON_CFLAGS := $(shell python3-config --cflags)
PYTHON_LDFLAGS := $(shell python3-config --ldflags --embed 2>/dev/null || python3-config --ldflags)

NUMPY_INCLUDE := $(shell python3 -c "import numpy; print(numpy.get_include())" 2>/dev/null)

INCLUDES = -I.
ifneq ($(NUMPY_INCLUDE),)
    INCLUDES += -I$(NUMPY_INCLUDE)
endif

CXXFLAGS += $(PYTHON_CFLAGS) $(INCLUDES)
LDFLAGS  += $(PYTHON_LDFLAGS)

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
