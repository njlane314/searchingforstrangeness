CXX = g++
CXXFLAGS = -Wall -O2 `root-config --cflags`
LDFLAGS = `root-config --libs`

SOURCES = analyser.C effpur.C

OBJECTS = $(SOURCES:.C=.o)

EXECUTABLES = analyser effpur

all: $(EXECUTABLES)

analyser: analyser.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

effpur: effpur.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.C
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLES)