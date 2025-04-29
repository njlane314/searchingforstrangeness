CXX = g++
CXXFLAGS = -Wall -Wextra -O2 -std=c++11
ROOTFLAGS = $(shell root-config --cflags --libs)
SRC = SelectionAnalyser.cpp
EXE = analyse

all: $(EXE)

$(EXE): $(SRC)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $(EXE) $(SRC)

clean:
	rm -f $(EXE)