CXX=g++
CXXFLAGS= -O2 -std=c++11 -g
OBJFILES= main_steane.o random.o circuit.o extractor.o toric_code.o decoder.o
TARGET = simulator

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJFILES) 

clean:
	rm -f $(OBJFILES) $(TARGET) *~
