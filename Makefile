CXX=g++
CXXFLAGS= -O2 -std=c++11
OBJFILES= main.o random.o circuit.o extractor.o toric_code.o
TARGET = simulator

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJFILES) 

clean:
	rm -f $(OBJFILES) $(TARGET) *~
