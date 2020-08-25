CXX=g++
CXXFLAGS= -O2 -std=c++11
OBJFILES= main.o
TARGET = ToricCode

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJFILES) 

clean:
	rm -f $(OBJFILES) $(TARGET) *~
