CXX = g++ 
CXXFLAGS = -Wall -g

main: main.o network.o node.o
	$(CXX) $(CXXFLAGS) -o run main.o network.o node.o

main.o: network.h node.h main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp

network.o: network.h node.h network.cpp
	$(CXX) $(CXXFLAGS) -c network.cpp

node.o: node.h node.cpp
	$(CXX) $(CXXFLAGS) -c node.cpp

clean:
	rm -rf *o run
