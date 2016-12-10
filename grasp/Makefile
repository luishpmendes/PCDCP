CPP = g++
CARGS = -Wall -m64 -g -fno-inline -std=c++11

clean: 
	rm -rf grasp solutionAggregator

grasp: grasp.cpp
	$(CPP) $(CARGS) -o grasp grasp.cpp -lpthread -lm

solutionAggregator: solutionAggregator.cpp
	$(CPP) $(CARGS) -o solutionAggregator solutionAggregator.cpp