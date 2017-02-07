CPP = g++
CARGS = -Wall -m64 -g -fno-inline -std=c++11

clean: 
	rm -rf geneticAlgorithm solutionAggregator

geneticAlgorithm: geneticAlgorithm.cpp
	$(CPP) $(CARGS) -o geneticAlgorithm geneticAlgorithm.cpp -lpthread -lm

solutionAggregator: solutionAggregator.cpp
	$(CPP) $(CARGS) -o solutionAggregator solutionAggregator.cpp