INC = /home/lmendes/gurobi652/linux64/include/
CPP = g++
CARGS = -Wall -m64 -O3 -fno-inline -std=c++11
CPPLIB = -L/home/lmendes/gurobi652/linux64/lib/ -lgurobi_c++ -lgurobi65

clean: 
	rm -rf linearProgram solutionAggregator

linearProgram: linearProgram.cpp
	$(CPP) $(CARGS) -o linearProgram linearProgram.cpp -I$(INC) $(CPPLIB) -lpthread -lm

solutionAggregator: solutionAggregator.cpp
	$(CPP) $(CARGS) solutionAggregator.cpp -o solutionAggregator
