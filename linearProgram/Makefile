INC = /opt/gurobi652/linux64/include/
CPP = g++
CARGS = -Wall -m64 -g -O0 -fno-inline
CPPLIB = -L/opt/gurobi652/linux64/lib/ -lgurobi_c++ -lgurobi65

clean: 
	rm -rf linearProgram

linearProgram: linearProgram.cpp
	$(CPP) $(CARGS) -o linearProgram linearProgram.cpp -I$(INC) $(CPPLIB) -lpthread -lm
