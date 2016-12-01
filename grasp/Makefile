CPP = g++
CARGS = -Wall -m64 -g -fno-inline -std=c++11

clean: 
	rm -rf grasp

grasp: grasp.cpp
	$(CPP) $(CARGS) -o grasp grasp.cpp -lpthread -lm
