CPPC = g++
CFLAGS = -Wall
DEBUGFLAGS = -O3 -fno-inline -std=c++11

COMPILER = $(CPPC)
CFLAGS += $(DEBUGFLAGS)

clean: 
	rm -rf plotGenerator

plotGenerator: plotGenerator.cpp
	$(COMPILER) $(CFLAGS) plotGenerator.cpp -o plotGenerator
