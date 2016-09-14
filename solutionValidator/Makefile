CPPC = g++
CFLAGS = -Wall
DEBUGFLAGS = -O3 -fno-inline -std=c++11

COMPILER = $(CPPC)
CFLAGS += $(DEBUGFLAGS)

clean: 
	rm -rf solutionValidator validationAggregator

solutionValidator: solutionValidator.cpp
	$(COMPILER) $(CFLAGS) solutionValidator.cpp -o solutionValidator

validationAggregator: validationAggregator.cpp
	$(COMPILER) $(CFLAGS) validationAggregator.cpp -o validationAggregator
