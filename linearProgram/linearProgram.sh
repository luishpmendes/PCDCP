#!/bin/bash
rm -rf output;
mkdir -p output;
make linearProgram;
timeLimit=360.0;
for n in 50 100 200
do
    for d in 0.3 0.5 0.7
    do
        for k in 0 10 20
        do
            for t in 0 1
            do
                for i in 0 1 2
                do
                    echo "N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i;
                    mkdir -p "output/N"$n"D"${d//.}"K"$k"T"$t"I"$i;
                    ./linearProgram $i $timeLimit < "../input/instanceN"$n"D"${d//.}"K"$k"T"$t"I"$i".in" > "output/N"$n"D"${d//.}"K"$k"T"$t"I"$i"/result.out";
                done
            done
        done
    done
done
make solutionAggregator;
./solutionAggregator > "output/solution.csv";
make clean;
