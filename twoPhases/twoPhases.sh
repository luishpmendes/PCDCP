#!/bin/bash
rm -rf output;
mkdir -p output;
make phase1;
make phase2;
timeLimit1=30;
timeLimit2=30;
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
                    echo "phase1 - N"$n"D"${d//.}"K"$k"T"$t"I"$i;
                    mkdir -p "output/N"$n"D"${d//.}"K"$k"T"$t"I"$i;
                    ./phase1 $i $timeLimit1 < "../input/instanceN"$n"D"${d//.}"K"$k"T"$t"I"$i".in" > "output/N"$n"D"${d//.}"K"$k"T"$t"I"$i"/result.tmp";
                done
            done
        done
    done
done
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
                    echo "phase2 - N"$n"D"${d//.}"K"$k"T"$t"I"$i;
                    mkdir -p "output/N"$n"D"${d//.}"K"$k"T"$t"I"$i;
                    ./phase2 $i $timeLimit2 < "output/N"$n"D"${d//.}"K"$k"T"$t"I"$i"/result.tmp" > "output/N"$n"D"${d//.}"K"$k"T"$t"I"$i"/result.out";
                done
            done
        done
    done
done
make solutionAggregator;
./solutionAggregator > "output/solution.csv";
make clean;
