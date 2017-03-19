#!/bin/bash

cd linearProgram;

rm -rf output;
mkdir -p output;
make linearProgram;
timeLimit=60;
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

cd ..;

cd twoPhases;

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

cd ..;

cd grasp;

rm -rf output;
mkdir -p output;
make grasp;
iterationLimit=60;
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
                    for alpha in 0.3 0.5 0.7
                    do
                        echo "N"$n"D"${d//.}"K"$k"T"$t"I"$i"A"${alpha//.};
                        mkdir -p "output/N"$n"D"${d//.}"K"$k"T"$t"I"$i"A"${alpha//.};
                        ./grasp $i $iterationLimit $alpha < "../input/instanceN"$n"D"${d//.}"K"$k"T"$t"I"$i".in" > "output/N"$n"D"${d//.}"K"$k"T"$t"I"$i"A"${alpha//.}"/result.out";
                    done
                done
            done
        done
    done
done
make solutionAggregator;
./solutionAggregator > "output/solution.csv";
make clean;

cd ..;

cd geneticAlgorithm;

rm -rf output;
mkdir -p output;
make geneticAlgorithm;
timeLimit=60;
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
                    for ps in 10 50 100
                    do
                        for mr in 0.1 0.2 0.3
                        do
                            echo "N"$n"D"${d//.}"K"$k"T"$t"I"$i"PS"$ps"MR"${mr//.};
                            mkdir -p "output/N"$n"D"${d//.}"K"$k"T"$t"I"$i"PS"$ps"MR"${mr//.};
                            ./geneticAlgorithm $i $timeLimit $ps $mr < "../input/instanceN"$n"D"${d//.}"K"$k"T"$t"I"$i".in" > "output/N"$n"D"${d//.}"K"$k"T"$t"I"$i"PS"$ps"MR"${mr//.}"/result.out";
                        done
                    done
                done
            done
        done
    done
done
make solutionAggregator;
./solutionAggregator > "output/solution.csv";
make clean;
