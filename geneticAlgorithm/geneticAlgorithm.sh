#!/bin/bash
rm -rf output;
mkdir -p output;
make geneticAlgorithm;
timeLimit=10;
for n in 10 20 50 100
do
    for d in 0.3 0.5 0.7
    do
        for k in 0 10 20
        do
            for t in 0 1
            do
                for p in 0.2
                do
                    for i in 0
                    do
                        for populationSize in 10 50 100
                        do
                            for mutationRate in 0.1 0.2 0.3
                            do
                                echo "N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i;
                                mkdir -p "output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"PS"$populationSize"MR"${mutationRate//.};
                                ./geneticAlgorithm $i $timeLimit $populationSize $mutationRate < "../input/instanceN"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i".in" > "output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"PS"$populationSize"MR"${mutationRate//.}"/result.out";
                            done
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
