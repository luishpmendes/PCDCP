#!/bin/bash
rm -rf output;
mkdir -p output;
make grasp;
iterationLimit=10000;
for n in 10 20 50 100 250
do
    for d in 0.3 0.5 0.7
    do
        for k in 0 5 10 20
        do
            for t in 0 1
            do
                for p in 0.1 0.5 1.0
                do
                    for i in {0..9}
                    do
                        for alpha in 0.3 0.5 0.7
                        do
                            echo "N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"A"${alpha//.};
                            mkdir -p "output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"A"${alpha//.};
                            ./grasp $iterationLimit $alpha < "../input/instanceN"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i".in" > "output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"A"${alpha//.}"/result.out";
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
