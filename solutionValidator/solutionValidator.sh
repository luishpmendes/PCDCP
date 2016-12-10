#!/bin/bash
make solutionValidator;
for path in "linearProgram" "twoPhases" "grasp"
do
    for n in 10 20 50 100
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
                            echo $path" - N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i;
                            rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"/validation.txt";
                            ./solutionValidator $path $n ${d//.} $k $t ${p//.} $i > "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"/validation.txt";
                        done
                    done
                done
            done
        done
    done
done
make validationAggregator;
for path in "linearProgram" "twoPhases"
do
    ./validationAggregator $path > "../"$path"/output/validation.txt";
done
make clean;
