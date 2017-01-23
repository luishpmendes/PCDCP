#!/bin/bash
make solutionValidator;
for path in "linearProgram" "twoPhases"
do
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
                        for alpha in 0.3 0.5 0.7
                        do
                            echo "grasp - N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"A"${alpha//.};
                            rm -f "../grasp/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"A"${alpha//.}"/validation.txt";
                            ./solutionValidator $path $n ${d//.} $k $t ${p//.} $i ${alpha//.} > "../grasp/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"A"${alpha//.}"/validation.txt";
                        done
                    done
                done
            done
        done
    done
done
make validationAggregator;
for path in "linearProgram" "twoPhases" "grasp"
do
    ./validationAggregator $path > "../"$path"/output/validation.csv";
done
make clean;
