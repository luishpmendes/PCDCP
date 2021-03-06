#!/bin/bash
make solutionValidator;
for path in "linearProgram" "twoPhases"
do
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
                        echo $path" - N"$n"D"${d//.}"K"$k"T"$t"I"$i;
                        rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"I"$i"/validation.txt";
                        ./solutionValidator $path $n ${d//.} $k $t $i > "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"I"$i"/validation.txt";
                    done
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
                    for alpha in 0.3 0.5 0.7
                    do
                        echo "grasp - N"$n"D"${d//.}"K"$k"T"$t"I"$i"A"${alpha//.};
                        rm -f "../grasp/output/N"$n"D"${d//.}"K"$k"T"$t"I"$i"A"${alpha//.}"/validation.txt";
                        ./solutionValidator grasp $n ${d//.} $k $t $i ${alpha//.} > "../grasp/output/N"$n"D"${d//.}"K"$k"T"$t"I"$i"A"${alpha//.}"/validation.txt";
                    done
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
                    for ps in 10 50 100
                    do
                        for mr in 0.1 0.2 0.3
                        do
                            echo "geneticAlgorithm - N"$n"D"${d//.}"K"$k"T"$t"I"$i"PS"$ps"MR"${mr//.};
                            rm -f "../geneticAlgorithm/output/N"$n"D"${d//.}"K"$k"T"$t"I"$i"PS"$ps"MR"${mr//.}"/validation.txt";
                            ./solutionValidator geneticAlgorithm $n ${d//.} $k $t $i $ps ${mr//.} > "../geneticAlgorithm/output/N"$n"D"${d//.}"K"$k"T"$t"I"$i"PS"$ps"MR"${mr//.}"/validation.txt";
                        done
                    done
                done
            done
        done
    done
done
make validationAggregator;
for path in "linearProgram" "twoPhases" "grasp" "geneticAlgorithm"
do
    ./validationAggregator $path > "../"$path"/output/validation.csv";
done
make clean;
