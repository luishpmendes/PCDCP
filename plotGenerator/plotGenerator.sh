#!/bin/bash
make plotGenerator;
for path in "linearProgram" "twoPhases" "grasp"
do
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
                            echo $path" - N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i;
                            rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"/vertices.txt";
                            rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"/edges.txt";
                            rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"/solutionVertices.txt";
                            rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"/solutionEdges.txt";
                            rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"/result.pdf";
                            ./plotGenerator $path $n ${d//.} $k $t ${p//.} $i;
                            gnuplot -e "path = '${path}'; n='${n}'; d='${d}'; k='${k}'; t='${t}'; p='${p}'; D='${d//.}'; P='${p//.}'; i='$i';" generator.plt;
                        done
                    done
                done
            done
        done
    done
done
make clean;
