#!/bin/bash
make plotGenerator;
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
                        echo $path" - N"$n"D"${d//.}"K"$k"T"$t"P"${p//.};
                        rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"/vertices.txt";
                        rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"/edges.txt";
                        rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"/solutionVertices.txt";
                        rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"/solutionEdges.txt";
                        rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"/graph.pdf";
                        ./plotGenerator $path $n ${d//.} $k $t ${p//.};
                        gnuplot -e "path = '${path}'; n='${n}'; d='${d}'; k='${k}'; t='${t}'; p='${p}'; D='${d//.}'; P='${p//.}';" generator2.plt;
                    done
                done
            done
        done
    done
done
make clean;
