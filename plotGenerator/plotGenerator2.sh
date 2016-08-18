#!/bin/bash
make plotGenerator;
path="linearProgram";
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
                    echo "N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.};
                    rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.}"/vertices.txt";
                    rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.}"/edges.txt";
                    rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.}"/solutionVertices.txt";
                    rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.}"/solutionEdges.txt";
                    rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.}"/graph.pdf";
                    ./plotGenerator $path $n ${d//.} $k $t $r ${p//.};
                    gnuplot -e "path = '${path}'; n='${n}'; d='${d}'; k='${k}'; t='${t}'; r='${r}'; p='${p}'; D='${d//.}'; P='${p//.}';" generator2.plt;
                done
            done
        done
    done
done
make clean;
