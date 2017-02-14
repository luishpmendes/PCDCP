#! /usr/bin/gnuplot
set terminal pdf color font "times"
unset key
set size ratio -1 1,1
set xrange [0:1]
set yrange [0:1]
if (path eq "grasp") {
    set output "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."A".A."/result.pdf"
} else if (path eq "geneticAlgorithm") {
    set output "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."PS".ps."MR".MR."/result.pdf"
} else {
    set output "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."/result.pdf"
}
set title "Instancia do PCDCP (N = ".n."; D = ".d."; K = ".k."; T = ".t."; P = ".p.")"
if (path eq "grasp") {
    vertices = "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."A".A."/vertices.txt"
} else if (path eq "geneticAlgorithm") {
    vertices = "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."PS".ps."MR".MR."/vertices.txt"
} else {
    vertices = "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."/vertices.txt"
}
if (path eq "grasp") {
    solutionVertices = "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."A".A."/solutionVertices.txt"
} else if (path eq "geneticAlgorithm") {
    solutionVertices = "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."PS".ps."MR".MR."/solutionVertices.txt"
} else {
    solutionVertices = "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."/solutionVertices.txt"
}
if (path eq "grasp") {
    solutionEdges = "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."A".A."/solutionEdges.txt"
} else if (path eq "geneticAlgorithm") {
    solutionEdges = "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."PS".ps."MR".MR."/solutionEdges.txt"
} else {
    solutionEdges = "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."/solutionEdges.txt"
}
plot vertices using ($1):($2):($3) with points pointtype 7 pointsize variable linecolor rgb "black", \
	 solutionEdges with lines linetype 1 linewidth 1 linecolor rgb "forest-green", \
	 solutionVertices with circles linetype 1 linewidth 0.5 linecolor rgb "red"
