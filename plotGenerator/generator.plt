set terminal postscript eps color font "times"
set output 'plot.eps'
unset key
set size ratio -1 1,1
set title "Instancia do PCDCP"
set xrange [0:1]
set yrange [0:1]
plot 'vertices.txt' u ($1):($2):($3) with points pt 7 ps variable lc rgb "black", 'edges.txt' with lines lt 0.5 lw 0.5 lc rgb "light-gray", 'solutionEdges.txt' with lines lt 1 lw 1 lc rgb "forest-green", 'solutionVertices.txt' with circles lt 1 lw 0.5 rgb "red";
