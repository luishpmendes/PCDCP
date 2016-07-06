if (!exists("n")) n='5'
if (!exists("d")) d='03'
if (!exists("k")) k='0'
if (!exists("t")) t='0'
if (!exists("path")) path='linearProgram'
set terminal postscript eps color font "times"
#set terminal postscript eps color font "times" fontscale 1.0 size 500, 300 
set output '../'.path.'/N'.n.'D'.d.'K'.k.'T'.t.'/plot.eps'
unset key
set size ratio -1 1,1
set title "Instancia do PCDCP (n=".n."; d=".d."; k=".k."; t=".t.")"
#set title sprintf("Instancia do PCDCP (n=%d; d=%f; k=%d; t=%d)", n, d, k, t)
#set title "Instancia do DVRP (n=250; k=0,15)"
set xrange [0:1]
set yrange [0:1]
plot 'vertices.txt' u ($1):($2):($3) with points pt 7 ps variable lc rgb "black", 'edges.txt' with lines lt 0.5 lw 0.5 lc rgb "light-gray", 'solutionEdges.txt' with lines lt 1 lw 1 lc rgb "forest-green", 'solutionVertices.txt' with circles lt 1 lw 0.5 rgb "red"

#plot 'dom_set.txt' with circles lt 1 lw 0.5 lc rgb "red", 'edges.txt' with lines lt 1 lw 2 lc rgb "forest-green", 'points.txt' u ($1):($2) with points pt 7 ps 0.5 lc rgb "black", 'dom_set.txt' u ($1):($2) with points pt 7 lc rgb "red", 'depot.txt' u ($1):($2) with points pt 7 lc rgb "blue"
