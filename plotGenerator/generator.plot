set terminal postscript eps color font "times"
#set terminal postscript eps color font "times" fontscale 1.0 size 500, 300 
set output 'circles.eps'
unset key
set size ratio -1 1,1
set title "Instancia do DVRP (n=250; k=0,15)" 
set xrange [0:1]
set yrange [0:1]
plot 'points.txt' u ($1):($2):($3) with points pt 7 ps variable lc rgb "black", 'edges.txt' with lines lt 1 lw 1 lc rgb "forest-green"

#plot 'dom_set.txt' with circles lt 1 lw 0.5 lc rgb "red", 'edges.txt' with lines lt 1 lw 2 lc rgb "forest-green", 'points.txt' u ($1):($2) with points pt 7 ps 0.5 lc rgb "black", 'dom_set.txt' u ($1):($2) with points pt 7 lc rgb "red", 'depot.txt' u ($1):($2) with points pt 7 lc rgb "blue"
