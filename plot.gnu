set title "Energy Vs Iteration Number For System"
set xlabel "Iteration"
set ylabel "Energy"
set key left
plot "energy.txt" using 2:1 title "Energy" with linespoints
