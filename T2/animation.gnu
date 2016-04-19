reset
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 2 # --- blue
unset key
set border 0
unset tics
set view 342,0
box = `head -n1 out.txt | awk '{print $6}'`
set xrange [0:box]
set yrange [0:box]
stats "out.txt"
do for [ii=1:STATS_blocks-1:450] {
	    plot "out.txt" index (ii-1) using 1:2 w p ls 1
		pause 0.1	
} 
