reset
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 2 # --- blue
unset key
set border 0
unset tics
set view 342,0
set xrange [-2:2]
set yrange [-2:2]

stats "out.txt"
do for [ii=1:STATS_blocks-1:100] {
	    plot "out.txt" index (ii-1) using 1:2 w p ls 1
		pause 0.1	
} 
