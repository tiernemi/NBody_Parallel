reset
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 2 # --- blue
unset key
set border 0
unset tics
set view 342,0
box = `head -n1 out.txt0 | awk '{print $6}'`
cmdStr = "ls -l out* | wc -l"
numFiles = system(cmdStr)
stats "out.txt0"
set xrange [0:box]
set yrange [0:box]
do for [ii=1:STATS_blocks-1:15] {
	    plot for [i=0:numFiles-1] 'out.txt'.i index (ii-1) using 1:2 w p pt 7 lc i ps 2
		pause 0.1	
} 



