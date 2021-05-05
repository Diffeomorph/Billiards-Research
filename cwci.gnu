
set term pos color
set out 'cwc.eps'
set xrange [0.34:0.36]
set yrange [0.938:0.942]
set pointsize 0.001
plot 'cwc.out' using 1:2 with dots
