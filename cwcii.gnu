
set term pos color
set out 'cwc.eps'
plot 'cwc.out' using 3:4 with lines, 'cwc.plot' with dots
