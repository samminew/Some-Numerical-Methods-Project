#TASK 1
set title "Iteration vrs Root Iterative Methods"
set xlabel "Iteration"
set ylabel "f(x)"
set yrange [-10:15]
plot "data.dat" using 1:2 title "Bisection" with linespoints lc "blue" lw 2 pt 7, "data.dat" using 1:3 title "Secant" with linespoints lc "orange" lw 2 pt 7, "data.dat" using 1:4 title "NewtonR" with linespoints lc "green" lw 2 pt 7, "data.dat" using 1:5 title "Fixed point" with linespoints lc "red" lw 2 pt 7

pause -1

#TASK 2
set title "Estimated temp and Exact temp Comparison"
set xlabel "Time"
set xrange [0:500]
set ylabel "Temperature"
set yrange [-1100:1700]
plot "euler.dat" using 1:2 title "Theta(480)" with linespoints lc "red" lw 2 pt 6, "euler.dat" using 1:3 title "Exact" with linespoints lc "purple" lw 2 pt 6

pause -1

