GNUPLOT BOOK:

http://www.gnuplot.info/gnuplot_cvs.pdf


FILLED CONTOUR PLOT: 

set pm3d
unset surface  
set view map
set key outside
set pm3d interpolate 20,20 # interpolate the color
set palette rgbformulae 22,13,10
set xlabel 'x, m'
set ylabel 'y, m'
set format x '%.2f'
set format y '%.2f'
set format z '%.2f'
splot 'post.01' using 1:2:3 notitle with lines lt 1




CONTOUR PLOT WITH LABELS

unset surface  
set view map
set contour
set key outside
set cntrparam cubicspline  # smooth out the lines
set cntrparam levels 20    # sets the num of contour lines
set cntrlabel start 100 interval -1 font ",7" # set spacing between labels along contour line
set xlabel 'x, m'
set ylabel 'y, m'
set format x '%.2f'
set format y '%.2f'
set format z '%.2f'
#set xrange [0:1]
#set yrange [0:0.1]
splot 'post.01' using 1:2:7 notitle with labels, 'post.01' using 1:2:7 notitle with lines lt 1
