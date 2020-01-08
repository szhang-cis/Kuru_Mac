#GNUPLOT
set encoding utf8
set terminal postscript eps enhanced 'Helvetica' 20 color
set output 'thickness.eps'
#set multiplot
#set size 0.33,0.33
#
#set origin 0.0,0.5
set key bottom left
set title 'Thickness'
set xlabel 'Time[Months]'
set ylabel 'Thickness Wall[mm]'
plot [0:*][*:*] 'thickn.dat' u 1:($2*1000) t '' w l
#
#
#
#
#
