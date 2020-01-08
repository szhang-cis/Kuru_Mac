#GNUPLOT
set encoding utf8
set terminal postscript eps enhanced 'Helvetica' 20 color
set output 'outeradio.eps'
#set multiplot
#set size 0.33,0.33
#
#set origin 0.0,0.5
set key top left
set title 'Outer Radius'
set xlabel 'Time[Months]'
set ylabel 'Radius[mm]'
plot [0:*][*:*] 'coordo.dat' u 1:(($2**2+$3**2)**0.5*1000) t '' w l
#
#
#
#
#
