#GNUPLOT
set encoding utf8
set terminal postscript eps enhanced 'Helvetica' 10 color
set output 'analitic.eps'
set multiplot
set size 0.5

set origin 0.0,0.5
set title 'Circumferential Stress'
set xlabel 'Steps'
set ylabel 'Cauchy Stress[kPa]'
plot 'analitic.dat' u 1:($5/1000) t 'Fit' w l

set origin 0.0,0.0
set title 'Longitudinal Stress'
set xlabel 'Steps'
set ylabel 'Cauchy Stress[kPa]'
plot 'analitic.dat' u 1:($6/1000) t 'Fit' w l

set origin 0.5,0.5
set title 'Radius'
set xlabel 'Steps'
set ylabel 'Radius[mm]'
plot 'analitic.dat' u 1:(23.0*$3) t '' w l

set origin 0.5,0.0
set title 'Thickness'
set xlabel 'Steps'
set ylabel 'Thickness'
plot 'analitic.dat' u 1:(2.0/$3) t '' w l
