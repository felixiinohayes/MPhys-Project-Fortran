ef= 17.22439820
#set xrange [ -0.12 : 0.12]
set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 0.5 size 5.00in, 5.00in
set output "band.pdf"
set title "B =       0.500T (along y-axis) \n beta_0 =       3.000 (along y-axis)"
set border
set xtics
set ytics
set encoding iso_8859_1
set xlabel "k_y"
set ylabel "E"
set size ratio 0 1.0,1.0
unset key
set mytics 2
set parametric
set trange [-10:10]
set multiplot
plot "band.dat" u 1:($2-ef) with l lt 1 lw 1 lc rgb "black"
unset multiplot
