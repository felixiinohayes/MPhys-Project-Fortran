ef=  4.19601857
set xrange [ -1 : 1]
set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 1 size 5.00in, 5.00in
set output "band.pdf"
set border
set xtics
set ytics
set encoding iso_8859_1
set size ratio 0 1.0,1.0
set yrange [3.5: 8.5]
unset key
set mytics 2
set parametric
set trange [-10:10]
set multiplot
#plot "band.dat" every 4 u 1:($2-ef):(column(3)*4) with points pt 7 ps variable lc rgb "royalblue"
#plot "band.dat" every 4 u 1:($2-ef):(column(4)*4) with points pt 7 ps variable lc rgb "light-red"
#plot "band.dat" every 4 u 1:($2-ef):(column(5)*4) with points pt 7 ps variable lc rgb "forest-green"
plot "band.dat" u 1:($2-ef) with l lt 1 lw 3.5 lc rgb "black"
unset multiplot

