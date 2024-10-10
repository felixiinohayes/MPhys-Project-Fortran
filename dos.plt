# Set the output file type
set terminal pngcairo enhanced font "DejaVu,30" fontscale 1 size 2000, 2000
# Set the output file name
set output 'dos.png'

# Now plot the data with lines and points
plot '< sort -nk1 dos_layer_test.dat' using 1:2 smooth csplines w lines ls 1 title 'y1', \
     '' using 1:3 smooth csplines w lines title 'y2', \
     '' using 1:4 smooth csplines w lines title 'y3', \
    #  '' using 1:5 smooth csplines w lines title 'y4', \
    #  '' using 1:6 smooth csplines w lines title 'y5', \
    #  '' using 1:7 smooth csplines w lines title 'y6'