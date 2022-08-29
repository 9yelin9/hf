#set data style dots
set nokey
set xrange [0: 3.69219]
set yrange [  3.94973 :  7.08587]
set arrow from  0.67371,   3.94973 to  0.67371,   7.08587 nohead
set arrow from  1.45164,   3.94973 to  1.45164,   7.08587 nohead
set arrow from  1.84060,   3.94973 to  1.84060,   7.08587 nohead
set arrow from  2.39069,   3.94973 to  2.39069,   7.08587 nohead
set arrow from  2.86707,   3.94973 to  2.86707,   7.08587 nohead
set xtics (" L "  0.00000," G "  0.67371," X "  1.45164," W "  1.84060," L "  2.39069," K "  2.86707," G "  3.69219)
plot "wannier90_band.dat"
pause mouse key
