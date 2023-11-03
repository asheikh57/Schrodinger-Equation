set datafile separator ","

set xrange [-5:5]
set yrange [0:10]
plot "energy1.csv" title 'ground state' with lines linestyle 1, \
     "energy2.csv" title 'first excited state' with lines linestyle 2, \
     "energy3.csv" title 'second excited state' with lines linestyle 3, \
     "energy4.csv" title 'third excited state' with lines linestyle 4, \
     "energy5.csv" title 'fourth excited state' with lines linestyle 5, \
     "energy6.csv" title 'fifth excited state' with lines linestyle 6, \
     0.25*x*x*x*x title 'potential' with lines linestyle 8

