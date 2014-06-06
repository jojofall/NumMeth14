###gnuplot commands:
unset key
set size 0.26,0.5
set pm3d map
#set dgrid3d 254,333
set palette defined (1 'white', 2 'blue', 3 'red',4 'yellow' ,5 'pink' )
unset xtics
unset ytics
unset colorbox
#set xrange [1:254]
#set yrange [1:333]
#splot "CorrectSegImage.dat"  u 1:(-$2):3 t ""
splot "saSweep.dat"  u 1:(-$2):3 t ""


 se term postscript eps enhanced color
 se output "new.ps"
 replot
