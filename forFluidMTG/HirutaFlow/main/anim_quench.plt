set term gif animate delay 10
set output "quench.gif"

set pm3d
set pm3d map

set xr [0:8*pi]
set yr [0:2*pi]

set cbr [-8:8]

set palette defined ( 0 '#0fffee',1 '#0090ff', 2 '#000fff',3 '#000090',4 '#ffffff',5 '#7f0000', 6 '#ee0000', 7 '#ff7000', 8 '#ffee00')

n = 0
title(n) = sprintf("Vorticity(t = %d)", 10*n)

do for [n = 0:250] {
    set title title(n) font "Times,20"
    splot "vorticity.dat" index n u 1:2:3 ti "" w pm3d
}

set output