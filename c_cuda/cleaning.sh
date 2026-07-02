######### rotina para limpar o código #########
set +e

#apaga parte transiente
rm -rfv transient/data/*.dat 
rm -rfv transient/*.png

#apaga os restarts transientes
rm -rfv data/results/PTZH/*.dat 
rm -rfv data/results/V/*.dat
rm -rfv data/results/U/*.dat

#rm -rfv data/restart/*.dat

#apaga graficos gerados pelo pos_graphics
rm -rfv output/*.png

#apaga dados gerados pelo pos_process.f90, preservando output_variables.dat
find data -maxdepth 1 -type f -name '*.dat' ! -name 'output_variables.dat' -print -delete
rm -rfv data/error.dat
rm -rfv data/flametip.dat
rm -rfv data/probe.dat
rm -rfv data/probe_point.dat
rm -rfv data/*.vtk
shopt -s nullglob
for f in *.log *.gif; do
    rm -vf "$f"
done
shopt -u nullglob
