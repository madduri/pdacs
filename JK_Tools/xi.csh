#PBS -q debug
#PBS -A hacc
#PBS -l nodes=2:ppn=8
#PBS -l walltime=00:30:00 
#PBS -N xi
#PBS -o xi.log
#PBS -e xi.err
#PBS -j oe  

cd $PBS_O_WORKDIR  

#Example arguments: 
#mpirun ... xi.out ngrid inputfile nfiles outputfile boxsize [line of sight w omega_m]
#last 3 args are needed for redshift space power spectra only

mpirun -n 16 /global/u1/j/jkwan/powerspectra/PDACS/xi.out 128 /project/projectdirs/hacc/Coyote/M000/G004/Gadget_1.0000 32 /project/projectdirs/cosmosim/anl/jkwan/testing/xi_M000_pdacs_test 936  



