#PBS -q debug
#PBS -A hacc
#PBS -l nodes=2:ppn=8
#PBS -l walltime=00:30:00 
#PBS -N pow
#PBS -o pow.log
#PBS -e pow.err
#PBS -j oe  

cd $PBS_O_WORKDIR  

module load fftw
module load gsl 

#echo $PWD

cd /project/projectdirs/hacc/PDACS/JK_Tools/
#Example arguments: 
#mpirun ... powerspec.out ngrid inputfile nfiles outputfile boxsize [line of sight w omega_m]
#last 3 args are needed for redshift space power spectra only

mpirun -n 16 ./powerspec.out 128 /project/projectdirs/hacc/PDACS/Coyote/Grid/M001/L180/G001/snapshot_011 8 m001.pk.dat 107.586 



