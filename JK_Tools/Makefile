powerspec_coyote.out: powerspec_coyote.o
	cc powerspec_coyote.o -L/opt/fftw/2.1.5.2/lib -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw -o powerspec_coyote.out

powerspec_coyote.o: powerspec_coyote.c
	cc powerspec_coyote.c -c -g -I/opt/fftw/2.1.5.2/include

powerspec_haloes.out: powerspec_haloes.o
	cc powerspec_haloes.o -L/opt/fftw/2.1.5.2/lib -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw -L/usr/common/usg/gsl/1.14/lib -lgsl -lgslcblas -lm -o powerspec_haloes.out

powerspec_haloes.o: powerspec_haloes.c
	cc powerspec_haloes.c -c -g -I/opt/fftw/2.1.5.2/include -I/usr/common/usg/gsl/1.14/include

powerspec_haloes_combined.out: powerspec_haloes_combined.o
	cc powerspec_haloes_combined.o -L/opt/fftw/2.1.5.2/lib -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw -L/usr/common/usg/gsl/1.14/lib -lgsl -lgslcblas -lm -o powerspec_haloes_combined.out

powerspec_haloes_combined.o: powerspec_haloes_combined.c
	cc powerspec_haloes_combined.c -c -g -I/opt/fftw/2.1.5.2/include -I/usr/common/usg/gsl/1.14/include

powerspec_haloes_binned.out: powerspec_haloes_binned.o
	cc powerspec_haloes_binned.o -L/opt/fftw/2.1.5.2/lib -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw -L/usr/common/usg/gsl/1.14/lib -lgsl -lgslcblas -lm -o powerspec_haloes_binned.out

powerspec_haloes_binned.o: powerspec_haloes_binned.c
	cc powerspec_haloes_binned.c -c -g -I/opt/fftw/2.1.5.2/include -I/usr/common/usg/gsl/1.14/include

powerspec_haloes_averaged.out: powerspec_haloes_averaged.o
	cc powerspec_haloes_averaged.o -L/opt/fftw/2.1.5.2/lib -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw -L/usr/common/usg/gsl/1.14/lib -lgsl -lgslcblas -lm -o powerspec_haloes_averaged.out

powerspec_haloes_averaged.o: powerspec_haloes_averaged.c
	cc powerspec_haloes_averaged.c -c -g -I/opt/fftw/2.1.5.2/include -I/usr/common/usg/gsl/1.14/include


powerspec.out: powerspec.o
	mpicc powerspec.o ${FFTW} -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw -L/usr/common/usg/gsl/1.14/lib -lgsl -lgslcblas -lm -o powerspec.out

powerspec.o: powerspec.c
	mpicc powerspec.c -c -g ${FFTW}

power_cross_auto.out: power_cross_auto.o
	cc power_cross_auto.o -L/opt/fftw/2.1.5.2/lib -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw -o power_cross_auto.out

power_cross_auto.o: power_cross_auto.c
	cc power_cross_auto.c -c -g -I/opt/fftw/2.1.5.2/include

powerspec_cross_haloes.out: powerspec_cross_haloes.o
	cc powerspec_cross_haloes.o -L/opt/fftw/2.1.5.2/lib -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw -L/usr/common/usg/gsl/1.14/lib -lgsl -lgslcblas -lm -o powerspec_cross_haloes.out

powerspec_cross_haloes.o: powerspec_cross_haloes.c
	cc powerspec_cross_haloes.c -c -g -I/opt/fftw/2.1.5.2/include -I/usr/common/usg/gsl/1.14/include


xi.out: xi.o
	mpicc xi.o ${FFTW} -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw -lm -o xi.out

xi.o: xi.c
	mpicc xi.c -c -g ${FFTW}


subsample.out: subsample.o
	cc subsample.o -L/usr/common/usg/gsl/1.14/lib -lgsl -lgslcblas -lm -o subsample.out

subsample.o: subsample.c
	cc subsample.c -c -g -I/usr/common/usg/gsl/1.14/include

hmf.out: hmf.o
	cc hmf.o -lm -o hmf.out

hmf.o: hmf.c
	cc hmf.c -c -g 

residuals.out: residuals.o
	cc residuals.o -L/usr/common/usg/gsl/1.14/lib -lgsl -lgslcblas  -lm -o residuals.out

residuals.o: residuals.c
	cc residuals.c -c -g -I/usr/common/usg/gsl/1.14/include


powerspec_avg.out: powerspec_avg.o
	cc powerspec_avg.o -lm -o powerspec_avg.out

powerspec_avg.o: powerspec_avg.c
	cc powerspec_avg.c -c -g 
 

butterworth.out: butterworth.o
	cc butterworth.o -L/opt/fftw/2.1.5.2/lib -lrfftw -lfftw -L/usr/common/usg/gsl/1.14/lib -lgsl -lgslcblas  -lm -o butterworth.out

butterworth.o: butterworth.c
	cc butterworth.c -c -g -I/opt/fftw/2.1.5.2/include -I/usr/common/usg/gsl/1.14/include

fft_test.out: fft_test.o
	cc fft_test.o -L/opt/fftw/2.1.5.2/lib -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw -o fft_test.out

fft_test.o: fft_test.c
	cc fft_test.c -c -g -I/opt/fftw/2.1.5.2/include