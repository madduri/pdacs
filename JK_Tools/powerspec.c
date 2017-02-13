#include <stdio.h>
#include <stdlib.h>
#include <rfftw_mpi.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <sys/time.h>


#define SKIP(fp) fread(&dummy, sizeof(dummy), 1, (fp));
#define index(i,j,k,ngrid)  ((k) + (ngrid) * ((j) + (ngrid) * (i)))
#define fftw_index_before(i,j,k,ngrid)  ((k) + (2*((ngrid)/2+1)) * ((j) + (ngrid) * (i)))
#define fftw_index_after(i,j,k,ngrid)  ((k) + ((ngrid)/2+1)* ((j) + (ngrid) * (i)))
#define TSC     


struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
};

struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
/*   float Mass; */
  int    Id;
} ;

struct random_particles 
{
  double  Pos[3];
  double Radius;
} ;



struct extracted_haloes
{
  /* npart,nvpart,X,Y,Z,Vx,Vy,Vz,Mvir,Rvir,sigV */

  int   Any;
  int Mass;
  float  Pos[3];
  float  Vel[3];
  float    Mvir;
  float Rvir;
  float sigV;

} ;

int Ngas, TotPartInHalo, ntot_withmasses;
int NumHaloes,pc =0,pc_new,pc_sph; 

int ReadByteSwappedFloat(FILE *fptr,float *n); 
int ReadByteSwappedInt(FILE *fptr, int *n); 
void swaparrays(void** p1, void** p2); 
int impose_window_function(struct particle_data **P, int *NumPart, float BoxSize, int los); 
int load_lightcone(struct particle_data **P, struct io_header_1 *header1, int *NumPart);
int load_random_cat(struct particle_data **P, int *NumPart, int *nlocal, int random_no);
int load_harlequin(struct particle_data **P, int *NumPart, int *nlocal, int random_no); 
int load_horizon(struct particle_data **P, int *NumPart, int horizon_no);
int load_catalog(struct particle_data **P, int *NumPart, int *TotNumPart, char input_fname[]);
int load_density_field(char input_fname[], int ngrid, fftw_real *dens_in); 
int allocate_memory(int n, struct particle_data **P);
int get_random_dist(int m, struct random_particles **Prand, float survey_depth, float *survey_size);
int grid_dist_ngp(int n, int NumPart, int TotNumPart, float BoxSize, struct particle_data *P, fftw_real *dens_in);
int grid_dist_cic_mpi(int rank, int nproc, int n, int NumPart, int TotNumPart, float BoxSize, float mass, struct particle_data *P, 
		      int local_nx, int local_x_start,fftw_real *dens_in); 
int grid_dist_cic_mpi_2(int rank, int nproc, int n, int NumPart, int TotNumPart, float BoxSize, float mass, struct particle_data *P, 
			int local_nx, int local_x_start,fftw_real *dens_in); 
int grid_dist_tsc(int n, int NumPart, int TotNumPart, struct io_header_1 *header1, struct particle_data *P, fftw_real *dens_in);
int grid_dist_tsc_mpi(int rank, int nproc, int n, int NumPart, int TotNumPart, float BoxSize, float mass, struct particle_data *P, 
		      int local_nx, int local_x_start,fftw_real *dens_in);
int grid_dist_tsc_mpi_2(int rank, int nproc, int n, int NumPart, int TotNumPart, float BoxSize, float mass, struct particle_data *P, 
			int local_nx, int local_x_start,fftw_real *dens_in);
int grid_dist_D12(int rank, int nproc, int n, int NumPart, int TotNumPart, float BoxSize, float mass, struct particle_data *P, 
		  int local_nx, int local_x_start,fftw_real *dens_in); 
int swap_densities(int rank, int nproc, int ngrid, int local_nx, fftw_real *dens_m, fftw_real *dens_p, fftw_real *dens_in);
int get_pk(int ngrid, int NumPart, int local_ny, int local_y_start, float BoxSize, fftw_complex *dens_out, double *dpower,int *ipower);
int get_pk_new(int me, int Ng, int Np, int local_y_start_after_transpose, int local_ny_after_transpose, float BoxSize, fftw_complex *dgrid, double *pkall, char filename[]); 
int get_pk_binned(int ngrid, int NumPart, int local_ny, int local_y_start, float BoxSize, fftw_complex *dens_out, double *dpower, int *ipower);
int get_pk_redshift(int ngrid, int Numpart, int nmu, float *theta, int local_ny, int local_y_start, float BoxSize,
		    fftw_complex *dens_out, int los,double *dpower, int *ipower);
int get_pk_redshift_binned(int ngrid, int NumPart, int nmu, float *theta, int local_ny, int local_y_start, float BoxSize, 
		    fftw_complex *dens_out, int los, double *dpower, int *ipower); 
int get_pk_new_redshift(int rank, int nproc, int n, int nmu, int NumPart, int local_y_start_after_transpose, int local_ny_after_transpose, 
			float BoxSize, int los, fftw_complex*dgrid, char *output_fname1, char *output_fname2); 
int get_pk_para_perp(int ngrid, int NumPart, int nmu, float *theta, int local_ny, int local_y_start, float BoxSize, 
		     fftw_complex *dens_out, int los, double *dpower, int *ipower);
int redshift_space(int NumPart, float BoxSize, int los, float redshift, float w, float Omega0, struct particle_data *P);
int redshift_space_with_errors(int NumPart, float BoxSize, int los, float redshift, float w, float Omega0, struct particle_data *P);
int integrate_comoving(double *Omega_Lambda, int npts, double chi[], double redshift[]); 
int func (double t, const double y[], double f[], void *params); 
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);
int load_snapshot(int nfiles, int fileno, char input_fname[], int *NumPart, int *TotNumPart, struct particle_data **P, struct io_header_1 *header1); 
int load_hacc(int nfiles, int fileno, char input_fname[], int *NumPart, double *TotNumPart, float BoxSize, struct particle_data **P, int byteswapped); 
int my_modulo(int ngrid, int index);
float wrap_around(float hx, float x0, float x1, float BoxSize); 
int angular_projection(int k, int ks, int nmu, float *theta);
int equal_angle(int k, int i, int j, int ks, int nbin); 
int reallocate_memory(int n, struct particle_data **P);
int subdivide_lightcone(int NumPart,  float BoxSize, struct particle_data **P); 
int write_particles(char filename[], int NumPart, struct particle_data *P); 
int recenter_particles(int NumPart, struct particle_data *P, float *BoxSize, float *NewBoxSize); 
int randomise_coords(int *NumPart, struct particle_data *P, float BoxSize, int rank); 
unsigned long int get_random_seed(int verbose);
int apply_fkp_weighting(int NumPart, float BoxSize, int ngrid, float *W, fftw_real *dens); 
int get_weight_function(int ngrid, float BoxSize, float **W); 
double get_norm(float x, float y, float z); 
int setup_grid(int n, int local_nx, fftw_real *dens_in); 
int finish_grid(int n, int local_nx, int TotNumPart, float BoxSize, fftw_real *dens_in); 
int logtransform(int n, int local_nx, int TotNumPart, float BoxSize, fftw_real *dens_in);
int get_sigma_8_gal(int halfngrid, int nmu, float BoxSize, int *ipower, double *dpower, double *sig8);
int jackknife(int rank, int n, int local_nx, int local_x_start, float BoxSize, fftw_real *dens_in); 

int main(int argc, char **argv)
{
    int numprocs, rank, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int local_nx, local_x_start, local_ny_after_transpose,
	local_y_start_after_transpose, total_local_size;

    int ngrid = (int)atoi(argv[1]); /*set size of grid*/
    char input_fname[200], output_fname[200];
    sprintf(input_fname,"%s",argv[2]); 
    int files = (int)atoi(argv[3]);
    sprintf(output_fname,"%s",argv[4]); 
    float BoxSize = (float)atof(argv[5]); 

    int NumPart=0;
    int TotNumPart=0;
    FILE *fp, *fout;

    struct particle_data *Psurvey;
    struct io_header_1 header1;
    double *dpower, *dpower_local, *spower, *spower_local, *Pk, FFTnorm, dk; 
    int *ipower, *ipower_local;
    int i,j,k; 

    fftw_complex *dens_out;
    fftw_real *dens_in, *work;
    rfftwnd_mpi_plan p;



    MPI_Status *status;


    int mass = 1.;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &namelen);

    /* Setup arrays for MPI version of FFTW */

    p = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,
				ngrid, ngrid, ngrid,
				FFTW_REAL_TO_COMPLEX,
				FFTW_ESTIMATE);

    /* This is to ensure appropriate chunks of the data are parceled out to the processors */
    /*Remember: FFTW works with slabs of data so there is only a local version of
      x and no other coordinate*/

    rfftwnd_mpi_local_sizes(p, &local_nx, &local_x_start,
			    &local_ny_after_transpose,
			    &local_y_start_after_transpose,
			    &total_local_size);

    if (!(dens_in = (fftw_real*) malloc(sizeof(fftw_real) * total_local_size)))
    {
	fprintf(stderr, "Unable to allocate memory to density grid \n"); 
	exit(1); 
    }


    if (argc < 4)
    {
	fprintf(stderr, "Incorrect number of arguments entered: no of files, randnum, output filename, line of sight(optional)\n");
	return (0);
    }


    int byteswapped = 0; 
/*     float BoxSize = 1000; */


    int los;

/*     sprintf(input_fname, "/project/projectdirs/hacc/Coyote/M%03d/G004/Gadget_%.4f", model_no, scalefactor); */



    for (i = 0; i < files; i++)
    {
	Psurvey = NULL; 
/* 	load_density_field(input_fname, ngrid, dens_in); */
/* 	load_catalog(&Psurvey, &NumPart, &TotNumPart, input_fname); */
	load_snapshot(files, i, input_fname, &NumPart, &TotNumPart, &Psurvey, &header1);
/* 	load_hacc(files, fileno, input_fname, &NumPart, &TotNumPart, BoxSize, &Psurvey, byteswapped); */


	if (argc > 6)
	  {


	    los = (int)atoi(argv[6]);
	    float redshift = header1.redshift;
	    float w0 = (float)atoi(argv[7]);
	    float Omega0 = (float)atoi(argv[8]);
	    	    redshift_space(NumPart, BoxSize, los, redshift, w0, Omega0, Psurvey); /* convert to redshift space */
	    impose_window_function(&Psurvey, &NumPart, BoxSize, los);

	  }

	
/* 		grid_dist_tsc_mpi_2(rank, numprocs, ngrid, NumPart, TotNumPart, BoxSize, mass, Psurvey, local_nx, local_x_start, dens_in); */
	grid_dist_cic_mpi_2(rank, numprocs, ngrid, NumPart, TotNumPart, BoxSize, mass, Psurvey, local_nx, local_x_start, dens_in);
/* 	grid_dist_cic_mpi_2(rank, numprocs, ngrid, NumPart, TotNumPart, BoxSize, mass, Psub, local_nx, local_x_start, dens_in); */

/* 	if (i == files -1) /\*last file, finish up gridding*\/ */
/* 	  finish_grid(ngrid, local_nx, TotNumPart, BoxSize, dens_in); */
 
	free(Psurvey);

    }

/*     fprintf(stderr, "Finished CIC density deposition\n");  */


    /*  workspace is the same size as the data: */

    if (!(work = (fftw_real*) malloc(sizeof(fftw_real) * total_local_size)))
      {
	fprintf(stderr, "Unable to allocate memory to density grid \n"); 
	exit(1); 
      }

    rfftwnd_mpi(p, 1, dens_in, work, FFTW_TRANSPOSED_ORDER);

    /* the data is now complex, so typecast a pointer: */
    dens_out = (fftw_complex*) dens_in;

    /* Need to calculate squared modulus and then average dens_out in spherical shells in k-space */

    fftw_free(work); 
    rfftwnd_mpi_destroy_plan(p);

    /*     fprintf(stderr, "Finished FFTW. Now squaring and averaging over k-space. Also normalising. \n"); */
    FFTnorm         = (double)(ngrid*ngrid*ngrid);
    dk = 2*M_PI/BoxSize;   

    if (argc > 6)
      {
	int nmu = 21;
/* 	sprintf(output_fname, "/project/projectdirs/cosmosim/anl/heitmann/Coyote_xxs/powerspectra/pk_dd_M%03d_%d_cic_%.4f_z%d.dat",model_no,ngrid,header1.time,los); */
	char output_fname1[256], output_fname2[256]; 
	sprintf(output_fname1, "%s.multipole",output_fname);
	sprintf(output_fname2, "%s.2D",output_fname);
	get_pk_new_redshift(rank, numprocs, ngrid, nmu, TotNumPart, local_y_start_after_transpose, local_ny_after_transpose, BoxSize, los, dens_out, output_fname1, output_fname2); 

      }
    else
      {
	dpower_local = calloc(ngrid,sizeof(double)); dpower = calloc(ngrid,sizeof(double)); 
	ipower_local = calloc(ngrid,sizeof(int)); ipower = calloc(ngrid,sizeof(int));
	Pk = calloc(2*ngrid,sizeof(double));

	get_pk_new(rank, ngrid, TotNumPart, local_y_start_after_transpose, local_ny_after_transpose, BoxSize, dens_out, Pk, output_fname);


      }





  

  
    /*     fprintf(stderr, "Finished \n"); */
    fftw_free(dens_in);
/*     free(dpower); */
/*     free(ipower); */
/*     free(dpower_local);       */
/*     free(ipower_local); */
/*     free(Pk); */
    /*     free(W);  */
    MPI_Finalize();
    exit(0);
}



int get_weight_function(int ngrid, float BoxSize, float **W)
{
  int i, j, k; 
  float norm=0; 
  float nbar = 3e-4;     /*Mean number density of galaxies*/


  int index_cen = rint(0.5*BoxSize*ngrid/BoxSize+0.5); /*This is the coordinate of the center*/

  *W = malloc(ngrid*ngrid*ngrid*sizeof(float)); /*W is the weight function W(i,j,k) is the number of galaxies expected without clustering*/

  for (i =0; i < ngrid; i++)
    for (j =0; j < ngrid; j++)
      for (k =0; k < ngrid; k++)
	{

	  /*Work out the "radius" at each grid point */
	  int radius = floor(sqrt((i-index_cen) * (i-index_cen) + (j-index_cen)*(j-index_cen) + (k-index_cen)*(k-index_cen))); 

	  /* Distance between successive radii is BoxSize/ngrid */
	  float radius_inner = (float)radius - 0.5;/* *BoxSize/(float)ngrid;  */
	  float radius_outer = (float)radius + 0.5;/* *BoxSize/(float)ngrid;  */

	  float volume = 4*M_PI/3 * (pow(BoxSize/ngrid * radius_outer, 3)-pow(BoxSize/ngrid * radius_inner,3)) ;  

	  /*The expected number density in each cell is nbar*/

	  if (radius > index_cen)
	    (*W)[index(i,j,k,ngrid)] = 0;
	  else
	    (*W)[index(i,j,k,ngrid)] = volume * nbar; 
	  norm += (*W)[index(i,j,k,ngrid)]; 
	  /* 		fprintf(stderr, "%f %f %f %d %f\n", i,j,k, radius, W[index(i,j,k,ngrid)]);  */

	}

  FILE*fp = fopen("check_w.dat", "w"); 

  for (i =0; i < ngrid; i++)
    for (j =0; j < ngrid; j++)
      for (k =0; k < ngrid; k++)
	{
	  /* 		(*W)[index(i,j,k,ngrid)]/=norm;  */
	  fprintf(fp, "%d %d %d %f\n", i, j, k, (*W)[index(i,j,k,ngrid)]); 
	}

  fclose(fp); 

  return(0);
}

int apply_fkp_weighting(int NumPart, float BoxSize, int ngrid, float *W, fftw_real *dens)
{
  int i, j, k; 
  float ngrid3 = ngrid*ngrid*ngrid; 
  float nbar = 3e-4;     /*Mean number density of galaxies*/
  float P0 = 5000;
    


  for (i =0; i < ngrid; i++)
    for (j =0; j < ngrid; j++)
      for (k =0; k < ngrid; k++)
	dens[fftw_index_before(i,j,k,ngrid)] /= (1+ W[index(i,j,k,ngrid)]* ngrid3 * nbar * P0);


  return(0); 
}

int impose_window_function(struct particle_data **P, int *NumPart, float BoxSize, int los)
{
  int i, j, NewCount;
     

  /*     struct particle_data *P_dum = calloc(*NumPart, sizeof(struct particle_data)); */

  /*     for (i = 0, NewCount = 0; i < *NumPart ; i++) */
  /*     { */

  /* 	float cen = 0.5 * BoxSize; */

  /* 	float limits[7] = {0., 293.79, 575.124, 843.494, 1098.73, 1340.9, 1570.29}; */

  /* /\* 	float r_inner = 1098.73; *\/ */
  /* /\* 	float r_outer = 1340.9; *\/ */
  /* 	float r_inner = limits[5]; */
  /* 	float r_outer = limits[6]; */

  /* 	float xdist = (*P)[i].Pos[0]-cen; */
  /* 	float ydist = (*P)[i].Pos[1]-cen; */
  /* 	float zdist = (*P)[i].Pos[2]-cen; */

  /* /\* 	float radius = 0.5 * BoxSize; *\/ */
  /* /\* 	float check_eqn = xdist*xdist + ydist*ydist + zdist*zdist - radius*radius; *\/ */

  /* 	float check_eqn_inner = xdist*xdist + ydist*ydist + zdist*zdist-r_inner*r_inner; */
  /* 	float check_eqn_outer = xdist*xdist + ydist*ydist + zdist*zdist-r_outer*r_outer; */

  /* 	if (check_eqn_inner > 0 && check_eqn_outer < 0) */
  /* 	{ */
  /* 	    P_dum[NewCount].Pos[0] = (*P)[i].Pos[0]; */
  /* 	    P_dum[NewCount].Pos[1] = (*P)[i].Pos[1]; */
  /* 	    P_dum[NewCount].Pos[2] = (*P)[i].Pos[2]; */
  /* 	    P_dum[NewCount].Vel[0] = (*P)[i].Vel[0]; */
  /* 	    P_dum[NewCount].Vel[1] = (*P)[i].Vel[1]; */
  /* 	    P_dum[NewCount].Vel[2] = (*P)[i].Vel[2]; */
  /* 	    NewCount++; */
  /* 	} */
  /*     } */


  for (i = 0, NewCount = 0; i < *NumPart; i++)
    {
      if ((*P)[i].Pos[los] > 1.) 
	(*P)[i].Pos[los] -= 1.;
      else if((*P)[i].Pos[los] < 0)
	(*P)[i].Pos[los] += 1.;
    }

  /*     swaparrays(P, &P_dum);  */

  /*     free(P_dum);      */

  /*     *NumPart = NewCount; */


  return(0);
}

void swaparrays(void** p1, void** p2) 
{
  void* p3;
  p3 = *p1;
  *p1 = *p2;
  *p2 = p3;
}


/* int impose_mass_cutoff(int *NumPart, float mass_cutoff, struct particle_data **P) */
/* { */
/*     struct particle_data *P_dum;  */
/*     int i, NewNumPart;  */
/*     float err = 0.01*mass_cutoff; /\* 1% error *\/ */
    
/*     for (i = 0, NewNumPart = 0; i < *NumPart; i++) */
/* 	if (0.5*fabs((*P)[i].mass - mass_cutoff )< err) */
/* 	{ */
/* 	    P_dum[NewNumPart].Pos[0] = (*P)[i].Pos[0]; */
/* 	    P_dum[NewNumPart].Pos[1] = (*P)[i].Pos[1]; */
/* 	    P_dum[NewNumPart].Pos[2] = (*P)[i].Pos[2]; */
/* 	    NewNumPart++; */
/* 	} */

/*     swaparrays(P, &P_dum);  */

/*     free(P_dum);      */

/*     *NumPart = NewNumPart; */


/*     return(0);  */
/* } */

int subdivide_lightcone(int NumPart,  float BoxSize, struct particle_data **P)
{
  int i,j,k; 

  struct particle_data *Psub; 
  if (!(Psub = calloc(NumPart, sizeof(struct particle_data))))
    {
      fprintf(stderr, "Could not allocate memory to Psub\n"); 
      exit(1);
    }

  float HalfBox = 0.5*BoxSize; 

  for(i = 0; i < NumPart; i++)
    {

      if ((*P)[i].Pos[0] >= HalfBox && (*P)[i].Pos[1] >= HalfBox && (*P)[i].Pos[2] >= HalfBox) 
	{
	  for (j = 0; j < 3; j++)
	    (Psub)[i].Pos[j] = (*P)[i].Pos[j]-0.5*BoxSize; 

	}
      else if ((*P)[i].Pos[0] <= HalfBox && (*P)[i].Pos[1] >= HalfBox && (*P)[i].Pos[2] >= HalfBox) 
	{		
	  (Psub)[i].Pos[0] = (*P)[i].Pos[0]+0.5*BoxSize; 
	  (Psub)[i].Pos[1] = (*P)[i].Pos[1]-0.5*BoxSize; 
	  (Psub)[i].Pos[2] = (*P)[i].Pos[2]-0.5*BoxSize; 

	}

      else if ((*P)[i].Pos[0] <= HalfBox && (*P)[i].Pos[1] <= HalfBox && (*P)[i].Pos[2] >= HalfBox) 
	{
	  (Psub)[i].Pos[0] = (*P)[i].Pos[0]+0.5*BoxSize; 
	  (Psub)[i].Pos[1] = (*P)[i].Pos[1]+0.5*BoxSize; 
	  (Psub)[i].Pos[2] = (*P)[i].Pos[2]-0.5*BoxSize; 

	}

      else if ((*P)[i].Pos[0] >= HalfBox && (*P)[i].Pos[1] <= HalfBox && (*P)[i].Pos[2] >= HalfBox)
	{
	  (Psub)[i].Pos[0] = (*P)[i].Pos[0]-0.5*BoxSize; 
	  (Psub)[i].Pos[1] = (*P)[i].Pos[1]+0.5*BoxSize; 
	  (Psub)[i].Pos[2] = (*P)[i].Pos[2]-0.5*BoxSize; 

	}

      else if ((*P)[i].Pos[0] >= HalfBox && (*P)[i].Pos[1] >= HalfBox && (*P)[i].Pos[2] <= HalfBox) 
	{
	  (Psub)[i].Pos[0] = (*P)[i].Pos[0]-0.5*BoxSize; 
	  (Psub)[i].Pos[1] = (*P)[i].Pos[1]-0.5*BoxSize; 
	  (Psub)[i].Pos[2] = (*P)[i].Pos[2]+0.5*BoxSize; 

	}

	
      else if ((*P)[i].Pos[0] <= HalfBox && (*P)[i].Pos[1] >= HalfBox && (*P)[i].Pos[2] <= HalfBox) 
	{
	  (Psub)[i].Pos[0] = (*P)[i].Pos[0]+0.5*BoxSize; 
	  (Psub)[i].Pos[1] = (*P)[i].Pos[1]-0.5*BoxSize; 
	  (Psub)[i].Pos[2] = (*P)[i].Pos[2]+0.5*BoxSize; 

	}

      else if ((*P)[i].Pos[0] <= HalfBox && (*P)[i].Pos[1] <= HalfBox && (*P)[i].Pos[2] <= HalfBox) 
	{
	  for (j = 0; j < 3; j++)
	    (Psub)[i].Pos[j] = (*P)[i].Pos[j]+0.5*BoxSize; 

	}
      else if ((*P)[i].Pos[0] >= HalfBox && (*P)[i].Pos[1] <= HalfBox && (*P)[i].Pos[2] <= HalfBox) 
	{
	  (Psub)[i].Pos[0] = (*P)[i].Pos[0]-0.5*BoxSize; 
	  (Psub)[i].Pos[1] = (*P)[i].Pos[1]+0.5*BoxSize; 
	  (Psub)[i].Pos[2] = (*P)[i].Pos[2]+0.5*BoxSize; 

	}

      else 
	{
	  fprintf(stderr, "Incorrect case i.e. rank selected \n"); 
	  fprintf(stderr, "%f %f %f \n", (*P)[i].Pos[0], (*P)[i].Pos[1], (*P)[i].Pos[2]); 	  
	  exit(1);
	}
    }

  swaparrays(P, &Psub); 

  free(Psub);     


  return(0);
}


int load_random_cat(struct particle_data **P, int *NumPart, int *nlocal, int random_no)
{
  FILE *fp;
  char   fname[200];
  int    i,j,k,NumHaloes;
  struct particle_data *Pdummy;

  /* Open file containing lightcone */

  /*     sprintf(fname,"/nfs/cluster/cosmic/jkwan/HorizonRun/windowfunctions/window.gaussian.2.dat"); */
  sprintf(fname,"/home/ssi/jkwan/testwinbox_0.dat");

  if(!(fp = fopen(fname,"rb")))
    {
      fprintf(stderr, "Unable to open file %s \n", fname); 
      exit(1);
    }


  fread(NumPart,sizeof(int),1,fp);
  if(!random_no)
    {
      allocate_memory(*NumPart, P);
      *nlocal = 0;
    }

  int n =*nlocal, nreads = 1;

  for(n=*nlocal;nreads!=0;n++)
    nreads = fread(&((*P)[n].Pos[0]),sizeof(float),3, fp);


  *nlocal = n-1;


  fclose(fp);

  return (0);
}


int load_harlequin(struct particle_data **P, int *NumPart, int *nlocal, int random_no)
{
  FILE *fp;
  char   fname[200];
  int    i,j,k,NumHaloes;
  struct particle_data *Pdummy;

  /* Open file containing lightcone */

  /*     sprintf(fname,"/nfs/cluster/cosmic/jkwan/HorizonRun/windowfunctions/window.gaussian.2.dat"); */
  sprintf(fname,"/nfs/cluster/cosmic/jkwan/HorizonRun/extractedhaloes_%d", random_no);

  if(!(fp = fopen(fname,"r")))
    {
      fprintf(stderr, "Unable to open file %s \n", fname);
      exit(1);
    }


  /*     fread(NumPart,sizeof(int),1,fp); */
  /*     if(!random_no) */
  /*     { */
  /* 	allocate_memory(*NumPart, P); */
  /* 	*nlocal = 0; */
  /*     } */
  /*     else */
  /* 	reallocate_memory(*NumPart+*nlocal, P); */

  int n =*nlocal, nreads = 1;
  float any; 

  /*     for(n=*nlocal;nreads!=0;n++) */
  /* 	nreads = fread(&((*P)[n].Pos[0]),sizeof(float),3, fp); */

  /*     for(n=*nlocal;nreads!=0;n++) */
  /*       { */
  /* 	nreads = fread(&((*P)[n].Pos[0]),sizeof(float),3, fp); */
  /* 	fread(&((*P)[n].Vel[0]),sizeof(float),3, fp); */

  /* 	if ((*P)[n].Pos[0] < 25) */
  /* 	  (*P)[n].Pos[0] += (*P)[n].Vel[0]*0.0909091/100.; */
  /* 	else */
  /* 	  (*P)[n].Pos[0] += (*P)[n].Vel[0]/100.; */
  /*       } */
  /*     *nlocal = n-1; */



  fscanf(fp,"%d",nlocal);

  fprintf(stderr, "%d %d\n", *nlocal, *NumPart);  
  if(!random_no)
    allocate_memory(*nlocal, P);
  else
    reallocate_memory(*NumPart+*nlocal, P);


  for(n=*NumPart;n< (*NumPart+*nlocal) ;n++)
    {
      fscanf(fp,"%f %f %f %f %f %f %f ",&((*P)[n].Pos[0]),&((*P)[n].Pos[1]) , 
	     &((*P)[n].Pos[2]), &((*P)[n].Vel[0]),&((*P)[n].Vel[1]) , 
	     &((*P)[n].Vel[2]), &any); 
  
    }

  *NumPart += *nlocal;



  fclose(fp);

  return (0);
}


int load_horizon(struct particle_data **P, int *NumPart, int horizon_no)
{
  int num_categories = 6;
  FILE *fp[num_categories]; 
  char filenames[num_categories][256];
  int num_haloes_in_file,i, j;
  struct info *haloes_dummy;
  float h; 


  sprintf(filenames[0], "/nfs/cluster/cosmic/jkwan/HorizonRun/mock_sdss_%d/x_sdss3_0%d.dat", horizon_no, horizon_no);
  sprintf(filenames[1], "/nfs/cluster/cosmic/jkwan/HorizonRun/mock_sdss_%d/y_sdss3_0%d.dat", horizon_no, horizon_no);
  sprintf(filenames[2], "/nfs/cluster/cosmic/jkwan/HorizonRun/mock_sdss_%d/z_sdss3_0%d.dat", horizon_no, horizon_no);
  sprintf(filenames[3], "/nfs/cluster/cosmic/jkwan/HorizonRun/mock_sdss_%d/vx_sdss3_0%d.dat", horizon_no, horizon_no);
  sprintf(filenames[4], "/nfs/cluster/cosmic/jkwan/HorizonRun/mock_sdss_%d/vy_sdss3_0%d.dat", horizon_no, horizon_no);
  sprintf(filenames[5], "/nfs/cluster/cosmic/jkwan/HorizonRun/mock_sdss_%d/vz_sdss3_0%d.dat", horizon_no, horizon_no);
  /*   sprintf(filenames[6], "/nfs/cluster/cosmic/jkwan/HorizonRun/mock_sdss_%d/mass_sdss3_0%d.dat", horizon_no, horizon_no); */

  /*   sprintf(filenames[0], "/nfs/cluster/cosmic/jkwan/HorizonRun/x_cubic_neal00401.dat");  */
  /*   sprintf(filenames[1], "/nfs/cluster/cosmic/jkwan/HorizonRun/y_cubic_neal00401.dat");  */
  /*   sprintf(filenames[2], "/nfs/cluster/cosmic/jkwan/HorizonRun/z_cubic_neal00401.dat");  */
  /*   sprintf(filenames[3], "/nfs/cluster/cosmic/jkwan/HorizonRun/vx_cubic_neal00401.dat");  */
  /*   sprintf(filenames[4], "/nfs/cluster/cosmic/jkwan/HorizonRun/vy_cubic_neal00401.dat");  */
  /*   sprintf(filenames[5], "/nfs/cluster/cosmic/jkwan/HorizonRun/vz_cubic_neal00401.dat");  */
  /*   sprintf(filenames[6], "/nfs/cluster/cosmic/jkwan/HorizonRun/mock_sdss_%d/mass_sdss3_0%d.dat", horizon_no, horizon_no); */

  float **x;



  for (i = 0; i < num_categories ; i++)
    {

      if (!(fp[i] = fopen(filenames[i], "rb")))
	{
	  fprintf(stderr, "could not open file %s \n", filenames[i]);
	  exit(1);
	}


      fread(NumPart, sizeof(int), 1, fp[i]);

      if (i == 0)
	x = malloc(num_categories*sizeof(float * ));

      *(x+i) = malloc((*NumPart)*sizeof(float)); 

      fread(*(x+i), sizeof(float), *NumPart, fp[i]); 


    }

  *P = calloc(*NumPart, sizeof(struct particle_data));

  float offset = 1570.;
  /*   float offset = 0.; */

  for(j=0; j < *NumPart; j++)
    {
      (*P)[j].Pos[0] = *(*x+j) + offset;
      (*P)[j].Pos[1] = *(*(x+1)+j) + offset;
      (*P)[j].Pos[2] = *(*(x+2)+j) + offset;
      (*P)[j].Vel[0] = *(*(x+3)+j);
      (*P)[j].Vel[1] = *(*(x+4)+j);
      (*P)[j].Vel[2] = *(*(x+5)+j);
      /*         (*P)[j].Mass = *(*(x+6)+j); */
    }
  


  for (i = 0; i < num_categories; i++)
    free(*(x+i));


  free(x);

  return(0);

}

int load_density_field(char input_fname[], int ngrid, fftw_real *dens_in)

{
  FILE *fp; 
  fp = fopen(input_fname,"r");
  int count=0; 
/*   int nread = fscanf(fp, "%lf", &(dens_in[count])); */

  for (count = 0; count < ngrid*ngrid*ngrid; count++)
      fscanf(fp, "%lf", &(dens_in[count]));
		     
 /*  while (nread != EOF) */
/*     { */
/*       count++;  */
/*       nread = fscanf(fp, "%lf", &(dens_in[count])); */

/*     } */

  return(0); 
}

int load_catalog(struct particle_data **P, int *NumPart, int *TotNumPart, char input_fname[])
{
  FILE *fp;
  char   fname[200];
  int    i,j,k,n,NumHaloes, chk_eof;
  float any;
  struct particle_data Pdummy;

  /* Open file containing lightcone */

 
  fp = fopen(input_fname,"r");

  /*     fscanf(fp,"%f %f %f %f %f %f %f",&any, &any, &any, &any, &any,&any,&any);  */

/*   fscanf(fp,"%d",NumPart); */

  /*     fprintf(stderr, "%d\n",*NumPart); */

/*   *TotNumPart = *NumPart; */
  int maxpart = 1e8;
  (*P) = (struct particle_data*)malloc(maxpart*sizeof(struct particle_data)); 
/*   allocate_memory(*NumPart, P); */

  n = *TotNumPart; 
  int nread = fscanf(fp,"%f %f %f",&((*P)[n].Pos[0]),&((*P)[n].Pos[1]) , 
		 &((*P)[n].Pos[2]));  

  while(nread!=EOF)
    {
      n++; 
      nread = fscanf(fp,"%f %f %f",&((*P)[n].Pos[0]),&((*P)[n].Pos[1]) , 
	     &((*P)[n].Pos[2]));  


      /*       (*P)[n].Pos[0] /= 1000.; */
      /*       (*P)[n].Pos[1] /= 1000.; */
      /*       (*P)[n].Pos[2] /= 1000.; */
    }
  *NumPart = n;
  (*TotNumPart) += n; 

  fclose(fp);

  return (0);
}



int load_lightcone(struct particle_data **P, struct io_header_1 *header1, int *NumPart)
{
  FILE *fp;
  char   fname[200];
  int    i,j,k,dummy;
  int    t,n,off,pc_sph;
  struct particle_data *Pdummy;

  /* Open file containing lightcone */

  sprintf(fname,"/nfs/cluster/cosmic/jkwan/lightcone_0.248/SDSS_blake08.dat");
  fp = fopen(fname,"r");
  
  fread(&dummy, sizeof(dummy), 1, fp);
  fread(header1, sizeof(struct io_header_1), 1, fp);
  fread(&dummy, sizeof(dummy), 1, fp);


  {
    for(k=0, *NumPart=0, ntot_withmasses=0; k<5; k++)
      *NumPart+= (*header1).npart[k];
    Ngas= (*header1).npart[0];
  }

  for(k=0, ntot_withmasses=0; k<5; k++)
    {
      if((*header1).mass[k]==0)
	ntot_withmasses+= (*header1).npart[k];
    }

  allocate_memory(*NumPart,P);
  allocate_memory(*NumPart,&Pdummy);

  SKIP(fp);
  for(k=0,pc_new=pc;k<6;k++)
    {
      for(n=0;n<(*header1).npart[k];n++)
	{
	  fread(&Pdummy[pc_new].Pos[0], sizeof(float), 3, fp);
	  Pdummy[pc_new].Pos[0] /= 1000.;       /*Convert to Mpc/h */
	  Pdummy[pc_new].Pos[1] /= 1000.;
	  Pdummy[pc_new].Pos[2] /= 1000.;
	  pc_new++;
	}
    }
  SKIP(fp);

    

  SKIP(fp);
  for(k=0,pc_new=pc;k<6;k++)
    {
      for(n=0;n<(*header1).npart[k];n++)
	{
	  fread(&Pdummy[pc_new].Vel[0], sizeof(float), 3, fp);
	  pc_new++;
	}
    }
  SKIP(fp);

  /*   SKIP(fp); */
  /*   for(k=0,pc_new=pc;k<6;k++) */
  /*     { */
  /*       for(n=0;n<(*header1).npart[k];n++) */
  /* 	{ */
  /* 	  fread(&Pdummy[pc_new].Radius, sizeof(float), 1, fp); */
  /* 	  pc_new++; */
  /* 	} */
  /*     } */
  /*   SKIP(fp); */

  /*  SKIP;
      for(k=0,pc_new=pc;k<6;k++)
      {
      for(n=0;n<(*header1).npart[k];n++)
      {
      fread(&Pdummy[pc_new].Mass, sizeof(float), 1, fp);
      pc_new++;
      }
      }
      SKIP;*/
    

  /*   SKIP(fp); */

  /*   for(k=0,pc_new=pc;k<6;k++) */
  /*     { */
  /*       for(n=0;n<(*header1).npart[k];n++) */
  /* 	{ */
  /* 	  fread(&Pdummy[pc_new].Id, sizeof(int), 1, fp); */
  /* 	  pc_new++; */
  /* 	} */
  /*     } */
  /*   SKIP(fp); */

  for (k = 0; k < *NumPart; k++)
    (*P)[k] = Pdummy[k];

  free(Pdummy);
  fclose(fp);
}


/* this routine allocates the memory for the 
 * particle data.
 */


int allocate_memory(int n, struct particle_data **P)  /*Need to pass &P as args */
{
  printf("allocating memory...\n");

  if(!(*P=malloc(n*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }

  return (0);
}



int grid_dist_ngp(int n, int NumPart, int TotNumPart, float BoxSize, struct particle_data *P, fftw_real *dens_in)
{
  int i, j, k,l, index, flag = 0;
  float rhobar;
  double xnorm,ynorm,znorm, deltax, deltay, deltaz;

  /*    BoxSize = 256.;*/

  /* Need to clear contents of dens_in */

  for (i=0;i<n;i++) 
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
	dens_in[index(i,j,k,n)] = 0.0;     

  /* Need to place galaxies and random sample on to a grid */ 

  deltax = BoxSize/n;
  deltay = BoxSize/n;
  deltaz = BoxSize/n;

  for (l=0; l < NumPart; l++)
    {
      /* coordinates in subcube units [0,NGRID] */
      xnorm = P[l].Pos[0]/deltax;
      ynorm = P[l].Pos[1]/deltay;
      znorm = P[l].Pos[2]/deltaz;

      /* index of nearest grid point [0,NGRID] */
      i  = rint(xnorm);
      j  = rint(ynorm);
      k  = rint(znorm);

      /* distance to nearest grid point (TSC ONLY) */
      /*hx  = rrx - (double)ix;
	hy  = rry - (double)iy;
	hz  = rrz - (double)iz;*/
      
      /* keep track of peridoc boundaries -> [0,NGRID-1] ; NGRID -> 0   */
      i = my_modulo(n,i);
      j = my_modulo(n,j);
      k = my_modulo(n,k);
      /*ix=(int)fmod(ix,NGRID);
	iy=(int)fmod(iy,NGRID);
	iz=(int)fmod(iz,NGRID);*/

      dens_in[index(i,j,k,n)] += 1;

    }


  /*   rhobar = TotNumPart/(BoxSize*BoxSize*BoxSize); */
  rhobar = 1;  

  /* Calculate the density contrast */

  for (i=0;i<n;i++) 
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
	dens_in[index(i,j,k,n)] = dens_in[index(i,j,k,n)]/rhobar-1;  

  FILE*fp=fopen("check_ngp.dat", "w"); 

  for (i=0;i<n;i++) 
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
	fprintf(fp, "%d %d %d %f\n", i, j, k, dens_in[index(i,j,k,n)]); 

  fclose(fp); 
  return(0);

}
int jackknife(int rank, int n, int local_nx, int local_x_start, float BoxSize, fftw_real *dens_in)
{
  int i,j,k; 

  int ncubes = 3; 

  float fftw_grid_length = BoxSize/n;    // the length of one fftw cube

  float jackknife_grid_length  = BoxSize/ncubes; // the length of one jackknife cube

  int gridpts_per_cube = (int)(jackknife_grid_length/fftw_grid_length+0.5); 

  int cube_coords[3] = {1, 2, 0}; 

  int xmin = cube_coords[0]*gridpts_per_cube; int xmax = gridpts_per_cube*(cube_coords[0]+1); 
  int ymin = cube_coords[1]*gridpts_per_cube; int ymax = gridpts_per_cube*(cube_coords[1]+1); 
  int zmin = cube_coords[2]*gridpts_per_cube; int zmax = gridpts_per_cube*(cube_coords[2]+1); 

  /*   fprintf(stderr, "fftw grid length = %f, jackknife_grid_length = %f, gridpts_per_cube = %d \n",  */
  /* 	  fftw_grid_length, jackknife_grid_length, gridpts_per_cube);  */

  FILE *fptest; 
  char filetest[250];
  /*   sprintf(filetest, "jackknife_%d.dat", rank); */
  
  /*   fptest=fopen(filetest,"w");  */


  for (i=0;i<local_nx;i++) 
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
	if(i+local_x_start >= xmin && i+local_x_start < xmax && j >= ymin && j < ymax && k >= zmin && k < zmax)
	  {
	    dens_in[fftw_index_before(i,j,k,n)] = 0.0;     
	    /* 	    fprintf(fptest, "%d %d %d \n", i+local_x_start, j, k);  */
	  }

  /*   if (local_nx < gridpts_per_cube && local_x_start < (gridpts_per_cube) */
  /*     { */
  /*       for (i=0;i<local_nx;i++)  */
  /* 	for (j = 0; j < gridpts_per_cube; j++) */
  /* 	  for (k = 0; k < gridpts_per_cube; k++) */
  /* 	    { */
  /* 	    dens_in[fftw_index_before(i,j,k,n)] = 0.0;      */
  /* 	    fprintf(fptest, "%d %d %d \n", i+local_x_start, j, k);  */
  /* 	    } */
  /*     } */
  /*   else */
  /*     { */
  /*       for (i=0;i<(gridpts_per_cube-local_x_start);i++)  */
  /* 	for (j = 0; j < gridpts_per_cube; j++) */
  /* 	  for (k = 0; k < gridpts_per_cube; k++) */
  /* 	    { */
  /* 	    dens_in[fftw_index_before(i,j,k,n)] = 0.0;      */
  /* 	    fprintf(fptest, "%d %d %d \n", i+local_x_start, j, k);  */
  /* 	    } */
 /*       fprintf(stderr, "HERE!! Rank: %d\n", rank);  */
  /*     } */





  /*     fclose(fptest); */


  return(0); 
}

int grid_dist_cic_mpi_2(int rank, int nproc, int n, int NumPart, int TotNumPart, float BoxSize, float mass, struct particle_data *P, 
		      int local_nx, int local_x_start,fftw_real *dens_in)
{
    int i, j, k,l, index, flag = 0;
    float rhobar;
    double xnorm,ynorm,znorm, deltax, deltay, deltaz, dx, dy, dz;
    double hx, hxp, hxm, hy, hyp, hym, hz, hzp, hzm;
    int ip, im, jp, jm, kp, km;
    FILE *fp;
    int check_in_box; 


    double ng = 0.9999999*n; 
    int Ng = n; 
    int Ng2 = n+2; 
    int count=0,overload[2]; 
    int ii[3];

    overload[0] = 0;
    overload[1] = 0;

    for (l=0; l < NumPart; l++)
    {

	i  = ng*P[l].Pos[0]; 
/* 	dx = ng*P[l].Pos[0]-i; */

/* 	if (i==(local_x_start-1) && fabs(ng*P[l].Pos[0]-local_x_start)<0.001) */
/* 	{ */
/* 	    i=local_x_start; */
/* /\* 	    dx = ng*P[l].Pos[0]-i; *\/ */
/* 	} */


	
    
	if (i>=(local_x_start-1) && i<(local_nx+local_x_start))
	{
	    if (i < (local_nx+local_x_start-1) && i > local_x_start-1)
	    {
		i     -= local_x_start;
		ii[0] = (i+1);
		dx    = ng*P[l].Pos[0]-(i+local_x_start);
		j  = ng*P[l].Pos[1]; 
		ii[1] = (j+1)%Ng;
		dy = ng*P[l].Pos[1]-j; 
		k  = ng*P[l].Pos[2]; 
		ii[2] = (k+1)%Ng;
		dz = ng*P[l].Pos[2]-k;

/* 		ix   -= xlim[0]; */
/* 		ii[0] = (ix+1); */
/* 		dx    = ng*P[nn].Pos[0]-(ix+xlim[0]); */
/* 		iy    = ng*P[nn].Pos[1]; */
/* 		ii[1] = (iy+1)%Ng; */
/* 		dy    = ng*P[nn].Pos[1]-iy; */
/* 		iz    = ng*P[nn].Pos[2]; */
/* 		ii[2] = (iz+1)%Ng; */
/* 		dz    = ng*P[nn].Pos[2]-iz; */

		dens_in[Ng*Ng2* i    + Ng2* j    + k   ]   += (1.0-dx)*(1.0-dy)*(1.0-dz);
		dens_in[Ng*Ng2*ii[0] + Ng2* j    + k   ]   +=      dx *(1.0-dy)*(1.0-dz);
		dens_in[Ng*Ng2* i    + Ng2*ii[1] + k   ]   += (1.0-dx)*     dy *(1.0-dz);
		dens_in[Ng*Ng2* i    + Ng2* j    +ii[2]]   += (1.0-dx)*(1.0-dy)*     dz ;
		dens_in[Ng*Ng2*ii[0] + Ng2*ii[1] + k   ]   +=      dx *     dy *(1.0-dz);
		dens_in[Ng*Ng2*ii[0] + Ng2* j    +ii[2]]   +=      dx *(1.0-dy)*     dz ;
		dens_in[Ng*Ng2* i    + Ng2*ii[1] +ii[2]]   += (1.0-dx)*     dy *     dz ;
		dens_in[Ng*Ng2*ii[0] + Ng2*ii[1] +ii[2]]   +=      dx *     dy *     dz ;
		
	    }
	    else if (i == (local_nx+local_x_start-1))
	    {
		i = local_nx-1;
		ii[0] = (i+1);
		dx    = ng*P[l].Pos[0]-(i+local_x_start);
		j  = ng*P[l].Pos[1]; 
		ii[1] = (j+1)%Ng;
		dy = ng*P[l].Pos[1]-j; 
		k  = ng*P[l].Pos[2]; 
		ii[2] = (k+1)%Ng;
		dz = ng*P[l].Pos[2]-k;

		dens_in[Ng*Ng2* i    + Ng2* j    + k   ]   += (1.0-dx)*(1.0-dy)*(1.0-dz);
		dens_in[Ng*Ng2* i    + Ng2*ii[1] + k   ]   += (1.0-dx)*     dy *(1.0-dz);
		dens_in[Ng*Ng2* i    + Ng2* j    +ii[2]]   += (1.0-dx)*(1.0-dy)*     dz ;
		dens_in[Ng*Ng2* i    + Ng2*ii[1] +ii[2]]   += (1.0-dx)*     dy *     dz ;

		overload[0]++;
/* 		count++;  */

	    }
	    else if (i == (local_x_start-1) && nproc > 1)
	    {
		i    -= local_x_start;
/* 	    ii[0] = local_x_start; // so ii[0] is now local_nx */
		ii[0] = 0;
		dx    = ng*P[l].Pos[0]-(i+local_x_start);
/* 		dx = ng*(P[l].Pos[0]/(double)BoxSize)-(double)(i+local_x_start); */
		j  = ng*P[l].Pos[1];
		ii[1] = (j+1)%Ng;
		dy = ng*P[l].Pos[1]-j;
		k  = ng*P[l].Pos[2];
		ii[2] = (k+1)%Ng;
		dz = ng*P[l].Pos[2]-k;

		overload[1]++;


		dens_in[Ng*Ng2*ii[0] + Ng2* j    + k   ]   +=      dx *(1.0-dy)*(1.0-dz);
		dens_in[Ng*Ng2*ii[0] + Ng2*ii[1] + k   ]   +=      dx *     dy *(1.0-dz);
		dens_in[Ng*Ng2*ii[0] + Ng2* j    +ii[2]]   +=      dx *(1.0-dy)*     dz ;
		dens_in[Ng*Ng2*ii[0] + Ng2*ii[1] +ii[2]]   +=      dx *     dy *     dz ;


	    }

	}
	if (i == n-1 && rank == 0)
	{
	    i    -= local_x_start;
	    dx    = ng*P[l].Pos[0]-(i+local_x_start);
	    ii[0] = 0; 
	    j  = ng*P[l].Pos[1];
	    ii[1] = (j+1)%Ng;
	    dy = ng*P[l].Pos[1]-j;
	    k  = ng*P[l].Pos[2];
	    ii[2] = (k+1)%Ng;
	    dz = ng*P[l].Pos[2]-k;

	    dens_in[Ng*Ng2*ii[0] + Ng2* j    + k   ]   +=      dx *(1.0-dy)*(1.0-dz);
	    dens_in[Ng*Ng2*ii[0] + Ng2*ii[1] + k   ]   +=      dx *     dy *(1.0-dz);
	    dens_in[Ng*Ng2*ii[0] + Ng2* j    +ii[2]]   +=      dx *(1.0-dy)*     dz ;
	    dens_in[Ng*Ng2*ii[0] + Ng2*ii[1] +ii[2]]   +=      dx *     dy *     dz ;


	    if (nproc==1)
	    {
		i = n-1;
		dx    = ng*P[l].Pos[0]-(i+local_x_start);
		dens_in[Ng*Ng2* i    + Ng2* j    + k   ]   += (1.0-dx)*(1.0-dy)*(1.0-dz);
		dens_in[Ng*Ng2* i    + Ng2*ii[1] + k   ]   += (1.0-dx)*     dy *(1.0-dz);
		dens_in[Ng*Ng2* i    + Ng2* j    +ii[2]]   += (1.0-dx)*(1.0-dy)*     dz ;
		dens_in[Ng*Ng2* i    + Ng2*ii[1] +ii[2]]   += (1.0-dx)*     dy *     dz ;

	    }

	    overload[1]++;


	}
    }

    return(0); 
}


int grid_dist_cic_mpi(int rank, int nproc, int n, int NumPart, int TotNumPart, float BoxSize, float mass, struct particle_data *P, 
			int local_nx, int local_x_start,fftw_real *dens_in)
{
  int i, j, k,l, index, flag = 0;
  float rhobar;
  double xnorm,ynorm,znorm, deltax, deltay, deltaz, dx, dy, dz;
  double hx, hxp, hxm, hy, hyp, hym, hz, hzp, hzm;
  int ip, im, jp, jm, kp, km;
  FILE *fp;
  int check_in_box; 


  int n2 = 2*(n/2+1);
  fftw_real *dens_p; 
  fftw_real *dens_m; 




  /* Need to place galaxies and random sample on to a grid */ 

  deltax = BoxSize/n;
  deltay = BoxSize/n;
  deltaz = BoxSize/n;



  for (l=0; l < NumPart; l++)
    {
      /* coordinates in subcube units [0,NGRID] */
      xnorm = P[l].Pos[0]/deltax;
      ynorm = P[l].Pos[1]/deltay;
      znorm = P[l].Pos[2]/deltaz;


      /* index of nearest grid point [0,NGRID] */
      i  = (int)(xnorm+0.5);
      j  = (int)(ynorm+0.5);
      k  = (int)(znorm+0.5);

      /* 	/\*Distance of particle from cell center*\/ */
/*       dx=xnorm-floor(xnorm)-0.5; */
/*       dy=ynorm-floor(ynorm)-0.5; */
/*       dz=znorm-floor(znorm)-0.5; */

      dx = xnorm-(double)i;
      dy = ynorm-(double)j;
      dz = znorm-(double)k;

/*       im = floor(xnorm);  */
/*       ip = im+1;  */

/*       jm = floor(ynorm);  */
/*       jp = jm+1;  */

/*       km = floor(znorm);  */
/*       kp = km+1;  */

      ip = i; hxp = 1-fabs(dx); hxm = fabs(dx); 
      jp = j; hyp = 1-fabs(dy); hym = fabs(dy); 
      kp = k; hzp = 1-fabs(dz); hzm = fabs(dz); 


/*       jm = j-1;  */
/*       jp = j;  */

/*       km = k-1; */
/*       kp = k;  */


      if (dx > 0)
	{
/* 	  hxp = 1-dx;  */
/* 	  hxm = dx; */
	  im = i+1;
 
	}
      else
	{
/* 	  hxp = fabs(dx);  */
/* 	  hxm = 1-fabs(dx);  */
	  im = i-1;

	}

      if (dy > 0)
	{
/* 	  hyp = 1-dy;  */
/* 	  hym = dy;  */
	  jm = j+1;
	}
      else
	{
/* 	  hyp = fabs(dy);  */
/* 	  hym = 1-fabs(dy);  */
	  jm = j-1; 
	}

      if (dz > 0)
	{
/* 	  hzp = 1-dz;  */
/* 	  hzm = dz;  */
	  km = k+1;
	}
      else
	{
/* 	  hzp = fabs(dz);  */
/* 	  hzm = 1-fabs(dz);  */
	  km=k-1;
	}



      /* FFTW is split along x direction */
      /* local_nx is the width of the slice */
      /* local_x_start is the starting xnorm of local slab*/

      if(ip >= n) ip=ip-n;  if(im >= n) im=im-n;            /*Wrap around for period boundary conditions*/
      if(jp >= n) jp=jp-n;  if(jm >= n) jm=jm-n;
      if(kp >= n) kp=kp-n;  if(km >= n) km=km-n;     

      if(ip < 0) ip=ip+n;   if(im < 0) im=im+n;             
      if(jp < 0) jp=jp+n;   if(jm < 0) jm=jm+n;             
      if(kp < 0) kp=kp+n;   if(km < 0) km=km+n;             

      /*Determine a new i coordinate i.e. i needs to be mapped from [0,NGRID] to [-1,local_nx]*/
      check_in_box = im-local_x_start;  /*check_in_box checks that particle is inside slab or in adjacent slab*/

      im = im - local_x_start; 
      ip = ip - local_x_start; 



      if (ip >= -1 && ip <= local_nx && im >= -1 && im <= local_nx)  
	{
     
	  /* Indicies of neighbouring grid point, p(lus) (ngp+1) and m(inus) (ngp-1), taking care of
	     periodic boundaries*/


	  if (im == -1 && ip == 0)
	    {
	      /*Particle inside, neightbour outside, lower edge  i.e. ip = 0*/
	      dens_in[fftw_index_before(ip,jp,kp,n)]  += hxp*hyp*hzp*mass;
	      dens_in[fftw_index_before(ip,jp,km,n)]  += hxp*hyp*hzm*mass;
	      dens_in[fftw_index_before(ip,jm,kp,n)]  += hxp*hym*hzp*mass;
	      dens_in[fftw_index_before(ip,jm,km,n)]  += hxp*hym*hzm*mass;

	    }
	  else if (ip == -1 && im == 0)
	    {
	      /*Particle outside, neighbour inside, lower edge  i.e. ip = 0*/
	      dens_in[fftw_index_before(im,jp,kp,n)]  += hxm*hyp*hzp*mass;
	      dens_in[fftw_index_before(im,jp,km,n)]  += hxm*hyp*hzm*mass;
	      dens_in[fftw_index_before(im,jm,kp,n)]  += hxm*hym*hzp*mass;
	      dens_in[fftw_index_before(im,jm,km,n)]  += hxm*hym*hzm*mass;

	    }

	  else if (im == local_nx && ip == (local_nx-1) )
	    {
	      /*Particle inside, neigbour outside, upper edge i.e. ip = local_nx-1*/
	      dens_in[fftw_index_before(ip,jp,kp,n)]  += hxp*hyp*hzp*mass;
	      dens_in[fftw_index_before(ip,jp,km,n)]  += hxp*hyp*hzm*mass;
	      dens_in[fftw_index_before(ip,jm,kp,n)]  += hxp*hym*hzp*mass;
	      dens_in[fftw_index_before(ip,jm,km,n)]  += hxp*hym*hzm*mass;

	    }
	  else if (im == (local_nx-1) && ip == local_nx)
	    {
	      /*Particle outside, neighbour inside, upper edge i.e. im = local_nx-1*/
	      dens_in[fftw_index_before(im,jp,kp,n)]  += hxm*hyp*hzp*mass;
	      dens_in[fftw_index_before(im,jp,km,n)]  += hxm*hyp*hzm*mass;
	      dens_in[fftw_index_before(im,jm,kp,n)]  += hxm*hym*hzp*mass;
	      dens_in[fftw_index_before(im,jm,km,n)]  += hxm*hym*hzm*mass;

	    }
	  else 
	    {
	      /*Bin everything*/
	      /* 		if (im < 0)  */
	      /* 		    fprintf(stderr, "last bin: im < 0, %f \n", xnorm);  */
	      /* 		if (ip >= local_nx) */
	      /* 		    fprintf(stderr, "last bin: ip >= local_nx, %f \n", xnorm);  */

	      dens_in[fftw_index_before(ip,jp,kp,n)]  += hxp*hyp*hzp*mass;
	      dens_in[fftw_index_before(ip,jp,km,n)]  += hxp*hyp*hzm*mass;
	      dens_in[fftw_index_before(ip,jm,kp,n)]  += hxp*hym*hzp*mass;
	      dens_in[fftw_index_before(ip,jm,km,n)]  += hxp*hym*hzm*mass;

	      dens_in[fftw_index_before(im,jp,kp,n)]  += hxm*hyp*hzp*mass;
	      dens_in[fftw_index_before(im,jp,km,n)]  += hxm*hyp*hzm*mass;
	      dens_in[fftw_index_before(im,jm,kp,n)]  += hxm*hym*hzp*mass;
	      dens_in[fftw_index_before(im,jm,km,n)]  += hxm*hym*hzm*mass;

	    }


	}
/*       else if(rank == 0 && im== n && nproc > 1) */
/* 	{ */
/* 	  /\*wrap around to top of box*\/ */
/* 	  dens_in[fftw_index_before(0 ,jp,kp,n)]  += hxp*hyp*hzp*mass; */
/* 	  dens_in[fftw_index_before(0 ,jp,km,n)]  += hxp*hyp*hzm*mass; */
/* 	  dens_in[fftw_index_before(0 ,jm,kp,n)]  += hxp*hym*hzp*mass; */
/* 	  dens_in[fftw_index_before(0 ,jm,km,n)]  += hxp*hym*hzm*mass; */
/* 	} */
/*       else if(rank == nproc-1 && im==0 && nproc > 1) */
/* 	{ */

/* 	  /\*wrap around to bottom of box*\/ */
/* 	  dens_in[fftw_index_before(local_nx-1,jp,kp,n)]  += hxm*hyp*hzp*mass; */
/* 	  dens_in[fftw_index_before(local_nx-1,jp,km,n)]  += hxm*hyp*hzm*mass; */
/* 	  dens_in[fftw_index_before(local_nx-1,jm,kp,n)]  += hxm*hym*hzp*mass; */
/* 	  dens_in[fftw_index_before(local_nx-1,jm,km,n)]  += hxm*hym*hzm*mass; */
/* 	} */


	    
	
    }


  return(0);

}


int grid_dist_tsc(int n, int NumPart, int TotNumPart, struct io_header_1 *header1, struct particle_data *P, fftw_real *dens_in)
{
  int i, j, k,l, index, flag = 0;
  float BoxSize, rhobar;
  double xnorm,ynorm,znorm, deltax, deltay, deltaz, dx, dy, dz;
  double hx, hxp, hxm, hy, hyp, hym, hz, hzp, hzm;
  int ip, im, jp, jm, kp, km;
  FILE *fp;
  double mass;


  for (i=0;i<n;i++) 
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
	dens_in[index(i,j,k,n)] = 0.0;     

  BoxSize = (*header1).BoxSize;
  mass = (*header1).mass[1];

  /*    BoxSize = 256.;*/

  /* Need to clear contents of dens_in */


  /* Need to place galaxies and random sample on to a grid */ 

  deltax = BoxSize/n;
  deltay = BoxSize/n;
  deltaz = BoxSize/n;


  for (l=0; l < NumPart; l++)
    {
      /* coordinates in subcube units [0,NGRID] */
      xnorm = P[l].Pos[0]/deltax;
      ynorm = P[l].Pos[1]/deltay;
      znorm = P[l].Pos[2]/deltaz;


      /* index of nearest grid point [0,NGRID] */
      i  = floor(xnorm);
      j  = floor(ynorm);
      k  = floor(znorm);

      if(i >= n) i=i-n;             /*Wrap around for period boundary conditions*/
      if(j >= n) j=j-n;
      if(k >= n) k=k-n;


      /*Distance of particle from cell center*/
      /*       dx=xnorm-floor(xnorm)-0.5;               /\*ASSUMES SUBCUBES ARE OF UNIT LENGTH, IF NOT MULTIPLY 0.5 BY DELTAX*\/ */
      /*       dy=ynorm-floor(ynorm)-0.5; */
      /*       dz=znorm-floor(znorm)-0.5; */

      dx = xnorm-(double)i;
      dy = ynorm-(double)j;
      dz = znorm-(double)k;
     
      /* Indicies of neighbouring grid point, p(lus) (ngp+1) and m(inus) (ngp-1), taking care of
	 periodic boundaries*/
      ip = my_modulo(n,i+1); 
      im = my_modulo(n,i-1); 

      jp = my_modulo(n,j+1); 
      jm = my_modulo(n,j-1); 

      kp = my_modulo(n,k+1); 
      km = my_modulo(n,k-1); 
      
      /*Weighting for TSC scheme */
      hx  =0.75 - dx*dx;
      hxp =0.5  * (0.5+dx)*(0.5+dx) ;
      hxm =0.5  * (0.5-dx)*(0.5-dx); 
      hy  =0.75 - dy*dy;
      hyp =0.5  * (0.5+dy)*(0.5+dy) ;
      hym =0.5  * (0.5-dy)*(0.5-dy); 
      hz  =0.75 - dz*dz;
      hzp =0.5  * (0.5+dz)*(0.5+dz) ;
      hzm =0.5  * (0.5-dz)*(0.5-dz);

      /* Sweep over ix+1 first...*/
      dens_in[index(ip,jp,kp,n)]  += hxp*hyp*hzp*mass;
      dens_in[index(ip,jp,k ,n)]  += hxp*hyp*hz *mass;
      dens_in[index(ip,jp,km,n)]  += hxp*hyp*hzm*mass;
      dens_in[index(ip,j ,kp,n)]  += hxp*hy *hzp*mass;
      dens_in[index(ip,j ,k ,n)]  += hxp*hy *hz *mass;
      dens_in[index(ip,j ,km,n)]  += hxp*hy *hzm*mass;
      dens_in[index(ip,jm,kp,n)]  += hxp*hym*hzp*mass;
      dens_in[index(ip,jm,k ,n)]  += hxp*hym*hz *mass;
      dens_in[index(ip,jm,km,n)]  += hxp*hym*hzm*mass;

      /* Then over ix...*/
      dens_in[index(i ,jp,kp,n)]  += hx *hyp*hzp*mass;
      dens_in[index(i ,jp,k ,n)]  += hx *hyp*hz *mass;
      dens_in[index(i ,jp,km,n)]  += hx *hyp*hzm*mass;
      dens_in[index(i ,j ,kp,n)]  += hx *hy *hzp*mass;
      dens_in[index(i ,j ,k ,n)]  += hx *hy *hz *mass;
      dens_in[index(i ,j ,km,n)]  += hx *hy *hzm*mass;
      dens_in[index(i ,jm,kp,n)]  += hx *hym*hzp*mass;
      dens_in[index(i ,jm,k ,n)]  += hx *hym*hz *mass;
      dens_in[index(i ,jm,km,n)]  += hx *hym*hzm*mass;

      /* Finally over ix-1.*/
      dens_in[index(im,jp,kp,n)]  += hxm*hyp*hzp*mass;
      dens_in[index(im,jp,k ,n)]  += hxm*hyp*hz *mass;
      dens_in[index(im,jp,km,n)]  += hxm*hyp*hzm*mass;
      dens_in[index(im,j ,kp,n)]  += hxm*hy *hzp*mass;
      dens_in[index(im,j ,k ,n)]  += hxm*hy *hz *mass;
      dens_in[index(im,j ,km,n)]  += hxm*hy *hzm*mass;
      dens_in[index(im,jm,kp,n)]  += hxm*hym*hzp*mass;
      dens_in[index(im,jm,k ,n)]  += hxm*hym*hz *mass;
      dens_in[index(im,jm,km,n)]  += hxm*hym*hzm*mass;


    }


  rhobar = TotNumPart/(BoxSize*BoxSize*BoxSize);   

  

  /* Calculate the density contrast */

  for (i=0;i<n;i++) 
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
	dens_in[index(i,j,k,n)] = dens_in[index(i,j,k,n)]/rhobar-1;

  return(0);

}
int setup_grid(int n, int local_nx, fftw_real *dens_in)
{
  int i, j, k;     
  for (i=0;i<local_nx;i++) 
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
	dens_in[fftw_index_before(i,j,k,n)] = 0.0;     




  return(0);
}

int grid_dist_tsc_mpi(int rank, int nproc, int n, int NumPart, int TotNumPart, float BoxSize, float mass, struct particle_data *P, 
		      int local_nx, int local_x_start,fftw_real *dens_in)
{
  int i, j, k,l, index, flag = 0;
  float rhobar;
  double xnorm,ynorm,znorm, deltax, deltay, deltaz, dx, dy, dz;
  double hx, hxp, hxm, hy, hyp, hym, hz, hzp, hzm;
  int ip, im, jp, jm, kp, km;
  FILE *fp;
  int check_in_box; 


  int n2 = 2*(n/2+1);
  fftw_real *dens_p; 
  fftw_real *dens_m; 

  if(!(dens_p = calloc(n*n2, sizeof(fftw_real))))
    {
      fprintf(stderr, "Unable to allocate memory to dens_p \n"); 
      exit(1);
    }
  
  if(!(dens_m = calloc(n*n2, sizeof(fftw_real))))
    {
      fprintf(stderr, "Unable to allocate memory to dens_m \n"); 
      exit(1);
    }
  



  /* Need to place galaxies and random sample on to a grid */ 

  deltax = BoxSize/n;
  deltay = BoxSize/n;
  deltaz = BoxSize/n;



  for (l=0; l < NumPart; l++)
    {
      /* coordinates in subcube units [0,NGRID] */
      xnorm = P[l].Pos[0]/deltax;
      ynorm = P[l].Pos[1]/deltay;
      znorm = P[l].Pos[2]/deltaz;


      /* index of nearest grid point [0,NGRID] */
      i  = (int)(xnorm+0.5);
      j  = (int)(ynorm+0.5);
      k  = (int)(znorm+0.5);

      /* FFTW is split along x direction */
      /* local_nx is the width of the slice */
      /* local_x_start is the starting xnorm of local slab*/

      if(i >= n) i=i-n;             /*Wrap around for period boundary conditions*/
      if(j >= n) j=j-n;
      if(k >= n) k=k-n;

      if(i < 0) i=i+n;             
      if(j < 0) j=j+n;
      if(k < 0) k=k+n;

      check_in_box = i-local_x_start;  /*check_in_box checks that particle is inside slab*/

      if (check_in_box >= 0 && check_in_box < local_nx )  
	{

	  /*Determine a new i coordinate i.e. i needs to be mapped from [0,NGRID] to [0,local_nx)*/
	  i -=local_x_start;

	  /*Distance of particle from cell center*/
	  /* 	  dx=xnorm-floor(xnorm)-0.5;               /\*ASSUMES SUBCUBES ARE OF UNIT LENGTH, IF NOT MULTIPLY 0.5 BY DELTAX*\/ */
	  /* 	  dy=ynorm-floor(ynorm)-0.5; */
	  /* 	  dz=znorm-floor(znorm)-0.5; */


	  dx = xnorm-(double)i;
	  dy = ynorm-(double)j;
	  dz = znorm-(double)k;
     
	  /* Indicies of neighbouring grid point, p(lus) (ngp+1) and m(inus) (ngp-1), taking care of
	     periodic boundaries*/


	  /* 	  ip = my_modulo(n,i+1);  */
	  /* 	  im = my_modulo(n,i-1);  */

	  ip = i+1; 
	  im = i-1; 

	  jp = my_modulo(n,j+1); 
	  jm = my_modulo(n,j-1); 

	  kp = my_modulo(n,k+1); 
	  km = my_modulo(n,k-1); 
      
	  /*Weighting for TSC scheme */
	  hx  =0.75 - dx*dx;
	  hxp =0.5  * (0.5+dx)*(0.5+dx) ;
	  hxm =0.5  * (0.5-dx)*(0.5-dx); 
	  hy  =0.75 - dy*dy;
	  hyp =0.5  * (0.5+dy)*(0.5+dy) ;
	  hym =0.5  * (0.5-dy)*(0.5-dy); 
	  hz  =0.75 - dz*dz;
	  hzp =0.5  * (0.5+dz)*(0.5+dz) ;
	  hzm =0.5  * (0.5-dz)*(0.5-dz);

	  /* 	  if (check_in_box == (local_nx + 1 )) /\*Only i is affected by slicing of cube *\/ */
	  if (ip == local_nx)
	    {
	      /*These go the next processor with rankup ...*/
	      dens_p[kp+(n2*jp)] += hxp*hyp*hzp*mass;
	      dens_p[k +(n2*jp)] += hxp*hyp*hz *mass;
	      dens_p[km+(n2*jp)] += hxp*hyp*hzm*mass;
	      dens_p[kp+(n2* j)] += hxp*hy *hzp*mass;
	      dens_p[k +(n2* j)] += hxp*hy *hz *mass;
	      dens_p[km+(n2* j)] += hxp*hy *hzm*mass;
	      dens_p[kp+(n2*jm)] += hxp*hym*hzp*mass;
	      dens_p[k +(n2*jm)] += hxp*hym*hz *mass;
	      dens_p[km+(n2*jm)] += hxp*hym*hzm*mass;
	    }
	  else
	    {
	      /* Sweep over ix+1 first...*/
	      dens_in[fftw_index_before(ip,jp,kp,n)]  += hxp*hyp*hzp*mass;
	      dens_in[fftw_index_before(ip,jp,k ,n)]  += hxp*hyp*hz *mass;
	      dens_in[fftw_index_before(ip,jp,km,n)]  += hxp*hyp*hzm*mass;
	      dens_in[fftw_index_before(ip,j ,kp,n)]  += hxp*hy *hzp*mass;
	      dens_in[fftw_index_before(ip,j ,k ,n)]  += hxp*hy *hz *mass;
	      dens_in[fftw_index_before(ip,j ,km,n)]  += hxp*hy *hzm*mass;
	      dens_in[fftw_index_before(ip,jm,kp,n)]  += hxp*hym*hzp*mass;
	      dens_in[fftw_index_before(ip,jm,k ,n)]  += hxp*hym*hz *mass;
	      dens_in[fftw_index_before(ip,jm,km,n)]  += hxp*hym*hzm*mass;
	    }

	  /* Then over ix...*/
	  dens_in[fftw_index_before(i ,jp,kp,n)]  += hx *hyp*hzp*mass;
	  dens_in[fftw_index_before(i ,jp,k ,n)]  += hx *hyp*hz *mass;
	  dens_in[fftw_index_before(i ,jp,km,n)]  += hx *hyp*hzm*mass;
	  dens_in[fftw_index_before(i ,j ,kp,n)]  += hx *hy *hzp*mass;
	  dens_in[fftw_index_before(i ,j ,k ,n)]  += hx *hy *hz *mass;
	  dens_in[fftw_index_before(i ,j ,km,n)]  += hx *hy *hzm*mass;
	  dens_in[fftw_index_before(i ,jm,kp,n)]  += hx *hym*hzp*mass;
	  dens_in[fftw_index_before(i ,jm,k ,n)]  += hx *hym*hz *mass;
	  dens_in[fftw_index_before(i ,jm,km,n)]  += hx *hym*hzm*mass;

	  /* 	  if (check_in_box == 0)  */
	  if (im == -1)
	    {
	      /* These go to the next processor with rank-1 i.e rankdown */
	      dens_m[kp+(n2*jp)] += hxm*hyp*hzp*mass;
	      dens_m[k +(n2*jp)] += hxm*hyp*hz *mass;
	      dens_m[km+(n2*jp)] += hxm*hyp*hzm*mass;
	      dens_m[kp+(n2* j)] += hxm*hy *hzp*mass;
	      dens_m[k +(n2* j)] += hxm*hy *hz *mass;
	      dens_m[km+(n2* j)] += hxm*hy *hzm*mass;
	      dens_m[kp+(n2*jm)] += hxm*hym*hzp*mass;
	      dens_m[k +(n2*jm)] += hxm*hym*hz *mass;
	      dens_m[km+(n2*jm)] += hxm*hym*hzm*mass;
	    }

	  else
	    {
	      /* Finally over ix-1.*/
	      dens_in[fftw_index_before(im,jp,kp,n)]  += hxm*hyp*hzp*mass;
	      dens_in[fftw_index_before(im,jp,k ,n)]  += hxm*hyp*hz *mass;
	      dens_in[fftw_index_before(im,jp,km,n)]  += hxm*hyp*hzm*mass;
	      dens_in[fftw_index_before(im,j ,kp,n)]  += hxm*hy *hzp*mass;
	      dens_in[fftw_index_before(im,j ,k ,n)]  += hxm*hy *hz *mass;
	      dens_in[fftw_index_before(im,j ,km,n)]  += hxm*hy *hzm*mass;
	      dens_in[fftw_index_before(im,jm,kp,n)]  += hxm*hym*hzp*mass;
	      dens_in[fftw_index_before(im,jm,k ,n)]  += hxm*hym*hz *mass;
	      dens_in[fftw_index_before(im,jm,km,n)]  += hxm*hym*hzm*mass;
	    }
	}
    }


  swap_densities(rank, nproc, n, local_nx, dens_m, dens_p, dens_in);




  free(dens_m); free(dens_p); 

  return(0);

}

int grid_dist_tsc_mpi_2(int rank, int nproc, int n, int NumPart, int TotNumPart, float BoxSize, float mass, struct particle_data *P, 
			int local_nx, int local_x_start,fftw_real *dens_in)
{
  int i, j, k,l, index, flag = 0;
  float rhobar;
  double xnorm,ynorm,znorm, deltax, deltay, deltaz, dx, dy, dz;
  double hx, hxp, hxm, hy, hyp, hym, hz, hzp, hzm;
  int ip, im, jp, jm, kp, km;
  FILE *fp;
  int check_in_box; 


  int n2 = 2*(n/2+1);
  fftw_real *dens_p; 
  fftw_real *dens_m; 

  /*     if(!(dens_p = calloc(n*n2, sizeof(fftw_real)))) */
  /*     { */
  /* 	fprintf(stderr, "Unable to allocate memory to dens_p \n");  */
  /* 	exit(1); */
  /*     } */
  
  /*     if(!(dens_m = calloc(n*n2, sizeof(fftw_real)))) */
  /*     { */
  /* 	fprintf(stderr, "Unable to allocate memory to dens_m \n");  */
  /* 	exit(1); */
  /*     } */
  



  /* Need to place galaxies and random sample on to a grid */ 

  deltax = BoxSize/n;
  deltay = BoxSize/n;
  deltaz = BoxSize/n;



  for (l=0; l < NumPart; l++)
    {
      /* coordinates in subcube units [0,NGRID] */
      xnorm = P[l].Pos[0]/deltax;
      ynorm = P[l].Pos[1]/deltay;
      znorm = P[l].Pos[2]/deltaz;


      /* index of nearest grid point [0,NGRID] */
      i  = (int)(xnorm+0.5);
      j  = (int)(ynorm+0.5);
      k  = (int)(znorm+0.5);
      /* 	i  = floor(xnorm); */
      /* 	j  = floor(ynorm); */
      /* 	k  = floor(znorm); */


      /* 	/\*Distance of particle from cell center*\/ */
      /* 	dx=xnorm-floor(xnorm)-0.5;             */
      /* 	dy=ynorm-floor(ynorm)-0.5; */
      /* 	dz=znorm-floor(znorm)-0.5; */

      dx = xnorm-(double)i;
      dy = ynorm-(double)j;
      dz = znorm-(double)k;

      hx  =0.75 - dx*dx;
      hxp =0.5  * (0.5+dx)*(0.5+dx) ;
      hxm =0.5  * (0.5-dx)*(0.5-dx); 
      hy  =0.75 - dy*dy;
      hyp =0.5  * (0.5+dy)*(0.5+dy) ;
      hym =0.5  * (0.5-dy)*(0.5-dy); 
      hz  =0.75 - dz*dz;
      hzp =0.5  * (0.5+dz)*(0.5+dz) ;
      hzm =0.5  * (0.5-dz)*(0.5-dz);


      /* FFTW is split along x direction */
      /* local_nx is the width of the slice */
      /* local_x_start is the starting xnorm of local slab*/

      if(i >= n) i=i-n;             /*Wrap around for period boundary conditions*/
      if(j >= n) j=j-n;
      if(k >= n) k=k-n;

      if(i < 0) i=i+n;             
      if(j < 0) j=j+n;
      if(k < 0) k=k+n;


      check_in_box = i-local_x_start;  /*check_in_box checks that particle is inside slab or in adjacent slab*/

      if (check_in_box >= -1 && check_in_box < local_nx+1 )  
	{

	  /*Determine a new i coordinate i.e. i needs to be mapped from [0,NGRID] to [-1,local_nx]*/
	  /* 	    i -=local_x_start; */
	  i = check_in_box;

     
	  /* Indicies of neighbouring grid point, p(lus) (ngp+1) and m(inus) (ngp-1), taking care of
	     periodic boundaries*/


	  /* 	  ip = my_modulo(n,i+1);  */
	  /* 	  im = my_modulo(n,i-1);  */

	  ip = i+1; 
	  im = i-1; 


	  jp = my_modulo(n,j+1); 
	  jm = my_modulo(n,j-1); 

	  kp = my_modulo(n,k+1); 
	  km = my_modulo(n,k-1); 
      
	  /*Weighting for TSC scheme */

	  if (i == -1)
	    {
	      /*Only bin the ones whose interpolation overlap into this slab i.e. ip = 0*/
	      dens_in[fftw_index_before(ip,jp,kp,n)]  += hxp*hyp*hzp*mass;
	      dens_in[fftw_index_before(ip,jp,k ,n)]  += hxp*hyp*hz *mass;
	      dens_in[fftw_index_before(ip,jp,km,n)]  += hxp*hyp*hzm*mass;
	      dens_in[fftw_index_before(ip,j ,kp,n)]  += hxp*hy *hzp*mass;
	      dens_in[fftw_index_before(ip,j ,k ,n)]  += hxp*hy *hz *mass;
	      dens_in[fftw_index_before(ip,j ,km,n)]  += hxp*hy *hzm*mass;
	      dens_in[fftw_index_before(ip,jm,kp,n)]  += hxp*hym*hzp*mass;
	      dens_in[fftw_index_before(ip,jm,k ,n)]  += hxp*hym*hz *mass;
	      dens_in[fftw_index_before(ip,jm,km,n)]  += hxp*hym*hzm*mass;

	    }
	  else if (i == local_nx) 
	    {
	      /*Only bin the ones whose interpolation region overlap into this slab i.e. im = local_nx-1*/
	      dens_in[fftw_index_before(im,jp,kp,n)]  += hxm*hyp*hzp*mass;
	      dens_in[fftw_index_before(im,jp,k ,n)]  += hxm*hyp*hz *mass;
	      dens_in[fftw_index_before(im,jp,km,n)]  += hxm*hyp*hzm*mass;
	      dens_in[fftw_index_before(im,j ,kp,n)]  += hxm*hy *hzp*mass;
	      dens_in[fftw_index_before(im,j ,k ,n)]  += hxm*hy *hz *mass;
	      dens_in[fftw_index_before(im,j ,km,n)]  += hxm*hy *hzm*mass;
	      dens_in[fftw_index_before(im,jm,kp,n)]  += hxm*hym*hzp*mass;
	      dens_in[fftw_index_before(im,jm,k ,n)]  += hxm*hym*hz *mass;
	      dens_in[fftw_index_before(im,jm,km,n)]  += hxm*hym*hzm*mass;

	    }
	  else if (i == 0)

	    {

	      /* 		if (ip >= local_nx) */
	      /* 		    fprintf(stderr, "ip >= local_nx, %f\n", xnorm);  */


	      dens_in[fftw_index_before(ip,jp,kp,n)]  += hxp*hyp*hzp*mass;
	      dens_in[fftw_index_before(ip,jp,k ,n)]  += hxp*hyp*hz *mass;
	      dens_in[fftw_index_before(ip,jp,km,n)]  += hxp*hyp*hzm*mass;
	      dens_in[fftw_index_before(ip,j ,kp,n)]  += hxp*hy *hzp*mass;
	      dens_in[fftw_index_before(ip,j ,k ,n)]  += hxp*hy *hz *mass;
	      dens_in[fftw_index_before(ip,j ,km,n)]  += hxp*hy *hzm*mass;
	      dens_in[fftw_index_before(ip,jm,kp,n)]  += hxp*hym*hzp*mass;
	      dens_in[fftw_index_before(ip,jm,k ,n)]  += hxp*hym*hz *mass;
	      dens_in[fftw_index_before(ip,jm,km,n)]  += hxp*hym*hzm*mass;

	      dens_in[fftw_index_before(i ,jp,kp,n)]  += hx *hyp*hzp*mass;
	      dens_in[fftw_index_before(i ,jp,k ,n)]  += hx *hyp*hz *mass;
	      dens_in[fftw_index_before(i ,jp,km,n)]  += hx *hyp*hzm*mass;
	      dens_in[fftw_index_before(i ,j ,kp,n)]  += hx *hy *hzp*mass;
	      dens_in[fftw_index_before(i ,j ,k ,n)]  += hx *hy *hz *mass;
	      dens_in[fftw_index_before(i ,j ,km,n)]  += hx *hy *hzm*mass;
	      dens_in[fftw_index_before(i ,jm,kp,n)]  += hx *hym*hzp*mass;
	      dens_in[fftw_index_before(i ,jm,k ,n)]  += hx *hym*hz *mass;
	      dens_in[fftw_index_before(i ,jm,km,n)]  += hx *hym*hzm*mass;

	    }

	  else if (i == local_nx-1)
	    {

	      /* 		if (im < 0)  */
	      /* 		    fprintf(stderr, "im < 0, %f \n", xnorm);  */


	      dens_in[fftw_index_before(i ,jp,kp,n)]  += hx *hyp*hzp*mass;
	      dens_in[fftw_index_before(i ,jp,k ,n)]  += hx *hyp*hz *mass;
	      dens_in[fftw_index_before(i ,jp,km,n)]  += hx *hyp*hzm*mass;
	      dens_in[fftw_index_before(i ,j ,kp,n)]  += hx *hy *hzp*mass;
	      dens_in[fftw_index_before(i ,j ,k ,n)]  += hx *hy *hz *mass;
	      dens_in[fftw_index_before(i ,j ,km,n)]  += hx *hy *hzm*mass;
	      dens_in[fftw_index_before(i ,jm,kp,n)]  += hx *hym*hzp*mass;
	      dens_in[fftw_index_before(i ,jm,k ,n)]  += hx *hym*hz *mass;
	      dens_in[fftw_index_before(i ,jm,km,n)]  += hx *hym*hzm*mass;

	      dens_in[fftw_index_before(im,jp,kp,n)]  += hxm*hyp*hzp*mass;
	      dens_in[fftw_index_before(im,jp,k ,n)]  += hxm*hyp*hz *mass;
	      dens_in[fftw_index_before(im,jp,km,n)]  += hxm*hyp*hzm*mass;
	      dens_in[fftw_index_before(im,j ,kp,n)]  += hxm*hy *hzp*mass;
	      dens_in[fftw_index_before(im,j ,k ,n)]  += hxm*hy *hz *mass;
	      dens_in[fftw_index_before(im,j ,km,n)]  += hxm*hy *hzm*mass;
	      dens_in[fftw_index_before(im,jm,kp,n)]  += hxm*hym*hzp*mass;
	      dens_in[fftw_index_before(im,jm,k ,n)]  += hxm*hym*hz *mass;
	      dens_in[fftw_index_before(im,jm,km,n)]  += hxm*hym*hzm*mass;

	    }
	  else
	    {
	      /*Bin everything*/
	      /* 		if (im < 0)  */
	      /* 		    fprintf(stderr, "last bin: im < 0, %f \n", xnorm);  */
	      /* 		if (ip >= local_nx) */
	      /* 		    fprintf(stderr, "last bin: ip >= local_nx, %f \n", xnorm);  */

	      dens_in[fftw_index_before(ip,jp,kp,n)]  += hxp*hyp*hzp*mass;
	      dens_in[fftw_index_before(ip,jp,k ,n)]  += hxp*hyp*hz *mass;
	      dens_in[fftw_index_before(ip,jp,km,n)]  += hxp*hyp*hzm*mass;
	      dens_in[fftw_index_before(ip,j ,kp,n)]  += hxp*hy *hzp*mass;
	      dens_in[fftw_index_before(ip,j ,k ,n)]  += hxp*hy *hz *mass;
	      dens_in[fftw_index_before(ip,j ,km,n)]  += hxp*hy *hzm*mass;
	      dens_in[fftw_index_before(ip,jm,kp,n)]  += hxp*hym*hzp*mass;
	      dens_in[fftw_index_before(ip,jm,k ,n)]  += hxp*hym*hz *mass;
	      dens_in[fftw_index_before(ip,jm,km,n)]  += hxp*hym*hzm*mass;

	      dens_in[fftw_index_before(i ,jp,kp,n)]  += hx *hyp*hzp*mass;
	      dens_in[fftw_index_before(i ,jp,k ,n)]  += hx *hyp*hz *mass;
	      dens_in[fftw_index_before(i ,jp,km,n)]  += hx *hyp*hzm*mass;
	      dens_in[fftw_index_before(i ,j ,kp,n)]  += hx *hy *hzp*mass;
	      dens_in[fftw_index_before(i ,j ,k ,n)]  += hx *hy *hz *mass;
	      dens_in[fftw_index_before(i ,j ,km,n)]  += hx *hy *hzm*mass;
	      dens_in[fftw_index_before(i ,jm,kp,n)]  += hx *hym*hzp*mass;
	      dens_in[fftw_index_before(i ,jm,k ,n)]  += hx *hym*hz *mass;
	      dens_in[fftw_index_before(i ,jm,km,n)]  += hx *hym*hzm*mass;

	      dens_in[fftw_index_before(im,jp,kp,n)]  += hxm*hyp*hzp*mass;
	      dens_in[fftw_index_before(im,jp,k ,n)]  += hxm*hyp*hz *mass;
	      dens_in[fftw_index_before(im,jp,km,n)]  += hxm*hyp*hzm*mass;
	      dens_in[fftw_index_before(im,j ,kp,n)]  += hxm*hy *hzp*mass;
	      dens_in[fftw_index_before(im,j ,k ,n)]  += hxm*hy *hz *mass;
	      dens_in[fftw_index_before(im,j ,km,n)]  += hxm*hy *hzm*mass;
	      dens_in[fftw_index_before(im,jm,kp,n)]  += hxm*hym*hzp*mass;
	      dens_in[fftw_index_before(im,jm,k ,n)]  += hxm*hym*hz *mass;
	      dens_in[fftw_index_before(im,jm,km,n)]  += hxm*hym*hzm*mass;

	    }


	}
      else if(rank == 0 && i== n-1 && nproc > 1)
	{
	  jp = my_modulo(n,j+1);
	  jm = my_modulo(n,j-1);

	  kp = my_modulo(n,k+1);
	  km = my_modulo(n,k-1);


	  /*wrap around to top of box*/
	  dens_in[fftw_index_before(0 ,jp,kp,n)]  += hx *hyp*hzp*mass;
	  dens_in[fftw_index_before(0 ,jp,k ,n)]  += hx *hyp*hz *mass;
	  dens_in[fftw_index_before(0 ,jp,km,n)]  += hx *hyp*hzm*mass;
	  dens_in[fftw_index_before(0 ,j ,kp,n)]  += hx *hy *hzp*mass;
	  dens_in[fftw_index_before(0 ,j ,k ,n)]  += hx *hy *hz *mass;
	  dens_in[fftw_index_before(0 ,j ,km,n)]  += hx *hy *hzm*mass;
	  dens_in[fftw_index_before(0 ,jm,kp,n)]  += hx *hym*hzp*mass;
	  dens_in[fftw_index_before(0 ,jm,k ,n)]  += hx *hym*hz *mass;
	  dens_in[fftw_index_before(0 ,jm,km,n)]  += hx *hym*hzm*mass;
	}
      else if(rank == nproc-1 && i==0 && nproc > 1)
	{
	  jp = my_modulo(n,j+1);
	  jm = my_modulo(n,j-1);

	  kp = my_modulo(n,k+1);
	  km = my_modulo(n,k-1);


	  /*wrap around to bottom of box*/
	  dens_in[fftw_index_before(local_nx-1,jp,kp,n)]  += hx *hyp*hzp*mass;
	  dens_in[fftw_index_before(local_nx-1,jp,k ,n)]  += hx *hyp*hz *mass;
	  dens_in[fftw_index_before(local_nx-1,jp,km,n)]  += hx *hyp*hzm*mass;
	  dens_in[fftw_index_before(local_nx-1,j ,kp,n)]  += hx *hy *hzp*mass;
	  dens_in[fftw_index_before(local_nx-1,j ,k ,n)]  += hx *hy *hz *mass;
	  dens_in[fftw_index_before(local_nx-1,j ,km,n)]  += hx *hy *hzm*mass;
	  dens_in[fftw_index_before(local_nx-1,jm,kp,n)]  += hx *hym*hzp*mass;
	  dens_in[fftw_index_before(local_nx-1,jm,k ,n)]  += hx *hym*hz *mass;
	  dens_in[fftw_index_before(local_nx-1,jm,km,n)]  += hx *hym*hzm*mass;
	}


	    
	
    }



  return(0);

}


/* int grid_dist_D12(int rank, int nproc, int n, int NumPart, int TotNumPart, float BoxSize, float mass, struct particle_data *P,  */
/* 		  int local_nx, int local_x_start,fftw_real *dens_in) */
/* { */
/*   int i, j, k,l, index, flag = 0; */
/*   float rhobar; */
/*   double xnorm,ynorm,znorm, deltax, deltay, deltaz, dx, dy, dz; */
/*   double hx, hxp, hxm, hy, hyp, hym, hz, hzp, hzm; */
/*   int ip, im, jp, jm, kp, km; */
/*   FILE *fp; */
/*   int check_in_box;  */


/*   int n2 = 2*(n/2+1); */
/*   fftw_real *dens_p;  */
/*   fftw_real *dens_m;  */

/*   int npts = 694; */
/*   double *x12 = malloc(npts*sizeof(double));  */
/*   double *s12 = malloc(npts*sizeof(double));  */

/*   fp = fopen("/home/562/jxk562/powerspectra/D12.txt","r"); */
/*   for (i = 0; i < npts; i++) */
/*     fscanf(fp, "%lf %lf", &(x12[i]), &(s12[i]));  */
/*   fclose(fp);  */

/*   gsl_interp_accel *acc  = gsl_interp_accel_alloc (); */
/*   gsl_spline *spline     = gsl_spline_alloc (gsl_interp_cspline, npts); */

/*   gsl_spline_init (spline, x12, s12, npts); */


/*   /\* Need to place galaxies and random sample on to a grid *\/  */

/*   deltax = BoxSize/n; */
/*   deltay = BoxSize/n; */
/*   deltaz = BoxSize/n; */



/*   for (l=0; l < NumPart; l++) */
/*     { */
/*       /\* coordinates in subcube units [0,NGRID] *\/ */
/*       xnorm = P[l].Pos[0]/deltax; */
/*       ynorm = P[l].Pos[1]/deltay; */
/*       znorm = P[l].Pos[2]/deltaz; */


/*       int ni, nj, nk; 	/\* indices of nearest grid point [0,NGRID] *\/ */
/*       ni  = (int)(xnorm+0.5); */
/*       nj  = (int)(ynorm+0.5); */
/*       nk  = (int)(znorm+0.5); */


/*       dx = xnorm-(double)i;               */
/*       dy = ynorm-(double)j; */
/*       dz = znorm-(double)k; */

	
/*       float  xc, yc, zc; /\* 	/\\*Distance of particle from cell center*\\/ *\/ */
/*       xc=xnorm-floor(xnorm)-0.5; */
/*       yc=ynorm-floor(ynorm)-0.5; */
/*       zc=znorm-floor(znorm)-0.5;  */


/*       int nni, nnj, nnk; /\*next nearest grid points i.e. vertices of cell surrounding particle *\/ */
/*       int xp, xm, yp, ym, zp, zm; /\*Limits on indicies*\/ */
/*       /\*This would give you the surrounding cell but we want the surrounding 6^3 grid points*\/ */
/*       /\* 	if (xc > 0) *\/                    */
/*       /\* 	  nni = ni+1;  *\/ */
/*       /\* 	else *\/ */
/*       /\* 	  nni = ni-1;  *\/ */

/*       /\* 	if (yc > 0) *\/ */
/*       /\* 	  nnj = nj+1;  *\/ */
/*       /\* 	else *\/ */
/*       /\* 	  nnj = nj-1;  *\/ */

/*       /\* 	if (xz > 0) *\/ */
/*       /\* 	  nnk = nk+1;  *\/ */
/*       /\* 	else *\/ */
/*       /\* 	  nnk = nk-1;  *\/ */

/*       if (xc > 0) */
/* 	{ */
/* 	  nni = ni+3;  */
/* 	  ni  = ni-2;  */
/* 	  xp = nni;  */
/* 	  xm = ni;  */
/* 	} */
/*       else */
/* 	{ */
/* 	  nni = ni-3;  */
/* 	  ni = ni +2;  */
/* 	  xp = ni;  */
/* 	  xm = nni;  */
/* 	} */
/*       if (yc > 0) */
/* 	{ */
/* 	  nnj = nj+3;  */
/* 	  nj = nj -2;  */
/* 	  yp = nnj;  */
/* 	  ym = nj;  */
/* 	} */
/*       else */
/* 	{ */
/* 	  nnj = nj-3;  */
/* 	  nj = nj+2;  */
/* 	  yp = nj;  */
/* 	  ym = nnj;  */
/* 	} */

/*       if (zc > 0) */
/* 	{ */
/* 	  nnk = nk+3;  */
/* 	  nk = nk -2;  */
/* 	  zp = nnk;  */
/* 	  zm = nk;  */
/* 	} */
/*       else */
/* 	{ */
/* 	  nnk = nk-3;  */
/* 	  nk = nk+2;  */
/* 	  zp = nk;  */
/* 	  zm = nnk;  */
/* 	} */
	  


/*       for (i = xm; i < xp; i++) */
/* 	{ */
/* 	  int check_in_box = i-local_x_start;  */
/* 	  if (check_in_box >= 0 && check_in_box < local_nx)  */
/* 	    { */
/* 	      for (j = ym; j < yp; j++) */
/* 		for (k = zm;  k < zp; k++) */
/* 		  { */


/* 		    /\*Need to account for periodic boundary conditions...*\/ */
/* 		    /\*MPI: Accept all particles in a slab of thickness 6 grid points in the x direction*\/ */


/* 		    if(i >= n) i=i-n;             /\*Wrap around for period boundary conditions*\/ */
/* 		    if(j >= n) j=j-n; */
/* 		    if(k >= n) k=k-n; */

/* 		    if(i < 0) i=i+n;              */
/* 		    if(j < 0) j=j+n; */
/* 		    if(k < 0) k=k+n; */

/* 		    hx = fabs(my_modulo(n,xnorm - i));   /\*Distance  between the current grid point (i,j,k) and the particle (xnorm, ynorm, znorm) *\/ */
/* 		    hy = fabs(my_modulo(n,ynorm - j)); */
/* 		    hz = fabs(my_modulo(n,znorm - k)); */

/* 		    hx = wrap_around(hx, (float)xnorm, (float)i, BoxSize);  /\*Shortest distance between grid point and particle*\/ */
/* 		    hy = wrap_around(hy, (float)ynorm, (float)j, BoxSize);   */
/* 		    hz = wrap_around(hz, (float)znorm, (float)k, BoxSize);  */

/* 		    double wx,wy,wz;  /\*weights in each direction*\/ */
	
/* 		    wx = gsl_spline_eval(spline, hx, acc); */
/* 		    wy = gsl_spline_eval(spline, hy, acc); */
/* 		    wz = gsl_spline_eval(spline, hz, acc); */


/* 		    dens_in[fftw_index_before(i,j,k,n)]  += wx*wy*wz*mass; */
/* 		  } */
/* 	    } */
/* 	} */
/*     } */
/*   gsl_spline_free(spline);  */
/*   gsl_interp_accel_free (acc); */
/*   free(x12); free(s12);  */

/*   return(0); */

/* } */



int finish_grid(int n, int local_nx, int TotNumPart, float BoxSize, fftw_real *dens_in)
{

  /*   float rhobar = TotNumPart/(BoxSize*BoxSize*BoxSize); */
  float rhobar = TotNumPart/((float)n*(float)n*(float)n);


  /* Calculate the density contrast */
  int i, j, k; 

  for (i=0; i<local_nx; i++) 
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
	dens_in[fftw_index_before(i,j,k,n)] = dens_in[fftw_index_before(i,j,k,n)]/rhobar-1;  


  return(0); 
}

int logtransform(int n, int local_nx, int TotNumPart, float BoxSize, fftw_real *dens_in)
{

  /*    FILE *fp = fopen("/short/u95/jkwan/snapshots/loggedtest_zspace.dat", "w"); */

  /*   float rhobar = TotNumPart/(BoxSize*BoxSize*BoxSize); */
  float rhobar = TotNumPart/((float)n*(float)n*(float)n);

  /*   fwrite(&n, sizeof(int), 3, fp); */

  /* Calculate the density contrast */
  int i, j, k; 


  for (i=0; i<local_nx; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
	{
	  dens_in[fftw_index_before(i,j,k,n)] = dens_in[fftw_index_before(i,j,k,n)]/rhobar-1;
	  if (dens_in[fftw_index_before(i,j,k,n)]!=0)
	    dens_in[fftw_index_before(i,j,k,n)] = log((double)dens_in[fftw_index_before(i,j,k,n)]+1);

	  /* 	  fwrite(&(dens_in[fftw_index_before(i,j,k,n)]), sizeof(float), 1, fp); */
	  /* 	  fprintf(fp, "%f \n", dens_in[fftw_index_before(i,j,k,n)]); */
	}

  /*     fclose(fp); */


  return(0); 
}


int swap_densities(int rank, int nproc, int ngrid, int local_nx, fftw_real *dens_m, fftw_real *dens_p, fftw_real *dens_in)
{
  int rankup, rankdown; 


  /* lowest rank -> highest rank */

  /* 0 -> 1 -> 2 ....-> nproc-1   */
  
  if (rank == 0)
    rankdown = nproc - 1;  /*Wrap onto last processor if you're rank 0 */
  else
    rankdown = rank - 1; 

  if (rank == nproc - 1) /* Wrap onto first processor if you're the highest ranked processor */ 
    rankup = 0; 
  else 
    rankup = rank + 1; 



  fftw_real *dens_buf_p, *dens_buf_m, *dens_buf; 

  MPI_Request request; 
    
  MPI_Status status; 

  int ngrid2 = 2*(ngrid/2+1);

  if(!(dens_buf = calloc(ngrid*ngrid2, sizeof(fftw_real))))
    {
      fprintf(stderr, "Unable to allocate memory to dens_buf_p \n");
      exit(1);
    }
  
 
  int j, k;
  double any; 


  /*Sending densities to next processor ... one row at a time*/
  MPI_Isend(dens_p, ngrid*ngrid2*sizeof(fftw_real), MPI_BYTE, rankup, 1, MPI_COMM_WORLD, &request);

  MPI_Recv(dens_buf, ngrid*ngrid2*sizeof(fftw_real), MPI_BYTE, rankdown, 1, MPI_COMM_WORLD, &status);

  MPI_Wait (&request, &status);

  for (j=0; j < ngrid; j++)
    for (k=0; k < ngrid; k++)
      dens_in[fftw_index_before(0,j,k,ngrid)] += dens_buf[fftw_index_before(0,j,k,ngrid)];


  int count; 
  MPI_Get_count(&status, MPI_BYTE, &count);

  /* 	/\*Sending densities to previous processor ... *\/ */
  MPI_Isend(dens_m, ngrid*ngrid2*sizeof(fftw_real), MPI_BYTE, rankdown, 0, MPI_COMM_WORLD, &request);

  /*Receiving (blocking) densities from next processor ... */
  MPI_Recv(dens_buf, ngrid*ngrid2*sizeof(fftw_real), MPI_BYTE, rankup, 0, MPI_COMM_WORLD, &status);

  MPI_Wait (&request, &status);
  /*     } */
  MPI_Get_count(&status, MPI_BYTE, &count);

  for (j=0; j < ngrid; j++)
    for (k=0; k < ngrid; k++)
      dens_in[fftw_index_before(local_nx-1,j,k,ngrid)] += dens_buf[fftw_index_before(0,j,k,ngrid)];


  free(dens_buf); 
    

  return(0); 
}
int get_pk_new(int me, int Ng, int Np, int local_y_start_after_transpose, int local_ny_after_transpose, float BoxSize, fftw_complex *dgrid, double *pkall, char filename[])
{
    int		xlim[2],ylim[2],fftsize,ibin=0;
    int		leftpart,rightpart,*lcnt,*acnt;
    float		k2,d2,kmax2,scale,dx,dy,dz; //Nbin = 1000
    double		ng,mdens,totm,B,tPi2,wk[3],*pk;
    int		Ng2,nn,ix,iy,iz,ii[3],ip;
    int             Nbin = Ng; 
    MPI_Status	status;

    kmax2 = 2*M_PI*0.8*Ng/2; 
    kmax2 = kmax2*kmax2; 
  mdens = (double)Np/pow((double)Ng,3);
  scale = 1.0/pow((double)Ng,3.0)/mdens;
  scale = scale*scale;
  ylim[0] = local_y_start_after_transpose; 
  ylim[1] = local_ny_after_transpose;
  tPi2 = 2.0*M_PI*2.0*M_PI; 
  B = 1./Ng; 
  
  lcnt = malloc(  Nbin*sizeof(int));
  acnt = malloc(  Nbin*sizeof(int));
  pk   = malloc(2*Nbin*sizeof(double));
  for (ix=0; ix<  Nbin; ix++) lcnt[ix]=acnt[ix]=0;
  for (ix=0; ix<2*Nbin; ix++) pk[ix]=pkall[ix]=0;


  for (iy=ylim[0]; iy<ylim[0]+ylim[1]; iy++)
    for (ix=0; ix<Ng; ix++)
      for (iz=0; iz<Ng/2+1; iz++) {
        if (ix>Ng/2)
          ii[0] = ix - Ng;
        else
          ii[0] = ix;
        if (iy>Ng/2)
          ii[1] = iy - Ng;
        else
          ii[1] = iy;
        if (iz>Ng/2)
          ii[2] = iz - Ng;
        else
          ii[2] = iz;
        k2  = ii[0]*ii[0]+
              ii[1]*ii[1]+
              ii[2]*ii[2];
        k2 *= tPi2;
        if (k2>0 && k2<kmax2) {
          ip   = Ng*(Ng/2+1)*(iy-ylim[0])+(Ng/2+1)*ix+iz;
          d2   = dgrid[ip].re*dgrid[ip].re+dgrid[ip].im*dgrid[ip].im;
          d2  *= scale;
          wk[0]= j0(2*M_PI*ii[0]*B/2);
          wk[1]= j0(2*M_PI*ii[1]*B/2);
          wk[2]= j0(2*M_PI*ii[2]*B/2);
          d2  /= wk[0]*wk[0]*wk[1]*wk[1]*wk[2]*wk[2];
          d2  /= wk[0]*wk[0]*wk[1]*wk[1]*wk[2]*wk[2];
/*           d2  -= 1.0/totm; */

//          d2  *= k2*sqrt(k2)/2/M_PI/M_PI;
/*           ibin = Nbin*log(k2/tPi2)/(log(kmax2/tPi2))+0.49999; */
          ibin = sqrt((float)k2/tPi2)/1.0;
	  //	  fprintf(stderr, "%f %f \n", d2, k2); 
          if(ibin>=0 && ibin<Nbin) {
            if (iz>0 || (iz==0 && ii[0]>0) || (iz==0 && ii[0]==0 && ii[1]>0)) {
              pk[2*ibin+0] += sqrt(k2);
              pk[2*ibin+1] += d2;
              lcnt[ibin]++;
            }
          }
        }
      }

  MPI_Reduce(lcnt,acnt, Nbin,MPI_INT  ,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(pk,pkall,2*Nbin,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if (me==0) {
    for (ibin=0; ibin<Nbin; ibin++)
      //      fprintf(stderr, "%f %f \n", pkall[2*ibin+0], pkall[2*ibin+1]);  
      if (acnt[ibin]>0) {
        pkall[2*ibin+0] /= acnt[ibin];
        pkall[2*ibin+1] /= acnt[ibin];
	//	fprintf(stderr, "%f %f \n", pk[2*ibin+1], pkall[2*ibin+1]); 
      }
      else {
        pkall[2*ibin+0] = -1;
        pkall[2*ibin+1] = -1;
      }

/*       FILE *fpoutnew = fopen("/global/u1/j/jkwan/powerspectra/pk_dd_M000_256_test_new_1.0000.dat","w"); */
      FILE *fpoutnew = fopen(filename,"w");
      for (ibin=0; ibin<Nbin; ibin++)
	if (pkall[2*ibin+0]>0)
	  fprintf(fpoutnew, "%12.4e %12.4e\n",pkall[2*ibin+0]/BoxSize,pkall[2*ibin+1]*pow(BoxSize,3.));
      fflush(fpoutnew);
      fclose(fpoutnew); 

  }

/* 	if (rank == 0) */
/* 	  { */
/* 	    int Nbin = 1000;  */
/* 	    int ibin;  */
/* 	    FILE *fpoutnew = fopen("/global/u1/j/jkwan/powerspectra/pk_dd_M000_256_test_new_1.0000.dat","w"); */
/* 	    for (ibin=0; ibin<Nbin; ibin++) */
/* 	      if (pkall[2*ibin+0]>0) */
/* 		fprintf(fpoutnew, "%12.4e %12.4e\n",pkall[2*ibin+0]/BoxSize,pkall[2*ibin+1]); */
/* 	    fflush(fpoutnew); */
/* 	    fclose(fpoutnew);  */
/* 	  } */

  return(0); 
}

double	j0(double x)
/* The zeroth order spherical Bessel function. */
{
  if (x>0.05)
    return( sin(x)/x );
  else
    return( 1 + x*x*(-1./6.+x*x*(1./120.-x*x/5040)) );
}


int get_pk(int ngrid, int NumPart, int local_ny, int local_y_start, float BoxSize, fftw_complex *dens_out, double *dpower, int *ipower)
{
  int i,j,k, l,ks, *kmod;
  float ks_min, ks_max, dk;
  double kx, ky, kz, kR, kmodulus, fcorrect,scorrect, window;

  /* kmax is given by the Nyquist frequency: k < k_n = 2pi/L n/2 */
  /* kmin is limited by the size of the box: k > 2pi/L */
  /* and the spacing of the grid needs to be L/n */

  ks_max = M_PI/BoxSize * ngrid;
  ks_min = 2.*M_PI/BoxSize;
  dk = 2*M_PI/BoxSize;   
  /*  /FFTnorm * pow3(io.header.boxsize)/FFTnorm;  */
  float number_dens = NumPart/powf(BoxSize,3);
  /*   float fftw_norm = powf(ngrid/BoxSize,3)/powf(ngrid,3.);   /\* original *\/ */
  float fftw_norm = powf(BoxSize/ngrid,3)/powf(ngrid,3.);
  float FFTnorm = (float)ngrid*(float)ngrid*(float)ngrid;
  /*The data is now transposed*/
  
  for (j = 0; j < local_ny; ++j)
    for (i = 0; i < ngrid; ++i)
      for (k = 0; k < (ngrid/2+1); ++k)
	{ 
	  /* Assign each shell an index, ks, the spherical wave number*/
	  ks = (int) floor( sqrt( (float)pow(i,2) +
				  (float)pow(j+local_y_start,2) +
				  (float)pow(k,2)));

/* 	  double power = pow(c_re(dens_out[fftw_index_after(j,i,k,ngrid)]),2)+ pow(c_im(dens_out[fftw_index_after(j,i,k,ngrid)]),2); */

	  c_re(dens_out[fftw_index_after(j,i,k,ngrid)]) = pow(c_re(dens_out[fftw_index_after(j,i,k,ngrid)]),2)+
	    pow(c_im(dens_out[fftw_index_after(j,i,k,ngrid)]),2);

	  c_im(dens_out[fftw_index_after(j,i,k,ngrid)]) = 0.0;
	  fcorrect = 1.; scorrect = 1.;

	  kx  = M_PI*(double)i/(double)ngrid; /*Missing a factor of 2?*/
	  ky  = M_PI*(double)(j+local_y_start)/(double)ngrid;
	  kz  = M_PI*(double)k/(double)ngrid;

	  if (kx > 0) fcorrect *= sin(kx)/kx;
	  if (ky > 0) fcorrect *= sin(ky)/ky;
	  if (kz > 0) fcorrect *= sin(kz)/kz;

	  scorrect *= (1. - pow(sin(kx),2.) + 2./15.*pow(sin(kx),4.));
	  scorrect *= (1. - pow(sin(ky),2.) + 2./15.*pow(sin(ky),4.));
	  scorrect *= (1. - pow(sin(kz),2.) + 2./15.*pow(sin(kz),4.));

	  fcorrect = 1./pow(fcorrect,4.);


	  if(ks >= 0 && ks < ngrid/2.)  /*A factor of dk has been taken out */
	    {
	      dpower[ks] += fcorrect*(double)c_re(dens_out[fftw_index_after(j,i,k,ngrid)])*fftw_norm;
	      ipower[ks] += 1;
	    }
	}  


    
  return(0);

}

int get_pk_binned(int ngrid, int NumPart, int local_ny, int local_y_start, float BoxSize, fftw_complex *dens_out, double *dpower, int *ipower)
{
  int i,j,k, l,ks, *kmod, nk;
  float ks_min, ks_max, dk, ks_step;
  double kx, ky, kz, kR, kmodulus, fcorrect,scorrect, window;

  /* kmax is given by the Nyquist frequency: k < k_n = 2pi/L n/2 */
  /* kmin is limited by the size of the box: k > 2pi/L */
  /* and the spacing of the grid needs to be L/n */
  nk = 20; 
  ks_max = log10(0.1); // log10(0.1)
  ks_min = -2; // log10(0.01)
  ks_step = (ks_max-ks_min)/(float)(nk-1);
  //  ks_step = 0.1; // ten bins between k = 0.01 and k = 0.1
  dk = 2*M_PI/BoxSize;   
  //dk = 2*M_PI/ngrid;   

  float number_dens = NumPart/powf(BoxSize,3);
  float fftw_norm = powf(ngrid/BoxSize,3)/powf(ngrid,3.);
  /*The data are now transposed*/
  
  for (j = 0; j < local_ny; ++j)
    for (i = 0; i < ngrid; ++i)
      for (k = 0; k < (ngrid/2+1); ++k)
	{ 
	  /* Assign each shell an index, ks, the spherical wave number*/
	  ks = (int) floor( sqrt( (float)pow(i,2) +
				  (float)pow(j+local_y_start,2) +
				  (float)pow(k,2)));

	  c_re(dens_out[fftw_index_after(j,i,k,ngrid)]) = pow(c_re(dens_out[fftw_index_after(j,i,k,ngrid)]),2)+
	    pow(c_im(dens_out[fftw_index_after(j,i,k,ngrid)]),2);

	  c_im(dens_out[fftw_index_after(j,i,k,ngrid)]) = 0.0;
	  fcorrect = 1.; scorrect = 1.;

	  kx  = M_PI*(double)i/(double)ngrid; /*Missing a factor of 2?*/
	  ky  = M_PI*(double)(j+local_y_start)/(double)ngrid;
	  kz  = M_PI*(double)k/(double)ngrid;

	  if (kx > 0) fcorrect *= sin(kx)/kx;
	  if (ky > 0) fcorrect *= sin(ky)/ky;
	  if (kz > 0) fcorrect *= sin(kz)/kz;

#ifdef TSC

	  fcorrect = 1./pow(fcorrect,6.);  /* Check that the index is correct */

	  scorrect *= (1. - pow(sin(kx),2.) + 2./15.*pow(sin(kx),4.));
	  scorrect *= (1. - pow(sin(ky),2.) + 2./15.*pow(sin(ky),4.));
	  scorrect *= (1. - pow(sin(kz),2.) + 2./15.*pow(sin(kz),4.));
#else 
	  fcorrect = 1./pow(fcorrect,2.);

#endif

	  int index = (int)((log10(ks*dk)-ks_min)/ks_step+0.5);
	  if(index >= 0 && index < (int)((ks_max-ks_min)/ks_step+0.5))  /*A factor of dk has been taken out */
	    {
	      dpower[index] += fcorrect*(double)c_re(dens_out[fftw_index_after(j,i,k,ngrid)])*fftw_norm;/*     - fcorrect*scorrect/number_dens;     */
	      ipower[index] += 1;
	    }
	}  


    
  return(0);

}


int get_pk_para_perp(int ngrid, int NumPart, int nmu, float *theta, int local_ny, int local_y_start, float BoxSize, 
		     fftw_complex *dens_out, int los, double *dpower, int *ipower)
{
  int i,j,k, l,ks, *kmod;
  float ks_min, ks_max, dk, idk;
  double kx, ky, kz, fcorrect,scorrect;
  int mu;


  /* kmax is given by the Nyquist frequency: k < k_n = 2pi/L n/2 */
  /* kmin is limited by the size of the box: k > 2pi/L */
  /* and the spacing of the grid needs to be L/n */

  ks_max = M_PI/BoxSize * ngrid;
  ks_min = 2.*M_PI/BoxSize;
  dk = 2*M_PI/BoxSize;
  idk = 0.01;    

  int nk = 50; 
  float k_max = 0.1; // log10(0.1)
  float k_min = 0; // log10(0.01)
  float k_step = (k_max-k_min)/(float)(nk-1);

  float number_dens = NumPart/powf(BoxSize,3);
  /*   float fftw_norm = powf(ngrid/BoxSize,3)/powf(ngrid,3.); */
  float fftw_norm = powf(BoxSize/ngrid,3.)/powf(ngrid,3.);

  /*The data are now transposed*/
  
  for (j = 0; j < local_ny; ++j)
    for (i = 0; i < ngrid; ++i)
      for (k = 0; k < (ngrid/2+1); ++k)
	{ 
	  /* Assign each shell an index, ks, the spherical wave number*/
	  ks = (int) floor( sqrt( (float)pow(i,2) +
				  (float)pow(j+local_y_start,2) +
				  (float)pow(k,2)));

	  c_re(dens_out[fftw_index_after(j,i,k,ngrid)]) = pow(c_re(dens_out[fftw_index_after(j,i,k,ngrid)]),2)+
	    pow(c_im(dens_out[fftw_index_after(j,i,k,ngrid)]),2);
	  c_im(dens_out[fftw_index_after(j,i,k,ngrid)]) = 0.0;
	  fcorrect = 1.; scorrect = 1.;

	  kx  = M_PI*(double)i/(double)ngrid;
	  ky  = M_PI*(double)(j+local_y_start)/(double)ngrid;
	  kz  = M_PI*(double)k/(double)ngrid;

	  if (kx > 0) fcorrect *= sin(kx)/kx;
	  if (ky > 0) fcorrect *= sin(ky)/ky;
	  if (kz > 0) fcorrect *= sin(kz)/kz;

#ifdef TSC

	  fcorrect = 1./pow(fcorrect,6.);  

	  scorrect *= (1. - pow(sin(kx),2.) + 2./15.*pow(sin(kx),4.));
	  scorrect *= (1. - pow(sin(ky),2.) + 2./15.*pow(sin(ky),4.));
	  scorrect *= (1. - pow(sin(kz),2.) + 2./15.*pow(sin(kz),4.));
#else 
	  fcorrect = 1./pow(fcorrect,2.);
#endif


	  /*mu needs to be binned*/
	  /*dpower and ipower need to be made two dimensional*/


	   

	  int k_perp, k_para;
	  

	  if(ks > 0 && ks < ngrid/2.)  /*A factor of dk has been taken out */
	    {

	      float root2 = 1./sqrt(2.);
	      float tmp; 
	      int n;


	      if (los == 0)
		{
		  /* 		  k_para = (float)i; */
		  /* 		  k_perp =(float)j*(float)j + (float)k*(float)k; */
		  /* 		  mu = angular_projection(i, ks, nmu, theta); */
		  k_para = i;
		  k_perp = (int)sqrt((float)j*(float)j + (float)k*(float)k);
		}
	      
	      else if (los == 1)
		{
		  /* 		  k_para = (float)j; */
		  /* 		  k_perp =(float)i*(float)i + (float)k*(float)k; */
		  /* 		  mu = angular_projection(j+local_y_start, ks, nmu, theta);  */
		  k_para = j;
		  k_perp = (int)sqrt((float)i*(float)i + (float)k*(float)k);
		}
	      
	      else
		{
		  /* 		  k_para =(float)k; */
		  /* 		  k_perp =(float)i*(float)i + (float)j*(float)j; */
		  /* 		  mu = angular_projection(k, ks, nmu, theta); */

		  k_para = k;
		  k_perp = (int)sqrt((float)i*(float)i + (float)j*(float)j);
		}
	      


	      /* 		    if ( tmp > 0.951) */
	      /* 			mu = 4; */
	      /* 		    else if (tmp > 0.809) */
	      /* 			mu = 3; */
	      /* 		    else if (tmp > 0.588) */
	      /* 			mu = 2; */
	      /* 		    else if (tmp > 0.309) */
	      /* 			mu = 1; */
	      /* 		    else */
	      /* 			mu = 0; */

	      /* 	      int index_perp = (int)((log10(k_perp*dk)-k_min)/k_step+0.5); */
	      /* 	      int index_para = (int)((log10(k_para*dk)-k_min)/k_step+0.5); */
	      int index_perp = k_perp;/* (int)(k_perp*dk/k_step); */
	      int index_para = k_para; /* (int)(k_para*dk/k_step); */
	      if(index_perp >= 0 && index_perp < nk && index_para >= 0 && index_para < nk)
		/* 	      if (index_perp+index_para*nk < (nk-1)+(nk-1)*nk && index_perp >=0 && index_para >=0) */
		{
		  dpower[index_perp+index_para*nk] += fcorrect*(double)c_re(dens_out[fftw_index_after(j,i,k,ngrid)])*fftw_norm;/*     - fcorrect*scorrect/number_dens;     */
		  ipower[index_perp+index_para*nk] += 1;
		}

	      /* 	      dpower[k_perp+k_para*ngrid] += fcorrect*(double)c_re(dens_out[fftw_index_after(j,i,k,ngrid)])*fftw_norm - fcorrect*scorrect/number_dens;     */
	      /* 	      ipower[k_perp+k_para*ngrid] += 1; */
	    }
	}  


    
  return(0);

}

int get_pk_new_redshift(int rank, int nproc, int n, int nmu, int NumPart, int local_y_start_after_transpose, int local_ny_after_transpose, float BoxSize, int los, fftw_complex*dgrid, char *output_fname1, char *output_fname2) 
{
    int i, j, k,l, index, flag = 0;
    float rhobar;
    double xnorm,ynorm,znorm, deltax, deltay, deltaz, dx, dy, dz;
    double hx, hxp, hxm, hy, hyp, hym, hz, hzp, hzm;
    int ip, im, jp, jm, kp, km;
    int check_in_box; 

    double		k2,k_para,k_perp,d2,kmax2;
    double		ng,mdens,totm,B,tPi2,wk[3],legendre2;
    int		Ng,Ng2,nn,ix,iy,iz,ii[3];

    
    ng = 0.9999999*n; 
    Ng = n; 
    Ng2 = n+2; 
    int count=0,overload[2]; 

    B    = 1.0/Ng;
    kmax2= 2*M_PI*(Ng/2) * 0.8;
    kmax2= kmax2*kmax2;
    tPi2 = (2.*M_PI)*(2.*M_PI);

    /* Calculate mean density. */
    totm = NumPart;
    mdens= totm/pow((double)n,3.0);


    float scale = 1.0/pow((double)Ng,3.0)/mdens;
    scale = scale*scale;
    
    int Nbin = Ng; 
/*     int *lcnt = malloc(  Nbin*sizeof(int)); */
/*     int *acnt = malloc(  Nbin*sizeof(int)); */
/*     float *pk   = malloc(2*Nbin*sizeof(float));   */
/*     float *pkall= malloc(2*Nbin*sizeof(float)); */

    double *kmodes = malloc(Nbin*sizeof(double));
    double *kmodes_para = malloc(Nbin*sizeof(double));
    double *kmodes_perp = malloc(Nbin*sizeof(double));
    double *pk = malloc(Nbin*nmu*sizeof(double));
    double *pk_mono = malloc(Nbin*sizeof(double));  
    double *pk_quad = malloc(Nbin*sizeof(double));
    int *lcnt = malloc(Nbin*nmu*sizeof(int));
    int *lcnt_para= malloc(Nbin*sizeof(int));
    int *lcnt_perp= malloc(Nbin*sizeof(int));
    int *lcnt_multipole= malloc(Nbin*sizeof(int)); 

    double *kmodes_all = malloc(Nbin*sizeof(double));
    double *kmodes_para_all = malloc(Nbin*sizeof(double));
    double *kmodes_perp_all = malloc(Nbin*sizeof(double));
    double *pk_all = malloc(Nbin*nmu*sizeof(double));
    double *pk_mono_all = malloc(Nbin*sizeof(double));  
    double *pk_quad_all = malloc(Nbin*sizeof(double));
    int *lcnt_all = malloc(Nbin*nmu*sizeof(int));
    int *lcnt_para_all= malloc(Nbin*sizeof(int));
    int *lcnt_perp_all= malloc(Nbin*sizeof(int));
    int *lcnt_multipole_all= malloc(Nbin*sizeof(int)); 


    int ibin=0, ibin_para=0, ibin_perp=0;
    int mu; 

    for (ix=0; ix<  Nbin*nmu; ix++) 
    {
	lcnt[ix]=0;
        pk[ix]=0;
    }

    for (ix=0; ix<  Nbin; ix++) 
    {
	lcnt_multipole[ix]=0;
	lcnt_para[ix]=0;
	lcnt_perp[ix]=0;
	kmodes[ix]=0;
	kmodes_para[ix]=0;
	kmodes_perp[ix]=0;
	pk_mono[ix]=0;
	pk_quad[ix]=0;
    }

    for (iy=local_y_start_after_transpose; iy<local_y_start_after_transpose+local_ny_after_transpose; iy++)
	for (ix=0; ix<Ng; ix++)
	    for (iz=0; iz<Ng/2+1; iz++) {
		if (ix>Ng/2)
		    ii[0] = ix - Ng;
		else
		    ii[0] = ix;
		if (iy>Ng/2)
		    ii[1] = iy - Ng;
		else
		    ii[1] = iy;
		if (iz>Ng/2)
		    ii[2] = iz - Ng;
		else
		    ii[2] = iz;
		k2  = ii[0]*ii[0]+
		    ii[1]*ii[1]+
		    ii[2]*ii[2];
		k2 *= tPi2;
/* 	fprintf(stderr, "%d %d %d %f \n", ii[0], ii[1], ii[2], k2);  */
		if (k2>0 && k2<kmax2) {
		    ip   = Ng*(Ng/2+1)*(iy-local_y_start_after_transpose)+(Ng/2+1)*ix+iz;
		    d2   = dgrid[ip].re*dgrid[ip].re+dgrid[ip].im*dgrid[ip].im;
		    d2  *= scale;

		    double fcorrect = 1; 

		    double kx  = M_PI*(double)ii[0]/(double)n; 
		    double ky  = M_PI*(double)ii[1]/(double)n;
		    double kz  = M_PI*(double)ii[2]/(double)n;

		    if (kx != 0) fcorrect *= sin(kx)/kx;
		    if (ky != 0) fcorrect *= sin(ky)/ky;
		    if (kz != 0) fcorrect *= sin(kz)/kz;
		

		    wk[0]= j0(2*M_PI*ii[0]*B/2);
		    wk[1]= j0(2*M_PI*ii[1]*B/2);
		    wk[2]= j0(2*M_PI*ii[2]*B/2);
		    d2  -= 1.0/totm; 
		    d2  /= wk[0]*wk[0]*wk[1]*wk[1]*wk[2]*wk[2];
		    d2  /= wk[0]*wk[0]*wk[1]*wk[1]*wk[2]*wk[2];


//		    d2  *= k2*sqrt(k2)/2/M_PI/M_PI;
/* 		    ibin = Nbin*log(k2/tPi2)/(log(kmax2/tPi2))+0.49999; */
/* 		    ibin = 2*sqrt(k2/tPi2)/1.0; //original */
		    ibin = sqrt(k2/tPi2)/1.0;
/* 		    ibin = 4.0*sqrt((float)k2/tPi2)/1.0; */
/* 		    if(ibin>=0 && ibin<2*n) { */
		    if(ibin>=0 && ibin<Nbin) {
			if (iz>0 || (iz==0 && ii[0]>0) || (iz==0 && ii[0]==0 && ii[1]>0)) {

			    switch (los)
			    {
				case 0:
/* 				    ibin_para = (int)floor(sqrt(ii[0]*ii[0]+0.49999)); */
/* 				    ibin_perp = (int)floor(sqrt(ii[1]*ii[1]+ii[2]*ii[2])); */

				    ibin_para = sqrt(ii[0]*ii[0]*1.0)/1.0;
				    ibin_perp = sqrt(ii[1]*ii[1]*1.0+ii[2]*ii[2]*1.0)/1.0;
				    k_para = sqrt(tPi2*ii[0]*ii[0]); 
				    k_perp = 2.*M_PI*sqrt(ii[1]*ii[1]+ii[2]*ii[2]);
				    legendre2 = 2.5*(3.*(float)ii[0]*(float)ii[0]/(k2/tPi2)-1.);
   				    mu = (int)((float)(nmu-1)*(float)ibin_para/sqrt(k2/tPi2)+0.5);
				    break;
				case 1:
/* 				    ibin_para = (int)floor(sqrt(ii[1]*ii[1]+0.49999)); */
/* 				    ibin_perp = (int)floor(sqrt(ii[0]*ii[0]+ii[2]*ii[2])); */

				    ibin_para = sqrt(ii[1]*ii[1]*1.0)/1.0;
				    ibin_perp = sqrt(ii[0]*ii[0]*1.0+ii[2]*ii[2]*1.0)/1.0;
				    k_para = sqrt(tPi2*ii[1]*ii[1]); 
				    k_perp = 2.*M_PI*sqrt(ii[0]*ii[0]+ii[2]*ii[2]);
				    legendre2 = 2.5*(3.*(float)ii[1]*(float)ii[1]/(k2/tPi2)-1.);
   				    mu = (int)((float)(nmu-1)*(float)ibin_para/sqrt(k2/tPi2)+0.5);
				    break;
				case 2:
/* 				    ibin_para = (int)floor(sqrt(ii[2]*ii[2]+0.49999)); */
/* 				    ibin_perp = (int)floor(sqrt(ii[0]*ii[0]+ii[1]*ii[1])); */

				    ibin_para = sqrt(ii[2]*ii[2]*1.0)/1.0;
				    ibin_perp = sqrt(ii[0]*ii[0]*1.0+ii[1]*ii[1]*1.0)/1.0;
				    k_para = sqrt(tPi2*ii[2]*ii[2]); 
				    k_perp = 2.*M_PI*sqrt(ii[0]*ii[0]+ii[1]*ii[1]);
				    legendre2 = 2.5*(3.*(float)ii[2]*(float)ii[2]/(k2/tPi2)-1.);
   				    mu = (int)((float)(nmu-1)*(float)ibin_para/sqrt(k2/tPi2)+0.5);
				    break;
				default:
				    ibin_para = 0; ibin_perp = 0; 
				    break;
			    }


			    pk[mu+ibin*nmu] += d2;
			    kmodes_para[ibin_para]+=k_para; 
			    kmodes_perp[ibin_perp]+=k_perp; 
			    lcnt_para[ibin_para]++;
			    lcnt_perp[ibin_perp]++;
			    lcnt[mu+ibin*nmu]++;

			    kmodes[ibin] += sqrt(k2);
			    pk_mono[ibin] += d2;
			    pk_quad[ibin] += d2*legendre2;
			    lcnt_multipole[ibin]++; 



/* 			    pk[ibin_para+ibin_perp*nmu] += d2; */
/* 			    kmodes_para[ibin_para]+=k_para;  */
/* 			    kmodes_perp[ibin_perp]+=k_perp;  */
/* 			    lcnt_para[ibin_para]++; */
/* 			    lcnt_perp[ibin_perp]++; */
/* 			    lcnt[ibin_para+ibin_perp*nmu]++; */

/* 			    kmodes[ibin] += sqrt(k2); */
/* 			    pk_mono[ibin] += d2; */
/* 			    pk_quad[ibin] += d2*legendre2; */
/* 			    lcnt_multipole[ibin]++;  */
/* 			    fprintf(stderr, "Rank %d mine: %d %d %d %d %d %f %lf\n", rank, ix, iy, iz, ibin, ip, k2, dgrid[ip].re*dgrid[ip].re+dgrid[ip].im*dgrid[ip].im); */
			}
		    }
		}
	    }


    

    MPI_Reduce(kmodes_perp, kmodes_perp_all, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(kmodes_para, kmodes_para_all, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(pk, pk_all, n*nmu, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(lcnt, lcnt_all, n*nmu, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(kmodes, kmodes_all, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(pk_mono, pk_mono_all, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(pk_quad, pk_quad_all, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(lcnt_multipole, lcnt_multipole_all, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(lcnt_para, lcnt_para_all, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(lcnt_perp, lcnt_perp_all, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);



/*     sprintf(output_fname1, "/home/jkwan/powerspectra/pk_M%03d_%d_multipole_%.4f_z%d.dat",model_no,ngrid,header1.time,los); */
/*     sprintf(output_fname2, "/home/jkwan/powerspectra/pk_M%03d_%d_2D_%.4f_z%d.dat",model_no,ngrid,header1.time,los); */

	if (rank == 0)
	{
	  FILE *fp_multi = fopen(output_fname1, "w"); 
	  FILE *fp_2D = fopen(output_fname2, "w"); 


//	    fout = fopen(output_fname, "w");
	    if(!(fp_multi = fopen(output_fname1, "w")))
	    {
		printf("can't open file `%s`\n",output_fname1);
		exit(0);
	    }

	    if(!(fp_2D = fopen(output_fname2, "w")))
	    {
		printf("can't open file `%s`\n",output_fname2);
		exit(0);
	    }



	    for (i = 0; i < n; i++) /*Don't include DC mode*/
	    {

		kmodes_perp_all[i] /= lcnt_perp_all[i];
		kmodes_para_all[i] /= lcnt_para_all[i];
		pk_mono_all[i]/=lcnt_multipole_all[i];
		pk_quad_all[i]/=lcnt_multipole_all[i];
		kmodes_all[i] /= lcnt_multipole_all[i];
		if (lcnt_multipole_all[i] != 0)
		    fprintf(fp_multi, "%f %f %f \n", kmodes_all[i]/BoxSize, pk_mono_all[i]*pow(BoxSize,3.), pk_quad_all[i]*pow(BoxSize,3.));
	    }

	    for (i = 0; i < n; i++) /*Don't include DC mode*/
	      //	      if(lcnt_perp_all[i] != 0 && lcnt_para_all[i]!=0)
		if (lcnt_multipole_all[i] != 0)
		{
		  for (j = 0; j< nmu; j++)
		    {
		      if (lcnt_all[j+nmu*i]!=0)
			{
			  pk_all[j+nmu*i]/=lcnt_all[j+nmu*i];
			  float mu = kmodes_para_all[j]/sqrt(kmodes_perp_all[i]*kmodes_perp_all[i]+kmodes_para_all[j]*kmodes_para_all[j]);
			  fprintf(fp_2D, "%f %f %f \n", kmodes_all[i]/BoxSize, (float)j/(float)(nmu-1), pk_all[j+nmu*i]*pow(BoxSize,3.));
			  //			  fprintf(fp_2D, "%f %f %f %f \n", kmodes_para_all[j]/BoxSize, kmodes_perp_all[i]/BoxSize, mu, pk_all[j+nmu*i]*pow(BoxSize,3.));
			}
		    }
		}

		fclose(fp_multi);
		fclose(fp_2D);

/* 	    fflush(fout); */
/* 	    fclose(fout); */

	}





    free(kmodes);    free(kmodes_para);    free(kmodes_perp);
    free(pk);    free(pk_mono);    free(pk_quad);
    free(lcnt);    free(lcnt_para);    free(lcnt_perp);    free(lcnt_multipole); 


    free(kmodes_all);    free(kmodes_para_all);    free(kmodes_perp_all);
    free(pk_all);    free(pk_mono_all);    free(pk_quad_all);
    free(lcnt_all);    free(lcnt_para_all);    free(lcnt_perp_all);    free(lcnt_multipole_all); 



    return(0); 
}

int get_pk_redshift(int ngrid, int NumPart, int nmu, float *theta, int local_ny, int local_y_start, float BoxSize, 
		    fftw_complex *dens_out, int los, double *dpower, int *ipower)
{
  int i,j,k, l,ks, *kmod;
  float ks_min, ks_max, dk, idk;
  double kx, ky, kz, fcorrect,scorrect;
  int mu;


  /* kmax is given by the Nyquist frequency: k < k_n = 2pi/L n/2 */
  /* kmin is limited by the size of the box: k > 2pi/L */
  /* and the spacing of the grid needs to be L/n */

  ks_max = M_PI/BoxSize * ngrid;
  ks_min = 2.*M_PI/BoxSize;
  dk = 2*M_PI/BoxSize;
  idk = 0.01;    

  float number_dens = NumPart/powf(BoxSize,3);
  float fftw_norm = powf(BoxSize/ngrid,3.)/powf(ngrid,3.);

  /*The data are now transposed*/
  
  for (j = 0; j < local_ny; ++j)
    for (i = 0; i < ngrid; ++i)
      for (k = 0; k < (ngrid/2+1); ++k)
	{ 
	  /* Assign each shell an index, ks, the spherical wave number*/
	  ks = (int) floor( sqrt( (float)pow(i,2) +
				  (float)pow(j+local_y_start,2) +
				  (float)pow(k,2)));

	  c_re(dens_out[fftw_index_after(j,i,k,ngrid)]) = pow(c_re(dens_out[fftw_index_after(j,i,k,ngrid)]),2)+
	    pow(c_im(dens_out[fftw_index_after(j,i,k,ngrid)]),2);
	  c_im(dens_out[fftw_index_after(j,i,k,ngrid)]) = 0.0;
	  fcorrect = 1.; scorrect = 1.;

	  kx  = M_PI*(double)i/(double)ngrid;
	  ky  = M_PI*(double)(j+local_y_start)/(double)ngrid;
	  kz  = M_PI*(double)k/(double)ngrid;

	  if (kx > 0) fcorrect *= sin(kx)/kx;
	  if (ky > 0) fcorrect *= sin(ky)/ky;
	  if (kz > 0) fcorrect *= sin(kz)/kz;

#ifdef TSC

	  fcorrect = 1./pow(fcorrect,6.);  

	  scorrect *= (1. - pow(sin(kx),2.) + 2./15.*pow(sin(kx),4.));
	  scorrect *= (1. - pow(sin(ky),2.) + 2./15.*pow(sin(ky),4.));
	  scorrect *= (1. - pow(sin(kz),2.) + 2./15.*pow(sin(kz),4.));
#else 
	  fcorrect = 1./pow(fcorrect,2.);
#endif


	  /*mu needs to be binned*/
	  /*dpower and ipower need to be made two dimensional*/

	  if(ks > 0 && ks < ngrid/2.)  /*A factor of dk has been taken out */
	    {

	      float root2 = 1./sqrt(2.);
	      float tmp; 
	      int n;


	      if (los == 0)
		mu = angular_projection(i,ks, nmu, theta);
	      else if (los == 1)
		mu = angular_projection(j+local_y_start, ks, nmu, theta);
	      else
		mu = angular_projection(k, ks, nmu, theta);
		    


	      /* 	      if (los == 0) */
	      /* 		mu = equal_angle(i, j, k, ks, nmu); */
	      /* 	      else if (los == 1) */
	      /* 		mu = equal_angle(j+local_y_start, i, k, ks, nmu);/\* angular_projection(j+local_y_start, ks, nmu, theta);  *\/ */
	      /* 	      else */
	      /* 		mu = equal_angle(k,i,j,ks,nmu); /\* angular_projection(k, ks, nmu, theta); *\/ */
		    


	      /* 		    if ( tmp > 0.951) */
	      /* 			mu = 4; */
	      /* 		    else if (tmp > 0.809) */
	      /* 			mu = 3; */
	      /* 		    else if (tmp > 0.588) */
	      /* 			mu = 2; */
	      /* 		    else if (tmp > 0.309) */
	      /* 			mu = 1; */
	      /* 		    else */
	      /* 			mu = 0; */


	      dpower[mu+ks*nmu] += fcorrect*(double)c_re(dens_out[fftw_index_after(j,i,k,ngrid)])*fftw_norm;// - fcorrect*scorrect/number_dens;    
	      ipower[mu+ks*nmu] += 1;
	    }
	}  


    
  return(0);

}

int get_pk_redshift_binned(int ngrid, int NumPart, int nmu, float *theta, int local_ny, int local_y_start, float BoxSize, 
			   fftw_complex *dens_out, int los, double *dpower, int *ipower) 
{
  int i,j,k, l,ks, *kmod, nk;
  float ks_min, ks_max, dk, ks_step;
  double kx, ky, kz, kR, kmodulus, fcorrect,scorrect, window;
  int mu; 

  /* kmax is given by the Nyquist frequency: k < k_n = 2pi/L n/2 */
  /* kmin is limited by the size of the box: k > 2pi/L */
  /* and the spacing of the grid needs to be L/n */

  nk = 20; 
  ks_max =  log10(0.1); // log10(0.1)
  ks_min = log10(0.01); // log10(0.01)
  ks_step = (ks_max-ks_min)/(float)(nk-1);
  //ks_step = 0.1; // ten bins between k = 0.01 and k = 0.1
  dk = 2*M_PI/BoxSize;   
  //dk = 2*M_PI/ngrid;   

  float number_dens = NumPart/powf(BoxSize,3);
  float fftw_norm = powf(ngrid/BoxSize,3)/powf(ngrid,3.);
  /*The data are now transposed*/
  
  for (j = 0; j < local_ny; ++j)
    for (i = 0; i < ngrid; ++i)
      for (k = 0; k < (ngrid/2+1); ++k)
	{ 
	  /* Assign each shell an index, ks, the spherical wave number*/
	  ks = (int) floor( sqrt( (float)pow(i,2) +
				  (float)pow(j+local_y_start,2) +
				  (float)pow(k,2)));

	  c_re(dens_out[fftw_index_after(j,i,k,ngrid)]) = pow(c_re(dens_out[fftw_index_after(j,i,k,ngrid)]),2)+
	    pow(c_im(dens_out[fftw_index_after(j,i,k,ngrid)]),2);

	  c_im(dens_out[fftw_index_after(j,i,k,ngrid)]) = 0.0;
	  fcorrect = 1.; scorrect = 1.;

	  kx  = M_PI*(double)i/(double)ngrid; /*Missing a factor of 2?*/
	  ky  = M_PI*(double)(j+local_y_start)/(double)ngrid;
	  kz  = M_PI*(double)k/(double)ngrid;

	  if (kx > 0) fcorrect *= sin(kx)/kx;
	  if (ky > 0) fcorrect *= sin(ky)/ky;
	  if (kz > 0) fcorrect *= sin(kz)/kz;

#ifdef TSC

	  fcorrect = 1./pow(fcorrect,6.);  /* Check that the index is correct */

	  scorrect *= (1. - pow(sin(kx),2.) + 2./15.*pow(sin(kx),4.));
	  scorrect *= (1. - pow(sin(ky),2.) + 2./15.*pow(sin(ky),4.));
	  scorrect *= (1. - pow(sin(kz),2.) + 2./15.*pow(sin(kz),4.));
#else 
	  fcorrect = 1./pow(fcorrect,2.);

#endif

	  int index = (int)((log10(ks*dk)-ks_min)/ks_step+0.5);
	  if(index >= 0 && index < (int)((ks_max-ks_min)/ks_step+0.5))  /*A factor of dk has been taken out */
	    {


	      /* 	      if (los == 0) */
	      /* 		mu = equal_angle(i, j, k, ks, nmu); */
	      /* 	      else if (los == 1) */
	      /* 		mu = equal_angle(j+local_y_start, i, k, ks, nmu);/\* angular_projection(j+local_y_start, ks, nmu, theta);  *\/ */
	      /* 	      else */
	      /* 		mu = equal_angle(k,i,j,ks,nmu); /\* angular_projection(k, ks, nmu, theta); *\/ */


	      if (los == 0)
		mu = angular_projection(i,ks, nmu, theta);
	      else if (los == 1)
		mu = angular_projection(j+local_y_start, ks, nmu, theta);
	      else
		mu = angular_projection(k, ks, nmu, theta);
		    

	      dpower[mu+index*nmu] += fcorrect*(double)c_re(dens_out[fftw_index_after(j,i,k,ngrid)])*fftw_norm - fcorrect*scorrect/number_dens;    
	      ipower[mu+index*nmu] += 1;

	    }
	}  


    
  return(0);

}


int get_sigma_8_gal(int ngrid, int nmu, float BoxSize, int *ipower, double *dpower, double *sig8)
{
  int i; 

  int ik, itf, in;
  double kh, k, h, x, win, delta;
  double lnk, dlnk, lnko;
  double dsig8, dsig8o, sig8o, powers;
  double sigr8 = 8., dummy = 0;
         
  /*           !Calculate MTrans%sigma_8^2 = int dk/k win**2 T_k**2 P(k), where win is the FT of a spherical top hat */
  /*           !of radius sigr8 h^{-1} Mpc */
          
  /*   double H=h0/100.; */

  int halfngrid = ngrid/2;
  for (in = 0; in < nmu; in++)
    {

      lnko=0;
      dsig8o=0;
      sig8o=0;
      for (ik = 1; ik < halfngrid; ik++) 
	{
	  k = 2.*M_PI/BoxSize * ik; /* (ik+1); */
	  delta = k*k; /* *MTrans%TransferData(transfer_power_var,ik,itf); */
	  /*                !if (CP%NonLinear/=NonLinear_None) delta= delta* MTrans%NonLinearScaling(ik,itf) */
	  /*                !sigma_8 defined "as though it were linear" */

	  x= k *sigr8;
	  win =3*(sin(x)-x*cos(x))/pow(x,3.);
	  /* 	  lnk=log(k); */
	  lnk = k;
	  if (ik==1) 
	    dlnk=0.5;
	  /*                  !Approx for 2._dl/(CP%InitPower%an(in)+3)  [From int_0^k_1 dk/k k^4 P(k)] */
	  /*                  !Contribution should be very small in any case  */
	  else
	    dlnk=lnk-lnko;
	  if (ipower[ik*nmu+in]!=0)
	    powers = dpower[ik*nmu+in]/(double)ipower[ik*nmu+in];
	  else 
	    powers = 0.; 
	  dsig8=pow(win*delta,2.)*powers;
	  sig8[in]+=(dsig8+dsig8o)*dlnk/2.;
	  dsig8o=dsig8;
	  lnko=lnk;
	  
	}
      /*       sig8[in]*=1./(2.*M_PI*M_PI); */
      fprintf(stderr, "sigma8 = %d %g \n", in, sig8[in]);

    }

  return(0);
}

int equal_angle(int k, int i, int j, int ks, int nbin)
{
  int theta;
  float theta_tmp;

  theta_tmp = acos((float)k/(float)ks)*180/M_PI; /*angle betwen 0 and 90*/

  /*   fprintf(stderr, "%d %d %d %d %g ", i,j,k, ks, theta_tmp, theta); */

  theta_tmp /= (float)(nbin);

  theta = rint(theta_tmp);

  /*   fprintf(stderr, "%d \n", theta);  */

  return theta;
}

int angular_projection(int k, int ks, int nmu, float *theta)
{

  /*Projects spherical wave number ks along the line of sight*/

  int mu, n;
  float mu_tmp;

  int nk = 10; 
  float ks_max =  log10(0.1); // log10(0.1)
  float ks_min = log10(0.01); // log10(0.01)
  float ks_step = (ks_max-ks_min)/(float)(nk-1);
  float k_step = 

    mu_tmp = (float)(nmu-1)*(float)k/(float)ks;
  /*     mu_tmp = (float)(nmu-1)*(float)k/pow(10.,(log10(ks_min)+ks*ks_step)); */
  mu = rint(mu_tmp);

  /*     for (n = 0; n < nmu; n++) */
  /* 	if (mu_tmp >= theta[n] ) */
  /* 	{ */
  /* 	    mu = (nmu-1)-n; */
  /* 	    break; */
  /* 	} */

  /*     if (mu_tmp == 1.) */
  /* 	mu = 1; */
  /*     else  */
  /* 	mu = 0; */

  return mu;
}

int measure_window(int ngrid, int TotNumPart, int local_ny, int local_y_start, struct io_header_1 *header1, fftw_complex *dens_out, double *kpower, int *ipower)
{
  int i,j,k, l,ks, *kmod;
  float ks_min, ks_max, dk;
  double kx, ky, kz, kR, kmodulus, fcorrect,scorrect, window;

  /* kmax is given by the Nyquist frequency: k < k_n = 2pi/L n/2 */
  /* kmin is limited by the size of the box: k > 2pi/L */
  /* and the spacing of the grid needs to be L/n */

  ks_max = M_PI/(*header1).BoxSize * ngrid;
  ks_min = 2.*M_PI/(*header1).BoxSize;
  dk = 2*M_PI/(*header1).BoxSize;   

  float number_dens = TotNumPart/powf((*header1).BoxSize,3);
  float fftw_norm = powf(ngrid/(*header1).BoxSize,3)/powf(ngrid,3.);
  /*The data is now transposed*/

  int set_index = 20; 

  int ks_set = (int) floor( sqrt( (float)pow(set_index,2) +
				  (float)pow(set_index,2) +
				  (float)pow(set_index,2)));

  
  for (j = 0; j < local_ny; ++j)
    for (i = 0; i < ngrid; ++i)
      for (k = 0; k < (ngrid/2+1); ++k)
	{ 
	  /* Assign each shell an index, ks, the spherical wave number*/

	  ks = (int) floor( sqrt( (float)pow(i,2) +
				  (float)pow(j+local_y_start,2) +
				  (float)pow(k,2)));


	  c_re(dens_out[fftw_index_after(j,i,k,ngrid)]) = pow(c_re(dens_out[fftw_index_after(j,i,k,ngrid)]),2)+
	    pow(c_im(dens_out[fftw_index_after(j,i,k,ngrid)]),2);

	  c_im(dens_out[fftw_index_after(j,i,k,ngrid)]) = 0.0;
	  fcorrect = 1.; scorrect = 1.;

	  kx  = M_PI*(double)i/(double)ngrid; /*Missing a factor of 2?*/
	  ky  = M_PI*(double)(j+local_y_start)/(double)ngrid;
	  kz  = M_PI*(double)k/(double)ngrid;

	  if (kx > 0) fcorrect *= sin(kx)/kx;
	  if (ky > 0) fcorrect *= sin(ky)/ky;
	  if (kz > 0) fcorrect *= sin(kz)/kz;

#ifdef TSC

	  fcorrect = 1./pow(fcorrect,6.);  /*Check that the index is correct*/

	  scorrect *= (1. - pow(sin(kx),2.) + 2./15.*pow(sin(kx),4.));
	  scorrect *= (1. - pow(sin(ky),2.) + 2./15.*pow(sin(ky),4.));
	  scorrect *= (1. - pow(sin(kz),2.) + 2./15.*pow(sin(kz),4.));
#else 
	  fcorrect = 1./pow(fcorrect,2.);

#endif
	  /* 		kmodulus = sqrt(kx*kx+ky*ky+kz*kz); */
	  /* 		kR = kmodulus*(*header1).BoxSize/2.; */
	  /* 		window = 8.*M_PI/pow(kmodulus,3.) * (sin(kR) + (kR)*cos(kR)); */
	  /* 		fcorrect /= pow(window,2.); */

	  if(ks == ks_set)  /*A factor of dk has been taken out */
	    {
	      /* 		    dpower[ks] += fcorrect*(double)c_re(dens_out[fftw_index_after(j,i,k,ngrid)])*fftw_norm; - fcorrect*scorrect/number_dens;  */
	      ipower[ks] += 1;
	      kpower[index(i,j,k,ngrid)] = fcorrect*(double)c_re(dens_out[fftw_index_after(j,i,k,ngrid)])*fftw_norm; - fcorrect*scorrect/number_dens;
	    }
	}  


    
  return(0);

}


int read_header(FILE *fd, int files, int *NumPart, struct io_header_1 *header1)
{
  char   buf[200], fname[200];
  int    i,j,k,dummy;
  int    t,n,off,pc_sph;
  struct io_header_1 dummyheader;


  fread(&dummy, sizeof(dummy), 1, fd);
  fread(header1, sizeof(dummyheader), 1, fd);
  fread(&dummy, sizeof(dummy), 1, fd);

  /*  if(files==1)
      {*/
  for(k=0, (*NumPart)=0, ntot_withmasses=0; k<5; k++)
    (*NumPart)+= (*header1).npartTotal[k];
  Ngas= (*header1).npart[0];
  /*    }
	else
	{
	for(k=0, (*NumPart)=0, ntot_withmasses=0; k<5; k++)
	(*NumPart)+= (*header1).npartTotal[k];
	Ngas= (*header1).npartTotal[0];
	}

	for(k=0, ntot_withmasses=0; k<5; k++)
	{
	if((*header1).mass[k]==0)
	ntot_withmasses+= (*header1).npart[k];
	}
  */

  return (0);
}

int load_hacc(int nfiles, int fileno, char input_fname[], int *NumPart, double *TotNumPart, float BoxSize, struct particle_data **P, int byteswapped)
{
    char buf[200];
    FILE *fd; 
    if (nfiles == 1)
	sprintf(buf,"%s",input_fname);
    else
	sprintf(buf,"%s.%d",input_fname,fileno);

/* 	fprintf(stderr, "Now reading %s\n", buf);  */
    if(!(fd=fopen(buf,"r")))
    {
	printf("can't open file `%s`\n",buf);
	exit(0);
    }


    /*Need to allocate particle data P*/
    /*Guess a size for P*/
    int maxpart = 2.5e7; 


    if(!(*P = malloc(maxpart*sizeof(struct particle_data))))
    {
	fprintf(stderr, "Couldn't allocate %d for particle data\n", maxpart);
	MPI_Abort(MPI_COMM_WORLD, -1); 
    }

    
/*     if (!(*P = malloc(2048.*2048.*2048./500.*sizeof(struct particle_data)))) */
/*     { */
/* 	sprintf(stderr, "Couldn't allocate memory for particles \n");  */
/* 	return(1); */
/*     }  */

    int nread = 1; 
    float any[7]; 

    *NumPart = 0; 

    while (*NumPart < maxpart)
    {
	if (byteswapped)
	{
	    int ID; 
	    nread = ReadByteSwappedFloat(fd,&(any[0])); 
	    if (!nread) break;
	    nread = ReadByteSwappedFloat(fd,&(any[1])); 
	    nread = ReadByteSwappedFloat(fd,&(any[2])); 
	    nread = ReadByteSwappedFloat(fd,&(any[3])); 
	    nread = ReadByteSwappedFloat(fd,&(any[4])); 
	    nread = ReadByteSwappedFloat(fd,&(any[5])); 
	    nread = ReadByteSwappedFloat(fd,&(any[6])); 
	    nread = ReadByteSwappedInt(fd,&ID); // On the Blue Gene systems, sizeof(long int) = 4

	}
	else
	{
	    long int ID; 
	    nread = fread(&(any[0]), sizeof(float), 7, fd); 
	    if (!nread) break;
	    nread = fread(&ID, sizeof(long int), 1, fd); 
/* 	    if (*NumPart < 10) */
/* 		fprintf(stderr, "%f %f %f %f %f %f %ld\n", any[0], any[2], any[4], any[1], any[3], any[5], ID);  */
	}


	(*P)[*NumPart].Pos[0] = any[0];
	(*P)[*NumPart].Pos[1] = any[2];
	(*P)[*NumPart].Pos[2] = any[4];
	(*P)[*NumPart].Vel[0] = any[1];
	(*P)[*NumPart].Vel[1] = any[3];
	(*P)[*NumPart].Vel[2] = any[5];

	(*NumPart)++; 
    }
    *TotNumPart += (double)*NumPart;

    return(0); 
}

int ReadByteSwappedFloat(FILE *fptr,float *n)
{
  unsigned char *cptr,tmp;

  if(sizeof(float) != 4)
    {
      fprintf(stderr,"ReadFloat: sizeof(float)=%ld and not 4\n",sizeof(float));
      exit(0);
    }

  if (fread(n,4,1,fptr) != 1)
    return(0);

  cptr = (unsigned char *)n;
  tmp     = cptr[0];
  cptr[0] = cptr[3];
  cptr[3] = tmp;
  tmp     = cptr[1];
  cptr[1] = cptr[2];
  cptr[2] = tmp;

  return(1);
}

int ReadByteSwappedInt(FILE *fptr, int *n)
{
  unsigned char *cptr,tmp;
// A long int on BGP is actually 4 bytes
  if (fread(n,4,1,fptr) != 1)
      return(0);

  cptr = (unsigned char *)n;
  tmp     = cptr[0];
  cptr[0] = cptr[3];
  cptr[3] = tmp;
  tmp     = cptr[1];
  cptr[1] = cptr[2];
  cptr[2] = tmp;

  return(1);
}


int load_snapshot(int nfiles, int fileno, char input_fname[], int *NumPart, int *TotNumPart, struct particle_data **P, struct io_header_1 *header1)
{
  FILE *fd;
  char   buf[200];
  int    i,j,k,dummy, pc, pc_new;
  int    t,n,off;

  pc = 0; 

  if (nfiles == 1)
    sprintf(buf,"%s",input_fname);
  else
    sprintf(buf,"%s.%d",input_fname,fileno);

  
  if(!(fd=fopen(buf,"r")))
    {
      printf("can't open file `%s`\n",buf);
      exit(0);
    }
	
/*       printf("reading `%s' ...\n ",buf); fflush(stdout); */


  fread(&dummy, sizeof(dummy), 1, fd);
  fread(header1, sizeof(struct io_header_1), 1, fd);
  fread(&dummy, sizeof(dummy), 1, fd);
      
  fread(&dummy, sizeof(dummy), 1, fd);

  (*NumPart) = 0; (*TotNumPart) = 0;

  for (k = 0; k < 6; k++)
    {
      (*NumPart) += (*header1).npart[k];
      (*TotNumPart) += (*header1).npartTotal[k];
    }

/*       fprintf(stderr, "NumPart in this file = %d \n", (*NumPart)); */
/*       fprintf(stderr, "TotNumPart in this file = %d \n", (*TotNumPart)); */

  if(!(*P = malloc(*NumPart * sizeof(struct particle_data))))
    {
      fprintf(stderr, "Couldn't allocate %d for particle data\n", (*NumPart));
      MPI_Abort(MPI_COMM_WORLD, -1); 
    }



  for(k=0,pc_new=pc;k<6;k++)
    {
      for(n=0;n<(*header1).npart[k];n++)
	{
	  fread(&((*P)[pc_new].Pos[0]), sizeof(float), 3, fd);
/* 	  fprintf(stderr, "%d %f %f %f\n", pc_new, (*P)[pc_new].Pos[0]/1000.,(*P)[pc_new].Pos[1]/1000.,(*P)[pc_new].Pos[2]/1000.); */
	  (*P)[pc_new].Pos[0]/=(*header1).BoxSize;
	  (*P)[pc_new].Pos[1]/=(*header1).BoxSize;
	  (*P)[pc_new].Pos[2]/=(*header1).BoxSize;

	  /* 	  (*P)[pc_new].Pos[0]/=header1.HubbleParam; */  //Change this to box size
/* 	  (*P)[pc_new].Pos[1]/=header1.HubbleParam; */
/* 	  (*P)[pc_new].Pos[2]/=header1.HubbleParam; */
	  pc_new++;
	}
    }
  SKIP(fd);
     

  SKIP(fd);
  for(k=0,pc_new=pc;k<6;k++)
    {
      for(n=0;n<(*header1).npart[k];n++)
	{
	  fread(&((*P)[pc_new].Vel[0]), sizeof(float), 3, fd);
	  (*P)[pc_new].Vel[0] *= sqrt((*header1).time);
	  (*P)[pc_new].Vel[1] *= sqrt((*header1).time);
	  (*P)[pc_new].Vel[2] *= sqrt((*header1).time);
	  pc_new++;
	}
    }
  SKIP(fd);
    
  fclose(fd);

  return(0);

}

int redshift_space(int NumPart, float BoxSize, int los, float redshift, float w, float Omega0, struct particle_data *P)
{
  int i ; /*norm is the normalisation of line of sight vector*/
  float H_z, z, x_los, y_los, z_los; /*x_los is the difference between point and the center*/
  

  /*Old: P[i].Pos[los] += P[i].Vel[los]/(100.*HubbleParam)/(1.+z);*/

  float omega_lambda0 = 1-Omega0; 
  float omega_m0 = Omega0; 

  for (i = 0; i < NumPart; i++)
    {

      z =redshift; 

      H_z = sqrt(omega_lambda0*pow((1.+z),3.*(1.+w))+omega_m0*pow(1+z,3));

      /*Full spherical distortions, hopefully*/

      /* 	x_los = P[i].Pos[0] - cen; */
      /* 	y_los = P[i].Pos[1] - cen; */
      /* 	z_los = P[i].Pos[2] - cen;  */

      /* 	double norm = get_norm(x_los, y_los, z_los); */

      /* 	x_los /= norm;  */
      /* 	y_los /= norm;  */
      /* 	z_los /= norm;  */

      /* 	double dotproduct = P[i].Vel[0] * x_los + P[i].Vel[1] * y_los + P[i].Vel[2] * z_los; */

      /* 	P[i].Pos[0] += (1.+z)*dotproduct*x_los/H_z/100.; */
      /* 	P[i].Pos[1] += (1.+z)*dotproduct*y_los/H_z/100.; */
      /* 	P[i].Pos[2] += (1.+z)*dotproduct*z_los/H_z/100.; */

      /*Use plain parallel approx. and line of sight to be across boxes*/
      P[i].Pos[los] += (1.+z)*P[i].Vel[los]/H_z/100/BoxSize;

    }


  return (0);
}



double get_angle(float x1, float y1, float z1, float x2, float y2, float z2)
{
  double dotproduct, angle, mag1, mag2;

  dotproduct = x1 * x2 + y1 * y2 + z1 * z2;

  mag1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);

  mag2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);

  angle = dotproduct/mag1/mag2;  /* actually returns cos of angular separation */

  return angle;

}

double get_norm(float x, float y, float z)
{
  double norm; 

  norm = x*x + y*y + z*z;
  norm = sqrt(norm); 

  return norm;

}

/* double get_dotproduct(float x1, float y1, float z1, float x2, float y2, float z2) */
/* { */
/*     double dotproduct;  */

/*     dotproduct = x1 * x2 + y1 * y2 + z1 * z2; */

/*     return dotproduct; */

/* } */

int my_modulo(int ngrid, int index)
{
  if (index > (ngrid-1)) 
    index = index-ngrid;
  if (index < 0) 
    index = index+ngrid;

  
  return index;
}

float wrap_around(float hx, float x0, float x1, float BoxSize) 
{
  if (hx > 0.5*BoxSize)
    if (x1 > x0)
      return x0+BoxSize-x1; 
    else 
      return x1+BoxSize-x0; 
  else
    return hx; 

}

int reallocate_memory(int n, struct particle_data **P)  /*Need to pass &P as args */
{

  struct particle_data *Ptmp;
  printf("reallocating memory...\n");


  if(!( Ptmp = realloc(*P, n*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  fprintf(stderr, "reallocation successful \n");

  *P = Ptmp;

  return (0);
}



int write_particles(char filename[], int NumPart, struct particle_data *P)
{
    FILE *fp; 


    /*BINARY*/
    fp = fopen(filename, "wb");

    fwrite(&NumPart,sizeof(int),1,fp);        /*total number of particles in file*/

    int i;

    for (i = 0; i <NumPart; i++)
	fwrite(&(P[i].Pos[0]),sizeof(float),3,fp);

    fclose(fp);

    /*FORMATTED*/
/*     fp = fopen(filename, "w"); */

/*     fprintf(fp, "%d\n", NumPart);        /\*total number of particles in file*\/ */

/*     int i;  */

/*     for (i = 0; i <NumPart; i++) */
/* 	fprintf(fp, "%f %f %f \n", (P[i].Pos[0]),(P[i].Pos[1]), (P[i].Pos[2])); */

/*     fclose(fp); */



    return(0); 
}

int recenter_particles(int NumPart, struct particle_data *P, float *BoxSize, float *NewBoxSize)
{
    int i; 

    for (i = 0; i < NumPart; i++)
    {
	P[i].Pos[0] += 0.5*(*NewBoxSize - *BoxSize); 
	P[i].Pos[1] += 0.5*(*NewBoxSize - *BoxSize); 
	P[i].Pos[2] += 0.5*(*NewBoxSize - *BoxSize); 
    }

    *BoxSize = *NewBoxSize; 


    return(0); 
}



unsigned long int get_random_seed(int verbose)
{

 unsigned int seed;
 struct timeval tv;
 FILE *devrandom;

 if ((devrandom = fopen("/dev/random","r")) == NULL) {
   gettimeofday(&tv,0);
   seed = tv.tv_sec + tv.tv_usec;
   if(verbose == 1) printf("Got seed %u from gettimeofday()\n",seed);
 } else {
   fread(&seed,sizeof(seed),1,devrandom);
   if(verbose == 1) printf("Got seed %u from /dev/random\n",seed);
   fclose(devrandom);
 }

 return(seed);

}
