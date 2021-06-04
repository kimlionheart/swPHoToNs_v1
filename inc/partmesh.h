#ifndef PARTMESH_H
#define PARTMESH_H


void partition_fft();
void partition_fft2();

void partmesh_thread();
void partmesh();

void cic();

#ifdef PMTHREAD
void* pm_thread(void *arg);
#endif

void powerspectrum(char powname[]) ;

//int NP_TEST;


#define convolution convolution_
#define conv_pmonly conv_pmonly_
#define densitykspace densitykspace_


//#define N 1024
//extern void convolution(double *data, int *nside);
//extern void decomp_2d_init(int *global_size, int* vproc, int *local_start, int *local_end, int *local_size);
//


#endif
