#include "photoNs.h"
#include <sys/time.h>
#include <stdio.h>
#include <string.h>

static FILE *flog;

void LogMessage(int loop, double a, double time_short, double time_pm, double time_total, double imbalance) {
	if (0== PROC_RANK) {
		if (0 == flog) {
			printf(" error log files!\n");
			exit(0);
		}

		fprintf(flog, "%5d %3d  %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", loop_step, adaptive_level_maximum,  imbalance, a, time_pm, time_short, dtime_p2p, dtime_m2l, dtime_fmm, dtime_fmm_remote , time_total );
		//    fprintf(flog, "\n");
	}
}

void Logfile_flush(){
	if (0== PROC_RANK) {
		fflush(flog);
	}
}

void initializeLogfile() {
	if (0== PROC_RANK) {
		char fname[256];
//		sprintf(fname, "%s/LOG%s.TXT", OutputPath, CodeProj);
		sprintf(fname, "%s/LOG%s.TXT", OutputPath, CodeProj);

		flog = fopen(fname, "w");
		if (0 == flog) {
			printf(" error log files!\n");
			exit(0);
		}
		fprintf(flog, "### n lvl  imbalance a_t      dTpm      dTshort  dTp2p    dTm2l  dTfmm    dText     Ttot \n");

	}
}



void finalizeLogfile() {

	if (0== PROC_RANK) {
		if (0 != flog) {
			fclose(flog);
		}
	}
}

void reset_mem() {
	int n;
	for (n=0; n<MEMBLOCK; n++) {
        if (n!= 6)
		memused[n] = 0;
	}

}

size_t total_mem_used() {
	int n;
	size_t total = 0;
	for (n=0; n<MEMBLOCK; n++)
		total += memused[n];	

	TotalMemory = total;

	if (total > MaxMemory)
		MaxMemory = total;

	return total;
} 

void* pmalloc (size_t size, int idx) {
	void *alloc = malloc(size);

	if (alloc == 0) {
		printf("[%d]out of memory (malloc) idx=%d\n", PROC_RANK, idx);
		exit(1000);
	}
	memset(alloc, 0, size);

	if (idx >= MEMBLOCK) {
		printf("uncorrect idx in xmallox!\n");
		exit(0);
	}

	memused[idx] = size ;
	TotalMemory += memused[idx];

	if (TotalMemory > MaxMemory)
		MaxMemory = TotalMemory;

	return alloc;
}

void pfree(void* ptr, int idx) {

	if (idx >= MEMBLOCK) {
		printf("uncorrect idx in xmallox!\n");
		exit(0);
	}

	free(ptr);

	TotalMemory -= memused[idx];

	memused[idx] = 0 ;
}

void mem_shift(int source, int target){
	memused[target] = memused[source];
}


void* xmalloc (size_t size) {
	void *alloc = malloc(size);

	if (alloc == 0) {
		printf("out of memory (malloc)");
		exit(1000);
	}
	memset(alloc, 0, size);

	return alloc;
}


double dtime() // Return double style timestamp using gettimeofday().
{
	double tseconds;
	struct timeval mytime;

	gettimeofday(&mytime, NULL);
	tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);

	return (tseconds);

} /* dtime() */



void ic_pmfmm(){
	//	NPART = 10000;

	long seed = 135463;
	int n;
	double r, tha, phi;
	printf(" box = %lf\n", BOXSIZE );
	double	RADIUS = BOXSIZE *0.5;

	for (n=0; n<NPART; n++) 
	{
		//		r     = pow( RADIUS * RADIUS * RADIUS * ran3(&seed),0.333333);
		r     =  RADIUS * ran3(&seed);
		tha   = asin(2*ran3(&seed)-1);
		phi   = 2.0*M_PI*ran3(&seed);


		//		part[n].pos[0] = r*cos(tha)*cos(phi) + 0.5*BOXSIZE;
		//		part[n].pos[1] = r*cos(tha)*sin(phi) + 0.5*BOXSIZE;
		//		part[n].pos[2] = r*sin(tha)+ 0.5*BOXSIZE;

		//		part[n].pos[0] = 0.0;
		//		part[n].pos[1] = 0.0;
		//		part[n].pos[2] = 0.0;

		part[n].vel[0] = 0.0;
		part[n].vel[1] = 0.0;
		part[n].vel[2] = 0.0;

		part[n].acc[0] = 0.0;
		part[n].acc[1] = 0.0;
		part[n].acc[2] = 0.0;

		part[n].mass = 0.0;


	}
	part[NPART-1].pos[0] = 0.5*BOXSIZE;
	part[NPART-1].pos[1] = 0.5*BOXSIZE;
	part[NPART-1].pos[2] = 0.5*BOXSIZE;

	part[NPART-1].mass = 1.0;//BOXSIZE;
}
/*
#ifdef DIRECT_NBODY
void direct_nbody1(int NPART, int step) {
	int n, m;
	printf(" direct - - -\n");

	double rs = splitRadius;
	double coeff = 2.0/sqrt(M_PI);
	int idx;

	for (n=0; n<NPART; n++) {
		for (m=0; m<NPART; m++) {
			if ( n == m)
				continue;

			double dx[3];
			dx[0] = part[m].pos[0] - part[n].pos[0];
			dx[1] = part[m].pos[1] - part[n].pos[1];
			dx[2] = part[m].pos[2] - part[n].pos[2];
			double x2;
			double ir3;

			x2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
			idx = (int) (x2 * invGravFunc);
			//printf(" pp idx = %d(%d)\n", idx, NumGravFunc);
			if (idx < NumGravFunc && idx >=0 )
				ir3 = part[m].mass * gravfunc[idx];
			else
				ir3 = 0.0;

			part[n].acc_direct[0] += dx[0] * ir3 ;
			part[n].acc_direct[1] += dx[1] * ir3 ;
			part[n].acc_direct[2] += dx[2] * ir3 ;
			//			part[n].acc_direct[3] -= part[m].mass/dr;
		}
	}

	FILE *fd = fopen("cmp.txt","w");
	for (n=0; n<NPART; n++) {
		fprintf(fd, "%e %e %e %e %e %e %e %e %e\n", part[n].acc[0], part[n].acc[1], part[n].acc[2], part[n].acc_direct[0], part[n].acc_direct[1], part[n].acc_direct[2], part[n].acc_pm[0], part[n].acc_pm[1],part[n].acc_pm[2]);
	}
	fclose(fd);



	fd = fopen("gravity.txt","w");

	for (n=0; n<NPART; n++) {

			double dx[3];
			dx[0] = part[n].pos[0];
			dx[1] = part[n].pos[1];
			dx[2] = part[n].pos[2];
			double x2;

			x2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
	
		double fpm = part[n].acc_pm[0]*part[n].acc_pm[0]+ part[n].acc_pm[1]*part[n].acc_pm[1]+ part[n].acc_pm[2]*part[n].acc_pm[2];
		double fs = part[n].acc[0]*part[n].acc[0]+ part[n].acc[1]*part[n].acc[1]+ part[n].acc[2]*part[n].acc[2];
		double ft = part[n].acc_direct[0]*part[n].acc_direct[0]+ part[n].acc_direct[1]*part[n].acc_direct[1]+ part[n].acc_direct[2]*part[n].acc_direct[2];

		fpm = sqrt(fpm);
		fs = sqrt(fs);
		ft = sqrt(ft);

		fprintf(fd, "%d %e %e %e %e %e %e\n", n, ft, fpm, fs, part[n].pos[0], part[n].pos[1], part[n].pos[2] );
	}

	fclose(fd);


}
#endif
*/






/*
void check_gravity1() {
#ifndef DIRECT_NBODY
	printf(" !!! MARCO : DIRECT_NBODY !!! \n");
	exit(1);
#endif

	int n,  step = 100;

	if (PROC_SIZE  != 1) {
		printf("\n !!!  error in check_gravity !!!\n");
		exit(0);
	}

	direct_nbody1(NPART, step);
}


*/




#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1; i<=54; i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1; k<=4; k++)
			for (i=1; i<=55; i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
void setup(int NPART, double BOXSIZE) {
	int n;
	seed = 378412 + PROC_RANK;

	//#define OUTPUTIC
#ifndef OUTPUTIC
	for (n=0; n<NPART; n++) {
		part[n].pos[0] = ran3(&seed)*BOXSIZE;
		part[n].pos[1] = ran3(&seed)*BOXSIZE;
		part[n].pos[2] = ran3(&seed)*BOXSIZE;
		part[n].mass = 1.0/NPART;
	}
#else
	FILE *fd = fopen("ic.txt","w");
	for (n=0; n<NPART; n++) {
		part[n].pos[0] = ran3(&seed)*BOXSIZE;
		part[n].pos[1] = ran3(&seed)*BOXSIZE;
		part[n].pos[2] = ran3(&seed)*BOXSIZE;
		part[n].mass = 1.0/NPART;

		fprintf(fd, "%lf %lf %lf %lf\n", part[n].pos[0], part[n].pos[1], part[n].pos[2], part[n].mass);
	}
	fclose(fd);
#endif

}


void setup_parallel_zoomin(int NPART, double minBox, double BOXSIZE) {
	int n;
	seed = 378412;

	int start = PROC_RANK * NPART;
	int total = 0;

	int th = (int)( 0.9*(double)(PROC_SIZE)*(double)NPART);
	//	double dummy;
	for (n=0; n<start; n++) {
		//		dummy = ran3(&seed) ;
		//		dummy = ran3(&seed) ;
		//		dummy = ran3(&seed) ;
		ran3(&seed) ;
		ran3(&seed) ;
		ran3(&seed) ;
		total ++;
	}

	//#define OUTPUTIC
#ifndef OUTPUTIC
	for (n=0; n<NPART; n++) {
		if (total > th) {
			part[n].pos[0] = ran3(&seed)*BOXSIZE;
			part[n].pos[1] = ran3(&seed)*BOXSIZE;
			part[n].pos[2] = ran3(&seed)*BOXSIZE;
			part[n].mass = 1.0/((double)(PROC_SIZE)*(double)NPART);
		}
		else {

			part[n].pos[0] = ran3(&seed)*minBox+0.5*BOXSIZE;
			part[n].pos[1] = ran3(&seed)*minBox+0.5*BOXSIZE;
			part[n].pos[2] = ran3(&seed)*minBox+0.5*BOXSIZE;
			part[n].mass = 0.00001/((double)(PROC_SIZE)*(double)NPART);
		}
		total ++;
	}
#else
	FILE *fd = fopen("ic.txt","w");
	for (n=0; n<NPART; n++) {
		part[n].pos[0] = ran3(&seed)*BOXSIZE;
		part[n].pos[1] = ran3(&seed)*BOXSIZE;
		part[n].pos[2] = ran3(&seed)*BOXSIZE;
		part[n].mass = 1.0/((double)(PROC_SIZE)*(double)NPART);

		fprintf(fd, "%lf %lf %lf %lf\n", part[n].pos[0], part[n].pos[1], part[n].pos[2], part[n].mass);
	}
	fclose(fd);
#endif
}

void setup_parallel(int NPART, double BOXSIZE) {
	int n;
	seed = 378412;

	int start = PROC_RANK * NPART;
	int total = 0;
	//	double dummy;

	for (n=0; n<start; n++) {
		//		dummy = ran3(&seed) ;
		//		dummy = ran3(&seed) ;
		//		dummy = ran3(&seed) ;

		ran3(&seed);
		ran3(&seed);
		ran3(&seed);
		ran3(&seed);
		ran3(&seed);
		ran3(&seed);

		total ++;
	}

	printf(" NPART = %d\n", NPART);
	//#define OUTPUTIC

	double pmass=100.0/((double)(PROC_SIZE)*(double)NPART)/GravConst;

	for (n=0; n<NPART; n++) {
		part[n].pos[0] = ran3(&seed)*BOXSIZE;
		part[n].pos[1] = ran3(&seed)*BOXSIZE;
		part[n].pos[2] = ran3(&seed)*BOXSIZE;

		part[n].vel[0] = 10.0* (ran3(&seed) - 0.5);
		part[n].vel[1] = 10.0* (ran3(&seed) - 0.5);
		part[n].vel[2] = 10.0* (ran3(&seed) - 0.5);

		part[n].acc[0] = 0.0;
		part[n].acc[1] = 0.0;
		part[n].acc[2] = 0.0;

		part[n].mass = pmass;
		total ++;
	}
	printf(" IC : total = %d, BOXSIZE = %lf mass = %e\n", total, BOXSIZE, pmass);

}

/*
   void setup_parallel_inplace(int NPART, double BOXSIZE) 
   {
   int n;
   seed = 378412;

   double pmass=100.0/((double)(PROC_SIZE)*(double)NPART)/GravConst;

   for (n=0; n<NPART; n++) {
   part[n].pos[0] = ran3(&seed)*BOXSIZE;
   part[n].pos[1] = ran3(&seed)*BOXSIZE;
   part[n].pos[2] = ran3(&seed)*BOXSIZE;

   part[n].vel[0] = 10.0* (ran3(&seed) - 0.5);
   part[n].vel[1] = 10.0* (ran3(&seed) - 0.5);
   part[n].vel[2] = 10.0* (ran3(&seed) - 0.5);

   part[n].acc[0] = 0.0;
   part[n].acc[1] = 0.0;
   part[n].acc[2] = 0.0;

   part[n].mass = pmass;
   total ++;
   }
   printf(" IC : total = %d, BOXSIZE = %lf mass = %e\n", total, BOXSIZE, pmass);

   }
   */

void setup_parallel_sphere(int NPART, double RADIUS) {
	int n;
	seed = 378412;

	//BOXSIZE = 2000000.0;


	int start = PROC_RANK * NPART;
	int total = 0;

	for (n=0; n<start; n++) {
		ran3(&seed) ;
		ran3(&seed) ;
		ran3(&seed) ;
		ran3(&seed) ;
		total ++;
	}
	printf(" NPART = %d\n", NPART);

	double r, tha, phi;

	for (n=0; n<NPART; n++) 
	{
		r     = pow( RADIUS * RADIUS * RADIUS * ran3(&seed),0.333333);
		tha   = asin(2*ran3(&seed)-1);
		phi   = 2.0*M_PI*ran3(&seed);


		part[n].pos[0] = r*cos(tha)*cos(phi) + 0.5*BOXSIZE;
		part[n].pos[1] = r*cos(tha)*sin(phi) + 0.5*BOXSIZE;
		part[n].pos[2] = r*sin(tha)+ 0.5*BOXSIZE;

		part[n].vel[0] = 0.0;
		part[n].vel[1] = 0.0;
		part[n].vel[2] = 0.0;

		part[n].acc[0] = 0.0;
		part[n].acc[1] = 0.0;
		part[n].acc[2] = 0.0;

		part[n].mass = 1.0/((double)(PROC_SIZE)*(double)NPART)/GravConst;
		total ++;

	}

	printf(" IC : total = %d, BOXSIZE = %lf \n", total, BOXSIZE);

}

