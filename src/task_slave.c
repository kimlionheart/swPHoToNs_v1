/*
 *     photoNs-2
 *
 *     2018 - 10 - 7
 *	qwang@nao.cas.cn
 */	 
#include <stdlib.h>
#include <stdio.h>
//#include <omp.h>

//#include "photoNs.h"
#include "task.h"
#include "slave.h"


extern int     tnp_core[NCORE][NSEG];
extern int     snp_core[NCORE][NSEG];

extern double   spx_core[NCORE][NSEG][NPACK];
extern double   spy_core[NCORE][NSEG][NPACK];
extern double   spz_core[NCORE][NSEG][NPACK];

extern double   tpx_core[NCORE][NSEG][NPACK];
extern double   tpy_core[NCORE][NSEG][NPACK];
extern double   tpz_core[NCORE][NSEG][NPACK];

extern double   tpx_core[NCORE][NSEG][NPACK];
extern double   tpy_core[NCORE][NSEG][NPACK];
extern double   tpz_core[NCORE][NSEG][NPACK];

extern double   task_param[4];

int LENDATA  = NSEG * NPACK * 8 ;

#ifdef SWATHREAD


__thread_local double   par[4];

__thread_local volatile unsigned long get_reply, put_reply;
__thread_local volatile unsigned long start, end;
__thread_local int  my_id;

__thread_local int  tnp[NSEG];
__thread_local int  snp[NSEG];

__thread_local double   spx[NSEG][NPACK];
__thread_local double   spy[NSEG][NPACK];
__thread_local double   spz[NSEG][NPACK];

__thread_local double   tpx[NSEG][NPACK];
__thread_local double   tpy[NSEG][NPACK];
__thread_local double   tpz[NSEG][NPACK];

__thread_local double   accx[NPACK];
__thread_local double   accy[NPACK];
__thread_local double   accz[NPACK];

__thread_local double   dx, dy, dz, dr, ir3, gf, coeff, mg, irs, r, soft;
//#define SWMATH

void fmm_task_parallel_kernel_athread() {

	int  n, s, s1, s2,i,j;
	int  ns, nt;
//	double   accx[NPACK];
//	double   accy[NPACK];
//	double   accz[NPACK];

//	double   dx, dy, dz, dr, ir3, gf, coeff, mg, irs, r, soft;

	/////////////// transport ///////////////////

	my_id = athread_get_id(-1);

//printf(" my id = %d, ", my_id);
//fflush(NULL);

	get_reply = 0;

	athread_get(PE_MODE, &task_param[0], &par[0], 4*8, &get_reply, 0, 0, 0);
	// int  NSEG * 2, numpart for target & source
	athread_get(PE_MODE, &tnp_core[my_id][0], &tnp[0], NSEG*4, &get_reply, 0, 0, 0);
	athread_get(PE_MODE, &snp_core[my_id][0], &snp[0], NSEG*4, &get_reply, 0, 0, 0);
	// double   NSEG * 3,  source
	athread_get(PE_MODE, &spx_core[my_id][0], &spx[0], LENDATA, &get_reply,0,0,0);
	athread_get(PE_MODE, &spy_core[my_id][0], &spy[0], LENDATA, &get_reply,0,0,0);
	athread_get(PE_MODE, &spz_core[my_id][0], &spz[0], LENDATA, &get_reply,0,0,0);
	// double   NSEG * 3  target
	athread_get(PE_MODE, &tpx_core[my_id][0], &tpx[0], LENDATA, &get_reply,0,0,0);
	athread_get(PE_MODE, &tpy_core[my_id][0], &tpy[0], LENDATA, &get_reply,0,0,0);
	athread_get(PE_MODE, &tpz_core[my_id][0], &tpz[0], LENDATA, &get_reply,0,0,0);

	while(get_reply !=9);

	/////////////// transport ///////////////////

	coeff = par[0];
	irs = par[1];
	mg  = par[2];
	soft = par[3];
//	printf(" myid = %d\n", my_id);
//gf = 1.0;	

	for (s=0; s<NSEG; s++) {

		nt = tnp[s];
		ns = snp[s];

		for (i=0; i<NPACK; i++) {
			accx[i] = 0.0;
			accy[i] = 0.0;
			accz[i] = 0.0;
		}

		for (i=0; i<nt; i++) {
			for (j=0; j<ns; j++) {

				dx = tpx[s][i] - spx[s][j];	
				dy = tpy[s][i] - spy[s][j];	
				dz = tpz[s][i] - spz[s][j];	
				dr =  sqrt(dx*dx + dy*dy + dz*dz);
			//	dr = dx;
				r = dr*irs; 
				if (dr < soft)
					dr = soft;

				ir3 = mg/(dr*dr*dr);
//	dr = dx;
//	ir3 = 1000000000.0;
#ifdef SWMATH
//		gf = ( slave__sw5cg_erfc(r) + slave__sw5cg_exp(-r*r) * coeff * r ) ;
#else
			//	gf =  ( erfc( r ) +  coeff * r * exp(-r*r) ) ;
#endif

				gf =  ( erfc( r ) +  coeff * r * exp(-r*r) ) ;
		//		gf = 1.0;

				gf *=  ir3 ;
				accx[i] -= gf*dx;
				accy[i] -= gf*dy;
				accz[i] -= gf*dz;
			}
		}
		for (i=0; i<nt; i++) {
			tpx[s][i] = accx[i];
			tpy[s][i] = accy[i];
			tpz[s][i] = accz[i];
		}
	}

	/////////////// transport ///////////////////

	put_reply = 0 ;
	athread_put(PE_MODE, &tpx[0], &tpx_core[my_id][0], LENDATA, &put_reply, 0,0);
	athread_put(PE_MODE, &tpy[0], &tpy_core[my_id][0], LENDATA, &put_reply, 0,0);
	athread_put(PE_MODE, &tpz[0], &tpz_core[my_id][0], LENDATA, &put_reply, 0,0);
	while(put_reply != 3);
	/////////////// transport ///////////////////

}

#endif ////  SWATHREAD  ////
 
