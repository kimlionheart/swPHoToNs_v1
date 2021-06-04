#include "photoNs.h"
#include "toptree.h"
#include "remotes.h"
#include "domains.h"
#include "fmm.h"
#include "snapshot.h"
#include "partmesh.h"
#include "initial.h"
//#include "gptl.h"
#include <mpi.h>
//#include <pthread.h>



void make_title(){
	if (0 == PROC_RANK) {
		printf(" \n\n\n\n\n");	
		printf("                    / photoNs 2 /\n\n");
		printf("\n");	

	}
}

void make_appendix(double total_time) {
	if (0 == PROC_RANK) {
		printf(" \n\n");	
		printf("                      complete !\n\n");
		printf(" Elapsed (wall-clock) Time %.1lf[sec]\n", total_time);
		printf(" Project: %s\n", CodeProj);
		printf(" output: %s\n", OutputPath);
		if (sizeof(VALIDATION) > 0)
			printf(" VALIDATION: %s\n", VALIDATION);
		printf(" \n\n");	
	}
}

double total_time_start; 
double total_time;


void driver(double ai, double af, int snap_idx, int nstep_fix) {
	//	double total_time_start = dtime();
	//double total_time;
	int n;
	double memproc;  
	double time0, time1, DTIME_DOMAIN;
	double DTIME_TOTAL_DOMAIN ;
	double DTIME_MAXIMUM;

	double DTIME_PM;
	double time_loop;
	double maxmem;  
	double maxgb;

	int loop;
	int Nstep = nstep_fix;


	double dKick_half;
	size_t mem;
	double dloga =  (log(af) - log(ai) ) /Nstep;
	double loga_i, loga_f, dk, dkh, dd;


	double imbalance ;


	for (n=0; n<NPART; n++) {
		part[n].acc[0] = 0.0;
		part[n].acc[1] = 0.0;
		part[n].acc[2] = 0.0;

		part[n].acc_pm[0] = 0.0;
		part[n].acc_pm[1] = 0.0;
		part[n].acc_pm[2] = 0.0;
	}
	DTIME_FRACTION = 1.0;

	reset_mem();

	//  use_sample_estimate();
	// use_time_estimate();

	MPI_Barrier(MPI_COMM_WORLD);//add by lxj
		if (0==PROC_RANK) {
			printf("begin domain_initialize\n");
		}
//	GPTLstart("domain");
	domain_initialize();
	MPI_Barrier(MPI_COMM_WORLD);//add by lxj
		if (0==PROC_RANK) {
			printf("begin domain_decomposition\n");
		}
	domain_decomposition() ;
	//GPTLstop("domain");

	MPI_Barrier(MPI_COMM_WORLD);

	maxmem = TotalMemory/(1024.0*1024.0*1024);  
	MPI_Allreduce(&maxmem, &maxgb, 1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);	 

	if (0==PROC_RANK)
		printf(" max memory %lf GB\n", maxgb);




#ifndef PMONLY
	MPI_Barrier(MPI_COMM_WORLD);//add by lxj
		if (0==PROC_RANK) {
			printf("begin fmm_construct\n");
		}
	//GPTLstart("fmm_construct");
	fmm_construct() ;
	//GPTLstop("fmm_construct");
		if (0==PROC_RANK) {
			printf("begin fmm_prepare\n");
		}
	//GPTLstart("fmm_prepare");
	fmm_prepare();
	//GPTLstop("fmm_prepare");
	MPI_Barrier(MPI_COMM_WORLD);//add by lxj
		if (0==PROC_RANK) {
			printf("end fmm_prepare\n");
		}
#endif


#ifdef PMTHREAD
	pthread_t tpm;
	pthread_create(&tpm, NULL, pm_thread, NULL);

#endif

#ifndef PMONLY
	MPI_Barrier(MPI_COMM_WORLD);//add by lxj
		if (0==PROC_RANK) {
			printf("begin fmm_task\n");
		}
	//GPTLstart("fmm_task");
	fmm_task();
	//GPTLstop("fmm_task");
	MPI_Barrier(MPI_COMM_WORLD);
		if (0==PROC_RANK) {
			printf("end fmm_task\n");
		}
#endif


#ifndef PMTHREAD
	dtime_pm = dtime();
	//GPTLstart("partmesh");
	partmesh();
	//GPTLstop("partmesh");
	dtime_pm = dtime() - dtime_pm;
	MPI_Barrier(MPI_COMM_WORLD);//add by lxj
		if (0==PROC_RANK) {
			printf("end partmesh\n");
		}
#endif


#ifdef PMTHREAD
	pthread_join(tpm, NULL);
#endif

#ifndef PMONLY
	MPI_Barrier(MPI_COMM_WORLD);//add by lxj
		if (0==PROC_RANK) {
			printf("begin fmm_ext\n");
		}
	//GPTLstart("fmm_ext");
	fmm_ext();
	//GPTLstop("fmm_ext");
#endif

	MPI_Barrier(MPI_COMM_WORLD); //add by lxj
		if (0==PROC_RANK) {
			printf("\n dtime_pm = %lf\n dtime_m2l = %lf, dtime_p2p = %lf [sec]\n", dtime_pm, dtime_m2l, dtime_p2p);
			printf(" dtime_fmm/task = %lf, task = %lf [sec]\n", dtime_fmm, dtime_task);
			if (dtime_pm > dtime_task)
				printf("    WARNING !  dTime PM > dTime Task\n");

			printf(" dtime_prep = %lf, ext = %lf [sec]\n", dtime_prep, dtime_ext); 

		}


	//	fmm_solver_total() ;

	//mem = total_mem_used();

	//	printf(" total memory %lu GB, %lu MB, %lu KB  at [%d]\n", mem/1024/1024/1024, mem/1024/1024, mem/1024, PROC_RANK );

	//	printf(" max memory %lu MB\n", MaxMemory/1024/1024);

	//////////////////////// LOOP //////////////////////// 

	//	int loop;
	//	int Nstep = nstep_fix;


	//	double dloga =  (log(af) - log(ai) ) /Nstep;
	//	double loga_i, loga_f, dk, dkh, dd;

	MPI_Barrier(MPI_COMM_WORLD);

	imbalance =  0.0;

	//    use_time_estimate();
	//   use_sample_estimate();

	for (loop = 0; loop < Nstep; loop++) 
//	for (loop = 0; loop < 10; loop++) 
//	for (loop = 0; loop < 1; loop++) 
	{
		time_loop = dtime();
		loop_step ++;

		loga_i = loop * dloga + log( ai);
		loga_f = (loop+1) * dloga + log( ai);

		dk = kick_loga(loga_i, loga_f);
		dd = drift_loga(loga_i, loga_f);

		if (PROC_RANK  == 0)
			printf("\n\nLOOP         a=( %lf to %lf )        %5d\n\n",  exp(loga_i), exp(loga_f), loop_step);

		dkh = 0.5*dk*GravConst;

		dKick_half = dkh;
		pack2pack_count = 0;
		pack2pack_inleaf_count = 0;



		for (n=0; n<NPART; n++) {
			part[n].vel[0] += part[n].acc_pm[0]*dkh;
			part[n].vel[1] += part[n].acc_pm[1]*dkh;
			part[n].vel[2] += part[n].acc_pm[2]*dkh;
		}


		for (n=0; n<NPART; n++) {
			part[n].vel[0] += part[n].acc[0]*dkh;
			part[n].vel[1] += part[n].acc[1]*dkh;
			part[n].vel[2] += part[n].acc[2]*dkh;
		}

		for (n=0; n<NPART; n++) {
			part[n].pos[0] += part[n].vel[0]*dd;
			part[n].pos[1] += part[n].vel[1]*dd;
			part[n].pos[2] += part[n].vel[2]*dd;
		}


		for (n=0; n<NPART; n++) {
			while (part[n].pos[0]  < 0.0 )
				part[n].pos[0] += BOXSIZE;

			while (part[n].pos[1]  < 0.0 )
				part[n].pos[1] += BOXSIZE;

			while (part[n].pos[2]  < 0.0 )
				part[n].pos[2] += BOXSIZE;

			while (part[n].pos[0]  >= BOXSIZE )
				part[n].pos[0] -= BOXSIZE;

			while (part[n].pos[1]  >= BOXSIZE )
				part[n].pos[1] -= BOXSIZE;

			while (part[n].pos[2]  >= BOXSIZE )
				part[n].pos[2] -= BOXSIZE;
		}



		//GPTLstart("domain");
#ifndef PMONLY
		fmm_deconstruct() ;
#endif



		time0 = dtime();

		//use_time_estimate();
		domain_decomposition() ;
		//GPTLstop("domain");

		DTIME_DOMAIN = dtime() - time0;

		if (0 ==PROC_RANK)
			printf(" * DOM reconstruct     (dTime = %lf sec)\n", DTIME_DOMAIN);

		//double DTIME_PM;
		time0 = dtime();

		for (n=0; n<NPART; n++) {
			part[n].acc_pm[0] = 0.0;
			part[n].acc_pm[1] = 0.0;
			part[n].acc_pm[2] = 0.0;
		}


		for (n=0; n<NPART; n++) {
			part[n].acc[0] = 0.0;
			part[n].acc[1] = 0.0;
			part[n].acc[2] = 0.0;
		}	

#ifndef PMTHREAD
//		partmesh();
#endif

#ifndef PMONLY
		//GPTLstart("fmm_construct");
		fmm_construct() ;
		//GPTLstop("fmm_construct");
		//GPTLstart("fmm_prepare");
		fmm_prepare();
		//GPTLstop("fmm_prepare");
#endif

		time1 = dtime();

		DTIME_PM = time1 - time0;

#ifdef PMTHREAD
		pthread_create(&tpm, NULL, pm_thread, NULL);

#endif


#ifndef PMONLY
		//GPTLstart("fmm_task");
		fmm_task();
		//GPTLstop("fmm_task");
#endif

#ifndef PMTHREAD
		//GPTLstart("partmesh");
		partmesh();
		//GPTLstop("partmesh");
#endif

#ifdef PMTHREAD
		pthread_join(tpm, NULL);
#endif

#ifndef PMONLY
		//GPTLstart("fmm_ext");
		fmm_ext();
		//GPTLstop("fmm_ext");
#endif



		//		if (0 == PROC_RANK)
		//		printf(" * PM  complete        (dTime = %lf sec)\n", DTIME_PARTMESH);

		maxmem = TotalMemory/(1024.0*1024.0*1024);  

		MPI_Allreduce(&maxmem, &maxgb, 1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);	    

		//		if (0==PROC_RANK)
		//		printf(" max memory %lf GB\n", maxgb);





		for (n=0; n<NPART; n++) {
			part[n].vel[0] += part[n].acc[0]*dkh;
			part[n].vel[1] += part[n].acc[1]*dkh;
			part[n].vel[2] += part[n].acc[2]*dkh;
		}


		for (n=0; n<NPART; n++) {
			part[n].vel[0] += part[n].acc_pm[0]*dkh;
			part[n].vel[1] += part[n].acc_pm[1]*dkh;
			part[n].vel[2] += part[n].acc_pm[2]*dkh;
		}


		DTIME_TOTAL_DOMAIN = 0.0;
		DTIME_MAXIMUM= 0.0;
		//		double DTIME_TOTAL_DOMAIN = 0.0;
		//	double DTIME_MAXIMUM= 0.0;
		//       DTIME_THIS_DOMAIN = pack2pack_count;

		MPI_Allreduce( &DTIME_THIS_DOMAIN, &DTIME_TOTAL_DOMAIN, 
				1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);

		MPI_Allreduce( &DTIME_THIS_DOMAIN, &DTIME_MAXIMUM, 
				1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);

		DTIME_FRACTION = DTIME_THIS_DOMAIN*PROC_SIZE/(DTIME_TOTAL_DOMAIN+0.0001);


		//		if (0 ==PROC_RANK)
		//		printf(" * FMM complete        (dTime = %lf sec)\n", DTIME_MAXIMUM);

		imbalance =  1.0 - DTIME_TOTAL_DOMAIN /((double) PROC_SIZE * DTIME_MAXIMUM ); 

		if (0 ==PROC_RANK)
			printf("\n work-load imbalance = %1.3f %%\n", imbalance*100.0);




		total_time = dtime() - total_time_start;



		//	if ( 0 == loop%5 )   Logfile_flush();
		Logfile_flush();

		//		MPI_Barrier(MPI_COMM_WORLD);
		//        size_t mem = total_mem_used();
		//maxmem = MaxMemory/(1024.0*1024.0*1024);  

		memproc =  total_mem_used()/(1024.0*1024.0*1024);  

		MPI_Allreduce(&memproc, &maxgb, 1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);

		if (0==PROC_RANK)
			printf("\n max memory %lf GB\n", maxgb);

		if (0==PROC_RANK) {
			printf("\n dtime_pm = %lf\n dtime_m2l = %lf, dtime_p2p = %lf [sec]\n", dtime_pm, dtime_m2l, dtime_p2p);
			printf(" dtime_fmm/task = %lf, task = %lf [sec]\n", dtime_fmm, dtime_task);
			if (dtime_pm > dtime_task)
				printf("    WARNING !  dTime PM > dTime Task\n");

			printf(" dtime_prep = %lf, ext = %lf [sec]\n", dtime_prep, dtime_ext); 

		}

		time_loop = dtime() - time_loop;
		if (0==PROC_RANK)
			printf(" dtime_loop = %lf, total_time= %lf[sec]\n", time_loop, total_time ); 


		LogMessage(loop_step, 0.5*(exp(loga_f)+exp(loga_i)), time_loop, DTIME_PM, total_time, imbalance ) ;

	}
	//////////////////////// LOOP //////////////////////// 


#ifndef PMONLY
	fmm_deconstruct() ;
#endif

	domain_finalize();

	Redshift_Time = 1.0/af - 1.0;
	//GPTLstart("write");

	write_snapshot(SnapFormat, snap_idx, PROC_RANK );
	//GPTLstop("write");
}

int main(int argc, char* argv[]) 
{

	MPI_Init(&argc, &argv);
        athread_init();
	//GPTLinitialize();
//	swlu_prof_init(); 
//        swlu_prof_start();
//	swlu_debug_init();
	MPI_Comm_size(MPI_COMM_WORLD, &PROC_SIZE);
	MPI_Comm_rank(MPI_COMM_WORLD, &PROC_RANK);
	//GPTLstart("total");

	my_gb=0.0;

#ifdef PMTHREAD
	MPI_Comm_dup(MPI_COMM_WORLD, &PM_COMM_WORLD);
#endif

	MPI_Barrier(MPI_COMM_WORLD);//add by lxj
	make_title();
	MPI_Barrier(MPI_COMM_WORLD);//add by lxj
	if (0 == PROC_RANK) printf("begin initialize_nbody!\n");

	total_time_start = dtime(); 
	//GPTLstart("init");
	initialize_nbody(argv[1]);//设置预处理的值
	//GPTLstop("init");
	MPI_Barrier(MPI_COMM_WORLD);//add by lxj
	if (0 == PROC_RANK) {
		printf("end init!\n");
	}

	loop_step = 0;
	count_crushing_step = 0;

#ifdef SWATHREAD
	athread_init();
#endif 
//	 driver(0.02, 0.05, 3, 3);
	driver(0.02, 1.0, 4, 360);

	//	int n;
	//	for (n=1; n<number_snaptime; n++) {			
	//		driver(snaptime[n-1], snaptime[n], n, 100);
	//	}

#ifdef SWATHREAD
	athread_halt();
#endif 

	finalize_nbody();

	total_time = dtime() - total_time_start;
	make_appendix(total_time);
	//GPTLstop("total");
	//GPTLpr(PROC_RANK);
	//GPTLfinalize();
//swlu_prof_stop();
//swlu_prof_print();
//	GPTLpr_summary_file(MPI_COMM_WORLD, "outfile");
	MPI_Finalize();
	return 0;
}


