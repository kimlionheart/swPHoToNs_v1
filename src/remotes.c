/*
 * photoNs-2
 *
 *      2018 - 8 - 12
 *	qwang@nao.cas.cn
 */	 
#include <pthread.h>
#include "photoNs.h"
#include "remotes.h"
#include <athread.h>

static int numbody;
static int numnode;
void p2p_kernel_ex(int inode, int jnode);
 extern SLAVE_FUN(p2p_kernel_ex_slave)(slave_data*);

static INT64 NTASKP2Pext;
static INT64 NTASKM2Lext;
void p2p_kernel_ex(int inode, int jnode)
{
	double dx[3], ir3, ir, dr,x2;
	int ip, jp;

	double rs = splitRadius;	
	double coeff = 2.0/sqrt(M_PI);

	int n,m, idx;
	int jstart = exrtree[jnode].son[0];
	int jend   = exrtree[jnode].son[0] + exrtree[jnode].npart;

	for (ip=leaf[inode].ipart; ip<leaf[inode].ipart+leaf[inode].npart; ip++)
	{
		for (jp=jstart; jp<jend; jp++)
		{
			dx[0] = exrbody[jp].pos[0] - part[ip].pos[0];
			dx[1] = exrbody[jp].pos[1] - part[ip].pos[1];
			dx[2] = exrbody[jp].pos[2] - part[ip].pos[2];


			dr = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2])+SoftenScale ;
//			if (dr < SoftenScale)
//				ir3 = part[jp].mass/(SoftenScale*SoftenScale*SoftenScale);
//			else
				ir3 = part[jp].mass/(dr*dr*dr);


#ifdef LONGSHORT
			double drs = 0.5*dr/rs;
			ir3 *= (erfc(drs) + coeff*drs*exp(-drs*drs));
#endif


			part[ip].acc[0] +=  dx[0] * ir3;
			part[ip].acc[1] +=  dx[1] * ir3;
			part[ip].acc[2] +=  dx[2] * ir3;

			p2p_count_remote++;

		}

	}
}
/*
void walk_m2l_remote2(int im, int jm)
{
	double dx, dy, dz, r2, dr, wi, wj;
	double dist[3];

	if ( im < first_leaf ) {
		printf(" error1 \n");
		exit(0);
	}
	if ( jm < 0 ) {
		printf(" error3 im = %d jm = %d\n", im, jm);
		exit(0);
	}

	// copy leaf
	if ( im < first_node && exrtree[jm].npart <= MAXLEAF ) 
	{

#ifndef P2POFF
		double t0 = dtime();
		p2p_kernel_ex(im, jm);			
		dtime_p2p_remote += dtime() - t0;
#endif

		return;
	}


	if ( im < first_node && exrtree[jm].npart > MAXLEAF ) 
	{
		dx = leaf[im].center[0] - exrtree[jm].center[0];
		dy = leaf[im].center[1] - exrtree[jm].center[1];
		dz = leaf[im].center[2] - exrtree[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[im].width, exrtree[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag || exrtree[jm].son[0]<0 || exrtree[jm].son[1] < 0) {
			m2l(dx, dy, dz, exrtree[jm].M, leaf[im].L);
		}
		else if (0 == flag) {
			walk_m2l_remote2(im, exrtree[jm].son[0]);
			walk_m2l_remote2(im, exrtree[jm].son[1]);
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}

		return ;

	}

	if ( im >= first_node && exrtree[jm].npart <= MAXLEAF ) {
		dx = btree[im].center[0] - exrtree[jm].center[0];
		dy = btree[im].center[1] - exrtree[jm].center[1];
		dz = btree[im].center[2] - exrtree[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;
		int flag = acceptance(btree[im].width, exrtree[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag) {
			m2l(dx, dy, dz, exrtree[jm].M, btree[im].L);
		}
		else if (0 == flag) {
			walk_m2l_remote2(btree[im].son[0], jm);
			walk_m2l_remote2(btree[im].son[1], jm);
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}


		return;
	}

	//printf(" %d %lf\n", exrtree[jm].npart, exrtree[jm].M[2]);
	if ( im >= first_node && exrtree[jm].npart > MAXLEAF ) 
	{
		dx = btree[im].center[0] - exrtree[jm].center[0];
		dy = btree[im].center[1] - exrtree[jm].center[1];
		dz = btree[im].center[2] - exrtree[jm].center[2];

		r2 = dx*dx + dy*dy + dz*dz ;

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;
		int flag = acceptance(btree[im].width, exrtree[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag ) {
			m2l(dx, dy, dz, exrtree[jm].M, btree[im].L);
		}
		else if (0 == flag) {
			if ( btree[im].width[0]+btree[im].width[1]+btree[im].width[2] 
					> exrtree[jm].width[0]+exrtree[jm].width[1]+exrtree[jm].width[2]
					|| exrtree[jm].son[0] < 0 || exrtree[jm].son[1] < 0 ) 
			{
				walk_m2l_remote2(btree[im].son[0], jm);
				walk_m2l_remote2(btree[im].son[1], jm);
			}
			else {
				walk_m2l_remote2(im, exrtree[jm].son[0]);
				walk_m2l_remote2(im, exrtree[jm].son[1]);
			}


		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return;

	}
}

*/

void walk_m2l_remote_task(int im, int jm)
{
	double dx, dy, dz, r2, dr, wi, wj;
	double dist[3];

	if ( im < first_leaf ) {
		printf(" error1 \n");
		exit(0);
	}
	if ( jm < 0 ) {
		printf(" error3 im = %d jm = %d\n", im, jm);
		exit(0);
	}

	// copy leaf
	if ( im < first_node && exrtree[jm].npart <= MAXLEAF ) 
	{

#ifndef P2POFF
		taskP2Pexts[idxP2Pext] = im;
		taskP2Pextt[idxP2Pext] = jm;
		idxP2Pext++ ;
		if(idxP2Pext>=NTASKP2Pext){
			printf("lxj-idxP2Pext=%d,NTASKP2Pext=%d\n",idxP2Pext,NTASKP2Pext);
			exit(0);
		}

		//		p2p_kernel_ex(im, jm);			

#endif

		return;
	}


	if ( im < first_node && exrtree[jm].npart > MAXLEAF ) 
	{
		dx = leaf[im].center[0] - exrtree[jm].center[0];
		dy = leaf[im].center[1] - exrtree[jm].center[1];
		dz = leaf[im].center[2] - exrtree[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[im].width, exrtree[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag || exrtree[jm].son[0]<0 || exrtree[jm].son[1] < 0) {


			taskM2Lexts[idxM2Lext] = jm;
			taskM2Lextt[idxM2Lext] = im;
			idxM2Lext ++;
			if(idxM2Lext>=NTASKM2Lext){
				printf("lxj-idxM2Lext=%d,NTASKM2Lext=%d\n",idxM2Lext,NTASKM2Lext);
				exit(0);
			}


			//        	m2l(dx, dy, dz, exrtree[jm].M, leaf[im].L);
		}
		else if (0 == flag) {
			walk_m2l_remote_task(im, exrtree[jm].son[0]);
			walk_m2l_remote_task(im, exrtree[jm].son[1]);
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}

		return ;

	}

	if ( im >= first_node && exrtree[jm].npart <= MAXLEAF ) {
		dx = btree[im].center[0] - exrtree[jm].center[0];
		dy = btree[im].center[1] - exrtree[jm].center[1];
		dz = btree[im].center[2] - exrtree[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;
		int flag = acceptance(btree[im].width, exrtree[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag) {
			taskM2Lexts[idxM2Lext] = jm;
			taskM2Lextt[idxM2Lext] = im;
			idxM2Lext ++;
			if(idxM2Lext>=NTASKM2Lext){
				printf("lxj-idxM2Lext=%d,NTASKM2Lext=%d\n",idxM2Lext,NTASKM2Lext);
				exit(0);
			}


			//			m2l(dx, dy, dz, exrtree[jm].M, btree[im].L);
		}
		else if (0 == flag) {
			walk_m2l_remote_task(btree[im].son[0], jm);
			walk_m2l_remote_task(btree[im].son[1], jm);
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}


		return;
	}

	//printf(" %d %lf\n", exrtree[jm].npart, exrtree[jm].M[2]);
	if ( im >= first_node && exrtree[jm].npart > MAXLEAF ) 
	{
		dx = btree[im].center[0] - exrtree[jm].center[0];
		dy = btree[im].center[1] - exrtree[jm].center[1];
		dz = btree[im].center[2] - exrtree[jm].center[2];

		r2 = dx*dx + dy*dy + dz*dz ;

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;
		int flag = acceptance(btree[im].width, exrtree[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag ) {

			taskM2Lexts[idxM2Lext] = jm;
			taskM2Lextt[idxM2Lext] = im;
			idxM2Lext ++;
			if(idxM2Lext>=NTASKM2Lext){
				printf("lxj-idxM2Lext=%d,NTASKM2Lext=%d\n",idxM2Lext,NTASKM2Lext);
				exit(0);
			}


			//			m2l(dx, dy, dz, exrtree[jm].M, btree[im].L);
		}
		else if (0 == flag) {
			if ( btree[im].width[0]+btree[im].width[1]+btree[im].width[2] 
					> exrtree[jm].width[0]+exrtree[jm].width[1]+exrtree[jm].width[2]
					|| exrtree[jm].son[0] < 0 || exrtree[jm].son[1] < 0 ) 
			{
				walk_m2l_remote_task(btree[im].son[0], jm);
				walk_m2l_remote_task(btree[im].son[1], jm);
			}
			else {
				walk_m2l_remote_task(im, exrtree[jm].son[0]);
				walk_m2l_remote_task(im, exrtree[jm].son[1]);
			}


		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return;

	}
}



void prepare_sendtree2(int isend, int ilocal, int tNode, int D, double displace[3])
{
	int n, idx, p;
	double dx, dy, dz, dr;

	if (ilocal < first_node && ilocal >= first_leaf ) {
		int ileaf = ilocal;

		exstree[isend].npart = leaf[ileaf].npart;
		exstree[isend].center[0] = leaf[ileaf].center[0] + displace[0];
		exstree[isend].center[1] = leaf[ileaf].center[1] + displace[1];
		exstree[isend].center[2] = leaf[ileaf].center[2] + displace[2];

		exstree[isend].width[0] = leaf[ileaf].width[0];
		exstree[isend].width[1] = leaf[ileaf].width[1];
		exstree[isend].width[2] = leaf[ileaf].width[2];

		for (p=0; p<NMULTI; p++) {
			exstree[isend].M[p] = leaf[ileaf].M[p];
		}

		exstree[isend].son[0] = numbody;

		for (p=leaf[ileaf].ipart; p<leaf[ileaf].ipart+leaf[ileaf].npart; p++) {
			exsbody[numbody].pos[0] = part[p].pos[0] + displace[0];
			exsbody[numbody].pos[1] = part[p].pos[1] + displace[1];
			exsbody[numbody].pos[2] = part[p].pos[2] + displace[2];
			exsbody[numbody].mass = part[p].mass;
			numbody++;
		}

		exstree[isend].son[1] = numbody;
		numnode ++;

		return;
	}

	dx = toptree[tNode].center[0] - btree[ilocal].center[0] - displace[0];
	dy = toptree[tNode].center[1] - btree[ilocal].center[1] - displace[1];
	dz = toptree[tNode].center[2] - btree[ilocal].center[2] - displace[2];

	if (dx < 0.0)
		dx = -dx;
	if (dy < 0.0)
		dy = -dy;
	if (dz < 0.0)
		dz = -dz;

	dx -= (toptree[tNode].width[0] + btree[ilocal].width[0])*0.5;
	dy -= (toptree[tNode].width[1] + btree[ilocal].width[1])*0.5;
	dz -= (toptree[tNode].width[2] + btree[ilocal].width[2])*0.5;

	dr = 0.0;
	if (dx > 0.0)
		dr += dx*dx;
	if (dy > 0.0)
		dr += dy*dy;
	if (dz > 0.0)
		dr += dz*dz;
	dr = sqrt(dr);

	exstree[isend].npart     = btree[ilocal].npart;

	exstree[isend].center[0] = btree[ilocal].center[0] + displace[0];
	exstree[isend].center[1] = btree[ilocal].center[1] + displace[1];
	exstree[isend].center[2] = btree[ilocal].center[2] + displace[2];

	exstree[isend].width[0]  = btree[ilocal].width[0];
	exstree[isend].width[1]  = btree[ilocal].width[1];
	exstree[isend].width[2]  = btree[ilocal].width[2];

	for (p=0; p<NMULTI; p++) {
		exstree[isend].M[p] = btree[ilocal].M[p];
	}

	numnode ++;

	double width_max;
	width_max = btree[ilocal].width[0];
	if (width_max < btree[ilocal].width[1])    
		width_max = btree[ilocal].width[1];
	if (width_max < btree[ilocal].width[2])    
		width_max = btree[ilocal].width[2];


#ifdef LONGSHORT
	if (dr >= cutoffRadius) {
		exstree[isend].son[0] = -1;
		exstree[isend].son[1] = -1;
		return;
	}
#endif


	if ( width_max < 0.95*open_angle*dr ) {
		exstree[isend].son[0] = -1;
		exstree[isend].son[1] = -1;
		return;
	}
	else {
		for (n=0; n<NSON; n++) {
			idx = btree[ilocal].son[n] ;
			if (idx >= first_leaf) {
				exstree[isend].son[n] = numnode;
				prepare_sendtree2(numnode, idx, tNode, (D+1)%3, displace);
			} 
		}
	}

}


//static INT64 NTASKP2Pext;
//static INT64 NTASKM2Lext;

void fmm_remote_task (int numrecv)
{

	//	NTASKP2Pext = (INT64)(NPART*400.0);
	//	NTASKM2Lext = (INT64)((last_node - first_node)*100.0 );
	NTASKP2Pext = (INT64)(NPART*400.0);
	NTASKM2Lext = (INT64)((last_node - first_node)*400.0 );

	//	printf(" Ntask %ld %ld\n",NTASKP2P, NTASKM2L );
	idxM2Lext = 0;
	idxP2Pext = 0;

	taskM2Lexts = (int*)malloc(sizeof(int)*NTASKM2Lext);
	taskM2Lextt = (int*)malloc(sizeof(int)*NTASKM2Lext);
	taskP2Pexts = (int*)malloc(sizeof(int)*NTASKP2Pext);
	taskP2Pextt = (int*)malloc(sizeof(int)*NTASKP2Pext);

	walk_m2l_remote_task(first_node, 0);
	//	printf("\n idxM2Lext = %d [%lu],  idxP2Pext = %d [%lu]\n", idxM2Lext , NTASKM2Lext, idxP2Pext, NTASKP2Pext);
	//	walk_m2l_remote2(first_node, 0);
	int n;


//#pragma omp parallel for num_threads(NumThread)
	for (n=0; n<idxM2Lext; n++) {

		int im = taskM2Lextt[n];	
		int jm = taskM2Lexts[n];	

		double dx, dy, dz;
		//   omp_set_lock(&writelock);
		if ( im < first_node) {
			dx = leaf[im].center[0] - exrtree[jm].center[0];
			dy = leaf[im].center[1] - exrtree[jm].center[1];
			dz = leaf[im].center[2] - exrtree[jm].center[2];

			m2l(dx, dy ,dz, exrtree[jm].M, leaf[im].L);
		}
		if ( im >= first_node ) {

			dx = btree[im].center[0] - exrtree[jm].center[0];
			dy = btree[im].center[1] - exrtree[jm].center[1];
			dz = btree[im].center[2] - exrtree[jm].center[2];

			m2l(dx, dy ,dz, exrtree[jm].M, btree[im].L);
		}

		//    omp_unset_lock(&writelock);
	}

//printf("ppppp\n");
//n = task_ex();
	int slave_task_num_ex = idxP2Pext/64;
	int slave_task_remainder_ex = idxP2Pext%64;
	int p2p_ex_slave_part[64][slave_task_num_ex+1][4]; 
	int i,j,k,im,jm,inode,jnode;
	double dtime_get_ex = 0.0;
	double dtime_all_ex = 0.0;
	double dtime_comput_ex = 0.0;
	slave_data slave_all_ex;
	slave_all_ex.s_i_part = &p2p_ex_slave_part;
	slave_all_ex.s_j_part = &p2p_ex_slave_part;
	slave_all_ex.slave_remote = exrbody;
	slave_all_ex.mass = part[0].mass;
	slave_all_ex.slave_part = part;
	slave_all_ex.rs =  splitRadius;
	slave_all_ex.coeff = 2.0/sqrt(M_PI);
	slave_all_ex.SoftenScale = SoftenScale;
	slave_all_ex.idxp2p = idxP2Pext;
n=0;
	for(i=0;i < slave_task_remainder_ex;i++){
		for(j =0;j<=slave_task_num_ex;j++)
		{
			im = taskP2Pexts[n];
			jm = taskP2Pextt[n];
			p2p_ex_slave_part[i][j][0] = leaf[im].npart;
		 	p2p_ex_slave_part[i][j][1] = leaf[im].ipart;
		 	p2p_ex_slave_part[i][j][2] = exrtree[jm].son[0];
		 	p2p_ex_slave_part[i][j][3] = exrtree[jm].npart;

n++;
		}
	
	}
	for(i=slave_task_remainder_ex;i < 64;i++){
		for(j =0;j<slave_task_num_ex;j++)
		{
			im = taskP2Pexts[n];
			jm = taskP2Pextt[n];

			p2p_ex_slave_part[i][j][0] = leaf[im].npart;
		 	p2p_ex_slave_part[i][j][1] = leaf[im].ipart;
		 	p2p_ex_slave_part[i][j][2] = exrtree[jm].son[0];
		 	p2p_ex_slave_part[i][j][3] = exrtree[jm].npart;
	
n++;
		}
	
	}


athread_spawn(p2p_kernel_ex_slave,&slave_all_ex);
athread_join();
/*	for (n=0; n<idxP2Pext; n++) {
		int im = taskP2Pexts[n];	
		int jm = taskP2Pextt[n];

		p2p_kernel_ex(im, jm);

	}
*/
//dtime_p2p_ex = dtime()-dtime_p2p_ex;

	free(taskM2Lexts);
	free(taskM2Lextt);
	free(taskP2Pexts);
	free(taskP2Pextt);
}

void fmm_remote(int idx, double displace[3])
{
	int srank, d, rrank;
	int n, mi, mj, mk;


	srank = (PROC_RANK + idx ) % PROC_SIZE;
	rrank = (PROC_RANK - idx + PROC_SIZE) % PROC_SIZE;

	//  printf("[%d]s = %d %d , %lf %lf %lf\n", PROC_RANK, srank ,rrank, displace[0], displace[1], displace[2]);

	double bdl[3], bdr[3];


	for (d=0; d<3; d++) {
		bdl[d] = toptree[this_domain].center[d] - 0.5*toptree[this_domain].width[d];
		bdr[d] = toptree[this_domain].center[d] + 0.5*toptree[this_domain].width[d];
	}

	width_this_domain = toptree[this_domain].width[0];

	if ( width_this_domain < toptree[this_domain].width[1])
		width_this_domain = toptree[this_domain].width[1];

	if ( width_this_domain < toptree[this_domain].width[2])
		width_this_domain = toptree[this_domain].width[2];

	int sendnumnode, sendnumbody;
	int recvnumnode, recvnumbody;

	int tNode = srank + mostleft ;
	if (tNode > last_domain)
		tNode -= PROC_SIZE;

	numbody = 0;
	numnode = 0;

	prepare_sendtree2(numnode, first_node, tNode, direct_local_start, displace);
//printf(" %d %d\n", numbody, numnode);
//fflush(stdout);

	sendnumnode = numnode ;
	sendnumbody = numbody ;

	MPI_Request request, request2;
	MPI_Status status, status2;

	MPI_Isend(&sendnumnode, 1, MPI_INT, srank, 101, MPI_COMM_WORLD, &request);
	MPI_Recv( &recvnumnode, 1, MPI_INT, rrank, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Wait(&request, &status);

	MPI_Isend(&sendnumbody, 1, MPI_INT, srank, 102, MPI_COMM_WORLD, &request2);
	MPI_Recv( &recvnumbody, 1, MPI_INT, rrank, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Wait(&request2, &status2);


	if ( 0 == recvnumbody && 0 == recvnumnode)
		return;

	MPI_Isend(exstree, sendnumnode, strReNode, srank, 111, MPI_COMM_WORLD, &request);
	MPI_Recv( exrtree, recvnumnode, strReNode, rrank, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Wait(&request, &status);    

	MPI_Isend(exsbody, sendnumbody, strReBody, srank, 112, MPI_COMM_WORLD, &request2);
	MPI_Recv( exrbody, recvnumbody, strReBody, rrank, 112, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Wait(&request2, &status2);    

	fmm_remote_task ( recvnumnode );
	//	walk_m2l_remote2(first_node, 0);

}




