/*
 *     photoNs-2
 *
 *     2017 - 8 - 21
 *	qwang@nao.cas.cn
 */	 

//#include <omp.h>
#include <math.h>
#include "photoNs.h"
#include "fmm.h"
#include "task.h"
#include <athread.h>

static int *task_s[2];
static int *task_t[2];
static int ntask[2];

static int stask;
static int ttask;
static int idxtask;
static int LEN_TASK;
extern SLAVE_FUN(p2p_kernel_slave)(slave_data*);

void bksort_inplace(int D, Body *p, int length, int npart[2], double *split)
{
	npart[0] = npart[1] = 0;
	if (length <2 ) {
		npart[1] = length;
		return;
	}
	Body temp;
	if (length == 2 )
	{
		npart[0] = npart[1] = 1;
		*split = 0.5*(p[0].pos[D] + p[1].pos[D]);
		if (p[0].pos[D] > p[1].pos[D]) {
			temp = p[0];
			p[0] = p[1];
			p[1] = temp;
		}
		return;
	}
	int n;

	double mean = 0.0;

	for (n=0; n<length; n++) {
		mean += p[n].pos[D];
	}
	mean /= (double)length;

	int but = length - 1;

	n=0;
	while ( n < but ) {
		if ( p[n].pos[D] > mean ) {
			temp = p[n];

			//	while ( p[but].pos[D] > mean )
			while ( p[but].pos[D] > mean && but > n )
				but--;

			p[n] = p[but];
			p[but] = temp;
		}
		n++;
	}
	//	npart[0] = but + 1;
	npart[0] = but;
	npart[1] = length - npart[0];

	*split = mean;
}


void build_kdtree(int direct, int iPart, int length, int iNode)
{
	if (length == 0)
		return;

	//	printf(" iNode = %d\n", iNode);

	btree[iNode].npart = length;
	Body *p = part + iPart;

	int npart[2];
	double split;

	bksort_inplace(direct, p, length, npart, &split);

	btree[iNode].split = split;

	//printf(" direct = %d | npart = %d %d, split = %e (%e ,%e) length = %d\n", direct, npart[0], npart[1], split,min, max, length);

	int n, ip;
	ip = iPart;
	int new_direct = (direct + 1)%3;

	//printf(" direct = %d | npart = %d %d, split = %lf  length = %d\n", direct, npart[0], npart[1], split - 0.5*BOXSIZE, length);

	for (n=0; n<NSON; n++) {

		if (npart[n]  <= MAXLEAF) {
			// create a new leaf
			leaf[last_leaf].npart = npart[n];
			leaf[last_leaf].ipart = ip;

			btree[iNode].son[n] = last_leaf;

			last_leaf++;
		}
		else {
			last_node++;
			btree[iNode].son[n] = last_node;
			build_kdtree(new_direct, ip, npart[n],  last_node);
		}

		ip += npart[n];
	}

}


void center_kdtree(int direct, int iNode, double left[3], double right[3])
{
	int n, idx;
	btree[iNode].width[0] = right[0]-left[0];
	btree[iNode].width[1] = right[1]-left[1];
	btree[iNode].width[2] = right[2]-left[2];
	btree[iNode].center[0] = 0.5*(right[0]+left[0]);
	btree[iNode].center[1] = 0.5*(right[1]+left[1]);
	btree[iNode].center[2] = 0.5*(right[2]+left[2]);

	int newd = (direct + 1)%3;
	double tmp;
	for (n=0; n<NSON; n++) {
		idx = btree[iNode].son[n];

		if (idx < last_leaf ) {

			leaf[idx].width[0] = btree[iNode].width[0];
			leaf[idx].width[1] = btree[iNode].width[1];
			leaf[idx].width[2] = btree[iNode].width[2];

			leaf[idx].center[0] = btree[iNode].center[0];
			leaf[idx].center[1] = btree[iNode].center[1];
			leaf[idx].center[2] = btree[iNode].center[2];

			if (0==n){
				leaf[idx].width[direct] = btree[iNode].split-left[direct];
				leaf[idx].center[direct] = 0.5*(left[direct]+btree[iNode].split);
			}

			if (1==n) {
				leaf[idx].width[direct] = right[direct] - btree[iNode].split;
				leaf[idx].center[direct] = 0.5*(right[direct]+btree[iNode].split);
			}


		} else {

			if (0==n) {
				tmp = right[direct];
				right[direct] = btree[iNode].split;
				center_kdtree(newd, idx, left, right);
				right[direct] =tmp;
			}
			if (1==n) {
				tmp = left[direct];
				left[direct] = btree[iNode].split;
				center_kdtree(newd, idx, left, right);
				left[direct] =tmp;

			}
		}
	}

}


void build_localtree() {

	dtime_p2p = 0.0;
	dtime_p2p_remote = 0.0;
	dtime_p2p_mirror = 0.0;

	dtime_p2p_adlc = 0.0;
	dtime_p2p_adex = 0.0;

double bdleft[3];
 double bdright[3];
    int d;
int direct = direct_local_start;

	for (d=0; d<3; d++) {
		bdleft[d] = toptree[this_domain].center[d] - 0.5*toptree[this_domain].width[d];
		bdright[d] = toptree[this_domain].center[d] + 0.5*toptree[this_domain].width[d];
	}

/*
	// 0...NPART, fist_leaf...last_leaf, first_node...last_node
*/

	int NNODE = (int)(2.0*((double)NPART)/((double)MAXLEAF));
	int NLEAF = (int)(2.0*((double)NPART)/((double)MAXLEAF));

	if (NNODE > NPART)
		NNODE = NPART +1 ;
	if (NLEAF > NPART)
		NLEAF = NPART +1;

	first_leaf = last_leaf = NPART;
	first_node = last_node = NPART+NLEAF;

	// 0...NPART, fist_leaf...last_leaf, first_node...last_node

	leaf  = (Pack*)pmalloc(sizeof(Pack)*NLEAF, 30);
	btree = (Node*)pmalloc(sizeof(Node)*NNODE, 31);

	int n, m;
	for (n=0; n<NNODE; n++) {
		btree[n].npart = 0;
		btree[n].son[0] = btree[n].son[1] = -1;
		btree[n].split = 0.0;
		btree[n].width[0] = 0.0;
		btree[n].width[1] = 0.0;
		btree[n].width[2] = 0.0;

		btree[n].center[0] = 0.0;
		btree[n].center[1] = 0.0;
		btree[n].center[2] = 0.0;

		for (m=0; m<NMULTI; m++) {
			btree[n].M[m] = 0.0;
			btree[n].L[m] = 0.0;
		}
	}


	for (n=0; n<NLEAF; n++) {
		leaf[n].npart = 0;
		leaf[n].ipart = 0;

		leaf[n].width[0] =0.0;
		leaf[n].width[1] =0.0;
		leaf[n].width[2] =0.0;

		leaf[n].center[0] = 0.0;
		leaf[n].center[1] = 0.0;
		leaf[n].center[2] = 0.0;

		for (m=0; m<NMULTI; m++) {
			leaf[n].M[m] = 0.0;
			leaf[n].L[m] = 0.0;
		}

	}

	btree -= first_node;
	leaf  -= first_leaf;

	build_kdtree(direct,0, NPART,  first_node);

	center_kdtree(direct, first_node, bdleft, bdright);
}


// 0==open, 1==accept, -1==abort
inline int acceptance(double wi[3], double wj[3], double dist[3]) 
{
	double min[3];
	double w[3];
    double dm2, comp = 1.0;

	w[0] = ( wi[0] + wj[0] ) * 0.5;	
	w[1] = ( wi[1] + wj[1] ) * 0.5;	
	w[2] = ( wi[2] + wj[2] ) * 0.5;	

	if (dist[0] < 0.0) {
		min[0] = -dist[0] - w[0];
	}
	else {
		min[0] = dist[0] - w[0];
	}

	if (dist[1] < 0.0) {
		min[1] = -dist[1] - w[1];
	}
	else {
		min[1] = dist[1] - w[1];
	}

	if (dist[2] < 0.0) {
		min[2] = -dist[2] - w[2];
	}
	else {
		min[2] = dist[2] - w[2];
	}

	// neighbour
	if (min[0] + min[1] + min[2] < 0.0001)
		return 0;

        dm2 =  min[0]*min[0] + min[1]*min[1] + min[2]*min[2] ;


        if ( dm2 <= (min[0]+min[1]+min[2])*(min[0]+min[1]+min[2])*0.5 )
                comp *= 0.1;


#ifdef LONGSHORT

        if ( dm2 > cutoffRadius * cutoffRadius) {
                return -1;

        }

        if (dm2 > splitRadius*splitRadius)
                comp = 0.1;

#endif

        if (  w[0]*w[0] + w[1]*w[1] + w[2]*w[2] < comp * open_angle * open_angle * dm2 )  {
                return 1;
        }
        else
                return 0;

}


void fmm_construct( ) {

	int rank, d;
	int n, mi, mj, mk;
	double t0, t1, t2, t3, t4, t5, t6, t7, t8;

	t0 = dtime();

	construct_toptree(PROC_SIZE);

	for (n=0; n<=last_domain; n++) {
		toptree[n].split = domtree[n].split;
	}

#ifdef PERIODIC_CONDITION

	double bdl[3] = {     0.0,     0.0,     0.0 };
	double bdr[3] = { BOXSIZE, BOXSIZE, BOXSIZE };

#else

	double bdl[3] = { BoxMinimum, BoxMinimum, BoxMinimum };
	double bdr[3] = { BoxMaximum, BoxMaximum, BoxMaximum };

#endif

	int direct = 0;


	center_toptree(direct, 0, bdl, bdr);

}



////////////////////////////////////////////////////////////////

///// need stort im, jm and for m2l or p2p list. meanwhile the extremely case should be considered.

int *ts;
int *tt;
//pthread_t tid;
int p1st = 1;
void *status;
int alt;

void *task_compute_p2p(void *arg);
int argt[2];

void turn2compute_p2p(){
	if (p1st == 1) {
		p1st = 0;
	} 
	else {
//		pthread_join(tid, &status);
		athread_join();
	}
	ntask[alt] = idxtask;

	argt[0] = alt;
	argt[1] = ntask[alt];
//	pthread_create(&tid, NULL, task_compute_p2p, (void*)argt);
	task_compute_p2p((void*)argt);
	//last join (if not the first)
	//fork current
	//exchange & reset
	
//	task_compute_p2p(NULL);
// printf(" turn %d, idxtas = %d\n", alt, idxtask);
	alt = (alt+1)%2;

	ts = task_s[alt];
	tt = task_t[alt];

	idxP2P +=idxtask;
	idxtask = 0;

}

void walk_task_p2p(int im, int jm)
{
	if ( -1 == im  || -1 == jm)
		return;

	if ( im < first_leaf || jm < first_leaf ) {
		printf("error\n");
		exit(0);
	}

	if ( im == jm ) {
		if ( im <  last_leaf ) {	
			*(ts + idxtask) = jm;
			*(tt + idxtask) = im;
			idxtask ++;

			if ( idxtask == LEN_TASK ) {
				//printf(" turn \n");
				turn2compute_p2p();

			}

		}
		if ( im >= first_node ) {
			walk_task_p2p(btree[im].son[0], btree[jm].son[0]);
			walk_task_p2p(btree[im].son[0], btree[jm].son[1]);
			walk_task_p2p(btree[im].son[1], btree[jm].son[0]);
			walk_task_p2p(btree[im].son[1], btree[jm].son[1]);
		}
		return;
	}
	double dx, dy, dz, r2, dr, wi, wj;
	double dist[3];

	if ( im < last_leaf && jm < last_leaf ) {
			*(ts + idxtask) = jm;
			*(tt + idxtask) = im;
			idxtask ++;

			if ( idxtask == LEN_TASK ) {

				turn2compute_p2p();

			}
		return;
	}

	if ( im < first_node && jm >= first_node ) 
	{
		dx = leaf[im].center[0] - btree[jm].center[0];
		dy = leaf[im].center[1] - btree[jm].center[1];
		dz = leaf[im].center[2] - btree[jm].center[2];
		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[im].width, btree[jm].width, dist);

		if ( 1 == flag ) {
		}
		else if (0 == flag) {
			walk_task_p2p(im, btree[jm].son[0]);
			walk_task_p2p(im, btree[jm].son[1]);
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return ;
	}	

	if ( im >= first_node && jm < first_node ) {
		dx = btree[im].center[0] - leaf[jm].center[0];
		dy = btree[im].center[1] - leaf[jm].center[1];
		dz = btree[im].center[2] - leaf[jm].center[2];
		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[jm].width, btree[im].width, dist);

		if ( 1 == flag ) {
		}
		else if (0 == flag) {
			walk_task_p2p(btree[im].son[0], jm);
			walk_task_p2p(btree[im].son[1], jm);
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return ;
	}

	if ( im >= first_node && jm >= first_node ) {
		dx = btree[im].center[0] - btree[jm].center[0];
		dy = btree[im].center[1] - btree[jm].center[1];
		dz = btree[im].center[2] - btree[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(btree[im].width, btree[jm].width, dist);

		if ( 1 == flag ) {

		}
		else if (0 == flag) {
			if ( btree[im].width[0]+btree[im].width[1] + btree[im].width[2] 
			   > btree[jm].width[0]+btree[jm].width[1] + btree[jm].width[2] ) 
			{
				walk_task_p2p(btree[im].son[0], jm);
				walk_task_p2p(btree[im].son[1], jm);
			}
			else {
				walk_task_p2p(im, btree[jm].son[0]);
				walk_task_p2p(im, btree[jm].son[1]);
			}
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return ;
	}
}


void *task_compute_m2l(void *arg);

void turn2compute_m2l(){
	if (p1st == 1) {
		p1st = 0;
	} 
	else {
//		pthread_join(tid, &status);
	}
	ntask[alt] = idxtask;

	argt[0] = alt;
	argt[1] = ntask[alt];
//	pthread_create(&tid, NULL, task_compute_m2l, (void*)argt);
	task_compute_m2l((void*)argt);
	//last join (if not the first)
	//fork current
	//exchange & reset
	
//	task_compute_p2p(NULL);
// printf(" turn %d, idxtas = %d\n", alt, idxtask);
	alt = (alt+1)%2;

	ts = task_s[alt];
	tt = task_t[alt];

	idxM2L +=idxtask;
	idxtask = 0;

}



void walk_task_m2l(int im, int jm)
{
	if ( -1 == im  || -1 == jm)
		return;

	if ( im < first_leaf || jm < first_leaf ) {
		printf("error\n");
		exit(0);
	}

	if ( im == jm ) {

		if ( im < last_leaf ) // leaf 
		{

		}

		if ( im >= first_node ) 
		{
			walk_task_m2l(btree[im].son[0], btree[jm].son[0]);
			walk_task_m2l(btree[im].son[0], btree[jm].son[1]);
			walk_task_m2l(btree[im].son[1], btree[jm].son[0]);
			walk_task_m2l(btree[im].son[1], btree[jm].son[1]);
		}

		return ;
	}
	double dx, dy, dz, r2, dr, wi, wj;
	double dist[3];

	if ( im < last_leaf && jm < last_leaf ) {

#ifndef P2POFF
//		double t0 = dtime();

	//		taskP2Ps[idxP2P] = jm;
	//		taskP2Pt[idxP2P] = im;

	//		idxP2P ++;
//		p2p_kernel(im, jm); // outer leaf
		
//		dtime_p2p += dtime() - t0;

#endif
		return;
	}

	if ( im < first_node && jm >= first_node ) 
	{
		dx = leaf[im].center[0] - btree[jm].center[0];
		dy = leaf[im].center[1] - btree[jm].center[1];
		dz = leaf[im].center[2] - btree[jm].center[2];
		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[im].width, btree[jm].width, dist);

		if ( 1 == flag ) {
//			m2l(dx, dy, dz, btree[jm].M, leaf[im].L);
			*(ts + idxtask) = jm;
			*(tt + idxtask) = im;
			idxtask ++;
			if ( idxtask == LEN_TASK ) {
				turn2compute_m2l();
			}
	
//			taskM2Ls[idxM2L] = jm;
//			taskM2Lt[idxM2L] = im;
//			idxM2L ++;

		//	printf(" m2l : %d -> %d\n", jm, im);
		}
		else if (0 == flag) {
			walk_task_m2l(im, btree[jm].son[0]);
			walk_task_m2l(im, btree[jm].son[1]);
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}

		return ;


	}	

	if ( im >= first_node && jm < first_node ) {
		dx = btree[im].center[0] - leaf[jm].center[0];
		dy = btree[im].center[1] - leaf[jm].center[1];
		dz = btree[im].center[2] - leaf[jm].center[2];
		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[jm].width, btree[im].width, dist);

		if ( 1 == flag ) {
			*(ts + idxtask) = jm;
			*(tt + idxtask) = im;
			idxtask ++;
			if ( idxtask == LEN_TASK ) {
				turn2compute_m2l();
			}
		}
		else if (0 == flag) {
			walk_task_m2l(btree[im].son[0], jm);
			walk_task_m2l(btree[im].son[1], jm);
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return ;

	}


	if ( im >= first_node && jm >= first_node ) {
		dx = btree[im].center[0] - btree[jm].center[0];
		dy = btree[im].center[1] - btree[jm].center[1];
		dz = btree[im].center[2] - btree[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(btree[im].width, btree[jm].width, dist);

		if ( 1 == flag ) {
			*(ts + idxtask) = jm;
			*(tt + idxtask) = im;
			idxtask ++;
			if ( idxtask == LEN_TASK ) {
				turn2compute_m2l();
			}

		}
		else if (0 == flag) {
			if ( btree[im].width[0]+btree[im].width[1] + btree[im].width[2] 
			   > btree[jm].width[0]+btree[jm].width[1] + btree[jm].width[2] ) 
			{
				walk_task_m2l(btree[im].son[0], jm);
				walk_task_m2l(btree[im].son[1], jm);
			}
			else {
				walk_task_m2l(im, btree[jm].son[0]);
				walk_task_m2l(im, btree[jm].son[1]);
			}
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return;
	}
}


static double dTlocal, dTremote, dTmirror;
static INT64 NTASKP2P;
static INT64 NTASKM2L;

void fmm_prepare() 
{
	int rank, d;
	int n, mi, mj, mk;
	double t0, t1, t2, t3, t4, t5, t6, t7, t8;
	double dTlocal, dTremote, dTmirror;
   	dtime_prep = dtime();

	dtime_p2p_remote = 0.0;
	dtime_p2p_mirror = 0.0;

	t0 = dtime();
	m2l_count = 0;
	double bdl[3], bdr[3];

	for (d=0; d<3; d++) {
		bdl[d] = toptree[this_domain].center[d] - 0.5*toptree[this_domain].width[d];
		bdr[d] = toptree[this_domain].center[d] + 0.5*toptree[this_domain].width[d];
	}

    	build_localtree();

	for (n=first_leaf; n<last_leaf; n++)
		p2m(leaf[n].ipart, leaf[n].npart, leaf[n].center, leaf[n].M);

	walk_m2m(first_node);

	int idx[8];

	connect_local_toptree();


	width_this_domain = toptree[this_domain].width[0];

	if ( width_this_domain < toptree[this_domain].width[1])
		width_this_domain = toptree[this_domain].width[1];

	if ( width_this_domain  < toptree[this_domain].width[2])
		width_this_domain = toptree[this_domain].width[2];

	ExtDomain = (int*)pmalloc(sizeof(ExtDomain)*PROC_SIZE, 32);

	for (n=0; n<PROC_SIZE; n++)
		ExtDomain[n] = 0;

	dtime_prep = dtime() - dtime_prep;

	t1 = dtime();
    	dtime_fmm = t1 - t0;
	idxM2L = 0; 
	idxP2P = 0;


}

void task_prepare_p2p() {
	alt = 0;
	p1st = 1;
	idxtask = 0;

	ts = task_s[alt];
	tt = task_t[alt];

	walk_task_p2p(first_node, first_node);


}


void task_prepare_m2l() {
	alt = 0;
	p1st = 1;
	idxtask = 0;

	ts = task_s[alt];
	tt = task_t[alt];

	walk_task_m2l(first_node, first_node);


}
int p2p_slave_part[64][256][4];
slave_data slave_all;
void *task_compute_p2p(void *argt) {
	int n;
	int *par = argt;
	int c = par[0];
	int nt =  par[1];
    	double t0 = dtime();
int slave_task_num = nt/64;
int slave_task_remainder = nt%64;
//int p2p_slave_part[64][slave_task_num][2];
	int p2p_slave_part1[64][slave_task_num][2];
	int i,j,k,m,l,ip,jp,inode,jnode;
	 n = 0;

k=0;
//dtime_compute = dtime();
/*for(i=0;i < slave_task_remainder;i++){
	for(j =0;j<=slave_task_num;j++){
	
	inode = task_t[c][n];
	jnode = task_s[c][n];
			     p2p_slave_part[i][j][0] = leaf[inode].npart;
 			     p2p_slave_part[i][j][1] = leaf[inode].ipart;
 			     p2p_slave_part[i][j][2] = leaf[jnode].npart;
 			     p2p_slave_part[i][j][3] = leaf[jnode].ipart;
	n++;
	}
}*/
//for(i=slave_task_remainder;i < 64;i++){
for(i=0;i < 64;i++){
for(j =0;j<slave_task_num;j++){
	inode = task_t[c][n];
	jnode = task_s[c][n];
			     p2p_slave_part[i][j][0] = leaf[inode].npart;
 			     p2p_slave_part[i][j][1] = leaf[inode].ipart;
 			     p2p_slave_part[i][j][2] = leaf[jnode].npart;
 			     p2p_slave_part[i][j][3] = leaf[jnode].ipart;
	n++;
if(leaf[jnode].npart==0||leaf[jnode].npart==0) printf("oooo");	
	}

}
//printf("ooo %lf %lf",part[].pos[leaf[].ipart],part[].pos[1]);
slave_all.s_i_part = &p2p_slave_part1[0][0][0];
slave_all.s_j_part = &p2p_slave_part[0][0][0]; 
slave_all.mass = part[0].mass;
slave_all.slave_part = part;
slave_all.rs = splitRadius;
slave_all.coeff = 2.0/sqrt(M_PI);
slave_all.SoftenScale = SoftenScale;
slave_all.idxp2p = nt;
//printf("iiii %p  %p \n",part,&p2p_slave_part[0][0][0]);
athread_spawn(p2p_kernel_slave,&slave_all);

	for (; n<nt; n++) {
		int im = task_t[c][n];
		int jm = task_s[c][n];
  		p2p_kernel(im, jm);
	}
    	
//athread_join();
dtime_p2p += dtime() - t0;
}

void *task_compute_m2l(void *argt) {
	int n;
	int *par = argt;
	int c = par[0];
	int nt =  par[1];
    	double t0 = dtime();

	for (n=0; n<nt; n++) {
		int im = task_t[c][n];
		int jm = task_s[c][n];
		double dx, dy, dz;
		//   omp_set_lock(&writelock);
		if ( im < first_node ) {
			dx = leaf[im].center[0] - btree[jm].center[0];
			dy = leaf[im].center[1] - btree[jm].center[1];
			dz = leaf[im].center[2] - btree[jm].center[2];
			m2l(dx, dy ,dz, btree[jm].M, leaf[im].L);
		}
		else if ( jm < first_node ) {
			dx = btree[im].center[0] - leaf[jm].center[0];
			dy = btree[im].center[1] - leaf[jm].center[1];
			dz = btree[im].center[2] - leaf[jm].center[2];
			m2l(dx, dy ,dz, leaf[jm].M, btree[im].L);
		}
		else {
			dx = btree[im].center[0] - btree[jm].center[0];
			dy = btree[im].center[1] - btree[jm].center[1];
			dz = btree[im].center[2] - btree[jm].center[2];
			m2l(dx, dy ,dz, btree[jm].M, btree[im].L);
		}	
	}
    	dtime_m2l += dtime() - t0;
}

void fmm_task() {
	int n;
	double t1  = dtime();
	LEN_TASK = 16384;

	task_s[0] = (int*)malloc(sizeof(int)*LEN_TASK);
	if(task_s[0]==0){
		printf("[%d]out of memory (malloc) task_s[0]\n", PROC_RANK);
		exit(1000);
	}
	task_t[0] = (int*)malloc(sizeof(int)*LEN_TASK);
	if(task_t[0]==0){
		printf("[%d]out of memory (malloc) task_t[0]\n", PROC_RANK);
		exit(1000);
	}

	task_s[1] = (int*)malloc(sizeof(int)*LEN_TASK);
	if(task_s[1]==0){
		printf("[%d]out of memory (malloc) task_s[1]\n", PROC_RANK);
		exit(1000);
	}
	task_t[1] = (int*)malloc(sizeof(int)*LEN_TASK);
	if(task_t[1]==0){
		printf("[%d]out of memory (malloc) task_t[1]\n", PROC_RANK);
		exit(1000);
	}
//double t2  = dtime();
	//GPTLstart("p2p");
	task_prepare_p2p();
	if(PROC_RANK==0){
		printf("end task_prepare_p2p\n");
	}
//printf(" build tree = %lf\n", dtime() - t2);
	if (p1st != 1) {
		athread_join();
//		pthread_join(tid, &status);
//		printf(" final join\n");
	}
	idxP2P +=idxtask;
	if(PROC_RANK==0) printf("idxP2P=%ld\n",idxP2P);
	ntask[alt] = idxtask;
//	int arg[2];
	argt[0] = alt;
	argt[1] = ntask[alt];

//	pthread_create(&tid, NULL, task_compute_p2p, argt);
	task_compute_p2p(argt);

	athread_join();
	//GPTLstop("p2p");
//	pthread_join(tid, &status);
	if(PROC_RANK==0){
		printf("end task_compute_p2p\n");
	}

	//GPTLstart("m2l");
    	double t0 = dtime();
	task_prepare_m2l();
//printf(" build tree = %lf\n", dtime() - t2);
	if (p1st != 1) {
//		pthread_join(tid, &status);
//		printf(" final join\n");
	}
	if(PROC_RANK==0){
		printf("end task_prepare_m2l\n");
	}
	idxM2L +=idxtask;
	if(PROC_RANK==0) printf("idxM2L=%ld\n",idxM2L);

	ntask[alt] = idxtask;
//	int arg[2];
	argt[0] = alt;
	argt[1] = ntask[alt];

//	pthread_create(&tid, NULL, task_compute_m2l, argt);
	task_compute_m2l(argt);

//	pthread_join(tid, &status);
	//GPTLstop("m2l");

	if(PROC_RANK==0){
		printf("end task_compute_m2l\n");
	}

	free(task_s[0]);
	free(task_t[0]);

	free(task_s[1]);
	free(task_t[1]);
	dtime_task=dtime() - t1;
	
	if (0==PROC_RANK) printf(" fmm time = %lf,task_m2l=%lf \n", dtime() - t1,dtime() - t0);
}


void fmm_ext(){
	int rank, d;
	int n, mi, mj, mk;
	double t0, t1, t2, t3, t4, t5, t6, t7, t8;
	double dTlocal, dTremote, dTmirror;

	t0 = dtime();

	dtime_ext = dtime();

	//printf("\n idxM2L = %d [%lu],  idxP2P = %d [%lu]\n", idxM2L , NTASKM2L, idxP2P, NTASKP2P);

//	free(taskM2Ls);
//	free(taskM2Lt);
//	free(taskP2Ps);
//	free(taskP2Pt);

	double lenExTree = 2.0/MAXLEAF;
	double lenExBody = 3.5;
	//NPART_MEAN = (int)((double)NPART_TOTAL)/((double)PROC_SIZE);

	lenExTree *= ((double)NPART_TOTAL)/((double)PROC_SIZE) * sizeof(RemoteNode) ;
	lenExBody *= ((double)NPART_TOTAL)/((double)PROC_SIZE) * sizeof(RemoteBody) ;

//	printf(" %lf %lf\n", lenExTree, lenExBody);
//	fflush(stdout);

	exstree = (RemoteNode*)malloc((size_t)lenExTree);
//	exstree = (RemoteNode*)malloc(sizeof(RemoteNode)*NPART_MEAN*lenExTree);
	if (exstree == 0) {
		printf("[%d]out of memory (malloc) exstree=%d\n",PROC_RANK);
		exit(1000);
	} //add by lxj
	exsbody = (RemoteBody*)malloc((size_t)lenExBody);
//	exsbody = (RemoteBody*)malloc(sizeof(RemoteBody)*NPART_MEAN*lenExBody);
	if (exsbody == 0) {
		printf("[%d]out of memory (malloc) exsbody=%d\n",PROC_RANK);
		exit(1000);
	} //add by lxj
	exrtree = (RemoteNode*)malloc((size_t)lenExTree);
//	exrtree = (RemoteNode*)malloc(sizeof(RemoteNode)*NPART_MEAN*lenExTree);
	if (exrtree == 0) {
		printf("[%d]out of memory (malloc) exrtree=%d\n",PROC_RANK);
		exit(1000);
	} //add by lxj
	exrbody = (RemoteBody*)malloc((size_t)lenExBody);
//	exrbody = (RemoteBody*)malloc(sizeof(RemoteBody)*NPART_MEAN*lenExBody);
	if (exrbody == 0) {
		printf("[%d]out of memory (malloc) exrbody=%d\n",PROC_RANK);
		exit(1000);
	} //add by lxj

	t4 = dtime();

	MPI_Barrier(MPI_COMM_WORLD);

	double shift[3] = {0.0, 0.0, 0.0};

	t5 = dtime();
//	GPTLstart("fmm_remote");

	for (n=1; n<PROC_SIZE; n++) {
		fmm_remote(n, shift);
	}

	t6 = dtime();

#ifdef PERIODIC_CONDITION
	for (mi=-1; mi<=1; mi++) {	
		for (mj=-1; mj<=1; mj++) {	
			for (mk=-1; mk<=1; mk++) {
				if ( 0== mi && 0== mj && 0 == mk)
					continue;

				shift[0] = (double)mi*BOXSIZE;	
				shift[1] = (double)mj*BOXSIZE;	
				shift[2] = (double)mk*BOXSIZE;	

				for (n=0; n<PROC_SIZE; n++) {
					fmm_remote(n, shift);
				}
			}
		}
	}
#endif
//	GPTLstop("fmm_remote");
	t7 = dtime();
//printf("ok\n");
	MPI_Barrier(MPI_COMM_WORLD);
	free( exstree );
	free( exsbody );
	free( exrtree );
	free( exrbody );
	//GPTLstart("walk_l2l");
	walk_l2l(first_node);
	//GPTLstop("walk_l2l");
	//GPTLstart("l2p");

	for (n=first_leaf; n<last_leaf; n++)
		l2p(n);
	//GPTLstop("l2p");

	t8 = dtime();

    	dtime_fmm += (t8 - t7);

//	dTlocal = (t8-t7) + (t1-t0);
	double dTm2l = (t3-t2);
	double dTp2p = (t4-t3);
	dTremote = (t6-t5);
	dTmirror = (t7-t6);

	//	DTIME_THIS_DOMAIN = dTlocal + dTremote + dTmirror;
//	DTIME_THIS_DOMAIN = dTlocal + dTm2l + dTp2p + dTremote;

	DTIME_THIS_DOMAIN = idxP2P + idxM2L;


	dtime_fmm_remote = t7 - t5;

//	printf(" dtime_fmm_remote %lf %lf\n", dTremote, dtime_p2p_remote);


    dtime_ext = dtime() - dtime_ext;

//    if (0==PROC_RANK)
  //      printf(" FMM ext dTIME = %lf sec\n", dtime() - t0);
}

void fmm_deconstruct() 
{


	btree += first_node;
	//printf(" free node %p\n", btree);
	//	free(btree);
	pfree(btree, 31);

	leaf  += first_leaf;
	//	free(leaf);
	pfree(leaf, 30);



	pfree(ExtDomain, 32);

	deconstruct_toptree();

}

