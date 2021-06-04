/*
 *     photoNs-2
 *
 *     2018 - 8 - 4
 *	qwang@nao.cas.cn
 */	 

#include "photoNs.h"
#include "adaptive.h"
#include "initial.h"
#include "remotes.h"


void active_particle (double ai, double af)
{
	active_count = 0;
	ErrTolIntAccuracy = 0.025;
	int n;
	int act, max;
	double ti, tf, dt;
	double a2, acc, tc;

	ti = t_flat_lcdm_a (ai);
	tf = t_flat_lcdm_a (af);
	dt = tf - ti;

	max = 0;

	for (n = 0; n < NPART; n++)
	{
		acc = part[n].acc[0] + part[n].acc_pm[0];
		a2 = acc * acc;

		acc = part[n].acc[1] + part[n].acc_pm[1];
		a2 += acc * acc;

		acc = part[n].acc[2] + part[n].acc_pm[2];
		a2 += acc * acc;

		acc = GravConst * sqrt (a2);

		tc = sqrt (2.0 * ErrTolIntAccuracy * SoftenScale / acc);

		act = 0;
		if (dt > tc)
		{
			do
			{
				tc *= 2;
				act++;
			}
			while (dt > tc);
			active_count++;
			//      printf(" active level = %d\n", act);
		}

		part[n].active = act;

		if (act > max)
			max = act;
	}
	adaptive_level_maximum = max;
	//adaptive_level_maximum = 0; 		
	//for (n=0; n<NPART; n++) 
	//			part[n].active = 0;

	//	printf (" [%d] active_count = %ld, max = %d\n", PROC_RANK, active_count,
	//			adaptive_level_maximum);
}


int update_leaf(int active_level, int iPart, int npart) {
	int n, idx;
	double dx[3];
	double m0, m1, m2, m3, m4, m5;

	int updated = 0;

	for (idx=iPart; idx<npart+iPart;  idx++) {
		if (active_level == part[idx].active) {
			updated = active_level;
			break;
		}
	}
	return updated;

}


void p2p_kernel_remote_adaptive(int inode, int jnode, int level)
{
	double dx[3], ir3, ir, dr,x2;
	int ip, jp;
pack2pack_count++;
	int n,m, idx;
	int jstart = recvtree[jnode].son[0];
	int jend   = recvtree[jnode].son[0] + recvtree[jnode].npart;

	for (ip=leaf[inode].ipart; ip<leaf[inode].ipart+leaf[inode].npart; ip++)
	{
		if (level != part[ip].active)
			continue;		

		for (jp=jstart; jp<jend; jp++)
		{
			dx[0] = recvbody[jp].pos[0] - part[ip].pos[0];
			dx[1] = recvbody[jp].pos[1] - part[ip].pos[1];
			dx[2] = recvbody[jp].pos[2] - part[ip].pos[2];


#ifdef TABLE_GRAVITY		
			x2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
			idx = (int) (x2 * invGravFunc);
			if (idx < NumGravFunc && idx > 0)
				ir3 = recvbody[jp].mass*gravfunc[idx];
			else
				ir3 = 0.0;
#endif

			part[ip].acc[0] +=  dx[0] * ir3;
			part[ip].acc[1] +=  dx[1] * ir3;
			part[ip].acc[2] +=  dx[2] * ir3;

			p2p_count_remote++;

		}

	}
}


void update_p2p_remote(int im, int jm, int level)
{
	double dx, dy, dz, r2, dr, wi, wj;
	double dist[3];

	if ( im < first_leaf ) {
		printf(" error1 \n");
		exit(0);
	}
	if ( jm < 0 ) {
		printf(" error7 im = %d jm = %d\n", im, jm);
		exit(0);
	}

	//	printf(" m2l = %d, %d\n", im , recvtree[jm].npart);


	// copy leaf
	if ( im < first_node && recvtree[jm].npart <= MAXLEAF ) {
		double t0 = dtime();

#ifndef P2POFF
		if (level == leaf[im].updated)
			p2p_kernel_remote_adaptive(im, jm, level);			
		//	p2p_kernel_remote(im, jm);			
#endif

		dtime_p2p_adex += dtime() - t0;

		return;
	}

	//	printf(" m2l = p2p , first_node=%d\n", first_node);
	//	fflush(stdout);
	if ( im < first_node && recvtree[jm].npart > MAXLEAF ) {
		dx = leaf[im].center[0] - recvtree[jm].center[0];
		dy = leaf[im].center[1] - recvtree[jm].center[1];
		dz = leaf[im].center[2] - recvtree[jm].center[2];


		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[im].width, recvtree[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag || recvtree[jm].son[0]<0 || recvtree[jm].son[1] < 0) {
			//			m2l(dx, dy, dz, recvtree[jm].M, leaf[im].L);
		}
		else if (0 == flag) {
			update_p2p_remote(im, recvtree[jm].son[0],level);
			update_p2p_remote(im, recvtree[jm].son[1],level);
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}


		return ;

	}

	//	printf(" m2l = m2p \n");
	//	fflush(stdout);
	if ( im >= first_node && recvtree[jm].npart <= MAXLEAF ) {
		dx = btree[im].center[0] - recvtree[jm].center[0];
		dy = btree[im].center[1] - recvtree[jm].center[1];
		dz = btree[im].center[2] - recvtree[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;
		int flag = acceptance(btree[im].width, recvtree[jm].width, dist);
		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag) {
			//			m2l(dx, dy, dz, recvtree[jm].M, btree[im].L);
		}
		else if (0 == flag) {
			update_p2p_remote(btree[im].son[0], jm, level);
			update_p2p_remote(btree[im].son[1], jm, level);
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}


		return;
	}

	//printf(" %d %lf\n", recvtree[jm].npart, recvtree[jm].M[2]);
	if ( im >= first_node && recvtree[jm].npart > MAXLEAF ) 
	{
		dx = btree[im].center[0] - recvtree[jm].center[0];
		dy = btree[im].center[1] - recvtree[jm].center[1];
		dz = btree[im].center[2] - recvtree[jm].center[2];

		r2 = dx*dx + dy*dy + dz*dz ;

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;
		int flag = acceptance(btree[im].width, recvtree[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag ) {
			//			m2l(dx, dy, dz, recvtree[jm].M, btree[im].L);
		}
		else if (0 == flag) {
			if (btree[im].width[0]+btree[im].width[1] + btree[im].width[2] 
					>  recvtree[jm].width[0]+recvtree[jm].width[1] + recvtree[jm].width[2]
					|| recvtree[jm].son[0] < 0 || recvtree[jm].son[1] < 0) {
				update_p2p_remote(btree[im].son[0], jm, level);
				update_p2p_remote(btree[im].son[1], jm, level);
			}
			else {
				update_p2p_remote(im, recvtree[jm].son[0], level);
				update_p2p_remote(im, recvtree[jm].son[1], level);
			}


		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}



		return;


	}
}



void p2p_kernel_mirror_adaptive(int inode, int jnode, int level)
{
	double dx[3], ir3, ir, dr,x2;
	int ip, jp;

	int n,m, idx;
	int jstart = recvtreem[jnode].son[0];
	int jend   = recvtreem[jnode].son[0] + recvtreem[jnode].npart;

	for (ip=leaf[inode].ipart; ip<leaf[inode].ipart+leaf[inode].npart; ip++)
	{
		if (part[ip].active != level)
			continue;

		for (jp=jstart; jp<jend; jp++)
		{
			dx[0] = recvbodym[jp].pos[0] - part[ip].pos[0];
			dx[1] = recvbodym[jp].pos[1] - part[ip].pos[1];
			dx[2] = recvbodym[jp].pos[2] - part[ip].pos[2];


#ifdef TABLE_GRAVITY		
			x2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
			idx = (int) (x2 * invGravFunc);
			if (idx < NumGravFunc && idx > 0)
				ir3 = recvbodym[jp].mass*gravfunc[idx];
			else
				ir3 = 0.0;
#endif

			part[ip].acc[0] +=  dx[0] * ir3;
			part[ip].acc[1] +=  dx[1] * ir3;
			part[ip].acc[2] +=  dx[2] * ir3;

			p2p_count_remote++;

		}

	}
}



void update_p2p_mirror(int im, int jm, int level)
{
	double dx, dy, dz, r2, dr, wi, wj;
	double dist[3];

	if ( im < first_leaf ) {
		printf(" error1 \n");
		exit(0);
	}
	if ( jm < 0 ) {
		printf(" error8 im = %d jm = %d\n", im, jm);
		exit(0);
	}

	//	printf(" m2l = %d, %d\n", im , recvtree[jm].npart);


	// copy leaf
	if ( im < first_node && recvtreem[jm].npart <= MAXLEAF ) {
		double t0 = dtime();

#ifndef P2POFF
		if (level == leaf[im].updated)
			p2p_kernel_mirror_adaptive(im, jm, level);			
#endif

		dtime_p2p_adex += dtime() - t0;

		return;
	}

	//	printf(" m2l = p2p , first_node=%d\n", first_node);
	//	fflush(stdout);
	if ( im < first_node && recvtreem[jm].npart > MAXLEAF ) {
		dx = leaf[im].center[0] - recvtreem[jm].center[0];
		dy = leaf[im].center[1] - recvtreem[jm].center[1];
		dz = leaf[im].center[2] - recvtreem[jm].center[2];


		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[im].width, recvtreem[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag || recvtreem[jm].son[0]<0 || recvtreem[jm].son[1] < 0) {
			//			m2l(dx, dy, dz, recvtreem[jm].M, leaf[im].L);
		}
		else if (0 == flag) {
			update_p2p_mirror(im, recvtreem[jm].son[0],level);
			update_p2p_mirror(im, recvtreem[jm].son[1],level);
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}


		return ;

	}

	//	printf(" m2l = m2p \n");
	//	fflush(stdout);
	if ( im >= first_node && recvtreem[jm].npart <= MAXLEAF ) {
		dx = btree[im].center[0] - recvtreem[jm].center[0];
		dy = btree[im].center[1] - recvtreem[jm].center[1];
		dz = btree[im].center[2] - recvtreem[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;
		int flag = acceptance(btree[im].width, recvtreem[jm].width, dist);
		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag) {
			//	m2l(dx, dy, dz, recvtree[jm].M, btree[im].L);
		}
		else if (0 == flag) {
			update_p2p_mirror(btree[im].son[0], jm, level);
			update_p2p_mirror(btree[im].son[1], jm, level);
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}


		return;
	}

	//printf(" %d %lf\n", recvtree[jm].npart, recvtree[jm].M[2]);
	if ( im >= first_node && recvtreem[jm].npart > MAXLEAF ) 
	{
		dx = btree[im].center[0] - recvtreem[jm].center[0];
		dy = btree[im].center[1] - recvtreem[jm].center[1];
		dz = btree[im].center[2] - recvtreem[jm].center[2];

		r2 = dx*dx + dy*dy + dz*dz ;

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;
		int flag = acceptance(btree[im].width, recvtreem[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag ) {
			//	m2l(dx, dy, dz, recvtree[jm].M, btree[im].L);
		}
		else if (0 == flag) {
			if (btree[im].width[0]+btree[im].width[1] + btree[im].width[2] 
					>  recvtreem[jm].width[0]+recvtreem[jm].width[1] + recvtreem[jm].width[2]
					|| recvtreem[jm].son[0] < 0 || recvtreem[jm].son[1] < 0) {
				update_p2p_mirror(btree[im].son[0], jm, level);
				update_p2p_mirror(btree[im].son[1], jm, level);
			}
			else {
				update_p2p_mirror(im, recvtreem[jm].son[0], level);
				update_p2p_mirror(im, recvtreem[jm].son[1], level);
			}
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return;
	}
}



//void update_p2p_remote(int im, int jm, int adaptive_level);

void fmm_adaptive_solver( int adaptive_level ){
	int n, idx;

	double t0, t1, t2, t3, t4, t5, t6, t7, t8;
	double dTlocal, dTremote, dTmirror;

	t0 = dtime();

	for (n=first_leaf; n<last_leaf; n++)
		leaf[n].updated = update_leaf(adaptive_level, leaf[n].ipart, leaf[n].npart);


	for (n=first_leaf; n<last_leaf; n++) {
		if ( leaf[n].updated == adaptive_level )
			for (idx=leaf[n].ipart; idx<leaf[n].npart+leaf[n].ipart; idx++){

				part[idx].acc[0] = 0.0;		
				part[idx].acc[1] = 0.0;		
				part[idx].acc[2] = 0.0;		

				l2p(n);
			}
	}

	t1 = dtime();
	update_p2p_local(first_node, first_node, adaptive_level);



	t2 = dtime();
	int rank;
	int flag_npart_remotenode;

	recvtree = bk_recvtree;
	recvbody = bk_recvbody;

	for (rank=0; rank<PROC_SIZE; rank++) {
		if (recvcnt_node[rank] > 0 ) {
			if (rank != PROC_RANK)
				update_p2p_remote(first_node, 0, adaptive_level);
			//	walk_m2l_remote(first_node, 0);

			recvtree += recvcnt_node[rank];
			recvbody += recvcnt_body[rank];
		}
	}

	for (rank=0; rank<PROC_SIZE; rank++) {
		if (recvcnt_node[rank] > 0 ) {
			recvtree -= recvcnt_node[rank];
			recvbody -= recvcnt_body[rank];
		}
	}


	t3 = dtime();
	//////////////////////////////////////////////
	
	t4 = dtime();

	recvtreem = bk_recvtreem;
	recvbodym = bk_recvbodym;

	for (rank=0; rank<PROC_SIZE; rank++) {
		if (recvcnt_nodem[rank] > 0 ) {

			flag_npart_remotenode = recvtreem[0].npart;
			for (n=0; n<recvcnt_nodem[rank]; n++) {
				if ( recvtreem[n].npart   ==  flag_npart_remotenode) {
					update_p2p_mirror(first_node, 0, adaptive_level);
				}
			}
			recvtreem += recvcnt_nodem[rank];
			recvbodym += recvcnt_bodym[rank];
		}
	}

	for (rank=0; rank<PROC_SIZE; rank++) {
		if (recvcnt_nodem[rank] > 0 ) {
			recvtreem -= recvcnt_nodem[rank];
			recvbodym -= recvcnt_bodym[rank];
		}
	}
	t5 = dtime();

	dTlocal =  (t2-t1);
	dTremote = (t5-t2);

	dtime_fmm_adlc = dTlocal  - dtime_p2p_adlc;
	dtime_fmm_adex = dTremote - dtime_p2p_adex;

}

void kdk2_level (double dkh, double dd, int level)
{
	if (level > adaptive_level_maximum)
		return;
	int n;

	kdk2_level(0.5*dkh, 0.5*dd, level+1);
	for (n = 0; n < NPART; n++)
	{
		if (part[n].active == level)
		{
			part[n].vel[0] += (part[n].acc_pm[0] + part[n].acc[0]) * dkh;
			part[n].vel[1] += (part[n].acc_pm[1] + part[n].acc[1]) * dkh;
			part[n].vel[2] += (part[n].acc_pm[2] + part[n].acc[2]) * dkh;	
		}
	}
	//printf(" cnt = %d, NPART = %d\n", cnt, NPART);
}


void kdk_level (double dkh, double dd, int level)
{
	if (level > adaptive_level_maximum)
		return;
	int n;

	for (n = 0; n < NPART; n++) {
		if (part[n].active == level) {
			part[n].vel[0] += (part[n].acc_pm[0] + part[n].acc[0]) * dkh;
			part[n].vel[1] += (part[n].acc_pm[1] + part[n].acc[1]) * dkh;
			part[n].vel[2] += (part[n].acc_pm[2] + part[n].acc[2]) * dkh;
		}
	}

	for (n = 0; n < NPART; n++) {
		if (part[n].active == level)	{
			part[n].pos[0] += part[n].vel[0] * dd * 0.5;
			part[n].pos[1] += part[n].vel[1] * dd * 0.5;
			part[n].pos[2] += part[n].vel[2] * dd * 0.5;
		}
	}


	kdk_level(0.5*dkh, 0.5*dd, level+1);

	for (n = 0; n < NPART; n++) {
		if (part[n].active == level)	{
			part[n].pos[0] += part[n].vel[0] * dd * 0.5;
			part[n].pos[1] += part[n].vel[1] * dd * 0.5;
			part[n].pos[2] += part[n].vel[2] * dd * 0.5;
		}
	}

	fmm_adaptive_solver( level );

	if (level == adaptive_level_maximum) {
		count_crushing_step ++;
		if ( 0 == PROC_RANK)
			printf("\n STEP : %5d  - LEVEL %d\n",count_crushing_step, adaptive_level_maximum);
	}

	for (n = 0; n < NPART; n++) {
		if (part[n].active == level) {
			part[n].vel[0] += (part[n].acc_pm[0] + part[n].acc[0]) * dkh;
			part[n].vel[1] += (part[n].acc_pm[1] + part[n].acc[1]) * dkh;
			part[n].vel[2] += (part[n].acc_pm[2] + part[n].acc[2]) * dkh;
		}
	}

	for (n = 0; n < NPART; n++)
	{
		if (part[n].active == level)
		{
			part[n].vel[0] += (part[n].acc_pm[0] + part[n].acc[0]) * dkh;
			part[n].vel[1] += (part[n].acc_pm[1] + part[n].acc[1]) * dkh;
			part[n].vel[2] += (part[n].acc_pm[2] + part[n].acc[2]) * dkh;	
		}
	}

	for (n = 0; n < NPART; n++) {
		if (part[n].active == level)	{
			part[n].pos[0] += part[n].vel[0] * dd * 0.5;
			part[n].pos[1] += part[n].vel[1] * dd * 0.5;
			part[n].pos[2] += part[n].vel[2] * dd * 0.5;
		}
	}

	kdk_level(0.5*dkh, 0.5*dd, level+1);

	for (n = 0; n < NPART; n++) {
		if (part[n].active == level)	{
			part[n].pos[0] += part[n].vel[0] * dd * 0.5;
			part[n].pos[1] += part[n].vel[1] * dd * 0.5;
			part[n].pos[2] += part[n].vel[2] * dd * 0.5;
		}
	}

	fmm_adaptive_solver( level );

	if (level == adaptive_level_maximum) {
		count_crushing_step ++;
		if ( 0 == PROC_RANK)
			printf(" STEP : %5d  - LEVEL %d\n\n",count_crushing_step, adaptive_level_maximum);
	}

	for (n = 0; n < NPART; n++) {
		if (part[n].active == level) {
			part[n].vel[0] += (part[n].acc_pm[0] + part[n].acc[0]) * dkh;
			part[n].vel[1] += (part[n].acc_pm[1] + part[n].acc[1]) * dkh;
			part[n].vel[2] += (part[n].acc_pm[2] + part[n].acc[2]) * dkh;
		}
	}


	//	kdk2_level(0.5*dkh, 0.5*dd, level+1);

	//printf(" cnt = %d, NPART = %d\n", cnt, NPART);

}

/*
   void _kdk_level (double dkh, double dd, int level)
   {
   if (level > adaptive_level_maximum)
   return;
   int n;

   for (n = 0; n < NPART; n++) {
   if (part[n].active == level) {
   part[n].vel[0] += (part[n].acc_pm[0] + part[n].acc[0]) * dkh;
   part[n].vel[1] += (part[n].acc_pm[1] + part[n].acc[1]) * dkh;
   part[n].vel[2] += (part[n].acc_pm[2] + part[n].acc[2]) * dkh;
   }
   }

   if (level < adaptive_level_maximum) {
   kdk_level (0.5 * dkh, 0.5 * dd, level + 1);
   }

   for (n = 0; n < NPART; n++) {
   if (part[n].active == level)	{
   part[n].pos[0] += part[n].vel[0] * dd;
   part[n].pos[1] += part[n].vel[1] * dd;
   part[n].pos[2] += part[n].vel[2] * dd;
   }
   }

   fmm_adaptive_solver( level );

   if (level < adaptive_level_maximum) {
   kdk_level (0.5 * dkh, 0.5 * dd, level + 1);
   }

   for (n = 0; n < NPART; n++)
   {
   if (part[n].active == level)
   {
   part[n].vel[0] += (part[n].acc_pm[0] + part[n].acc[0]) * dkh;
   part[n].vel[1] += (part[n].acc_pm[1] + part[n].acc[1]) * dkh;
   part[n].vel[2] += (part[n].acc_pm[2] + part[n].acc[2]) * dkh;	
   }
   }
//printf(" cnt = %d, NPART = %d\n", cnt, NPART);
}
*/

int leaf_updated(int active_level, int iPart, int npart) 
{
	int idx;
	int updated = 0;

	for (idx=iPart; idx<npart+iPart;  idx++) {
		if (active_level == part[idx].active) {
			updated = active_level;
			break;
		}
	}

	return updated;
}

/*
   int update_p2m(int active_level, int iPart, int npart, double center[3], double Mold[], int *updated) 
   {
   int n, idx;
   double dx[3];
   double m0, m1, m2, m3, m4, m5;

   int active = 0;
 *updated = 0;
 for (idx=iPart; idx<npart+iPart;  idx++) {
 if (active_level == part[idx].active) {
 *updated = active_level;
 active = 1;
 break;
 }
 }

 if ( 0 == active ) {
 return 0;	
 }

 double cent[3];
 cent[0] = cent[1] = cent[2] = 0.0;
 for (idx=iPart; idx<npart+iPart;  idx++) {
 cent[0] += part[idx].pos[0];
 cent[1] += part[idx].pos[1];
 cent[2] += part[idx].pos[2];
 }
 cent[0] /= npart;
 cent[1] /= npart;
 cent[2] /= npart;

 double dm, dm2=0.0;
 dm = cent[0] - center[0];
 dm2 += dm*dm;

 dm = cent[1] - center[1];
 dm2 += dm*dm;

 dm = cent[2] - center[2];
 dm2 += dm*dm;

 double eps = 0.000001;

 if ( dm2 < eps )
 return 0;

 center[0] = cent[0];
 center[1] = cent[1];
 center[2] = cent[2];

 double M[NMULTI];
 for (n=0; n<NMULTI; n++) 
 M[n] = 0.0;

 for (idx=iPart; idx<npart+iPart;  idx++) 
 {
 dx[0] = part[idx].pos[0] - center[0];
 dx[1] = part[idx].pos[1] - center[1];	
 dx[2] = part[idx].pos[2] - center[2];	

 m0 = part[idx].mass;
 m1 = -m0;

 M[0] += m0;

 M[X] += m1*dx[0];
 M[Y] += m1*dx[1];
 M[Z] += m1*dx[2];

#ifdef QUADRUPOLE
m2 = m0;

M[XX] += m2*dx[0]*dx[0]/2;
M[XY] += m2*dx[0]*dx[1];
M[XZ] += m2*dx[0]*dx[2];
M[YY] += m2*dx[1]*dx[1]/2;
M[YZ] += m2*dx[1]*dx[2];
M[ZZ] += m2*dx[2]*dx[2]/2;

#endif

#ifdef  OCTUPOLE
m3 = -m0;

M[XXX] += m3 * dx[0] * dx[0] * dx[0]/6;
M[XXY] += m3 * dx[1] * dx[0] * dx[0]/2;
M[XXZ] += m3 * dx[2] * dx[0] * dx[0]/2;
M[XYY] += m3 * dx[1] * dx[1] * dx[0]/2;
M[XYZ] += m3 * dx[2] * dx[1] * dx[0];
M[XZZ] += m3 * dx[2] * dx[2] * dx[0]/2;
M[YYY] += m3 * dx[1] * dx[1] * dx[1]/6;
M[YYZ] += m3 * dx[1] * dx[2] * dx[1]/2;
M[YZZ] += m3 * dx[2] * dx[2] * dx[1]/2;
M[ZZZ] += m3 * dx[2] * dx[2] * dx[2]/6;	

#endif

#ifdef  HEXADECAPOLE
m4 = m0;

M[XXXX] += m4 * dx[0] * dx[0] * dx[0] * dx[0]/24;
M[XXXY] += m4 * dx[0] * dx[0] * dx[0] * dx[1]/6;
M[XXXZ] += m4 * dx[0] * dx[0] * dx[0] * dx[2]/6;
M[XXYY] += m4 * dx[0] * dx[0] * dx[1] * dx[1]/4;
M[XXYZ] += m4 * dx[0] * dx[0] * dx[1] * dx[2]/2;
M[XXZZ] += m4 * dx[0] * dx[0] * dx[2] * dx[2]/4;
M[XYYY] += m4 * dx[0] * dx[1] * dx[1] * dx[1]/6;
M[XYYZ] += m4 * dx[0] * dx[1] * dx[1] * dx[2]/2;
M[XYZZ] += m4 * dx[0] * dx[1] * dx[2] * dx[2]/2;
M[XZZZ] += m4 * dx[0] * dx[2] * dx[2] * dx[2]/6;
M[YYYY] += m4 * dx[1] * dx[1] * dx[1] * dx[1]/24;
M[YYYZ] += m4 * dx[1] * dx[1] * dx[1] * dx[2]/6;
M[YYZZ] += m4 * dx[1] * dx[1] * dx[2] * dx[2]/4;
M[YZZZ] += m4 * dx[1] * dx[2] * dx[2] * dx[2]/6;
M[ZZZZ] += m4 * dx[2] * dx[2] * dx[2] * dx[2]/24;
#endif
}

for (n=0; n<NMULTI; n++) 
Mold[n] = M[n];

return 1;

}


int node_updated(int iNode, int adaptive_level)
{
	int n, idx, ch[2];


	for (n=0; n<NSON; n++) {
		idx = btree[iNode].son[n];
		if (idx < 0) {
			printf(" hehe \n");
			return 0;
		} else
			if (idx < first_node) {
				ch[n] = leaf[idx].updated;
			} 
			else {
				ch[n] = node_updated(idx, adaptive_level);
			}
	}

	if ( adaptive_level == ch[0] || adaptive_level == ch[1] ) {
		btree[iNode].updated = adaptive_level;
	}
	else {
		btree[iNode].updated = 0;
	}

	return btree[iNode].updated;
}


void p2p_kernel_active(int inode, int jnode, int level) {
	double dx[3], dr, ir3, ir, mp, x2;
	int ip, jp;

	int n,m, idx;

	pack2pack_count++;

	double rs = splitRadius;

	double coeff = 2.0/sqrt(M_PI);

	leaf[inode].pack2pack_cnt++;

	for (ip=leaf[inode].ipart; ip<leaf[inode].ipart+leaf[inode].npart; ip++)
	{

		if (part[ip].active != level)
			continue;

		for (jp=leaf[jnode].ipart; jp<leaf[jnode].ipart+leaf[jnode].npart; jp++)
		{
			if (jp == ip)
				continue;

			dx[0] = part[jp].pos[0] - part[ip].pos[0];
			dx[1] = part[jp].pos[1] - part[ip].pos[1];
			dx[2] = part[jp].pos[2] - part[ip].pos[2];

#ifdef TABLE_GRAVITY	
			x2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
			idx = (int) (x2 * invGravFunc);

			if (idx < NumGravFunc && idx >=0 )
				ir3 = part[jp].mass * gravfunc[idx];
			else
				ir3 = 0.0;
#else

			dr = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
			if (dr < SoftenScale)
				ir3 = part[jp].mass/(SoftenScale*SoftenScale*SoftenScale);
			else
				ir3 = part[jp].mass/(dr*dr*dr);


#ifdef LONGSHORT
			double drs = 0.5*dr/rs;
			ir3 *= (erfc(drs) + coeff*drs*exp(-drs*drs));
#endif

#endif
			part[ip].acc[0] +=  dx[0] * ir3;
			part[ip].acc[1] +=  dx[1] * ir3;
			part[ip].acc[2] +=  dx[2] * ir3;
			//			part[ip].acc[3] -= mp * ir;

			p2p_count++;
		}
	}
}
*/

void update_p2p_local(int im, int jm, int adaptive_level)
{
	if ( -1 == im  || -1 == jm )
		return;

	if ( im < first_leaf || jm < first_leaf ) {
		printf("error\n");
		exit(0);
	}

	if ( im == jm ) {
		if ( im < last_leaf ) // leaf 
		{	

#ifndef P2POFF
			double t0 = dtime();

			if ( leaf[im].updated == adaptive_level )
				p2p_kernel(im, jm); // inner leaf
			//	p2p_kernel_inleaf(im, jm); // inner leaf

			dtime_p2p_adlc += dtime() - t0;
#endif

		}

		if ( im >= first_node ) {
			update_p2p_local(btree[im].son[0], btree[jm].son[0], adaptive_level);
			update_p2p_local(btree[im].son[0], btree[jm].son[1], adaptive_level);
			update_p2p_local(btree[im].son[1], btree[jm].son[0], adaptive_level);
			update_p2p_local(btree[im].son[1], btree[jm].son[1], adaptive_level);
		}

		return ;
	}

	double dx, dy, dz, r2, dr, wi, wj;
	double dist[3];

	if ( im < last_leaf && jm < last_leaf ) {


#ifndef P2POFF
		double t0 = dtime();

		if (leaf[im].updated == adaptive_level)
			p2p_kernel(im, jm); // outer leaf

		dtime_p2p_adlc += dtime() - t0;
#endif


		return;
	}


	if ( im < first_node && jm >= first_node ) {
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
			update_p2p_local(im, btree[jm].son[0], adaptive_level);
			update_p2p_local(im, btree[jm].son[1], adaptive_level);
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
			update_p2p_local(btree[im].son[0], jm, adaptive_level);
			update_p2p_local(btree[im].son[1], jm, adaptive_level);
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
			//			m2l(dx, dy, dz, btree[jm].M, btree[im].L);
		}
		else if (0 == flag) {
			if (btree[im].width[0]+btree[im].width[1] + btree[im].width[2] > btree[jm].width[0]+btree[jm].width[1] + btree[jm].width[2] ) {
				update_p2p_local(btree[im].son[0], jm, adaptive_level);
				update_p2p_local(btree[im].son[1], jm, adaptive_level);
			}
			else {
				update_p2p_local(im, btree[jm].son[0], adaptive_level);
				update_p2p_local(im, btree[jm].son[1], adaptive_level);
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

/*
   void update_p2p_remote(int im, int jm, int adaptive_level)
   {
   double dx, dy, dz, r2, dr, wi, wj;
   double dist[3];

   if ( im < first_leaf ) {
   printf(" error1 \n");
   exit(0);
   }
   if ( jm < 0 ) {
   printf(" error2 im = %d jm = %d\n", im, jm);
   exit(0);
   }

//	printf(" m2l = %d, %d\n", im , recvtree[jm].npart);


// copy leaf
if ( im < first_node && recvtree_remote[jm].npart <= MAXLEAF ) {
double t0 = dtime();

#ifndef P2POFF
if ( leaf[im].updated == adaptive_level )
p2p_kernel_remote_adaptive(im, jm, adaptive_level);			
#endif

dtime_p2p_remote += dtime() - t0;

return;
}

//	printf(" m2l = p2p , first_node=%d\n", first_node);
//	fflush(stdout);
if ( im < first_node && recvtree_remote[jm].npart > MAXLEAF ) {
dx = leaf[im].center[0] - recvtree_remote[jm].center[0];
dy = leaf[im].center[1] - recvtree_remote[jm].center[1];
dz = leaf[im].center[2] - recvtree_remote[jm].center[2];


dist[0] = dx;
dist[1] = dy;
dist[2] = dz;

int flag = acceptance(leaf[im].width, recvtree_remote[jm].width, dist);

if (-1 == flag) {
return;
}
else if ( 1 == flag || recvtree_remote[jm].son[0]<0 || recvtree_remote[jm].son[1] < 0) {
//			m2l(dx, dy, dz, recvtree[jm].M, leaf[im].L);
}
else if (0 == flag) {
update_p2p_remote(im, recvtree_remote[jm].son[0], adaptive_level);
update_p2p_remote(im, recvtree_remote[jm].son[1], adaptive_level);
}
else {
printf(" error acceptance \n");
exit(0);
}


return ;


}

//	printf(" m2l = m2p \n");
//	fflush(stdout);
if ( im >= first_node && recvtree_remote[jm].npart <= MAXLEAF ) {
dx = btree[im].center[0] - recvtree_remote[jm].center[0];
dy = btree[im].center[1] - recvtree_remote[jm].center[1];
dz = btree[im].center[2] - recvtree_remote[jm].center[2];

dist[0] = dx;
dist[1] = dy;
dist[2] = dz;
int flag = acceptance(btree[im].width, recvtree_remote[jm].width, dist);
if (-1 == flag) {
	return;
}
else if ( 1 == flag) {
	//			m2l(dx, dy, dz, recvtree_remote[jm].M, btree[im].L);
}
else if (0 == flag) {
	update_p2p_remote(btree[im].son[0], jm, adaptive_level);
	update_p2p_remote(btree[im].son[1], jm, adaptive_level);
}
else {
	printf(" error acceptance \n");
	exit(0);
}


return;
}

//printf(" %d %lf\n", recvtree_remote[jm].npart, recvtree[jm].M[2]);
if ( im >= first_node && recvtree_remote[jm].npart > MAXLEAF ) 
{
	dx = btree[im].center[0] - recvtree_remote[jm].center[0];
	dy = btree[im].center[1] - recvtree_remote[jm].center[1];
	dz = btree[im].center[2] - recvtree_remote[jm].center[2];

	r2 = dx*dx + dy*dy + dz*dz ;
	//		dr = sqrt(r2) ;

	dist[0] = dx;
	dist[1] = dy;
	dist[2] = dz;
	int flag = acceptance(btree[im].width, recvtree_remote[jm].width, dist);

	if (-1 == flag) {
		return;
	}
	else if ( 1 == flag ) {
		//			m2l(dx, dy, dz, recvtree_remote[jm].M, btree[im].L);
	}
	else if (0 == flag) {
		if (btree[im].width[0]+btree[im].width[1] + btree[im].width[2] 
				>  recvtree_remote[jm].width[0] + recvtree_remote[jm].width[1] + recvtree_remote[jm].width[2]
				|| recvtree[jm].son[0] < 0 || recvtree[jm].son[1] < 0) 
		{
			update_p2p_remote(btree[im].son[0], jm, adaptive_level);
			update_p2p_remote(btree[im].son[1], jm, adaptive_level);
		}
		else {
			update_p2p_remote(im, recvtree_remote[jm].son[0], adaptive_level);
			update_p2p_remote(im, recvtree_remote[jm].son[1], adaptive_level);
		}
	}
	else {
		printf(" error acceptance \n");
		exit(0);
	}

	return;
}
}


void p2p_kernel_periodic(int inode, int jnode,int level)
{
	double dx[3], ir3, ir, dr,x2;
	int ip, jp;

	int n,m, idx;
	int jstart = recvtree[jnode].son[0];
	int jend   = recvtree[jnode].son[0] + recvtree[jnode].npart;

	//	double rs = splitRadius ;
	//	double coeff = 2.0/sqrt(M_PI);

	for (ip=leaf[inode].ipart; ip<leaf[inode].ipart+leaf[inode].npart; ip++)
	{
		if (part[ip].active != level)
			continue;
		for (jp=jstart; jp<jend; jp++)
		{
			dx[0] = recvbody[jp].pos[0] - part[ip].pos[0];
			dx[1] = recvbody[jp].pos[1] - part[ip].pos[1];
			dx[2] = recvbody[jp].pos[2] - part[ip].pos[2];


#ifdef TABLE_GRAVITY		
			x2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
			idx = (int) (x2 * invGravFunc);
			if (idx < NumGravFunc && idx > 0)
				ir3 = recvbody[jp].mass*gravfunc[idx];
			else
				ir3 = 0.0;
#else


			dr = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
			if (dr < SoftenScale)
				ir3 = recvbody[jp].mass/(SoftenScale * SoftenScale * SoftenScale);
			else
				ir3 = recvbody[jp].mass/(dr*dr*dr);

#ifdef LONGSHORT
			double drs = 0.5*dr/rs;
			ir3 *= (erfc(drs) + coeff*drs*exp(-drs*drs));
#endif

#endif

			part[ip].acc[0] +=  dx[0] * ir3;
			part[ip].acc[1] +=  dx[1] * ir3;
			part[ip].acc[2] +=  dx[2] * ir3;

			p2p_count_remote++;

		}
	}
}


void update_p2p_periodic(int im, int jm, int adaptive_level)
{
	double dx, dy, dz, r2, dr, wi, wj;
	double dist[3];

	if ( im < first_leaf ) {
		printf(" error1 \n");
		exit(0);
	}
	if ( jm < 0 ) {
		printf(" error6 im = %d jm = %d\n", im, jm);
		exit(0);
	}

	// copy leaf
	if ( im < first_node && recvtree[jm].npart <= MAXLEAF ) {
		double t0 = dtime();

#ifndef P2POFF
		if ( leaf[im].updated == adaptive_level )
			p2p_kernel_periodic(im, jm, adaptive_level);			
#endif

		dtime_p2p_remote += dtime() - t0;

		return;
	}

	if ( im < first_node && recvtree[jm].npart > MAXLEAF ) {
		dx = leaf[im].center[0] - recvtree[jm].center[0];
		dy = leaf[im].center[1] - recvtree[jm].center[1];
		dz = leaf[im].center[2] - recvtree[jm].center[2];


		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[im].width, recvtree[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag || recvtree[jm].son[0]<0 || recvtree[jm].son[1] < 0) {
			//			m2l(dx, dy, dz, recvtree[jm].M, leaf[im].L);
		}
		else if (0 == flag) {
			update_p2p_periodic(im, recvtree[jm].son[0], adaptive_level);
			update_p2p_periodic(im, recvtree[jm].son[1], adaptive_level);
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}


		return ;


	}

	if ( im >= first_node && recvtree[jm].npart <= MAXLEAF ) {
		dx = btree[im].center[0] - recvtree[jm].center[0];
		dy = btree[im].center[1] - recvtree[jm].center[1];
		dz = btree[im].center[2] - recvtree[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;
		int flag = acceptance(btree[im].width, recvtree[jm].width, dist);
		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag) {
			//			m2l(dx, dy, dz, recvtree[jm].M, btree[im].L);
		}
		else if (0 == flag) {
			update_p2p_periodic(btree[im].son[0], jm, adaptive_level);
			update_p2p_periodic(btree[im].son[1], jm, adaptive_level);
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}


		return;
	}

	//printf(" %d %lf\n", recvtree[jm].npart, recvtree[jm].M[2]);
	if ( im >= first_node && recvtree[jm].npart > MAXLEAF ) 
	{
		dx = btree[im].center[0] - recvtree[jm].center[0];
		dy = btree[im].center[1] - recvtree[jm].center[1];
		dz = btree[im].center[2] - recvtree[jm].center[2];

		r2 = dx*dx + dy*dy + dz*dz ;
		//		dr = sqrt(r2) ;

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;
		int flag = acceptance(btree[im].width, recvtree[jm].width, dist);

		if (-1 == flag) {
			return;
		}
		else if ( 1 == flag ) {
			//			m2l(dx, dy, dz, recvtree[jm].M, btree[im].L);
		}
		else if (0 == flag) {
			if (btree[im].width[0]+btree[im].width[1] + btree[im].width[2] 
					>  recvtree[jm].width[0]+recvtree[jm].width[1] + recvtree[jm].width[2]
					|| recvtree[jm].son[0] < 0 || recvtree[jm].son[1] < 0) {
				update_p2p_periodic(btree[im].son[0], jm, adaptive_level);
				update_p2p_periodic(btree[im].son[1], jm, adaptive_level);
			}
			else {
				update_p2p_periodic(im, recvtree[jm].son[0], adaptive_level);
				update_p2p_periodic(im, recvtree[jm].son[1], adaptive_level);
			}
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}

		return;

	}
}




void update_local( int adaptive_level )
{
	int n, idx;
	int up;
	int rank;
	int flag_npart_remotenode ;

	for (n=first_leaf; n<last_leaf; n++) 
		leaf[n].updated = 0;

	int active_count =0;
	for (n=first_leaf; n<last_leaf; n++) {
		leaf[n].updated = leaf_updated(adaptive_level, leaf[n].ipart, leaf[n].npart);
		if ( adaptive_level ==  leaf[n].updated )
			active_count ++;
	}

	node_updated(first_node, adaptive_level);

	for (n=first_leaf; n<last_leaf; n++) {
		if ( leaf[n].updated == adaptive_level )
			for (idx=leaf[n].ipart; idx<leaf[n].npart+leaf[n].ipart;  idx++) {
				part[idx].acc[0] = 0.0;		
				part[idx].acc[1] = 0.0;		
				part[idx].acc[2] = 0.0;		
			}
	}

	for (n=first_leaf; n<last_leaf; n++)
		if ( leaf[n].updated == adaptive_level )
			l2p(n);

	update_p2p_local(first_node, first_node, adaptive_level);

	////// remote //////

	for (rank=0; rank<PROC_SIZE; rank++) {
		if (recvcnt_node_remote[rank] > 0 ) {

			flag_npart_remotenode = recvtree[0].npart;
			for (n=0; n<recvcnt_node_remote[rank]; n++) {
				if ( recvtree_remote[n].npart   ==  flag_npart_remotenode) {
					update_p2p_remote(first_node, n, adaptive_level);
				}
			}
			recvtree_remote += recvcnt_node_remote[rank];
			recvbody_remote += recvcnt_body_remote[rank];
		}
	}


	for (rank=0; rank<PROC_SIZE; rank++) {
		if (recvcnt_node_remote[rank] > 0 ) {
			recvtree_remote -= recvcnt_node_remote[rank];
			recvbody_remote -= recvcnt_body_remote[rank];
		}
	}

	////// remote //////


	////// periodic //////

	for (rank=0; rank<PROC_SIZE; rank++) {
		if (recvcnt_node[rank] > 0 ) {

			flag_npart_remotenode = recvtree[0].npart;
			for (n=0; n<recvcnt_node[rank]; n++) {
				if ( recvtree[n].npart   ==  flag_npart_remotenode) {
					update_p2p_periodic(first_node, n, adaptive_level);
				}
			}
			recvtree += recvcnt_node[rank];
			recvbody += recvcnt_body[rank];
		}
	}

	for (rank=0; rank<PROC_SIZE; rank++) {
		if (recvcnt_node[rank] > 0 ) {
			recvtree -= recvcnt_node[rank];
			recvbody -= recvcnt_body[rank];
		}
	}

	////// periodic //////


}

*/


