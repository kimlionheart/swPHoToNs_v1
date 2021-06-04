#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>
#include "partmesh.h"
#include "photoNs.h"

#ifdef PMTHREAD


void* pm_thread(void *arg) {
    dtime_pm = dtime();

    partmesh_thread();

    dtime_pm = dtime() - dtime_pm;

	return NULL;
}


void partmesh_thread()
{
	int nproc, rank;
	MPI_Comm_size(PM_COMM_WORLD, &nproc);	
	MPI_Comm_rank(PM_COMM_WORLD, &rank );

	double time_0, time_1, t0, t1;
	time_0 = dtime();

	int n, c ;
	double domain_min[3] = {BOXSIZE,BOXSIZE,BOXSIZE};
	double domain_max[3] = {0.0, 0.0, 0.0};
	int local_xmin[3], local_xmax[3], isize[3];

	domain_min[0] = domain_min[1] = domain_min[2] = BOXSIZE;
	domain_max[0] = domain_max[1] = domain_max[2] = 0.0;

	for (n=0; n<NPART; n++ ) {
		if (domain_min[0] > part[n].pos[0] )
			domain_min[0] = part[n].pos[0];

		if (domain_min[1] > part[n].pos[1] )
			domain_min[1] = part[n].pos[1];

		if (domain_min[2] > part[n].pos[2] )
			domain_min[2] = part[n].pos[2];


		if (domain_max[0] < part[n].pos[0] )
			domain_max[0] = part[n].pos[0];

		if (domain_max[1] < part[n].pos[1] )
			domain_max[1] = part[n].pos[1];

		if (domain_max[2] < part[n].pos[2] )
			domain_max[2] = part[n].pos[2];
	}

//printf("domain max - [%d] %lf %lf %lf\n", PROC_RANK, domain_max[0], domain_max[1], domain_max[2]);
//printf("domain min - [%d] %lf %lf %lf\n", PROC_RANK, domain_min[0], domain_min[1], domain_min[2]);

	local_xmin[0] = (int) ((double)domain_min[0]*NSIDE / BOXSIZE);
	local_xmin[1] = (int) ((double)domain_min[1]*NSIDE / BOXSIZE);
	local_xmin[2] = (int) ((double)domain_min[2]*NSIDE / BOXSIZE);

	local_xmin[0] -= 1;
	local_xmin[1] -= 1;
	local_xmin[2] -= 1;

	local_xmin[0] -= 1;
	local_xmin[1] -= 1;
	local_xmin[2] -= 1;

	local_xmin[0] -= 1;
	local_xmin[1] -= 1;
	local_xmin[2] -= 1;


	local_xmax[0] = (int) ((double)domain_max[0]*NSIDE / BOXSIZE);
	local_xmax[1] = (int) ((double)domain_max[1]*NSIDE / BOXSIZE);
	local_xmax[2] = (int) ((double)domain_max[2]*NSIDE / BOXSIZE);

	isize[0] = local_xmax[0] - local_xmin[0]  +2 +1 +2;
	isize[1] = local_xmax[1] - local_xmin[1]  +2 +1 +2;
	isize[2] = local_xmax[2] - local_xmin[2]  +2 +1 +2;

//	printf(" p2 - %d %d %d\n", isize[0], isize[1], isize[2]);
//	printf(" p2[%d] - min: %d %d %d\n",PROC_RANK, local_xmin[0], local_xmin[1], local_xmin[2]);
//	fflush(stdout);

//MPI_Barrier(PM_COMM_WORLD);

	long meshsize = isize[0] * isize[1] * isize[2];
	double* mesh = (double*)pmalloc(sizeof(double) * meshsize, 50);


	for (n=0; n<meshsize; n++) 
		mesh[n] = 0.0;

	int i,j,k,idx;
	int ii, jj, kk;
	double norm = NSIDE/BOXSIZE;
	double delta = 1.0/norm;

	double wi, wj , wk, win, wjn, wkn;
	c=0;	
	for (n=0; n<NPART; n++) {
		//		i = (int) (part[n].pos[0]*norm) - local_xmin[0] ;	
		//		j = (int) (part[n].pos[1]*norm) - local_xmin[1] ;	
		//		k = (int) (part[n].pos[2]*norm) - local_xmin[2] ;	
		//		idx = (i*isize[1] + j)*isize[2] + k;
		//		mesh[idx] += 1.0;

		i = (int) (part[n].pos[0] *norm) ;	
		j = (int) (part[n].pos[1] *norm) ;	
		k = (int) (part[n].pos[2] *norm) ;	

		wi = (part[n].pos[0] - (i+0.5)*delta)*norm;

		if (wi > 0) {
			ii = i + 1;
		}
		else {
			wi = -wi;
			ii = i - 1;
		}
		win =  1.0 - wi;

		wj = (part[n].pos[1] - (j+0.5)*delta)*norm;
		if (wj > 0) {
			jj = j + 1;
		}
		else {
			wj = -wj;
			jj = j - 1;
		}
		wjn =  1.0 - wj;

		wk = (part[n].pos[2] - (k+0.5)*delta)*norm;
		if (wk > 0) {
			kk = k + 1;
		}
		else {
			wk = -wk;
			kk = k - 1;
		}
		wkn =  1.0 - wk;

		i -= local_xmin[0] ;	
		j -= local_xmin[1] ;	
		k -= local_xmin[2] ;	

		ii -= local_xmin[0] ;	
		jj -= local_xmin[1] ;	
		kk -= local_xmin[2] ;	


		idx = (i*isize[1] + j)*isize[2] + k;
		mesh[idx] += part[n].mass*win*wjn*wkn;	

		idx = (ii*isize[1] + j)*isize[2] + k;
		mesh[idx] += part[n].mass*wi *wjn*wkn;	

		idx = (i*isize[1] + jj)*isize[2] + k;
		mesh[idx] += part[n].mass*win*wj *wkn;	

		idx = (i*isize[1] + j)*isize[2] + kk;
		mesh[idx] += part[n].mass*win*wjn*wk ;	

		idx = (ii*isize[1] + jj)*isize[2] + k;
		mesh[idx] += part[n].mass*wi *wj *wkn;	

		idx = (ii*isize[1] + j)*isize[2] + kk;
		mesh[idx] += part[n].mass*wi *wjn*wk ;	

		idx = (i*isize[1] + jj)*isize[2] + kk;
		mesh[idx] += part[n].mass*win*wj *wk ;	

		idx = (ii*isize[1] + jj)*isize[2] + kk;
		mesh[idx] += part[n].mass* wi *wj *wk ;	

		c++;
	}

	double renormal = (NSIDE/BOXSIZE);
	renormal = renormal*renormal*renormal;

	c=0;
	for (i=0; i<isize[0]; i++)
		for (j=0; j<isize[1]; j++)
			for (k=0; k<isize[2]; k++) {

				mesh[c] *= renormal;
				c++;
			}

	/////////////////////////////////////////

	int *nsendkey = (int*)pmalloc(sizeof(int) * nproc ,51);
	int *nsendisp = (int*)pmalloc(sizeof(int) * nproc ,52);
	int *nrecvkey = (int*)pmalloc(sizeof(int) * nproc ,53); 
	int *nrecdisp = (int*)pmalloc(sizeof(int) * nproc ,54);


	for (n=0; n<nproc; n++) {
		nsendkey[n] = 0;
		nsendisp[n] = 0;
		nrecvkey[n] = 0;
		nrecdisp[n] = 0;
	}	

	int l,m,q, d=0;
	c =0;
	int rk ,rj;

	for (m=0; m<isize[1]; m++) { 
		j = m + local_xmin[1] ;
		if (j >= NSIDE)
			j -= NSIDE;
		if (j < 0)
			j += NSIDE;

		for (q=0; q<isize[2]; q++) {
			k = q + local_xmin[2] ;
			if ( k >= NSIDE)				
				k -= NSIDE;
			if ( k < 0)
				k += NSIDE;

			rj = 0;
			while (j > MSIZE0[rj]) {
				rj++;	
			}

			rk = 0;
			while (k > MSIZE1[rk]) {
				rk++;	
			}
			d = rj + rk;

			nsendkey[d] += isize[0];
		}
	}

	nsendisp[0] = 0;
	for (n=1; n<nproc; n++) {
		nsendisp[n] = nsendkey[n-1] + nsendisp[n-1];
	}

	if (0 == rank)
		printf("meshsize = %d %d\n", meshsize,sizeof(MKey) );
	fflush(stdout);

	MKey *sendbuff = (MKey*) pmalloc(sizeof(MKey) * meshsize, 55);

	for (n=0; n<meshsize; n++) {
		sendbuff[n].x =  local_xmin[0];
		sendbuff[n].y =  local_xmin[1];
		sendbuff[n].z =  local_xmin[2];
		sendbuff[n].v = 0.0;
	}

	c = 0;
	for (l=0; l<isize[0]; l++) 
		for (m=0; m<isize[1]; m++) 
			for (q=0; q<isize[2]; q++) {

				i = l + local_xmin[0];
				j = m + local_xmin[1];
				k = q + local_xmin[2];				

				ii = i;
				jj = j;
				kk = k;

				if (i >= NSIDE)
					i -= NSIDE;
				if (i < 0)
					i += NSIDE;

				if (j >= NSIDE)
					j -= NSIDE;
				if (j < 0)
					j += NSIDE;

				if (k >= NSIDE)
					k -= NSIDE;
				if (k < 0) 
					k += NSIDE;
			
			rk = k/pside[1];
			if (rk > vproc[1])
				rk = vproc[1];

			rj = j/pside[0];
			if (rj > vproc[0])
				rj = vproc[0];
	
	
			rj = 0;
			while (j > MSIZE0[rj]) {
				rj++;	
			}
			rk = 0;
			while (k > MSIZE1[rk]) {
				rk++;	
			}


			d = rj + rk;


				sendbuff[nsendisp[d]].x = ii;
				sendbuff[nsendisp[d]].y = jj;
				sendbuff[nsendisp[d]].z = kk;
				sendbuff[nsendisp[d]].v = mesh[c];			
				nsendisp[d]++;	

				c++;

			}

	int local_sendcnt = c;
	////////////////////////////////////
	nsendisp[0] = 0;
	for (n=1; n<nproc; n++) {
		nsendisp[n] = nsendkey[n-1] + nsendisp[n-1];
	}

	MPI_Alltoall(nsendkey,1,MPI_INT,nrecvkey,1,MPI_INT,PM_COMM_WORLD);


	nrecdisp[0] = 0;
	int tot_recvkey = nrecvkey[0];	
	for (n=1; n<nproc; n++) {
		nrecdisp[n] = nrecdisp[n-1] + nrecvkey[n-1];
		tot_recvkey += nrecvkey[n];
	}



	MKey *recvbuff = (MKey*)pmalloc(sizeof(MKey)*tot_recvkey, 56);

	MPI_Status status;
//	MPI_Status *vstatus = (MPI_Status*)malloc(sizeof(MPI_Status)  * nproc);
//	MPI_Request *vreq  = (MPI_Request*)malloc(sizeof(MPI_Request) * nproc);

	MPI_Status *vstatus = (MPI_Status*)pmalloc(sizeof(MPI_Status)  * nproc, 57);
	MPI_Request *vreq  = (MPI_Request*) pmalloc(sizeof(MPI_Request) * nproc, 58);

	time_0 = dtime();

//	#define MYALLTOALL
#ifndef MYALLTOALL
/*	for (n=0; n<nproc; n++) {
		nsendkey[n] *= sizeof(MKey);
		nrecvkey[n] *= sizeof(MKey);
		nsendisp[n] *= sizeof(MKey);
		nrecdisp[n] *= sizeof(MKey);
	}

	MPI_Alltoallv(sendbuff, nsendkey, nsendisp, MPI_BYTE, recvbuff, nrecvkey, nrecdisp, MPI_BYTE, PM_COMM_WORLD);
*/
MPI_Alltoallv(sendbuff, nsendkey, nsendisp, strMKey, recvbuff, nrecvkey, nrecdisp, strMKey, PM_COMM_WORLD);
#else

	for (n=0; n<nsendkey[rank]; n++) {
		*(recvbuff+nrecdisp[rank] + n) = *(sendbuff+nsendisp[rank] + n);
	}
	for (n=0;n<nproc;n++) {
		if (n != rank && nsendkey[n] > 0) {
			//	MPI_Send(sendbuff+nsendisp[n], nsendkey[n]*sizeof(MKey), MPI_BYTE, n,n, PM_COMM_WORLD);
//			MPI_Isend(sendbuff+nsendisp[n], nsendkey[n]*sizeof(MKey), MPI_BYTE, n,n, PM_COMM_WORLD, &vreq[n]);
MPI_Isend(sendbuff+nsendisp[n], nsendkey[n], strMKey, n,n, PM_COMM_WORLD, &vreq[n]);
		}
	}

	for (n=0;n<nproc;n++) {
		if (n != rank && nrecvkey[n] > 0)
			//		MPI_Recv(recvbuff, nrecvkey[n], MPI_BYTE, n,rank, PM_COMM_WORLD, &status);
//			MPI_Recv(recvbuff+nrecdisp[n], nrecvkey[n]*sizeof(MKey), MPI_BYTE, n,rank, PM_COMM_WORLD, &status);
MPI_Recv(recvbuff+nrecdisp[n], nrecvkey[n], strMKey, n,rank, PM_COMM_WORLD, &status);
	}	
	for (n=0;n<nproc;n++) {
		if (n != rank && nsendkey[n] > 0)
			MPI_Wait(&vreq[n], &vstatus[n]);
	}
	MPI_Barrier(PM_COMM_WORLD);
#endif

	for (m=0; m<local_xsize[0]*local_xsize[1]*local_xsize[2]; m++) {
		data[m] = 0.0;
	}	

	for (n=0; n<tot_recvkey; n++) {

		i = recvbuff[n].x ;
		j = recvbuff[n].y ;
		k = recvbuff[n].z ;

		if (i < 0)
			i += NSIDE ;
		if (i >= NSIDE)
			i -= NSIDE;
		if (j < 0)
			j += NSIDE;
		if (j >= NSIDE)
			j -= NSIDE;
		if (k < 0)
			k += NSIDE;
		if (k >= NSIDE)
			k -= NSIDE;

		i -= local_xstart[0];
		j -= local_xstart[1];
		k -= local_xstart[2];

		data[(i*local_xsize[1]+j)*local_xsize[2]+k] += recvbuff[n].v;

	}


	int nside[3] = {NSIDE, NSIDE, NSIDE};
	double param[2] = { splitRadius, BOXSIZE};

#ifndef PMONLY
	convolution(data,nside,param);
#else
	conv_pmonly(data,nside,param);
#endif

	for (n=0; n<tot_recvkey; n++) {
		i = recvbuff[n].x ;
		j = recvbuff[n].y ;
		k = recvbuff[n].z ;

		if (i < 0)
			i += NSIDE ;
		if (i >= NSIDE)
			i -= NSIDE;

		if (j < 0)
			j += NSIDE;
		if (j >= NSIDE)
			j -= NSIDE;

		if (k < 0)
			k += NSIDE;
		if (k >= NSIDE)
			k -= NSIDE;

		i -= local_xstart[0];
		j -= local_xstart[1];
		k -= local_xstart[2];

		idx = (i*local_xsize[1]+j)*local_xsize[2]+k;

		recvbuff[n].v = data[(i*local_xsize[1]+j)*local_xsize[2]+k];

	}
	FILE *fd;
	char fname[80];


#ifndef MYALLTOALL

//	MPI_Alltoallv(recvbuff, nrecvkey, nrecdisp, MPI_BYTE, sendbuff, nsendkey, nsendisp, MPI_BYTE, PM_COMM_WORLD);
MPI_Alltoallv(recvbuff, nrecvkey, nrecdisp, strMKey, sendbuff, nsendkey, nsendisp, strMKey, PM_COMM_WORLD);
#else

	for (n=0; n<nsendkey[rank]; n++) {
		*(sendbuff+nsendisp[rank] + n) = *(recvbuff+nrecdisp[rank] + n)  ;
	}

	for (n=0;n<nproc;n++) {
		if (n != rank && nrecvkey[n] > 0)
			MPI_Isend(recvbuff+nrecdisp[n], nrecvkey[n]*sizeof(MKey), MPI_BYTE, n,n, PM_COMM_WORLD, &vreq[n]);
	}
	for (n=0;n<nproc;n++) {

		if (n != rank && nsendkey[n] > 0)
			MPI_Recv(sendbuff+nsendisp[n], nsendkey[n]*sizeof(MKey), MPI_BYTE, n,rank, PM_COMM_WORLD, &status);
	}
	for (n=0;n<nproc;n++) {
		if (n != rank && nrecvkey[n] >0)
			MPI_Wait(&vreq[n], &vstatus[n]);
	}

#endif

	pfree(vstatus, 57);
	pfree(vreq, 58);

	for (c=0; c<local_sendcnt; c++) {
		i = sendbuff[c].x;
		j = sendbuff[c].y;
		k = sendbuff[c].z;

		i -= local_xmin[0] ;
		j -= local_xmin[1] ;
		k -= local_xmin[2] ;

		idx = (i*isize[1] + j)*isize[2] + k;

		mesh[idx] = sendbuff[c].v ;
	}


	time_1 = dtime();

	double invx = 0.5*NSIDE/BOXSIZE;
	int nx, ny, nz;

	nx = isize[1]*isize[2];
	ny = isize[2];
	nz = 1;

	int idx1, idx2; 
	double dpx, dpxn, dpy, dpyn, dpz, dpzn;
	double dp[8];

	for (n=0; n<NPART; n++) {
		i = (int) (part[n].pos[0] *norm) ;	
		j = (int) (part[n].pos[1] *norm) ;	
		k = (int) (part[n].pos[2] *norm) ;	

		wi = (part[n].pos[0] - (i+0.5)*delta)*norm;
		if (wi > 0) {
			ii = i + 1;
		}
		else {
			wi = -wi;
			ii = i - 1;
		}
		win =  1.0 - wi;

		wj = (part[n].pos[1] - (j+0.5)*delta)*norm;
		if (wj > 0) {
			jj = j + 1;
		}
		else {
			wj = -wj;
			jj = j - 1;
		}
		wjn =  1.0 - wj;

		wk = (part[n].pos[2] - (k+0.5)*delta)*norm;
		if (wk > 0) {
			kk = k + 1;
		}
		else {
			wk = -wk;
			kk = k - 1;
		}
		wkn =  1.0 - wk;

		i -= local_xmin[0] ;	
		j -= local_xmin[1] ;	
		k -= local_xmin[2] ;	

		ii -= local_xmin[0] ;	
		jj -= local_xmin[1] ;	
		kk -= local_xmin[2] ;	

		idx = (i*isize[1] + j)*isize[2] + k;

		if (idx + nx > meshsize)
			printf(" error : %d %d %d\n", i, j, k);

		if (idx - nx < 0)
			printf(" error : %d %d %d\n", i, j, k);

			double f1 = 4.0/3.0;
			double f2 = 1.0/6.0;
			idx1 = ((i-1)*isize[1] + j)*isize[2] + k;
			idx2 = ((i+1)*isize[1] + j)*isize[2] + k;
			dp[0] =  f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-2)*isize[1] + j)*isize[2] + k;
			idx2 = ((i+2)*isize[1] + j)*isize[2] + k;
			dp[0] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-1)*isize[1] + j)*isize[2] + k;
			idx2 = ((ii+1)*isize[1] + j)*isize[2] + k;
			dp[1] = f1 *invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-2)*isize[1] + j)*isize[2] + k;
			idx2 = ((ii+2)*isize[1] + j)*isize[2] + k;
			dp[1] -= f2 * invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-1)*isize[1] + jj)*isize[2] + k;
			idx2 = ((i+1)*isize[1] + jj)*isize[2] + k;
			dp[2] = f1 * invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-2)*isize[1] + jj)*isize[2] + k;
			idx2 = ((i+2)*isize[1] + jj)*isize[2] + k;
			dp[2] -= f2 * invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-1)*isize[1] + jj)*isize[2] + k;
			idx2 = ((ii+1)*isize[1] + jj)*isize[2] + k;
			dp[3] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-2)*isize[1] + jj)*isize[2] + k;
			idx2 = ((ii+2)*isize[1] + jj)*isize[2] + k;
			dp[3] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-1)*isize[1] + j)*isize[2] + kk;
			idx2 = ((i+1)*isize[1] + j)*isize[2] + kk;
			dp[4] = f1* invx*(mesh[idx2] - mesh[idx1]);


			idx1 = ((i-2)*isize[1] + j)*isize[2] + kk;
			idx2 = ((i+2)*isize[1] + j)*isize[2] + kk;
			dp[4] -= f2 *invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-1)*isize[1] + j)*isize[2] + kk;
			idx2 = ((ii+1)*isize[1] + j)*isize[2] + kk;
			dp[5] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-2)*isize[1] + j)*isize[2] + kk;
			idx2 = ((ii+2)*isize[1] + j)*isize[2] + kk;
			dp[5] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-1)*isize[1] + jj)*isize[2] + kk;
			idx2 = ((i+1)*isize[1] + jj)*isize[2] + kk;
			dp[6] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-2)*isize[1] + jj)*isize[2] + kk;
			idx2 = ((i+2)*isize[1] + jj)*isize[2] + kk;
			dp[6]-= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-1)*isize[1] + jj)*isize[2] + kk;
			idx2 = ((ii+1)*isize[1] + jj)*isize[2] + kk;
			dp[7] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-2)*isize[1] + jj)*isize[2] + kk;
			idx2 = ((ii+2)*isize[1] + jj)*isize[2] + kk;
			dp[7]-= f2*invx*(mesh[idx2] - mesh[idx1]);

			part[n].acc_pm[0] = win*wjn*wkn*dp[0] 
				+ wi *wjn*wkn*dp[1]
				+ win*wj *wkn*dp[2]
				+ wi *wj *wkn*dp[3]
				+ win*wjn*wk *dp[4]
				+ wi *wjn*wk *dp[5]
				+ win*wj *wk *dp[6]
				+ wi *wj *wk *dp[7];

			idx1 = (i*isize[1] + j-1)*isize[2] + k;
			idx2 = (i*isize[1] + j+1)*isize[2] + k;
			dp[0] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + j-2)*isize[2] + k;
			idx2 = (i*isize[1] + j+2)*isize[2] + k;
			dp[0] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j-1)*isize[2] + k;
			idx2 = (ii*isize[1] + j+1)*isize[2] + k;
			dp[1] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j-2)*isize[2] + k;
			idx2 = (ii*isize[1] + j+2)*isize[2] + k;
			dp[1] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj-1)*isize[2] + k;
			idx2 = (i*isize[1] + jj+1)*isize[2] + k;
			dp[2] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj-2)*isize[2] + k;
			idx2 = (i*isize[1] + jj+2)*isize[2] + k;
			dp[2] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj-1)*isize[2] + k;
			idx2 = (ii*isize[1] + jj+1)*isize[2] + k;
			dp[3] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj-2)*isize[2] + k;
			idx2 = (ii*isize[1] + jj+2)*isize[2] + k;
			dp[3] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + j-1)*isize[2] + kk;
			idx2 = (i*isize[1] + j+1)*isize[2] + kk;
			dp[4] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + j-2)*isize[2] + kk;
			idx2 = (i*isize[1] + j+2)*isize[2] + kk;
			dp[4] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j-1)*isize[2] + kk;
			idx2 = (ii*isize[1] + j+1)*isize[2] + kk;
			dp[5] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j-2)*isize[2] + kk;
			idx2 = (ii*isize[1] + j+2)*isize[2] + kk;
			dp[5] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj-1)*isize[2] + kk;
			idx2 = (i*isize[1] + jj+1)*isize[2] + kk;
			dp[6] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj-2)*isize[2] + kk;
			idx2 = (i*isize[1] + jj+2)*isize[2] + kk;
			dp[6] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj-1)*isize[2] + kk;
			idx2 = (ii*isize[1] + jj+1)*isize[2] + kk;
			dp[7] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj-2)*isize[2] + kk;
			idx2 = (ii*isize[1] + jj+2)*isize[2] + kk;
			dp[7] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			part[n].acc_pm[1] = win*wjn*wkn*dp[0] 
				+ wi *wjn*wkn*dp[1]
				+ win*wj *wkn*dp[2]
				+ wi *wj *wkn*dp[3]
				+ win*wjn*wk *dp[4]
				+ wi *wjn*wk *dp[5]
				+ win*wj *wk *dp[6]
				+ wi *wj *wk *dp[7];



			idx1 = (i*isize[1] + j)*isize[2] + k-1;
			idx2 = (i*isize[1] + j)*isize[2] + k+1;
			dp[0] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + j)*isize[2] + k-2;
			idx2 = (i*isize[1] + j)*isize[2] + k+2;
			dp[0] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j)*isize[2] + k-1;
			idx2 = (ii*isize[1] + j)*isize[2] + k+1;
			dp[1] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j)*isize[2] + k-2;
			idx2 = (ii*isize[1] + j)*isize[2] + k+2;
			dp[1] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + jj)*isize[2] + k-1;
			idx2 = (i*isize[1] + jj)*isize[2] + k+1;
			dp[2] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj)*isize[2] + k-2;
			idx2 = (i*isize[1] + jj)*isize[2] + k+2;
			dp[2] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj)*isize[2] + k-1;
			idx2 = (ii*isize[1] + jj)*isize[2] + k+1;
			dp[3] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj)*isize[2] + k-2;
			idx2 = (ii*isize[1] + jj)*isize[2] + k+2;
			dp[3] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + j)*isize[2] + kk-1;
			idx2 = (i*isize[1] + j)*isize[2] + kk+1;
			dp[4] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + j)*isize[2] + kk-2;
			idx2 = (i*isize[1] + j)*isize[2] + kk+2;
			dp[4] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j)*isize[2] + kk-1;
			idx2 = (ii*isize[1] + j)*isize[2] + kk+1;
			dp[5] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j)*isize[2] + kk-2;
			idx2 = (ii*isize[1] + j)*isize[2] + kk+2;
			dp[5] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + jj)*isize[2] + kk-1;
			idx2 = (i*isize[1] + jj)*isize[2] + kk+1;
			dp[6] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj)*isize[2] + kk-2;
			idx2 = (i*isize[1] + jj)*isize[2] + kk+2;
			dp[6] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj)*isize[2] + kk-1;
			idx2 = (ii*isize[1] + jj)*isize[2] + kk+1;
			dp[7] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj)*isize[2] + kk-2;
			idx2 = (ii*isize[1] + jj)*isize[2] + kk+2;
			dp[7] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			part[n].acc_pm[2] = win*wjn*wkn*dp[0] 
				+ wi *wjn*wkn*dp[1]
				+ win*wj *wkn*dp[2]
				+ wi *wj *wkn*dp[3]
				+ win*wjn*wk *dp[4]
				+ wi *wjn*wk *dp[5]
				+ win*wj *wk *dp[6]
				+ wi *wj *wk *dp[7];

			//////////////////

	}


	if (mesh != NULL)
		pfree(mesh, 50);

	if (sendbuff != NULL )
		pfree(sendbuff, 55);

	if (nsendkey != NULL)
		pfree(nsendkey, 51);

	if (nsendisp != NULL)
		pfree(nsendisp, 52);

	if (nrecdisp != NULL)
		pfree(nrecdisp, 54); 

	if (nrecvkey != NULL )
		pfree(nrecvkey, 53); 

	if (recvbuff != NULL)
		pfree(recvbuff, 56);

}

#else

void partmesh()
{
	int nproc, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank );


	double time_0, time_1, t0, t1;
	time_0 = dtime();

	int n, c ;
	double domain_min[3] = {BOXSIZE,BOXSIZE,BOXSIZE};
	double domain_max[3] = {0.0, 0.0, 0.0};
	int local_xmin[3], local_xmax[3], isize[3];

	//GPTLstart("partmesh_MPI1");
	domain_min[0] = domain_min[1] = domain_min[2] = BOXSIZE;
	domain_max[0] = domain_max[1] = domain_max[2] = 0.0;

	for (n=0; n<NPART; n++ ) {
		if (domain_min[0] > part[n].pos[0] )
			domain_min[0] = part[n].pos[0];

		if (domain_min[1] > part[n].pos[1] )
			domain_min[1] = part[n].pos[1];

		if (domain_min[2] > part[n].pos[2] )
			domain_min[2] = part[n].pos[2];


		if (domain_max[0] < part[n].pos[0] )
			domain_max[0] = part[n].pos[0];

		if (domain_max[1] < part[n].pos[1] )
			domain_max[1] = part[n].pos[1];

		if (domain_max[2] < part[n].pos[2] )
			domain_max[2] = part[n].pos[2];
	}


//printf("domain max - [%d] %lf %lf %lf\n", PROC_RANK, domain_max[0], domain_max[1], domain_max[2]);
//printf("domain min - [%d] %lf %lf %lf\n", PROC_RANK, domain_min[0], domain_min[1], domain_min[2]);
//fflush(stdout);

	local_xmin[0] = (int) ((double)domain_min[0]*NSIDE / BOXSIZE);
	local_xmin[1] = (int) ((double)domain_min[1]*NSIDE / BOXSIZE);
	local_xmin[2] = (int) ((double)domain_min[2]*NSIDE / BOXSIZE);

	local_xmin[0] -= 1;
	local_xmin[1] -= 1;
	local_xmin[2] -= 1;

	local_xmin[0] -= 1;
	local_xmin[1] -= 1;
	local_xmin[2] -= 1;

	local_xmin[0] -= 1;
	local_xmin[1] -= 1;
	local_xmin[2] -= 1;


	local_xmax[0] = (int) ((double)domain_max[0]*NSIDE / BOXSIZE);
	local_xmax[1] = (int) ((double)domain_max[1]*NSIDE / BOXSIZE);
	local_xmax[2] = (int) ((double)domain_max[2]*NSIDE / BOXSIZE);

	isize[0] = local_xmax[0] - local_xmin[0]  +2 +1 +2;
	isize[1] = local_xmax[1] - local_xmin[1]  +2 +1 +2;
	isize[2] = local_xmax[2] - local_xmin[2]  +2 +1 +2;

//	printf(" p2 - %d %d %d\n", isize[0], isize[1], isize[2]);
//	printf(" p2[%d] - min: %d %d %d\n",PROC_RANK, local_xmin[0], local_xmin[1], local_xmin[2]);
//	fflush(stdout);

MPI_Barrier(MPI_COMM_WORLD);

	long meshsize = isize[0] * isize[1] * isize[2];
	double* mesh = (double*)pmalloc(sizeof(double) * meshsize, 50);

	for (n=0; n<meshsize; n++) 
		mesh[n] = 0.0;

	int i,j,k,idx;
	int ii, jj, kk;
	double norm = NSIDE/BOXSIZE;
	double delta = 1.0/norm;

	double wi, wj , wk, win, wjn, wkn;
	c=0;	
	for (n=0; n<NPART; n++) {
		//		i = (int) (part[n].pos[0]*norm) - local_xmin[0] ;	
		//		j = (int) (part[n].pos[1]*norm) - local_xmin[1] ;	
		//		k = (int) (part[n].pos[2]*norm) - local_xmin[2] ;	
		//		idx = (i*isize[1] + j)*isize[2] + k;
		//		mesh[idx] += 1.0;

		i = (int) (part[n].pos[0] *norm) ;	
		j = (int) (part[n].pos[1] *norm) ;	
		k = (int) (part[n].pos[2] *norm) ;	

		wi = (part[n].pos[0] - (i+0.5)*delta)*norm;

		if (wi > 0) {
			ii = i + 1;
		}
		else {
			wi = -wi;
			ii = i - 1;
		}
		win =  1.0 - wi;

		wj = (part[n].pos[1] - (j+0.5)*delta)*norm;
		if (wj > 0) {
			jj = j + 1;
		}
		else {
			wj = -wj;
			jj = j - 1;
		}
		wjn =  1.0 - wj;

		wk = (part[n].pos[2] - (k+0.5)*delta)*norm;
		if (wk > 0) {
			kk = k + 1;
		}
		else {
			wk = -wk;
			kk = k - 1;
		}
		wkn =  1.0 - wk;

		i -= local_xmin[0] ;	
		j -= local_xmin[1] ;	
		k -= local_xmin[2] ;	

		ii -= local_xmin[0] ;	
		jj -= local_xmin[1] ;	
		kk -= local_xmin[2] ;	


		idx = (i*isize[1] + j)*isize[2] + k;
		mesh[idx] += part[n].mass*win*wjn*wkn;	

		idx = (ii*isize[1] + j)*isize[2] + k;
		mesh[idx] += part[n].mass*wi *wjn*wkn;	

		idx = (i*isize[1] + jj)*isize[2] + k;
		mesh[idx] += part[n].mass*win*wj *wkn;	

		idx = (i*isize[1] + j)*isize[2] + kk;
		mesh[idx] += part[n].mass*win*wjn*wk ;	

		idx = (ii*isize[1] + jj)*isize[2] + k;
		mesh[idx] += part[n].mass*wi *wj *wkn;	

		idx = (ii*isize[1] + j)*isize[2] + kk;
		mesh[idx] += part[n].mass*wi *wjn*wk ;	

		idx = (i*isize[1] + jj)*isize[2] + kk;
		mesh[idx] += part[n].mass*win*wj *wk ;	

		idx = (ii*isize[1] + jj)*isize[2] + kk;
		mesh[idx] += part[n].mass* wi *wj *wk ;	

		c++;

		//		printf(" mesh[%d,%d,%d]=%lf\n", i,j,k, mesh[idx]);

	}

	//	printf(" c = %d npart = %d\n", c, NPART);
	double renormal = (NSIDE/BOXSIZE);
	renormal = renormal*renormal*renormal;

	c=0;
	for (i=0; i<isize[0]; i++)
		for (j=0; j<isize[1]; j++)
			for (k=0; k<isize[2]; k++) {

				mesh[c] *= renormal;
				c++;
			}

	////////////////////////////////////////////////////////////////////
	//
	/*
	   double pos[3] = { 0.07, 0.07, 0.07  };
	   for (n=0; n<meshsize; n++) 
	   mesh[n] = 0.0;

	   if (rank ==0 ){
	   i = (int) (pos[0]*norm)-local_xmin[0];	
	   j = (int) (pos[1]*norm)-local_xmin[1];	
	   k = (int) (pos[2]*norm)-local_xmin[2];	

	   idx = (i*isize[1] + j)*isize[2] + k;

	   mesh[idx] = 1.0;	
	   }
	   */
	//	
	//////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////

#ifdef CHECK_PARTMESH
	FILE *fd;
	char fname[90];
	sprintf(fname, "map.%d-%d", PROC_SIZE, PROC_RANK);
	fd = fopen(fname, "w");

	c=0;
	for (i=0; i<isize[0]; i++)
		for (j=0; j<isize[1]; j++)
			for (k=0; k<isize[2]; k++) {
				fprintf(fd, "%d %d %d %lf\n", i+local_xmin[0], j+local_xmin[1], k+local_xmin[2] ,mesh[c]);
				//                        printf("%d %d %d %e\n", i+local_xmin[0], j+local_xmin[1], k+local_xmin[2] ,mesh[c]);
				c++;
			}

	fclose(fd);
#endif
	/////////////////////////////////////////

	int *nsendkey = (int*)pmalloc(sizeof(int) * nproc ,51);
	int *nsendisp = (int*)pmalloc(sizeof(int) * nproc ,52);
	int *nrecvkey = (int*)pmalloc(sizeof(int) * nproc ,53); 
	int *nrecdisp = (int*)pmalloc(sizeof(int) * nproc ,54);

	for (n=0; n<nproc; n++) {
		nsendkey[n] = 0;
		nsendisp[n] = 0;
		nrecvkey[n] = 0;
		nrecdisp[n] = 0;
	}	

	int l,m,q, d=0;
	c =0;
	int rk ,rj;

	for (m=0; m<isize[1]; m++) { 
		j = m + local_xmin[1] ;
		if (j >= NSIDE)
			j -= NSIDE;
		if (j < 0)
			j += NSIDE;

		for (q=0; q<isize[2]; q++) {
			k = q + local_xmin[2] ;
			if ( k >= NSIDE)				
				k -= NSIDE;
			if ( k < 0)
				k += NSIDE;

			
	//		rk = k/pside[1];
	//		if (rk > vproc[1])
	//			rk = vproc[1];

	//		rj = j/pside[0];
	//		if (rj > vproc[0])
	//			rj = vproc[0];

	
			rj = 0;
			while (j > MSIZE0[rj]) {
//		printf(" j = %d %d - %d %d\n", j, rj, MSIZE0[0], MSIZE0[1]);
				rj++;	
			}

			rk = 0;
			while (k > MSIZE1[rk]) {
				rk++;	
			}
			d = rj + rk;

//	printf(" %d %d - %d %d = %d %d- %d\n", j, k,rj, rk, MSIZE0[rj], MSIZE1[rk] , d);
//			d = rj  + rk * vproc[0] ;

//			d = ( k/local_xsize[1] ) + (j/local_xsize[0]) * vproc[1];
	
//			d = ( k/local_xsize[2] ) + j/local_xsize[1] * vproc[1] ;	
		//	d = ( j/pside[0] )* vproc[1] + k/pside[1];	

			nsendkey[d] += isize[0];
		}
	}

	nsendisp[0] = 0;
	for (n=1; n<nproc; n++) {
		nsendisp[n] = nsendkey[n-1] + nsendisp[n-1];
	}

	MKey *sendbuff = (MKey*) pmalloc(sizeof(MKey) * meshsize, 55);
//printf(" mesh size = %ld \n", meshsize);
	for (n=0; n<meshsize; n++) {
		sendbuff[n].x =  local_xmin[0];
		sendbuff[n].y =  local_xmin[1];
		sendbuff[n].z =  local_xmin[2];
		sendbuff[n].v = 0.0;
	}

	c = 0;
	for (l=0; l<isize[0]; l++) 
		for (m=0; m<isize[1]; m++) 
			for (q=0; q<isize[2]; q++) {

				i = l + local_xmin[0];
				j = m + local_xmin[1];
				k = q + local_xmin[2];				

				ii = i;
				jj = j;
				kk = k;

				if (i >= NSIDE)
					i -= NSIDE;
				if (i < 0)
					i += NSIDE;

				if (j >= NSIDE)
					j -= NSIDE;
				if (j < 0)
					j += NSIDE;

				if (k >= NSIDE)
					k -= NSIDE;
				if (k < 0) 
					k += NSIDE;
			
			rk = k/pside[1];
			if (rk > vproc[1])
				rk = vproc[1];

			rj = j/pside[0];
			if (rj > vproc[0])
				rj = vproc[0];
	
	
			rj = 0;
			while (j > MSIZE0[rj]) {
				rj++;	
			}
			rk = 0;
			while (k > MSIZE1[rk]) {
				rk++;	
			}


			d = rj + rk;


//			d = rj  + rk * vproc[0] ;
			//	d = ( j/pside[0] )* vproc[1] + k/pside[1];

				sendbuff[nsendisp[d]].x = ii;
				sendbuff[nsendisp[d]].y = jj;
				sendbuff[nsendisp[d]].z = kk;
				sendbuff[nsendisp[d]].v = mesh[c];			
				nsendisp[d]++;	

				c++;

			}

	int local_sendcnt = c;
	////////////////////////////////////
	nsendisp[0] = 0;
	for (n=1; n<nproc; n++) {
		nsendisp[n] = nsendkey[n-1] + nsendisp[n-1];
	}

	//	for (n=0; n<nproc; n++) {
	//		nsendisp[n] -= nsendkey[n];
	//	}

	/*
	   time_1 = dtime();
	   if (0==rank) {
	   printf(" dt prep mesh : %lf, %d\n", time_1 - time_0, NPART);
	   fflush(stdout);
	   }
	   time_0 = dtime();
	   */
	MPI_Alltoall(nsendkey,1,MPI_INT,nrecvkey,1,MPI_INT,MPI_COMM_WORLD);

	/*
	   time_1 = dtime();

	   if (0==rank) {
	   printf(" dt mpi_alltoall info : %lf\n", time_1 - time_0);
	   fflush(stdout);
	   }

*/


	nrecdisp[0] = 0;
	int tot_recvkey = nrecvkey[0];	
	for (n=1; n<nproc; n++) {
		nrecdisp[n] = nrecdisp[n-1] + nrecvkey[n-1];
//printf("[%d] tot recvkey = %d %d \n", rank, nrecdisp[n], tot_recvkey);
//fflush(stdout);

		tot_recvkey += nrecvkey[n];
	}



	MKey *recvbuff = (MKey*)pmalloc(sizeof(MKey)*tot_recvkey, 56);

	MPI_Status status;
//	MPI_Status *vstatus = (MPI_Status*)malloc(sizeof(MPI_Status)  * nproc);
//	MPI_Request *vreq  = (MPI_Request*)malloc(sizeof(MPI_Request) * nproc);

	MPI_Status *vstatus = (MPI_Status*)pmalloc(sizeof(MPI_Status)  * nproc, 57);
	MPI_Request *vreq  = (MPI_Request*) pmalloc(sizeof(MPI_Request) * nproc, 58);

	time_0 = dtime();
	//#define MYALLTOALL
#ifndef MYALLTOALL
/*
	for (n=0; n<nproc; n++) {
		nsendkey[n] *= sizeof(MKey);
		nrecvkey[n] *= sizeof(MKey);
		nsendisp[n] *= sizeof(MKey);
		nrecdisp[n] *= sizeof(MKey);
	}

	MPI_Alltoallv(sendbuff, nsendkey, nsendisp, MPI_BYTE, recvbuff, nrecvkey, nrecdisp, MPI_BYTE, MPI_COMM_WORLD);
*/
MPI_Alltoallv(sendbuff, nsendkey, nsendisp, strMKey, recvbuff, nrecvkey, nrecdisp, strMKey, MPI_COMM_WORLD);
#else

	for (n=0; n<nsendkey[rank]; n++) {
		*(recvbuff+nrecdisp[rank] + n) = *(sendbuff+nsendisp[rank] + n);
	}
	for (n=0;n<nproc;n++) {
		if (n != rank && nsendkey[n] > 0) {
			//	MPI_Send(sendbuff+nsendisp[n], nsendkey[n]*sizeof(MKey), MPI_BYTE, n,n, MPI_COMM_WORLD);
			MPI_Isend(sendbuff+nsendisp[n], nsendkey[n]*sizeof(MKey), MPI_BYTE, n,n, MPI_COMM_WORLD, &vreq[n]);
		}
	}

	for (n=0;n<nproc;n++) {
		if (n != rank && nrecvkey[n] > 0)
			//		MPI_Recv(recvbuff, nrecvkey[n], MPI_BYTE, n,rank, MPI_COMM_WORLD, &status);
			MPI_Recv(recvbuff+nrecdisp[n], nrecvkey[n]*sizeof(MKey), MPI_BYTE, n,rank, MPI_COMM_WORLD, &status);
	}	
	for (n=0;n<nproc;n++) {
		if (n != rank && nsendkey[n] > 0)
			MPI_Wait(&vreq[n], &vstatus[n]);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	for (m=0; m<local_xsize[0]*local_xsize[1]*local_xsize[2]; m++) {
		data[m] = 0.0;
	}	

	for (n=0; n<tot_recvkey; n++) {

		i = recvbuff[n].x ;
		j = recvbuff[n].y ;
		k = recvbuff[n].z ;

		if (i < 0)
			i += NSIDE ;
		if (i >= NSIDE)
			i -= NSIDE;
		if (j < 0)
			j += NSIDE;
		if (j >= NSIDE)
			j -= NSIDE;
		if (k < 0)
			k += NSIDE;
		if (k >= NSIDE)
			k -= NSIDE;

		i -= local_xstart[0];
		j -= local_xstart[1];
		k -= local_xstart[2];

		data[(i*local_xsize[1]+j)*local_xsize[2]+k] += recvbuff[n].v;

	}

	//GPTLstop("partmesh_MPI1");

	int nside[3] = {NSIDE, NSIDE, NSIDE};
	double param[2] = { splitRadius, BOXSIZE};

	//GPTLstart("partmesh_conv");
#ifndef PMONLY
	convolution(data,nside,param);
#else
	conv_pmonly(data,nside,param);
#endif
	//GPTLstop("partmesh_conv");

	//GPTLstart("partmesh_MPI2");
	for (n=0; n<tot_recvkey; n++) {
		i = recvbuff[n].x ;
		j = recvbuff[n].y ;
		k = recvbuff[n].z ;


		if (i < 0)
			i += NSIDE ;
		if (i >= NSIDE)
			i -= NSIDE;

		if (j < 0)
			j += NSIDE;
		if (j >= NSIDE)
			j -= NSIDE;

		if (k < 0)
			k += NSIDE;
		if (k >= NSIDE)
			k -= NSIDE;

		i -= local_xstart[0];
		j -= local_xstart[1];
		k -= local_xstart[2];

		idx = (i*local_xsize[1]+j)*local_xsize[2]+k;

		recvbuff[n].v = data[(i*local_xsize[1]+j)*local_xsize[2]+k];

	}
	FILE *fd;
	char fname[80];


#ifndef MYALLTOALL

	MPI_Alltoallv(recvbuff, nrecvkey, nrecdisp, strMKey, sendbuff, nsendkey, nsendisp, strMKey, MPI_COMM_WORLD);
//	MPI_Alltoallv(recvbuff, nrecvkey, nrecdisp, MPI_BYTE, sendbuff, nsendkey, nsendisp, MPI_BYTE, MPI_COMM_WORLD);
#else

	for (n=0; n<nsendkey[rank]; n++) {
		*(sendbuff+nsendisp[rank] + n) = *(recvbuff+nrecdisp[rank] + n)  ;
	}

	for (n=0;n<nproc;n++) {
		if (n != rank && nrecvkey[n] > 0)
			MPI_Isend(recvbuff+nrecdisp[n], nrecvkey[n]*sizeof(MKey), MPI_BYTE, n,n, MPI_COMM_WORLD, &vreq[n]);
	}
	for (n=0;n<nproc;n++) {
		if (n != rank && nsendkey[n] > 0)
			MPI_Recv(sendbuff+nsendisp[n], nsendkey[n]*sizeof(MKey), MPI_BYTE, n,rank, MPI_COMM_WORLD, &status);
	}
	for (n=0;n<nproc;n++) {
		if (n != rank && nrecvkey[n] >0)
			MPI_Wait(&vreq[n], &vstatus[n]);
	}
//	pfree(vstatus, 57);
//	pfree(vreq, 58);

#endif
	pfree(vstatus, 57);
	pfree(vreq, 58);


	for (c=0; c<local_sendcnt; c++) {
		i = sendbuff[c].x;
		j = sendbuff[c].y;
		k = sendbuff[c].z;

		i -= local_xmin[0] ;
		j -= local_xmin[1] ;
		k -= local_xmin[2] ;

		idx = (i*isize[1] + j)*isize[2] + k;

		mesh[idx] = sendbuff[c].v ;
	}


	time_1 = dtime();

	//	FILE *fd;
	//	char fname[90];

	//#define CHECK_PARTMESH
#ifdef CHECK_PARTMESH
	sprintf(fname, "pot2.%d-%d", PROC_SIZE, PROC_RANK); 
	fd = fopen(fname, "w");

	c=0;
	for (i=0; i<isize[0]; i++)
		for (j=0; j<isize[1]; j++)
			for (k=0; k<isize[2]; k++) {
				fprintf(fd, "%d %d %d %e\n", i+local_xmin[0],
						j+local_xmin[1], k+local_xmin[2] ,mesh[c]);
				c++;
			}

	fclose(fd);
#endif
	//#undef CHECK_PARTMESH

	//GPTLstop("partmesh_MPI2");
	//GPTLstart("partmesh_accpm");

	double invx = 0.5*NSIDE/BOXSIZE;
	int nx, ny, nz;

	nx = isize[1]*isize[2];
	ny = isize[2];
	nz = 1;

	int idx1, idx2; 
	double dpx, dpxn, dpy, dpyn, dpz, dpzn;
	double dp[8];

	for (n=0; n<NPART; n++) {
		i = (int) (part[n].pos[0] *norm) ;	
		j = (int) (part[n].pos[1] *norm) ;	
		k = (int) (part[n].pos[2] *norm) ;	

		wi = (part[n].pos[0] - (i+0.5)*delta)*norm;
		if (wi > 0) {
			ii = i + 1;
		}
		else {
			wi = -wi;
			ii = i - 1;
		}
		win =  1.0 - wi;

		wj = (part[n].pos[1] - (j+0.5)*delta)*norm;
		if (wj > 0) {
			jj = j + 1;
		}
		else {
			wj = -wj;
			jj = j - 1;
		}
		wjn =  1.0 - wj;

		wk = (part[n].pos[2] - (k+0.5)*delta)*norm;
		if (wk > 0) {
			kk = k + 1;
		}
		else {
			wk = -wk;
			kk = k - 1;
		}
		wkn =  1.0 - wk;

		i -= local_xmin[0] ;	
		j -= local_xmin[1] ;	
		k -= local_xmin[2] ;	

		ii -= local_xmin[0] ;	
		jj -= local_xmin[1] ;	
		kk -= local_xmin[2] ;	

		idx = (i*isize[1] + j)*isize[2] + k;

		if (idx + nx > meshsize)
			printf(" error : %d %d %d\n", i, j, k);

		if (idx - nx < 0)
			printf(" error : %d %d %d\n", i, j, k);

	


			double f1 = 4.0/3.0;
			double f2 = 1.0/6.0;
			idx1 = ((i-1)*isize[1] + j)*isize[2] + k;
			idx2 = ((i+1)*isize[1] + j)*isize[2] + k;
			dp[0] =  f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-2)*isize[1] + j)*isize[2] + k;
			idx2 = ((i+2)*isize[1] + j)*isize[2] + k;
			dp[0] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-1)*isize[1] + j)*isize[2] + k;
			idx2 = ((ii+1)*isize[1] + j)*isize[2] + k;
			dp[1] = f1 *invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-2)*isize[1] + j)*isize[2] + k;
			idx2 = ((ii+2)*isize[1] + j)*isize[2] + k;
			dp[1] -= f2 * invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-1)*isize[1] + jj)*isize[2] + k;
			idx2 = ((i+1)*isize[1] + jj)*isize[2] + k;
			dp[2] = f1 * invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-2)*isize[1] + jj)*isize[2] + k;
			idx2 = ((i+2)*isize[1] + jj)*isize[2] + k;
			dp[2] -= f2 * invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-1)*isize[1] + jj)*isize[2] + k;
			idx2 = ((ii+1)*isize[1] + jj)*isize[2] + k;
			dp[3] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-2)*isize[1] + jj)*isize[2] + k;
			idx2 = ((ii+2)*isize[1] + jj)*isize[2] + k;
			dp[3] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-1)*isize[1] + j)*isize[2] + kk;
			idx2 = ((i+1)*isize[1] + j)*isize[2] + kk;
			dp[4] = f1* invx*(mesh[idx2] - mesh[idx1]);


			idx1 = ((i-2)*isize[1] + j)*isize[2] + kk;
			idx2 = ((i+2)*isize[1] + j)*isize[2] + kk;
			dp[4] -= f2 *invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-1)*isize[1] + j)*isize[2] + kk;
			idx2 = ((ii+1)*isize[1] + j)*isize[2] + kk;
			dp[5] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-2)*isize[1] + j)*isize[2] + kk;
			idx2 = ((ii+2)*isize[1] + j)*isize[2] + kk;
			dp[5] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-1)*isize[1] + jj)*isize[2] + kk;
			idx2 = ((i+1)*isize[1] + jj)*isize[2] + kk;
			dp[6] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-2)*isize[1] + jj)*isize[2] + kk;
			idx2 = ((i+2)*isize[1] + jj)*isize[2] + kk;
			dp[6]-= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-1)*isize[1] + jj)*isize[2] + kk;
			idx2 = ((ii+1)*isize[1] + jj)*isize[2] + kk;
			dp[7] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-2)*isize[1] + jj)*isize[2] + kk;
			idx2 = ((ii+2)*isize[1] + jj)*isize[2] + kk;
			dp[7]-= f2*invx*(mesh[idx2] - mesh[idx1]);

			part[n].acc_pm[0] = win*wjn*wkn*dp[0] 
				+ wi *wjn*wkn*dp[1]
				+ win*wj *wkn*dp[2]
				+ wi *wj *wkn*dp[3]
				+ win*wjn*wk *dp[4]
				+ wi *wjn*wk *dp[5]
				+ win*wj *wk *dp[6]
				+ wi *wj *wk *dp[7];

			idx1 = (i*isize[1] + j-1)*isize[2] + k;
			idx2 = (i*isize[1] + j+1)*isize[2] + k;
			dp[0] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + j-2)*isize[2] + k;
			idx2 = (i*isize[1] + j+2)*isize[2] + k;
			dp[0] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j-1)*isize[2] + k;
			idx2 = (ii*isize[1] + j+1)*isize[2] + k;
			dp[1] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j-2)*isize[2] + k;
			idx2 = (ii*isize[1] + j+2)*isize[2] + k;
			dp[1] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj-1)*isize[2] + k;
			idx2 = (i*isize[1] + jj+1)*isize[2] + k;
			dp[2] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj-2)*isize[2] + k;
			idx2 = (i*isize[1] + jj+2)*isize[2] + k;
			dp[2] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj-1)*isize[2] + k;
			idx2 = (ii*isize[1] + jj+1)*isize[2] + k;
			dp[3] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj-2)*isize[2] + k;
			idx2 = (ii*isize[1] + jj+2)*isize[2] + k;
			dp[3] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + j-1)*isize[2] + kk;
			idx2 = (i*isize[1] + j+1)*isize[2] + kk;
			dp[4] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + j-2)*isize[2] + kk;
			idx2 = (i*isize[1] + j+2)*isize[2] + kk;
			dp[4] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j-1)*isize[2] + kk;
			idx2 = (ii*isize[1] + j+1)*isize[2] + kk;
			dp[5] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j-2)*isize[2] + kk;
			idx2 = (ii*isize[1] + j+2)*isize[2] + kk;
			dp[5] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj-1)*isize[2] + kk;
			idx2 = (i*isize[1] + jj+1)*isize[2] + kk;
			dp[6] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj-2)*isize[2] + kk;
			idx2 = (i*isize[1] + jj+2)*isize[2] + kk;
			dp[6] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj-1)*isize[2] + kk;
			idx2 = (ii*isize[1] + jj+1)*isize[2] + kk;
			dp[7] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj-2)*isize[2] + kk;
			idx2 = (ii*isize[1] + jj+2)*isize[2] + kk;
			dp[7] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			part[n].acc_pm[1] = win*wjn*wkn*dp[0] 
				+ wi *wjn*wkn*dp[1]
				+ win*wj *wkn*dp[2]
				+ wi *wj *wkn*dp[3]
				+ win*wjn*wk *dp[4]
				+ wi *wjn*wk *dp[5]
				+ win*wj *wk *dp[6]
				+ wi *wj *wk *dp[7];



			idx1 = (i*isize[1] + j)*isize[2] + k-1;
			idx2 = (i*isize[1] + j)*isize[2] + k+1;
			dp[0] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + j)*isize[2] + k-2;
			idx2 = (i*isize[1] + j)*isize[2] + k+2;
			dp[0] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j)*isize[2] + k-1;
			idx2 = (ii*isize[1] + j)*isize[2] + k+1;
			dp[1] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j)*isize[2] + k-2;
			idx2 = (ii*isize[1] + j)*isize[2] + k+2;
			dp[1] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + jj)*isize[2] + k-1;
			idx2 = (i*isize[1] + jj)*isize[2] + k+1;
			dp[2] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj)*isize[2] + k-2;
			idx2 = (i*isize[1] + jj)*isize[2] + k+2;
			dp[2] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj)*isize[2] + k-1;
			idx2 = (ii*isize[1] + jj)*isize[2] + k+1;
			dp[3] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj)*isize[2] + k-2;
			idx2 = (ii*isize[1] + jj)*isize[2] + k+2;
			dp[3] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + j)*isize[2] + kk-1;
			idx2 = (i*isize[1] + j)*isize[2] + kk+1;
			dp[4] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + j)*isize[2] + kk-2;
			idx2 = (i*isize[1] + j)*isize[2] + kk+2;
			dp[4] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j)*isize[2] + kk-1;
			idx2 = (ii*isize[1] + j)*isize[2] + kk+1;
			dp[5] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j)*isize[2] + kk-2;
			idx2 = (ii*isize[1] + j)*isize[2] + kk+2;
			dp[5] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + jj)*isize[2] + kk-1;
			idx2 = (i*isize[1] + jj)*isize[2] + kk+1;
			dp[6] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj)*isize[2] + kk-2;
			idx2 = (i*isize[1] + jj)*isize[2] + kk+2;
			dp[6] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj)*isize[2] + kk-1;
			idx2 = (ii*isize[1] + jj)*isize[2] + kk+1;
			dp[7] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj)*isize[2] + kk-2;
			idx2 = (ii*isize[1] + jj)*isize[2] + kk+2;
			dp[7] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			part[n].acc_pm[2] = win*wjn*wkn*dp[0] 
				+ wi *wjn*wkn*dp[1]
				+ win*wj *wkn*dp[2]
				+ wi *wj *wkn*dp[3]
				+ win*wjn*wk *dp[4]
				+ wi *wjn*wk *dp[5]
				+ win*wj *wk *dp[6]
				+ wi *wj *wk *dp[7];

			//////////////////

		

	}

	//GPTLstop("partmesh_accpm");

	if (mesh != NULL)
		pfree(mesh, 50);

	if (sendbuff != NULL )
		pfree(sendbuff, 55);

	if (nsendkey != NULL)
		pfree(nsendkey, 51);

	if (nsendisp != NULL)
		pfree(nsendisp, 52);

	if (nrecdisp != NULL)
		pfree(nrecdisp, 54); 

	if (nrecvkey != NULL )
		pfree(nrecvkey, 53); 

	if (recvbuff != NULL)
		pfree(recvbuff, 56);

}

#endif

void powerspectrum(char powname[]) 
{
	int nproc, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank );


	double time_0, time_1, t0, t1;
	time_0 = dtime();

	int n, c ;
	double domain_min[3] = {BOXSIZE,BOXSIZE,BOXSIZE};
	double domain_max[3] = {0.0, 0.0, 0.0};
	int local_xmin[3], local_xmax[3], isize[3];

	domain_min[0] = domain_min[1] = domain_min[2] = BOXSIZE;
	domain_max[0] = domain_max[1] = domain_max[3] = 0.0;

	for (n=0; n<NPART; n++ ) {
		if (domain_min[0] > part[n].pos[0] )
			domain_min[0] = part[n].pos[0];

		if (domain_min[1] > part[n].pos[1] )
			domain_min[1] = part[n].pos[1];

		if (domain_min[2] > part[n].pos[2] )
			domain_min[2] = part[n].pos[2];


		if (domain_max[0] < part[n].pos[0] )
			domain_max[0] = part[n].pos[0];

		if (domain_max[1] < part[n].pos[1] )
			domain_max[1] = part[n].pos[1];

		if (domain_max[2] < part[n].pos[2] )
			domain_max[2] = part[n].pos[2];
	}

	//printf("domain max - [%d] %lf %lf %lf\n", PROC_RANK, domain_max[0], domain_max[1], domain_max[2]);
	//printf("domain min - [%d] %lf %lf %lf\n", PROC_RANK, domain_min[0], domain_min[1], domain_min[2]);

	local_xmin[0] = (int) ((double)domain_min[0]*NSIDE / BOXSIZE);
	local_xmin[1] = (int) ((double)domain_min[1]*NSIDE / BOXSIZE);
	local_xmin[2] = (int) ((double)domain_min[2]*NSIDE / BOXSIZE);

	local_xmin[0] -= 1;
	local_xmin[1] -= 1;
	local_xmin[2] -= 1;

	local_xmin[0] -= 1;
	local_xmin[1] -= 1;
	local_xmin[2] -= 1;


	local_xmax[0] = (int) ((double)domain_max[0]*NSIDE / BOXSIZE);
	local_xmax[1] = (int) ((double)domain_max[1]*NSIDE / BOXSIZE);
	local_xmax[2] = (int) ((double)domain_max[2]*NSIDE / BOXSIZE);

	isize[0] = local_xmax[0] - local_xmin[0]  + 2 +1 ;
	isize[1] = local_xmax[1] - local_xmin[1]  + 2 +1 ;
	isize[2] = local_xmax[2] - local_xmin[2]  + 2 +1 ;

	//	printf(" p2 - %d %d %d\n", isize[0], isize[1], isize[2]);
	//	printf(" p2[%d] - min: %d %d %d\n",PROC_RANK, local_xmin[0], local_xmin[1], local_xmin[2]);

	long meshsize = isize[0] * isize[1] * isize[2];
	double* mesh = (double*)pmalloc(sizeof(double) * meshsize, 59);

	for (n=0; n<meshsize; n++) 
		mesh[n] = 0.0;

	int i,j,k,idx;
	int ii, jj, kk;
	double norm = NSIDE/BOXSIZE;
	double delta = 1.0/norm;

	double wi, wj , wk, win, wjn, wkn;
	c=0;	
	for (n=0; n<NPART; n++) {
		//		i = (int) (part[n].pos[0]*norm) - local_xmin[0] ;	
		//		j = (int) (part[n].pos[1]*norm) - local_xmin[1] ;	
		//		k = (int) (part[n].pos[2]*norm) - local_xmin[2] ;	
		//		idx = (i*isize[1] + j)*isize[2] + k;
		//		mesh[idx] += 1.0;

		i = (int) (part[n].pos[0] *norm) ;	
		j = (int) (part[n].pos[1] *norm) ;	
		k = (int) (part[n].pos[2] *norm) ;	

		wi = (part[n].pos[0] - (i+0.5)*delta)*norm;
		if (wi > 0) {
			ii = i + 1;
		}
		else {
			wi = -wi;
			ii = i - 1;
		}
		win =  1.0 - wi;

		wj = (part[n].pos[1] - (j+0.5)*delta)*norm;
		if (wj > 0) {
			jj = j + 1;
		}
		else {
			wj = -wj;
			jj = j - 1;
		}
		wjn =  1.0 - wj;

		wk = (part[n].pos[2] - (k+0.5)*delta)*norm;
		if (wk > 0) {
			kk = k + 1;
		}
		else {
			wk = -wk;
			kk = k - 1;
		}
		wkn =  1.0 - wk;

		i -= local_xmin[0] ;	
		j -= local_xmin[1] ;	
		k -= local_xmin[2] ;	

		ii -= local_xmin[0] ;	
		jj -= local_xmin[1] ;	
		kk -= local_xmin[2] ;	


		idx = (i*isize[1] + j)*isize[2] + k;
		mesh[idx] += part[n].mass*win*wjn*wkn;	


		idx = (ii*isize[1] + j)*isize[2] + k;
		mesh[idx] += part[n].mass*wi *wjn*wkn;	

		idx = (i*isize[1] + jj)*isize[2] + k;
		mesh[idx] += part[n].mass*win*wj *wkn;	

		idx = (i*isize[1] + j)*isize[2] + kk;
		mesh[idx] += part[n].mass*win*wjn*wk ;	

		idx = (ii*isize[1] + jj)*isize[2] + k;
		mesh[idx] += part[n].mass*wi *wj *wkn;	

		idx = (ii*isize[1] + j)*isize[2] + kk;
		mesh[idx] += part[n].mass*wi *wjn*wk ;	

		idx = (i*isize[1] + jj)*isize[2] + kk;
		mesh[idx] += part[n].mass*win*wj *wk ;	

		idx = (ii*isize[1] + jj)*isize[2] + kk;
		mesh[idx] += part[n].mass* wi *wj *wk ;	

		c++;

		//		printf(" mesh[%d,%d,%d]=%lf\n", i,j,k, mesh[idx]);

	}

	//	printf(" c = %d npart = %d\n", c, NPART);
	double renormal = (NSIDE/BOXSIZE);
	renormal = renormal*renormal*renormal;

	c=0;
	for (i=0; i<isize[0]; i++)
		for (j=0; j<isize[1]; j++)
			for (k=0; k<isize[2]; k++) {

				mesh[c] *= renormal;
				c++;
			}

	////////////////////////////////////////////////////////////////////
	//
	/*
	   double pos[3] = { 0.07, 0.07, 0.07  };
	   for (n=0; n<meshsize; n++) 
	   mesh[n] = 0.0;

	   if (rank ==0 ){
	   i = (int) (pos[0]*norm)-local_xmin[0];	
	   j = (int) (pos[1]*norm)-local_xmin[1];	
	   k = (int) (pos[2]*norm)-local_xmin[2];	

	   idx = (i*isize[1] + j)*isize[2] + k;

	   mesh[idx] = 1.0;	
	   }
	   */
	//	
	//////////////////////////////////////////////////////////////////////



	/////////////////////////////////////////

#ifdef CHECK_PARTMESH
	FILE *fd;
	char fname[80];
	sprintf(fname, "map.%d-%d", PROC_SIZE, PROC_RANK);
	fd = fopen(fname, "w");

	c=0;
	for (i=0; i<isize[0]; i++)
		for (j=0; j<isize[1]; j++)
			for (k=0; k<isize[2]; k++) {
				fprintf(fd, "%d %d %d %e\n", i+local_xmin[0], j+local_xmin[1], k+local_xmin[2] ,mesh[c]);

				c++;
			}

	fclose(fd);
#endif	



	/////////////////////////////////////////

	int *nsendkey = (int*)pmalloc(sizeof(int) * nproc ,60);
	int *nsendisp = (int*)pmalloc(sizeof(int) * nproc ,61);
	int *nrecvkey = (int*)pmalloc(sizeof(int) * nproc ,62); 
	int *nrecdisp = (int*)pmalloc(sizeof(int) * nproc ,63);

	for (n=0; n<nproc; n++) {
		nsendkey[n] = 0;
		nsendisp[n] = 0;
		nrecvkey[n] = 0;
		nrecdisp[n] = 0;
	}	

	int l,m,q, d=0;
	c =0;

	for (m=0; m<isize[1]; m++) { 
		j = m + local_xmin[1] ;
		if (j >= NSIDE)
			j -= NSIDE;
		if (j < 0)
			j += NSIDE;

		for (q=0; q<isize[2]; q++) {
			k = q + local_xmin[2] ;
			if ( k >= NSIDE)				
				k -= NSIDE;
			if ( k < 0)
				k += NSIDE;


			d = ( j/pside[0] )* vproc[1] + k/pside[1];	

			nsendkey[d] += isize[0];
		}
	}

	nsendisp[0] = 0;
	for (n=1; n<nproc; n++) {
		nsendisp[n] = nsendkey[n-1] + nsendisp[n-1];
	}

	MKey *sendbuff = (MKey*) pmalloc(sizeof(MKey) * meshsize, 64);


	for (n=0; n<meshsize; n++) {
		sendbuff[n].x =  local_xmin[0];
		sendbuff[n].y =  local_xmin[1];
		sendbuff[n].z =  local_xmin[2];
		sendbuff[n].v = 0.0;
	}

	c = 0;
	for (l=0; l<isize[0]; l++) 
		for (m=0; m<isize[1]; m++) 
			for (q=0; q<isize[2]; q++) {

				i = l + local_xmin[0];
				j = m + local_xmin[1];
				k = q + local_xmin[2];				

				ii = i;
				jj = j;
				kk = k;

				if (i >= NSIDE)
					i -= NSIDE;
				if (i < 0)
					i += NSIDE;

				if (j >= NSIDE)
					j -= NSIDE;
				if (j < 0)
					j += NSIDE;

				if (k >= NSIDE)
					k -= NSIDE;
				if (k < 0) 
					k += NSIDE;

				d = ( j/pside[0] )* vproc[1] + k/pside[1];

				sendbuff[nsendisp[d]].x = ii;
				sendbuff[nsendisp[d]].y = jj;
				sendbuff[nsendisp[d]].z = kk;
				sendbuff[nsendisp[d]].v = mesh[c];			
				nsendisp[d]++;	

				c++;

			}

	int local_sendcnt = c;
	////////////////////////////////////
	nsendisp[0] = 0;
	for (n=1; n<nproc; n++) {
		nsendisp[n] = nsendkey[n-1] + nsendisp[n-1];
	}

	//	for (n=0; n<nproc; n++) {
	//		nsendisp[n] -= nsendkey[n];
	//	}
	/*
	   time_1 = dtime();

	   if (0==rank) {
	   printf(" dt prep mesh : %lf, %d\n", time_1 - time_0, NPART);
	   fflush(stdout);
	   }
	   time_0 = dtime();
	   */
	MPI_Alltoall(nsendkey,1,MPI_INT,nrecvkey,1,MPI_INT,MPI_COMM_WORLD);

	/*
	   time_1 = dtime();

	   if (0==rank) {
	   printf(" dt mpi_alltoall info : %lf\n", time_1 - time_0);
	   fflush(stdout);
	   }
	   */

	nrecdisp[0] = 0;
	int tot_recvkey = nrecvkey[0];	
	for (n=1; n<nproc; n++) {
		nrecdisp[n] = nrecdisp[n-1] + nrecvkey[n-1];
		tot_recvkey += nrecvkey[n];
	}



	MKey *recvbuff = (MKey*)pmalloc(sizeof(MKey)*tot_recvkey, 65);

	MPI_Status status;
	MPI_Status *vstatus = (MPI_Status*)pmalloc(sizeof(MPI_Status)  * nproc, 66);
	MPI_Request *vreq  = (MPI_Request*) pmalloc(sizeof(MPI_Request) * nproc,67);

	time_0 = dtime();
	//#define MYALLTOALL
#ifndef MYALLTOALL

	for (n=0; n<nproc; n++) {
		nsendkey[n] *= sizeof(MKey);
		nrecvkey[n] *= sizeof(MKey);
		nsendisp[n] *= sizeof(MKey);
		nrecdisp[n] *= sizeof(MKey);
	}

	MPI_Alltoallv(sendbuff, nsendkey, nsendisp, MPI_BYTE, recvbuff, nrecvkey, nrecdisp, MPI_BYTE, MPI_COMM_WORLD);\

#else

		for (n=0; n<nsendkey[rank]; n++) {
			*(recvbuff+nrecdisp[rank] + n) = *(sendbuff+nsendisp[rank] + n);
		}
	for (n=0;n<nproc;n++) {
		if (n != rank && nsendkey[n] > 0) {
			//	MPI_Send(sendbuff+nsendisp[n], nsendkey[n]*sizeof(MKey), MPI_BYTE, n,n, MPI_COMM_WORLD);
			MPI_Isend(sendbuff+nsendisp[n], nsendkey[n]*sizeof(MKey), MPI_BYTE, n,n, MPI_COMM_WORLD, &vreq[n]);
		}
	}

	for (n=0;n<nproc;n++) {
		if (n != rank && nrecvkey[n] > 0)
			//		MPI_Recv(recvbuff, nrecvkey[n], MPI_BYTE, n,rank, MPI_COMM_WORLD, &status);
			MPI_Recv(recvbuff+nrecdisp[n], nrecvkey[n]*sizeof(MKey), MPI_BYTE, n,rank, MPI_COMM_WORLD, &status);
	}	
	for (n=0;n<nproc;n++) {
		if (n != rank && nsendkey[n] > 0)
			MPI_Wait(&vreq[n], &vstatus[n]);
	}

	MPI_Barrier(MPI_COMM_WORLD);

#endif

	for (m=0; m<local_xsize[0]*local_xsize[1]*local_xsize[2]; m++) {
		data[m] = 0.0;
	}	

	for (n=0; n<tot_recvkey; n++) {

		i = recvbuff[n].x ;
		j = recvbuff[n].y ;
		k = recvbuff[n].z ;

		if (i < 0)
			i += NSIDE ;
		if (i >= NSIDE)

			i -= NSIDE;
		if (j < 0)
			j += NSIDE;
		if (j >= NSIDE)
			j -= NSIDE;
		if (k < 0)
			k += NSIDE;
		if (k >= NSIDE)
			k -= NSIDE;

		i -= local_xstart[0];
		j -= local_xstart[1];
		k -= local_xstart[2];

		data[(i*local_xsize[1]+j)*local_xsize[2]+k] += recvbuff[n].v;

	}


	int nside[3] = {NSIDE, NSIDE, NSIDE};
	double param[2] = {splitRadius, BOXSIZE};

#define NPOW 64

	double power[NPOW];
	int powcnt[NPOW];

	double power_total[NPOW];
	int powcnt_total[NPOW];

	for (n=0; n<NPOW; n++) {
		power[n] = 0.0;
		power_total[n] = 0.0;

		powcnt[n] = 0;
		powcnt_total[n] = 0;
	}

	densitykspace(data, nside, power, powcnt, param);

	MPI_Allreduce(powcnt,powcnt_total,NPOW,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

	MPI_Reduce(power,power_total,NPOW,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

	if ( 0 == PROC_RANK ) {
		FILE *fd;
		char fname[80];

		sprintf(fname, "power_%s.txt", powname);
		fd = fopen(fname, "w");

		printf(" output powerspectrum for '%s' !\n", powname);

		for (n=0; n<NPOW/2; n++) {
			fprintf( fd, "%d %e\n", n, power_total[n]/powcnt_total[n] );
		}

		fclose(fd);
	}

#undef NPOW

	if (mesh != NULL)
		pfree(mesh, 50);

	if (sendbuff != NULL )
		pfree(sendbuff, 55);

	if (nsendkey != NULL)
		pfree(nsendkey, 51);

	if (nsendisp != NULL)
		pfree(nsendisp, 52);

	if (nrecdisp != NULL)
		pfree(nrecdisp, 54);

	if (nrecvkey != NULL )
		pfree(nrecvkey, 53);

	if (recvbuff != NULL)
		pfree(recvbuff, 56);

}



