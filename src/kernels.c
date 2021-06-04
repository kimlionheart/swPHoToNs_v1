#include "photoNs.h"
#include <stdio.h>
#include <stdlib.h>

/*
void p2p_kernel_inleaf(int inode, int jnode) {
	double dx[3], dr, ir3, ir, x2;
	int ip, jp;
	int n,m, idx;
	int start = leaf[inode].ipart;
	int end = leaf[inode].ipart+leaf[inode].npart;
	int endi = end - 1;

	pack2pack_inleaf_count ++;


	leaf[inode].pack2pack_cnt++;

	double rs = splitRadius;

	double coeff = 2.0/sqrt(M_PI);
//printf(" split - -  %lf - \n", splitRadius);
	for (ip=start; ip<endi; ip++)
	{
		for (jp=ip+1; jp<end; jp++)
		{
			dx[0] = part[jp].pos[0] - part[ip].pos[0];
			dx[1] = part[jp].pos[1] - part[ip].pos[1];
			dx[2] = part[jp].pos[2] - part[ip].pos[2];

#ifdef TABLE_GRAVITY
			x2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
			idx = (int) (x2 * invGravFunc);

			if (idx < NumGravFunc && idx > 0)
				ir3 = gravfunc[idx];
			else
				ir3 = 0.0;
#else

			dr = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] );
			if (dr < SoftenScale)
				ir3 = 1.0/(SoftenScale * SoftenScale * SoftenScale);
			else
				ir3 = 1.0/(dr*dr*dr);


#ifdef LONGSHORT
			double drs = 0.5*dr/rs;
			ir3 *= (erfc(drs) + coeff*drs*exp(-drs*drs));
#endif            

#endif // TABLE_GRAVITY

			part[ip].acc[0] +=  part[jp].mass * dx[0] * ir3;
			part[ip].acc[1] +=  part[jp].mass * dx[1] * ir3;
			part[ip].acc[2] +=  part[jp].mass * dx[2] * ir3;
			
			part[jp].acc[0] -=  part[ip].mass * dx[0] * ir3;
			part[jp].acc[1] -=  part[ip].mass * dx[1] * ir3;
			part[jp].acc[2] -=  part[ip].mass * dx[2] * ir3;

			p2p_count_inleaf++;



		}
	}
}

*/

void p2p_kernel(int inode, int jnode) {
	double dx[3], dr, ir3, ir, mp, x2;
	int ip, jp;

	int n,m, idx;

	pack2pack_count++;

	double rs = splitRadius;

	double coeff = 2.0/sqrt(M_PI);



	leaf[inode].pack2pack_cnt++;
//printf(" %d %d \n", inode, jnode);
	for (ip=leaf[inode].ipart; ip<leaf[inode].ipart+leaf[inode].npart; ip++)
	{
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
			//printf(" pp idx = %d(%d)\n", idx, NumGravFunc);
			if (idx < NumGravFunc && idx >=0 )
				ir3 = part[jp].mass * gravfunc[idx];
			else
				ir3 = 0.0;
                    
#else

//			dr = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] + SoftenScale * SoftenScale);
//			ir3 = 1.0/(dr*dr*dr);

			dr = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]) + SoftenScale ;

//			if (dr < SoftenScale)
//				ir3 = part[jp].mass/(SoftenScale*SoftenScale*SoftenScale);
//			else
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


/*

float rinvsqrt( float number )
{
	long i;
	float x2, y;
	const float threehalfs = 1.5F;

	x2 = number * 0.5F;
	y  = number;
	i  = * ( long * ) &y;                      
	i  = 0x5f3759df - ( i >> 1 );               
	y  = * ( float * ) &i;
	y  = y * ( threehalfs - ( x2 * y * y ) ); 
	y  = y * ( threehalfs - ( x2 * y * y ) );   

	return y;
}

void p2p_kernel_rinvsqrt(int inode, int jnode) {
	float dx[3], dr, ir3, ir, mp;
	int ip, jp;

	for (ip=leaf[inode].ipart; ip<leaf[inode].ipart+leaf[inode].npart; ip++)
	{
		for (jp=leaf[jnode].ipart; jp<leaf[jnode].ipart+leaf[jnode].npart; jp++)
		{
			if (jp == ip)
				continue;

			dx[0] = part[jp].pos[0] - part[ip].pos[0];
			dx[1] = part[jp].pos[1] - part[ip].pos[1];
			dx[2] = part[jp].pos[2] - part[ip].pos[2];

			ir = rinvsqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]+SoftSquare);
			ir3 = part[jp].mass*ir*ir*ir;


			part[ip].acc[0] +=  dx[0] * ir3;
			part[ip].acc[1] +=  dx[1] * ir3;
			part[ip].acc[2] +=  dx[2] * ir3;
			//			part[ip].acc[3] -= mp * ir;
		}
	}
}

*/



