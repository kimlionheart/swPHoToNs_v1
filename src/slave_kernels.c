#include <slave.h>
//#include <athread.h>
#include <dma.h>
#include <simd.h>
#include <math.h>
#include "math_dldm.h"
#define CPE
#define LWPF_KERNELS K(A) K(B) K(C)
#define LWPF_UNIT U(TEST)
#include "lwpf2.h"
#define MAX_PACKAGE 16
#define slave_task_alone  256
#define FAST_ERF
#ifndef _ADMA_H_INCLUDE
#define _ADMA_H_INCLUDE

#define DMA_GET(da,mode,src,dest,len,re_addr,n)    \
	({   				\
	 dma_set_op(&da, DMA_GET);	\
	 dma_set_mode(&da, mode);	\
	 dma_set_size(&da, len);	\
	 dma_set_reply(&da, re_addr);	\
	 dma(da, src, dest);		\
	 dma_wait((re_addr), n);	\
	 	 })				

#define DMA_PUT(da,mode,src,dest,len,re_addr,n)	   \
	({				\
	 dma_set_op(&da, DMA_PUT);	\
	 dma_set_mode(&da, mode);	\
	 dma_set_size(&da, len);	\
	 dma_set_reply(&da, re_addr);	\
	 dma(da, src, dest);		\
	 dma_wait((re_addr), n);	\
	 	 })

#endif
typedef struct {
	double pos[3];
	double replenish1;
	double acc[3];
	double replenish2;
	double vel[3];
	double replenish3;
	double acc_pm[3];
	double replenish4;
	double mass;
	double len;
	long id;
	int active;

} Body;
typedef struct {
	double pos[3];
	double mass;
} RemoteBody;

typedef struct{
	RemoteBody *slave_remote;
	Body *slave_body;
	double rs;
	double coeff;
	double SoftenScale;
	double mass;
	int  *s_i_part;
	int  *s_j_part;
	int idxp2p;
}slave_data;
void p2p_kernel_slave(slave_data *slave_all){
	//lwpf_enter(TEST);
	slave_data total;
	math_init();
	total.idxp2p = 0;
//	int xyx[128][2];
	int i, k, j, ip, jp, l, num = 0, count = 0;;
	double dx[4], dr, ir3, ir, mp, x2;
	int slave_task_num = 0;
	dma_desc get_desc = 0;
	dma_desc put_desc = 0;
	dma_desc get_leaf = 0;
	dma_desc get_part = 0;
	volatile unsigned long get_reply = 0;
	//volatile unsigned long put_reply = 0 ;
	volatile int put_reply = 0;
	int my_id;
	my_id = athread_get_id(-1);
	//if(my_id == 0){
	get_reply = 0;
	athread_get(PE_MODE, slave_all, &total, sizeof(slave_data), &get_reply , 0, 0, 0);
//	 while (get_reply != 1);
	while (!wait_get(get_reply));
	double rs = total.rs;
	double coeff = total.coeff;
	double SoftenScale = total.SoftenScale;

	double rs2 = rs*2;
	double mass = total.mass;
	slave_task_num = total.idxp2p / 64;
	int slave_limit = total.idxp2p%64;
	
//	if(my_id<slave_limit) {
//		slave_task_num++;
//	}
	
	
	int slave_base = 0;
	
	int part[256][4];	
	doublev4 node_i_data_v4[16][5];
	double(*node_i_data)[5][4] = node_i_data_v4;
	doublev4 node_j_data_v4[16][5];
	double(*node_j_data)[5][4] = node_j_data_v4;

doublev4 ir3_v4;
	doublev4 dx_v4;
	doublev4 dx1_v4;
	doublev4 dx2_v4;
	doublev4 dx3_v4;
	doublev4 d2x_v4;
	doublev4 d2x1_v4;
	doublev4 d2x2_v4;
	doublev4 d2x3_v4;
	doublev4 dr_v4;
	doublev4 p12l,p12h,p34l,p34h;
	doublev4 half =simd_set_doublev4(0.5,0.5,0.5,0.5);
	doublev4 two = simd_set_doublev4(-2.0,-2.0,-2.0,-2.0);
	doublev4 one = simd_set_doublev4(-1.0,-1.0,-1.0,-1.0);
	doublev4 zone = simd_set_doublev4(1.0,1.0,1.0,1.0);
	doublev4 ztwo = simd_set_doublev4(2.0,2.0,2.0,2.0);
	doublev4 zero = simd_set_doublev4(0.0,0.0,0.0,0.0);
	doublev4 SoftenScale_v4 = simd_set_doublev4(SoftenScale,SoftenScale,SoftenScale,SoftenScale);//转入
	doublev4 coeff_v4 = simd_set_doublev4(coeff,coeff,coeff,coeff);
	doublev4 rs2_v4 = simd_set_doublev4(rs2,rs2,rs2,rs2);
	doublev4 mass_v4 = simd_set_doublev4(mass,mass,mass,mass);
doublev4 invln2 = simd_set_doublev4(1.44269504088896338700e+00,1.44269504088896338700e+00,1.44269504088896338700e+00,1.44269504088896338700e+00);
doublev4 ln2hi = simd_set_doublev4(6.93147180369123816490e-01,6.93147180369123816490e-01,6.93147180369123816490e-01,6.93147180369123816490e-01);
doublev4 ln2lo = simd_set_doublev4(1.90821492927058770002e-10,1.90821492927058770002e-10,1.90821492927058770002e-10,1.90821492927058770002e-10);
doublev4 p1 = simd_set_doublev4(1.66666666666666019037e-01,1.66666666666666019037e-01,1.66666666666666019037e-01,1.66666666666666019037e-01);
doublev4 p2 = simd_set_doublev4(-2.77777777770155933842e-03,-2.77777777770155933842e-03,-2.77777777770155933842e-03,-2.77777777770155933842e-03);
doublev4 p3 = simd_set_doublev4(6.61375632143793436117e-05,6.61375632143793436117e-05,6.61375632143793436117e-05,6.61375632143793436117e-05);
doublev4 p4 = simd_set_doublev4(-1.65339022054652515390e-06,-1.65339022054652515390e-06,-1.65339022054652515390e-06,-1.65339022054652515390e-06);
doublev4 p5 = simd_set_doublev4(4.13813679705723846039e-08,4.13813679705723846039e-08,4.13813679705723846039e-08,4.13813679705723846039e-08);
  

doublev4 a0 = simd_set_doublev4(0.3275911,0.3275911,0.3275911,0.3275911);
doublev4 a1 = simd_set_doublev4(0.254829592,0.254829592,0.254829592,0.254829592);
doublev4 a2 = simd_set_doublev4(-0.284496736,-0.284496736,-0.284496736,-0.284496736);
doublev4 a3 = simd_set_doublev4(1.421413741,1.421413741,1.421413741,1.421413741);
doublev4 a4 = simd_set_doublev4(-1.453152027,-1.453152027,-1.453152027,-1.453152027);
doublev4 a5 = simd_set_doublev4(1.061405429,1.061405429,1.061405429,1.061405429);
doublev4 t1,t2,t3,t4,t5;

union {double f; uint64_t i;} u1,u2,u3,u4;
    double ln2h[4];
    doublev4 ln2v4;
    int n;
    doublev4 uv4;
    doublev4 hi;
    doublev4 lo;
    doublev4 x_v4;
    doublev4 xx_v4;
    doublev4 c;
    doublev4 y;
      doublev4 k_v4 ;
	doublev4 drs_v4;
	doublev4 drs2_v4;
	doublev4 dxx_v4;
	double drs_erfc[4];
	double drs_exp[4];
	double drs2[4];
	double drs;
	double drs4[4];
	double dxx[4];
	double ir34[4];
	doublev4 dxx1_v4,dxx2_v4,dxx3_v4,dxx4_v4;
	doublev4 drs_erfc_v4;
	doublev4 drs_exp_v4;
	int ip1,ip2,ip3,ip4;
int jp1,jp2,jp3,jp4;
int parti; 
int partj; 
int numm;


	get_reply =0;
	 athread_get(PE_MODE, total.s_j_part+my_id*slave_task_alone * 4,part,sizeof(int)*(slave_task_alone)*4, &get_reply,0,0,0);
	 while (get_reply != 1);


	for (i = 0; i<slave_task_num; i++)
	{

		get_reply = 0;
	athread_get(PE_MODE,total.slave_body+part[i][1],&node_i_data_v4[0][0],sizeof(Body)*part[i][0], &get_reply , 0, 0, 0);
	while (get_reply != 1);
	get_reply = 0;
	athread_get(PE_MODE,total.slave_body+part[i][3],&node_j_data_v4[0][0],sizeof(Body)*part[i][2], &get_reply , 0, 0, 0);
	while (get_reply != 1);

	
		parti = part[i][0]-1;
			partj = part[i][2]-1;
			//limit = partj - 3;
			count = part[i][0]*part[i][2];
			count = count -3;
			ip4=0;
			jp4 = -1;
			//if(my_id == 0) printf("xxx%d  %d\n",parti,partj);
			for(numm = 0;numm < count;numm+=4){
				ip1=ip4;jp1=jp4+1;if(jp1>partj) {ip1++;jp1=0;}
				ip2 = ip1;jp2 = jp1+1;if(jp2>partj) {ip2++;jp2=0;}
				ip3 = ip2;jp3 = jp2+1;if(jp3>partj) {ip3++;jp3=0;}
				ip4 = ip3;jp4 = jp3+1;if(jp4>partj) {ip4++;jp4=0;}
					
					dx_v4 = simd_vsubd(node_j_data_v4[jp1][0], node_i_data_v4[ip1][0]);
					dx1_v4 = simd_vsubd(node_j_data_v4[jp2][0], node_i_data_v4[ip2][0]);
					dx2_v4 = simd_vsubd(node_j_data_v4[jp3][0], node_i_data_v4[ip3][0]);
					dx3_v4 = simd_vsubd(node_j_data_v4[jp4][0], node_i_data_v4[ip4][0]);
					
					d2x_v4 =dx_v4*dx_v4;
					d2x1_v4 =dx1_v4*dx1_v4;
					d2x2_v4 =dx2_v4*dx2_v4;
					d2x3_v4 =dx3_v4*dx3_v4;
					p12l=simd_vshff(d2x1_v4,d2x_v4,0x44);
					p12h=simd_vshff(d2x1_v4,d2x_v4,0xee);
					p34l=simd_vshff(d2x3_v4,d2x2_v4,0x44);
					p34h=simd_vshff(d2x3_v4,d2x2_v4,0xee);
					d2x_v4 = simd_vshff(p34l,p12l,0x88);
					d2x1_v4=simd_vshff(p34l,p12l,0xdd);
					d2x2_v4=simd_vshff(p34h,p12h,0x88);//转置
					dr_v4 = d2x_v4+d2x1_v4+d2x2_v4;
					dr_v4 = simd_vsqrtd(dr_v4)+SoftenScale_v4;

					ir3_v4 = mass_v4/(dr_v4*dr_v4*dr_v4);//质量相除；;
					drs_v4 = dr_v4 / rs2_v4;
					//拆分，分别求解
//					simd_store(drs_v4,&drs4);
//					drs_erfc[0] = DERF(drs4[0]);
//					drs_erfc[1] = DERF(drs4[1]);
//					drs_erfc[2] = DERF(drs4[2]);
//					drs_erfc[3] = DERF(drs4[3]);//尝试更换库
//					simd_load(drs_erfc_v4,&drs_erfc);//向量装入
					
					drs2_v4 = zero - drs_v4*drs_v4;
			//		simd_store(drs2_v4, &drs2);//先乘在0-在拆分
					ln2v4 = invln2 *drs2_v4+half;
					simd_store(ln2v4,&ln2h);
				    ln2h[0] = floor(ln2h[0]);
				    ln2h[1] = floor(ln2h[1]);
					ln2h[2] = floor(ln2h[2]);
					ln2h[3] = floor(ln2h[3]);
					simd_load(k_v4,&ln2h);
					n = (int)(ln2h[0]);
					u1.i = (uint64_t)(0x3ff+n)<<52;
					n = (int)(ln2h[1]);
					u2.i = (uint64_t)(0x3ff+n)<<52;
					n = (int)(ln2h[2]);
					u3.i = (uint64_t)(0x3ff+n)<<52;
					n = (int)(ln2h[3]);
					u4.i = (uint64_t)(0x3ff+n)<<52;
					uv4 = simd_set_doublev4(u1.f,u2.f,u3.f,u4.f);
					hi = drs2_v4 - k_v4*ln2hi;
					lo = k_v4*ln2lo;
					x_v4 = hi - lo;
					xx_v4 = x_v4*x_v4;
					c = x_v4 - xx_v4*(p1+xx_v4*(p2+xx_v4*(p3+xx_v4*(p4+xx_v4*p5))));
					y = zone + (x_v4*c/(ztwo-c) - lo + hi);
				
					drs_exp_v4 = y*uv4;
				
//				simd_store(drs_exp_v4,&dxx);
//					simd_store(drs2_v4,&drs2);
//					drs_exp[0] = DEXP(drs2[0]);
//					drs_exp[1] = DEXP(drs2[1]);
//					drs_exp[2] = DEXP(drs2[2]);
//					drs_exp[3] = DEXP(drs2[3]);
//					simd_load(drs_exp_v4,&drs_exp);
						   
					
//if(fabs(drs_exp[0]-dxx[0])>1e-10) printf("   %.17lf  %.17lf   %.17lf   \n",drs_exp[0],dxx[0],drs2[0]);
//if(fabs(drs_exp[1]-dxx[1])>1e-10) printf("   %d  %.17lf   %.17lf   \n",my_id,drs_exp[1]-dxx[1],drs2[1]);
//if(fabs(drs_exp[2]-dxx[2])>1e-10) printf("   %d  %.17lf   %.17lf   \n",my_id,drs_exp[2]-dxx[2],drs2[2]);
//if(fabs(drs_exp[3]-dxx[3])>1e-10) printf("   %d  %.17lf   %.17lf   \n",my_id,drs_exp[3]-dxx[3],drs2[3]);
#ifdef FAST_ERF
				t1 = zone/(zone+a0*drs_v4);
				t2 = t1*t1;
				t3 = t2*t1;
				t4 = t3*t1;
				t5 = t4*t1;
				drs_erfc_v4 = zone-(a1*t1+a2*t2+a3*t3+a4*t4+a5*t5)*drs_exp_v4;
#else
					simd_store(drs_v4,&drs4);
					drs_erfc[0] = DERF(drs4[0]);
					drs_erfc[1] = DERF(drs4[1]);
					drs_erfc[2] = DERF(drs4[2]);
					drs_erfc[3] = DERF(drs4[3]);//尝试更换库
					simd_load(drs_erfc_v4,&drs_erfc);//向量装入
#endif

				ir3_v4 = ir3_v4 * (zone - drs_erfc_v4 + coeff_v4*drs_v4*drs_exp_v4);
			

					//dxx_v4 = dx_v4*ir3_v4;
					//拆解和装入
					simd_store(ir3_v4,&ir34);
					dxx1_v4 = simd_set_doublev4(ir34[0],ir34[0],ir34[0],0);
					dxx2_v4 = simd_set_doublev4(ir34[1],ir34[1],ir34[1],0);
					dxx3_v4 = simd_set_doublev4(ir34[2],ir34[2],ir34[2],0);
					dxx4_v4 = simd_set_doublev4(ir34[3],ir34[3],ir34[3],0);

					node_i_data_v4[ip1][1] = dxx1_v4*dx_v4 +node_i_data_v4[ip1][1];
					node_i_data_v4[ip2][1] = dxx2_v4*dx1_v4 +node_i_data_v4[ip2][1];
					node_i_data_v4[ip3][1] = dxx3_v4*dx2_v4 +node_i_data_v4[ip3][1];
				 	node_i_data_v4[ip4][1] = dxx4_v4*dx3_v4 +node_i_data_v4[ip4][1];
			}
		
			if (ip4 == parti){
			jp4++;
				for(;jp4<=partj;jp4++){
//			if(my_id == 0) printf("xxx\n");
					dx_v4 = simd_vsubd(node_j_data_v4[jp4][0], node_i_data_v4[ip4][0]);
					simd_store(dx_v4, &dx);
					if (dx[0] == 0 && dx[1] == 0 && dx[2] == 0)  continue;
					dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]) + SoftenScale;
					ir3 = mass / (dr*dr*dr);

					drs = dr / rs2;
					ir3 = ir3 * (1- DERF(drs) + coeff*drs*DEXP(-drs*drs));
					ir3_v4 = simd_set_doublev4(ir3, ir3, ir3, 0);


					dxx_v4 = dx_v4*ir3_v4;
					node_i_data_v4[ip4][1] = dxx_v4 + node_i_data_v4[ip4][1];
				}
			}
			else {
			jp4++;
				for(;jp4<=partj;jp4++){
//			if(my_id == 0) printf("xxx\n");
					dx_v4 = simd_vsubd(node_j_data_v4[jp4][0], node_i_data_v4[ip4][0]);
					simd_store(dx_v4, &dx);
					if (dx[0] == 0 && dx[1] == 0 && dx[2] == 0)  continue;
					dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]) + SoftenScale;
					ir3 = mass / (dr*dr*dr);

					drs = dr / rs2;
					ir3 = ir3 * (1- DERF(drs) + coeff*drs*DEXP(-drs*drs));
					ir3_v4 = simd_set_doublev4(ir3, ir3, ir3, 0);


					dxx_v4 = dx_v4*ir3_v4;
					node_i_data_v4[ip4][1] = dxx_v4 + node_i_data_v4[ip4][1];
				}
			ip4++;
				for(;ip4<=parti;ip4++)
					for(jp4 = 0;jp4<=partj;jp4++){
					dx_v4 = simd_vsubd(node_j_data_v4[jp4][0], node_i_data_v4[ip4][0]);
					simd_store(dx_v4, &dx);
					if (dx[0] == 0 && dx[1] == 0 && dx[2] == 0)  continue;
					dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]) + SoftenScale;
					ir3 = mass / (dr*dr*dr);

					drs = dr / rs2;
					ir3 = ir3 * (1- DERF(drs) + coeff*drs*DEXP(-drs*drs));
					ir3_v4 = simd_set_doublev4(ir3, ir3, ir3, 0);


					dxx_v4 = dx_v4*ir3_v4;
					node_i_data_v4[ip4][1] = dxx_v4 + node_i_data_v4[ip4][1];
				}
			}
		put_reply = 0;

		athread_put(PE_MODE, &node_i_data_v4[0][0], total.slave_body + part[i][1], sizeof(Body)*part[i][0], &put_reply, 0, 0);
		while (put_reply != 1);


	}
}

void p2p_kernel_ex_slave(slave_data *slave_all){
	slave_data total;
	math_init();
	int my_id;
	volatile unsigned long get_reply = 0;
	volatile unsigned long get_reply_ex = 0;
	dma_desc get_desc = 0;
	get_reply=0;
	DMA_GET(get_desc,PE_MODE, slave_all, &total, sizeof(slave_data), &get_reply,1);
//	while(!wait_get(get_reply));
	dma_wait((&get_reply),1);
	int i,k,j,ip,jp,l,num=0,count=0;
	double dx[4], dr, ir3, ir, mp, x2;
	double rs = total.rs;
	double coeff = total.coeff;
	double mass = total.mass;
	double SoftenScale = total.SoftenScale;	
	double rs2 = rs*2;
	int slave_task_num = 0;
	dma_desc put_desc = 0;
	dma_desc get_leaf = 0;
	dma_desc get_part = 0;
	dma_desc get_remote = 0;
	volatile unsigned long put_reply=0;
	my_id = athread_get_id(-1);
	slave_task_num = total.idxp2p/64;
	int slave_limit = total.idxp2p%64;
	
	if(my_id<slave_limit) {
		slave_task_num++;

	}
	int part[256][4];
	doublev4 node_i_data_v4[16][5];
	double (*node_i_data)[5][4]=node_i_data_v4;
	doublev4 node_j_data_v4[16];
	double (*node_j_data)[4]=node_j_data_v4;
	
	doublev4 ir3_v4;
	doublev4 dx_v4;
	doublev4 dx1_v4;
	doublev4 dx2_v4;
	doublev4 dx3_v4;
	doublev4 d2x_v4;
	doublev4 d2x1_v4;
	doublev4 d2x2_v4;
	doublev4 d2x3_v4;
	doublev4 dr_v4;
	doublev4 p12l,p12h,p34l,p34h;
	doublev4 half =simd_set_doublev4(0.5,0.5,0.5,0.5);
	doublev4 two = simd_set_doublev4(-2.0,-2.0,-2.0,-2.0);
	doublev4 one = simd_set_doublev4(-1.0,-1.0,-1.0,-1.0);
	doublev4 zone = simd_set_doublev4(1.0,1.0,1.0,1.0);
	doublev4 ztwo = simd_set_doublev4(2.0,2.0,2.0,2.0);
	doublev4 zero = simd_set_doublev4(0.0,0.0,0.0,0.0);
	doublev4 SoftenScale_v4 = simd_set_doublev4(SoftenScale,SoftenScale,SoftenScale,SoftenScale);//转入
	doublev4 coeff_v4 = simd_set_doublev4(coeff,coeff,coeff,coeff);
	doublev4 rs2_v4 = simd_set_doublev4(rs2,rs2,rs2,rs2);
	doublev4 mass_v4 = simd_set_doublev4(mass,mass,mass,mass);
doublev4 invln2 = simd_set_doublev4(1.44269504088896338700e+00,1.44269504088896338700e+00,1.44269504088896338700e+00,1.44269504088896338700e+00);
doublev4 ln2hi = simd_set_doublev4(6.93147180369123816490e-01,6.93147180369123816490e-01,6.93147180369123816490e-01,6.93147180369123816490e-01);
doublev4 ln2lo = simd_set_doublev4(1.90821492927058770002e-10,1.90821492927058770002e-10,1.90821492927058770002e-10,1.90821492927058770002e-10);
doublev4 p1 = simd_set_doublev4(1.66666666666666019037e-01,1.66666666666666019037e-01,1.66666666666666019037e-01,1.66666666666666019037e-01);
doublev4 p2 = simd_set_doublev4(-2.77777777770155933842e-03,-2.77777777770155933842e-03,-2.77777777770155933842e-03,-2.77777777770155933842e-03);
doublev4 p3 = simd_set_doublev4(6.61375632143793436117e-05,6.61375632143793436117e-05,6.61375632143793436117e-05,6.61375632143793436117e-05);
doublev4 p4 = simd_set_doublev4(-1.65339022054652515390e-06,-1.65339022054652515390e-06,-1.65339022054652515390e-06,-1.65339022054652515390e-06);
doublev4 p5 = simd_set_doublev4(4.13813679705723846039e-08,4.13813679705723846039e-08,4.13813679705723846039e-08,4.13813679705723846039e-08);
doublev4 a0 = simd_set_doublev4(0.3275911,0.3275911,0.3275911,0.3275911);
doublev4 a1 = simd_set_doublev4(0.254829592,0.254829592,0.254829592,0.254829592);
doublev4 a2 = simd_set_doublev4(-0.284496736,-0.284496736,-0.284496736,-0.284496736);
doublev4 a3 = simd_set_doublev4(1.421413741,1.421413741,1.421413741,1.421413741);
doublev4 a4 = simd_set_doublev4(-1.453152027,-1.453152027,-1.453152027,-1.453152027);
doublev4 a5 = simd_set_doublev4(1.061405429,1.061405429,1.061405429,1.061405429);
doublev4 t1,t2,t3,t4,t5;



union {double f; uint64_t i;} u1,u2,u3,u4;
    double ln2h[4];
    doublev4 ln2v4;
    int n;
    doublev4 uv4;
    doublev4 hi;
    doublev4 lo;
    doublev4 x_v4;
    doublev4 xx_v4;
    doublev4 c;
    doublev4 y;
      doublev4 k_v4 ;
	doublev4 drs_v4;
	doublev4 drs2_v4;
	doublev4 dxx_v4;
	double drs_erfc[4];
	double drs_exp[4];
//	double drs2[4];
	double drs;
	double drs4[4];
	double dxx[4];
	double ir34[4];
	doublev4 dxx1_v4,dxx2_v4,dxx3_v4,dxx4_v4;
	doublev4 drs_erfc_v4;
	doublev4 drs_exp_v4;
	int ip1,ip2,ip3,ip4;
int jp1,jp2,jp3,jp4;
int parti; 
int partj; 
int numm;
	


	get_reply =0;
	 athread_get(PE_MODE, total.s_j_part+my_id*slave_task_alone * 4,part,sizeof(int)*(slave_task_alone)*4, &get_reply,0,0,0);
	 while (get_reply != 1);
	 
	for(i = 0;i<slave_task_num;i++)
	{

get_reply = 0;
athread_get(PE_MODE,total.slave_body+part[i][1],node_i_data,sizeof(Body)*part[i][0], &get_reply,0,0,0);	
while (get_reply != 1);

get_reply_ex = 0;
athread_get(PE_MODE,total.slave_remote+part[i][2],node_j_data, sizeof(RemoteBody)*part[i][3],&get_reply,0,0,0);
while (get_reply != 1);
	parti = part[i][0]-1;
			partj = part[i][3]-1;
			//limit = partj - 3;
			count = part[i][0]*part[i][3];
			count = count -3;
			ip4=0;
			jp4 = -1;
			//if(my_id == 0) printf("xxx%d  %d\n",parti,partj);
			for(numm = 0;numm < count;numm+=4){
				ip1=ip4;jp1=jp4+1;if(jp1>partj) {ip1++;jp1=0;}
				ip2 = ip1;jp2 = jp1+1;if(jp2>partj) {ip2++;jp2=0;}
				ip3 = ip2;jp3 = jp2+1;if(jp3>partj) {ip3++;jp3=0;}
				ip4 = ip3;jp4 = jp3+1;if(jp4>partj) {ip4++;jp4=0;}
					
					dx_v4  = simd_vsubd(node_j_data_v4[jp1], node_i_data_v4[ip1][0]);
				        dx1_v4 = simd_vsubd(node_j_data_v4[jp2], node_i_data_v4[ip2][0]);
					dx2_v4 = simd_vsubd(node_j_data_v4[jp3], node_i_data_v4[ip3][0]);
					dx3_v4 = simd_vsubd(node_j_data_v4[jp4], node_i_data_v4[ip4][0]);
					
					d2x_v4 =dx_v4*dx_v4;
					d2x1_v4 =dx1_v4*dx1_v4;
					d2x2_v4 =dx2_v4*dx2_v4;
					d2x3_v4 =dx3_v4*dx3_v4;
					p12l=simd_vshff(d2x1_v4,d2x_v4,0x44);
					p12h=simd_vshff(d2x1_v4,d2x_v4,0xee);
					p34l=simd_vshff(d2x3_v4,d2x2_v4,0x44);
					p34h=simd_vshff(d2x3_v4,d2x2_v4,0xee);
					d2x_v4 = simd_vshff(p34l,p12l,0x88);
					d2x1_v4=simd_vshff(p34l,p12l,0xdd);
					d2x2_v4=simd_vshff(p34h,p12h,0x88);//转置
					dr_v4 = d2x_v4+d2x1_v4+d2x2_v4;
					dr_v4 = simd_vsqrtd(dr_v4)+SoftenScale_v4;

					ir3_v4 = mass_v4/(dr_v4*dr_v4*dr_v4);//质量相除；;
					drs_v4 = dr_v4 / rs2_v4;
					//拆分，分别求解
	//				simd_store(drs_v4,&drs4);
	//				drs_erfc[0] = DERF(drs4[0]);
	//				drs_erfc[1] = DERF(drs4[1]);
	//				drs_erfc[2] = DERF(drs4[2]);
	//				drs_erfc[3] = DERF(drs4[3]);//尝试更换库
	//				simd_load(drs_erfc_v4,&drs_erfc);//向量装入
					
					drs2_v4 = zero - drs_v4*drs_v4;
//				simd_store(drs2_v4, &drs2);//先乘在0-在拆分*					ln2v4 = invln2 *drs2_v4+half;
					simd_store(ln2v4,&ln2h);
				    ln2h[0] = floor(ln2h[0]);
				    ln2h[1] = floor(ln2h[1]);
					ln2h[2] = floor(ln2h[2]);
					ln2h[3] = floor(ln2h[3]);
					simd_load(k_v4,&ln2h);
					n = (int)(ln2h[0]);
					u1.i = (uint64_t)(0x3ff+n)<<52;
					n = (int)(ln2h[1]);
					u2.i = (uint64_t)(0x3ff+n)<<52;
					n = (int)(ln2h[2]);
					u3.i = (uint64_t)(0x3ff+n)<<52;
					n = (int)(ln2h[3]);
					u4.i = (uint64_t)(0x3ff+n)<<52;
					uv4 = simd_set_doublev4(u1.f,u2.f,u3.f,u4.f);
					hi = drs2_v4 - k_v4*ln2hi;
					lo = k_v4*ln2lo;
					x_v4 = hi - lo;
					xx_v4 = x_v4*x_v4;
					c = x_v4 - xx_v4*(p1+xx_v4*(p2+xx_v4*(p3+xx_v4*(p4+xx_v4*p5))));
					y = zone + (x_v4*c/(ztwo-c) - lo + hi);
					drs_exp_v4 = y*uv4;
				//simd_store(drs_exp_v4,&dxx);
		      			   
//				simd_store(drs2_v4,&drs4);
//				drs_exp[0] = DEXP(drs4[0]);
//				drs_exp[1] = DEXP(drs4[1]);
//				drs_exp[2] = DEXP(drs4[2]);
//				drs_exp[3] = DEXP(drs4[3]);
//				simd_load(drs_exp_v4,&drs_exp);
#ifdef FAST_ERF
				t1 = zone/(zone+a0*drs_v4);
				t2 = t1*t1;
				t3 = t2*t1;
				t4 = t3*t1;
				t5 = t4*t1;
				drs_erfc_v4 = zone-(a1*t1+a2*t2+a3*t3+a4*t4+a5*t5)*drs_exp_v4;
#else
					simd_store(drs_v4,&drs4);
					drs_erfc[0] = DERF(drs4[0]);
					drs_erfc[1] = DERF(drs4[1]);
					drs_erfc[2] = DERF(drs4[2]);
					drs_erfc[3] = DERF(drs4[3]);//尝试更换库
					simd_load(drs_erfc_v4,&drs_erfc);//向量装入
#endif

					ir3_v4 = ir3_v4 * (zone - drs_erfc_v4 + coeff_v4*drs_v4*drs_exp_v4);
			

					//dxx_v4 = dx_v4*ir3_v4;
					//拆解和装入
					simd_store(ir3_v4,&ir34);
					dxx1_v4 = simd_set_doublev4(ir34[0],ir34[0],ir34[0],0);
					dxx2_v4 = simd_set_doublev4(ir34[1],ir34[1],ir34[1],0);
					dxx3_v4 = simd_set_doublev4(ir34[2],ir34[2],ir34[2],0);
					dxx4_v4 = simd_set_doublev4(ir34[3],ir34[3],ir34[3],0);

					node_i_data_v4[ip1][1] = dxx1_v4*dx_v4 +node_i_data_v4[ip1][1];
					node_i_data_v4[ip2][1] = dxx2_v4*dx1_v4 +node_i_data_v4[ip2][1];
					node_i_data_v4[ip3][1] = dxx3_v4*dx2_v4 +node_i_data_v4[ip3][1];
				 	node_i_data_v4[ip4][1] = dxx4_v4*dx3_v4 +node_i_data_v4[ip4][1];
			}
		
			if (ip4 == parti){
			jp4++;
				for(;jp4<=partj;jp4++){
//			if(my_id == 0) printf("xxx\n");
					dx_v4 = simd_vsubd(node_j_data_v4[jp4], node_i_data_v4[ip4][0]);
					simd_store(dx_v4, &dx);
					if (dx[0] == 0 && dx[1] == 0 && dx[2] == 0)  continue;
					dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]) + SoftenScale;
					ir3 = mass / (dr*dr*dr);

					drs = dr / rs2;
					ir3 = ir3 * (1- DERF(drs) + coeff*drs*DEXP(-drs*drs));
					ir3_v4 = simd_set_doublev4(ir3, ir3, ir3, 0);


					dxx_v4 = dx_v4*ir3_v4;
					node_i_data_v4[ip4][1] = dxx_v4 + node_i_data_v4[ip4][1];
				}
			}
			else {
			jp4++;
				for(;jp4<=partj;jp4++){
//			if(my_id == 0) printf("xxx\n");
					dx_v4 = simd_vsubd(node_j_data_v4[jp4], node_i_data_v4[ip4][0]);
					simd_store(dx_v4, &dx);
					if (dx[0] == 0 && dx[1] == 0 && dx[2] == 0)  continue;
					dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]) + SoftenScale;
					ir3 = mass / (dr*dr*dr);

					drs = dr / rs2;
					ir3 = ir3 * (1- DERF(drs) + coeff*drs*DEXP(-drs*drs));
					ir3_v4 = simd_set_doublev4(ir3, ir3, ir3, 0);


					dxx_v4 = dx_v4*ir3_v4;
					node_i_data_v4[ip4][1] = dxx_v4 + node_i_data_v4[ip4][1];
				}
			ip4++;
				for(;ip4<=parti;ip4++)
					for(jp4 = 0;jp4<=partj;jp4++){
//			if(my_id == 0) printf("xxx\n");
					dx_v4 = simd_vsubd(node_j_data_v4[jp4], node_i_data_v4[ip4][0]);
					simd_store(dx_v4, &dx);
					if (dx[0] == 0 && dx[1] == 0 && dx[2] == 0)  continue;
					dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]) + SoftenScale;
					ir3 = mass / (dr*dr*dr);

					drs = dr / rs2;
					ir3 = ir3 * (1- DERF(drs) + coeff*drs*DEXP(-drs*drs));
					ir3_v4 = simd_set_doublev4(ir3, ir3, ir3, 0);


					dxx_v4 = dx_v4*ir3_v4;
					node_i_data_v4[ip4][1] = dxx_v4 + node_i_data_v4[ip4][1];
				}
			}

	put_reply = 0;
	athread_put(PE_MODE,node_i_data,total.slave_body+part[i][1],sizeof(Body)*part[i][0],&put_reply, 0, 0);
	while (put_reply != 1);

	}
	
	
}
