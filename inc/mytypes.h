#ifndef MYTYPES_H
#define MYTYPES_H

//#define HYDRODYN


//#define DIRECT_NBODY
//#define NMULTI 4
#define QUADRUPOLE
//#define NMULTI 10
#define OCTUPOLE
#define NMULTI 20
//#define HEXADECAPOLE
//#define NMULTI 35

double masspart;
double IRS;
double COEFF;
double softlen;
double param[4];
#define NSON 2
typedef unsigned long INT64;

typedef struct {
	int  updated;
	int  changed;
	int  npart;
	int  son[NSON];
	double width[3];
	double split;
	//	double tmp[2];
	double center[3];
	double M[NMULTI];
	double L[NMULTI];
	int m2l_count;
	int m2l_disp;
} Node;

typedef struct {
	int  updated;
	int  changed;
	int npart;
	int ipart;
	double width[3];
	//	double tmp[2];
	double center[3];
	double M[NMULTI];
	double L[NMULTI];
	int pack2pack_count;
	int m2l_count;
	int m2l_disp;
	int pack2pack_cnt;
} Pack;
/*
typedef struct {
	double pos[3];
	double replenish1;
	double acc[3];
	double replenish2;
	double vel[3];
	double replenish3;
    	double acc_pm[3];
	double replenish4;
} Body;
*/


typedef struct {
	double pos[3];
	double replenish1;
	double acc[3];
	double replenish2;
	double vel[3];
	double replenish3;
    	double acc_pm[3];
	double replenish4;
#ifdef TESTGRAV
	double acc_split[3];
#endif

#ifdef TREE_CODE
	double acc_tree[4];
#endif

#ifdef DIRECT_NBODY
	double acc_direct[4];
#endif

	double mass;
	double len;
    	long id;
	int active;

} Body;

Body *part;

Node *btree;
Pack *leaf;

double my_gb;

#endif /// MYTYPES_H ////


