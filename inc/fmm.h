#ifndef FMM_H
#define FMM_H

#include "kernels.h"
#include "operator.h"
#include "utility.h"
#include "toptree.h"
#include "remotes.h"
#include "adaptive.h"


void check_gravity() ;
//void build_localtree(int NPART, double bdleft[], double bdright[], int max_leaf, double open, int direct);
void build_localtree();
void fmm_solver() ;
int intergra();
void fmm_construct( ) ;
void fmm_deconstruct( ) ;


void fmm_prepare();
void fmm_task();
void fmm_ext();

void fmm_solver_total();


#endif

