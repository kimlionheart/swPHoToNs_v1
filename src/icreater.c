/* use Zel'dovich approximation to generate the IC, modified from 2LPTic */


#include "photoNs.h"
#include "initial.h"
#include "icreater.h"
#include <math.h>

double PrimordialIndex = 1.0;
double InputSpectrum_UnitLength_in_cm =  3.085678e21;
double	UnitLength_in_cm = 3.085678e21;

typedef double FLOAT;

int direct_ic;
long seedi;

void rand3_(FLOAT* rn) {
	*rn = ran3(&seedi);
}

double Dplus, Norm;
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran31(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1; i<=54; i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1; k<=4; k++)
			for (i=1; i<=55; i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC



#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran32(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1; i<=54; i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1; k<=4; k++)
			for (i=1; i<=55; i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran33(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	long mj,mk;
	int i,ii,k;
	static int iff=0;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1; i<=54; i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1; k<=4; k++)
			for (i=1; i<=55; i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
//////// NR subrutine //////
void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=(double *)malloc((size_t) ((n+1)*sizeof(double)));
	d=(double *)malloc((size_t) ((n+1)*sizeof(double)));

	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) {printf("Error in routine polint");exit(111);}
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free(d);
	free(c);
}

#define FUNC(x) ((*func)(x))
double trapzd(double (*func)(double), double a, double b, int n)
{
	double x,tnm,sum,del;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC

#define EPS 1.0e-8
#define JMAX 40
#define JMAXP (JMAX+1)
#define K 5

double qromb(double (*func)(double), double a, double b)
{
	void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
	double trapzd(double (*func)(double), double a, double b, int n);
	//        void nrerror(char error_text[]);
	double ss,dss;
	double s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	printf("Too many steps in routine qromb\n");
	exit(222);
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K
//////// NR subrutine //////


static  int NPowerTable;

static struct pow_table
{
	double logk, logD;
}
*PowerTable;


//////// 2LPTic ///////
double growth_int(double a) {
	return pow(a / (OmegaM0 + (1 - OmegaM0 - OmegaX0) * a + OmegaX0 * a * a * a), 1.5);
}

double growth(double a) {
	double hubble_a;

	hubble_a = sqrt(OmegaM0 / (a * a * a) + (1 - OmegaM0 - OmegaX0) / (a * a) + OmegaX0);

	return hubble_a * qromb(growth_int, 0, a);
}


void read_power_table(void)
{
	FILE *fd;
	char buf[500];
	double k, p, p1, p2;

	sprintf(buf, "inp.dat");

//	sprintf(buf, "spec.dat");
	
	printf("reading spectrum\n");

	if(!(fd = fopen(buf, "r")))
	{
		printf("can't read input spectrum in file '%s'\n", buf);
		exit(0);
		// FatalError(17);
	}

	NPowerTable = 0;
	do
	{
//		if(fscanf(fd, " %lg %lg ", &k, &p) == 2)
		if(fscanf(fd, " %lg %lg %lg %lg ", &k, &p, &p1, &p2) == 4)
			NPowerTable++;
		else
			break;
	}
	while(1);

	fclose(fd);

	//if(ThisTask == 0)
	{
		printf("found %d pairs of values in input spectrum table\n", NPowerTable);
		fflush(stdout);
	}


	PowerTable = malloc(NPowerTable * sizeof(struct pow_table));

	//  sprintf(buf, FileWithInputSpectrum);

	if(!(fd = fopen(buf, "r")))
	{
		printf("can't read input spectrum in file '%s'\n", buf);
		//    FatalError(18);
		exit(0);
	}

	NPowerTable = 0;
	do
	{
//		if(fscanf(fd, " %lg %lg ", &k, &p) == 2)
		if(fscanf(fd, " %lg %lg %lg %lg ", &k, &p, &p1, &p2) == 4)
		{
			PowerTable[NPowerTable].logk = log10(k);
			PowerTable[NPowerTable].logD = log10(p);
			NPowerTable++;
		}
		else
			break;
	}
	while(1);

	fclose(fd);

	//  qsort(PowerTable, NPowerTable, sizeof(struct pow_table), compare_logk);
}



double tk_eh(double k)          /* from Martin White */
{
  double q, theta, ommh2, a, s, gamma, L0, C0;
  double tmp;
  double omegam, ombh2, hubble;

  /* other input parameters */
  hubble = Hubble0;

  omegam = OmegaM0;

  ombh2 = 0.04 * Hubble0 * Hubble0;

  k *= (3.085678e24 / UnitLength_in_cm);        /* convert to h/Mpc */

  theta = 2.728 / 2.7;
  ommh2 = omegam * hubble * hubble;
  s = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * hubble;
  a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
    + 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
  gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * k * s)));
  gamma *= omegam * hubble;
  q = k * theta * theta / gamma;
  L0 = log(2. * exp(1.) + 1.8 * q);
  C0 = 14.2 + 731. / (1. + 62.5 * q);
  tmp = L0 / (L0 + C0 * q * q);
  return (tmp);
}


double PowerSpec_EH(double k)   /* Eisenstein & Hu */
{
  return Norm * k * pow(tk_eh(k), 2);
}





double PowerSpec_Tabulated(double k)
{
	double logk, logD, P, kold, u, dlogk, Delta2;
	int binlow, binhigh, binmid;

	kold = k;

	k *= (3.085678e24 / UnitLength_in_cm);     /* convert to h/Mpc */

	logk = log10(k);

	if(logk < PowerTable[0].logk || logk > PowerTable[NPowerTable - 1].logk)
		return 0;

	binlow = 0;
	binhigh = NPowerTable - 1;

	while(binhigh - binlow > 1)
	{
		binmid = (binhigh + binlow) / 2;
		if(logk < PowerTable[binmid].logk)
			binhigh = binmid;
		else
			binlow = binmid;
	}

	dlogk = PowerTable[binhigh].logk - PowerTable[binlow].logk;

	if(dlogk == 0){
		printf(" input spect dlog = 0\n");
		exit(0);
	}
	u = (logk - PowerTable[binlow].logk) / dlogk;

	logD = (1 - u) * PowerTable[binlow].logD + u * PowerTable[binhigh].logD;

	Delta2 = pow(10.0, logD);

	P = Norm * Delta2 / (4 * M_PI * kold * kold * kold);
//	P = Norm * Delta2 / (4 * M_PI * kold * kold);

	return P;
}

double PowerSpec(double k)
{
	double power, alpha, Tf;

	power =  PowerSpec_EH(k);

//	power = PowerSpec_Tabulated(k);
//	power *= pow(k, PrimordialIndex - 1.0);

	return power;
}

double r_tophat;

double sigma2_int(double k)
{
	double kr, kr3, kr2, w, x;

	kr = r_tophat * k;
	kr2 = kr * kr;
	kr3 = kr2 * kr;

	if(kr < 1e-8)
		return 0;

	w = 3 * (sin(kr) / kr3 - cos(kr) / kr2);
	x = 4 * M_PI * k * k * w * w * PowerSpec(k);

	return x;
}



double TopHatSigma2(double R)
{
	r_tophat = R;

	return qromb(sigma2_int, 0, 500.0 * 1 / R);   /* note: 500/R is here chosen as 
							 integration boundary (infinity) */
}



void initialize_powerspectrum(void)
{
	double res;

	int nproc, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank );


	double R8 = 8 * (3.085678e24 / UnitLength_in_cm);    /* 8 Mpc/h */

//	read_power_table();

	double Sigma8 = 0.809;
	Norm = 1.0;
	res = TopHatSigma2(R8);


	Norm = Sigma8 * Sigma8 / res;

//	printf("Normalization adjusted to  Sigma8=%g   (Normfac=%g)\n\n", Sigma8, Norm);

//	printf(" Omega = %lf %lf %lf\n", OmegaM0, OmegaX0, InitialTime);	
	Dplus = growth(1.0)/growth(InitialTime); // flat lcdm
#ifdef MYRAND
	seedi = 1237842 + rank;
#else
	srand((unsigned) rank);
#endif		
	//  Dplus = GrowthFactor(InitTime, 1.0);
}


void ps_ic_k_(FLOAT *kn, FLOAT *dre, FLOAT *dim) {
	double fac = pow(2*M_PI/BOXSIZE,1.5);
	double k2 = (*kn);
	double phase, k = sqrt(k2);
	double p_of_k = 2* PowerSpec(k);
	double amp = fac*sqrt(p_of_k)/Dplus ;
	double rn;
#ifdef MYRAND
	if ( 0 == direct_ic ) {

		rn = ran31(&seedi)/2147483647.0;
		while (rn == 0.0)
			rn = ran31(&seedi)/2147483647.0;
		amp *= sqrt(-log(rn)) ;

		rn = ran31(&seedi)/2147483647.0;
		while (rn == 0.0)
			rn = ran31(&seedi)/2147483647.0;

		phase = (FLOAT)( 2.0*M_PI*rn );

	}

	if ( 1 == direct_ic ) {

		rn = ran32(&seedi)/2147483647.0;
		while (rn == 0.0)
			rn = ran32(&seedi)/2147483647.0;
		amp *= sqrt(-log(rn)) ;

		rn = ran32(&seedi)/2147483647.0;
		while (rn == 0.0)
			rn = ran32(&seedi)/2147483647.0;

		phase = (FLOAT)( 2.0*M_PI*rn );


	}
	if ( 2 == direct_ic ) {

		rn = ran33(&seedi)/2147483647.0;
		while (rn == 0.0)
			rn = ran33(&seedi)/2147483647.0;
		amp *= sqrt(-log(rn)) ;

		rn = ran33(&seedi)/2147483647.0;
		while (rn == 0.0)
			rn = ran33(&seedi)/2147483647.0;

		phase = (FLOAT)( 2.0*M_PI*rn );


	}
#else
	rn = rand()/2147483647.0;
	while (rn == 0.0)
		rn = rand()/2147483647.0;
	amp *= sqrt(-log(rn)) ;

	rn = rand()/2147483647.0;
	while (rn == 0.0)
		rn = rand()/2147483647.0;
	phase = (FLOAT)( 2.0*M_PI*rn );
#endif
	 *dre = (FLOAT)amp*sin(phase);
	*dim = (FLOAT)amp*cos(phase);
}



void icgradient(int dir)
{
	double time_0, time_1, t0, t1;
	time_0 = dtime();
	int n, c ;
	int nproc, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank );

	if (0 == rank)
	printf(" ic gradient direction [%d]\n", dir);

	////// setup /////
	direct_ic = dir;
	initialize_powerspectrum();


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


	MPI_Barrier(MPI_COMM_WORLD);

	long meshsize = isize[0] * isize[1] * isize[2];
	if (0==rank) {
printf(" mesh size = %ld, %d %d %d\n", meshsize, isize[0], isize[1], isize[2]);
fflush(stdout);

	}
	double* mesh = (double*)pmalloc(sizeof(double) * meshsize, 50);

	for (n=0; n<meshsize; n++) 
		mesh[n] = 0.0;

	int i,j,k,idx;
	int ii, jj, kk;
	double norm = NSIDE/BOXSIZE;
	double delta = 1.0/norm;

	double wi, wj , wk, win, wjn, wkn;
/*
	c=0;	
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
//	double renormal = (NSIDE/BOXSIZE);
	double renormal = (NSIDE);
	renormal = renormal*renormal*renormal;
*/
	c=0;
	for (i=0; i<isize[0]; i++)
		for (j=0; j<isize[1]; j++)
			for (k=0; k<isize[2]; k++) {
				mesh[c] = 0;
			}

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
				//		printf(" j = %d %d - %d %d\n", j, rj, MSIZE0[0], MSIZE0[1]);
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
	MPI_Alltoall(nsendkey,1,MPI_INT,nrecvkey,1,MPI_INT,MPI_COMM_WORLD);

	nrecdisp[0] = 0;
	long tot_recvkey = nrecvkey[0];	
	for (n=1; n<nproc; n++) {
		nrecdisp[n] = nrecdisp[n-1] + nrecvkey[n-1];
		tot_recvkey += nrecvkey[n];
	}

	if (0==rank) {
printf(" tot recv = %ld\n", tot_recvkey);
fflush(stdout);

	}


	MKey *recvbuff = (MKey*)pmalloc(sizeof(MKey)*tot_recvkey, 56);

	MPI_Status status;
	//	MPI_Status *vstatus = (MPI_Status*)malloc(sizeof(MPI_Status)  * nproc);
	//	MPI_Request *vreq  = (MPI_Request*)malloc(sizeof(MPI_Request) * nproc);

	MPI_Status *vstatus = (MPI_Status*)pmalloc(sizeof(MPI_Status)  * nproc, 57);
	MPI_Request *vreq  = (MPI_Request*) pmalloc(sizeof(MPI_Request) * nproc, 58);

	time_0 = dtime();
	//#define MYALLTOALL
//#ifndef MYALLTOALL
//	MPI_Alltoallv(sendbuff, nsendkey, nsendisp, strMKey, recvbuff, nrecvkey, nrecdisp, strMKey, MPI_COMM_WORLD);
//#else

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
//#endif

/*
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
*/
	int nside[3] = {NSIDE, NSIDE, NSIDE};
	double param[3] = { splitRadius, BOXSIZE, dir+1};

	genicgrad2(data,nside,param);

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
//	FILE *fd;
//	char fname[80];


//#ifndef MYALLTOALL

//	MPI_Alltoallv(recvbuff, nrecvkey, nrecdisp, strMKey, sendbuff, nsendkey, nsendisp, strMKey, MPI_COMM_WORLD);
	//	MPI_Alltoallv(recvbuff, nrecvkey, nrecdisp, MPI_BYTE, sendbuff, nsendkey, nsendisp, MPI_BYTE, MPI_COMM_WORLD);
//#else

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

//#endif

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
		idx1 = ((i)*isize[1] + j)*isize[2] + k;
dp[0] = mesh[idx1];

		idx1 = ((ii)*isize[1] + j)*isize[2] + k;
dp[1] = mesh[idx1];
		idx1 = ((i)*isize[1] + jj)*isize[2] + k;
dp[2] = mesh[idx1];
		idx1 = ((ii)*isize[1] + jj)*isize[2] + k;
dp[3] = mesh[idx1];
		idx1 = ((i)*isize[1] + j)*isize[2] + kk;
dp[4] = mesh[idx1];
		idx1 = ((ii)*isize[1] + j)*isize[2] + kk;

dp[5] = mesh[idx1];
		idx1 = ((i)*isize[1] + jj)*isize[2] + kk;


dp[6] = mesh[idx1];
		idx1 = ((ii)*isize[1] + jj)*isize[2] + kk;


dp[7] = mesh[idx1];
		part[n].acc_pm[dir] = win*wjn*wkn*dp[0] 
			+ wi *wjn*wkn*dp[1]
			+ win*wj *wkn*dp[2]
			+ wi *wj *wkn*dp[3]
			+ win*wjn*wk *dp[4]
			+ wi *wjn*wk *dp[5]
			+ win*wj *wk *dp[6]
			+ wi *wj *wk *dp[7];

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


void ic_lcdm0(double box, int npartside) {
	int proc = this_domain;
	double xl[3], xr[3];
	xl[0] = xl[1] = xl[2] = 0;
	xr[0] = xr[1] = xr[2] = box;

	int left, right;
	int mid = PROC_SIZE;
	double frac;
	int i, d=0;

	proc = PROC_RANK;

	while (mid>1)
	{
		left = mid/2 + mid%2;
		right= mid - left;
		frac = (double)left/((double)mid);


		if (proc < left){        
			xr[d] = frac * (xr[d]-xl[d]) + xl[d];
			mid = left;
			//    printf(" [%d] left = %d \n", PROC_RANK, left);
		}
		else{
			xl[d] = frac * (xr[d]-xl[d]) + xl[d];
			proc-=left;
			mid = right;
			//    printf(" [%d] right= %d \n", PROC_RANK, right);
		}
		d = (d+1)%3;
	}


	long n;
	seed = 378412+PROC_RANK;


	int xistart[3];
	int xiend[3];
	double xi, dxi = BOXSIZE / npartside ;

	int npart = 0;

	for (d=0; d<3; d++) {
		i = 0;
		while (i<npartside) {
			xi = (0.5+i)*dxi;
			if (xi > xl[d])
				break;
			i++;
		}
		xistart[d] = i;
		while (i<npartside) {
			xi = (0.5+i)*dxi;
			if (xi > xr[d])
				break;
			i++;
		}
		xiend[d] = i-1;
	}


	npart = (xiend[0] - xistart[0] + 1) * (xiend[1]-xistart[1]+1) * (xiend[2]-xistart[2]+1);


	NPART = npart;

}



void ic_lcdm1(double box, int npartside) {
	int proc = this_domain;
	double xl[3], xr[3];
	xl[0] = xl[1] = xl[2] = 0;
	xr[0] = xr[1] = xr[2] = box;

	int left, right;
	int mid = PROC_SIZE;
	double frac;
	int i,j,k, d=0;

	proc = PROC_RANK;

	while (mid>1)
	{
		left = mid/2 + mid%2;
		right= mid - left;
		frac = (double)left/((double)mid);


		if (proc < left){        
			xr[d] = frac * (xr[d]-xl[d]) + xl[d];
			mid = left;
			//    printf(" [%d] left = %d \n", PROC_RANK, left);
		}
		else{
			xl[d] = frac * (xr[d]-xl[d]) + xl[d];
			proc-=left;
			mid = right;
			//    printf(" [%d] right= %d \n", PROC_RANK, right);
		}
		d = (d+1)%3;
	}


	long n;
	seed = 378412+PROC_RANK;


	int xistart[3];
	int xiend[3];
	double xi, dxi = BOXSIZE / npartside ;

	int npart = 0;

	for (d=0; d<3; d++) {
		i = 0;
		while (i<npartside) {
			xi = (0.5+i)*dxi;
			if (xi > xl[d])
				break;
			i++;
		}
		xistart[d] = i;
		while (i<npartside) {
			xi = (0.5+i)*dxi;
			if (xi > xr[d])
				break;
			i++;
		}
		xiend[d] = i;
	}

	npart = (xiend[0] - xistart[0] ) * (xiend[1]-xistart[1]) * (xiend[2]-xistart[2]);

	n = 0;
	for (i=xistart[0]; i<xiend[0]; i++)
		for (j=xistart[1]; j<xiend[1]; j++)
			for (k=xistart[2]; k<xiend[2]; k++) {

				part[n].pos[0] = (i+0.5)*dxi;
				part[n].pos[1] = (j+0.5)*dxi;
				part[n].pos[2] = (k+0.5)*dxi;

				part[n].vel[0] = 0.0;
				part[n].vel[1] = 0.0;
				part[n].vel[2] = 0.0;

				part[n].acc[0] = 0.0;
				part[n].acc[1] = 0.0;
				part[n].acc[2] = 0.0;


				part[n].acc_pm[0] = 0.0;
				part[n].acc_pm[1] = 0.0;
				part[n].acc_pm[2] = 0.0;
			
				n++;
			}

	double pmass=(OmegaM0*3*0.01)/(8.0*M_PI*GravConst)*pow(box/npartside,3.0);
	MASSPART = pmass;
//	for (n=0; n<NPART; n++) {
//		part[n].mass = pmass;
//	}

}


void ic_lcdm2(double box, int npart) {
	// printf(" zeldovich \n"); 
	int n, d;
	double ai = InitialTime; 
	double hubble_ai =  0.1*sqrt(OmegaM0/(ai*ai*ai)+OmegaX0);
	double F_Omega_ai = pow(OmegaM0/(OmegaM0+ai*ai*ai*OmegaX0), 0.6);
	double vel_prefac = ai*hubble_ai*F_Omega_ai/sqrt(ai);
	vel_prefac *= pow(ai,1.5);

	double dis[3];

	for (d=0; d<3; d++) {

		for (n=0; n<npart; n++) {
			part[n].acc_pm[d] = 0.0;
		}
		icgradient(d);
/*
		for (n=0; n<npart; n++) {

			dis[d] = part[n].acc_pm[d];

			part[n].pos[d] += dis[d];

			while(part[n].pos[d] >= BOXSIZE) {
				part[n].pos[d] -= BOXSIZE;
			}
			while(part[n].pos[d] < 0.0) {
				part[n].pos[d] += BOXSIZE;
			}

			part[n].vel[d] = dis[d]*vel_prefac;
		}
*/

	}


	for (d=0; d<3; d++) {
/*
		for (n=0; n<npart; n++) {
			part[n].acc_pm[d] = 0.0;
		}
		icgradient(d);
*/
		for (n=0; n<npart; n++) {

			dis[d] = part[n].acc_pm[d];

			part[n].pos[d] += dis[d];

			while(part[n].pos[d] >= BOXSIZE) {
				part[n].pos[d] -= BOXSIZE;
			}
			while(part[n].pos[d] < 0.0) {
				part[n].pos[d] += BOXSIZE;
			}

			part[n].vel[d] = dis[d]*vel_prefac;
		}


	}

}


