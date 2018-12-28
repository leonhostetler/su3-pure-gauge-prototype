/**************************************************************
* main.c -- This is the main simulation script
* 
* Author: Leon Hostetler, Jul 16, 2018
*
**************************************************************/
#include <stdio.h>
#include <stdlib.h>              /* Random */
#include <math.h>                /* sqrt */ /* Remember to add -lm to compilation */
#include "su3matrices.h"
#include "complex.h"

const double EPS = 0.25;         /* Small number used in building SU2 matrices */
const int DIM = 4;               /* The number of dimensions */
const int N[4] = {12,12,12,12};  /* The lattice size in each direction */
const int NEQUI = 10;            /* The number of equilibrium sweeps to make */
const int NMEAS = 10;            /* The number of configurations to save */
const int NDISC = 1;             /* The number of configurations to discard between saved ones */
const double BETA = 3.0;
int start_type = 1;              /* 0 for cold start, 1 for hot start */

/* Must manually update these preprocessor directives */
#define SU3ELEMENTS 1000
#define LATTICESITES (12*12*12*12)

/* Function Prototypes */
void rsu2_near1();
su3 rsu3_near1();
void build_SU3_set();
void init_conf();
void print_links();
void update_conf();
void update_link(int n, int mu);
double action_diff(int n, int mu, su3 cand_link);
int find_site(int n, int direction[4]);
void n_to_xyzt(int n, int xyzt[DIM]);
int xyzt_to_n(int xyzt[4]);
void print_sites();
void run_tests();
double plaq_avg();

/* Global Variables */
su3 SU3_set[2*SU3ELEMENTS];  /* The set (i.e. vector) of random SU(3) matrices */
su3 links[LATTICESITES][4];  /* The array of link variables U[n][mu] */
int acpt_count = 0;          /* Count the Metropolis acceptances */
int metr_count = 0;          /* Count the Metropolis decisions */
int file_flag = 0;           /* Print action differences to file after equilibration */

/* Main Program */
int main()
{
    /* Run a few simple validation tests */
    run_tests();

    printf("\nSimulation Information:");
    printf("\nLattice Size: %d x %d x %d x %d", N[0],N[1],N[2],N[3]);
    printf("\nBETA  = %lf", BETA);
    printf("\n\nStatistics Setup:");
    printf("\nEPS   = %lf", EPS);
    printf("\nNEQUI = %d", NEQUI);
    printf("\nNMEAS = %d", NMEAS);
    printf("\nNDISC = %d", NDISC);

    printf("\n\nInitialization:");
    fflush(stdout);
    
    /* Initialize */
    init_conf();

    /*print_links();*/

    /* Calculate plaquette average after initialization */
    double plaqs;
    plaqs = plaq_avg();
    if(start_type==1) {
        printf("\nHot start plaquette average: %lf", plaqs);
        fflush(stdout);
    } else {
        printf("\nCold start plaquette average: %lf", plaqs);
        fflush(stdout);
    }

    /* Equilibrate */
    build_SU3_set();
    printf("\n\nEquilibration:");
    fflush(stdout);

    int i;
    for(i=0; i<NEQUI; ++i) {
        update_conf();
        build_SU3_set();
    }
    printf("\nEquilibration complete");
    fflush(stdout);

    /*print_links();*/

    printf("\n\nData:");

    /* Reset accept/reject counters for after equilibration */
    acpt_count = 0;
    metr_count = 0;
    file_flag = 1;

    /* Loop over iterations */
    double plaquette_avg = 0.0;
    double plaquette;
    int j,k;
    for(j=0; j<NMEAS; ++j) {
        /*printf("\nStarting measurement loop j = %d", j);*/
        fflush(stdout);
        for(k=0; k<NDISC; ++k) {
            update_conf();
            build_SU3_set();
        }
        /* Measurements or save configurations */
        plaquette = plaq_avg();
        plaquette_avg += plaquette;
        printf("\nMeasurement %d, plaq_avg = %lf", j, plaquette);
    }
    plaquette_avg = plaquette_avg/((double) NMEAS);
    printf("\nAverage plaquette average over NMEAS measurements: %lf", plaquette_avg);
    printf("\nMetropolis decisions: %d", metr_count);
    printf("\nAcceptances: %d", acpt_count);
    printf("\nAcceptance rate: %lf\n", acpt_count/((double) metr_count));

    /*print_links();*/
    /*print_sites();*/
}


/***********************************************************
init_conf -- initialize the lattice

parameters
  none

returns
  void

NOTE: Hot/Cold start choice is set in global variable start_type.
  Cold start sets initial link variables to identity matrix. Hot
  start sets initial link variables to random SU(3) matrices near
  unity.
***********************************************************/
void init_conf()
{
  
    if(start_type==1){  /* Hot start */
        int random_index;
        int M = 0;
        int N = 2*SU3ELEMENTS-1;
        int k,m;

        /* Fill initial configuration with random SU(3) matrices near unity */
        for(k=0; k<LATTICESITES; ++k) {
            build_SU3_set();
            for(m=0; m<DIM; ++m) {
                random_index = M + rand() / (RAND_MAX/(N-M+1)+1);
                su3 X = SU3_set[random_index];
                links[k][m] = X;
            }
        }

    } else {  /* Cold start */

        /* Create the SU(3) identity matrix */
        su3 identity;
        int i,j;
        for(i=0; i<3; ++i) {
            for(j=0; j<3; ++j) {
                identity.a[i][j].real = 0.0;
                identity.a[i][j].imag = 0.0;
            }
        }
        identity.a[0][0].real = 1.0;
        identity.a[1][1].real = 1.0;
        identity.a[2][2].real = 1.0;

        int k,m;
        for(k=0; k<LATTICESITES; ++k) {
            for(m=0; m<DIM; ++m) {
                links[k][m] = identity;
            }
        }
    }

    return;
}

/***********************************************************
print_links -- print all of the links in the lattice. May be
	useful for debugging.

parameters
  none

returns
  void
***********************************************************/
void print_links()
{
    int k,m;
    for(k=0; k<LATTICESITES; ++k) {
        for(m=0; m<DIM; ++m) {
            printf("\nn = %d, mu = %d\n", k, m);
            su3_pri(links[k][m]);
        }
        printf("\n-----------------------\n");
    }
    return;
}


/***********************************************************
update_conf -- make an entire sweep of the lattice by
	updating every link

parameters
  none

returns
  void
***********************************************************/
void update_conf()
{
    int n, mu;
    for(n=0; n<LATTICESITES; ++n) {
        for(mu=0; mu<DIM; ++mu) {
            update_link(n, mu);
        }
    }
    return;
}


/***********************************************************
update_link -- update a single link

parameters
  n  -- the lattice site associated with the link
  mu -- the direction associated with the link

returns
  void

NOTE: Naive implementation of random select. Will be biased.
NOTE: Naive (and wrong) implementation of random double.
***********************************************************/
void update_link(int n, int mu)
{
    /*
      Propose a new link by multiplying it by a randomly
      selected matrix from the SU3_set
    */
    int random_index;
    int M = 0;
    int N = 2*SU3ELEMENTS-1;
    random_index = M + rand() / (RAND_MAX/(N-M+1)+1);

    su3 X = SU3_set[random_index];
    su3 new_link;   /* The new link candidate */
    new_link = su3_mul(X, links[n][mu]);

    /* Calculate the action difference if new_link is accepted */
    double actdif;
    actdif = action_diff(n, mu, new_link);

    /* Metropolis accept/reject */
    double r;     /* A random number in [0,1) */
    r = rand()/((double) RAND_MAX);
    metr_count += 1;
    if(r <= exp(-actdif)) {
        links[n][mu] = new_link;
        acpt_count += 1;
    }

    return;
}


/***********************************************************
action_diff -- calculate the action difference if a single
	link is substituted

parameters
  n  -- the lattice site associated with the link
  mu -- the direction associated with the link
  cand_link -- the candidate link to replace the current link

returns
  the action difference

NOTE: Only works for DIM = 4
***********************************************************/
double action_diff(int n, int mu, su3 cand_link)
{
    /* Initialize A with zeros */
    su3 A;
    int i,j;
    for(i=0; i<3; ++i) {
        for(j=0; j<3; ++j) {
            A.a[i][j].real = 0.0;
            A.a[i][j].imag = 0.0;
        }
    }

    /* Initialize temporary SU3 matrices */
    su3 U1, U2, U3, U4, U5, U6;

    /* This loop is the implementation of the sum Eq. (4.20)
       given in Gattringer/Lang */
    int nu;
    for(nu=0; nu<DIM; ++nu) {
        if(nu == mu)
            continue;  /* Skip the case nu = mu */

        /*printf("\naction_diff loop n, mu, nu = %d, %d, %d", n, mu, nu);*/

        /* Initialize dim-dimensional "direction of movement" vectors */
        int d1[4] = {};
        int d2[4] = {};
        int d3[4] = {};
        int d4[4] = {};

        /* There are four distinct movements involved */
        d1[mu] = 1;   /* Move +1 in the direction mu */
        d2[nu] = 1;   /* Move +1 in the direction nu */
        d3[mu] = 1;   /* Move +1 in the direction mu, then -1 in direction nu */
        d3[nu] = -1;
        d4[nu] = -1;  /* Move -1 in the direction nu */
        
        /*printf(" d1 = %d %d %d %d", d1[0],d1[1],d1[2],d1[3]);
        printf(" d2 = %d %d %d %d", d2[0],d2[1],d2[2],d2[3]);
        printf(" d3 = %d %d %d %d", d3[0],d3[1],d3[2],d3[3]);
        printf(" d4 = %d %d %d %d", d4[0],d4[1],d4[2],d4[3]);*/

        /* Determine the links involved in this pair of staples */
        U1 = links[find_site(n,d1)][nu];
        U2 = links[find_site(n,d2)][mu];
        U3 = links[n][nu];
        U4 = links[find_site(n,d3)][nu];
        U5 = links[find_site(n,d4)][mu];
        U6 = links[find_site(n,d4)][nu];

        /* Add this pair of staples to A */
        su3 temp1, temp2;
        temp1 = su3_mul(su3_mul(U1, su3_hco(U2)), su3_hco(U3));
        temp2 = su3_mul(su3_mul(su3_hco(U4), su3_hco(U5)), U6);
        A = su3_add(su3_add(A, temp1), temp2);

    }

    /* Calculate the change in action */
    double deltaS;
    deltaS = -(BETA/3.0)*su3_rtr(su3_mul(su3_sub(cand_link, links[n][mu]),A));

    return(deltaS);
}


/***********************************************************
plaq_avg -- calculates the average plaquette value for the
	current lattice configuration

parameters
  none

returns
  the average value of the plaquettes over one lattice configuration

NOTE: Only works for DIM = 4
***********************************************************/
double plaq_avg()
{
    /* Initialize plaquette sum with zeros */
    su3 plaq;
    int i,j;
    for(i=0; i<3; ++i) {
        for(j=0; j<3; ++j) {
            plaq.a[i][j].real = 0.0;
            plaq.a[i][j].imag = 0.0;
        }
    }

    /* Temporary SU3 matrices */
    su3 U1, U2, U3, U4;

    /* Loop over plaquettes */
    int n, mu, nu;
    double sum = 0.0;
    for(n=0; n<LATTICESITES; ++n) {
        for(nu=0; nu<DIM; ++nu) {
            for(mu=0; mu<nu; ++mu) {

                /* Initialize dim-dimensional "direction of movement" vectors */
                int d1[4] = {};
                int d2[4] = {};

                /*printf("\naction_diff loop n, nu, mu = %d, %d, %d", n, nu, mu);*/

                /* There are two distinct movements involved */
                d1[mu] = 1;   /* Move +1 in the direction mu */
                d2[nu] = 1;   /* Move +1 in the direction nu */

                /*printf(" d1 = %d %d %d %d", d1[0],d1[1],d1[2],d1[3]);
                printf(" d2 = %d %d %d %d", d2[0],d2[1],d2[2],d2[3]);*/

                /* The links involved in this plaquette */
                U1 = links[n][mu];
                U2 = links[find_site(n,d1)][nu];
                U3 = links[find_site(n,d2)][mu];
                U4 = links[n][nu];

                /* Multiply the four link variables to get the plaquette */
                plaq = su3_mul(su3_mul(su3_mul(U1,U2),su3_hco(U3)),su3_hco(U4));

                /* Sum the real trace of all plaquettes */
                sum += su3_rtr(plaq); 
            }
        }
    }
    return (sum/((double) 18*LATTICESITES));
}


/***********************************************************
find_site -- given the address n of a lattice site and
	a movement vector, this function returns the lattice
	site at the new location

parameters
  n  -- the lattice site
  direction -- the vector specifying a movement from site n
	to a new site

returns
  the n associated with the new lattice site

NOTE: This is the function that enforces the boundary conditions
	(periodic). A move to a site outside the lattice is
	mapped back to a site inside the lattice.
***********************************************************/
int find_site(int n, int direction[DIM])
{
    /* Input validation */
    if((n<0) || (n>LATTICESITES-1))
        printf("find_site: (WARNING) Input outside lattice!");

    /* A vector to hold the x,y,z,t coordinates of old site */
    /* Convert n --> x,y,z,t */
    int old_site[DIM];
    n_to_xyzt(n, old_site);

    /* Add the direction vector to the old site to get new site */
    int new_site[DIM];
    int i;
    for(i=0; i<DIM; ++i) {
        new_site[i] = old_site[i] + direction[i];
    }

    /* Now we enforce the boundary conditions */
    for(i=0; i<DIM; ++i) {
        if(new_site[i] < 0)
            new_site[i] = new_site[i] + N[i];
        if(new_site[i] > N[i]-1)
            new_site[i] = new_site[i] - N[i];
    }

    /* Convert back: x,y,z,t --> n */
    int new_n;
    new_n = xyzt_to_n(new_site);

    /* Output validation */
    if((new_n<0) || (new_n>LATTICESITES-1))
        printf("find_site: (WARNING) Output outside lattice!");
    
    return (new_n);
}


/***********************************************************
n_to_xyzt -- this function converts a lattice address of the
	form (n) to one of the form (x,y,z,t)

parameters
  n -- the lattice site
  xyzt -- a pointer to a 4D vector to store the result

returns
  void (function updates nonlocal vector)

NOTE: Only works for DIM=4
***********************************************************/
void n_to_xyzt(int n, int xyzt[4])
{
    /* Input Validation */
    if((n<0)||(n>LATTICESITES-1))
        printf("\nn_to_xyzt: (WARNING) requested site outside lattice!"); 
    
    int temp = n;
    xyzt[0] = temp % N[0];
    temp = (temp-xyzt[0])/N[0];
    xyzt[1] = temp % N[1];
    temp = (temp-xyzt[1])/N[1];
    xyzt[2] = temp % N[2];
    xyzt[3] = (temp-xyzt[2])/N[2];

    return;
}


/***********************************************************
xyzt_to_n -- this function converts a lattice address of the
	form (x,y,z,t) to one of the form (n)

parameters
  xyzt -- a pointer to the old form in x,y,z,t form

returns
  the new lattice site in form (n)

NOTE: Only works for DIM=4
***********************************************************/
int xyzt_to_n(int xyzt[4])
{
    /* Input validation */
    int k;
    for(k=0; k<DIM; ++k) {
        if((xyzt[k]<0) || (xyzt[k]>N[k]-1))
            printf("\nxyzt_to_n: (WARNING) Site outside lattice!");
    }

    int x,y,z,t,n;

    x = xyzt[0];
    y = xyzt[1];
    z = xyzt[2];
    t = xyzt[3];

    n = t*N[0]*N[1]*N[2] + z*N[0]*N[1] + y*N[0] + x;

    return (n);
}


/***********************************************************
rsu2_near1 -- given the pointer to a complex 2x2 matrix, this
	function fills that matrix as a random SU(2) matrix
	near unity

parameters
  su2 -- the pointer to the matrix 

returns
  void

NOTE: The random number generator used here may not be appropriate
***********************************************************/
void rsu2_near1(complex su2[2][2])
{
    /* Select four random numbers */
    double r[4];
    int i = 0;   
    for(i=0; i<4; ++i)
        r[i] = (((double)rand()+1)/((double)RAND_MAX+2))-0.5;

    double magr;
    magr = sqrt(r[1]*r[1] + r[2]*r[2] + r[3]*r[3]); 
    
    int sign;
    if(r[0] < 0) {
        sign = -1;
    } else if (r[0] > 0) {
        sign = 1;
    } else {
        sign = 0;
    }

    double x[4];
    x[0] = sign*sqrt(1-EPS*EPS);
    
    for(i=1; i<4; ++i)
        x[i] = EPS*r[i]/magr;

    su2[0][0].real = x[0];
    su2[0][0].imag = x[3];
    su2[0][1].real = x[2];
    su2[0][1].imag = x[1];
    su2[1][0].real = -x[2];
    su2[1][0].imag = x[1];
    su2[1][1].real = x[0];
    su2[1][1].imag = -x[3];

    return;
}


/***********************************************************
rsu3_near1 -- this function returns a random SU(3) matrix
	near unity

parameters
  none

returns
  an SU(3) matrix structure
***********************************************************/
su3 rsu3_near1()
{
    /* Declare 3 SU(3) Matrices and initialize to zeros */
    su3 R, S, T;
    int i = 0;
    int j = 0;
    complex zero = {0,0};
    for(i=0; i<3; ++i) {
        for(j=0; j<3; ++j) {
            R.a[i][j] = zero;
            S.a[i][j] = zero;
            T.a[i][j] = zero;
        }
    }

    /* Get three random SU(2) matrices near unity */
    complex Q[2][2], V[2][2], W[2][2];
    rsu2_near1(Q);
    rsu2_near1(V);
    rsu2_near1(W);

    /* Fill the SU(3) matrices according to Gattringer/Lang Eq. (4.31) */
    R.a[0][0] = Q[0][0];
    R.a[0][1] = Q[0][1];
    R.a[1][0] = Q[1][0];
    R.a[1][1] = Q[1][1];
    R.a[2][2].real = 1.0;

    S.a[0][0] = V[0][0];
    S.a[0][2] = V[0][1];
    S.a[2][0] = V[1][0];
    S.a[2][2] = V[1][1];
    S.a[1][1].real = 1.0;

    T.a[1][1] = W[0][0];
    T.a[1][2] = W[0][1];
    T.a[2][1] = W[1][0];
    T.a[2][2] = W[1][1];
    T.a[0][0].real = 1.0;

    /* Create the SU(3) matrix by multiplying R*S*T */
    su3 rsu3matrix;
    rsu3matrix = su3_mul(su3_mul(R, S), T);

    return (rsu3matrix);

}


/***********************************************************
build_su3_set -- this function builds a set of random SU(3)
	matrices. The Hermitian conjugate of each element
	is also in the set.

parameters
  none

returns
  void
***********************************************************/
void build_SU3_set()
{
    int i=0;
    su3 random_su3;
    for(i=0; i<SU3ELEMENTS; ++i) {
        random_su3 = rsu3_near1();
        SU3_set[i] = random_su3;
        SU3_set[SU3ELEMENTS+i] = su3_hco(random_su3);
    }
}


/***********************************************************
print_sites -- this function prints all the lattice sites in
	both the xyzt and n representation. May be useful
	for debugging

parameters
  none

returns
  void

NOTE: Only works for DIM=4
***********************************************************/
void print_sites()
{
    printf("\nList all lattice sites");
    printf("\n n    x,y,z,t");
    int n;
    int xyzt[4];
    for(n=0; n<LATTICESITES; ++n) {
        n_to_xyzt(n, xyzt);
        printf("\n %d    %d %d %d %d", n, xyzt[0], xyzt[1], xyzt[2], xyzt[3]);
    }

    return;
}


/***********************************************************
run_tests -- this function runs various tests which may be 
	useful when debugging

parameters
  none

returns
  voids

Tests the following functions:
  -- xyzt_to_n
  -- n_to_xyzt
  -- find_site
***********************************************************/
void run_tests()
{

    /* Test xyzt_to_n and n_xyzt */
    int myxyzt[4];
    int t,q;
    for(t=0; t<LATTICESITES; ++t) {
        n_to_xyzt(t, myxyzt);
        q = xyzt_to_n(myxyzt);
        if(t != q)
            printf("WARNING: Incompatibility between xyzt_to_n and n_to_xyzt");
        myxyzt[0] = 0;
        myxyzt[1] = 0;
        myxyzt[2] = 0;
        myxyzt[3] = 0;
    }


    /* Validate find_site */
    int move_vector[4] = {0};
    int r, newsite;
    /* Not moving should return the same site */
    for(r=0; r<LATTICESITES; ++r) {
        newsite = find_site(r, move_vector);
        if(r != newsite)
            printf("\nfind_site: Error1");        
    }

    /* Moving by [Nx,Ny,Nz,Nt] should return to same site */
    for(r=0; r<DIM; ++r)
        move_vector[r] = N[r];
    for(r=0; r<LATTICESITES; ++r) {
        newsite = find_site(r, move_vector);
        if(r != newsite)
            printf("\nfind_site: Error2");        
    }

    /* Moving forward then backward in any direction and from
     any site should return one to the same site */
    int s, forw, backw;
    s = 0;
    for(s=0; s<DIM; ++s) {
        move_vector[0] = 0;
        move_vector[1] = 0;
        move_vector[2] = 0;
        move_vector[3] = 0;
        for(r=0; r<LATTICESITES; ++r) {
            move_vector[s] = 1; /* Move forward one unit in s direction */
            forw = find_site(r, move_vector);
           move_vector[s] = -1; /* Move back one unit in s direction */
            backw = find_site(forw, move_vector);
            /*printf("\nr, forw, backw, s: %d %d %d %d", r, forw, backw, s);*/
            if(backw != r)
                printf("\nfind_site: Error3");        
        }
    }
    return;
}

