/*************************************************************
 * Copyright (c) 2009 Steven McKay Curtis and Garritt L. Page
 *
 * I give permission for anyone to use and/or modify these
 * programs provided the following conditions are met:
 *
 * 1) This copyright message is retained.
 * 2) Any changes to the code are noted.
 *
 * These programs are provided WITHOUT ANY WARRNTY.
 *
 *************************************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>

#include <R_ext/Lapack.h>


//============================================================
// Allocating and Freeing Vectors and Matrices
//
// Note that the Matrix allocation functions allocate the
// memory for the entire matrix in a contiguous block.

//==================================================
// Vectors
double* R_Vector(int n){
    double* v;
    v = (double *) R_alloc(n,sizeof(double));
    return(v);
}

double* R_VectorInit(int n,double init){
    double* v;
    v = (double *) R_alloc(n,sizeof(double));
    for(int i=0; i<n; i++)
	v[i]=init;
    return(v);
}

/* For possible Calloc version
void R_FreeVector(double *v, int n){
    Free(v);
}
*/

//==================================================
// Matrices
double** R_Matrix(int nr, int nc){
    double** m;
    m = (double **) R_alloc(nr,sizeof(double*));
    *m = (double *) R_alloc(nr*nc,sizeof(double));
    for (int i=1; i<nr; i++)
	*(m+i) = *m + i*nc;
    return(m);
}

double** R_MatrixInit(int nr, int nc, double init){
    double** m;
    int i, j;
    m = (double **) R_alloc(nr,sizeof(double*));
    *m = (double *) R_alloc(nr*nc,sizeof(double));
    for (i=1; i<nr; i++)
	*(m+i) = *m + i*nc;
    for (i=0; i<nr; i++)
	for (j=0; j<nc; j++)
	    m[i][j] = init;
    return(m);
}

/* for possible Calloc versions
void R_FreeMatrix(double** m, int nr, int nc){
   Free(*m); // because it was allocated in a contiguous block
   Free(m);
//for (int i=nr-1; i>=nr; i--)
//Free(m[i]);
//Free(m);
}
*/

double** R_Data2Matrix(double* d, int nr, int nc){
    double** m;
    m = (double **) R_alloc(nr,sizeof(double*));
    for (int i=0; i<nr; i++)
	*(m+i) = d + i*nc;
    return(m);
}

int** R_Data2iMatrix(int* d, int nr, int nc){
    int** m;
    m = (int **) R_alloc(nr,sizeof(int*));
    for (int i=0; i<nr; i++)
	*(m+i) = d + i*nc;
    return(m);
}

//============================================================
// Matrix Operations
//


double quform(double *x, double* A, int dim){
    int i,j;
    double sm=0.0;
    for (i=1; i<dim; i++)
	for (j=0; j<i; j++)
	    sm += x[i]*x[j]*A[i*dim+j];
    sm*=2;
    for(i=0; i<dim; i++)
	sm += x[i]*x[i]*A[i*dim+i];
    return(sm);
}


double biform(double* x, double* A, double* y, int l){
    int i,j;
    double sm=0.0;
    for (i=0; i<l; i++)
	for (j=0; j<l; j++)
	    sm += x[i]*A[i*l+j]*y[j];
    return(sm);
}



//************************************************************
// Printing Matrices and Vectors
//
// The following functions
//
// Rprintvec, Rprintmat, RprintIvec, RprintImat,
//
// are modified versions of functions
// provided by Howard Seltman at the following web page:
// http://www.stat.cmu.edu/~hseltman/files/vmr.c
//
// I have modified the functions to work with R and to
// provide slightly modified output to suit my tastes.
//

// for doubles
void Rprintvec(char* title, double *v, int l){
    if (title!=NULL)
	Rprintf("%s\n",title);
    for (int i=0; i<l; i++)
	Rprintf("%f\n", v[i]);
    Rprintf("\n");
    return;
}

void Rprintmat(char* title, double **m, int nr, int nc){
    if (title!=NULL)
	Rprintf("%s\n", title);
    for (int i=0; i<nr; i++) {
	for (int j=0; j<nc; j++)
	    Rprintf("%f ", m[i][j]);
	Rprintf("\n");
    }
    Rprintf("\n");
    return;
}



// for integers
void RprintIvec(char* title, int* v, int n){
    if (title!=NULL)
	Rprintf("%s\n", title);
    for (int i=0; i<n; i++)
	Rprintf("%i\n", v[i]);
    Rprintf("\n");
    return;
}

void RprintImat(char* title, int** m, int nr, int nc){
    if (title!=NULL)
	Rprintf("%s\n", title);
    for (int i=0; i<nr; i++) {
	for (int j=0; j<nc; j++)
	    Rprintf("%i ", m[i][j]);
	Rprintf("\n");
    }
    Rprintf("\n");
    return;
}


// Print a vector as a matrix
void RprintVecAsMat(char* title, double *v, int nr, int nc){
    if (title!=NULL)
	Rprintf("%s\n", title);
    for (int i=0; i<nr; i++) {
	for (int j=0; j<nc; j++)
	    Rprintf("%f ", v[i*nc + j]);
	Rprintf("\n");
    }
    Rprintf("\n");
    return;
}


// Print a vector as a matrix
void RprintIVecAsMat(char* title, int *v, int nr, int nc){
    if (title!=NULL)
	Rprintf("%s\n", title);
    for (int i=0; i<nr; i++) {
	for (int j=0; j<nc; j++)
	    Rprintf("%d ", v[i*nc + j]);
	Rprintf("\n");
    }
    Rprintf("\n");
    return;
}


// Matrix tranpose
void mat_transpose(double *mat, double *tmat, int nr, int nc)
	{

//		RprintVecAsMat("Uk", mat, 1, (nr)*(nc));
		int i, j;
		for(i = 0; i < nr; i++)
			{
				for(j = 0; j < nc; j++)
					{

						tmat[j*nr + i] = mat[i*nc + j];


					}

			}
//		RprintVecAsMat("Uk", mat, 1, (nr)*(nc));

	}
// Matrix tranpose that destroys original matrix.
void mat_transpose2(double *m, int w, int h)
{
	int start, next, i;
	double tmp;

	for (start = 0; start <= w * h - 1; start++) {
		next = start;
		i = 0;
		do {	i++;
			next = (next % h) * w + next / h;
		} while (next > start);
		if (next < start || i == 1) continue;

		tmp = m[next = start];
		do {
			i = (next % h) * w + next / h;
			m[next] = (i == start) ? tmp : m[i];
			next = i;
		} while (next > start);
	}
}

// trace of a square matrix.
double mat_trace(double *matrix, int w)
{
	int i, ii;
	double sum;

 	sum = 0.0;
	for (i = 0; i < w; i++){
		for(ii = 0; ii < w; ii++){
			if(i == ii) sum = sum + matrix[i*w+ii];
		}
	}

	return(sum);
}


// factorial function using recursive arguments
int factorial(int n){
//	Rprintf("n = %d\n", n);
	int fact = 1;
	if(n < 0){ Rprintf("Cannot compute factorial of negative number "); return(0);}
	if(n <= 1 ){
		return(1);
	}else{
		fact = n*factorial(n-1);
	}
	return(fact);
}


// function that computes the mean of an contiguous array of memory
double mean(double* x, int n){
	int i;
	double mn;
	mn = 0.0;
	for(i = 0; i < n; i++){
		mn = mn + x[i]/((double) n);
	}

	return(mn);

}

//  Function that computes the sample variance (using n as the denominator)
double var(double* x, int n){

	int i;
	double mn, mn2, vr;
	mn = 0.0, mn2=0.0;
	for(i = 0; i < n; i++){
		mn = mn + x[i]/((double) n);
		mn2 = mn2 + (x[i]*x[i])/((double) n);
	}

	vr = mn2 - mn*mn;

	return(vr);

}






//============================================================
// Random Number Generation
//

void ran_mvnorm(double* m, double* cholV, int dim, double* z, double* out){
/*************************************************************
 * PURPOSE:
 * Generate a draw from a multivariate normal distribution.
 *
 * INPUTS:
 * m     = an array for the mean
 * cholV = an array for the cholesky decomposition
 *         of the variance (note this is a 1D array that
 *         that will be treated as 2D).
 * dim   = dimension of the multivariate normal distribution
 * z     = a scratch array of length dim
 *         to be filled with N(0,1) draws
 *
 * OUTPUTS:
 * out   = final output array to be filled with a random
 *         multivariate normal draw
 *
 *************************************************************/
    int i,j;
    for (i=0; i<dim; i++){
	z[i] = rnorm(0,1);
	out[i] = m[i];
	for (j=0; j<=i; j++)
	    out[i] += cholV[i*dim + j]*z[j];
    }
}


void ran_wish(int nu, double* cholS, int dim, double* z, double* x, double* zeros, double* out){
/*************************************************************
 * PURPOSE:
 * Generate a random draw from Wishart distribution with
 * degrees of freedom nu and parameter matrix S
 *
 * INPUTS:
 * nu    = degrees of freedom
 * cholS = cholesky decomposition of matrix parameter S
 *         This is a 1D array but is accessed
 *         as a 2D array as in cholS[i*dim + j]
 * dim   = the dimension of the Wishart distribution
 * z     = scratch vector of length dim to be passed
 *         to the function ran_mvnorm
 * x     = scratch vector of length dim to be passed
 *         to the function ran_mvnorm
 * zeros = vector of zeros of length dim to be passed
 *         to the function ran_mvnorm as the mean
 * out   = matrix to contain the random draw from the
 *         Wishart distribution
 *
 *************************************************************/

    int i, j, k;

    /* Zero the "out" matrix */
    for (j=0; j<dim; j++)
	for (k=0; k<dim; k++)
	    out[j*dim + k] = 0.0;

    for (i=0; i<nu; i++){
	ran_mvnorm(zeros,cholS,dim,z,x);
	for (j=0; j<dim; j++)
	    for (k=0; k<=j; k++)
		out[j*dim + k] += x[j]*x[k];
    }

    /* fill the upper triangular part with lower triangular part */
    for (j=0; j<dim; j++)
	for (k=0; k<j; k++)
	    out[k*dim + j] = out[j*dim + k];
}


void ran_dirich(double *alpha, int dim, double* scratch, double* out){
/*************************************************************
 * PURPOSE:
 * Generate a random draw from a Dirichlet distribution with
 * parameters alpha
 *
 * INPUTS:
 * alpha   = parameter vector associated with dirichlet
 * dim     = dimension of alpha vector
 * scratch = vector of length dim to be passed
 * out     = array that holds the random values
 *
 *************************************************************/

    int h;
	double sg = 0;
    /* Zero the "out" matrix */
//	RprintVecAsMat("alpha", alpha, 1, dim);
	for(h = 0; h < dim; h++)
		{

			scratch[h] = rgamma(alpha[h], 1);
//			Rprintf("rgamma = %f\n", scratch[h]);
			sg = sg + scratch[h];
//			Rprintf("sg = %f\n", sg);
		}

//	Rprintf("sg = %f\n", sg);
	for(h = 0; h < dim; h++) out[h] = scratch[h]/sg;

//	RprintVecAsMat("out", out, 1, dim);
}





/* The following provides a function that allows me to sample from a
   truncated normal distribution with the following arguments.  This
   function relies heavly on the Rmath library.  As an example here
   are the arguments for the pnorm function
   double pnorm(double x, double mu, double sigma, int lower_tail,
				int give_log );


  	m - mean of the normal prior to truncation;
  	s - standard deviation of normal distribution prior to truncation;
  	a - Lower bound of the truncation
  	b - Upper bound of the truncation */

double rtnorm(double m, double s, double a, double b){

  int jj;
  double a_term, b_term, Un, p, rtdraw, tmp, z;

//	Rprintf("a = %f\n", a);
//	Rprintf("b = %f\n", b);

//	Rprintf("(a-m)/s = %f\n", (a-m)/s);
	a_term = pnorm((a-m)/s, 0.0, 1.0, 1, 0);

//	Rprintf("(b-m)/s = %f\n", (b-m)/s);
	b_term = pnorm((b-m)/s, 0.0, 1.0, 1, 0);

//	Rprintf("a_term = %f\n", a_term);
//	Rprintf("b_term = %f\n", b_term);

	Un = runif(0,1);

	p = a_term + Un*(b_term - a_term);

	rtdraw = m + s*qnorm(p, 0.0, 1.0, 1, 0);

//	Rprintf("rtdraw = %f\n", rtdraw);
//	Rprintf("p = %f\n", p);

	if(p == 1.0){
		jj = 1;

		while(jj != 0){

			a = (a - m)/s;
    		b = (b - m)/s;

//			Rprintf("a = %f\n", a);
			tmp = (a + sqrt(a*a + 4))/2;
//			Rprintf("tmp = %f\n", tmp);
        	z = rexp(1/tmp) + a;
//			Rprintf("z = %f\n", z);
        	Un = runif(0,1);
//			Rprintf("Un = %f\n", Un);
//			Rprintf("exp(-(z - tmp)*(z - tmp)/2) = %f\n", exp(-(z - tmp)*(z - tmp)/2));

        	if((Un <= exp(-(z - tmp)*(z - tmp)/2)) & (z <= b)) jj = 0;
//			Rprintf("jj = %d\n", jj);

		}

//		Rprintf("z = %f\n", z);

		rtdraw = z*s + m;


	}


//	Rprintf("rtdraw = %f\n", rtdraw);

	return(rtdraw);

}


//============================================================
// Multivariate density functions
//
/* The following provides a density function for the multivariate
   normal distribution.  This function relies heavily on matrix.c
   To use it I must create an matrix.o(object file)

	*y -  is the observation vector for which the density will be computed
  	*mu - mean vector
  	*iSig - the inverse covariance matrix as a contiguos vector in row-major form
  	dim - is the dimension of the multivariate distribution
	ld - the log of the determinant of Sigma
	scr - scratch vector to hold
	logout - logical determines if log density is returned
*/

//Multivariate normal
double dmvnorm(double* y, double* mu, double* iSig, int dim, double ld, double* scr, int logout){
    int i;
    double qf, out;
    for(i=0; i<dim; i++) scr[i] = y[i] - mu[i];
    qf = quform(scr,iSig,dim);
    out = -(double) dim*M_LN_SQRT_2PI - 0.5*(ld + qf);
    if (logout)	return out;
    return exp(out);
}



//============================================================
// univariate scaled and shifted t-distribution
//

//
double dsst(double y, double mu, double s, double nu, int logout){

    double qf, out, lc;

	lc = lgamma(0.5*(nu + 1)) - (0.5*log(nu*M_PI) + lgamma(0.5*nu));
	qf = ((y-mu)/s)*((y-mu)/s);

    out = lc - 0.5*(nu+1)*log(1 + (1/nu)*qf);

    if (logout)	return out;
    return exp(out);
}



/* The following provides a density function for the Inverse Wishart
   distribution.  This function relies heavily on matrix.c
   To use it I must create an matrix.o(object file)
	Sig is argument and A is parameter.
	*ASigInv - is ASigma^{-1} that is found in the trace of the density.  It is a
	           contiguos vector of memory that is in row-major form
	detSig - is the determinant of the argument matrix
	detA - is the determinant of A
  	df - degrees of freedom of the inverse-wishart function
  	dim - is the dimension of the scale matrix distribution */

double dinvwish(double *ASigInv, double detSig, double detA, int df, int dim, int logout){

	int i;
	double out;
	double lgammasum;
	double trace;
	double p1, p2, p3, p4;

//	Rprintf("dim = %d\n", dim);

	lgammasum = 0.0;
	for(i = 0; i < dim; i++){
		lgammasum = lgammasum + lgammafn(0.5*(df + 1 - (i+1))); //I have to add 1 since C starts at zero
	}
//	Rprintf("lgammasum = %f\n", lgammasum);
	trace = 0.0;
	for(i = 0; i < dim*dim; i++){

		if(i % (dim+1) == 0) trace = trace + ASigInv[i];

	}
//	Rprintf("trace = %f\n", trace);
	p1 = 0.5*df*dim;
	p2 = 0.25*dim*(dim - 1);
	p3 = 0.5*df;
	p4 = 0.5*(df + dim + 1);
//	Rprintf("p1 = %f\n", p1);Rprintf("p2 = %f\n", p2);Rprintf("p3 = %f\n", p3);Rprintf("p4 = %f\n", p4);

	out = -p1*log(2) - (p2*log(M_PI) + lgammasum) + p3*log(detA) - p4*log(detSig) - 0.5*trace;
//	Rprintf("detA = %f\n", detA);	Rprintf("detSig = %f\n", detSig);
	if(logout) return out;
//	Rprintf("out = %f\n", out);

	return exp(out);

}
/*
   k <- nrow(S)
    gammapart <- 1
    for (i in 1:k) {
        gammapart <- gammapart * gamma((v + 1 - i)/2)
    }
    denom <- gammapart * 2^(v * k/2) * pi^(k * (k - 1)/4)
    detS <- det(S)
    detW <- det(W)
    hold <- S %*% solve(W)
    tracehold <- sum(hold[row(hold) == col(hold)])
    num <- detS^(v/2) * detW^(-(v + k + 1)/2) * exp(-1/2 * tracehold)
*/

// Inverse Gamma density
// the parameterization here provides a mean of beta/(alpha - 1)
double dinvgamma(double y, double alpha, double beta, int logout){

//	Rprintf("alpha = %f\n", alpha);
//	Rprintf("beta = %f\n", beta);

	double ldens;

	ldens = alpha*log(beta) - lgamma(alpha) - (alpha + 1)*log(y) - (beta/y);

	if(logout) return ldens;
	return exp(ldens);

}


// Normal-Inverse Gamma density
// the parameterization here is such that mu0 is prior mean value with k0 prior "observations"
// and s20 is the prior guess for sig2 with nu0 prior "observations"
double dN_IG(double mu, double sig2, double mu0, double k0, double a0, double b0, int logout){

//	Rprintf("alpha = %f\n", alpha);
//	Rprintf("beta = %f\n", beta);

	double ldens;

//	ldens =  0.5*(log(k0) - log(2*M_PI*sig2)) + a0*log(b0) -
//	         lgammafn(a0) - (a0 + 1)*log(sig2) -
//	         0.5*(1/sig2)*(k0*(mu-mu0)*(mu-mu0) + 2*b0);

	ldens = dnorm(mu, mu0, sqrt(sig2/k0),logout) +
	        dinvgamma(sig2, a0, b0, logout);
//	Rprintf("ldens = %f\n", ldens);
	if(logout){ return ldens;
	}else{return exp(ldens);}

}

// density for multinomial distribution
double dmultinom(int *y, double *prob, int n, int dim, int logout){

	int ii;
	double ldensity =0.0;

	ldensity = log(factorial(n));

	for(ii = 0; ii < dim; ii++){
		ldensity = ldensity + y[ii]*log(prob[ii]) - log(factorial(y[ii]));
	}

	if(logout) return ldensity;
	return exp(ldensity);
}



double ddirich(double *pi, double *alpha, int C, int logout){
/*************************************************************
 * PURPOSE:
 * Evaluate Dirichlet density with parameters alpha
 *
 * INPUTS:
 * pi   = C - dimensional simplex values
 * alpha     = C - dimensional parameter vector
 * C = dimension of probability vector
 * logout = logical indicating whether log density should be returned
 *
 *************************************************************/

    int ii;
	double suma = 0.0, sumla = 0.0, ldensity=0.0;
    /* Zero the "out" matrix */
//	RprintVecAsMat("alpha", alpha, 1, C);
//	RprintVecAsMat("pi", pi, 1, C);
	for(ii = 0; ii < C; ii++){

		suma = suma + alpha[ii];
		sumla = sumla + lgammafn(alpha[ii]);

	}

	for(ii = 0; ii < C; ii++){

		ldensity = ldensity + (alpha[ii]-1)*log(pi[ii]);

	}

	ldensity = ldensity + lgammafn(suma) - sumla;

	if(logout) return ldensity;
	return exp(ldensity);

}


// density for univariate truncated normal
double dtnorm(double x, double mu, double sigma, double l, double u, int logout){

	double den, num, ldensity;

	den = pnorm5(u,mu,sigma, 1, 0) - pnorm5(l,mu,sigma, 1, 0);
	num = dnorm(x, mu, sigma,1);

//	Rprintf("den = %f\n", den);

	ldensity = num - log(den);
	if(logout) return ldensity;
	return exp(ldensity);
}

// Determine how many unique elements are in an array
int unique_elements(int *arr, int len)
{
     if (len <= 0) return 0;
     int unique = 1;

     for (int outer = 1; outer < len; ++outer)
     {
        int is_unique = 1;
        for (int inner = 0; is_unique && inner < outer; ++inner)
        {
             if (arr[inner] == arr[outer]) is_unique = 0;
        }
        if (is_unique) ++unique;
     }
     return unique;
}


void cluster_indx(int *rho, int n, int nclus, int* clus_indx){
/*******************************************************************************
 * PURPOSE:
 * Take input vector of cluster labels and output indexes associated with
 * cluster membership.  For example input 1 2 1 2 3 provides output of
 * 1 3 2 4 5
 *
 * INPUTS:
 * nh - n dimensional integer array holding cluster labels
 * n - integer indicating the length of vector
 * nclus - integer indicating the number of clusters
 * clus_indx - n dimensional integer array holding cluster membership indexes
 *
 *******************************************************************************/
	int i, k, indx, min;

//	RprintIVecAsMat("rho", rho, 1, n);
	min=rho[0];
	for(i=1; i<n; i++){
		if(rho[i] < min) min = rho[i];
	}
//	Rprintf("min = %d\n", min);
//	Rprintf("nclus = %d\n", nclus);

	indx=0;
	for(k=0;k<nclus;k++){
		for(i=0;i<n;i++){
			if(rho[i] == k+min){
				clus_indx[indx] = i;
				indx = indx + 1;
			}
		}
	}
}


void remove_zero(int *rho, int n, int *rho_no_zero){
/*******************************************************************************
 * PURPOSE:
 * Take input vector and remove entries that are zero moving other entries down
 *
 * INPUTS:
 * nh - n dimensional integer array holding cluster sizes with possibly zero
 * n - integer indicating the length of vector
 *
 *******************************************************************************/
	int i, indx, nzero;

	nzero=0;
	for(i=0;i<n;i++){
		if(rho[i] != 0) ++nzero;
		rho_no_zero[i] = 0;
	}

	indx=0;
	for(i=0;i<n;i++){
		if(rho[i] != 0){
			rho_no_zero[indx] = rho[i];
			indx = indx + 1;
		}
		if(indx > nzero) rho_no_zero[i] = 0;
	}

}

void relabel(int *Si, int n, int *Sirelab, int *nhrelab, int *oldLab){
/*******************************************************************************
 * PURPOSE:
 * Relabel groups so that group with first unit is always 1, then group with
 * closest unit to the first that is not a part of the first group is 2, etc.
 *
 *
 * INPUTS:
 * Si - n dimensional integer array holding cluster labels, this will be d
 * n - integer indicating length of rho1 and rho2
 *
 * OUTPUTS:
 * Sirelab - n dimensional integer that will hold the relabeled group labels
 * nhrelab - n dimensional integer that will contain reorder cluster sizes
 *
 *******************************************************************************/
	int j;
	int shuffle = n, loc = 0, lab = 1;

	for(j = 0; j < n; j++){
		nhrelab[j] = 0, Sirelab[j] = 0;
	}

	while(shuffle > 0){
		for(j = 0; j < n; j++){
			if(Si[j] == Si[loc]){
				Sirelab[j] = lab;
				shuffle = shuffle - 1;
				nhrelab[lab-1] = nhrelab[lab-1] + 1;
			}
		}

		oldLab[lab-1] = Si[loc];

		lab = lab+1;
		for(j = 0; j < n; j++){
			if(Sirelab[j] == 0){
				loc = j;
				break;
			}
		}
	}
}

int compatibility(int *rho1, int *rho2, int n){
/*******************************************************************************
 * PURPOSE:
 * Compare rhot | gammat = 1 and rhot_1 | gammat = 1 to determine compatibility
 * I will relabel before comparison so that the first unit is always assigned to
 * cluster 1, the first unit not in cluster 1 will be assigned cluster 2, etc.
 *
 * INPUTS:
 * rho1 - n dimensional integer array holding cluster labels from reduced rho_t
 * rho2 - n dimensional integer array holding cluster labels from reduced rho_t-1
 * n - integer indicating length of rho1 and rho2
 *
 * OUTPUTS:
 * comp - binary indicating if two partitions are comatible
 * 		0 - not compatible
 * 		1 - compatible
 *******************************************************************************/
	int i, comp=1, scr1[n], scr2[n], scr3[n], scr4[n];
    //int nc1, nc2;
	//nc1 = unique_elements(rho1, n);
	//nc2 = unique_elements(rho2, n);

//	Rprintf("nc1 = %d\n", nc1);
//	Rprintf("nc2 = %d\n", nc2);


	relabel(rho1, n, scr1, scr3, scr4);
	relabel(rho2, n, scr2, scr3, scr4);

//	RprintIVecAsMat("rho1", rho1, 1, n);
//	RprintIVecAsMat("rho2", rho2, 1, n);


//	RprintIVecAsMat("relabelRho1", scr1, 1, n);
//	RprintIVecAsMat("relabelRho1", scr2, 1, n);


	for(i=0;i<n;i++){
		if((scr1[i] != scr2[i])){
			comp = 0;
			break;
		}
	}
	return(comp);
}



double partition_prob_crp(int *nh, int k, double M, int n, int logout){
/*******************************************************************************
 * PURPOSE:
 * Evaluate probability of partition based on CRP
 *
 * INPUTS:
 * nh - k dimensional integer array holding cluster sizes
 * k - integer indicating the number of clusters
 * M - precision parameter of DP
 * n - number of units that are being clustered
 * nh - scratch vector of length n which is used to determine cluster size
 * logout = logical indicating whether log density should be returned
 *
 *******************************************************************************/

    int i,j;
	double den = 0.0, num = 0.0, ldensity=0.0;

//	Rprintf("k = %d\n", k);
//	Rprintf("M = %f\n", M);
//	Rprintf("n = %d\n", n);
//	RprintIVecAsMat("nh = ", nh, 1, k);

	for(i = 0; i < n; i++){
		den = den + log(M + (i+1) - 1);
	}
	for(j = 0; j < k; j++){
		num = num + lgammafn(nh[j]);
	}


	ldensity = ((double) k)*log(M) + num - den;
	if(k == 0) ldensity = log(1); // This is an artifact of the MCMC algorithm

	if(logout) return ldensity;
	return exp(ldensity);

}





// Similarity function where the variance of x is computed only
double gsimconEV(double sumx, double sumx2, int n, double alpha, int logout){

//	Rprintf("sumx = %f\n", sumx);
//	Rprintf("sumx2 = %f\n", sumx2);
//	Rprintf("n = %d\n", n);
//	Rprintf("alpha = %f\n", alpha);

	double xbar, out, var;
	xbar = (1/(double) n)*sumx;
	var = ((1/(double) n)*sumx2 - xbar*xbar);
//	Rprintf("var = %f\n", var);
	out = -alpha*var;
	if(!logout) out = exp(out);

//	Rprintf("out = %f\n", out);


	return(out);

}




// normal-normal Similarity function with x following normal and m a normal as well (v is fixed).
// This is the result after integrating over mj.  One could also simply find the
// ratio between likelihood multiplied by prior divided by posterior

// The double dipper is included as an argument.

double gsimconNN(double m0, double v2, double s20, double sumx, double sumx2, double mle,
                 int n,  int DD, int cal, int logout){

//	Rprintf("m0 = %f\n", m0);
//	Rprintf("v2 = %f\n", v2);
//	Rprintf("s20 = %f\n", s20);
//	Rprintf("sumx = %f\n", sumx);
//	Rprintf("sumx2 = %f\n", sumx2);
//	Rprintf("n = %d\n", n);

	double mus, muss, s2s, s2ss;
	double ld1, ld2, ld3, ld4, ld5;
	double out;
	//double out1, out2, ld6;

	s2s = 1/((n/v2) + (1/s20));
	mus = s2s*((1/v2)*sumx + (1/s20)*m0);

//	Rprintf("mus = %f\n", mus);
//	Rprintf("s2s = %f\n", s2s);

	s2ss = 1/((n/v2) + (1/s2s));
	muss = s2ss*((1/v2)*sumx + (1/s2s)*mus);

//	Rprintf("muss = %f\n", muss);
//	Rprintf("s2ss = %f\n", s2ss);

	ld1 = -0.5*n*log(2*M_PI*v2) - 0.5*(1/v2)*sumx2;
	ld2 = dnorm(m0, 0, sqrt(s20),1);
	ld3 = dnorm(mus, 0, sqrt(s2s),1);
	ld4 = dnorm(muss, 0, sqrt(s2ss),1);

	ld5 = dnorm(mle, m0, sqrt(s20),1);
//	ld6 = dnorm(mle, mus, sqrt(s2s),1);
//	Rprintf("ld1 = %f\n", ld1);
//	Rprintf("ld2 = %f\n", ld2);
//	Rprintf("ld3 = %f\n", ld3);
//	Rprintf("ld4 = %f\n", ld4);


	out = ld1 + ld2 - ld3;
//	out1 = ld1 + ld3 - ld4;
//	out2 = ld5 - ld6;
//	Rprintf("out = %f\n", out);
//	Rprintf("out1 = %f\n", out1);
	if(DD==1) out = ld1 + ld3 - ld4;
	if(cal==1) out = ld5 - ld4;
	if(!logout) out = exp(out);
	return(out);

}






// normal-normal-IG Similarity function with x following normal and m,v a normal-IG.
// I didn't carry out integration explicitly over m and v. I simply used the fact that
// marginal likelihood (similarity) is equal to likelihood x prior / posterior.  This
// requires inserting a value for m and v which I use 0 and 1.

// The double dipper is included as an argument.


double gsimconNNIG(double m0, double k0, double nu0, double s20, double sumx, double sumx2,
				   double mnmle, double s2mle, int n, int DD, int cal, int logout){

//	Rprintf("sumx = %f\n", sumx);
//	Rprintf("sumx2 = %f\n", sumx2);
//	Rprintf("n = %d\n", n);

	double a0, b0, m0s, m0ss, k0s, k0ss, a0s, a0ss, b0s, b0ss;
	double ld1, ld2, ld3, ld4, ld5, ld6, out;
	double mu=10, v2=0.1;
	double xbar = sumx*(1/ (double) n);

	a0 = 0.5*nu0;
	b0 = 0.5*nu0*s20;


	m0s = (k0*m0 + n*xbar)/(k0 + n);
	k0s = k0 + (double) n;
	a0s = a0 + 0.5*n;
	b0s = b0 + 0.5*(sumx2 - n*xbar*xbar) + 0.5*n*k0*(xbar - m0)*(xbar - m0)/(k0+n);

	m0ss = (k0s*m0s + n*xbar)/(k0s + n);
	k0ss = k0s + (double) n;
	a0ss = a0s + 0.5*n;
	b0ss = b0s + 0.5*(sumx2 - n*xbar*xbar) + 0.5*n*k0s*(xbar - m0s)*(xbar - m0s)/(k0s+n);

	ld1 = -0.5*n*log(2*M_PI*v2) - 0.5*(1/v2)*(sumx2 - 2*sumx*mu + n*mu*mu);
	ld2 = dN_IG(mu, v2, m0, k0, a0, b0, 1);
	ld3 = dN_IG(mu, v2, m0s, k0s, a0s, b0s, 1);
	ld4 = dN_IG(mu, v2, m0ss, k0ss, a0ss, b0ss, 1);

	ld5 = dN_IG(mnmle, s2mle, m0, k0, a0, b0, 1);
	ld6 = dN_IG(mnmle, s2mle, m0s, k0s, a0s, b0s, 1);

	out = ld1 + ld2 - ld3;


	if(DD==1) out = ld1 + ld3 - ld4;
	if(cal==1) out = ld5 - ld6;
	if(!logout) out = exp(out);

//	Rprintf("out = %f\n", out);
	return(out);

}



// Similarity function with for a categorical x  dirichlet-multinomial with out
// where only on object is allocated (x_i is basically univariate that identifies which
// category ith individual has).  The integral in reality is a product of two ratios,
// but one of the ratios is constant in terms of $x$ and so disappears in the ratio
// of similarity functions when updating cluster labels and so is ignored in the
// function that follows.
double gsimcatDM(int* nobsj, double* dirweights, int C, int DD, int logout){

	int ii, sumc;
	double tmp1=0.0,tmp2=0.0,tmp3=0.0,tmp4=0.0,tmp5=0.0,tmp6=0.0;
	double out;
//	double ld1, ld2, ld3, ld4, out2, out3;
//	double *pi_vec = R_Vector(C);
//	double *dirweightsDD = R_Vector(C);
//	double *dirweightsDDs = R_Vector(C);

//	RprintIVecAsMat("nobsj", nobsj, 1, C);

	sumc=0;
	for(ii = 0; ii < C; ii++){
		sumc = sumc+nobsj[ii];

		tmp1 = tmp1 + dirweights[ii];
		tmp2 = tmp2 + lgamma(dirweights[ii]);

		tmp3 = tmp3 + (double) nobsj[ii] + dirweights[ii];
		tmp4 = tmp4 + lgamma( (double) nobsj[ii] + dirweights[ii]);

		tmp5 = tmp5 + 2*((double) nobsj[ii]) + dirweights[ii];
		tmp6 = tmp6 + lgamma( 2*((double) nobsj[ii]) + dirweights[ii]);
	}

//	Rprintf("tmp3 = %f\n", tmp3);
//	Rprintf("tmp5 = %f\n", tmp5);

	// The next piece of code is for the similarity likelihoodxprior/posterior

/*	for(ii=0;ii<C;ii++){
		sumc = sumc+nobsj[ii];
		pi_vec[ii] = 1/(double) C;
		dirweightsDD[ii] = (double) nobsj[ii] + dirweights[ii];
		dirweightsDDs[ii] = 2*((double) nobsj[ii]) + dirweights[ii];
	}



	ld1 = 0.0;
	for(ii=0;ii<C;ii++){
		ld1 = ld1 + nobsj[ii]*log(pi_vec[ii]);
	}
	ld2 = ddirich(pi_vec, dirweights, C, 1);
	ld3 = ddirich(pi_vec, dirweightsDD, C, 1);
	ld4 = ddirich(pi_vec, dirweightsDDs, C, 1);


	out1 = ld1 + ld2 - ld3;
	out2 = ld1 + ld3 - ld4;
*/

	out = (lgammafn(tmp1) - tmp2) + (tmp4 - lgammafn(tmp3));


//	Rprintf("out = %f\n", out);

	if(DD==1) out = (lgammafn(tmp3) - tmp4) + (tmp6 - lgammafn(tmp5));
	if(sumc==0) out = log(1);
	if(!logout) out = exp(out);
	return(out);

}




void PositivePart(double *x, double t,  int dim, int q, double* scratch){
/**************************************************************************************************
 * PURPOSE:
 * Evaluate the positive part function for a vector of values
 *
 * (x - t)^q_+	= (x-t)^qI[x > t]
 *
 * INPUTS:
 * x = input value
 * t = knot in a spline
 * dim = dimension of x
 * q = degree of positive part function
 * scratch = scratch vector to hold results
 *
 *************************************************************************************************/

    int ii;
	for(ii = 0; ii < dim; ii++){
		if(x[ii] > t){
			scratch[ii] = pow((x[ii] - t), q);
		} else {
			scratch[ii] = 0.0;
		}
	}

}




void psbf(double *x, double *knots, double* Bmat, int nobs,  int nknots, int q){
/**************************************************************************************************
 * PURPOSE:
 * Construct a q-degree spline basis
 *
 * (x - t)^q_+	= (x-t)^qI[x > t]
 *
 * INPUTS:
 * x = vector of data inputs (features
 * knots = vector containing location of knots in spline
 * Bmat = matrix that will hold the basis functions evaluated at values of x this
 *        is a contiguous array of memory whose dimension is nobs x ncols
 * nobs = dimension of x (number of observations)
 * ncols = number of basis functions (number of columns in basis matrix) = q + 1 + length(knots)
 * q = degree of positive part function
 *
 *************************************************************************************************/

    int ii,jj;
    int ncols = nknots+q+1;
	double *scr = R_Vector(nobs);

	for(ii = 0; ii < nobs; ii++){
//		Rprintf("ii = %d\n", ii);
		for(jj = 0; jj < ncols; jj++){
			if(jj < q + 1){
				Bmat[ii*ncols + jj] = pow(x[ii], jj);
			}
			if(jj >= (q+1)){
//				Rprintf("jj-q-1 = %d\n", jj-q-1);
//				Rprintf("knots[jj-q-1] = %f\n",knots[jj-q-1]);
				PositivePart(&x[ii], knots[jj-q-1], 1, q, scr);
//				RprintVecAsMat("scr", scr, 1, nobs);
				Bmat[ii*ncols + jj] = scr[0];
			}
		}
	}
}




double bs(int m, int bnum, int q, double xii, double *newknots){
/**************************************************************************************************
 * PURPOSE:
 * function called in bsb that carries out the recursive part of calculating a B-spline
 *
 *
 * INPUTS:
 * m = total number of knots (see next function m1 + 2*(q+1))
 * bnum = basis function number (this is actually iterated from 1 to m1 + q + 1)
 * q = degree of b-spine basis (up to cubic which is q=3)
 * xii = the iith entry of the x vector for which b-spline basis will be evaluated
 * nobs = number of observations in x
 * newknots = vector containing new knots
 *
 *************************************************************************************************/

	double result, result1, result2, temp1, temp2;

	if(q == 1){
		if(xii >= newknots[bnum] & xii <= newknots[bnum+1]){
			result = 1.0;
		}else{
			result = 0.0;
		}
	}else{
		result1 = bs(m, bnum, q-1, xii, newknots);
		temp1 = result1*((xii-newknots[bnum])/(newknots[bnum+q-1] - newknots[bnum]));


		if(temp1 != temp1) temp1=0.0;

		result2 = bs(m, bnum+1, q-1, xii, newknots);
		temp2 = result2*(newknots[bnum+q] - xii)/(newknots[bnum+q] - newknots[bnum+1]);


		if(temp2 != temp2) temp2=0.0;

		result = temp1 + temp2;



	}

	return(result);
}


void bsb(double *x, double *innerknots, double *Boundaryknots, double *basis,
         int nobs,  int m1, int q){
/**************************************************************************************************
 * PURPOSE:
 * Construct up to a cubic b-spline basis
 *
 * b_{j,1}(t) = (t - t_j)/(t_{j+1} - t_j) I[t_j <= t < t_{j+1}] +
 *              (t_{j+1} - t)/(t_{j+2} - t_{j+1}) i[t_{j+1} <= t < t_{j+2}]
 *
 *
 *
 * INPUTS:
 * x = vector of data inputs for which B-spline basis will be evaluated
 * inner.knots = vector containing location of inner knots in spline (length = m)
 * Boundar.knots = vector containing boundary knots (length = 2)
 * basis = a n x (m1+q+1) array of contiguous memory that will hold basis function evaluations
 * nobs = number of observations in x
 * m1 = number of inner knots
 * q = degree of b-spline (1 - linear, 2 - quadratic, 3 - cubic
 *
 *************************************************************************************************/

    int ii,jj;
    int m = m1 + 2*(q+1); /*total number of knots*/
	int k = m1 + q + 1;  /* number of spline basis */
	double *newknots = R_Vector(m);

//	RprintVecAsMat("innerknots =", innerknots, 1, m1);
//	RprintVecAsMat("Boundaryknots =", Boundaryknots, 1, 2);

	// create new vector of knots that contains boundary knots
//	Rprintf("q = %d\n", q);
	for(ii = 0; ii < m; ii++){
//		Rprintf("ii = %d\n", ii);
		if(ii < q+1) newknots[ii] = Boundaryknots[0];
		if(ii >= q+1 & ii < m1+q+1) newknots[ii] = innerknots[ii-(q+1)];
		if(ii >= m1 + q + 1) newknots[ii] = Boundaryknots[1];
	}
//	RprintVecAsMat("newknots =", newknots, 1, m);

	for(ii = 0; ii < nobs; ii++){

//		if(ii == 7 | ii == 8 | ii == 9) {
//			Rprintf("ii =================== %d\n", ii );
//			RprintVecAsMat("newknots =", newknots, 1, m);
//		}

		if(x[ii] == Boundaryknots[1]){

			basis[ii*k + k-1] = 1.0;

			for(jj = 0; jj < k-1; jj++){
				basis[ii*k + jj] = 0.0;
//				if(ii == 7 | ii == 8 | ii == 9) {Rprintf("basis = %f\n", basis[ii*k + jj]);}
			}

		}else{

			for(jj = 0; jj < k; jj++){

//				if(ii == 7 | ii == 8 | ii == 9) {
//					Rprintf("jj = %d\n", jj);
//					Rprintf("x[ii] = %f\n", x[ii]);
//					Rprintf("m = %d\n", m);
//					Rprintf("bs = %f\n", bs(m, jj, q+1, x[ii], newknots));
//				}

				basis[ii*k + jj] = bs(m, jj, q+1, x[ii], newknots);

//				if(ii == 7 | ii == 8 | ii == 9) {Rprintf("basis = %f\n", basis[ii*k + jj]);}
			}
		}
	}

//	RprintVecAsMat("basis", basis, nobs, k);
}





//////////////////////////////////////////////////////////////////////////////////
//
// Cohesion 1 (see Cohesion investigation file in sPPM latex file)
//
//////////////////////////////////////////////////////////////////////////////////

double Cohesion1(double *s1, double *s2, double epsilon, int dim, int lg){

	int ii;
	double cent1,cent2,sdist,out,maxdist,dist;
//	double sdist1, maxdist1, dist1;
	cent1=0, cent2=0;
//	Rprintf("dim = %d\n", dim);
	for(ii=0; ii<dim; ii++){

		cent1 = cent1 + s1[ii]/(double) dim;
		cent2 = cent2 + s2[ii]/(double) dim;

	}
//	Rprintf("cent1 = %f\n", cent1);
//	Rprintf("cent2 = %f\n", cent2);
	maxdist = sqrt((s1[0] - cent1)*(s1[0] - cent1) +
			       (s2[0] - cent2)*(s2[0] - cent2));
	sdist = 0;
	for(ii = 0; ii < dim; ii++){
		dist = sqrt((s1[ii] - cent1)*(s1[ii] - cent1) +
			        (s2[ii] - cent2)*(s2[ii] - cent2));

		sdist = sdist + dist;

		if(maxdist < dist) maxdist = dist;
	}

//	Rprintf("sdist = %f\n", sdist);
//	maxdist1 = 0;

//	Rprintf("maxdist1 = %f\n", maxdist1);
//	sdist1 = 0;
//	for(i=0; i<dim; i++){
//		Rprintf("i = %d\n", i);
//		for(ii=0; ii<i; ii++){
//			Rprintf("ii = %d\n", ii);
//			dist1 = sqrt((s1[i] - s1[ii])*(s1[i] - s1[ii]) +
//						 (s2[i] - s2[ii])*(s2[i] - s2[ii]));
//			sdist1 = sdist1 + dist1;
//			Rprintf("dist = %f\n", dist);
//			Rprintf("a = %f\n",a);
//			if(maxdist1 < dist1) maxdist1 = dist1;
//			Rprintf("maxdist1 = %f\n", maxdist1);

//		}
//	}
//	sdist=sdist;


//	Rprintf("Totaldist = %f\n", sdist);
//	Rprintf("Averagedist = %f\n", sdist/(double) dim);
//	Rprintf("maxdist = %f\n", maxdist);
//	Rprintf("maxdist1 = %f\n", maxdist1);

//	sdist = sdist/(double) dim;
//	sdist = maxdist1;
//	Rprintf("sdist = %f\n", sdist);

	if(sdist >= 1){
		out = -lgammafn(sdist*epsilon);
//		out = -lgammafn(10*sdist);
//		out = -log(sdist);
	}else if(sdist != 0){
//		out = -log(sdist + epsilon);
//		out = -log(1);
		out = -log(sdist);
	}else{
		out = log(1);
	}
//	Rprintf("out = %f\n", out);
	if(!lg) out = exp(out);
	return(out);

}

//////////////////////////////////////////////////////////////////////////////////
//
// Cohesion 2 (see Cohesion investigation file in sPPM latex file)
//
//////////////////////////////////////////////////////////////////////////////////

double Cohesion2(double *s1, double *s2, double a, int dim, int lg){
	int i,ii;
	double dist,out;
	for(i=0; i<dim; i++){
		for(ii=0; ii<dim; ii++){
			dist = sqrt((s1[i] - s1[ii])*(s1[i] - s1[ii]) +
						(s2[i] - s2[ii])*(s2[i] - s2[ii]));
//			Rprintf("dist = %f\n", dist);
			if(dist < a){
				out = 1;
			}else{
				out = 0.0;
				break;
			}

		}
		if(out<1) break;
	}

//	Rprintf("out = %f\n", out);
	if(lg) out = log(out);
	return(out);

}

//////////////////////////////////////////////////////////////////////////////////
//
// Cohesion 3 and 4 (see Cohesion investigation file in sPPM latex file)
//
//////////////////////////////////////////////////////////////////////////////////
double G2a(double a, int lg){
	double out;
	out = log(M_PI) + lgammafn(a) + lgammafn(a-0.5);
	if(!lg) out = exp(out);
	return(out);
}

double Cohesion3_4(double *s1, double *s2, double *mu0, double k0, double v0, double *L0,
					int dim, int Cohesion, int lg){

	int ii;
	double kn,knn,vn,vnn,out,sbar1,sbar2,dL0,dLn,dLnn;
	double s_sbar1,s_sbar2;
	double sbar_mu01, sbar_mu02, sbar_mun1, sbar_mun2, mun1, mun2, munn1, munn2;
	double Vs1, Vs2, Vs3, Vs4;
	double Vsbarmu01, Vsbarmu02, Vsbarmu03, Vsbarmu04;
	double Vsbarmun1, Vsbarmun2, Vsbarmun3, Vsbarmun4;
	double Ln1, Ln2, Ln3, Ln4, Lnn1, Lnn2, Lnn3, Lnn4;

	sbar1=0.0, sbar2=0.0;
	for(ii=0; ii<dim; ii++){

		sbar1 = sbar1 + s1[ii]/(double) dim;
		sbar2 = sbar2 + s2[ii]/(double) dim;

	}


	Vs1=0.0, Vs2=0.0, Vs3=0.0, Vs4=0.0;
	for(ii=0; ii<dim; ii++){

		s_sbar1 = s1[ii] - sbar1;
		s_sbar2 = s2[ii] - sbar2;


		Vs1 = Vs1 + s_sbar1*s_sbar1;
		Vs2 = Vs2 + s_sbar1*s_sbar2;
		Vs3 = Vs3 + s_sbar2*s_sbar1;
		Vs4 = Vs4 + s_sbar2*s_sbar2;

	}


	kn = k0 + dim; vn = v0 + dim;
	knn = kn + dim; vnn = vn + dim;

//	Rprintf("kn = %f\n", kn);
//	Rprintf("vn = %f\n", vn);

//	Rprintf("knn = %f\n", knn);
//	Rprintf("vnn = %f\n", vnn);

	mun1 = k0/(k0+dim)*mu0[0] + dim/(k0+dim)*sbar1;
	mun2 = k0/(k0+dim)*mu0[1] + dim/(k0+dim)*sbar2;
	munn1 = kn/(kn+dim)*mun1 + dim/(kn+dim)*sbar1;
	munn2 = kn/(kn+dim)*mun2 + dim/(kn+dim)*sbar2;


	sbar_mu01 = sbar1 - mu0[0];
	sbar_mu02 = sbar2 - mu0[1];
	sbar_mun1 = sbar1 - mun1;
	sbar_mun2 = sbar2 - mun2;


	Vsbarmu01 = sbar_mu01*sbar_mu01, Vsbarmun1 = sbar_mun1*sbar_mun1;
	Vsbarmu02 = sbar_mu01*sbar_mu02, Vsbarmun2 = sbar_mun1*sbar_mun2;
	Vsbarmu03 = sbar_mu02*sbar_mu01, Vsbarmun3 = sbar_mun2*sbar_mun1;
	Vsbarmu04 = sbar_mu02*sbar_mu02, Vsbarmun4 = sbar_mun2*sbar_mun2;

	Ln1 = L0[0] + Vs1 + k0*(dim)/(k0+dim)*Vsbarmu01;
	Ln2 = L0[1] + Vs2 + k0*(dim)/(k0+dim)*Vsbarmu02;
	Ln3 = L0[2] + Vs3 + k0*(dim)/(k0+dim)*Vsbarmu03;
	Ln4 = L0[3] + Vs4 + k0*(dim)/(k0+dim)*Vsbarmu04;

	Lnn1 = Ln1 + Vs1 + kn*(dim)/(kn+dim)*Vsbarmun1;
	Lnn2 = Ln2 + Vs2 + kn*(dim)/(kn+dim)*Vsbarmun2;
	Lnn3 = Ln3 + Vs3 + kn*(dim)/(kn+dim)*Vsbarmun3;
	Lnn4 = Ln4 + Vs4 + kn*(dim)/(kn+dim)*Vsbarmun4;


	dL0 = L0[0]*L0[3] - L0[1]*L0[2];
	dLn = Ln1*Ln4 - Ln2*Ln3;
	dLnn = Lnn1*Lnn4 - Lnn2*Lnn3;

//	RprintVecAsMat("L0", L0, 2, 2);
//	RprintVecAsMat("Ln", Ln, 2, 2);
//	RprintVecAsMat("Lnn", Lnn, 2, 2);



	if(Cohesion==3){
		out = -dim*log(M_PI) +
		      (G2a(0.5*vn, 1) - G2a(0.5*v0, 1)) +
		      (0.5*v0*log(dL0) - 0.5*vn*log(dLn)) +
		      (log(k0) - log(kn));
	}
	if(Cohesion==4){
//		Rprintf("G2a(0.5*vn) = %f\n", G2a(0.5*vn,1));
//		Rprintf("G2a(0.5*vnn) = %f\n", G2a(0.5*vnn,1));
		out = -dim*log(M_PI) +
		      (G2a(0.5*vnn,1) - G2a(0.5*vn,1)) +
		      (0.5*vn*log(dLn) - 0.5*vnn*log(dLnn)) +
		      (log(kn) - log(knn));
	}
//	Rprintf("out = %f\n", out);
	if(!lg) out = exp(out);
	return(out);

}

//////////////////////////////////////////////////////////////////////////////////
//
// Cohesion 5 (see Cohesion investigation file in sPPM latex file)
//
//////////////////////////////////////////////////////////////////////////////////

double Cohesion5(double *s1, double *s2, double phi, int dim, int lg){

	int ii;
	double cent1,cent2,sdist,out;
	cent1=0, cent2=0;
	for(ii=0; ii<dim; ii++){

		cent1 = cent1 + s1[ii]/(double) dim;
		cent2 = cent2 + s2[ii]/(double) dim;

	}
    //	Rprintf("cent1 = %f\n", cent1);
    //	Rprintf("cent2 = %f\n", cent2);
	sdist = 0;
	for(ii = 0; ii < dim; ii++){

        sdist = sdist + sqrt((s1[ii] - cent1)*(s1[ii] - cent1) +
                             (s2[ii] - cent2)*(s2[ii] - cent2));

	}
    //	Rprintf("sdist = %f\n", sdist);
	out = (-phi*sdist);
    //	Rprintf("out = %f\n", out);
	if(!lg) out = exp(out);
	return(out);

}


//////////////////////////////////////////////////////////////////////////////////
//
// Cohesion 6 (see Cohesion investigation file in sPPM latex file)
//
//////////////////////////////////////////////////////////////////////////////////

double Cohesion6(double *s1, double *s2, double phi, int dim, int lg){

	int ii;
	double cent1,cent2,sdist,out;
	cent1=0, cent2=0;
	for(ii=0; ii<dim; ii++){

		cent1 = cent1 + s1[ii]/(double) dim;
		cent2 = cent2 + s2[ii]/(double) dim;

	}
    //	Rprintf("cent1 = %f\n", cent1);
    //	Rprintf("cent2 = %f\n", cent2);
	sdist = 0;
	for(ii = 0; ii < dim; ii++){

        sdist = sdist + sqrt((s1[ii] - cent1)*(s1[ii] - cent1) +
                             (s2[ii] - cent2)*(s2[ii] - cent2));

	}

	out = -phi*log(sdist);
	if(dim==1) out = 0.0;

	if(!lg) out = exp(out);

	return(out);

}


void rPPMs(int cohesion, double M, int m, double *s1, double *s2, int *Si, int *nk, int *nh,
           double epsilon, double a, double *mu0, double k0, double v0, double *L0,
           double phi){
/**************************************************************************************************
 * Function that generates draws from a spatial Product Partition Model prior
 *
 * Inputs:
 *
 * cohesion - an integer indicating which cohesion function to use
 * M - DP scale parameter
 * m - integer that indicates dimension of s1 and s2
 * s1 - mx1 vector of 1st location coordinates (e.g., longitude)
 * s2 - mx1 vector of 2nd location coordinates (e.g., lattitude)
 * epsilon - double holding parameter for cohesion 1
 * a - double holding parameter for cohesion 2
 * mu0 - 2x1 double holding mean vector for cohesions 3,4
 * k0 - double holding variance scale parameter cohesions 3,4
 * v0 - double holding inverse-wishart scale parameter cohesions 3,4
 * L0 - 4x1 double holding matrix parameter for cohesions 3,4
 * phi - double holding tuning parameter for cohesions 5,6
 *
 * Outputs:
 * Si - nx1 scratch array of contiguous memory that holds partition of n objects
 * nk_aux - an integer indicating number of clusters
 * nh_aux - an integer indicating number of clusters
 *
 *************************************************************************************************/

	int i, ii, k;
	int  counter, iaux;
	double lCo, lCn, lC1, maxph, cprobh, denph, uu;
	double *s1tmp0 = R_Vector(m);
	double *s2tmp0 = R_Vector(m);
	double *s1tmp = R_Vector(m);
	double *s2tmp = R_Vector(m);
	double *ph = R_Vector(m);
	double *phtmp = R_Vector(m);
	double *probh = R_Vector(m);
	for(i=0;i<m;i++){
		Si[i]=0;
		nh[i]=0;
	}
	Si[0] = 1;
	nh[0] = 1;
	nk[0] = 1; // identifies the number of clusters
//	RprintVecAsMat("s1", s1, 1, m);
//	RprintVecAsMat("s2", s2, 1, m);
	for(i = 1; i < m; i++){
//		Rprintf("i = %d ====================\n", i);

		for(k=0; k < nk[0]; k++){
//			Rprintf("k = %d\n", k);
//			Rprintf("nh = %d\n", nh[k]);
//			Rprintf("nk[0] = %d\n", nk[0]);
			counter = 0;
			for(ii=0; ii<i; ii++){
//			 	Rprintf("ii = %d\n", ii);
//			 	Rprintf("Si[ii] = %d\n", Si[ii]);
				if(Si[ii] == k+1){
//					Rprintf("counter = %d\n", counter);
					s1tmp0[counter] = s1[ii];
					s2tmp0[counter] = s2[ii];
					s1tmp[counter] = s1[ii];
					s2tmp[counter] = s2[ii];
					counter = counter+1;

				}
			}

//			RprintVecAsMat("s1tmp0", s1tmp0, 1, nh[k]);
//			RprintVecAsMat("s2tmp0", s2tmp0, 1, nh[k]);

			s1tmp[nh[k]] = s1[i];
			s2tmp[nh[k]] = s2[i];

//			RprintVecAsMat("s1tmp", s1tmp, 1, nh[k]+1);
//			RprintVecAsMat("s2tmp", s2tmp, 1, nh[k]+1);

			if(cohesion==1){
				lCo = Cohesion1(s1tmp0, s2tmp0, epsilon, nh[k],1);
				lCn = Cohesion1(s1tmp, s2tmp, epsilon, nh[k]+1,1);
				lC1 = Cohesion1(&s1[i],&s2[i], epsilon, 1, 1);
			}

			if(cohesion==2){
				lCo = Cohesion2(s1tmp0, s2tmp0, a, nh[k],1);
				lCn = Cohesion2(s1tmp, s2tmp, a, nh[k]+1,1);
				lC1 = Cohesion2(&s1[i],&s2[i], a, 1, 1);
			}

			if(cohesion==3){
				lCo = Cohesion3_4(s1tmp0, s2tmp0, mu0, k0, v0, L0, nh[k], 3, 1);
				lCn = Cohesion3_4(s1tmp, s2tmp, mu0, k0, v0, L0, nh[k]+1,3, 1);
				lC1 = Cohesion3_4(&s1[i],&s2[i], mu0, k0, v0, L0, 1, 3, 1);

			}

			if(cohesion==4){

				lCo = Cohesion3_4(s1tmp0, s2tmp0, mu0, k0, v0, L0, nh[k], 4, 1);
				lCn = Cohesion3_4(s1tmp, s2tmp, mu0, k0, v0, L0, nh[k]+1,4, 1);
				lC1 = Cohesion3_4(&s1[i],&s2[i], mu0, k0, v0, L0, 1, 4, 1);

			}

			if(cohesion==5){

				lCo = Cohesion5(s1tmp0, s2tmp0, phi, nh[k],1);
				lCn = Cohesion5(s1tmp, s2tmp, phi, nh[k]+1,1);
				lC1 = Cohesion5(&s1[i],&s2[i], phi, 1, 1);

			}

			if(cohesion==6){

				lCo = Cohesion6(s1tmp0, s2tmp0, phi, nh[k],1);
				lCn = Cohesion6(s1tmp, s2tmp, phi, nh[k]+1,1);
				lC1 = Cohesion6(&s1[i],&s2[i], phi, 1, 1);

			}

			ph[k] = log(nh[k]) + lCn-lCo;
//			Rprintf("ph = %f\n", ph[k]);
		}

//		RprintVecAsMat("ph", ph, 1, nk[0]);


		ph[nk[0]] = log(M) + lC1;
		for(k = 0; k < nk[0]+1; k++) phtmp[k] = ph[k];
		R_rsort(phtmp,  nk[0]+1) ;

//		RprintVecAsMat("phtmp", phtmp, 1, nk[0]+1);

		maxph = phtmp[nk[0]];

//		Rprintf("maxph = %f\n", maxph);
//		RprintVecAsMat("ph", ph, 1, nk[0]+1);

		denph = 0.0;
		for(k = 0; k < nk[0]+1; k++){

			ph[k] = exp(ph[k] - maxph);
//					ph[k] = pow(exp(ph[k] - maxph), (1 - exp(-0.0001*(i+1))));
			denph = denph + ph[k];

		}

//		RprintVecAsMat("ph", ph, 1, nk[0]+1);

		for(k = 0; k < nk[0]+1; k++){

			probh[k] = ph[k]/denph;

		}
//		Rprintf("denph = %f\n", denph);

//		RprintVecAsMat("probh", probh, 1, nk[0]+1);

		uu = runif(0.0,1.0);
//		Rprintf("uu = %f\n", uu);

		cprobh= 0.0;

		for(k = 0; k < nk[0]+1; k++){

			cprobh = cprobh + probh[k];

			if (uu < cprobh){

				iaux = k+1;
				break;
			}
		}

		if(iaux <= nk[0]){

			Si[i] = iaux;
			nh[Si[i]-1] = nh[Si[i]-1] + 1;

		}else{

			nk[0] = nk[0] + 1;
			Si[i] = nk[0];
			nh[Si[i]-1] = 1;

		}
//		RprintIVecAsMat("Si", Si, 1, m);
//		Rprintf("nk[0] = %d\n", nk[0]);
//		RprintIVecAsMat("nh", nh, 1, nk[0]);
	}
}






