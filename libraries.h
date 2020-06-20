//Loading standard libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <malloc.h>
#include <omp.h>

/* GSL Libraries */

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "GLOBAL_CONSTANTS.h"

void Struc_params();

double Cpar_Ter(double x,double AC,double BC,double ABC);
double Cpar_QTer(double x,double y,double AD,double BD,double CD,double ABD,double BCD,double ACD);

double strain(double a_sub,double a_mat,double a_c,double C11,double C12,double C44);
double strain_bi(double a_sub,double a_mat,double a_c,double C11,double C12,double C44);

// GaInAs
double Eg_GaInAs(double x);
double VBO_GaInAs(double x);
double m_GaInAs(double x);
double a_GaInAs(double x);
double C11_GaInAs(double x);
double C12_GaInAs(double x);
double C44_GaInAs(double x);
double ac_GaInAs(double x);
double av_GaInAs(double x);

// AlInAs
double Eg_AlInAs(double x);
double VBO_AlInAs(double x);
double m_AlInAs(double x);
double a_AlInAs(double x);
double C11_AlInAs(double x);
double C12_AlInAs(double x);
double C44_AlInAs(double x);
double ac_AlInAs(double x);
double av_AlInAs(double x);

// AlGaAs
double bow_Eg_AlGaAs(double x);
double Eg_AlGaAs(double x);
double VBO_AlGaAs(double x);
double m_AlGaAs(double x);
double a_AlGaAs(double x);
double C11_AlGaAs(double x);
double C12_AlGaAs(double x);
double C44_AlGaAs(double x);
double ac_AlGaAs(double x);
double av_AlGaAs(double x);

// AlGaInAs
double Eg_AlGaInAs(double x);
double m_AlGaInAs(double x);
double a_AlGaInAs(double x);
double VBO_AlGaInAs(double x);
double ac_AlGaInAs(double x);
double av_AlGaInAs(double x); 
double C11_AlGaInAs(double x);
double C12_AlGaInAs(double x);
double C44_AlGaInAs(double x);

// QDot shape
double F_POT(double z);
double z_QDShape(double u);

// Matrix elements

struct f_params;

int Delta(int x,int y);
double kzero(int l,int m);

double Z_n(double z,int n,double L);
double Rfig_lm(double u, int l, int m);
double modFz_Onda(double z,double L,int l,int k);
double R_lm(double u, int l, int m);
double DR_lm(double u, int l, int m);
double Sn1n2(int n1,int n2,double s,double L);
double ASn1n2(int n1,int n2,double s,double L);
double Cn1n2(int n1,int n2,double s,double L);
double Bn1n2(int n1,int n2,double s,double L);

// QDOT Integrals
double QDOT_mass_fint1(double u,void *params);
double QDOT_mass_fint2(double u,void *params);
double QDOT_mass_fint3(double u,void *params);

double QDOT_mass_type1(int l1,int m1,int n1,int l2,int m2,int n2,double a,double b);
double QDOT_mass_type2(int l1,int m1,int n1,int l2,int m2,int n2,double a,double b);
double QDOT_mass_type3(int l1,int m1,int n1,int l2,int m2,int n2,double a,double b);

// QWell Integrals
double QWell_a_finteg(double u,void *params);
double QWell_a_integral(int l1,int m1,int n1,int l2,int m2,int n2,double a,double b);
double QWell_b_finteg(double u,void *params);
double QWell_b_integral(int l1,int m1,int n1,int l2,int m2,int n2,double a,double b);

double QDOT_bias_finteg(double u,void *params);
double QDOT_bias_integral(int l1,int m1,int n1,int l2,int m2,int n2,double a,double b);

double f1_Zn1Zn2(int n1,int n2,double u);
double f2_Zn1Zn2(int n1,int n2,double u);
double Integrando_RRB_lm(double u,void *params);
double INTEGRAL_RRBlm(int l1,int m1,int n1,int l2,int m2,int n2,double a,double b);

double VL(double z,double z_a,double z_b,double x_a,double x_b);
double VC(double z,double z_a,double z_b,double x_a,double x_b);

double linear_V(int l1,int m1,int n1,int l2,int m2,int n2,double z_a,double z_b,double L,double x_a,double x_b);
double Quadratic_V(int l1,int m1,int n1,int l2,int m2,int n2,double z_a,double z_b,double L,double x_a,double x_b);


double AELEMENT(int l1,int m1,int n1,int l2,int m2,int n2);

void map(double *z,int Nm);
double _Complex F_Onda(double r,double phi,double z,int l,int k);
void Lapackdiag(double _Complex *A, int dim,double _Complex *evals,double _Complex *evecs);

void arraysorted(double *a,int length,int *indx);
