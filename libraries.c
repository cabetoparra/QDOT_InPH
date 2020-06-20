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

/* initialization parameters */

#include "GLOBAL_CONSTANTS.h"
#include "libraries.h"

// Reading Structure
void Struc_params()
{
  FILE *al,*lg;
  int i,j,Nm;
  char MAT[100][100],Pos[100][100],dz[100][100],x[100][100],y[100][100],type[100][100],col[100],file[100];
  double mat[100][100],zz[50],u,z,CB;
  
  if(strcmp(HOLE,"si")==0){
    sprintf(col,"vb_");
    CB = 0.0;
  }
  else{
    sprintf(col,"cb_");
    CB = 1.0;
  }  
  
  al=fopen("../param.txt","r");
  i=0;
  do{
    fscanf(al,"%s %s %s %s %s %s",&Pos[i],&MAT[i],&dz[i],&x[i],&y[i],&type[i]);
    i++;
  }while(!feof(al));
  fclose(al);
  
  Nm=i-1;

  for(j=0;j<Nm;j++)
    {
      if(strcmp(MAT[j],"InAs")==0)
	{
	  mat[atoi(Pos[j])][0] = CB * Eg_InAs+VBO_InAs+strain(a_subs,a_InAs,ac_InAs,C11_InAs,C12_InAs,C44_InAs);
	  mat[atoi(Pos[j])][1] = m_InAs;
	  mat[atoi(Pos[j])][2] = mhh_100_InAs;
	  mat[atoi(Pos[j])][3] = mhh_110_InAs;
	  mat[atoi(Pos[j])][4] = atof(dz[j]);
	  
	  // Type of slab
	  if(strcmp(type[j],"Qdot")==0)
	    mat[atoi(Pos[j])][5] = 1;
	  else
	    mat[atoi(Pos[j])][5] = 0;
	}
      else if(strcmp(MAT[j],"InP")==0)
	{
	  mat[atoi(Pos[j])][0] = CB * Eg_InP+VBO_InP+strain(a_subs,a_InP,ac_InP,C11_InP,C12_InP,C44_InP);
	  mat[atoi(Pos[j])][1] = m_InP;
	  mat[atoi(Pos[j])][2] = mhh_100_InP;
	  mat[atoi(Pos[j])][3] = mhh_110_InP;
	  mat[atoi(Pos[j])][4] = atof(dz[j]); 	  

	  // Type of slab
	  if(strcmp(type[j],"Qdot")==0)
	    mat[atoi(Pos[j])][5] = 1;
	  else
	    mat[atoi(Pos[j])][5] = 0;
	  
	}
      else if(strcmp(MAT[j],"AlAs")==0)
	{
	  mat[atoi(Pos[j])][0] = CB * Eg_AlAs+VBO_AlAs+strain(a_subs,a_AlAs,ac_AlAs,C11_AlAs,C12_AlAs,C44_AlAs);
	  mat[atoi(Pos[j])][1] = m_AlAs;
	  mat[atoi(Pos[j])][2] = mhh_100_AlAs;
	  mat[atoi(Pos[j])][3] = mhh_110_AlAs;
	  mat[atoi(Pos[j])][4] = atof(dz[j]); 	  

	  // Type of slab
	  if(strcmp(type[j],"Qdot")==0)
	    mat[atoi(Pos[j])][5] = 1;
	  else
	    mat[atoi(Pos[j])][5] = 0;
	}
      else if(strcmp(MAT[j],"GaAs")==0)
	{
	  mat[atoi(Pos[j])][0] = CB * Eg_GaAs+VBO_GaAs+strain(a_subs,a_GaAs,ac_GaAs,C11_GaAs,C12_GaAs,C44_GaAs);
	  mat[atoi(Pos[j])][1] = m_GaAs;
	  mat[atoi(Pos[j])][2] = mhh_100_GaAs;
	  mat[atoi(Pos[j])][3] = mhh_110_GaAs;
	  mat[atoi(Pos[j])][4] = atof(dz[j]); 	  

	  // Type of slab
	  if(strcmp(type[j],"Qdot")==0)
	    mat[atoi(Pos[j])][5] = 1;
	  else
	    mat[atoi(Pos[j])][5] = 0;
	}
      else if(strcmp(MAT[j],"GaInAs")==0)
	{
	  mat[atoi(Pos[j])][0] = CB * Eg_GaInAs(atof(x[j]))+VBO_GaInAs(atof(x[j]))+strain(a_subs,a_GaInAs(atof(x[j])),ac_GaInAs(atof(x[j])),C11_GaInAs(atof(x[j])),C12_GaInAs(atof(x[j])),C44_GaInAs(atof(x[j])));
	  mat[atoi(Pos[j])][1]=m_GaInAs(atof(x[j]));
	  mat[atoi(Pos[j])][2]=m_GaInAs(atof(x[j]));
	  mat[atoi(Pos[j])][3]=m_GaInAs(atof(x[j]));
	  mat[atoi(Pos[j])][4]=atof(dz[j]); 	  

	  // Type of slab
	  if(strcmp(type[j],"Qdot")==0)
	    mat[atoi(Pos[j])][5] = 1;
	  else
	    mat[atoi(Pos[j])][5] = 0;
	}
      else if(strcmp(MAT[j],"AlInAs")==0)
	{
	  mat[atoi(Pos[j])][0] = CB * Eg_AlInAs(atof(x[j]))+VBO_AlInAs(atof(x[j]))+strain(a_subs,a_AlInAs(atof(x[j])),ac_AlInAs(atof(x[j])),C11_AlInAs(atof(x[j])),C12_AlInAs(atof(x[j])),C44_AlInAs(atof(x[j])));
	  mat[atoi(Pos[j])][1]=m_AlInAs(atof(x[j]));
	  mat[atoi(Pos[j])][2]=m_AlInAs(atof(x[j]));
	  mat[atoi(Pos[j])][3]=m_AlInAs(atof(x[j]));
	  mat[atoi(Pos[j])][4]=atof(dz[j]); 	 

	  // Type of slab
	  if(strcmp(type[j],"Qdot")==0)
	    mat[atoi(Pos[j])][5] = 1;
	  else
	    mat[atoi(Pos[j])][5] = 0;
	}
      else if(strcmp(MAT[j],"AlGaAs")==0)
	{
	  mat[atoi(Pos[j])][0] = CB * Eg_AlGaAs(atof(x[j]))+VBO_AlGaAs(atof(x[j]))+strain(a_subs,a_AlGaAs(atof(x[j])),ac_AlGaAs(atof(x[j])),C11_AlGaAs(atof(x[j])),C12_AlGaAs(atof(x[j])),C44_AlGaAs(atof(x[j])));
	  mat[atoi(Pos[j])][1]=m_AlGaAs(atof(x[j]));
	  mat[atoi(Pos[j])][2]=m_AlGaAs(atof(x[j]));
	  mat[atoi(Pos[j])][3]=m_AlGaAs(atof(x[j]));
	  mat[atoi(Pos[j])][4]=atof(dz[j]); 	  

	  // Type of slab
	  if(strcmp(type[j],"Qdot")==0)
	    mat[atoi(Pos[j])][5] = 1;
	  else
	    mat[atoi(Pos[j])][5] = 0;
	}
      else if(strcmp(MAT[j],"AlGaInAs")==0)
	{
	  mat[atoi(Pos[j])][0] = CB * Eg_AlGaInAs(atof(x[j]))+VBO_AlGaInAs(atof(x[j]))+strain(a_subs,a_AlGaInAs(atof(x[j])),ac_AlGaInAs(atof(x[j])),C11_AlGaInAs(atof(x[j])),C12_AlGaInAs(atof(x[j])),C44_AlGaInAs(atof(x[j])));
	  mat[atoi(Pos[j])][1]=m_AlGaInAs(atof(x[j]));
	  mat[atoi(Pos[j])][2]=m_AlGaInAs(atof(x[j]));
	  mat[atoi(Pos[j])][3]=m_AlGaInAs(atof(x[j]));
	  mat[atoi(Pos[j])][4]=atof(dz[j]); 	  

	  // Type of slab
	  if(strcmp(type[j],"Qdot")==0)
	    mat[atoi(Pos[j])][5] = 1;
	  else
	    mat[atoi(Pos[j])][5] = 0;
	}      
    }
  
  sprintf(file,"../datafile-%s/%sparameters.dat",Q_Shape,col);
  lg=fopen(file,"w");
  for(i=0;i<Nm-2;i++)
    fprintf(lg,"%7.4lf\t%7.4lf\t%7.4lf\t%7.4lf\t%7.1lf\t%d\n",mat[i][0],mat[i][1],mat[i][2],mat[i][3],mat[i][4],(int)rint(mat[i][5]));
  fclose(lg);
  
  zz[0]=0;
  for(i=0;i<Nm-2;i++)
    {
      zz[i+1]=zz[i]+mat[i][4];
    }

  if(strcmp(what2do,"diag") == 0){
    printf(">> Sample length = %0.2lf nm\n",(zz[Nm-2])/10);
    printf(">> Potential... done!\n");

    sprintf(file,"../datafile-%s/%sPotential.dat",Q_Shape,col);
    lg=fopen(file,"w");
    for(z=(zz[0]+0.00001);z<=(zz[Nm-2]);z+=zz[Nm-2]/10000)
      fprintf(lg,"%lf %lf\n",z/10,F_POT(z));
    fclose(lg);

    sprintf(file,"../datafile-%s/Qd_Shape.dat",Q_Shape);
    lg=fopen(file,"w");
    for(u=0.0;u <= 1;u+=1.0/1000.0)
      fprintf(lg,"%lf %lf\n",u*Rs/10,z_QDShape(u)/10);
    fclose(lg);
  }
  else; 
}

double strain(double a_sub,double a_mat,double a_c,double C11,double C12,double C44)
{
  double eparal,ezz,y;
  eparal=(a_sub-a_mat)/a_mat;
  ezz=-2*(C12/C11)*eparal;//Uniaxial strain
  return a_c*(2*eparal+ezz);
}

double strain_bi(double a_sub,double a_mat,double a_c,double C11,double C12,double C44)
{
  double eparal,ezz,y;
  eparal=(a_sub-a_mat)/a_mat;
  ezz=-((2*C11+4*C12-4*C44)/(C11+2*C12+4*C44))*eparal;//Biaxial strain
  
  return a_c*(2*eparal+ezz);
}

//******** STRUCTURE ROUTINES

double Cpar_Ter(double x,double AC,double BC,double ABC)
{
  double y;
  return y=(1-x)*AC+x*BC-x*(1-x)*ABC;
}

double Cpar_QTer(double x,double y,double AD,double BD,double CD,double ABD,double BCD,double ACD)
{
  double g,u,v,w;

  u=(1-x-y)/2;
  v=(2-x-2*y)/2;
  w=(2-2*x-y)/2;
  
  g=(x*y*Cpar_Ter(u,AD,BD,ABD)+y*(1-x-y)*Cpar_Ter(v,BD,CD,BCD)+(1-x-y)*x*Cpar_Ter(w,AD,CD,ACD))/(x*y+y*(1-x-y)+(1-x-y)*x);

  return g;
}


//----- TERNARIES PARAMETERS

//--- GaInAs

double Eg_GaInAs(double x){return Cpar_Ter(x,Eg_GaAs,Eg_InAs,bow_Eg_GaInAs);}
double VBO_GaInAs(double x){return Cpar_Ter(x,VBO_GaAs,VBO_InAs,bow_VBO_GaInAs);}
double m_GaInAs(double x){return Cpar_Ter(x,m_GaAs,m_InAs,bow_m_GaInAs);}
double a_GaInAs(double x){return Cpar_Ter(x,a_GaAs,a_InAs,bow_a_GaInAs);}
double C11_GaInAs(double x){return Cpar_Ter(x,C11_GaAs,C11_InAs,bow_C11_GaInAs);}
double C12_GaInAs(double x){return Cpar_Ter(x,C12_GaAs,C12_InAs,bow_C12_GaInAs);}
double C44_GaInAs(double x){return Cpar_Ter(x,C44_GaAs,C44_InAs,bow_C44_GaInAs);}
double ac_GaInAs(double x){return Cpar_Ter(x,ac_GaAs,ac_InAs,bow_ac_GaInAs);}
double av_GaInAs(double x){return Cpar_Ter(x,av_GaAs,av_InAs,bow_av_GaInAs);}

//--- AlInAs

double Eg_AlInAs(double x){return Cpar_Ter(x,Eg_InAs,Eg_AlAs,bow_Eg_AlInAs);}
double VBO_AlInAs(double x){return Cpar_Ter(x,VBO_InAs,VBO_AlAs,bow_VBO_AlInAs);}
double m_AlInAs(double x){return Cpar_Ter(x,m_InAs,m_AlAs,bow_m_AlInAs);}
double a_AlInAs(double x){return Cpar_Ter(x,a_InAs,a_AlAs,bow_a_AlInAs);}
double C11_AlInAs(double x){return Cpar_Ter(x,C11_InAs,C11_AlAs,bow_C11_AlInAs);}
double C12_AlInAs(double x){return Cpar_Ter(x,C12_InAs,C12_AlAs,bow_C12_AlInAs);}
double C44_AlInAs(double x){return Cpar_Ter(x,C44_InAs,C44_AlAs,bow_C44_AlInAs);}
double ac_AlInAs(double x){return Cpar_Ter(x,ac_InAs,ac_AlAs,bow_ac_AlInAs);}
double av_AlInAs(double x){return Cpar_Ter(x,av_InAs,av_AlAs,bow_av_AlInAs);}

//--- AlGaAs

double bow_Eg_AlGaAs(double x){return -0.127+1.310*x;}
double Eg_AlGaAs(double x){return Cpar_Ter(x,Eg_GaAs,Eg_AlAs,bow_Eg_AlGaAs(x));}
double VBO_AlGaAs(double x){return Cpar_Ter(x,VBO_GaAs,VBO_AlAs,bow_VBO_AlGaAs);}
double m_AlGaAs(double x){return Cpar_Ter(x,m_GaAs,m_AlAs,bow_m_AlGaAs);}
double a_AlGaAs(double x){return Cpar_Ter(x,a_GaAs,a_AlAs,bow_a_AlGaAs);}
double C11_AlGaAs(double x){return Cpar_Ter(x,C11_GaAs,C11_AlAs,bow_C11_AlGaAs);}
double C12_AlGaAs(double x){return Cpar_Ter(x,C12_GaAs,C12_AlAs,bow_C12_AlGaAs);}
double C44_AlGaAs(double x){return Cpar_Ter(x,C44_GaAs,C44_AlAs,bow_C44_AlGaAs);}
double ac_AlGaAs(double x){return Cpar_Ter(x,ac_GaAs,ac_AlAs,bow_ac_AlGaAs);}
double av_AlGaAs(double x){return Cpar_Ter(x,av_GaAs,av_AlAs,bow_av_AlGaAs);}


//----- QUATERNARIES PARAMETERS

//--- AlGaInAs

double Int_QTer(double x,double GA,double GB,double CAB)
{
  return (1-x/0.48)*GA+(x/0.48)*GB-(x/0.48)*(1-x/0.48)*CAB;
}

double Eg_AlGaInAs(double x){return Int_QTer(x,Eg_GaInAs(0.53),Eg_AlInAs(0.48),0.22);}
double m_AlGaInAs(double x){return Int_QTer(x,m_GaInAs(0.53),m_AlInAs(0.48),0.0);}
double a_AlGaInAs(double x){return Int_QTer(x,a_GaInAs(0.53),a_AlInAs(0.48),0.0);}
double VBO_AlGaInAs(double x){return Int_QTer(x,VBO_GaInAs(0.53),VBO_AlInAs(0.48),0.0);}
double ac_AlGaInAs(double x){return Int_QTer(x,ac_GaInAs(0.53),ac_AlInAs(0.48),0.0);}
double av_AlGaInAs(double x){return Int_QTer(x,av_GaInAs(0.53),av_AlInAs(0.48),0.0);}
double C11_AlGaInAs(double x){return Int_QTer(x,C11_GaInAs(0.53),C11_AlInAs(0.48),0.0);}
double C12_AlGaInAs(double x){return Int_QTer(x,C12_GaInAs(0.53),C12_AlInAs(0.48),0.0);}
double C44_AlGaInAs(double x){return Int_QTer(x,C44_GaInAs(0.53),C44_AlInAs(0.48),0.0);}


double F_POT(double z)
{
  double y,zz[100],V[100],dm[100],dhz[100],dhp[100],dz[100],tp[100];
  FILE *gl;
  int i,Nm,j;

  char col[100],file[100];
  
  if(strcmp(HOLE,"si")==0)
    sprintf(col,"vb_");
  else
    sprintf(col,"cb_");
  sprintf(file,"../datafile-%s/%sparameters.dat",Q_Shape,col);
  
  gl=fopen(file,"r");
  i=0;
  do{
    fscanf(gl,"%lf %lf %lf %lf %lf %lf",&V[i],&dm[i],&dhz[i],&dhp[i],&dz[i],&tp[i]);
    i++;
  }while(!feof(gl));
  fclose(gl);
  
  Nm=i-1;

  zz[0]=0.0;
  for(i=0;i<=Nm;i++)
    {
      zz[i+1]=(zz[i]+dz[i]);
    }
  
  for(j=0;j<=Nm;j++)
    {
      if(z>=zz[j] &&z<=zz[j+1])
	{
	  if(BIAS==0)
	    {
	      if(z>=zz[V1_pos]*1.0&&z<=zz[V1_pos+1])
		{	
		  if(type_V==1)
		    y=VL(z,zz[V1_pos],zz[V1_pos+1],x_con_a1,x_con_b1);
		  else if(type_V==2)
		    y=VC(z,zz[V1_pos],zz[V1_pos+1]*1.0,x_con_a1,x_con_b1);
		}
	      else if(z>=zz[V2_pos] &&z<=zz[V2_pos+1])
		{	
		  if(type_V==1)
		    y=VL(z,zz[V2_pos],zz[V2_pos+1],x_con_a2,x_con_b2);
		  else if(type_V==2)
		    y=VC(z,zz[V2_pos],zz[V2_pos+1],x_con_a2,x_con_b2);
		}
	      else
		y=V[j];
	    }
	  else if(BIAS==1)
	    {
	      if(z>=zz[V1_pos] && z<=zz[V1_pos+1])
		{	
		  if(type_V==1)
		    y=E_field*z+VL(z,zz[V1_pos],zz[V1_pos+1],x_con_a1,x_con_b1);
		  else if(type_V==2)
		    y=E_field*z+VC(z,zz[V1_pos],zz[V1_pos+1],x_con_a1,x_con_b1);
		}
	      else if(z>=zz[V2_pos]&&z<=zz[V2_pos+1])
		{	
		  if(type_V==1)
		    y=E_field*z+VL(z,zz[V2_pos]*1.0,zz[V2_pos+1]*1.0,x_con_a2,x_con_b2);
		  else if(type_V==2)
		    y=E_field*z+VC(z,zz[V1_pos]*1.0,zz[V2_pos+1]*1.0,x_con_a2,x_con_b2);
		}
	      else
		y=V[j]+E_field * z + depl_e * exp(-z/depl_L);
	    }
	}
    }
  return y;
}

//--
//---- ATIPICAL POTENTIALS...

//----------------------------------------------------------------------------------------------
//------- INTEGRALS ALONG Z-DIRECTION WITH Z^1 -> f1_Zn1n2  E Z^2 -> f2_Zn1n2

double f1_Zn1Zn2(int n1,int n2,double u)
{
  double y;

  if(n1==n2)
    y=2*(pow(u/2,2)-cos(2*n1*M_PI*u)/(8*pow(n1*M_PI,2))-u*sin(2*n1*M_PI*u)/(4*M_PI*n1));
  else 
    y=((cos((n2-n1)*M_PI*u)+(n2-n1)*u*M_PI*sin((n2-n1)*M_PI*u))/pow(n2-n1,2)-(cos((n2+n1)*M_PI*u)+(n2+n1)*u*M_PI*sin((n2+n1)*M_PI*u))/pow(n2+n1,2))/(pow(M_PI,2));
  
  return y;
}

double f2_Zn1Zn2(int n1,int n2,double u)//multiply by L^2 in Matrix Z element
{
  double y;

  if(n1==n2)
    y=2*(pow(u,3)/6-u*cos(2*n1*M_PI*u)/(4*pow(n1*M_PI,2))-(-1+2*pow(n1*M_PI*u,2))*sin(2*n1*M_PI*u)/(8*pow(M_PI*n1,3)));
  else 
    y=2*((2*(n2-n1)*M_PI*u*cos((n2-n1)*M_PI*u)+(-2+pow(n2*M_PI*u,2)-2*n1*n2*pow(M_PI*u,2)+pow(n1*M_PI*u,2))*sin((n2-n1)))/pow(n2-n1,3)-(2*(n2+n1)*M_PI*u*cos((n2+n1)*M_PI*u)+(-2+pow(n2*M_PI*u,2)-2*n1*n2*pow(M_PI*u,2)+pow(n1*M_PI*u,2))*sin((n2+n1)))/pow(n2+n1,3))/(2*pow(M_PI,3));
  
  return y;
}


double VL(double z,double z_a,double z_b,double x_a,double x_b)// LINEARS... TO RAMPS
{
  double yy,V_a,V_b;

  V_a=Eg_AlGaInAs(x_a)+VBO_AlGaInAs(x_a)+strain(a_InP,a_AlGaInAs(x_a),ac_AlGaInAs(x_a),C11_AlGaInAs(x_a),C12_AlGaInAs(x_a),C44_AlGaInAs(x_a));
  V_b=Eg_AlGaInAs(x_b)+VBO_AlGaInAs(x_b)+strain(a_InP,a_AlGaInAs(x_b),ac_AlGaInAs(x_b),C11_AlGaInAs(x_b),C12_AlGaInAs(x_b),C44_AlGaInAs(x_b));

 
  if(V_a<V_b)
    yy=((V_b-V_a)/(z_b-z_a))*(z-z_a)+V_a;
  else
    yy=((V_b-V_a)/(z_b-z_a))*(z-z_b)+V_b;

  return yy;
}

double VC(double z,double z_a,double z_b,double x_a,double x_b)// QUADRATICS... TO PARABOLS
{
  double y,V_a,V_b,a;

  V_a=Eg_AlGaInAs(x_a)+VBO_AlGaInAs(x_a)+strain(a_InP,a_AlGaInAs(x_a),ac_AlGaInAs(x_a),C11_AlGaInAs(x_a),C12_AlGaInAs(x_a),C44_AlGaInAs(x_a));
  V_b=Eg_AlGaInAs(x_b)+VBO_AlGaInAs(x_b)+strain(a_InP,a_AlGaInAs(x_b),ac_AlGaInAs(x_b),C11_AlGaInAs(x_b),C12_AlGaInAs(x_b),C44_AlGaInAs(x_b));
  
  if(V_a<V_b)
    {
      a=((V_b-V_a)/pow(z_b-z_a,2));
      y=((V_b-V_a)/pow(z_b-z_a,2))*pow(z-z_a,2)+V_a;
    }
  else
    {
      a=((V_a-V_b)/pow(z_a-z_b,2));
      y=((V_a-V_b)/pow(z_a-z_b,2))*pow(z-z_b,2)+V_b;
    }
  return y;
}

//--- INTEGRALS OF ATIPICAL POTENTIALS

double linear_V(int l1,int m1,int n1,int l2,int m2,int n2,double z_a,double z_b,double L,double x_a,double x_b)
{
  double y,V_a,V_b;

  V_a=Eg_AlGaInAs(x_a)+VBO_AlGaInAs(x_a)+strain(a_InP,a_AlGaInAs(x_a),ac_AlGaInAs(x_a),C11_AlGaInAs(x_a),C12_AlGaInAs(x_a),C44_AlGaInAs(x_a));
  V_b=Eg_AlGaInAs(x_b)+VBO_AlGaInAs(x_b)+strain(a_InP,a_AlGaInAs(x_b),ac_AlGaInAs(x_b),C11_AlGaInAs(x_b),C12_AlGaInAs(x_b),C44_AlGaInAs(x_b));


  if(V_a<V_b)
    y=(V_a-z_a*((V_b-V_a)/(z_b-z_a)))*(Sn1n2(n1,n2,z_b,L)-Sn1n2(n1,n2,z_a,L))+((V_b-V_a)/(z_b-z_a))*(f1_Zn1Zn2(n1,n2,z_b/L)-f1_Zn1Zn2(n1,n2,z_a/L))*(2*L);
  else
    y=(V_b-z_b*((V_b-V_a)/(z_b-z_a)))*(Sn1n2(n1,n2,z_b,L)-Sn1n2(n1,n2,z_a,L))+((V_b-V_a)/(z_b-z_a))*(f1_Zn1Zn2(n1,n2,z_b/L)-f1_Zn1Zn2(n1,n2,z_a/L))*(2*L);

  return y/Ry;
}

double Quadratic_V(int l1,int m1,int n1,int l2,int m2,int n2,double z_a,double z_b,double L,double x_a,double x_b)
{
  double y,V_a,V_b,a;

  V_a=Eg_AlGaInAs(x_a)+VBO_AlGaInAs(x_a)+strain(a_InP,a_AlGaInAs(x_a),ac_AlGaInAs(x_a),C11_AlGaInAs(x_a),C12_AlGaInAs(x_a),C44_AlGaInAs(x_a));
  V_b=Eg_AlGaInAs(x_b)+VBO_AlGaInAs(x_b)+strain(a_InP,a_AlGaInAs(x_b),ac_AlGaInAs(x_b),C11_AlGaInAs(x_b),C12_AlGaInAs(x_b),C44_AlGaInAs(x_b));

  if(V_a<V_b)
    {
      a=((V_b-V_a)/pow(z_b-z_a,2));
      y=(V_a+a*pow(z_a,2))*(Sn1n2(n1,n2,z_b,L)-Sn1n2(n1,n2,z_a,L))+2*a*z_a*(f1_Zn1Zn2(n1,n2,z_b/L)-f1_Zn1Zn2(n1,n2,z_a/L))*(L)+pow(a,2)*(f2_Zn1Zn2(n1,n2,z_b/L)-f2_Zn1Zn2(n1,n2,z_a/L))*pow(L,2);
    }
  else
    {
      a=((V_a-V_b)/pow(z_b-z_a,2));
      y=(V_b+a*pow(z_b,2))*(Sn1n2(n1,n2,z_b,L)-Sn1n2(n1,n2,z_a,L))+2*a*z_b*(f1_Zn1Zn2(n1,n2,z_b/L)-f1_Zn1Zn2(n1,n2,z_a/L))*(L)+pow(a,2)*(f2_Zn1Zn2(n1,n2,z_b/L)-f2_Zn1Zn2(n1,n2,z_a/L))*pow(L,2);
    }
  return y/Ry;
}


//--- QD - SHAPE-FUNCTION

double z_QDShape(double u)
{
  double y;

  // Cylinder
  if(strcmp(Q_Shape,"Cyl")==0){
    if(u >= 0.0 && u <= Ro/Rs)
      y=ho;
    else
      y=0;
  }
  // Trapezoide
  else if(strcmp(Q_Shape,"Trap")==0){
    
    if(u>=0.0 && u<=ro/Rs)
      y = ho;
    else if(u>=ro/Rs && u<Ro/Rs)
      y = (ho/(ro-Ro))*(u*Rs-ro)+ho;
    else
      y = 0;
  }
  // Cone
  else if(strcmp(Q_Shape,"Cone")==0){
    if(u>=0.0 && u<=Ro/Rs)
      y = ho * (1.0-u * (Rs/Ro));
    else
      y = 0;
  }
  else if(strcmp(Q_Shape,"Lens")==0){
    if(u>=0.0 && u<=Ro/Rs)
      y = ho * sqrt(1.0-pow(u * (Rs/Ro),2));
    else
      y = 0;
  }
  
  return y;
}


//--- ROUTINES FOR DIAGONALIZATION PROCESS

struct f_params
{
  int l1,m1,n1,l2,m2,n2;
};

int Delta(int x,int y)
{
  int g;

  if(x==y)
    g=1;
  else
    g=0;
  
  return g;
}

double kzero(int l,int m)
{
  double y=gsl_sf_bessel_zero_Jnu(l,m);
  return y;
}

//--- RADIAL FUNCTIONS AND ITS DERIVATES 

double R_lm(double u, int l,int m)
{
  double y;

  y=(sqrt(2.0)/(1.0*Rs*gsl_sf_bessel_Jn(l+1,kzero(l,m))))*gsl_sf_bessel_Jn(l,kzero(l,m)*(u));

  return y;
}

double DR_lm(double u, int l, int m)
{
  double y;

  y=(1.0/(sqrt(2.0)*Rs*gsl_sf_bessel_Jn(l+1,kzero(l,m))))*(gsl_sf_bessel_Jn(l-1,kzero(l,m)*(u))-gsl_sf_bessel_Jn(l+1,kzero(l,m)*(u)))*kzero(l,m);

  return y;
}

//--- ANTIDERIVATES ALONG Z-DIRECTION

double Sn1n2(int n1,int n2,double s,double L)
{
  double y;
  
  if(n1==n2)
    y = (2.0/L) * (s/2.0-sin(2*n1*M_PI*(s))/(4 * M_PI *n1));
  else
    y = (2.0/L) * (sin((n1-n2)*M_PI*(s))/(n1-n2)-sin((n1+n2)*M_PI*(s))/(n1+n2))/(2 * M_PI);
  
  return y;
}

double Cn1n2(int n1,int n2,double s,double L)
{
  double y;
  
  if(n1==n2)
    y = (2.0/L) * pow(n1*M_PI,2) *(s/2.0+sin(2*n1*M_PI*(s))/(4 * M_PI * n1));
  else
    y = (2.0/L) * (n1 * n2 * pow(M_PI,2))*(sin((n1-n2)*M_PI*(s))/(n1-n2)+sin((n1+n2)*M_PI*(s))/(n1+n2))/(2 * M_PI);
  
  return y;
}

double Bn1n2(int n1,int n2,double s,double L)
{
  double y;
  
  if(n1==n2)
    y = (2.0/L) * (pow(s/2.0,2)-cos(2 *n1 * M_PI * s)/(8*pow(n1 *M_PI,2))-s*sin(2 *n1 * M_PI * s)/(4*n1 *M_PI));
  else
    y = (2.0/L) * ((cos((n1-n2)*M_PI*s)+(n1-n2)*M_PI*s*sin((n1-n2)*M_PI*s))/pow(n1-n2,2)-(cos((n1+n2)*M_PI*s)+(n1+n2)*M_PI*s*sin((n1+n2)*M_PI*s))/pow(n1+n2,2))/(2*pow(M_PI,2));
  
  return y;
}


//-----------------------------------------------------------------------------------------------
//-- INTEGRAL ALONG R-DIRECTION FOR QDOTS

// QDOT Integrals

// type 1
double QDOT_mass_fint1(double u,void *params)//-- DERIVATES
{
  double y = 0.0;
  struct f_params *p=(struct f_params *) params;
  
  int    l1  = p->l1;
  int    m1  = p->m1;
  int    n1  = p->n1;
  int    l2  = p->l2;
  int    m2  = p->m2;
  int    n2  = p->n2;

  int i,Nm,j;
  FILE *gl;
  char col[100],file[100];
  double L,zz[100],V[100],dm[100],dhz[100],dhp[100],dz[100],tp[100];
  if(strcmp(HOLE,"si")==0)
    sprintf(col,"vb_");  
  else if(strcmp(HOLE,"no")==0)
    sprintf(col,"cb_");
  else;
  
  sprintf(file,"../datafile-%s/%sparameters.dat",Q_Shape,col);
  
  gl=fopen(file,"r");
  i=0;
  do{
    fscanf(gl,"%lf %lf %lf %lf %lf %lf",&V[i],&dm[i],&dhz[i],&dhp[i],&dz[i],&tp[i]);
    i++;
  }while(!feof(gl));
  fclose(gl);
  
  Nm = i-1;
  
  zz[0] = 0;
  for(i=0;i<Nm;i++) zz[i+1]=zz[i]+dz[i];
  
  L = zz[Nm];

  
  y += u * (DR_lm(u,l1,m1) * DR_lm(u,l2,m2)+(l1*l2/pow(u,2))*R_lm(u,l1,m1)*R_lm(u,l2,m2)) * (Sn1n2(n1,n2,(zz[P_pos]+z_QDShape(u))/L,L)-Sn1n2(n1,n2,zz[P_pos]/L,L));
  
  return y;
}

double QDOT_mass_type1(int l1,int m1,int n1,int l2,int m2,int n2,double a,double b)
{
  double y,erabs=1e-8,errel=1e-8,integral,error;
  FILE *gl;

  gsl_function F;
  struct f_params params={l1,m1,n1,l2,m2,n2};
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
  F.function=&QDOT_mass_fint1;
  F.params=&params;
  
  gsl_integration_qag(&F,a,b,erabs,errel,1000,1,w,&integral,&error);
  gsl_integration_workspace_free(w);
  
  y = integral;
  
  return y;
}

// type 2
double QDOT_mass_fint2(double u,void *params)//-- DERIVATES
{
  double y = 0.0;
  struct f_params *p=(struct f_params *) params;
  
  int    l1 = p->l1;
  int    m1 = p->m1;
  int    n1 = p->n1;
  int    l2 = p->l2;
  int    m2 = p->m2;
  int    n2 = p->n2;

  int i,Nm,j;
  FILE *gl;
  char col[100],file[100];
  double L,zz[100],V[100],dm[100],dhz[100],dhp[100],dz[100],tp[100];
  if(strcmp(HOLE,"si")==0)
    sprintf(col,"vb_");  
  else if(strcmp(HOLE,"no")==0)
    sprintf(col,"cb_");
  else;
  
  sprintf(file,"../datafile-%s/%sparameters.dat",Q_Shape,col);
  
  gl=fopen(file,"r");
  i=0;
  do{
    fscanf(gl,"%lf %lf %lf %lf %lf %lf",&V[i],&dm[i],&dhz[i],&dhp[i],&dz[i],&tp[i]);
    i++;
  }while(!feof(gl));
  fclose(gl);
  
  Nm = i-1;
  
  zz[0] = 0;
  for(i=0;i<Nm;i++) zz[i+1]=zz[i]+dz[i];
  
  L = zz[Nm];
  
  y += u * pow(Rs/L,2) * R_lm(u,l1,m1) * R_lm(u,l2,m2) * (Cn1n2(n1,n2,(zz[P_pos]+z_QDShape(u))/L,L)-Cn1n2(n1,n2,zz[P_pos]/L,L));
  
  return y;
}

double QDOT_mass_type2(int l1,int m1,int n1,int l2,int m2,int n2,double a,double b)
{
  double y,erabs=1e-8,errel=1e-8,integral,error;
  FILE *gl;

  gsl_function F;
  struct f_params params={l1,m1,n1,l2,m2,n2};
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
  F.function=&QDOT_mass_fint2;
  F.params=&params;

  gsl_integration_qag(&F,a,b,erabs,errel,1000,1,w,&integral,&error);
  gsl_integration_workspace_free(w);

  y = integral;
  
  return y;
}

// type 3
double QDOT_mass_fint3(double u,void *params)//-- DERIVATES
{
  double y = 0.0;
  struct f_params *p=(struct f_params *) params;
  
  int    l1 = p->l1;
  int    m1 = p->m1;
  int    n1 = p->n1;
  int    l2 = p->l2;
  int    m2 = p->m2;
  int    n2 = p->n2;

  int i,Nm,j;
  FILE *gl;
  char col[100],file[100];
  double L,zz[100],V[100],dm[100],dhz[100],dhp[100],dz[100],tp[100];
  if(strcmp(HOLE,"si")==0)
    sprintf(col,"vb_");  
  else if(strcmp(HOLE,"no")==0)
    sprintf(col,"cb_");
  else;
  
  sprintf(file,"../datafile-%s/%sparameters.dat",Q_Shape,col);
  
  gl=fopen(file,"r");
  i=0;
  do{
    fscanf(gl,"%lf %lf %lf %lf %lf %lf",&V[i],&dm[i],&dhz[i],&dhp[i],&dz[i],&tp[i]);
    i++;
  }while(!feof(gl));
  fclose(gl);
  
  Nm = i-1;
  
  zz[0] = 0;
  for(i=0;i<Nm;i++) zz[i+1]=zz[i]+dz[i];
  
  L = zz[Nm];
  
  y = u * R_lm(u,l1,m1)*R_lm(u,l2,m2)*(Sn1n2(n1,n2,(zz[P_pos]+z_QDShape(u))/L,L)-Sn1n2(n1,n2,zz[P_pos]/L,L));
  
  //printf("%lf\n",L);
  return y;
}

double QDOT_mass_type3(int l1,int m1,int n1,int l2,int m2,int n2,double a,double b)
{
  double y,erabs=1e-8,errel=1e-8,integral,error;
  FILE *gl;

  gsl_function F;
  struct f_params params={l1,m1,n1,l2,m2,n2};
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
  F.function=&QDOT_mass_fint3;
  F.params=&params;

  gsl_integration_qag(&F,a,b,erabs,errel,1000,1,w,&integral,&error);
  gsl_integration_workspace_free(w);

  y = integral;
  
  return y;
}


double QDOT_bias_finteg(double u,void *params)
{
  double y;
  struct f_params *p=(struct f_params *) params;

  int l1 = p->l1;
  int m1 = p->m1;
  int n1 = p->n1;
  int l2 = p->l2;
  int m2 = p->m2;
  int n2 = p->n2;

  int i,Nm,j;
  FILE *gl;
  char col[100],file[100];
  double L,zz[100],V[100],dm[100],dhz[100],dhp[100],dz[100],tp[100];
  if(strcmp(HOLE,"si")==0)
    sprintf(col,"vb_");  
  else if(strcmp(HOLE,"no")==0)
    sprintf(col,"cb_");
  else;
  
  sprintf(file,"../datafile-%s/%sparameters.dat",Q_Shape,col);
  
  gl=fopen(file,"r");
  i=0;
  do{
    fscanf(gl,"%lf %lf %lf %lf %lf %lf",&V[i],&dm[i],&dhz[i],&dhp[i],&dz[i],&tp[i]);
    i++;
  }while(!feof(gl));
  fclose(gl);
  
  Nm = i-1;
  
  zz[0] = 0;
  for(i=0;i<Nm;i++) zz[i+1]=zz[i]+dz[i];
  
  L = zz[Nm];
  
  y = u * R_lm(u,l1,m1)*R_lm(u,l2,m2)*(Bn1n2(n1,n2,(zz[P_pos]+z_QDShape(u))/L,L)-Bn1n2(n1,n2,zz[P_pos]/L,L));

  return y;
}

double QDOT_bias_integral(int l1,int m1,int n1,int l2,int m2,int n2,double a,double b)
{
  double y,erabs=1e-8,errel=1e-8,integral,error;
  FILE *gl;

  gsl_function F;
  struct f_params params={l1,m1,n1,l2,m2,n2};
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
  F.function=&QDOT_bias_finteg;
  F.params=&params;

  gsl_integration_qag(&F,a,b,erabs,errel,1000,1,w,&integral,&error);
  gsl_integration_workspace_free(w);

  y=integral;

  return y;
}


// QWell Integrals

double QWell_a_finteg(double u,void *params)//-- DERIVATES
{
  double y = 0.0;
  struct f_params *p=(struct f_params *) params;

  int l1=p->l1;
  int m1=p->m1;
  int n1=p->n1;
  int l2=p->l2;
  int m2=p->m2;
  int n2=p->n2;
  
  y += u * DR_lm(u,l1,m1) * DR_lm(u,l2,m2)+l1*l2*R_lm(u,l1,m1)*R_lm(u,l2,m2)/u;

  //* (Sn1n2(n1,n2,zz[P_pos]+z_QDShape(u),L)-Sn1n2(n1,n2,zz[P_pos],L));
  //y += u * pow(1.0*Rs,2) * R_lm(u,l1,m1) * R_lm(u,l2,m2) * (Cn1n2(n1,n2,zz[P_pos]+z_QDShape(u),L)-Cn1n2(n1,n2,zz[P_pos],L));
  
  return y;
}

double QWell_a_integral(int l1,int m1,int n1,int l2,int m2,int n2,double a,double b)
{
  double y,erabs=1e-8,errel=1e-8,integral,error;
  FILE *gl;

  gsl_function F;
  struct f_params params={l1,m1,n1,l2,m2,n2};
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
  F.function=&QWell_a_finteg;
  F.params=&params;

  gsl_integration_qag(&F,a,b,erabs,errel,1000,1,w,&integral,&error);
  gsl_integration_workspace_free(w);

  y=integral;

  return y;
}

double QWell_b_finteg(double u,void *params)
{
  double y = 0.0;
  
  struct f_params *p=(struct f_params *) params;

  int l1 = p->l1;
  int m1 = p->m1;
  int n1 = p->n1;
  int l2 = p->l2;
  int m2 = p->m2;
  int n2 = p->n2;
  
  y += u * R_lm(u,l1,m1) * R_lm(u,l2,m2);

  return y;
}

double QWell_b_integral(int l1,int m1,int n1,int l2,int m2,int n2,double a,double b)
{
  double y,erabs=1e-8,errel=1e-8,integral,error;
  FILE *gl;

  gsl_function F;
  struct f_params params={l1,m1,n1,l2,m2,n2};
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
  F.function=&QWell_b_finteg;
  F.params=&params;

  gsl_integration_qag(&F,a,b,erabs,errel,1000,1,w,&integral,&error);
  gsl_integration_workspace_free(w);

  y=integral;

  return y;
}

//-----------------------------------------------------------------------------------------------

//---- MATRIX ELEMENT 

double AELEMENT(int l1,int m1,int n1,int l2,int m2,int n2)
{
  double y,y_all,y_Qdot,y_Wells,y_Bias,L,zz[100],V[100],dm[100],dhz[100],dhp[100],dz[100],tp[100];
  FILE *gl;
  int i,Nm,j,hole;
  
  char col[100],file[100];
  if(strcmp(HOLE,"si")==0)
    sprintf(col,"vb_");  
  else if(strcmp(HOLE,"no")==0)
    sprintf(col,"cb_");
  else;
  
  sprintf(file,"../datafile-%s/%sparameters.dat",Q_Shape,col);
  
  gl=fopen(file,"r");
  i=0;
  do{
    fscanf(gl,"%lf %lf %lf %lf %lf %lf",&V[i],&dm[i],&dhz[i],&dhp[i],&dz[i],&tp[i]);
    i++;
  }while(!feof(gl));
  fclose(gl);
  
  Nm = i-1;
  
  zz[0] = 0;
  for(i=0;i<Nm;i++) zz[i+1]=zz[i]+dz[i];
  
  L = zz[Nm];
  
  //*** Initialization of the cylinder!
  if(strcmp(HOLE,"si")==0)
    y_all = Delta(n1,n2)*Delta(m1,m2)*Delta(l1,l2) * (V[0]/Ry-(pow(kzero(l1,m1)/(Rs/ao),2)/dhp[0]+pow(n1*M_PI/(zz[Nm]/ao),2)/dhz[0]));
  else if(strcmp(HOLE,"no")==0)
    y_all = Delta(n1,n2)*Delta(m1,m2)*Delta(l1,l2) * ((pow(n1*M_PI/(zz[Nm]/ao),2)+pow(kzero(l1,m1)/(Rs/ao),2))/dm[0]+V[0]/Ry);
  else;
  
  y_Wells = 0.0;
  y_Qdot  = 0.0;

  for(j=1;j<Nm-1;j++){
    if(j == P_pos){
      if(strcmp(HOLE,"si") == 0 ){
	y_Qdot += (-1.0) * pow(ao,3) * (L/ao) * Delta(l1,l2) *((1.0/dhp[P_pos]-1.0/dhp[0]) * QDOT_mass_type1(l1,m1,n1,l2,m2,n2,0.0,Ro/Rs)+(1.0/dhz[P_pos]-1.0/dhz[0])* QDOT_mass_type2(l1,m1,n1,l2,m2,n2,0.0,Ro/Rs));
	y_Qdot += pow(ao,3) * (L/ao) * Delta(l1,l2) * pow(Rs/ao,2) * ((V[P_pos]-V[0])/Ry) * QDOT_mass_type3(l1,m1,n1,l2,m2,n2,0.0,Ro/Rs);

	if(fabs(E_field) != 0)
	  y_Qdot += (E_field * L/Ry) * pow(ao,3) * (L/ao)* pow(Rs/ao,2) * Delta(l1,l2) * QDOT_bias_integral(l1,m1,n1,l1,m2,n2,0,Ro/Rs);
	
      }
      else if(strcmp(HOLE,"no") == 0 ){
	y_Qdot += pow(ao,3) * (L/ao) * Delta(l1,l2) *((1.0/dm[P_pos]-1.0/dm[0]) * (QDOT_mass_type1(l1,m1,n1,l2,m2,n2,0.0,Ro/Rs)+ QDOT_mass_type2(l1,m1,n1,l2,m2,n2,0.0,Ro/Rs)));
	y_Qdot += pow(ao,3) * (L/ao) * Delta(l1,l2) * pow(Rs/ao,2) * ((V[P_pos]-V[0])/Ry) * QDOT_mass_type3(l1,m1,n1,l2,m2,n2,0.0,Ro/Rs);

	if(fabs(E_field) != 0)
	  y_Qdot += (E_field * L/Ry) * pow(ao,3) * (L/ao)* pow(Rs/ao,2) * Delta(l1,l2) * QDOT_bias_integral(l1,m1,n1,l1,m2,n2,0,Ro/Rs);

      }
      else;
    }
    else{
      if(strcmp(HOLE,"si") == 0 ){
	y_Wells +=  (-1.0)*pow(ao,3) * (L/ao) * Delta(l1,l2) * ((1.0/dhp[j]-1.0/dhp[0]) * QWell_a_integral(l1,m1,n1,l1,m2,n2,0,1.0) * (Sn1n2(n1,n2,zz[j+1]/L,L)-Sn1n2(n1,n2,zz[j]/L,L)) + pow(Rs/L,2) * (1.0/dhz[j]-1.0/dhz[0]) * (Delta(m1,m2)/pow(Rs,2)) * (Cn1n2(n1,n2,zz[j+1]/L,L)-Cn1n2(n1,n2,zz[j]/L,L)));
	y_Wells += pow(ao,3) * (L/ao) * Delta(l1,l2) * pow(Rs/ao,2)  * ((V[j]-V[0])/Ry) * (Delta(m1,m2)/pow(Rs,2)) * (Sn1n2(n1,n2,zz[j+1]/L,L)-Sn1n2(n1,n2,zz[j]/L,L));

	if(fabs(E_field) != 0)
	  y_Wells += (E_field * L/Ry) * pow(ao,3) * (L/ao)* pow(Rs/ao,2) * Delta(l1,l2) * (Delta(m1,m2)/pow(Rs,2)) * (Bn1n2(n1,n2,zz[j+1]/L,L)-Bn1n2(n1,n2,zz[j]/L,L));

      }
      else if(strcmp(HOLE,"no") == 0 ){
	y_Wells += pow(ao,3) * (L/ao) * Delta(l1,l2) * (1.0/dm[j]-1.0/dm[0]) * ( QWell_a_integral(l1,m1,n1,l1,m2,n2,0,1.0) * (Sn1n2(n1,n2,zz[j+1]/L,L)-Sn1n2(n1,n2,zz[j]/L,L)) + pow(Rs/L,2) * (Delta(m1,m2)/pow(Rs,2)) * (Cn1n2(n1,n2,zz[j+1]/L,L)-Cn1n2(n1,n2,zz[j]/L,L)));
	y_Wells += pow(ao,3) * (L/ao) * Delta(l1,l2) * pow(Rs/ao,2)  * ((V[j]-V[0])/Ry) * (Delta(m1,m2)/pow(Rs,2)) * (Sn1n2(n1,n2,zz[j+1]/L,L)-Sn1n2(n1,n2,zz[j]/L,L));
	
	if(fabs(E_field) != 0)
	  y_Wells += (E_field * L/Ry) * pow(ao,3) * (L/ao)* pow(Rs/ao,2) * Delta(l1,l2) * (Delta(m1,m2)/pow(Rs,2)) * (Bn1n2(n1,n2,zz[j+1]/L,L)-Bn1n2(n1,n2,zz[j]/L,L));

      }
      else;
    }
  }
  y = y_all+y_Qdot+y_Wells;
  
  return y;
}


//.. To plot wave function
//--- Z-WAVEFUNCTION

double Z_n(double z,int n,double L)
{
  double y;

  y=sqrt(2.0/L)*sin(n*M_PI*(z/L));

  return y;
}

//--- R-WAVEFUNCTION

double Rfig_lm(double r,int l,int m)
{
  double y;

  y=(sqrt(2)/(Rs*gsl_sf_bessel_Jn(l+1,kzero(l,m))))*gsl_sf_bessel_Jn(l,kzero(l,m)*(r/Rs));

  return y;
}

//---- WAVE FUNCTIONS 

double modFz_Onda(double z,double L,int l,int k)
{
  double Coef[DIM];
  double _Complex  y;
  int ll[DIM],mm[DIM],nn[DIM],i,j,n;
  FILE *gl;
  char param[1000],param1[1000],file[1000];
  
  sprintf(file,"../datafile-%s/VECTOR_L%d_E%d.dat",Q_Shape,l,k);
  gl=fopen(file,"r");
  
  i=0;
  do{
    fscanf(gl,"%d %d %d %lf",&ll[i],&mm[i],&nn[i],&Coef[i]);
    i++;
  }while(!feof(gl));
  fclose(gl);
  n = i-1;
  
  y=0.0;
  for(i=0;i<n;++i)
    {
    for(j=0;j<n;++j)
      {
	y += conj(Coef[i])*Coef[j]*Z_n(z,nn[i],L)*Z_n(z,nn[j],L);
      }
    }
  
  return cabs(y);
}

double _Complex F_Onda(double r,double phi,double z,int l,int k)
{
  double Coefr[DIM],Coefi[DIM],zz[1000],V[1000],dm[1000],dhh[1000],dz[1000],tp[1000],L;
  double _Complex y;
  int ll[DIM],mm[DIM],nn[DIM],i,j,Nm,kk,m1,n1;
  FILE *gl;
  
  char col[100],file[100];

  if(strcmp(HOLE,"si")==0)
    sprintf(col,"vb_");
  else
    sprintf(col,"cb_");
  sprintf(file,"../datafile-%s/%sparameters.dat",Q_Shape,col);

  gl=fopen(file,"r");
  i=0;
  do{
    fscanf(gl,"%lf %lf %lf %lf %lf",&V[i],&dm[i],&dhh[i],&dz[i],&tp[i]);
    i++;
  }while(!feof(gl));
  fclose(gl);
  
  Nm = i-1;

  zz[0] = 0;
  for(i=0;i<Nm;i++) zz[i+1]=zz[i]+dz[i];

  L = zz[Nm];

  sprintf(file,"../datafile-%s/%sVECTOR_L%d_E%03d.dat",Q_Shape,col,l,k);
  gl=fopen(file,"r");
  
  i=0;
  do{
    fscanf(gl,"%d %d %d %lf %lf",&ll[i],&mm[i],&nn[i],&Coefr[i],&Coefi[i]);
    //    printf("%d %d %d %lf %lf\n",ll[i],mm[i],nn[i],Coefr[i],Coefi[i]);
    i++;
  }while(!feof(gl));

  fclose(gl);
  
  kk=0;
  for(m1=1;m1<=Mmax;m1++)
    {
      for(n1=1;n1<=Nmax;n1++)
	{
	  ll[kk]=l; mm[kk]=m1; nn[kk]=n1; kk++;
	}	  
    }
  
  
  y=0.0;
  for(j=0;j<i-1;j++)
    {
      y = y + (Coefr[j]+I*Coefi[j])*Rfig_lm(r,ll[j],mm[j])*Z_n(z,nn[j],L)*(1/sqrt(2*M_PI))*cexp(I*ll[j]*phi);
    }

  return y;
}


//*** Mapppings

void map(double *zz,int Nm)
{
  int i,j,k,l,Nz = 200,Nr = 100;
  FILE *lg;
  double zp,r,th,z1 = zz[0],z2=zz[Nm];
  double dr = (2.0*Rs/3.0)/Nr,dz = (z2-z1)/Nz;
  double cbevs[Lmax][Nlevels],vbevs[Lmax][Nlevels];
  double _Complex Fz;

  char file[100];

  if(strcmp(what2do,"energ")==0){
  for(l=0;l<Lmax;l++){
    sprintf(file,"../datafile-%s/%sENERG_L%d.dat",Q_Shape,"cb_",l);
    lg=fopen(file,"r");
    for(i = 0;i<Nlevels;i++){
      fscanf(lg,"%lf",&cbevs[l][i]);
    }
    fclose(lg);
    
    sprintf(file,"../datafile-%s/%sENERG_L%d.dat",Q_Shape,"vb_",l);
    lg=fopen(file,"r");
    for(i = 0;i<Nlevels;i++){
      fscanf(lg,"%lf",&vbevs[l][i]);
    }
    fclose(lg);
  }  
  
  lg = fopen("mappings.gpl","w");
  fprintf(lg,"set term postscript eps enhanced color\n");
  fprintf(lg,"set output 'fig.eps'\n set nokey\n\n");
  fprintf(lg,"set multiplot\n\n");
  fprintf(lg,"set xlabel 'Growth direction [nm]' font 'Helvetica,20'\n");
  fprintf(lg,"set ylabel 'Energy [eV]' font 'Helvetica,18' offset 0,0\n");
  fprintf(lg,"set style line 101 lc rgb '#808080' lt 1 lw 1\n");
  fprintf(lg,"set border 3 front ls 101\n");
  fprintf(lg,"set title 'Conduction Band : Electrons' font 'Helvetica,18'\n");
  fprintf(lg,"set tics nomirror out scale 0.75\n\n");

  //  fprintf(lg,"set xrange [%.2lf:%.2lf]\n",z1/10,z2/10.0);
  //  fprintf(lg,"set yrange [%.2lf:%.2lf]\n",0.1,0.9);
  fprintf(lg,"set xtics font 'Helvetica,16'\n");
  fprintf(lg,"set ytics font 'Helvetica,16'\n");
  fprintf(lg,"set label 'L = 0' font 'Helvetica,20'  at 92,0.15 tc rgb '#0000ff'\n");  
  fprintf(lg,"set label 'L = 1' font 'Helvetica,20'  at 109,0.15 tc rgb '#ff0000'\n");  
  fprintf(lg,"set label 'L = 2' font 'Helvetica,20'  at 125,0.15 tc rgb '#000000'\n");  
  fprintf(lg,"set size 0.5,0.95\n");
  fprintf(lg,"set origin 0.015,0.015\n");  
  fprintf(lg,"plot\t");
  for(l=0;l<Lmax;l++){
    for(i=0;i<Nlevels;i++){
      if(l==0)
	fprintf(lg,"%lf w l lw 3 dt 1 lc rgb '#0000ff',",cbevs[l][i]);
      else if(l==1)
	fprintf(lg,"%lf w l lw 3 dt 4 lc rgb '#ff0000',",cbevs[l][i]);
      else if(l==2)
	fprintf(lg,"%lf w l lw 4 dt 2 lc rgb '#000000',",cbevs[l][i]);
      else if(l==3)
	fprintf(lg,"%lf w l lw 4 dt 2 lc rgb '#00ff00',",cbevs[l][i]);
      else;
    }
  }
  fprintf(lg,"'../datafile-%s/%sPotential.dat' notitle w l lw 8 dt 1 lc rgb '#00ff00'\n\n",Q_Shape,"cb_");

  fprintf(lg,"unset yrange \n unset label\n");
  fprintf(lg,"set title 'Valence Band : Holes' font 'Helvetica,18'\n");
  fprintf(lg,"set ylabel ''\n");
  //  fprintf(lg,"set yrange [%.2lf:%.2lf]\n",-0.9,-0.25);
  fprintf(lg,"set size 0.5,0.95\n");
  fprintf(lg,"set origin 0.5,0.015\n");  
  fprintf(lg,"plot\t");
  for(l=0;l<Lmax;l++){
    for(i=0;i<Nlevels;i++){
      if(l==0)
	fprintf(lg,"%lf w l lw 3 dt 1 lc rgb '#0000ff',",vbevs[l][i]);
      else if(l==1)
	fprintf(lg,"%lf w l lw 3 dt 4 lc rgb '#ff0000',",vbevs[l][i]);
      else if(l==2)
	fprintf(lg,"%lf w l lw 4 dt 2 lc rgb '#000000',",vbevs[l][i]);
      else if(l==3)
	fprintf(lg,"%lf w l lw 4 dt 2 lc rgb '#00ff00',",vbevs[l][i]);
      else;
    }
  }
  fprintf(lg,"'../datafile-%s/%sPotential.dat' notitle w l lw 8 lc rgb '#00ff00'",Q_Shape,"vb_");

  fclose(lg);
  system("gnuplot mappings.gpl");
  }
  else if(strcmp(what2do,"maps")==0){
    lg = fopen("mappings.gpl","w");
    fprintf(lg,"set term postscript eps enhanced color\n");
    fprintf(lg,"set output 'maps.eps'\n\nset nokey\n");
    fprintf(lg,"set lmargin 0\n");
    fprintf(lg,"set rmargin 0\n");
    fprintf(lg,"set bmargin 0\n");
    fprintf(lg,"set tmargin 0\n");
    fprintf(lg,"set multiplot layout 2,3 margins 0.08,0.95,0.08,0.95 spacing 0.08,0.1\n\n");
    fprintf(lg,"set view map scale 1\nset border 0 lc rgb '#808080' lt 1 lw 1\n");
    fprintf(lg,"set nocolorbox\nset samples 500\nset isosamples 500\n");
    fprintf(lg,"set pm3d interpolate 0,0\nset xyplane relative 0\n");
    fprintf(lg,"set palette defined (0 0 0 0, 1 0 0 1, 3 0 1 0, 4 1 0 0, 6 1 1 1)\n");
    fprintf(lg,"unset tics\nset xtics 20 font 'Helvetica,12'\n");
    fprintf(lg,"set ytics font 'Helvetica,12'\n");
    fprintf(lg,"set yrange [%.2lf:%.2lf]\n",z1/10.0,z2/10.0);
    fprintf(lg,"set xrange [%.2lf:%.2lf]\n",-45.0,45.0);
    
    fprintf(lg,"splot '../datafile-%s/%sWAVEFUNCPlaneRZ_L%d_K%d.dat' u 1:2:3 notitle w pm3d,",Q_Shape,"vb_",0,0);
    fprintf(lg,"'../datafile-%s/Qd_Shape.dat' u 1:($2+%0.1lf):(0.0) w l lw 3 dt 3 lc rgb '#ffff00',",Q_Shape,zz[P_pos]/10.0);
    fprintf(lg,"'../datafile-%s/Qd_Shape.dat' u 1:(0.0*$2+%0.1lf):(0.0) w l lw 3 dt 3 lc rgb '#ffff00'\n",Q_Shape,zz[P_pos-1]/10.0);

    fprintf(lg,"splot '../datafile-%s/%sWAVEFUNCPlaneRZ_L%d_K%d.dat' u 1:2:3 notitle w pm3d,",Q_Shape,"vb_",1,0);
    fprintf(lg,"'../datafile-%s/Qd_Shape.dat' u 1:($2+%0.1lf):(0.0) w l lw 3 dt 3 lc rgb '#ffff00',",Q_Shape,zz[P_pos]/10.0);
    fprintf(lg,"'../datafile-%s/Qd_Shape.dat' u 1:(0.0*$2+%0.1lf):(0.0) w l lw 3 dt 3 lc rgb '#ffff00'\n",Q_Shape,zz[P_pos-1]/10.0);

    fprintf(lg,"splot '../datafile-%s/%sWAVEFUNCPlaneRZ_L%d_K%d.dat' u 1:2:3 notitle w pm3d,",Q_Shape,"vb_",2,0);
    fprintf(lg,"'../datafile-%s/Qd_Shape.dat' u 1:($2+%0.1lf):(0.0) w l lw 3 dt 3 lc rgb '#ffff00',",Q_Shape,zz[P_pos]/10.0);
    fprintf(lg,"'../datafile-%s/Qd_Shape.dat' u 1:(0.0*$2+%0.1lf):(0.0) w l lw 3 dt 3 lc rgb '#ffff00'\n",Q_Shape,zz[P_pos-1]/10.0);

    fprintf(lg,"splot '../datafile-%s/%sWAVEFUNCPlaneRZ_L%d_K%d.dat' u 1:2:3 notitle w pm3d,",Q_Shape,"vb_",0,1);
    fprintf(lg,"'../datafile-%s/Qd_Shape.dat' u 1:($2+%0.1lf):(0.0) w l lw 3 dt 3 lc rgb '#ffff00',",Q_Shape,zz[P_pos]/10.0);
    fprintf(lg,"'../datafile-%s/Qd_Shape.dat' u 1:(0.0*$2+%0.1lf):(0.0) w l lw 3 dt 3 lc rgb '#ffff00'\n",Q_Shape,zz[P_pos-1]/10.0);

    fprintf(lg,"splot '../datafile-%s/%sWAVEFUNCPlaneRZ_L%d_K%d.dat' u 1:2:3 notitle w pm3d,",Q_Shape,"vb_",1,1);
    fprintf(lg,"'../datafile-%s/Qd_Shape.dat' u 1:($2+%0.1lf):(0.0) w l lw 3 dt 3 lc rgb '#ffff00',",Q_Shape,zz[P_pos]/10.0);
    fprintf(lg,"'../datafile-%s/Qd_Shape.dat' u 1:(0.0*$2+%0.1lf):(0.0) w l lw 3 dt 3 lc rgb '#ffff00'\n",Q_Shape,zz[P_pos-1]/10.0);

    fprintf(lg,"splot '../datafile-%s/%sWAVEFUNCPlaneRZ_L%d_K%d.dat' u 1:2:3 notitle w pm3d,",Q_Shape,"vb_",2,1);
    fprintf(lg,"'../datafile-%s/Qd_Shape.dat' u 1:($2+%0.1lf):(0.0) w l lw 3 dt 3 lc rgb '#ffff00',",Q_Shape,zz[P_pos]/10.0);
    fprintf(lg,"'../datafile-%s/Qd_Shape.dat' u 1:(0.0*$2+%0.1lf):(0.0) w l lw 3 dt 3 lc rgb '#ffff00'\n",Q_Shape,zz[P_pos-1]/10.0);
    fprintf(lg,"unset multiplot\n");
    fclose(lg);
    system("gnuplot mappings.gpl");
    //    sprintf(file,"mv maps.eps %s\n",plot_file);
    //system(file);

  }
  else;
  
}


// ******** Diagonalization with lapack: For any matrix!
void Lapackdiag(double _Complex *A, int dim,double _Complex *evals,double _Complex *evecs){
  
  int i,j,N=dim,LDA=N,LDVL=N,LDVR=N;
  int n=N,lda=LDA,ldvl=LDVL,ldvr=LDVR,info;

  double _Complex wkopt;
  double _Complex  *work;
  int lwork = 48*n; // This should be tuned
  double rwork[2*N];
  double _Complex *vl;

  vl=malloc(LDVL*N*sizeof(double _Complex));
  work=malloc(lwork*sizeof(double _Complex));

  for(i=0;i<N;i++)evals[i]=0.0;
  for(i=0;i<N*N;i++){vl[i]=0.0;evecs[i]=0.0;}
  
  // TAKE IMMENSE CARE ===========>> matrix is overwritten after the process has taken place <<===========
  zgeev_("N", "V", &n, A, &n, evals, vl, &ldvl, evecs, &ldvr, work, &lwork, rwork, &info);

  if(((int)(work[0]) > lwork)||(lwork > 3*((int)(work[0])))){ // Warns if workspace is small or huge
    //printf("Eigenvalues workspace: actual dimension %d, optimal dimension %d\n", lwork, (int)(work[0]));
  }

  if(info != 0){
    printf("Eigenvalues: ERROR\n");
  }

  free(work);
  free(vl);

}


void arraysorted(double *a,int length,int *indx){
  int i,j,k;
  double x,atemp[length];//,outa[length];
  double temp;
  
  for(i=0; i<length; i++)atemp[i]=a[i];
  
  for(i=1; i<length; i++){
    for(j=0; j<i; j++){
      if(a[j] >= a[i]){
	temp = a[i];
	a[i] = a[j];
	a[j] = temp;
      }
    }
  }
  k=0;
  for(i=0;i<length;i++){
    for(j=0;j<length;j++){
      if(fabs(a[i]-atemp[j])<1e-22){
	indx[k]=j;
	k++;
	//	printf("%d\n",j);
      }
    }
  }
  //  for(i=0;i<length;i++)printf("%d\t%lf\t%lf\t%d\n",i,atemp[i],a[i],indx[i]);printf("\n");

  //exit(0);
  //for(i=0;i<length;i++)a[i]=atemp[i ];
}


/*
 for(j=1;j<Nm-1;j++){
    
    // Inserting the Qdot!
    if(j == P_pos){
      y_Qdot += hole * pow(ao,3) * (L/ao) * Delta(l1,l2) *(Delta(hole,1)*(1.0/dm[P_pos]-1.0/dm[0])+Delta(hole,-1)*(1.0/dhh[P_pos]-1.0/dhh[0])) * QDOT_mass_integral(l1,m1,n1,l1,m2,n2,0,Ro/Rs);
      y_Qdot += pow(ao,3) * (L/ao) * Delta(l1,l2) * pow(Rs/ao,2)  * ((V[P_pos]-V[0])/Ry) * QDOT_pot_integral(l1,m1,n1,l1,m2,n2,0,Ro/Rs);

      if(fabs(E_field) != 0)
	y_Qdot += (E_field * L/Ry) * pow(ao,3) * (L/ao)* pow(Rs/ao,2) * Delta(l1,l2) * QDOT_bias_integral(l1,m1,n1,l1,m2,n2,0,Ro/Rs);
    }
    else{
      y_Wells += hole * (pow(ao,3) * (L/ao) * Delta(l1,l2) * (Delta(hole,1)*(1.0/dm[j]-1.0/dm[0])+Delta(hole,-1)*(1.0/dhh[j]-1.0/dhh[0])) * (QWell_a_integral(l1,m1,n1,l1,m2,n2,0,1.0) * (Sn1n2(n1,n2,zz[j+1]/L,L)-Sn1n2(n1,n2,zz[j]/L,L)) + pow(Rs/L,2) * (Delta(m1,m2)/pow(Rs,2)) * (Cn1n2(n1,n2,zz[j+1]/L,L)-Cn1n2(n1,n2,zz[j]/L,L))));
      
      y_Wells += pow(ao,3) * (L/ao) * Delta(l1,l2) * pow(Rs/ao,2)  * ((V[j]-V[0])/Ry) * (Delta(m1,m2)/pow(Rs,2)) * (Sn1n2(n1,n2,zz[j+1]/L,L)-Sn1n2(n1,n2,zz[j]/L,L));
      
      if(fabs(E_field) != 0)
	y_Wells += (E_field * L/Ry) * pow(ao,3) * (L/ao)* pow(Rs/ao,2) * Delta(l1,l2) * (Delta(m1,m2)/pow(Rs,2)) * (Bn1n2(n1,n2,zz[j+1]/L,L)-Bn1n2(n1,n2,zz[j]/L,L));
      
    }
  }
 
 */

