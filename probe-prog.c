#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <malloc.h>
#include <omp.h>

/* OWN LIBRARIES */

#include "GLOBAL_CONSTANTS.h"
#include "libraries.h"

int m[DIM],n[DIM],l[DIM];
double _Complex integ[DIM*DIM];
double ms[DIM],ns[DIM],ls[DIM];
double Cr[DIM], Ci[DIM];
double eigvs[DIM];

//*************************************************************************************

int main(int argc, char *argv[])
{

  int k1,k2,kk,n1,n2,m1,m2,i,Nm,l1,l2,j,k,ll,lll,kkk,nkk,mk,p;
  double z,phi,zz[100],V[100],dm[100],dhz[100],dhp[100];
  double dz[100],tp[100],ENElk[1000][DIM],func,Osc_Sth,zp,r,dif_t,Fz,Fr1,Fr2,modFr,u;

  double _Complex *mij,*evec,*eval;
  FILE *lg,*lh;
  time_t start,end,timer;
  struct tm* tm_info;

  char col[100],file[100],nfile[100],buffer[26];

  printf("\n... Code started at:\t");
  time(&timer);
  tm_info = localtime(&timer);
  strftime(buffer, 26, "%Y-%m-%d %H:%M:%S", tm_info);
  puts(buffer);

  if(strcmp(HOLE,"si")==0)
    sprintf(col,"vb_");
  else
    sprintf(col,"cb_");

  if(strcmp(what2do,"diag")==0){
    printf("\n>> Computing %s Energies, please wait\n System dimension : %d states\n\n",col,DIM);
    printf("Shape of the QDoT : %s\n",Q_Shape);
    printf("w. height = %0.1lf nm\n", ho/10.0);
    printf("w. Base-Radius = %0.1lf nm\n", Ro/10.0);
  }
  else if(strcmp(what2do,"maps")==0)
    printf("\n>> Plotting Maps ... please wait\n\n");
  else if(strcmp(what2do,"energ")==0)
    printf("\n>> Plotting Energy Diagram ... please wait\n\n");
  else;

  /***************************
     COMPUTING THE POTENTIAL 
  ****************************/
  Struc_params();

  //*******************************************************
  sprintf(file,"../datafile-%s/%sparameters.dat",Q_Shape,col);
  lg=fopen(file,"r");
  i=0;
  do{
    fscanf(lg,"%lf %lf %lf %lf %lf %lf",&V[i],&dm[i],&dhz[i],&dhp[i],&dz[i],&tp[i]);
    i++;
  }while(!feof(lg));
  fclose(lg);
  Nm = i-1;
  zz[0]=0;  for(i=0;i<Nm;i++)zz[i+1]=(zz[i]+dz[i]);
  //*******************************************************

  
  if(strcmp(what2do,"diag") == 0){
    //  printf("Energy interval [ %lf , %lf]\n", V[P_pos] + E_field * zz[1],V[P_pos+1]+ E_field* zz[1]);
    for(l1=0;l1<Lmax;l1++)
      {
	printf("\nComputing ... L = %d\n\n",l1);
	time(&start);
	
	k=0;
	for(m1=1;m1<=Mmax;m1++)
	  {
	    for(n1=1;n1<=Nmax;n1++)
	      {
		l[k]=l1; m[k]=m1; n[k]=n1; k++;
	      }	  
	  }
	
	mij  = malloc(DIM*DIM*sizeof(double _Complex));
	evec = malloc(DIM*DIM*sizeof(double _Complex));
	eval = malloc( DIM *sizeof(double _Complex));
	for(i=0;i<DIM*DIM;i++){mij[i] = 0.0; evec[i] = 0.0;}
	for(i=0;i<DIM;i++)eval[i] = 0.0;
	
	// to parallelize...
	omp_set_dynamic(0);
#pragma omp parallel for private(j) num_threads(nthread)
	for(i=0;i<DIM;i++){
	  for(j=i;j<DIM;j++){
	    mij[i * DIM + j] = AELEMENT(l[i],m[i],n[i],l[j],m[j],n[j]);
	    mij[j * DIM + i] = conj(mij[i * DIM + j]);
	  }
	}
	time(&end);
	dif_t=difftime(end,start);
	printf("Matrix construction time ... in  %0.4lf segs...!\n",dif_t);
	
	time(&start);
	Lapackdiag(mij,DIM,eval,evec);
	time(&end);
	dif_t=difftime(end,start);
	printf("Diagonalization time ... in  %0.4lf segs...!\n",dif_t);
	
	int indx[DIM];
	double eig[DIM],evs;
	
	for(i=0;i<DIM;i++) eig[i] = creal(eval[i]);
	arraysorted(eig,DIM,indx);
	
	sprintf(file,"../datafile-%s/%sENERG_L%d.dat",Q_Shape,col,l1);
	lh=fopen(file,"w");
	
	//      printf("%lf %lf\n",V[P_pos],V[P_pos+1]);
	ll = 0;
	for(i=0;i<Nlevels;i++)
	  {
	    if(strcmp(HOLE,"si")==0)
	      k = DIM-1-i;
	    else
	      k = i;
	    
	    evs = creal(eval[indx[k]])*Ry;	  
	    
	    if( evs > 0.35 && evs < .9){
	      fprintf(lh,"%11.11lf\n",evs);
	      
	      sprintf(file,"../datafile-%s/%sVECTOR_L%d_E%03d.dat",Q_Shape,col,l1,ll);
	      lg=fopen(file,"w");
	      for(j=0;j<DIM;j++)
		{
		  double _Complex z = evec[ indx[k] * DIM + j];
		  fprintf(lg,"%d\t%d\t%d\t%11.11lf\t%11.11lf\n",l[j],m[j],n[j],creal(z),cimag(z));
		} 
	      fclose(lg);
	      ll = ll + 1;	      
	    }
	    
	  }
	printf("# states =  %d\n",ll);	  
	fclose(lh);
	printf("\n");
	
	free(evec); free(eval); free(mij);      
      }

    double evssel[100];
    
    ll = 0;
    printf("\nHellow ... Eigsel ON:\n");
    for(l1 = 0;l1<Lmax;l1++){
      
      k=0;
      for(m1=1;m1<=Mmax;m1++)
	{
	  for(n1=1;n1<=Nmax;n1++)
	    {
	      l[k]=l1; m[k]=m1; n[k]=n1; k++;
	    }	  
	}
      
      sprintf(file,"../datafile-%s/%sENERG_L%d.dat",Q_Shape,col,l1);
      lg=fopen(file,"r");
      i=0;
      do{
	fscanf(lg,"%lf",&eigvs[i]);
	i++;
      }while(!feof(lg));
      fclose(lg);
      
      int qts = i-1;
      printf("# of states to examinate = %d\n",qts);
      

      for(mk = 0;mk<qts;mk++){
	sprintf(file,"../datafile-%s/%sVECTOR_L%d_E%03d.dat",Q_Shape,col,l1,mk);
	lg=fopen(file,"r");
	i=0;
	do{
	  fscanf(lg,"%lf %lf %lf %lf %lf",&ls[i],&ms[i],&ns[i],&Cr[i],&Ci[i]);
	  i++;
	}while(!feof(lg));
	fclose(lg);
	
	double _Complex tot;
	double ntot;
	for(p=0;p<DIM*DIM;p++) integ[p] = 0.0;
	{// to parallelize...
	  const int ntr = 20;
	  omp_set_dynamic(0);
#pragma omp parallel for private(j) num_threads(nthread)
	  for(p=0;p<DIM;p++){
	    for(j=p;j<DIM;j++){
	      if(fabs(Cr[p])>1.0e-4 && fabs(Cr[j])>1.0e-4){
		integ[p*DIM+j] = conj(Cr[p]) * Cr[j] * QDOT_mass_type3(l[p],m[p],n[p],l[j],m[j],n[j],0.0,Ro/Rs); 
		integ[j*DIM+p] = integ[p*DIM+j];
	      }
	    }
	  }
	}
	tot = 0.0; for(p=0;p<DIM*DIM;p++) tot += integ[p];
	ntot = cabs(tot) * pow(ao,3) * pow(Rs/ao,2) * (zz[Nm]/ao);
	
	printf("E(%03d) = %0.11lf\tFilling => factor = %0.6lf\n",mk,eigvs[mk],ntot);
	if(ntot < 0.9){
	  sprintf(nfile,"rm %s",file);
	  system(nfile);
	}
	else{
	  //  fprintf(lh,"E(%03d) = %11.11lf\tFilling factor = %0.6lf\n",mk,eigvs[mk],ntot);
	  evssel[ll] = eigvs[mk];
	  //	  if(ll>3)   break;
	  ll++;
	}
      }
      //      fclose(lh);
    }
    
    if(ll>0){
      int idx[ll];
      arraysorted(evssel,ll,idx);
      
      sprintf(file,"../datafile-%s/new_%s_Ro%0.1lf_ENERG.dat",Q_Shape,col,Ro);
      lh=fopen(file,"w");
    
      for(i=0;i<ll;i++)
	fprintf(lh,"%11.11lf\n",mk,evssel[idx[ll-i-1]]);
      
      fclose(lh);
    }
    else
      printf("Zero confined levels!\n\nSchüss!\\");
      
    

  }  
  else if(strcmp(what2do,"maps") == 0){
    // Plotting mappings
    printf("mierda...\n");
    
    int i,j,k,l,Nz = 400,Nr = 100;
    FILE *lg;
    char file[100];
    double _Complex Fz;
    double zp,r,th,z1 = zz[0],z2=zz[Nm];
    double dr = (Rs/3.0)/Nr,dz = (z2-z1)/Nz;

    if(strcmp(HOLE,"si")==0)
      sprintf(col,"vb_");
    else
      sprintf(col,"cb_");
    
    l = 0;
    k = 0;
        
    sprintf(file,"../datafile-%s/%sWAVEFUNCPlaneRZ_L%d_K%d.dat",Q_Shape,col,l,k);
    lg=fopen(file,"w");
    
    for(r = 0; r< Rs/3.0; r += dr){
      for(zp = z1; zp <= z2; zp += dz){
	fprintf(lg,"%11.8lf\t%11.8lf\t%11.22lf\n",0.1*r,0.1*zp,pow(cabs(F_Onda(r,0.0,zp,l,k)),2));
      }
      fprintf(lg,"\n");
    }
    fclose(lg);

    
  }
  else;

  printf("\n...Code ended at:\t");
  time(&timer);
  tm_info = localtime(&timer);
  strftime(buffer, 26, "%Y-%m-%d %H:%M:%S", tm_info);
  puts(buffer);

  
  return 0;
}
