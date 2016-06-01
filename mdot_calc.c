#ifdef RAYTRACE_TG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "allvars.h"

#define pi 3.1415927
#define c 2.99792458e10
#define h_nu 6.6262e-27
#define k_B 1.3806e-16
#define pc 3.085678e18

double mdot_calc(int numtot)
  {
  FILE *sinkmasses;
  FILE *starmodel;
  FILE *starpres;
  int n, numstart, nID, nsteps_int = 20;
  double nsteps, mdot, t_cur, t, lum_acc;
  double sinkposx, sinkposy, sinkposz, dis, xpos, ypos, zpos, total, a, b;
  double ztime, time, time_avg, redshift, tacc, tacc_cur, sinkmass, sinkmass_read;
  double time0, min_fac=2.0;
  double lnx, lny, lnx_sqrd, lnxlny, Acoeff, Bcoeff, Ccoeff, Dcoeff;
  int ID, ID0, NumCurrentTimestep;
  char fsinkmasses[200];
  char fstarmodel[200], fstarpres[200];

  sprintf(fsinkmasses,"%s/sinkmasses_hiacc", All.OutputDir);
  sprintf(fstarmodel,"%s/starmodel_hiacc", All.OutputDir); 
  sprintf(fstarpres,"%s/starpres_hiacc", All.OutputDir);
  printf("mdot line 25\n");

  sinkmasses = fopen(fsinkmasses, "r"); 

  nID = 0;
  lnx = 0;
  lny = 0;
  lnxlny = 0;
  lnx_sqrd = 0;
  time_avg = 10.0;

  for(n = 0; n < numtot; n++)
   {
   fscanf(sinkmasses, "%lg %d %lg %d %lg %lg %lg", &ztime, &ID, &sinkmass_read, &NumCurrentTimestep, &sinkposx, &sinkposy, &sinkposz);
   redshift = (1.0/ztime) - 1.0;
   time=5.4e8/pow((1.e0+redshift)/10.e0,1.5e0);
   if(n==0)
     {
     time0 = time;
     ID0 = ID;
     }
   if(ztime <= All.Time)
     tacc_cur = time - time0;
 //  if((numtot - n) <= nsteps_int  && ID1 == 3755078)
 //    {
//     sinkmass =  sinkmass_read*1.e10/0.7;
//     lnx = lnx + log(tacc);
//     lny = lny + log(sinkmass);
//     lnxlny = lnxlny + log(tacc)*log(sinkmass);
 //    lnx_sqrd = lnx_sqrd + log(tacc)*log(tacc);
 //    nID++;
 //    }
   }

  fclose(sinkmasses);

  sinkmasses = fopen(fsinkmasses, "r");

//  for()
//{
  for(n = 0; n < numtot; n++)
    {
    fscanf(sinkmasses, "%lg %d %lg %d %lg %lg %lg", &ztime, &ID, &sinkmass_read, &NumCurrentTimestep, &sinkposx, &sinkposy, &sinkposz);
    redshift = (1.0/ztime) - 1.0;
    time=5.4e8/pow((1.e0+redshift)/10.e0,1.5e0);
    tacc = time - time0;
    //if(tacc_cur - tacc < 100.0  &&  tacc > 0 && ID == All.ray_center_ID)
    if(tacc_cur - tacc < time_avg  &&  tacc > 0 && ID == ID0 && ztime < All.Time)
      {
      sinkmass =  sinkmass_read*1.e10/0.7;
      lnx = lnx + log(tacc);
      lny = lny + log(sinkmass);
      lnxlny = lnxlny + log(tacc)*log(sinkmass);
      lnx_sqrd = lnx_sqrd + log(tacc)*log(tacc);
      nID++;
      }
     }

  fclose(sinkmasses);

  nsteps = (double)nID; 

  printf("mdot line 58\n");

  Bcoeff = (nsteps*(lnxlny) - lnx*lny)/(nsteps*lnx_sqrd - pow(lnx,2));
  Acoeff = (lny - Bcoeff*lnx)/nsteps;
  Acoeff = pow(2.71818,Acoeff);
 
  Ccoeff = Bcoeff*pow(Acoeff,(1.0/Bcoeff));
  Dcoeff = (Bcoeff-1.0)/Bcoeff;

  mdot= Ccoeff*pow(sinkmass, Dcoeff);

//  if(mdot < mdot_min)
//    time_avg = time_avg + 50.0*j;
//}

  if(mdot != mdot)
    mdot = 1.e-7;

  All.tacc = tacc_cur*3.15569e7;
  lum_acc = All.alpha*GRAVITY * (All.star_mass*1.98892e33)* (mdot*1.98892e33/3.1557e7) / All.star_rad;

  starmodel = fopen(fstarmodel, "a");
  fprintf(starmodel, "%17.13g %8d %15.6g %10d %15.11g  %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g\n", All.Time, ID0, All.tacc, All.NumCurrentTiStep, sinkmass, mdot, All.star_rad, lum_acc, All.lum_tot, All.Teff, ray.Q_H_ion, ray.Q_He_ion, ray.Qion_zams, All.alpha, All.Q_LW, All.Prad_avg);
  fclose(starmodel);

  starpres = fopen(fstarpres, "a");
  fprintf(starpres, "%17.13g  %10d %15.11g %15.11g  %15.11g %15.11g  %15.11g %15.11g\n", All.Time, All.NumCurrentTiStep, All.Fdir_avg, All.arad_avg, All.adir_avg, All.agrav_avg, All.Pres_avg, All.Prad_avg);
  fclose(starpres);

  printf("mdot_calc is working! A= %lg, B = %lg, C = %lg, D = %lg\n", Acoeff, Bcoeff, Ccoeff, Dcoeff);
  printf("mdot_calc time = %lg, mass = %lg mdot = %lg \n", tacc, sinkmass, mdot);

  return(mdot); 
  }

#endif

