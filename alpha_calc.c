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

double alpha_calc(int n_check)
  {
  FILE *sinkangmom;
  FILE *sinkangmom_near;
  int n, numstart, nID, nsteps_int = 20;
  double nsteps, alpha_val, t_cur, t, jfrac, jcent, jfrac_tot=0;
  double sinkposx, sinkposy, sinkposz, dis, xpos, ypos, zpos, total, a, b;
  double ztime, time, time_avg, redshift, tacc, tacc_cur, sinkmass, sinkmass_read;
  double time0;
  double d_sink_phys, vx_phys, vy_phys, vz_phys, jx, jy, jz, jtot;
  double n_jlow_doub, n_jhigh_doub, n_tot_doub;
  int ID, ID0, NumCurrentTimestepi, n_jlow=0, n_jhigh=0, n_tot=0, f_tot=0, f_check;
  char fsinkangmom[200], fsinkangmom_near[200];

  sprintf(fsinkangmom,"%s/sinkangmom_hiacc", All.OutputDir);
  sprintf(fsinkangmom_near,"%s/sinkangmom_near_hiacc", All.OutputDir);

  time_avg = 10.0;
  f_check = (int)(All.sinkmass_sum_tot*1.e10/.7/6.33e-4);
  f_tot = n_check;
  if(f_tot < f_check)
    f_tot = f_check;

  if(f_tot < 10)
    f_tot = 10;

  //printf("n_near = %d\n", All.n_near_max);

  sinkangmom = fopen(fsinkangmom, "r"); 

  for(n = 0; n < f_tot; n++)
   {
   fscanf(sinkangmom, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %d\n", &ztime, &d_sink_phys, &vx_phys, &vy_phys, &vz_phys, &jx, &jy, &jz, &jtot, &sinkmass_read, &ID);
   redshift = (1.0/ztime) - 1.0;
   time=5.4e8/pow((1.e0+redshift)/10.e0,1.5e0);
   if(n==0)
     {
     time0 = time;
     ID0 = ID;
     }
   if(ztime <= All.Time)
    tacc_cur = time - time0;
    }   
  fclose(sinkangmom);


  sinkangmom = fopen(fsinkangmom, "r");
  for(n = 0; n < f_tot; n++)
   {
   fscanf(sinkangmom, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %d\n", &ztime, &d_sink_phys, &vx_phys, &vy_phys, &vz_phys, &jx, &jy, &jz, &jtot, &sinkmass_read, &ID);
   redshift = (1.0/ztime) - 1.0;
   time=5.4e8/pow((1.e0+redshift)/10.e0,1.5e0);
   tacc = time - time0;

   jcent = pow((6.67e-8*sinkmass_read*(1.e10/0.7)*1.98892e33*d_sink_phys),0.5);
   jfrac = jtot/jcent;

   if(ID == ID0 && tacc_cur - tacc < time_avg  &&  tacc > 0 && tacc_cur - tacc > 0)
     jfrac_tot = jfrac_tot+jfrac;

   if(n == f_tot - 50)
     printf("f_tot = %d, n_jlow = %10d n_jhigh = %10d jtot = %10.5g jcent = %10.5g jfrac = %lg ,jfrac_tot = %lg, tacc = %lg\n", f_tot, n_jlow, n_jhigh, jtot, jcent, jfrac, jfrac_tot, tacc);

   if(ID == ID0 && tacc_cur - tacc < time_avg  &&  tacc > 0 && tacc_cur - tacc > 0)
     n_tot++;

   if(jfrac < 0.5  && ID == ID0 && tacc_cur - tacc < time_avg  &&  tacc > 0 && tacc_cur - tacc > 0)
     n_jlow++;

   if(jfrac >= 0.5  && ID == ID0 && tacc_cur - tacc < time_avg  &&  tacc > 0 && tacc_cur - tacc > 0)
     n_jhigh++;
   }
  fclose(sinkangmom);


  sinkangmom_near = fopen(fsinkangmom_near, "r");
  if(sinkangmom_near != NULL)
    {
    printf("Read in nearby particles!\n");
    for(n = 0; n < All.n_near_max; n++)
      {
      fscanf(sinkangmom_near, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %d\n", &ztime, &d_sink_phys, &vx_phys, &vy_phys, &vz_phys, &jx, &jy, &jz, &jtot, &sinkmass_read, &ID);
      redshift = (1.0/ztime) - 1.0;
      time=5.4e8/pow((1.e0+redshift)/10.e0,1.5e0);
      tacc = time - time0;

      jcent = pow((6.67e-8*sinkmass_read*(1.e10/0.7)*1.98892e33*d_sink_phys),0.5);
      jfrac = jtot/jcent;

      if(ID == ID0)
        jfrac_tot = jfrac_tot+jfrac;

      if(n == All.n_near_max - 20)
        printf("f_tot = %d, n_jlow = %10d n_jhigh = %10d jtot = %10.5g jcent = %10.5g jfrac = %lg, jfrac_tot = %lg, tacc = %lg\n", f_tot, n_jlow, n_jhigh, jtot, jcent, jfrac, jfrac_tot, tacc);

      if(ID == ID0)
        n_tot++;
 
      if(jfrac < 0.5  && ID == ID0)
        n_jlow++;
  
      if(jfrac >= 0.5  && ID == ID0)
        n_jhigh++;
      }
    fclose(sinkangmom_near);
    }


  n_jlow_doub = (double)n_jlow;
  n_tot_doub = (double)n_tot;
  //alpha_val = n_jlow_doub/n_tot_doub;
  alpha_val = 1. - (jfrac_tot/n_tot_doub);
  if(alpha_val > 1)
    alpha_val = 1;
  if(alpha_val < 0)
    alpha_val = 0;
  if(alpha_val != alpha_val)
    alpha_val = All.alpha;

  printf("n_jlow = %10d n_jhigh = %10d n_tot = %10d alpha_val = %10.5g tacc_cur = %lg, tacc = %lg\n", n_jlow, n_jhigh, n_tot, alpha_val, tacc_cur, tacc);

  return(alpha_val); 

  }

#endif

