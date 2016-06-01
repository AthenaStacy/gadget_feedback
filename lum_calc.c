#ifdef RAYTRACE_TG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"
#include "tags.h"

#define ls 3.e10
#define  PI  3.14159265358979323846

double lum_calc(int lum_type, double star_mass, double mdot, double nu_ion, double dt_raytrace) 
{
    FILE *startrack;
    double nu, nuint, numax, nupeak;
    double Ltot, Lacc, Lint, L_zams,  L_deut, L_hay, L_hen; 
    double rad, rad_AD, rad_KH, rad_zams, Teff4, Teff, Mmax, Mzams;
    double Fion, Qion, t_KH, t_KH_check, t_acc, rdot;
    double Fion_zams, Qion_zams, Teff4_zams, Teff_zams;
    double L_sphere, L_disk, rad_sphere, rad_sphere_old, rad_disk, rad_AD_disk, rad_disk_old, rad_KH_disk;
    double p1, p2, a1, a2, q1, q2, q3, b1, b2, b3, r2zams;
    double tcount1=0, tcount2=0, tcount3=0, tcount4=0, tcount5=0, tcount7=0;
    double Time, alpha, star_rad, Tint, r0, m0, r1, m1, mdot1, r2, m2, e2;
    int trans1, trans1a, trans2, set1;

    double sig=5.67e-5;
    double mdot0=4.41e-3;   //in sol masses per year
    double wien=5.88e10;
    double boltz=1.38e-16;
    double planck=6.626e-27;
    double nu_LW=2.7e15;

    double r4(double x, double y, double step, double a, double b, double c, double d, double f(double x, double y, double a, double b, double c, double d));
    double bbfunc(double x, double y, double a, double b, double c, double d);
    double bbfunc2(double x, double y, double a, double b, double c, double d);

    int i, imaxint, j, k, l, m, p, q, r, numint;
    int already_zams=0;
    char fstartrack[200];

    sprintf(fstartrack,"startrack_hiacc");
 
    //if(ThisTask == 0)
      printf("Recalculating Teff and L_tot\n");


    printf("ThisTask = %d, mass = %lg, rad = %lg, rad-rad_old = %lg, rad_KH = %lg, t_KH_check = %lg, Ltot = %lg\n", ThisTask, star_mass, All.star_rad, rad-All.star_rad, rad_KH*6.955e10, t_KH_check, Ltot);
  
    t_KH_check = GRAVITY*pow(star_mass*1.98892e33,2)/(All.star_rad*All.lum_tot);  //star_mass in sol.mass
                                                                                  //star_rad in cm


    //rad_AD, rad_KH, rad_AD_disk, rad_KH_disk in Rsol units
    rad_AD = 49.0 * pow((star_mass),(1.0/3.0)) * pow((mdot/mdot0),(1.0/3.0));
    rad_KH = 3.2e6 * (mdot) / pow((star_mass),2);            //mdot must also be in sol.mass. per year

    rad_AD_disk = 1.72*pow(star_mass,-.333333);
    rad_KH_disk = 1.1e6 * mdot / pow(star_mass,2.0); //just lower disk KH contraction radius by a factor of ~3 (?)

     Mmax = 7.*pow(mdot/1.e-3,.27);	 
     Mzams = 50.*pow(mdot/1.e-3,.645);
     r2zams = 0.28*pow(Mzams,0.61); 


    rad_sphere = rad_AD;
    if(rad_KH < rad_sphere)
         {
         rad_sphere = rad_KH;
         }

    rad_disk = rad_AD_disk;
    if(rad_KH_disk < rad_AD_disk)
         {
         rad_disk = rad_KH_disk;
         }


    rad = (All.alpha*rad_sphere + (1.-All.alpha)*rad_disk)*6.955e10;
    //"rad" is now in cm

    printf("rad_init = %lg, rad_sphere = %lg, rad_disk = %lg, rad_KH_disk = %lg, rad_ad_disk = %lg, rad_KH = %lg, rad_AD = %lg\n", rad, rad_sphere, rad_disk, rad_KH_disk, rad_AD_disk, rad_KH, rad_AD);

    All.Tint = 1.e6*(star_mass/0.3)*pow((rad/6.955e10)/2.4,-1);

    rad_zams = 0.28*pow(star_mass,0.61)*6.955e10;
 
    rdot = (rad - All.star_rad)/dt_raytrace;

    double khfac = 1.0, khfac_exp=1.0;

    //limit contraction rate to the KH contraction rate
    if(rad < rad_zams || (rdot < -khfac * (All.star_rad/t_KH_check) /*&& rad_KH < rad_AD*/))
     {
     rad = All.star_rad - khfac*(All.star_rad/t_KH_check)*dt_raytrace;
     //if(ThisTask == 0)
       printf("time = %lg, rdot = %lg, rad = %lg, rad-rad_old = %lg, rad_KH = %lg, t_KH_check = %lg, Ltot = %lg\n", dt_raytrace, rdot, rad, rad-All.star_rad, rad_KH*6.955e10, t_KH_check, Ltot);
     }

    //limit puff-up once star is greater than 10 solar masses.  At low masses puff-up is allowed...
    if(rdot > khfac_exp * (All.star_rad/t_KH_check) && star_mass > 10)
     {
     rad = All.star_rad + khfac_exp*(All.star_rad/t_KH_check)*dt_raytrace;
     printf("time = %lg, rdot = %lg, rad = %lg, rad-rad_old = %lg, rad_KH = %lg, t_KH_check = %lg, Ltot = %lg\n", dt_raytrace, rdot, rad, rad-All.star_rad, rad_KH*6.955e10, t_KH_check, Ltot);
     }


    if(rad < rad_zams || already_zams > 0)
      {
      rad = rad_zams;
      already_zams = 1;
      }

///////////////////////////////////////////////////////////////////////////////////////////////////////////

    Lacc = GRAVITY * (star_mass*1.98892e33)* mdot;  //lum_tot in cgs
    //Lacc = (-3./7.)*((1./3.) - (7.*All.alpha/3.))*(Lacc/rad)*1.98892e33/3.1557e7;
    Lacc = All.alpha*(Lacc/rad)*1.98892e33/3.1557e7;

    L_zams = 1.4e4*pow(star_mass/10.0,2)*3.839e33;

    L_hay = 4.*PI*pow(rad,2)*sig*pow(4500.,4);

    L_hen = pow(10,3.5)*pow(star_mass/9.,22/5)*pow(All.Teff/1.e4,4/5)*3.839e33;
   
    if(All.Tint >= 2.e6)	 
         L_deut = 1500*3.839e33*(mdot/1.e-3);

    if(L_hay > L_hen && rad > rad_zams)
	  Lint = L_hay;
	 
    if(L_hay <= L_hen && rad > rad_zams)
	   Lint = L_hen;
	 	 
    if(rad <= 1.05*rad_zams || already_zams > 0)
	  Lint = L_zams;

	 
     //Ltot = Lacc + L_KH + L_zams_add;
     Ltot = Lacc + Lint /*+ L_KH + L_deut*/;

    t_KH = GRAVITY*pow(star_mass*1.98892e33,2)/(rad_zams*L_zams);
    t_acc = All.tacc;
    All.t_KH = t_KH;
    All.t_acc = t_acc;

    Teff4 = Ltot / (4.e0 * PI * sig * pow(rad,2));
    All.Teff = pow(Teff4,0.25);

    Teff4_zams=L_zams/(4.e0*PI*sig*pow(rad_zams,2));
    Teff_zams=pow(Teff4_zams,0.25);

    //if(ThisTask == 0)
      printf("t_KH = %lg, t_acc = %lg \n", t_KH, t_acc);
    //if(ThisTask == 50)
      printf("star_mass = %lg, mdot = %lg, rad = %lg, lum_tot = %lg, Lacc = %lg, L_hay = %lg, L_hen = %lg, L_zams = %lg, Teff = %lg Tint = %lg\n", star_mass, mdot, rad, Ltot, Lacc, L_hay, L_hen, L_zams, All.Teff, All.Tint);

   // if(ThisTask == 0)
   //   {
      printf("mdot = %lg, rad = %lg, lum_tot = %lg, Teff = %lg\n", mdot, rad, Ltot, All.Teff);

      startrack = fopen(fstartrack, "w");
      fprintf(startrack, "%15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g  %8d %8d %8d %8d \n", 
            All.Time, All.alpha, All.star_rad, All.lum_tot, All.Tint, 
            All.r0, All.m0, All.r1, All.m1, All.mdot1, All.r2, All.m2, All.e2,
            All.trans1,  All.trans1a, All.trans2, All.set1);
      fclose(startrack);
    //  }

    nupeak=wien*All.Teff;
    numax=1.e2*nupeak;

    Fion=Fion_zams=0;
    Qion=Qion_zams=0;
    j=0;
    nuint=nu_ion;

    if(lum_type == 4)
      {
       numax =  nu_ion;
       nu_ion = nu_LW;
       nuint = nu_LW; 
       }

    while(nuint<= numax){
      nuint=nu_ion+j*(numax-nu_ion)/1.e5;
      Qion=r4(nuint, Qion, (numax-nu_ion)/1.e5, All.Teff, planck, nupeak, boltz, bbfunc2);
      j++;
    }

    Qion=PI*Qion;
    Qion=4.0*PI*pow(rad,2.e0)*Qion;

    j=0;
    nuint=nu_ion;
    while(nuint<= numax){
       nuint=nu_ion+j*(numax-nu_ion)/1.e5;
       Qion_zams=r4(nuint, Qion_zams, (numax-nu_ion)/1.e5, Teff_zams, planck, nupeak, boltz, bbfunc2);
       j++;
     }  

    Qion_zams=PI*Qion_zams;
    ray.Qion_zams=4*PI*pow(rad_zams,2.e0)*Qion_zams;

   //if(ThisTask == 0)
      printf("nu_ion = %lg, Qion = %lg\n", nu_ion, Qion);

    if(lum_type == 0 || lum_type == 4)
      return Qion;
    else if (lum_type ==1)
      {
      All.star_rad = rad;
      return Ltot;
      }
    else if (lum_type ==2)
      {
      return All.Teff;
      }
    else if (lum_type ==3)
      return rad;
    else
     {
     printf("lum_type not properly set!\n");
     exit(0);
     }
}

double r4(double x, double y, double step, double a, double b, double c, double d, double f(double x, double y, double a, double b, double c, double d))
{
double k1, k2, k3, k4;
k1 = step*f(x,y, a, b, c, d);
k2 = step*f(x+ step/2.0, y+k1/2.0, a, b, c, d);
k3 = step*f(x + step/2.0, y+k2/2.0, a, b, c, d);
k4 = step*f(x + step, y + k3, a, b, c, d);
return(y+(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0);
}

double bbfunc2(double nu, double Fion, double temp, double planck, double nupeak, double boltz)
{
double bb;
bb=2*pow(nu, 2.e0)*pow(ls, -2.e0)/(exp(planck*nu/(boltz*temp))-1.e0);
//bb=pow(nu, 2.e0);
return(bb);
}

#endif
