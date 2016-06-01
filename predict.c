#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"


/*! \file predict.c
 *  \brief drift particles by a small time interval
 *
 *  This function contains code to implement a drift operation on all the
 *  particles, which represents one part of the leapfrog integration scheme.
 */


/*! This function drifts all particles from the current time to the future:
 *  time0 - > time1
 *
 *  If there is no explicit tree construction in the following timestep, the
 *  tree nodes are also drifted and updated accordingly. Note: For periodic
 *  boundary conditions, the mapping of coordinates onto the interval
 *  [0,All.BoxSize] is only done before the domain decomposition, or for
 *  outputs to snapshot files.  This simplifies dynamic tree updates, and
 *  allows the domain decomposition to be carried out only every once in a
 *  while.
 */
void move_particles(long long int time0, long long int time1)
{
  FILE *AccData;
  char fAccData[200];
  int i, j, do_print = 0;
  double dt_drift, dt_gravkick, dt_hydrokick, dt_entr;
  double t0, t1;
#ifdef CHEMCOOL
  double t2, t3;
#endif
  double a3, a3inv, SinkCriticalDensity,hubble_param, hubble_param2; 
  double gam_fac, pressure_convert, nh_here;

  sprintf(fAccData,"%s/AccData_hiacc", All.OutputDir);

  t0 = second();

  if(All.ComovingIntegrationOn)
    {
      dt_drift = get_drift_factor(time0, time1);
      dt_gravkick = get_gravkick_factor(time0, time1);
      dt_hydrokick = get_hydrokick_factor(time0, time1);
      a3 = All.Time * All.Time * All.Time;
      a3inv = 1.e0/a3;
      if(ThisTask == 0)
        printf("dt_drift = %lg, dt_grav = %lg, dt_hydro = %lg\n", dt_drift, dt_gravkick, dt_hydrokick);
    }
  else
    {
      dt_drift = dt_gravkick = dt_hydrokick = (time1 - time0) * All.Timebase_interval;
      a3 = a3inv = 1.0;
    }

  hubble_param = All.HubbleParam;
  hubble_param2 = hubble_param*hubble_param;
  SinkCriticalDensity = All.SinkCriticalDens / All.UnitDensity_in_cgs * PROTONMASS / HYDROGEN_MASSFRAC * a3 / hubble_param2;

#ifdef CHEMCOOL
  t2 = second();
#endif

  for(i = 0; i < NumPart; i++)
    {

      if(P[i].Type == 0 && P[i].ID < 0) /*SINK*/
	continue;
      for(j = 0; j < 3; j++)
	P[i].Pos[j] += P[i].Vel[j] * dt_drift;

      if(P[i].Type == 0 && SphP[i].sink >0)
        {
         set_sink(a3, a3inv, hubble_param, hubble_param2, i);
         if(All.NumCurrentTiStep % 10000 == 0)
           printf("predict sink = %lg, acc_dir[0] = %lg, acc_dir[1] = %lg, acc_dir[2] = %lg, acc[0] = %lg, acc[1] = %lg, acc[2] = %lg\n", SphP[i].sink, SphP[i].HydroAccel[0], SphP[i].HydroAccel[1], SphP[i].HydroAccel[2], P[i].GravAccel[0], P[i].GravAccel[1], P[i].GravAccel[2]);
         continue;
        }

      if(P[i].Type == 0)
	{
          gam_fac = pow(All.Time,3*(5.0/3.0) -2)/pow(All.Time,3*(SphP[i].Gamma) - 2);
#ifdef PMGRID
	  for(j = 0; j < 3; j++)
	    SphP[i].VelPred[j] +=
	      (P[i].GravAccel[j] + P[i].GravPM[j]) * dt_gravkick + SphP[i].HydroAccel[j] * dt_hydrokick * gam_fac;
#else
	  for(j = 0; j < 3; j++)
	    SphP[i].VelPred[j] += P[i].GravAccel[j] * dt_gravkick + SphP[i].HydroAccel[j] * dt_hydrokick * gam_fac;
#endif
	  SphP[i].Density *= exp(-SphP[i].DivVel * dt_drift);
	  SphP[i].Hsml *= exp(0.333333333333 * SphP[i].DivVel * dt_drift);

	  if(SphP[i].Hsml < All.MinGasHsml)
	    SphP[i].Hsml = All.MinGasHsml;

          if(SphP[i].Prad > 0)
            {
            if((All.NumCurrentTiStep+10) % 10000 == 0)
              printf("predict acc_dir[0] = %lg, acc_dir[1] = %lg, acc_dir[2] = %lg, acc[0] = %lg, acc[1] = %lg, acc[2] = %lg\n", SphP[i].HydroAccel[0], SphP[i].HydroAccel[1], SphP[i].HydroAccel[2], P[i].GravAccel[0], P[i].GravAccel[1], P[i].GravAccel[2]); 
            }

      
         nh_here = SphP[i].Density*All.UnitDensity_in_cgs*All.HubbleParam*All.HubbleParam/a3*HYDROGEN_MASSFRAC/PROTONMASS;
         pressure_convert = pow(a3, SphP[i].Gamma) / All.UnitPressure_in_cgs / (All.HubbleParam*All.HubbleParam);
 
            

#ifdef POLYTROPE
          SphP[i].Pressure = get_pressure(SphP[i].Density);
#else
	  dt_entr = (double)(time1 - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;
#ifdef CHEMCOOL
          if (All.NeedAbundancesForOutput == 1) {
           /* If dt_entr is positive, we need to evolve the chemical network forward in time
            * to get the required values. On the other hand, if dt_entr is negative (i.e. this
            * particle has already evolved past the output time), then we need to go back in
            * time. Since we can't evolve the network backwards, we instead simply use the
            * values at the output time that we computed & stored earlier.
            * 
            * Note that the timestep is constrained such that we can only have overshot one
            * output time, so the values we stored are certain to be the values for that 
            * output time.
            */
           if (dt_entr > 0 && SphP[i].sink < 0.5) { /* Output is in this particle's future */
             do_chemcool_step(dt_entr, &SphP[i], 0, All.NumCurrentTiStep, ThisTask, P[i].ID);
           }
            SphP[i].Pressure = SphP[i].EntropyOut * pow(SphP[i].Density, SphP[i].Gamma);
         }
         else {
           SphP[i].Pressure = (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, SphP[i].Gamma);
         }
#else
	  SphP[i].Pressure = (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, GAMMA);
#endif /* CHEMCOOL */
#endif /* POLYTROPE */
	}
    }

#ifdef CHEMCOOL
  t3 = second();
#endif 

///////////////////////////////ARS printing out acceleration and velocity data
double fac1=1, fac2=1, vfac=1;

if(All.ComovingIntegrationOn)
    {
    fac1 = 1 / (All.Time * All.Time) / All.UnitTime_in_s;
    vfac = pow(All.Time, -1);
    }

if(All.NumCurrentTiStep == ThisTask+4 && do_print == 1)
  {
  AccData=fopen(fAccData,"a");
  for(i = 0; i < N_gas; i++)
    {
            pressure_convert = pow(a3, SphP[i].Gamma) / All.UnitPressure_in_cgs / (All.HubbleParam*All.HubbleParam);
            nh_here = SphP[i].Density*All.UnitDensity_in_cgs*All.HubbleParam*All.HubbleParam/a3*HYDROGEN_MASSFRAC/PROTONMASS; 
            if(All.ComovingIntegrationOn)
               fac2 = 1 / pow(All.Time, 3 * SphP[i].Gamma - 2) / All.UnitTime_in_s;
            if(nh_here > 1.e7 && SphP[i].sink < 0.5)
               {
               //AccData=fopen(fAccData,"a");
               fprintf(AccData, "%15.11lg %8d %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g\n" , 
                           All.Time, P[i].ID, nh_here, SphP[i].Ray_LW_coeff, 
                           P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], 
                           P[i].Vel[0]*vfac, P[i].Vel[1]*vfac, P[i].Vel[2]*vfac, 
                           SphP[i].HydroAccel[0]*fac2, SphP[i].HydroAccel[1]*fac2, SphP[i].HydroAccel[2]*fac2, 
                           P[i].GravAccel[0]*fac1, P[i].GravAccel[1]*fac1, P[i].GravAccel[2]*fac1);
               //fclose(AccData);
               }

    }
    fclose(AccData);
  }

if(All.NumCurrentTiStep == 3 && do_print == 1)
  {
  AccData=fopen(fAccData,"a");
  for(i = 0; i < N_gas; i++)
    {
            pressure_convert = pow(a3, SphP[i].Gamma) / All.UnitPressure_in_cgs / (All.HubbleParam*All.HubbleParam);
            nh_here = SphP[i].Density*All.UnitDensity_in_cgs*All.HubbleParam*All.HubbleParam/a3*HYDROGEN_MASSFRAC/PROTONMASS;
            if(All.ComovingIntegrationOn)
               fac2 = 1 / pow(All.Time, 3 * SphP[i].Gamma - 2) / All.UnitTime_in_s;
            if(SphP[i].sink > 0.5)
               {
               fprintf(AccData, "%15.11lg %8d %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g\n" ,
                           All.Time, P[i].ID, nh_here, SphP[i].Ray_LW_coeff,
                           P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],
                           P[i].Vel[0]*vfac, P[i].Vel[1]*vfac, P[i].Vel[2]*vfac,
                           SphP[i].HydroAccel[0]*fac2, SphP[i].HydroAccel[1]*fac2, SphP[i].HydroAccel[2]*fac2,
                           P[i].GravAccel[0]*fac1, P[i].GravAccel[1]*fac1, P[i].GravAccel[2]*fac1);
               }

    }
    fclose(AccData);
  }

////////////////////////////////////////////////////////////////////////////////////

  /* if domain-decomp and tree are not going to be reconstructed, update dynamically.  */
  if(All.NumForcesSinceLastDomainDecomp < All.TotNumPart * All.TreeDomainUpdateFrequency)
    {
      for(i = 0; i < Numnodestree; i++)
	for(j = 0; j < 3; j++)
	  Nodes[All.MaxPart + i].u.d.s[j] += Extnodes[All.MaxPart + i].vs[j] * dt_drift;

      force_update_len();

      force_update_pseudoparticles();
    }

  t1 = second();

#ifdef CHEMCOOL
  All.CPU_Chemcool += timediff(t2, t3);
#endif
  All.CPU_Predict += timediff(t0, t1);
}



/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 */
#ifdef PERIODIC
void do_box_wrapping(void)
{
  int i, j;
  double boxsize[3];

  for(j = 0; j < 3; j++)
    boxsize[j] = All.BoxSize;

#ifdef LONG_X
  boxsize[0] *= LONG_X;
#endif
#ifdef LONG_Y
  boxsize[1] *= LONG_Y;
#endif
#ifdef LONG_Z
  boxsize[2] *= LONG_Z;
#endif

  for(i = 0; i < NumPart; i++) 
    {
      if(P[i].Type == 0 && P[i].ID < 0) /*SINK*/
      continue;
    for(j = 0; j < 3; j++)
      {
	while(P[i].Pos[j] < 0)
	  P[i].Pos[j] += boxsize[j];

	while(P[i].Pos[j] >= boxsize[j])
	  P[i].Pos[j] -= boxsize[j];
      }
    }
}
#endif

#ifdef POLYTROPE
FLOAT get_pressure(FLOAT density) 
{
  FLOAT p, p1, p2, d1, d2, dd, logdens;
  int indx;

  /* Rescale from code units to CGS number density of H nuclei */
  density *= All.UnitDensity_in_cgs / ((1.0 + 4.0 *ABHE)*PROTONMASS);
  logdens = log10(density);
  if (logdens <= All.MinTabulatedDensity) {
    /* At densities below the minimum tabulated density, pressure scales as p ~ rho^All.PolyIndexLowDensity */
    p  = All.EOSPressure[0];
    p *= pow((pow(10, logdens) / pow(10, All.MinTabulatedDensity)), All.PolyIndexLowDensity);
  }
  else if (logdens >= All.MaxTabulatedDensity) {
    /* At densities above the maximum tabulated density, pressure scales as p ~ rho^All.PolyIndexHighDensity */
    p  = All.EOSPressure[All.EOSFullTableSize - 1];
    p *= pow((pow(10, logdens) / pow(10, All.MaxTabulatedDensity)), All.PolyIndexHighDensity);
  }
  else {
    indx = floor((logdens - All.MinTabulatedDensity) / All.EOSDensDel);
    d1 = All.EOSDensity[indx];
    d2 = All.EOSDensity[indx+1];
    p1 = All.EOSPressure[indx];
    p2 = All.EOSPressure[indx+1];
    dd = (logdens - d1) / (d2 - d1);
    p  = pow(10, log10(p1) + dd * (log10(p2) - log10(p1)));
  }

  /* Convert p back into code units */
  p /= All.UnitPressure_in_cgs;
  return p;
}

FLOAT get_energy(FLOAT density) 
{
  FLOAT e, e1, e2, d1, d2, dd, logdens;
  int indx;

  /* Rescale from code units to CGS number density of H nuclei */
  density *= All.UnitDensity_in_cgs / ((1.0 + 4.0 *ABHE)*PROTONMASS);
  logdens = log10(density);
  if (logdens <= All.MinTabulatedDensity) {
    /* At densities below the minimum tabulated density, energy scales as e ~ rho^All.PolyIndexLowDensity */
    e  = All.EOSEnergy[0];
    e *= pow((pow(10, logdens) / pow(10, All.MinTabulatedDensity)), All.PolyIndexLowDensity);
  }
  else if (logdens >= All.MaxTabulatedDensity) {
    /* At densities above the maximum tabulated density, energy scales as e ~ rho^All.PolyIndexHighDensity */
    e  = All.EOSEnergy[All.EOSFullTableSize - 1];
    e *= pow((pow(10, logdens) / pow(10, All.MaxTabulatedDensity)), All.PolyIndexHighDensity);
  }
  else {
    indx = floor((logdens - All.MinTabulatedDensity) / All.EOSDensDel);
    d1 = All.EOSDensity[indx];
    d2 = All.EOSDensity[indx+1];
    e1 = All.EOSEnergy[indx];
    e2 = All.EOSEnergy[indx+1];
    dd = (logdens - d1) / (d2 - d1);
    e  = pow(10, log10(e1) + dd * (log10(e2) - log10(e1)));
  }

  /* Convert e back into code units. Note that conversion factor
   * for the internal energy density is the same as that for the
   * pressure
   */
  e /= All.UnitPressure_in_cgs;
  return e;
}


#endif /* POLYTROPE */
