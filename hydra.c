#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include "allvars.h"
#include "proto.h"

/*! \file hydra.c
 *  \brief Computation of SPH forces and rate of entropy generation
 *
 *  This file contains the "second SPH loop", where the SPH forces are
 *  computed, and where the rate of change of entropy due to the shock heating
 *  (via artificial viscosity) is computed.
 */


static double hubble_a, atime, hubble_a2, fac_mu, fac_vsic_fix, a3inv;

#ifdef PERIODIC
static double boxSize, boxHalf;

#ifdef LONG_X
static double boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef LONG_Y
static double boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef LONG_Z
static double boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif
#endif



/*! This function is the driver routine for the calculation of hydrodynamical
 *  force and rate of change of entropy due to shock heating for all active
 *  particles .
 */
void hydro_force(void)
{
  long long ntot, ntotleft;
  int i, j, k, n, ngrp, maxfill, source, ndone;
  int *nbuffer, *noffset, *nsend_local, *nsend, *numlist, *ndonelist;
  int level, sendTask, recvTask, nexport, place;
  double soundspeed_i;
  double tstart, tend, sumt, sumcomm;
  double timecomp = 0, timecommsumm = 0, timeimbalance = 0, sumimbalance;
  MPI_Status status;
#ifdef CHEMCOOL
  FLOAT gamm1;
#endif

#ifdef PERIODIC
  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
#endif
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
#endif
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
#endif
#endif

  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
	+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) + All.OmegaLambda;

      hubble_a = All.Hubble * sqrt(hubble_a);
      hubble_a2 = All.Time * All.Time * hubble_a;

#ifndef CHEMCOOL
      fac_mu  = pow(All.Time, 3 * (GAMMA - 1) / 2) / All.Time;
      fac_vsic_fix = hubble_a * pow(All.Time, 3 * GAMMA_MINUS1);
#endif

      a3inv = 1 / (All.Time * All.Time * All.Time);
      atime = All.Time;
    }
  else
    hubble_a = hubble_a2 = atime = fac_mu = fac_vsic_fix = a3inv = 1.0;


  /* `NumSphUpdate' gives the number of particles on this processor that want a force update */
  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
      if(P[n].ID < 0 || SphP[n].sink > 0) /*SINK*/
	continue;
      if(P[n].Ti_endstep == All.Ti_Current)
	NumSphUpdate++;
    }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);


  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);


  i = 0;			/* first particle for this task */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */
      tstart = second();
      for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeHydro - NTask; i++) {
	if(P[i].ID < 0 || SphP[i].sink > 0) /*SINK*/
	  continue;
	if(P[i].Ti_endstep == All.Ti_Current)
	  {
	    ndone++;

	    for(j = 0; j < NTask; j++)
	      Exportflag[j] = 0;

	    hydro_evaluate(i, 0);

	    for(j = 0; j < NTask; j++)
	      {
		if(Exportflag[j])
		  {
		    for(k = 0; k < 3; k++)
		      {
			HydroDataIn[nexport].Pos[k] = P[i].Pos[k];
			HydroDataIn[nexport].Vel[k] = SphP[i].VelPred[k];
		      }
		    HydroDataIn[nexport].Hsml = SphP[i].Hsml;
		    HydroDataIn[nexport].Mass = P[i].Mass;
		    HydroDataIn[nexport].DhsmlDensityFactor = SphP[i].DhsmlDensityFactor;
		    HydroDataIn[nexport].Density = SphP[i].Density;
		    HydroDataIn[nexport].Pressure = SphP[i].Pressure;
                    HydroDataIn[nexport].Prad = SphP[i].Prad;
                    for(k = 0; k < 3; k++)
                       HydroDataIn[nexport].Prad_dir[k] = SphP[i].Prad_dir[k];
		    HydroDataIn[nexport].Timestep = P[i].Ti_endstep - P[i].Ti_begstep;
#ifdef CHEMCOOL
                    HydroDataIn[nexport].Gamma    = SphP[i].Gamma;
#endif /* CHEMCOOL */

		    /* calculation of F1 */
#ifdef CHEMCOOL
		    soundspeed_i = sqrt(SphP[i].Gamma * SphP[i].Pressure / SphP[i].Density);
                    if(All.ComovingIntegrationOn)
                      {
                        fac_mu  = pow(All.Time, 3 * (SphP[i].Gamma - 1.0) / 2) / All.Time;
		      }
#else
		    soundspeed_i = sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density);
#endif
		    HydroDataIn[nexport].F1 = fabs(SphP[i].DivVel) /
		      (fabs(SphP[i].DivVel) + SphP[i].CurlVel +
		       0.0001 * soundspeed_i / SphP[i].Hsml / fac_mu);

		    HydroDataIn[nexport].Index = i;
		    HydroDataIn[nexport].Task = j;
		    nexport++;
		    nsend_local[j]++;
		  }
	      }
	  }
      }
      tend = second();
      timecomp += timediff(tstart, tend);

      qsort(HydroDataIn, nexport, sizeof(struct hydrodata_in), hydro_compare_key);

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      tstart = second();

      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

      tend = second();
      timeimbalance += timediff(tstart, tend);



      /* now do the particles that need to be exported */

      for(level = 1; level < (1 << PTask); level++)
	{
	  tstart = second();
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeHydro)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&HydroDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A,
				   &HydroDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);

	  /* now do the imported particles */
	  tstart = second();
	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    hydro_evaluate(j, 1);
	  tend = second();
	  timecomp += timediff(tstart, tend);

	  /* do a block to measure imbalance */
	  tstart = second();
	  MPI_Barrier(MPI_COMM_WORLD);
	  tend = second();
	  timeimbalance += timediff(tstart, tend);

	  /* get the result */
	  tstart = second();
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeHydro)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&HydroDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B,
				   &HydroDataPartialResult[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, &status);

		      /* add the result to the particles */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  source = j + noffset[recvTask];
			  place = HydroDataIn[source].Index;

			  for(k = 0; k < 3; k++)
			    SphP[place].HydroAccel[k] += HydroDataPartialResult[source].Acc[k];

#ifndef POLYTROPE
#ifdef CHEMCOOL
			  SphP[place].DtEntropyVisc += HydroDataPartialResult[source].DtEntropy;
#else
			  SphP[place].DtEntropy += HydroDataPartialResult[source].DtEntropy;
#endif /* CHEMCOOL */
#endif /* POLYTROPE */

			  if(SphP[place].MaxSignalVel < HydroDataPartialResult[source].MaxSignalVel)
			    SphP[place].MaxSignalVel = HydroDataPartialResult[source].MaxSignalVel;
			}
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);

	  level = ngrp - 1;
	}

      MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
      for(j = 0; j < NTask; j++)
	ntotleft -= ndonelist[j];
    }

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);

#ifdef CHEMCOOL
  for(i=0; i<N_gas; i++) 
    {
      if(P[i].Ti_endstep == All.Ti_Current) 
        {
          gamm1 = SphP[i].Gamma - 1.0;
	  SphP[i].DtEntropyVisc *= gamm1 / (hubble_a2 * pow(SphP[i].Density, gamm1));
        }
    }
#endif /* CHEMCOOL */

  /* do final operations on results */
  tstart = second();

  for(i = 0; i < N_gas; i++) {
    if(P[i].ID < 0 || SphP[i].sink > 0) /*SINK*/
      continue;
    if(P[i].Ti_endstep == All.Ti_Current)
      {
#ifndef POLYTROPE
#ifndef CHEMCOOL
	SphP[i].DtEntropy *= GAMMA_MINUS1 / (hubble_a2 * pow(SphP[i].Density, GAMMA_MINUS1));
#endif 
#endif
#ifdef SPH_BND_PARTICLES
	if(P[i].ID == 0)
	  {
#ifndef POLYTROPE
#ifdef CHEMCOOL
	    SphP[i].DtEntropyVisc = 0;
#else
	    SphP[i].DtEntropy = 0;
#endif
#endif
	    for(k = 0; k < 3; k++)
	      SphP[i].HydroAccel[k] = 0;
	  }
#endif
      }
  }

  tend = second();
  timecomp += timediff(tstart, tend);

  /* collect some timing information */

  MPI_Reduce(&timecomp, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecommsumm, &sumcomm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.CPU_HydCompWalk += sumt / NTask;
      All.CPU_HydCommSumm += sumcomm / NTask;
      All.CPU_HydImbalance += sumimbalance / NTask;
    }
}


/*! This function is the 'core' of the SPH force computation. A target
 *  particle is specified which may either be local, or reside in the
 *  communication buffer.
 */
void hydro_evaluate(int target, int mode)
{
  int j, k, m, n, startnode, numngb; 
  int set_prad_zero = 0;
  long long int timestep;
  FLOAT *pos, *vel;
  FLOAT mass, h_i, dhsmlDensityFactor, rho, pressure, f1, f2,  prad;
  double acc[3], acc_dir[3], prad_dir[3], maxSignalVel;
  double dx, dy, dz, dvx, dvy, dvz;
  double h_i2, hinv, hinv4;
  double p_over_rho2_i, p_over_rho2_j, prad_rho2_i, prad_rho2_j, soundspeed_i, soundspeed_j;
  double hfc, dwk_i, vdotr, vdotr2, visc, mu_ij, rho_ij, vsig;
  double h_j, dwk_j, r, r2, u, hfc_visc;
#ifndef POLYTROPE
  double dtEntropy;
#endif

#ifndef NOVISCOSITYLIMITER
  double dt;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h_i = SphP[target].Hsml;
      mass = P[target].Mass;
      dhsmlDensityFactor = SphP[target].DhsmlDensityFactor;
      rho = SphP[target].Density;
      pressure = SphP[target].Pressure;
      prad = SphP[target].Prad;
      for(m = 0; m < 3; m++) prad_dir[m] = SphP[target].Prad_dir[m];
      timestep = P[target].Ti_endstep - P[target].Ti_begstep;
#ifdef CHEMCOOL
      soundspeed_i = sqrt(SphP[target].Gamma * pressure / rho);
      if(All.ComovingIntegrationOn)
        {
          fac_mu  = pow(All.Time, 3 * (SphP[target].Gamma - 1.0) / 2) / All.Time;
	}
#else
      soundspeed_i = sqrt(GAMMA * pressure / rho);
#endif
      f1 = fabs(SphP[target].DivVel) /
	(fabs(SphP[target].DivVel) + SphP[target].CurlVel +
	 0.0001 * soundspeed_i / SphP[target].Hsml / fac_mu);
    }
  else
    {
      pos = HydroDataGet[target].Pos;
      vel = HydroDataGet[target].Vel;
      h_i = HydroDataGet[target].Hsml;
      mass = HydroDataGet[target].Mass;
      dhsmlDensityFactor = HydroDataGet[target].DhsmlDensityFactor;
      rho = HydroDataGet[target].Density;
      pressure = HydroDataGet[target].Pressure;
      prad = HydroDataGet[target].Prad;
      for(m = 0; m < 3; m++) prad_dir[m] = HydroDataGet[target].Prad_dir[m];
      timestep = HydroDataGet[target].Timestep;
#ifdef CHEMCOOL
      soundspeed_i = sqrt(HydroDataGet[target].Gamma * pressure / rho);
#else
      soundspeed_i = sqrt(GAMMA * pressure / rho);
#endif
      f1 = HydroDataGet[target].F1;
    }

    //if(prad > 0)
    //  printf("hydra1 mode = %d, Prad = %lg, Prad_dir[0] = %lg, Prad_dir[1] = %lg, Prad_dir[2] = %lg, vel[0] = %lg, vel[1] = %lg, vel[2] = %lg\n", mode, prad, prad_dir[0], prad_dir[1], prad_dir[2], vel[0], vel[1], vel[2]);

  if(set_prad_zero == 1)
   {
   prad = 0;
   for(m = 0; m < 3; m++) prad_dir[m] = 0;
   }

  /* initialize variables before SPH loop is started */
  acc[0] = acc[1] = acc[2] = 0;
  acc_dir[0] = acc_dir[1] = acc_dir[2] = 0;
#ifndef POLYTROPE
  dtEntropy = 0;
#endif
  maxSignalVel = 0;

  p_over_rho2_i = pressure / (rho * rho) * dhsmlDensityFactor;
  prad_rho2_i = prad /  (rho * rho) * dhsmlDensityFactor;
  //prad_rho2_i /= All.UnitPressure_in_cgs;
  h_i2 = h_i * h_i;

  for(m = 0; m < 3; m++)
    //acc_dir[m] = prad_dir[m]*(6.e-18/pow(3.0857e21,2))/PROTONMASS
    //             *(All.UnitMass_in_g/.7)*(All.HubbleParam*All.HubbleParam)/(All.Time*All.Time);
    acc_dir[m] = prad_dir[m]*pow(h_i,2)/mass;
    //acc_dir[m] = prad_dir[m]/(rho*rho)/pow(h_i,4)*mass;

  /* Now start the actual SPH computation for this particle */
  startnode = All.MaxPart;
  do
    {
      numngb = ngb_treefind_pairs(&pos[0], h_i, &startnode);

      for(n = 0; n < numngb; n++)
	{
	  j = Ngblist[n];

          if(P[j].Mass > 1.e-9) 
            {
            printf("hydra.c SKIPPING SINK! mass = %lg\n", P[j].Mass);
            continue;  //ARS skipping the hugest sinks!!
            }   

	  dx = pos[0] - P[j].Pos[0];
	  dy = pos[1] - P[j].Pos[1];
	  dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  find the closest image in the given box size  */
	  if(dx > boxHalf_X)
	    dx -= boxSize_X;
	  if(dx < -boxHalf_X)
	    dx += boxSize_X;
	  if(dy > boxHalf_Y)
	    dy -= boxSize_Y;
	  if(dy < -boxHalf_Y)
	    dy += boxSize_Y;
	  if(dz > boxHalf_Z)
	    dz -= boxSize_Z;
	  if(dz < -boxHalf_Z)
	    dz += boxSize_Z;
#endif
	  r2 = dx * dx + dy * dy + dz * dz;
	  h_j = SphP[j].Hsml;
	  if(r2 < h_i2 || r2 < h_j * h_j)
	    {
	      r = sqrt(r2);
	      if(r > 0)
		{
		  p_over_rho2_j = SphP[j].Pressure / (SphP[j].Density * SphP[j].Density);
                  prad_rho2_j =  SphP[j].Prad / (SphP[j].Density * SphP[j].Density);

 
#ifdef CHEMCOOL
		  soundspeed_j = sqrt(SphP[j].Gamma * p_over_rho2_j * SphP[j].Density);
#else
		  soundspeed_j = sqrt(GAMMA * p_over_rho2_j * SphP[j].Density);
#endif
		  dvx = vel[0] - SphP[j].VelPred[0];
		  dvy = vel[1] - SphP[j].VelPred[1];
		  dvz = vel[2] - SphP[j].VelPred[2];
		  vdotr = dx * dvx + dy * dvy + dz * dvz;

		  if(All.ComovingIntegrationOn)
		    vdotr2 = vdotr + hubble_a2 * r2;
		  else
		    vdotr2 = vdotr;

		  if(r2 < h_i2)
		    {
		      hinv = 1.0 / h_i;
#ifndef  TWODIMS
		      hinv4 = hinv * hinv * hinv * hinv;
#else
		      hinv4 = hinv * hinv * hinv / boxSize_Z;
#endif
		      u = r * hinv;
		      if(u < 0.5)
			dwk_i = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		      else
			dwk_i = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		    }
		  else
		    {
		      dwk_i = 0;
		    }

		  if(r2 < h_j * h_j)
		    {
		      hinv = 1.0 / h_j;
#ifndef  TWODIMS
		      hinv4 = hinv * hinv * hinv * hinv;
#else
		      hinv4 = hinv * hinv * hinv / boxSize_Z;
#endif
		      u = r * hinv;
		      if(u < 0.5)
			dwk_j = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		      else
			dwk_j = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		    }
		  else
		    {
		      dwk_j = 0;
		    }

		  if(soundspeed_i + soundspeed_j > maxSignalVel)
		    maxSignalVel = soundspeed_i + soundspeed_j;

		  if(vdotr2 < 0)	/* ... artificial viscosity */
		    {
#ifdef CHEMCOOL
                    if(All.ComovingIntegrationOn)
                      {
                        fac_mu  = pow(All.Time, 3 * (SphP[j].Gamma - 1.0) / 2) / All.Time;
		      } 
#endif /* CHEMCOOL */
		      mu_ij = fac_mu * vdotr2 / r;	/* note: this is negative! */

                      vsig = soundspeed_i + soundspeed_j - 3 * mu_ij;

		      if(vsig > maxSignalVel)
			maxSignalVel = vsig;

		      rho_ij = 0.5 * (rho + SphP[j].Density);
		      f2 =
			fabs(SphP[j].DivVel) / (fabs(SphP[j].DivVel) + SphP[j].CurlVel +
						0.0001 * soundspeed_j / fac_mu / SphP[j].Hsml);

		      visc = 0.25 * All.ArtBulkViscConst * vsig * (-mu_ij) / rho_ij * (f1 + f2);

		      /* .... end artificial viscosity evaluation */
#ifndef NOVISCOSITYLIMITER
		      /* make sure that viscous acceleration is not too large */
		      dt = imax(timestep, (P[j].Ti_endstep - P[j].Ti_begstep)) * All.Timebase_interval;
		      if(dt > 0 && (dwk_i + dwk_j) < 0)
			{
#ifdef CHEMCOOL
                        if(All.ComovingIntegrationOn)
                          {
                            fac_vsic_fix = hubble_a * pow(All.Time, 3 * (SphP[j].Gamma - 1.0));
		          } 
#endif /* CHEMCOOL */
			  visc = dmin(visc, 0.5 * fac_vsic_fix * vdotr2 /
				      (0.5 * (mass + P[j].Mass) * (dwk_i + dwk_j) * r * dt));
			}
#endif
		    }
		  else
		    visc = 0;

                 // if(prad > 0)
                 //   printf("hydra2 Prad = %lg, Prad_i = %lg, Pres_i = %lg, Prad_j = %lg, Pres_j = %lg, UnitPress = %lg\n", prad, prad_rho2_i, p_over_rho2_i, prad_rho2_j, p_over_rho2_j, All.UnitPressure_in_cgs);


		  p_over_rho2_j *= SphP[j].DhsmlDensityFactor;
                  prad_rho2_j *= SphP[j].DhsmlDensityFactor;


		  hfc_visc = 0.5 * P[j].Mass * visc * (dwk_i + dwk_j) / r;

                  //ARS:  Add radiation pressure Prad here?
		  //hfc = hfc_visc + P[j].Mass * (p_over_rho2_i * dwk_i + p_over_rho2_j * dwk_j) / r;


                 // if(prad > 0)
                 //   printf("hydra4 Prad = %lg, Prad_i = %lg, Pres_i = %lg, Prad_j = %lg, Pres_j = %lg, UnitPress = %lg\n", prad, prad_rho2_i, p_over_rho2_i, prad_rho2_j, p_over_rho2_j, All.UnitPressure_in_cgs);
    
                  hfc = hfc_visc + P[j].Mass * ((p_over_rho2_i+prad_rho2_i) * dwk_i + (p_over_rho2_j+prad_rho2_i) * dwk_j) / r; 

		  acc[0] -= hfc * dx;
		  acc[1] -= hfc * dy;
		  acc[2] -= hfc * dz;

                  //if(fabs(acc_dir[0]) > 0)
                  //  printf("hydra4 acc_dir[0] = %lg, acc_dir[1] = %lg, acc_dir[2] = %lg, acc[0] = %lg, acc[1] = %lg, acc[2] = %lg\n", acc_dir[0], acc_dir[1], acc_dir[2], acc[0], acc[1], acc[2]);


#ifndef POLYTROPE
 		  dtEntropy += 0.5 * hfc_visc * vdotr2;
#endif
		}
	    }
	}
    }
  while(startnode >= 0);

  /* Now collect the result at the right place */
  if(mode == 0)
    {
      //if(prad > 0 && ((All.NumCurrentTiStep+10000) % 10 == 0))
      //  printf("hydra5 acc_dir[0] = %lg, acc_dir[1] = %lg, acc_dir[2] = %lg, acc[0] = %lg, acc[1] = %lg, acc[2] = %lg\n", acc_dir[0], acc_dir[1], acc_dir[2], acc[0], acc[1], acc[2]);

      for(k = 0; k < 3; k++)
	SphP[target].HydroAccel[k] = acc[k] + acc_dir[k];
#ifndef POLYTROPE
#ifdef CHEMCOOL
      SphP[target].DtEntropyVisc = dtEntropy;
#else
      SphP[target].DtEntropy     = dtEntropy;
#endif
#endif
      SphP[target].MaxSignalVel = maxSignalVel;
    }
  else
    {
      //if(prad > 0 && ((All.NumCurrentTiStep+10) % 10000 == 0))
      //  printf("hydra5 acc_dir[0] = %lg, acc_dir[1] = %lg, acc_dir[2] = %lg, acc[0] = %lg, acc[1] = %lg, acc[2] = %lg\n", acc_dir[0], acc_dir[1], acc_dir[2], acc[0], acc[1], acc[2]);

      for(k = 0; k < 3; k++)
	HydroDataResult[target].Acc[k] = acc[k] + acc_dir[k];
#ifndef POLYTROPE
      HydroDataResult[target].DtEntropy = dtEntropy;
#endif
      HydroDataResult[target].MaxSignalVel = maxSignalVel;
    }
}




/*! This is a comparison kernel for a sort routine, which is used to group
 *  particles that are going to be exported to the same CPU.
 */
int hydro_compare_key(const void *a, const void *b)
{
  if(((struct hydrodata_in *) a)->Task < (((struct hydrodata_in *) b)->Task))
    return -1;
  if(((struct hydrodata_in *) a)->Task > (((struct hydrodata_in *) b)->Task))
    return +1;
  return 0;
}
