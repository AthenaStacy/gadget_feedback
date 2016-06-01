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

void ghost()
{

double sim_r, r_here, a, a3, hubble, hubble2, hubble3, nh_here;
int i;

 printf("ray_x = %lg ray_y = %lg ray_z = %lg\n", ray.center_x, ray.center_y, ray.center_z);

 if(All.ComovingIntegrationOn)
   {
     a = All.Time;
     a3 = a*a*a;
     hubble = All.HubbleParam;
     hubble2 = hubble*hubble;
     hubble3 = hubble*hubble2;
   } 
 else
   a = a3 = hubble = hubble2 = hubble3 = 1.0;


 for(i = 0; i < N_gas; i++)
 {
    if(P[i].ID > 0 && P[i].ID != All.ray_center_ID)
    {
    sim_r = sqrt((P[i].Pos[0]-ray.center_x)*(P[i].Pos[0]-ray.center_x)+(P[i].Pos[1]-ray.center_y)*(P[i].Pos[1]-ray.center_y)+(P[i].Pos[2]-ray.center_z)*(P[i].Pos[2]-ray.center_z));
    r_here = sim_r * a / hubble * 1.0e3;
    nh_here = SphP[i].Density*All.UnitDensity_in_cgs*hubble2/a3*HYDROGEN_MASSFRAC/PROTONMASS;

    if(sim_r < 2.0*All.SofteningGas && SphP[i].sink < 0.5)
      {
      SphP[i].sink = -1;
      SphP[i].TracAbund[IHP] = 1.0 / (1.0 + ray.alpha_B_H / ray.HI_ion_rate * r_here * r_here * nh_here);
      SphP[i].TracAbund[IDP] = All.DeutAbund * SphP[i].TracAbund[IHP];
      SphP[i].TracAbund[IHEP] = ABHE / (1.0 + ray.alpha_B_H / ray.HeI_ion_rate * r_here * r_here * nh_here);
      SphP[i].TracAbund[IHEPP] = 0.0;

      printf("Found a ghost!\n");
      }

    }
  }
}

#endif
