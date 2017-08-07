/* Copyright (C) 2000, Dr. Antonio Munjiza
 *
 * This code is provided as part of the book entitled "The Combined
 * Finite Discrete Element Method". It is distributed WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE. Inclusion of a part or whole of this code into any other
 * commercial or research or other purpose code is not granted without author's
 * written explicit permission.
 * When results using whole or any part of this code
 * are published, Y code must be mentioned and acknowledgement to Dr Munjiza must be made.
 * Should you modify this source code, the Copyright (C) on the mdified code
 * as a whole belongs to Dr. A. Munjiza regardless of the extent or nature
 * of modifications.
 * Copyright (C) to whole of any code containing any part of this code
 * also belongs to Dr. A.Munjiza. 
 * Any code comprising any part of this source code
 * must be called Y program.
 * If you do not agree with this, you are not allowed to do 
 * any modifications to any part of this source code or included 
 * any part of it in any other program.
*/
#ifndef Y_NANOSYSTEM
#define Y_NANOSYSTEM

#include "Y_NANOframe.h"
#include "Y_NANOobj.h"
#include "Y_NANOatArDLJ.h"
#include "Y_NANOatArDBuf14_7.h"
#include "Y_NANOatAr.h"
#include "Y_NANOBnD.h"
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

#define NRADIST 100

class Y_NANOParameters
{ public:
   INT nsteps;          /* MaximumNumberOfTimeSteps             */
   INT ioutf;           /* OutputFrequency                      */
   INT iouti;           /* NumberOfTheLastSimulation (used only
                           if isfirst=0 or if isfirst=1         */
   INT isfirs;          /* 0=ContinueOldSimulationSameT,
                           1=ContinueOldSimulationDifferentT,
                           2=NewSimulationNewT                  */

   DBL a;               /* AtomSpacingInACubicalFaceCentered
                           UniformRaster                        */
   DBL dsdt;            /* SizeTimeStep                         */
   DBL dsizc;           /* SizeCoordinates                      */
   DBL dsizv;           /* SizeVelocities                       */
};

class Y_NANOsystem
{ public:
   Y_NANOsystem()
   { nbound=0;
     nprec=4;
     natoms=0;
     nx=0;
     ny=0;
     nz=0;
     nfreqfollow=10;
     ifquadratic=NO;
     
     dEPa=D0;
     dEPb=D0;
     dpress=D0;
     dEK=D0;
     dkbolt=0.1380658;
     dtime=0;

     v0ArDLJ=&vArDLJ;
     v0ArDBUF14_7=&vArDBUF14_7;
     v0BnD=&vBnD;
   }

   INT nbound;          /* NumberofBOUNDaries                   */

   INT nprec;           /* Precission                           */
   INT natoms;          /* TotalNumberOfAtoms                   */
   INT iconbnd;         /* NumberOfContactsWBoundaries          */
   INT iconato;         /* NumberOfContactsWAtoms               */
   INT nleft;           /* NumberOfAtomsOnTheLeftSphere         */
   INT nright;          /* NumberOfAtomsOnTheRightSphere        */
   INT i;               /* CurrentTimeStep                      */
   INT nx;              /* NumberOfCellsInXDirection            */
   INT ny;              /* NumberOfCellsInYDirection            */
   INT nz;              /* NumberOfCellsInZDirection            */
   INT nfreqfollow;     /* FrequencyToFollowAtoms               */
   INT nsetstofollow;
   INT ifquadratic;
   INT nprntfreq;
   INT iBNDModel;       /* iBNDModel=0 -> 3-9 Potential
                           iBNDModel=1 -> WCA Potential         */
   INT iARARModel;      /* ModelForArgonArgonInteraction
                           iARARModel=0 LennardJones 12-6
                           iARARModel=1 Buffer 14-7             */
  
   DBL dtime;           /* TotalSimulationTime                  */
   DBL deltat;          /* DELTAT                               */
   DBL dEPb;            /* POTentialEnergyWBoundary             */
   DBL dEPa;            /* POTentialEnergyWAtoms                */
   DBL dsphra;          /* RADiusSphere                         */
   DBL dcsizc;          /* SizeCoordinates                      */
   DBL dcsizv;          /* SizeVelocities                       */

/* Thermodynamic Data                                           */

   DBL dtemp;           /* Temperature                          */
   DBL dEK;             /* KINeticENergy                        */
   DBL dpress;          /* PRESSureOnTheBoundary                */
   DBL dkbolt;          /* Boltzmann'sConstant                  */
   DBL ddens;           /* SystemDensity                        */
   DBL dheatcurconvx;   /* HeatCurrentDueToConvectionAlongX     */
   DBL dheatcurconvy;   /* HeatCurrentDueToConvectionAlongY     */
   DBL dheatcurconvz;   /* HeatCurrentDueToConvectionAlongZ     */
   DBL dheatcurcontx;   /* HeatCurrentDueToContactAlongX        */
   DBL dheatcurconty;   /* HeatCurrentDueToContactAlongY        */
   DBL dheatcurcontz;   /* HeatCurrentDueToContactAlongZ        */
   DBL dheatcurx;       /* HeatCurrentTotalAlongX               */
   DBL dheatcury;       /* HeatCurrentTotalAlongY               */
   DBL dheatcurz;       /* HeatCurrentTotalAlongZ               */
   DBL dTxx;            /* StressTensorXX                       */
   DBL dTyy;            /* StressTensorYY                       */
   DBL dTzz;            /* StressTensorZZ                       */
   DBL dTxy;            /* StressTensorXY                       */
   DBL dTxz;            /* StressTensorXZ                       */
   DBL dTyz;            /* StressTensorYZ                       */
   DBL dTyx;            /* StressTensorYX                       */
   DBL dTzx;            /* StressTensorZX                       */
   DBL dTzy;            /* StressTensorZY                       */

/* Initial Minimum Distance Between Atoms                       */

   DBL delta;           /* AtomSpacing                          */

/* Cube Data                                                    */
   DBL dcubcc1x;        /* CubeFace1CoordinateX                 */
   DBL dcubcc1y;        /* CubeFace1CoordinateY                 */
   DBL dcubcc1z;        /* CubeFace1CoordinateZ                 */
   DBL dcubcc2x;        /* CubeFace1NormalX                     */
   DBL dcubcc2y;        /* CubeFace1NormalY                     */
   DBL dcubcc2z;        /* CubeFace1NormalZ                     */
   DBL dcubdist;        /* CubeDistanceBetweenAtoms             */
   DBL dvolume;         /* Volume                               */
   DBL darea;           /* SurfaceArea                          */
   DBL dwalltemp;       /* WallTemperature                      */
   DBL dwalltemp0;      /* ReferenceWallTemperature             */
/* End Cube Data                                                */

   Y_NANOobj *v0hdbnd;     /* Y_NANOpointertoHeaDBouNDary             */
   Y_NANOatArDLJ vArDLJ;   /* Y_NANOATomARgonDataLennardJones         */
   Y_NANOatArDLJ *v0ArDLJ; /* Y_NANOpointertoATomARgonDataLennardJones*/
   Y_NANOatArDBuf14_7  vArDBUF14_7;/* Y_NANOATomARgonDataBuffer14_7    */
   Y_NANOatArDBuf14_7  *v0ArDBUF14_7;/* Y_NANOpointertoATomARgonDataBuffer14_7*/
   Y_NANOParameters vpar;  /* Y_NANOparameters                        */
   Y_NANOBnD vBnD;         /* Y_NANOBoundaryData                      */
   Y_NANOBnD *v0BnD;       /* Y_NANOpointertoBoundaryData             */
   Y_NANOatAr *v0atar;

   Y_NANOobj *v0radist[NRADIST];
   INT  nradist;
   DBL  dmaxradist;
   DBL  dmaxradist2;
   DBL  dselectrad;
   DBL  dselectrad2;


/* Function Pointer for Contact                                 */
   INT (*pt2BNDArContactFunction)(DBL,DBL *,DBL *);
   INT (*pt2ArArContactFunction)(DBL,DBL *,DBL *);

};

#endif

