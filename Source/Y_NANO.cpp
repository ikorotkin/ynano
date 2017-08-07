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
#include "Y_NANOframe.h"
#include "Y_NANOin.h"
#include "Y_NANOout.h"
#include "Y_NANOobj.h"
#include "Y_NANOprim.h"
#include "Y_NANOcontact.h"
#include "Y_NANOmove.h"
#include "Y_NANOconvertfromPrim.h"
#include "Y_NANOconverttoPrim.h"
#include "Y_NANOtoolbox.h"
#include "Y_NANOcontactArAr.h"
#include "Y_NANOmoveAr.h"
#include "Y_NANOArfromPrim.h"
#include "Y_NANOArtoPrim.h"
#include "Y_NANOBndCubfromPrim.h"
#include "Y_NANOBndCubtoPrim.h"
#include "Y_NANObndCubical.h"
#include "Y_NANOcontactBnCubAr.h"
#include "Y_NANOsystem.h"
#include "Y_NANOrandvel.h"
#include "Y_NANOconvert.h"
#include "Y_NANOconvertArDbP.h"
#include "Y_NANOconvertArDIA.h"
#include "Y_NANOconvertArCoVeP.h"
#include "Y_NANOconvertCoVePAr.h"
#include "Y_NANOintpoint.h"
#include "Y_NANOdiapoint.h"
#include "Y_NANOdblpoint.h"
#include "Y_NANOcovepoint.h"
#include "Y_NANOContaiBinList.h"
#include "Y_NANOintegerize.h"
#include "Y_NANOintegerizeAr.h"
#include "Y_NANOSetEquilibrium.h"
#include "TMPPotEnerg.h"
#include "Y_NANOgraph.h"
#include "Y_NANOWriteGraph.h"
#include "TQuadratic.h"
#include <iostream>

/* Functions Registration   */

Y_NANOsystem vsys;                                   /* Y_NANOSYSTEMobject         */
Y_NANOsystem *v0sys=&vsys;                           /* Y_NANOSYSTEMptr            */

static Y_NANOContaiBinList vrcbl;                    /* Y_NANORubbContainerBinList */
Y_NANOContaiBinList *v0rcbl=&vrcbl;                  /* Y_NANORubbContainBinListptr*/
static Y_NANOatAr vratar;                            /* Y_NANORubbATomARgon        */
Y_NANOatAr *v0rar=&vratar;                           /* Y_NANORubbATomARgonptr     */
static 
Y_NANObndCubical vrbncub(Y_NANONULL,Y_NANONULL);           /* Y_NANORubbBouNdCUBical     */
Y_NANObndCubical *v0rcub=&vrbncub;                   /* Y_NANORubbBouNdCUBicptr    */
static Y_NANOintpoint vrintp;                        /* Y_NANORubbINTegerPnt       */
Y_NANOintpoint *v0rintp=&vrintp;                     /* Y_NANORubbINTegerPntptr    */
static Y_NANOdblpoint vrdblp;                        /* Y_NANORubbDouBLePnt        */
Y_NANOdblpoint *v0rdblp=&vrdblp;                     /* Y_NANORubbDouBLePntptr     */
static Y_NANOdiapoint vrdiap;                        /* Y_NANORubbDIAmeterPnt      */
Y_NANOdiapoint *v0rdiap=&vrdiap;                     /* Y_NANORubbDIAmeterPntptr   */
static Y_NANOcovepoint vrcovep;                      /* Y_NANORubbCoordVelPnt      */
Y_NANOcovepoint *v0rcovep=&vrcovep;                  /* Y_NANORubbCoordVelPntptr   */

Y_NANOtoolbox vtoolbox;                              /* TOOLBOXobj              */
Y_NANOcontact contact;                               /* CONTACTobj              */
Y_NANOmove move_;                                     /* MOVEobj                 */
Y_NANOconvertfromPrim convertfromprim;               /* CONVERTFROMPRIMitiveobj */
Y_NANOconverttoPrim converttoprim;                   /* CONVERTTOPRIMitiveobj   */
Y_NANOconvert convert;                               /* CONVERTobjects          */
Y_NANOintegerize integerize;                         /* INTEGERIZEobjects       */


static Y_NANOcontactArAr rcArAr(v0rar,v0rar,Y_NANONULL);/* CONTACTARAR             */
static Y_NANOmoveAr rmAr(v0rar,Y_NANONULL,Y_NANONULL);     /* MOVERARgon              */
static Y_NANOArfromPrim rafp(v0rar,v0rcbl,Y_NANONULL);  /* ARFROMPRIMitive         */
static Y_NANOArtoPrim ratp(v0rar,Y_NANONULL,Y_NANONULL);   /* ARTOPRIMitive           */
static 
Y_NANOBndCubfromPrim rbsfp(v0rcub,Y_NANONULL,Y_NANONULL);  /* BndCubFROMPRIMitive     */
static 
Y_NANOBndCubtoPrim rbstpr(Y_NANONULL,v0rcub,Y_NANONULL);   /* BndCubtoPRIMitive       */
static 
Y_NANOcontactBnCubAr rcSA(v0rcub,v0rar,Y_NANONULL);     /* Y_NANOcontactBnCubAr       */
static
 Y_NANOconvertArDbP rcADP(v0rar,v0rdblp,Y_NANONULL);    /* CONVERTARgonDBlPoint    */
static
 Y_NANOconvertArDIA rcADIP(v0rar,v0rdiap,Y_NANONULL);   /* CONVERTARgonDIAPoint    */
static
 Y_NANOconvertArCoVeP rcACVP(v0rar,v0rcovep,Y_NANONULL);/* CONVERTARgonCoordVelPnt */
static
 Y_NANOconvertCoVePAr rcCVPA(v0rcovep,v0rar,Y_NANONULL);/* CONVERTCoordVelPntARgon */
static
 Y_NANOintegerizeAr rinteAr(v0rar,v0rcbl,Y_NANONULL);   /* INTEgerizeARgon         */

void *v0heaphd;         /* HeadOfHeapSpace                      */
void *v0heaplt;         /* LastOfHeapSpace                      */


static void GenerateInitialRaster(CHR *c1fnam,Y_NANOContaiBinList *v0cbinlist)
{ INT iatoms;           /* NumberOfAtoms                        */
  INT nlsimu;           /* NumberOfLastSimulation               */
  INT isnewsim;         /* IsNewSimulationFlag                  */
  INT nx;               /* NumberOfAtomsXDirection              */
  INT ny;               /* NumberOfAtomsYDirection              */
  INT nz;               /* NumberOfAtomsZDirection              */
  INT i1;               /* LoopVariable                         */
  INT i2;               /* LoopVariable                         */
  INT i3;               /* LoopVariable                         */
  INT iradpos;          /* LoopVariable                         */
  DBL a;                /* HalfSideLengthOfUnitCell             */
  DBL deltax;           /* DeltaXOfTheCube                      */
  DBL deltay;           /* DeltaYOfTheCube                      */
  DBL deltaz;           /* DeltaZOfTheCube                      */
  DBL dminx;            /* MinimumCoordinateXInitialForAtoms    */
  DBL dminy;            /* MinimumCoordinateYInitialForAtoms    */
  DBL dminz;            /* MinimumCoordinateZInitialForAtoms    */
  DBL dmaxx;            /* MaximumCoordinateXInitialForAtoms    */
  DBL dmaxy;            /* MaximumCoordinateYInitialForAtoms    */
  DBL dmaxz;            /* MaximumCoordinateZInitialForAtoms    */
  DBL dpos2;            /* SquareOfThePositionOfTheAtom         */
  CHR c1simfil[100];    /* SimulationFilename                   */
  CHR c1nradist[100];   /* NumberOfRadialDistributionFileName   */
  FILE *fp;
  Y_NANOprim *v0prim;      /* Y_NANOpointertoPRIMitiveobject          */
  Y_NANOobj *v0con;        /* ContactorATom                        */
  Y_NANOrandvel vrandvel;  /* Randomvelocity generator;Replace the old velocity of atom Argon with new random velocities             */


  CHRcpy(c1simfil,c1fnam);
  CHRcpy(c1nradist,c1fnam);
  CHRcat(c1nradist,".radist");
  iatoms=0;
  nlsimu=v0sys->vpar.iouti;
  isnewsim=v0sys->vpar.isfirs;
  a=v0sys->vpar.a;
  v0prim=new Y_NANOprim();
 
 if(isnewsim==2)  /* Creates a completely new set of data
                           (new coordinates and new velocities)
                           using a face centered uniform rasted */
  { nlsimu=0;

/* Generate The parameteres for the Cube                        */
    v0prim->type=12;
    v0prim->nipar=0;
    v0prim->ndpar=7;
    v0prim->d1par[0]=v0sys->dcubcc1x;
    v0prim->d1par[1]=v0sys->dcubcc1y;
    v0prim->d1par[2]=v0sys->dcubcc1z;
    v0prim->d1par[3]=v0sys->dcubcc2x;
    v0prim->d1par[4]=v0sys->dcubcc2y;
    v0prim->d1par[5]=v0sys->dcubcc2z;
    v0prim->d1par[6]=v0sys->dcubdist;

    deltax=v0sys->dcubcc2x-v0sys->dcubcc1x;
    deltay=v0sys->dcubcc2y-v0sys->dcubcc1y;
    deltaz=v0sys->dcubcc2z-v0sys->dcubcc1z;

    vtoolbox.Y_NANOtool(&convertfromprim,v0prim);
    v0sys->dvolume=deltax*deltay*deltaz;
    v0sys->darea=D2*(deltax*deltay+deltax*deltaz+deltay*deltaz);

/* Generate Atoms Parameters                                    */
    v0prim->type=0;
    v0prim->nipar=1;
    v0prim->ndpar=11;
  
    nx=v0sys->nx;
    ny=v0sys->ny;
    nz=v0sys->nz;

    dminx=-((DBL)nx)*a;
    dmaxx=((DBL)nx)*a;
    dminy=-((DBL)ny)*a;
    dmaxy=((DBL)ny)*a;
    dminz=-((DBL)nz)*a;
    dmaxz=((DBL)nz)*a;

    v0prim->d1par[3]=D0;
    v0prim->d1par[4]=D0;
    v0prim->d1par[5]=D0;
    v0prim->d1par[6]=D0;
    v0prim->d1par[7]=D0;
    v0prim->d1par[8]=D0;
    v0prim->d1par[9]=D0;
    v0prim->d1par[10]=D0;

    iradpos=0;

/* Generate the Atoms inside the cube                           */
    for(i3=0;i3<=nz;i3++)
    { v0prim->d1par[2]=dminz+((DBL)i3)*D2*a;
      for(i2=0;i2<=ny;i2++)
      { v0prim->d1par[1]=dminy+((DBL)i2)*D2*a;
        for(i1=0;i1<=nx;i1++)
        { v0prim->d1par[0]=dminx+((DBL)i1)*D2*a;
          v0prim->i1par[0]=iatoms;
          vtoolbox.Y_NANOtool(&convertfromprim,v0prim,v0cbinlist);

          V3DLe2(dpos2,v0prim->d1par[0],v0prim->d1par[1],v0prim->d1par[2]);
          if((dpos2<v0sys->dselectrad2)&&(iradpos<NRADIST))
          { v0sys->v0radist[iradpos]=v0cbinlist->v0prev;
            iradpos=iradpos+1; 
          }
          iatoms=iatoms+1;
    } } }

    for(i3=0;i3<=nz;i3++)
    { v0prim->d1par[2]=dminz+((DBL)i3)*D2*a;
      for(i2=0;i2<ny;i2++)
      { v0prim->d1par[1]=dminy+a+((DBL)i2)*D2*a;
        for(i1=0;i1<nx;i1++)
        { v0prim->d1par[0]=dminx+a+((DBL)i1)*D2*a;
          v0prim->i1par[0]=iatoms;
          vtoolbox.Y_NANOtool(&convertfromprim,v0prim,v0cbinlist);

          V3DLe2(dpos2,v0prim->d1par[0],v0prim->d1par[1],v0prim->d1par[2]);
          if((dpos2<v0sys->dselectrad2)&&(iradpos<NRADIST))
          { v0sys->v0radist[iradpos]=v0cbinlist->v0prev;
            iradpos=iradpos+1; 
          }
          iatoms=iatoms+1;
    } } }

    for(i3=0;i3<nz;i3++)
    { v0prim->d1par[2]=dminz+a+((DBL)i3)*D2*a;
      for(i2=0;i2<=ny;i2++)
      { v0prim->d1par[1]=dminy+((DBL)i2)*D2*a;
        for(i1=0;i1<nx;i1++)
        { v0prim->d1par[0]=dminx+a+((DBL)i1)*D2*a;
          v0prim->i1par[0]=iatoms;
          vtoolbox.Y_NANOtool(&convertfromprim,v0prim,v0cbinlist);

          V3DLe2(dpos2,v0prim->d1par[0],v0prim->d1par[1],v0prim->d1par[2]);
          if((dpos2<v0sys->dselectrad2)&&(iradpos<NRADIST))
          { v0sys->v0radist[iradpos]=v0cbinlist->v0prev;
            iradpos=iradpos+1; 
          }
          iatoms=iatoms+1;
    } } }

    for(i3=0;i3<nz;i3++)
    { v0prim->d1par[2]=dminz+a+((DBL)i3)*D2*a;
      for(i2=0;i2<ny;i2++)
      { v0prim->d1par[1]=dminy+a+((DBL)i2)*D2*a;
        for(i1=0;i1<=nx;i1++)
        { v0prim->d1par[0]=dminx+((DBL)i1)*D2*a;
          v0prim->i1par[0]=iatoms;
          vtoolbox.Y_NANOtool(&convertfromprim,v0prim,v0cbinlist);

          V3DLe2(dpos2,v0prim->d1par[0],v0prim->d1par[1],v0prim->d1par[2]);
          if((dpos2<v0sys->dselectrad2)&&(iradpos<NRADIST))
          { v0sys->v0radist[iradpos]=v0cbinlist->v0prev;
            iradpos=iradpos+1; 
          }
          iatoms=iatoms+1;
   } } }

    v0sys->nradist=iradpos;
    fp=fopen(c1nradist,"w");
    CHRw(fp,"NATOMSREFERENCE ");
    INTw(fp,iradpos,5);
    CHRwcr(fp);
    fclose(fp);

/* Apply random velocity field                                  */
    v0con=v0cbinlist->v0next;
    while(v0con!=v0cbinlist)
    { vrandvel.getrandvel(v0con);
      v0con=v0con->v0next;
    }

/* Equilibrate System                                           */
    SetEquilibrium(v0cbinlist);
  }
  v0sys->natoms=iatoms;
}


static void CalculateRadialDistributionFunction(Y_NANOContaiBinList *v0cbinlist,
                                                INT istep)
{ INT i;
  INT iradpos;
  INT nradpos;
  INT icount;
  static DBL d1radpos[30000];
  static Y_NANOWriteGraph *v0wradpos=((Y_NANOWriteGraph *)NULL);
  static INT iout=0;
  static INT nbytes=0;
  DBL drelpos2;
  Y_NANOobj *v0con;
  Y_NANOatAr *v0atar;
  Y_NANOatAr *v0atarREF;
  Y_NANOgraph vgraph1;      /* GraphOutput                          */
  CHR c1tmp[100];

  iradpos=0;

  for(i=0;i<v0sys->nradist;i++)
  { v0atarREF=((Y_NANOatAr *)(v0sys->v0radist[i]));

    v0con=v0cbinlist->v0next;
    while(v0con!=v0cbinlist)
    { v0atar=((Y_NANOatAr *)v0con);
      if(v0atar!=v0atarREF)
      { V3DLe2(drelpos2,(v0atar->daccx-v0atarREF->daccx),
                        (v0atar->daccy-v0atarREF->daccy),
                        (v0atar->daccz-v0atarREF->daccz));
        if(drelpos2<v0sys->dmaxradist2)
        { d1radpos[iradpos]=DSQRT(drelpos2);
          iradpos=iradpos+1;
          if(iradpos>=30000)
          { CHRw(stderr,"Array dimension exceds maximum!");
            CHRwcr(stderr);
            exit(0);
          }
      } }
      v0con=v0con->v0next;
    }
  }
  nradpos=iradpos;

  if((v0wradpos==((Y_NANOWriteGraph *)NULL))||(nbytes>1000000000))
  { if(v0wradpos!=((Y_NANOWriteGraph *)NULL)) delete v0wradpos;

    CHRcpy(c1tmp,"radialpositions");
    ZCHRcat(c1tmp,iout);

    v0wradpos=new Y_NANOWriteGraph(c1tmp);

    iout=iout+1;
    nbytes=0;
  }

/* Integer Data                                                             */
  vgraph1.nipar=0;
/* Double Data                                                              */
  vgraph1.ndpar=nradpos+1;
  vgraph1.d1par[0]=v0sys->dtime;
  icount=1;
  for(iradpos=0;iradpos<nradpos;iradpos++)
  { vgraph1.d1par[icount]=d1radpos[iradpos];
    icount=icount+1;
  }
  v0wradpos->WriteGrph(&vgraph1);
  nbytes=nbytes+(v0sys->nprec+2)*vgraph1.ndpar+5;
}


static void CalculateStressTensor(Y_NANOContaiBinList *v0cbinlist)
{ Y_NANOobj *v0pars;
  Y_NANOatAr *v0atar;

  v0pars=v0cbinlist->v0next;
  while(v0pars!=v0cbinlist)
  { v0atar=(Y_NANOatAr *)v0pars;
    v0sys->dTxx=v0sys->dTxx+v0sys->v0ArDLJ->dama*v0atar->davcxo*v0atar->davcxo;
    v0sys->dTyy=v0sys->dTyy+v0sys->v0ArDLJ->dama*v0atar->davcyo*v0atar->davcyo;
    v0sys->dTzz=v0sys->dTzz+v0sys->v0ArDLJ->dama*v0atar->davczo*v0atar->davczo;
    v0sys->dTxy=v0sys->dTxy+v0sys->v0ArDLJ->dama*v0atar->davcxo*v0atar->davcyo;
    v0sys->dTxz=v0sys->dTxz+v0sys->v0ArDLJ->dama*v0atar->davcxo*v0atar->davczo;
    v0sys->dTyz=v0sys->dTyz+v0sys->v0ArDLJ->dama*v0atar->davcyo*v0atar->davczo;
    v0sys->dTyx=v0sys->dTyx+v0sys->v0ArDLJ->dama*v0atar->davcyo*v0atar->davcxo;
    v0sys->dTzx=v0sys->dTzx+v0sys->v0ArDLJ->dama*v0atar->davczo*v0atar->davcxo;
    v0sys->dTzy=v0sys->dTzy+v0sys->v0ArDLJ->dama*v0atar->davczo*v0atar->davcyo;

    v0pars=v0pars->v0next;
} }


static void CalculateHeatCurrent(Y_NANOContaiBinList *v0cbinlist)
{ Y_NANOobj *v0pars;
  Y_NANOatAr *v0atar;

  v0pars=v0cbinlist->v0next;
  while(v0pars!=v0cbinlist)
  { v0atar=(Y_NANOatAr *)v0pars;
    v0sys->dheatcurconvx=v0sys->dheatcurconvx+
                                    v0atar->davcxo*(v0atar->daeki+v0atar->daepo);
    v0sys->dheatcurconvy=v0sys->dheatcurconvy+
                                    v0atar->davcyo*(v0atar->daeki+v0atar->daepo);
    v0sys->dheatcurconvz=v0sys->dheatcurconvz+
                                    v0atar->davczo*(v0atar->daeki+v0atar->daepo);

    v0atar->davcxo=v0atar->davcx;
    v0atar->davcyo=v0atar->davcy;
    v0atar->davczo=v0atar->davcz;
    v0atar->daeki=D0;
    v0atar->daepo=D0;

    v0pars=v0pars->v0next;
  }

  v0sys->dheatcurx=v0sys->dheatcurconvx+v0sys->dheatcurcontx;
  v0sys->dheatcury=v0sys->dheatcurconvy+v0sys->dheatcurconty;
  v0sys->dheatcurz=v0sys->dheatcurconvz+v0sys->dheatcurcontz;
}

static void ReadDataFile(CHR *c1fnam)
{ CHR c1datfil[100];    /* DataFilename                         */
  CHR c1name[100];      /* VariableName                         */
  FILE *fpdat;          /* DataFile(DataFile)                   */

  CHRcpy(c1datfil,c1fnam);
  CHRcat(c1datfil,".ynano");
  fpdat=fopen(c1datfil,"r");

  CHRr(fpdat,c1name);
  while(FILEND(fpdat)==0)
  { if(CHRcmp(c1name,"IOUTI",5)==0)
    { INTr(fpdat,&v0sys->vpar.iouti);
    }
    else if(CHRcmp(c1name,"ISFIRST",7)==0)
    { INTr(fpdat,&v0sys->vpar.isfirs);
    }
    else if(CHRcmp(c1name,"IFQUADRATIC",11)==0)
    { INTr(fpdat,&v0sys->ifquadratic);
    }
    else if(CHRcmp(c1name,"NPRNTFREQ",9)==0)
    { INTr(fpdat,&v0sys->nprntfreq);
    }
    else if(CHRcmp(c1name,"BNDMODEL",8)==0)
    { INTr(fpdat,&v0sys->iBNDModel);
      if(v0sys->iBNDModel==0)
      { v0sys->pt2BNDArContactFunction=ContactBNDAr39;
      }
      else if(v0sys->iBNDModel==1)
      { v0sys->pt2BNDArContactFunction=ContactBNDArWCA;
      }
      else
      { CHRw(stderr,"Boundary Model Not Valid");
        CHRwcr(stderr);
        exit(1);
      }
    }
    else if(CHRcmp(c1name,"ARARMODEL",9)==0)
    { INTr(fpdat,&v0sys->iARARModel);
      if(v0sys->iARARModel==0)
      { v0sys->pt2ArArContactFunction=ContactArArLJ12_6;
      }
      else if(v0sys->iARARModel==1)
      { v0sys->pt2ArArContactFunction=ContactArArBUF14_7;
      }
      else
      { CHRw(stderr,"Argon-Argon Model Not Valid");
        CHRwcr(stderr);
        exit(1);
      }
    }
    else if(CHRcmp(c1name,"A",1)==0)
    { DBLr(fpdat,&v0sys->vpar.a);
    }
    else if(CHRcmp(c1name,"DMAXRADIST",10)==0)
    { DBLr(fpdat,&v0sys->dmaxradist);
      v0sys->dmaxradist2=v0sys->dmaxradist*v0sys->dmaxradist;
    }
    else if(CHRcmp(c1name,"DSELECTRAD",10)==0)
    { DBLr(fpdat,&v0sys->dselectrad);
      v0sys->dselectrad2=v0sys->dselectrad*v0sys->dselectrad;
    }
    else if(CHRcmp(c1name,"DTEMP",5)==0)
    { DBLr(fpdat,&v0sys->dtemp);
    }
    else if(CHRcmp(c1name,"DWALLTEMP",5)==0)
    { DBLr(fpdat,&v0sys->dwalltemp0);
    }
    else if(CHRcmp(c1name,"CUBCC1X",7)==0)
    { DBLr(fpdat,&v0sys->dcubcc1x);
    }
    else if(CHRcmp(c1name,"CUBCC1Y",7)==0)
    { DBLr(fpdat,&v0sys->dcubcc1y);
    }
    else if(CHRcmp(c1name,"CUBCC1Z",7)==0)
    { DBLr(fpdat,&v0sys->dcubcc1z);
    }
    else if(CHRcmp(c1name,"CUBCC2X",7)==0)
    { DBLr(fpdat,&v0sys->dcubcc2x);
    }
    else if(CHRcmp(c1name,"CUBCC2Y",7)==0)
    { DBLr(fpdat,&v0sys->dcubcc2y);
    }
    else if(CHRcmp(c1name,"CUBCC2Z",7)==0)
    { DBLr(fpdat,&v0sys->dcubcc2z);
    }
    else if(CHRcmp(c1name,"CUBDIST",7)==0)
    { DBLr(fpdat,&v0sys->dcubdist);
    }
    else if(CHRcmp(c1name,"NX",2)==0)
    { INTr(fpdat,&v0sys->nx);
    }
    else if(CHRcmp(c1name,"NY",2)==0)
    { INTr(fpdat,&v0sys->ny);
    }
    else if(CHRcmp(c1name,"NZ",2)==0)
    { INTr(fpdat,&v0sys->nz);
    }
    else if(CHRcmp(c1name,"NFREQFOLLOW",11)==0)
    { INTr(fpdat,&v0sys->nfreqfollow);
    }
    else if(CHRcmp(c1name,"NSETSTOFOLLOW",13)==0)
    { INTr(fpdat,&v0sys->nsetstofollow);
    }
    else
    { CHRw(stderr,"Wrong variable name!!\n");
    }
    CHRr(fpdat,c1name);
  }

  fclose(fpdat);
}


static void ResetCounters()
{ v0sys->iconato=0;
  v0sys->iconbnd=0;
  v0sys->dEPa=D0;
  v0sys->dEPb=D0;
  v0sys->dEK=D0;
  v0sys->dpress=D0;
  v0sys->dheatcurcontx=D0;
  v0sys->dheatcurconty=D0;
  v0sys->dheatcurcontz=D0;
  v0sys->dheatcurconvx=D0;
  v0sys->dheatcurconvy=D0;
  v0sys->dheatcurconvz=D0;
  v0sys->dheatcurx=D0;
  v0sys->dheatcury=D0;
  v0sys->dheatcurz=D0;
  v0sys->dTxx=D0;
  v0sys->dTyy=D0;
  v0sys->dTzz=D0;
  v0sys->dTxy=D0;
  v0sys->dTxz=D0;
  v0sys->dTyz=D0;
  v0sys->dTyx=D0;
  v0sys->dTzx=D0;
  v0sys->dTzy=D0;
}




static void ResolveContact(Y_NANOContaiBinList *v0cbinlist)
{ Y_NANOobj *v0con;        /* ContactorATom                                    */
  Y_NANOobj *v0tbnd;       /* pointertoTargetBouNDary                          */

/* Resolve Contact Between Boundaries and Atoms                             */
  v0tbnd=v0sys->v0hdbnd;
  while(v0tbnd!=Y_NANONULL)
  { v0con=v0cbinlist->v0next;
    while(v0con!=v0cbinlist)
    { vtoolbox.Y_NANOtool(&contact,v0tbnd,v0con);//ContactWithCube(v0bnd,v0at)
      v0con=v0con->v0next;
    }
    v0tbnd=v0tbnd->v0next;
  }

/* Resolve Contact Among Atoms Themselfs                                    */
  if(v0sys->ifquadratic==YES)
  { TQuadratic(v0cbinlist);
  }
  else
  { v0cbinlist->TLContact();
} }


static void PrintOutputParaviewEdgeConn(FILE *ftemp, INT inode,INT jnode)
{ 
  INTw(ftemp,2,19);
  CHRwsp(ftemp);
  CHRwsp(ftemp);
  INTw(ftemp,inode,19);
  CHRwsp(ftemp);
  CHRwsp(ftemp);
  INTw(ftemp,jnode,19);
  CHRwcr(ftemp);
}



static void PrintOutputParaview(INT nout,INT i,INT istart,INT nfrout,
                                CHR *c1fnam,Y_NANOContaiBinList *v0cbinlist)
{ 
  CHR c1worda[100];
  Y_NANOobj *v0con;     /* ContactorATom                                    */
  FILE *ftemp;
  INT iatom;
  INT natoms;
  INT iedge;
  INT ibndvertex;
  Y_NANOatAr *vatar;
  Y_NANObndCubical *v0bnd;
  static INT isfirst=YES;


  CHRw(stdout,"Output :");
  INTw(stdout,nout,4);
  CHRwcr(stdout);

  if(isfirst) /* PLOTTING THE CONTAINER                                     */
  { CHRcpy(c1worda,c1fnam);
    CHRcat(c1worda,"Container.vtk");
    ftemp=fopen(c1worda,"w");
    CHRw(ftemp,"# vtk DataFile Version 2.0");
    CHRwcr(ftemp);
    CHRw(ftemp,"Output for problem: ");
    CHRw(ftemp,c1fnam);
    CHRwcr(ftemp);
    CHRw(ftemp,"ASCII");
    CHRwcr(ftemp);
    CHRwcr(ftemp);
    CHRw(ftemp,"DATASET UNSTRUCTURED_GRID");
    CHRwcr(ftemp);
    CHRw(ftemp,"POINTS ");
    INTw(ftemp,8,19);
    CHRw(ftemp," float");
    CHRwcr(ftemp);

    v0bnd=(Y_NANObndCubical *)v0sys->v0hdbnd;
    DBLw(ftemp,v0bnd->dcc1x,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,v0bnd->dcc1y,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,v0bnd->dcc1z,19);
    CHRwcr(ftemp);
  
    DBLw(ftemp,v0bnd->dcc2x,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,v0bnd->dcc1y,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,v0bnd->dcc1z,19);
    CHRwcr(ftemp);
  
    DBLw(ftemp,v0bnd->dcc2x,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,v0bnd->dcc2y,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,v0bnd->dcc1z,19);
    CHRwcr(ftemp);
   
    DBLw(ftemp,v0bnd->dcc1x,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,v0bnd->dcc2y,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,v0bnd->dcc1z,19);
    CHRwcr(ftemp);
  
    DBLw(ftemp,v0bnd->dcc1x,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,v0bnd->dcc1y,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,v0bnd->dcc2z,19);
    CHRwcr(ftemp);

    DBLw(ftemp,v0bnd->dcc2x,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,v0bnd->dcc1y,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,v0bnd->dcc2z,19);
    CHRwcr(ftemp);

    DBLw(ftemp,v0bnd->dcc2x,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,v0bnd->dcc2y,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,v0bnd->dcc2z,19);
    CHRwcr(ftemp);

    DBLw(ftemp,v0bnd->dcc1x,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,v0bnd->dcc2y,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,v0bnd->dcc2z,19);
    CHRwcr(ftemp);

/* WRITING THE CONNECTIVITIES OF THE CELLS.                                 */
    CHRw(ftemp,"CELLS ");
    INTw(ftemp,12,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    INTw(ftemp,36,19);
    CHRwcr(ftemp);
/* CUBICAL CONTAINER'S EDGES                                                */
    PrintOutputParaviewEdgeConn(ftemp,0,1);
    PrintOutputParaviewEdgeConn(ftemp,1,2);
    PrintOutputParaviewEdgeConn(ftemp,2,3);
    PrintOutputParaviewEdgeConn(ftemp,3,0);
    PrintOutputParaviewEdgeConn(ftemp,4,5);
    PrintOutputParaviewEdgeConn(ftemp,5,6);
    PrintOutputParaviewEdgeConn(ftemp,6,7);
    PrintOutputParaviewEdgeConn(ftemp,7,4);
    PrintOutputParaviewEdgeConn(ftemp,0,4);
    PrintOutputParaviewEdgeConn(ftemp,1,5);
    PrintOutputParaviewEdgeConn(ftemp,2,6);
    PrintOutputParaviewEdgeConn(ftemp,3,7);

/* WRITING THE CELLTYPES OF THE CELLS.                                      */
    CHRw(ftemp,"CELL_TYPES ");
    INTw(ftemp,12,19);
    CHRwcr(ftemp);
    for(iedge=0;iedge<12;iedge++)
    { INTw(ftemp,3,19);
      CHRwcr(ftemp);
    }

/* WRITING THE RADIUS FOR EACH OF THE VERTICES.                             */
    CHRw(ftemp,"POINT_DATA");
    CHRwsp(ftemp);
    INTw(ftemp,8,19);
    CHRwcr(ftemp);
/* WRITING THE VELOCITY FOR EACH OF THE VERTICES. FOR THE CUBICAL CONTAINER
   THE VELOCITY OF THE VERTICES IS SET TO ZERO.                             */
    CHRw(ftemp,"VECTORS Velocity float");
    CHRwcr(ftemp);

    for(ibndvertex=0;ibndvertex<8;ibndvertex++)
    { DBLw(ftemp,D0,19);
      CHRwsp(ftemp);
      CHRwsp(ftemp);
      DBLw(ftemp,D0,19);
      CHRwsp(ftemp);
      CHRwsp(ftemp);
      DBLw(ftemp,D0,19);
      CHRwcr(ftemp);
    }
  
    if(ftemp!=FILENULL) fclose(ftemp);    
    isfirst=NO;
  }


  CHRcpy(c1worda,c1fnam);
  ZCHRcat(c1worda,nout);
  CHRcat(c1worda,".vtk");
  ftemp=fopen(c1worda,"w");
  CHRw(ftemp,"# vtk DataFile Version 2.0");
  CHRwcr(ftemp);
  CHRw(ftemp,"Output for problem: ");
  CHRw(ftemp,c1fnam);
  CHRwcr(ftemp);
  CHRw(ftemp,"ASCII");
  CHRwcr(ftemp);
  CHRwcr(ftemp);
  CHRw(ftemp,"DATASET UNSTRUCTURED_GRID");
  CHRwcr(ftemp);
  CHRw(ftemp,"POINTS ");
  INTw(ftemp,v0sys->natoms,19);
  CHRw(ftemp," float");
  CHRwcr(ftemp);

/* WRITING THE COORDINATES OF THE POINTS.                                   */
  iatom=0;
  v0con=v0cbinlist->v0next;
  while(v0con!=v0cbinlist)
  { vatar=(Y_NANOatAr *)v0con;

    DBLw(ftemp,vatar->daccx,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,vatar->daccy,19);
    CHRwsp(ftemp);  
    CHRwsp(ftemp);
    DBLw(ftemp,vatar->daccz,19);
    CHRwcr(ftemp);

    v0con=v0con->v0next;
    iatom=iatom+1;
  }

  CHRwcr(ftemp);
  natoms=iatom;

/* WRITING THE CONNECTIVITIES OF THE CELLS.                                 */
  CHRw(ftemp,"CELLS ");
  INTw(ftemp,natoms,19);
  CHRwsp(ftemp);
  CHRwsp(ftemp);
  INTw(ftemp,2*natoms,19);
  CHRwcr(ftemp);

  for(iatom=0;iatom<natoms;iatom++)
  { INTw(ftemp,1,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    INTw(ftemp,iatom,19);
    CHRwcr(ftemp);
  }
  CHRwcr(ftemp);

/* WRITING THE CELLTYPES OF THE CELLS.                                      */
  CHRw(ftemp,"CELL_TYPES ");
  INTw(ftemp,natoms,19);
  CHRwcr(ftemp);
  for(iatom=0;iatom<natoms;iatom++)
  { INTw(ftemp,1,19);
    CHRwcr(ftemp);
  }
  CHRwcr(ftemp);

/* WRITING THE RADIUS FOR EACH OF THE VERTICES.                             */
  CHRw(ftemp,"POINT_DATA");
  CHRwsp(ftemp);
  INTw(ftemp,natoms,19);
  CHRwcr(ftemp);
  CHRw(ftemp,"SCALARS Radius float 1");
  CHRwcr(ftemp);
  CHRw(ftemp,"LOOKUP_TABLE default");
  CHRwcr(ftemp);

  v0con=v0cbinlist->v0next;
  while(v0con!=v0cbinlist)
  { vatar=(Y_NANOatAr *)v0con;

    DBLw(ftemp,vatar->drad,19);
    CHRwcr(ftemp);

    v0con=v0con->v0next;
  }
  CHRwcr(ftemp);
  CHRwcr(ftemp);

/* WRITING THE VELOCITY FOR EACH OF THE VERTICES.                           */
  CHRw(ftemp,"VECTORS Velocity float");
  CHRwcr(ftemp);

  v0con=v0cbinlist->v0next;
  while(v0con!=v0cbinlist)
  { vatar=(Y_NANOatAr *)v0con;

    DBLw(ftemp,vatar->davcx,19);
    CHRwsp(ftemp);
    CHRwsp(ftemp);
    DBLw(ftemp,vatar->davcy,19);
    CHRwsp(ftemp);  
    CHRwsp(ftemp);
    DBLw(ftemp,vatar->davcz,19);
    CHRwcr(ftemp);

    v0con=v0con->v0next;
  }

  if(ftemp!=FILENULL) fclose(ftemp);
}

static void ReadAtomsIDToFollow(CHR *c1fnam,INT **i1ids,INT *natoms)
{ INT *i1idsL;          /* ArrayOfAtomsIDsLocal                 */
  INT natomsL;          /* NumberOfAtomsLocal                   */
  INT iatom;            /* LoopVariable                         */
  CHR c1datfil[100];    /* DataFilename                         */
  FILE *fpdat;          /* DataFile(DataFile)                   */

  CHRcpy(c1datfil,c1fnam);
  CHRcat(c1datfil,".ato");
  fpdat=fopen(c1datfil,"r");

  INTr(fpdat,&natomsL);
  i1idsL=((INT *)MALLOC(natomsL*sizeof(INT)));

  for(iatom=0;iatom<natomsL;iatom++)
  { INTr(fpdat,&i1idsL[iatom]);
  }
  fclose(fpdat);

  *i1ids=i1idsL;
  *natoms=natomsL;
}



static void OpenAtomsFiles(CHR *c1fnam,INT natoms)
{ INT iatom;        /* LoopVariable                           */
  CHR c0tmp[100]; 
  FILE *fptmp;      /* TMPVariable                            */

  for(iatom=0;iatom<natoms;iatom++)
  { CHRcpy(c0tmp,c1fnam);
    CHRcat(c0tmp,"atom");
    ZCHRcat(c0tmp,iatom);
    CHRcat(c0tmp,".xyz");
    fptmp=fopen(c0tmp,"w");
    INTw(fptmp,(v0sys->nsetstofollow*500),9);
    CHRwcr(fptmp);
    fprintf(fptmp,CHR_S,"x, y, z, vx, vy, vz, v");
    CHRwcr(fptmp);
    fclose(fptmp);
  }
}

#define NBUFFER 500
#define GETBUFFERCOORDX(d1cfax,ibuffer,iatom,natoms)\
                                   (d1cfax[(ibuffer)*natoms+(iatom)])
#define GETBUFFERCOORDY(d1cfay,ibuffer,iatom,natoms)\
                                   (d1cfay[(ibuffer)*natoms+(iatom)])
#define GETBUFFERCOORDZ(d1cfaz,ibuffer,iatom,natoms)\
                                   (d1cfaz[(ibuffer)*natoms+(iatom)])
#define GETBUFFERVELOCX(d1vfax,ibuffer,iatom,natoms)\
                                   (d1vfax[(ibuffer)*natoms+(iatom)])
#define GETBUFFERVELOCY(d1vfay,ibuffer,iatom,natoms)\
                                   (d1vfay[(ibuffer)*natoms+(iatom)])
#define GETBUFFERVELOCZ(d1vfaz,ibuffer,iatom,natoms)\
                                   (d1vfaz[(ibuffer)*natoms+(iatom)])
#define GETBUFFERVELOC(d1vefa,ibuffer,iatom,natoms)\
                                   (d1vefa[(ibuffer)*natoms+(iatom)])


static void FollowAtoms(DBL *d1cfax,DBL *d1cfay,DBL *d1cfaz,
                        DBL *d1vfax,DBL *d1vfay,DBL *d1vfaz,DBL *d1vefa,
                        INT *i1atomids,INT natomfollow,CHR *c1fnam,
                        Y_NANOContaiBinList *v0cbinlist)
{ INT iatom;        /* LoopVariable                           */
  static INT ibuffer=0;/* BufferCounter                       */
  INT jbuffer;      /* LoopVariable                           */
  DBL dacvel;       /* AtomCurrentVelocity                    */
  CHR c0tmp[100];
  FILE *fptmp;      /* TMPVariable                            */
  Y_NANOobj *v0con;    /* ContactorATom                        */
  Y_NANOatAr *v0atar;  /* ContactorATom                        */ 
  static INT i=0;

  if(i>=(v0sys->nsetstofollow))
  { return;
  }

  for(iatom=0;iatom<natomfollow;iatom++)
  { v0con=v0cbinlist->v0next;
    while(v0con!=v0cbinlist)
    { v0atar=((Y_NANOatAr *)v0con);
      if(i1atomids[iatom]==(v0atar->id))
      { GETBUFFERCOORDX(d1cfax,ibuffer,iatom,natomfollow)=v0atar->daccx;
        GETBUFFERCOORDY(d1cfay,ibuffer,iatom,natomfollow)=v0atar->daccy;
        GETBUFFERCOORDZ(d1cfaz,ibuffer,iatom,natomfollow)=v0atar->daccz;
        GETBUFFERVELOCX(d1vfax,ibuffer,iatom,natomfollow)=v0atar->davcx*1E2;
        GETBUFFERVELOCY(d1vfay,ibuffer,iatom,natomfollow)=v0atar->davcy*1E2;
        GETBUFFERVELOCZ(d1vfaz,ibuffer,iatom,natomfollow)=v0atar->davcz*1E2;

        V3DLen(dacvel,v0atar->davcx*1E2,v0atar->davcy*1E2,
                      v0atar->davcz*1E2);

        GETBUFFERVELOC(d1vefa,ibuffer,iatom,natomfollow)=dacvel;
      }
      v0con=v0con->v0next;
    }
  }
  ibuffer=ibuffer+1;

  if(ibuffer==NBUFFER)
  { for(iatom=0;iatom<natomfollow;iatom++)
    { CHRcpy(c0tmp,c1fnam);
      CHRcat(c0tmp,"atom");
      ZCHRcat(c0tmp,iatom);
      CHRcat(c0tmp,".xyz");
      fptmp=fopen(c0tmp,"a");

      for(jbuffer=0;jbuffer<NBUFFER;jbuffer++)
      {	INTw(fptmp,iatom,5);  
        CHRwsp(fptmp);
        DBLw(fptmp,GETBUFFERCOORDX(d1cfax,jbuffer,iatom,natomfollow),19);
        CHRwsp(fptmp);
        DBLw(fptmp,GETBUFFERCOORDY(d1cfay,jbuffer,iatom,natomfollow),19);
        CHRwsp(fptmp);
        DBLw(fptmp,GETBUFFERCOORDZ(d1cfaz,jbuffer,iatom,natomfollow),19);
        CHRwsp(fptmp);
        DBLw(fptmp,GETBUFFERVELOCX(d1vfax,jbuffer,iatom,natomfollow),19);
        CHRwsp(fptmp);
        DBLw(fptmp,GETBUFFERVELOCY(d1vfay,jbuffer,iatom,natomfollow),19);
        CHRwsp(fptmp);
        DBLw(fptmp,GETBUFFERVELOCZ(d1vfaz,jbuffer,iatom,natomfollow),19);
        CHRwsp(fptmp);
        DBLw(fptmp,GETBUFFERVELOC(d1vefa,jbuffer,iatom,natomfollow),19); 
        CHRwcr(fptmp); 
      }
      fclose(fptmp);
    }
    ibuffer=0;
	i=i+1;
  }
}

static void CheckPotentialAndForce()
{ FILE *fp;
  INT npoint;
  INT ifcontact;
  DBL dr;
  DBL dr2;
  DBL dr3;
  DBL dr4;
  DBL dr6;
  DBL dr7;
  DBL dr9;
  DBL dr10;
  DBL dr12;
  DBL dr13;
  DBL dpot;
  DBL dforc;
  DBL dus;
  DBL dfs;
  DBL dmsqerrU;
  DBL dmsqerrF;

  fp=fopen("check6_12","r");
  CHRw(stderr,"Checking for Correctness of Lennard Jones 12-6 Potential for Ar-Ar");
  CHRwcr(stderr);
  if(fp!=FILENULL)
  { dmsqerrU=D0;
    dmsqerrF=D0;
    npoint=0;

    DBLr(fp,&dr);
    DBLr(fp,&dus);
    DBLr(fp,&dfs);
    while(FILEND(fp)==0)
    { dr2=dr*dr;
      dr6=dr2*dr2*dr2;
      dr7=dr6*dr;
      dr12=dr6*dr6;
      dr13=dr12*dr;

      dpot=PotArAr(dr,dr6,dr12);

      dforc=v0sys->v0ArDLJ->dconstD/dr13-v0sys->v0ArDLJ->dconstE/dr7-
            v0sys->v0ArDLJ->dfcoff;

      dmsqerrU=dmsqerrU+(dpot-dus)*(dpot-dus);
      dmsqerrF=dmsqerrF+(dforc-dfs)*(dforc-dfs);
      npoint=npoint+1;

      DBLr(fp,&dr);
      DBLr(fp,&dus);
      DBLr(fp,&dfs);
    }
    fclose(fp);

    dmsqerrU=dmsqerrU/((DBL)npoint);
    dmsqerrF=dmsqerrF/((DBL)npoint);

    CHRw(stderr,"Mean square error for Potential ");
    DBLw(stderr,dmsqerrU,19);
    CHRwcr(stderr);
    CHRw(stderr,"Mean square error for Force ");
    DBLw(stderr,dmsqerrF,19);
    CHRwcr(stderr);
    fclose(fp);
  }
  else
  { CHRw(stderr,"Could not open check file\n");
  }

  fp=fopen("check3_9","r");
  CHRw(stderr,"Checking for Correctness of Lennard Jones 9-3 Potential for Bnd-Ar");
  CHRwcr(stderr);
  if(fp!=FILENULL)
  { dmsqerrU=D0;
    dmsqerrF=D0;
    npoint=0;

    DBLr(fp,&dr);
    DBLr(fp,&dus);
    DBLr(fp,&dfs);
    while(FILEND(fp)==0)
    { dr2=dr*dr;
      dr3=dr2*dr;
      dr4=dr3*dr;
      dr9=dr3*dr3*dr3;
      dr10=dr9*dr;

      dpot=PotBAr(dr,dr3,dr9);

      FORCEBOUND(dforc,dr10,dr4);

      dmsqerrU=dmsqerrU+(dpot-dus)*(dpot-dus);
      dmsqerrF=dmsqerrF+(dforc-dfs)*(dforc-dfs);
      npoint=npoint+1;

      DBLr(fp,&dr);
      DBLr(fp,&dus);
      DBLr(fp,&dfs);
    }
    fclose(fp);

    dmsqerrU=dmsqerrU/((DBL)npoint);
    dmsqerrF=dmsqerrF/((DBL)npoint);

    CHRw(stderr,"Mean square error for Potential ");
    DBLw(stderr,dmsqerrU,19);
    CHRwcr(stderr);
    CHRw(stderr,"Mean square error for Force ");
    DBLw(stderr,dmsqerrF,19);
    CHRwcr(stderr);
    fclose(fp);
  }
  else
  { CHRw(stderr,"Could not open check file\n");
  }

  fp=fopen("checkbuff14_7","r");
  CHRw(stderr,"Checking for Correctness of Buffer 14-7 Potential for Ar-Ar");
  CHRwcr(stderr);
  if(fp!=FILENULL)
  { dmsqerrU=D0;
    dmsqerrF=D0;
    npoint=0;

    DBLr(fp,&dr);
    DBLr(fp,&dus);
    DBLr(fp,&dfs);
    while(FILEND(fp)==0)
    { dr2=dr*dr;

      ifcontact=ContactArArBUF14_7(dr2,&dpot,&dforc);

      dmsqerrU=dmsqerrU+(dpot-dus)*(dpot-dus);
      dmsqerrF=dmsqerrF+(dforc-dfs)*(dforc-dfs);
      npoint=npoint+1;

      DBLr(fp,&dr);
      DBLr(fp,&dus);
      DBLr(fp,&dfs);
    }
    fclose(fp);

    dmsqerrU=dmsqerrU/((DBL)npoint);
    dmsqerrF=dmsqerrF/((DBL)npoint);

    CHRw(stderr,"Mean square error for Potential ");
    DBLw(stderr,dmsqerrU,19);
    CHRwcr(stderr);
    CHRw(stderr,"Mean square error for Force ");
    DBLw(stderr,dmsqerrF,19);
    CHRwcr(stderr);
    fclose(fp);
  }
  else
  { CHRw(stderr,"Could not open check file\n");
  }

  exit(1);
}



INT main(INT argc, char **argv)
{ INT isfirst=1;        /* Flag                                 */
  INT i;                /* dummyvariableforloops                */
  INT nout;             /* numberofoutputinstances              */
  INT nfrout;           /* OutputFrequency                      */
  INT iseq;             /* ActualSequence                       */
  INT nseq;             /* LastSequence                         */
  INT istart;           /* StartingTimeStep                     */
  INT iend;             /* EndingTimeStep                       */
  INT ntout;            /* TotalNumberOfOutsForSequence         */
  INT icount;           /* Counter                              */
  INT jj;               /* LoopVariable                         */
  INT i1nconta[NBUFFER];/* NumberOfAtomsInContactWithAtoms      */
  INT i1ncontb[NBUFFER];/* NumberOfAtomsInContactWithBoundary   */
  DBL d1EPa[NBUFFER];   /* PotentialEnergyAmongAtoms            */
  DBL d1EPb[NBUFFER];   /* PotentialEnergyAmongAtomsAndBoundary */
  DBL d1EK[NBUFFER];    /* TotalKineticEnergy                   */
  DBL d1press[NBUFFER]; /* TotalPressure                        */
  DBL d1time[NBUFFER];  /* SimulationTime                       */
  DBL d1heatcurconvx[NBUFFER];/* HeatCurrentDueToConvectionX    */
  DBL d1heatcurconvy[NBUFFER];/* HeatCurrentDueToConvectionY    */
  DBL d1heatcurconvz[NBUFFER];/* HeatCurrentDueToConvectionZ    */
  DBL d1heatcurcontx[NBUFFER];/* HeatCurrentDueToContactX       */
  DBL d1heatcurconty[NBUFFER];/* HeatCurrentDueToContactY       */
  DBL d1heatcurcontz[NBUFFER];/* HeatCurrentDueToContactZ       */
  DBL d1heatcurx[NBUFFER];/* HeatCurrentTotalX                  */
  DBL d1heatcury[NBUFFER];/* HeatCurrentTotalY                  */
  DBL d1heatcurz[NBUFFER];/* HeatCurrentTotalZ                  */
  DBL d1Txx[NBUFFER];   /* StressTensorTxx                      */
  DBL d1Tyy[NBUFFER];   /* StressTensorTyy                      */
  DBL d1Tzz[NBUFFER];   /* StressTensorTzz                      */
  DBL d1Txy[NBUFFER];   /* StressTensorTxy                      */
  DBL d1Txz[NBUFFER];   /* StressTensorTxz                      */
  DBL d1Tyz[NBUFFER];   /* StressTensorTyz                      */
  DBL d1Tyx[NBUFFER];   /* StressTensorTyx                      */
  DBL d1Tzx[NBUFFER];   /* StressTensorTzx                      */
  DBL d1Tzy[NBUFFER];   /* StressTensorTzy                      */
  DBL dTemp;            /* Temperature                          */
  DBL dPress;           /* Pressure                             */
  DBL dMDPress;         /* MolecularDynamicsPressure            */

/* Follow Atoms Variables                                       */
  DBL *d1cfax;          /* CoordinateFollowAtomsX               */
  DBL *d1cfay;          /* CoordinateFollowAtomsY               */
  DBL *d1cfaz;          /* CoordinateFollowAtomsZ               */
  DBL *d1vfax;          /* VelocityFollowAtomsX                 */
  DBL *d1vfay;          /* VelocityFollowAtomsY                 */
  DBL *d1vfaz;          /* VelocityFollowAtomsZ                 */
  DBL *d1vefa;          /* VelocityFollowAtoms                  */
  INT *i1atomids;       /* FollowAtomsIds                       */
  INT natomfollow;      /* NumberOfAtomsToFollow                */
/* End Follow Atoms Variables                                   */

  DBL dETOTALinitial=D1;
  DBL dETOTALactual;

  CHR c1fnam[100];      /* problemfilename                      */
  CHR c1seqfil[100];    /* SequenceFilename                     */
  CHR c1codfil[100];    /* CodedFilename                        */
  CHR c1enefil[100];    /* EnergyFilename                       */
  CHR c1fatoms[100];    /* FollowAtomsFilename                  */
  CHR c1tmp[100];

  Y_NANOContaiBinList vcbinlist;   /* ContainerBinalyList          */
  Y_NANOgraph vgraph;      /* GraphOutput                          */
  Y_NANOWriteGraph *v0wgr; /* WriteGraphOutput                     */

  v0heaphd=MALLOC(4000000000);  /* OMG!!! */
  v0heaplt=v0heaphd;

  FILE *fpene;          /* EnergyFile                           */
  FILE *fpseq;          /* SequenceFile(DataFile)               */

/*  CheckPotentialAndForce();*/

/*  CHRcpy(c1fnam,"TEST");*/
 /*  CHRcpy(c1fnam,"TEST");*/
  if(argc>1)
  { CHRcpy(c1fnam,argv[1]);
  }
  else
  { printf("double clik the ynanosetupfile to run the program");
  }
 
  i=0;
  while ((c1fnam[i]!='.')&&(c1fnam[i]!='\0'))
  { c1tmp[i]=c1fnam[i];
     i=i+1;
  }
  c1tmp[i]='\0';
  CHRcpy(c1fnam,c1tmp);

  printf("%s",c1fnam);

  CHRcpy(c1seqfil,c1fnam);
  CHRcpy(c1enefil,c1fnam);
  CHRcpy(c1fatoms,c1fnam);
  CHRcat(c1seqfil,".seq");
  CHRcat(c1enefil,".ene");

  fpene=fopen(c1enefil,"w");
  fclose(fpene);

/* Read The Data File                                                       */
  ReadDataFile(c1fnam);//Read the problem file;intialize system parameters<system.h>

/* Read The File With The IDs Of The Atoms To Follow                        */
  ReadAtomsIDToFollow(c1fatoms,&i1atomids,&natomfollow);//read atoms id to the array ilatomids[]

/* Allocate Memory For Follow Atom Process                                  */
  d1cfax=((DBL *)MALLOC(natomfollow*NBUFFER*sizeof(DBL)));
  d1cfay=((DBL *)MALLOC(natomfollow*NBUFFER*sizeof(DBL)));
  d1cfaz=((DBL *)MALLOC(natomfollow*NBUFFER*sizeof(DBL)));
  d1vfax=((DBL *)MALLOC(natomfollow*NBUFFER*sizeof(DBL)));
  d1vfay=((DBL *)MALLOC(natomfollow*NBUFFER*sizeof(DBL)));
  d1vfaz=((DBL *)MALLOC(natomfollow*NBUFFER*sizeof(DBL)));
  d1vefa=((DBL *)MALLOC(natomfollow*NBUFFER*sizeof(DBL)));//The maxmium number of atoms to follow is 500

/* Initialize The Atoms Files By Opening Them And Closing Them              */
  OpenAtomsFiles(c1fatoms,natomfollow);//Create the follow atoms file

/* Generate Atom Raster                                                     */
  GenerateInitialRaster(c1fnam,&vcbinlist);

  nout=v0sys->vpar.iouti;
  CHRw(stderr,"START"); CHRwcr(stderr);
  CHRw(stderr,"Number of Atoms: ");
  INTw(stderr,v0sys->natoms,4); CHRwcr(stderr);
  CHRw(stderr,"Radius of Sphere: ");
  DBLw(stderr,v0sys->dsphra,15); CHRwcr(stderr);

  nseq=0;
  iseq=-2;
  while(1)
  { fpseq=fopen(c1seqfil,"r");
    while((iseq!=nseq) && (iseq!=-1))
    { INTr(fpseq,&iseq);
      INTr(fpseq,&istart);
      INTr(fpseq,&iend);
      INTr(fpseq,&nfrout);
      DBLr(fpseq,&v0sys->deltat);
    }
    fclose(fpseq);

    if(iseq!=-1)
    { nseq=nseq+1;
      ntout=((iend-istart)/nfrout)-1;
      CHRcpy(c1codfil,c1fnam);
      ZCHRcat(c1codfil, iseq);
      CHRcat(c1codfil,".p");

      v0wgr=new Y_NANOWriteGraph(c1codfil);

      CHRw(stderr,"Running Sequence: "); INTw(stderr,iseq,4);
      CHRw(stderr," Delta t: "); DBLw(stderr,v0sys->deltat,15);
      CHRw(stderr,CHRRETURN);

      icount=0;

      v0sys->nprec=4;
      for(i=istart;i<iend;i++)
      { v0sys->i=i;
/* Print OUTPUT                                                             */
        if(((i-istart)%nfrout)==0)
        { if(((i-istart)/nfrout)==ntout)
          { v0sys->nprec=6;
            PrintOutputParaview(nout,i,istart,nfrout,c1fnam,&vcbinlist);
            v0sys->nprec=4;
          }
          else
          { PrintOutputParaview(nout,i,istart,nfrout,c1fnam,&vcbinlist);
          } 
          nout=nout+1;
        }

/* Follow Atoms                                                             */
        if(((i-istart)%(v0sys->nfreqfollow))==0)
        { FollowAtoms(d1cfax,d1cfay,d1cfaz,d1vfax,d1vfay,d1vfaz,d1vefa,
                    i1atomids,natomfollow,c1fnam,&vcbinlist);
        }

/* Setting al Counters and Acumulators to Zero                              */
        ResetCounters();

/* Resolve Contact Between Boundaries and Atoms                             */
        ResolveContact(&vcbinlist);

/* Calculate Radial Distribution Function                                   */
//        CalculateRadialDistributionFunction(&vcbinlist,i);

/* Calculate Stress Tensor                                                  */
        CalculateStressTensor(&vcbinlist);

/* Calculate Heat Current                                                   */
//        CalculateHeatCurrent(&vcbinlist);

/* Calculate the Initial Total Energy                                       */
        if(isfirst==1)
        { dETOTALinitial=v0sys->dEK+v0sys->dEPa+v0sys->dEPb;
          isfirst=0;
        }

/* Acumulate data in the buffer arrays                                      */
        i1nconta[icount]=v0sys->iconato;
        i1ncontb[icount]=v0sys->iconbnd;
        d1EPa[icount]=v0sys->dEPa;
        d1EPb[icount]=v0sys->dEPb;
        d1EK[icount]=v0sys->dEK;
        d1press[icount]=v0sys->dpress;
        d1time[icount]=v0sys->dtime;
        d1heatcurconvx[icount]=v0sys->dheatcurconvx;
        d1heatcurconvy[icount]=v0sys->dheatcurconvy;
        d1heatcurconvz[icount]=v0sys->dheatcurconvz;
        d1heatcurcontx[icount]=v0sys->dheatcurcontx;
        d1heatcurconty[icount]=v0sys->dheatcurconty;
        d1heatcurcontz[icount]=v0sys->dheatcurcontz;
        d1heatcurx[icount]=v0sys->dheatcurx;
        d1heatcury[icount]=v0sys->dheatcury;
        d1heatcurz[icount]=v0sys->dheatcurz;
        d1Txx[icount]=v0sys->dTxx;
        d1Tyy[icount]=v0sys->dTyy;
        d1Tzz[icount]=v0sys->dTzz;
        d1Txy[icount]=v0sys->dTxy;
        d1Txz[icount]=v0sys->dTxz;
        d1Tyz[icount]=v0sys->dTyz;
        d1Tyx[icount]=v0sys->dTyx;
        d1Tzx[icount]=v0sys->dTzx;
        d1Tzy[icount]=v0sys->dTzy;
        icount=icount+1;

/*  Pressure & Energy Output                                                */
        if((icount%NBUFFER)==0)
        { for(jj=0;jj<icount;jj++)
          { 
/* Integer Data                                                             */
            vgraph.nipar=3;
            vgraph.i1par[0]=v0sys->natoms;
            vgraph.i1par[1]=i1nconta[jj];
            vgraph.i1par[2]=i1ncontb[jj];

/* Double Data                                                              */
            vgraph.ndpar=23;
            vgraph.d1par[0]=d1time[jj];
            vgraph.d1par[1]=d1EPa[jj];
            vgraph.d1par[2]=d1EPb[jj];
            vgraph.d1par[3]=d1EK[jj];
            vgraph.d1par[4]=d1press[jj];
            vgraph.d1par[5]=d1heatcurconvx[jj];
            vgraph.d1par[6]=d1heatcurconvy[jj];
            vgraph.d1par[7]=d1heatcurconvz[jj];
            vgraph.d1par[8]=d1heatcurcontx[jj];
            vgraph.d1par[9]=d1heatcurconty[jj];
            vgraph.d1par[10]=d1heatcurcontz[jj];
            vgraph.d1par[11]=d1heatcurx[jj];
            vgraph.d1par[12]=d1heatcury[jj];
            vgraph.d1par[13]=d1heatcurz[jj];
            vgraph.d1par[14]=d1Txx[jj];
            vgraph.d1par[15]=d1Tyy[jj];
            vgraph.d1par[16]=d1Tzz[jj];
            vgraph.d1par[17]=d1Txy[jj];
            vgraph.d1par[18]=d1Txz[jj];
            vgraph.d1par[19]=d1Tyz[jj];
            vgraph.d1par[20]=d1Tyx[jj];
            vgraph.d1par[21]=d1Tzx[jj];
            vgraph.d1par[22]=d1Tzy[jj];
            v0wgr->WriteGrph(&vgraph);
          }
          icount=0;
        }
        dETOTALactual=v0sys->dEK+v0sys->dEPa+v0sys->dEPb;

        if(((i-istart)%v0sys->nprntfreq)==0)
        { dTemp=D2*v0sys->dEK/(D3*((DBL)v0sys->natoms)*v0sys->dkbolt);
          dPress=v0sys->dpress*1e8/v0sys->darea;
          dMDPress=1e8*(v0sys->dTxx+v0sys->dTyy+v0sys->dTzz)/(D3*v0sys->dvolume);
          CHRw(stdout,"RUNNING STEP: ");
          INTw(stdout,i,5);
          CHRwcr(stdout);
          CHRw(stdout,"WALL TEMPERATURE: ");
          DBLw(stdout,v0sys->dwalltemp,19);
          CHRwsp(stdout);
          CHRw(stdout,"Energy Difference: ");
          DBLw(stdout,((dETOTALactual-dETOTALinitial)*D100)/dETOTALinitial,19);
          CHRwcr(stdout);
          CHRw(stdout,"TEMPERATURE: ");
          DBLw(stdout,dTemp,19);
          CHRwcr(stdout);
          CHRw(stdout,"DEM PRESSURE: ");
          DBLw(stdout,dPress,19);
          CHRw(stdout,"  MD PRESSURE: ");
          DBLw(stdout,dMDPress,19);
          CHRwcr(stdout);
        }

        v0sys->dtime=v0sys->dtime+v0sys->deltat;
      }

      dETOTALactual=v0sys->dEK+v0sys->dEPa+v0sys->dEPb;
/*
      fpene=fopen(c1enefil,"a");
      CHRw(fpene,"Percentage Difference in energy: ");
      DBLw(fpene,((dETOTALactual-dETOTALinitial)*D100)/dETOTALinitial,15);
      CHRw(fpene,CHRRETURN);
      fclose(fpene);
*/
      delete v0wgr;
    }
    else if(iseq==(-1))
    { return(1);
      exit(0);
  } }

  return(1);
}

