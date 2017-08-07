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

#define MBIG 10000000
#define MSEED 161803398

CHR *DBL_S[20]=
{ "%le"     ,"%+1.0le" ,"%+2.0le" ,"%+3.0le" ,"%+4.0le" ,"%+5.0le" ,
  "%+6.0le" ,"%+7.0le" ,"%+8.1le" ,"%+9.2le" ,"%+10.3le",
  "%+11.4le","%+12.5le","%+13.6le","%+14.7le","%+15.8le",
  "%+16.9le","%+17.10le","%+18.11le","%+19.12le" 
};

CHR *INT_S[20]=
{"%ld","%1ld","%2ld","%3ld","%4ld","%5ld","%6ld","%7ld","%8ld","%9ld",
 "%10ld","%11ld","%12ld","%13ld","%14ld",
 "%15ld","%16ld" ,"%17ld","%18ld","%19ld"
};


DBL TranNumUniDistGen(INT *i0dum) /* returns uniform distributed random
                                     numbers over interval 0 to 1       */
{ static INT inext;
  static INT inextp;
  static INT ma[56];
  static INT iff=0;
  INT mj;
  INT mk;
  INT i;
  INT ii;
  INT k;
  DBL dranun;
	
  if(*i0dum<0 || iff==0)
  { iff=1;
    mj=MSEED-(*i0dum<0 ? -*i0dum : *i0dum);
    mj%=MBIG;
    ma[55]=mj;
    mk=1;
    for(i=1;i<=54;i++)
    { ii=(21*i)%55;
      ma[ii]=mk;
      mk=mj-mk;
      if(mk<0) mk=mk+MBIG;
      mj=ma[ii];
    }
    for(k=1;k<=4;k++)
    { for(i=1;i<=55;i++)
      { ma[i]=ma[i]-ma[1+(i+30)%55];
        if(ma[i]<0) ma[i]=ma[i]+MBIG;
    } }
    inext=0;
    inextp=31;
    *i0dum=1;
  }
  if(++inext==56) inext=1;
  if(++inextp==56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if(mj<0) mj=mj+MBIG;
  ma[inext]=mj;
  dranun=((DBL)mj)/MBIG;
  return dranun;
}

#define MDIV 1000
#define DDIV 1000.0

DBL TranNumGaussDisGen(INT *i0dum) /* returns Normal distrition (sigma =1/3;) 
                                      random numbers over interval -1 to 1    */
{ INT idiv;
  INT k;
  INT iter;
  INT idelk;
  static INT ipass=0;

  static DBL dxref[MDIV+1];
  static DBL dsref[MDIV+1];
  static DBL dtarea;
  DBL ddiv;
  DBL dsumcorr;
  DBL dsigci;
  DBL dsum;
  DBL dsumce;
  DBL deltx;
  DBL dcx;
  DBL dxmin;
  DBL dxmax;
  DBL dymin;
  DBL dymax;
  DBL dxcen;
  DBL dcxp;
  DBL dsump;
  DBL dsk;
  DBL dleni;
  DBL dlenc;
  DBL dxsec;
  DBL dysec;
  
  idelk=3;
  ddiv=DDIV;
  if(ipass==0)
  { dsigci=D1;
    deltx=D2*dsigci/ddiv;
    dcx=-dsigci;
    dtarea=D0;
    for(idiv=0;idiv<MDIV;idiv++) /* Scale confidence interval   */
    { dtarea=dtarea+(TDBLNormalDist(dcx+(D1-DGP2)*DP5*deltx)+
                     TDBLNormalDist(dcx+(D1+DGP2)*DP5*deltx));
      dcx=dcx+deltx;
    }
    dtarea=dtarea*DP5*deltx;
    dxref[0]=-dsigci; /* For scaled cumulative get corresponding x */
    dsref[0]=D0;
    dcxp=-dsigci;
    dsump=D0;
    for(idiv=1;idiv<(MDIV+1);idiv++)
    { dsumcorr=TDBLNormalise(idiv,ddiv);
      dsref[idiv]=dsumcorr;
      dcx=dcxp;
      dsum=dsump;
      for(k=0;k<MDIV;k++)
      { if(dsum<dsumcorr)
        { dsum=dsum+DP5*deltx*(TDBLNormalDist(dcx+(D1-DGP2)*DP5*deltx)+
                               TDBLNormalDist(dcx+(D1+DGP2)*DP5*deltx))/dtarea;
          dcx=dcx+deltx;
          dsump=dsum;
          dcxp=dcx;
      } }
      dxmin=dcx-deltx;
      dxmax=dcx;
      for(iter=0;iter<32;iter++)
      { dxcen=DP5*(dxmax+dxmin);
        dsumce=dsum-DP5*(dxmax-dxcen)*(TDBLNormalDist(dxcen+(D1-DGP2)*DP5*(dxmax-dxcen))+
                                       TDBLNormalDist(dxcen+(D1+DGP2)*DP5*(dxmax-dxcen)))/dtarea;
        if(dsumce<dsumcorr)
        { dxmin=dxcen;
        }
        else
        { dxmax=dxcen;
          dsum=dsumce;
      } }
      dxref[idiv]=dxcen;
    }
    ipass=1; 
  }

  dsumcorr=TranNumUniDistGen(i0dum);/* Get UNI number           */
  k=TINTDenormalise(dsumcorr,ddiv);/* Find x for UNI number     */
  k=k+idelk;
  dsk=TDBLNormalise(k,ddiv);
  while(dsk>dsumcorr)
  { k=k-1;
    dsk=TDBLNormalise(k,ddiv);
  }
  dxmin=dxref[k];
  dxmax=dxref[k+1];
  dymin=TDBLNormalise(k,ddiv);
  dymax=TDBLNormalise(k+1,ddiv);
  dleni=dxmax-dxmin;
  dlenc=dleni;

  dxsec=D0;
  while((dlenc/dleni)>0.001)
  { dxsec=dxmin+(dxmax-dxmin)*((dsumcorr-dymin)/(dymax-dymin+EPSILON));
    dlenc=MINIM(DABS(dxsec-dxmin),DABS(dxsec-dxmax));
    dxcen=DP5*(dxsec+dxmin);
    dysec=dymin+DP5*(dxsec-dxmin)*(TDBLNormalDist(dxcen+(D1-DGP2)*DP5*(dxsec-dxmin))+
                                   TDBLNormalDist(dxcen+(D1+DGP2)*DP5*(dxsec-dxmin)))/dtarea;
    if(dysec<dsumcorr)
    { dxmin=dxsec;
      dymin=dysec;
    }
    else
    { dxmax=dxsec;
      dymax=dysec;
    }
  }
  return dxsec;
}

/**************SAVING SPACE BY CODED ARRAYS**************/
#define ICODEBASE 90
static CHR *c1diga="0123456789-=qwertyuiop[]^asdfg";
static CHR *c1digb="hjkl;zxcvbnm,./~!@#$%&*()_+QWE";
static CHR *c1digc="RTYUIOP{}|ASDFGHJKL:ZXCVBNM<>?";
static CHR  c1dig[ICODEBASE+1];
static INT  i1dig[1000];
static INT  i1max[10];
static INT iffirst=0;

static void initcode()
{ INT i,j;

  iffirst=1;
  i1max[0]=1;
  for(i=1;i<10;i++)
  { i1max[i]=i1max[i-1]*ICODEBASE;
  }
  for(j=0;j<1000;j++)
  { i1dig[j]=0;
  }    
  CHRcpy(c1dig,c1diga);
  CHRcat(c1dig,c1digb);
  CHRcat(c1dig,c1digc);
  for(i=0;i<ICODEBASE;i++)
  { j=(INT)c1dig[i];
    if(j<0)j=j+500;
    if(j>1000)
    { CHRw(stderr,"Unexpected digit value");
      exit(1);
    }
    i1dig[j]=i;
} }
  
void codeCHRtoINT(CHR *c1code, INT *i1num)
{ INT inum,idig,ndigit,nnum,icar;

  if(iffirst==0)initcode();
  icar=(INT)c1code[0];
  if(icar<0)icar=icar+500;
  ndigit=i1dig[icar];
  icar=(INT)c1code[1];
  if(icar<0)icar=icar+500;
  nnum=i1dig[icar];
  i1num[0]=ndigit;
  i1num[1]=nnum; 
  for(inum=2;inum<nnum;inum++)
  { i1num[inum]=0;
    for(idig=ndigit-1;idig>=0;idig--)
    { icar=(INT)c1code[(inum-2)*ndigit+2+idig];
      if(icar<0)icar=icar+500;
      i1num[inum]=i1num[inum]+i1dig[icar]*i1max[idig];
} } }
 
void codeINTtoCHR(CHR *c1code, INT *i1num)
{ INT inum,idig,ival,ndigit,nnum;

  if(iffirst==0)initcode();
  ndigit=i1num[0];
  nnum=i1num[1];
  c1code[0]=c1dig[ndigit];
  c1code[1]=c1dig[nnum];
  c1code[(nnum-2)*ndigit+2]=CHRTERMINATE; 
  for(inum=2;inum<nnum;inum++)
  { ival=i1num[inum];
    for(idig=ndigit-1;idig>=0;idig--)
    { c1code[(inum-2)*ndigit+2+idig]=c1dig[(ival/i1max[idig])];
      ival=ival%i1max[idig];
} } }

void ZCHRcat(CHR *c1, INT inum)
{ INT i;
  INT j;
  INT iparse;
  CHR c1tmp[5];
  CHR c1num[5];

  for(i=0;i<5;i++)
  { c1num[i]='0';
  }

  SINTw(c1tmp,inum,0);
  iparse=0;
  while(c1tmp[iparse]!='\0')
  { iparse=iparse+1;
  }

  j=4;
  for(i=iparse;i>=0;i--)
  { c1num[j]=c1tmp[i];
    j=j-1;
  }
  CHRcat(c1,c1num);
}

INT Gcd(INT m, INT n)
{ 
  if(n==0)
  {
	  return m;
  }
  else 
  {
	  return Gcd(n, m%n);
  }
}

