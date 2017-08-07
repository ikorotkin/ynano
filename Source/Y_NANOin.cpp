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
#include "Y_NANOin.h"

Y_NANOin::Y_NANOin(CHR *f0name)
{ INT i;                /* Integerdummyforloops                 */
  INT j;                /* Integerdummyforcode                  */
  
  finpu.open(f0name);

  c1dig[0]='0';   c1dig[1]='1';   c1dig[2]='2';   c1dig[3]='3';   c1dig[4]='4';
  c1dig[5]='5';   c1dig[6]='6';   c1dig[7]='7';   c1dig[8]='8';   c1dig[9]='9';
  c1dig[10]='a';  c1dig[11]='b';  c1dig[12]='c';  c1dig[13]='d';  c1dig[14]='e';
  c1dig[15]='f';  c1dig[16]='g';  c1dig[17]='h';  c1dig[18]='i';  c1dig[19]='j';
  c1dig[20]='k';  c1dig[21]='l';  c1dig[22]='m';  c1dig[23]='n';  c1dig[24]='o';
  c1dig[25]='p';  c1dig[26]='q';  c1dig[27]='r';  c1dig[28]='s';  c1dig[29]='t';
  c1dig[30]='u';  c1dig[31]='v';  c1dig[32]='w';  c1dig[33]='x';  c1dig[34]='y';
  c1dig[35]='z';  c1dig[36]='A';  c1dig[37]='B';  c1dig[38]='C';  c1dig[39]='D';
  c1dig[40]='E';  c1dig[41]='F';  c1dig[42]='G';  c1dig[43]='H';  c1dig[44]='I';
  c1dig[45]='J';  c1dig[46]='K';  c1dig[47]='L';  c1dig[48]='M';  c1dig[49]='N';
  c1dig[50]='O';  c1dig[51]='P';  c1dig[52]='Q';  c1dig[53]='R';  c1dig[54]='S';
  c1dig[55]='T';  c1dig[56]='U';  c1dig[57]='V';  c1dig[58]='W';  c1dig[59]='X';
  c1dig[60]='Y';  c1dig[61]='Z';  c1dig[62]='!';  c1dig[63]='(';  c1dig[64]=')';
  c1dig[65]='[';  c1dig[66]=']';  c1dig[67]='{';  c1dig[68]='}';  c1dig[69]='^';
  c1dig[70]='_';  c1dig[71]='~';  c1dig[72]='&';  c1dig[73]='/';  c1dig[74]='|';
  c1dig[75]='=';  c1dig[76]=';';  c1dig[77]=':';  c1dig[78]='?';  c1dig[79]='`';

  nprec=3;

  for(i=0;i<MHASH;i++)
  { i1hat[i]=-1;
  }
  for(i=0;i<MBASE;i++)
  { j=(INT)c1dig[i];
    if(j<0)j=j+MHASH/2;
    i1hat[j]=i;
} }

INT Y_NANOin::IsTag(CHR *c1word, INT *i0)
{ INT itag;             /* IntegerTAGflag                       */
  CHR cdummy;           /* CharacterDUMMY                       */

  cdummy=c1word[*i0];

  if((cdummy!='+')&&(cdummy!='-')&&(cdummy!='.')&&
     (cdummy!=',')&&(cdummy!='#')&&(cdummy!='>')&&
     (cdummy!='"')&&(cdummy!='$')&&(cdummy!='%')&&
     (cdummy!='*')&&(cdummy!='@'))
  { itag=1;
  }
  else
  { itag=0;
  }
  return itag;
}

INT Y_NANOin::CHRtoTYPE(CHR *c1word, INT *i0)
{ INT i;                /* Integerdummyforpointer               */
  INT ibase;            /* Integerdummyforbase                  */
  INT itype;            /* IntegerTYPE                          */
  INT j;                /* Integerdummyforcode                  */

  i=*i0;
  ibase=1;
  itype=0;
  i=i+1;
  while(IsTag(c1word, &i))
  { j=(INT)c1word[i];
    if(j<0)j=j+MHASH/2;
    j=i1hat[j];
    itype=itype+j*ibase;
    ibase=ibase*MBASE;
    i=i+1;
  }
  *i0=i;
  return itype;
}   

INT Y_NANOin::CHRtoINT(CHR *c1word, INT *i0)
{ INT i;                /* Integerdummyforpointer               */
  INT ibase;            /* Integerdummyforbase                  */
  INT j;                /* Integerdummyforcode                  */
  INT ipar;             /* IntegerPARameter                     */
  

  i=*i0;
  ibase=1;
  i=i+1;
  ipar=0;
  while(IsTag(c1word, &i))
  { j=(INT)c1word[i];
    if(j<0)j=j+MHASH/2;
    j=i1hat[j];
    ipar=ipar+j*ibase;
    ibase=ibase*MBASE;
    i=i+1;
  }
  *i0=i;
  return ipar;
}

DBL Y_NANOin::CHRtoDBL(CHR *c1word, INT *i0)
{ INT i;                /* Integerdummyforpointer               */
  LONG ibase;           /* Integerdummyforbase                  */
  INT iexp;             /* IntegerEXPonent                      */
  LONG j;               /* Integerdummyforcode                  */
  INT k;                /* Integerdummyforloops                 */
  DBL dpar;             /* DoublePARameter                      */

  i=*i0;
  ibase=1;
  i=i+1;
  j=(INT)c1word[i];
  if(j<0)j=j+MHASH/2;
  j=i1hat[j];
  iexp=j-MBASE/2;
  i=i+1;
  dpar=0;
  while(IsTag(c1word, &i))
  { j=(INT)c1word[i];
    if(j<0)j=j+MHASH/2;
    j=i1hat[j];
    dpar=dpar+(DBL)(j*ibase);
    ibase=ibase*MBASE;
    i=i+1;
  }
  if(iexp>0)
  { for(k=0;k<iexp;k++)
    { dpar=dpar*(DBL)MBASE;
  } }
  else if(iexp<0)
  { for(k=0;k<(-iexp);k++)
    { dpar=dpar/(DBL)MBASE;
  } }
  *i0=i;
  return dpar;
}


void Y_NANOin::ReadObj(Y_NANOprim *v0prim)
{ INT j;                /* Integerdummyforloops                 */
    
  j=0;
  nipar=0;
  ndpar=0;
  ieoff=ReadFil(c1word);
  if(!ieoff)
  { while(c1word[j]!='>')
    { if(c1word[j]=='<')/* tpye                                 */
      { v0prim->type=CHRtoTYPE(c1word,&j);
      }
      else if(c1word[j]=='+')/* positive integer                */
      { v0prim->i1par[nipar]=CHRtoINT(c1word,&j);
        nipar=nipar+1;
      }
      else if(c1word[j]=='-')/* negative integer                */
      { v0prim->i1par[nipar]=-CHRtoINT(c1word,&j);
        nipar=nipar+1;
      }
      else if(c1word[j]=='.')/* positive double                 */
      { v0prim->d1par[ndpar]=CHRtoDBL(c1word,&j);
        ndpar=ndpar+1;
      }
      else if(c1word[j]==',')/* negative double                 */
      { v0prim->d1par[ndpar]=-CHRtoDBL(c1word,&j);
        ndpar=ndpar+1;
  } } }
  v0prim->nipar=nipar;
  v0prim->ndpar=ndpar;
}

INT Y_NANOin::ReadFil(CHR *c1word)
{ INT i;                /* Integerdummyforloops                 */
  INT iEOF;             /* IntegerEndOfFileflag                 */
  CHR cdummy;           /* CharacterDUMMY                       */
  ofstream ferr;        /* ERRorfileinstance                    */

  finpu.get(cdummy);
  i=0;
  while((cdummy!='<')&&(!finpu.eof()))
  { finpu.get(cdummy);
  }
  if(!finpu.eof())
  { c1word[i]=cdummy;
    i=i+1;
    finpu.get(cdummy);
    while((cdummy!='>')&&(!finpu.eof()))
    { if(cdummy=='<')
      { ferr.open("errorlog.tmp",ios::out);
        ferr.write("Incomplete Object!\n",20);
        ferr.close();
        exit(1);
      }
      if((cdummy!=' ') &&(cdummy!='\n')&&(cdummy!='\f')&&(cdummy!='\r')
       &&(cdummy!='\t')&&(cdummy!='\v')&&(cdummy!='\\')&&(cdummy!='\''))
      { c1word[i]=cdummy;
        i=i+1;
        finpu.get(cdummy);
      }
      else
      { finpu.get(cdummy);
    } }
    if(cdummy=='>')
    { c1word[i]=cdummy;
    }
    else
    { ferr.open("errorlog.tmp",ios::out);
      ferr.write("Incomplete Object!\n",20);
      ferr.close();
      exit(1);
    } 
    iEOF=0;
  }
  else
  { iEOF=1;
  }
  return iEOF;
}

