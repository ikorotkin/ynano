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
#include "Y_NANOWriteGraph.h"

Y_NANOWriteGraph::Y_NANOWriteGraph(CHR *f0name)
{ INT i;                /* dummyvariableforloops                */
  INT j;                /* dummyvariableforhashtable            */

  foutp.open(f0name);

  i1hat=TalINT1(MHASH);
  c1dig=TalCHR1(MBASE);

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

  nprec=4;

  for(i=0;i<MHASH;i++)
  { i1hat[i]=-1;
  }
  for(i=0;i<MBASE;i++)
  { j=(INT)c1dig[i];
    if(j<0)j=j+MHASH/2;
    i1hat[j]=i;
  }
}

Y_NANOWriteGraph::~Y_NANOWriteGraph()
{ FREE(i1hat);
  FREE(c1dig);
  foutp.close();
}

void Y_NANOWriteGraph::INTtoCHR(CHR *c1word, INT inum, INT *i0)
{ INT i;                /* positionindexinsideC1WORD            */
  INT idig;             /*                                      */

  i=*i0;
  if(inum<0)
  { c1word[i]='-';
    i=i+1;
    inum=-inum;
  }
  else if(inum==0)
  { c1word[i]='+';
    i=i+1;
    c1word[i]=c1dig[0];
    i=i+1;
  }
  else
  { c1word[i]='+';
    i=i+1;
  }
  while(inum>0)
  { idig=inum%MBASE;
    inum=inum/MBASE;
    c1word[i]=c1dig[idig];
    i=i+1;
  }
  *i0=i;
}

void Y_NANOWriteGraph::DBLtoCHR(CHR *c1word, DBL dnum, INT *i0)
{ INT i;                /*                                      */
  INT idig;             /*                                      */
  INT iexp;             /*                                      */
  INT j;                /*                                      */
  CHR c1dumm[MCHPO];    /*                                      */
  
  i=*i0;
  if(dnum<D0)
  { c1word[i]=',';
    i=i+1;
    dnum=-dnum;
  }
  else
  { c1word[i]='.';
    i=i+1;
  }
  
  iexp=0;
  while((dnum<D1)&&(iexp>-(MBASE/2-nprec)))
  { dnum=dnum*(DBL)MBASE;
    iexp=iexp-1;
  }
  while((dnum>=MBASE)&&(iexp<(MBASE/2-nprec)))
  { dnum=dnum/(DBL)MBASE;
    iexp=iexp+1;
  }
  for(j=0;j<nprec;j++)
  { idig=(INT)dnum;
    c1dumm[j]=c1dig[idig];
    dnum=dnum-idig;
    dnum=dnum*(DBL)MBASE;
  }
  iexp=iexp-(nprec-1);

  c1word[i]=c1dig[iexp+MBASE/2];
  i=i+1;
  for(j=(nprec-1);j>=0;j--)
  { c1word[i]=c1dumm[j];
    i=i+1;
  }
  *i0=i;
}

void Y_NANOWriteGraph::WriteGrph(Y_NANOgraph *v0graph)
{ INT i;                /*                                      */
  INT ipar;             /*                                      */
  INT j;                /*                                      */
  INT nipar;            /*                                      */
  INT ndpar;            /*                                      */
  DBL dpar;             /*                                      */

  i=0;
  c1word[i]='<';
  i=i+1;
  nipar=v0graph->nipar;
  for(j=0;j<nipar;j++)  /* Write INT parameters                 */
  { ipar=v0graph->i1par[j];
    INTtoCHR(c1word,ipar,&i);
  }
  ndpar=v0graph->ndpar;
  for(j=0;j<ndpar;j++)  /* Write DBL parameters                 */
  { dpar=v0graph->d1par[j];
    DBLtoCHR(c1word,dpar,&i);
  }
  c1word[i]='>';
  i=i+1;
  c1word[i]='\0';
  i=i+1;

  for(j=0;j<i;j++)
  { foutp.write(&c1word[j],1);
  }
  foutp.write("\n",1);
}

