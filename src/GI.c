/////////////////////////////////////////////////////////////////////////
//  Gauss Integrals - A tool for describing protein backbone geometry  //
//                                                                     //
//  Authors: Peter R{\o}gen & Boris Fain                               //
/////////////////////////////////////////////////////////////////////////

/* 
Copyright (c)   The researchers and professors of the Technical 
University of Denmark. All Rights Reserved
 
Permission to use, copy, modify and distribute any part of this GI 
software for educational, research and non-profit purposes, without fee, 
and without a written agreement is hereby granted, provided that the 
above copyright notice, this paragraph and the following four paragraphs 
appear in all copies.
 
Those desiring to incorporate this Software into commercial products or 
use for commercial purposes should contact the Technology Transfer 
Office, The Secretariat of Research and Innovation, Technical University 
of Denmark, Anker Engelundsvej, Building 101 A, 2800 Kgs. Lyngby, 
DENMARK. Fax (+45) 45888040.
 
IN NO EVENT SHALL THE TECHNICAL UNIVERSITY OF DENMARK BE LIABLE TO ANY 
PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL 
DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE,
EVEN IF THE TECHNICAL UNIVERSITY OF DENMARK  HAS BEEN ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.
 
THE  SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE TECHNICAL
UNIVERSITY OF DENMARK HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, 
UPDATES, ENHANCEMENTS OR MODIFICATIONS. THE TECHNICAL UNIVERSITY OF 
DENMARK MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND,
EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR 
THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR
OTHER RIGHTS.

CITATIONS
For publication of results, please cite:
P. R{\o}gen & B. Fain, Automatic classification of protein structure 
by using Gauss integrals, PNAS, 100(1), 119-124, 2003.

*/

#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
#define MAXLENGTH 1300
#define dipeptide 0  /*if 1 dipeptide ferquences of chains are printed*/
#define peptidesequence 0 /*if 1 peptide sequence of chains are printed*/
#define CylinderTransform 1 /*The Cylindertransform is normalized such 
                              that chain length and each structural measure 
                              has derivation one on a set of representatives 
                              of the homology (H) classes of CATH 2.4. 
   
                              Under this transform the Scaled Gauss Metric, 
                              SGM,  introduced in
   
                              P. R{\o}gen & B. Fain, Automatic classification 
                              of protein structure by using Gauss integrals, 
                              PNAS, 100(1), 119-124, 2003.
   
                              is given by the usual scalar product in R^30 of the 
                              last 30 columns of each output line.

 
                              If Cylindertransform is 0, GI writes the raw Gauss
                              integrals as output.
                           */

FILE *farrayo, *ferroro;

double twoPI=2.0*M_PI;
double polyl[MAXLENGTH][3];                     /*C-alpha coordinates*/
double unitvectors[MAXLENGTH][MAXLENGTH][3];
double omega[MAXLENGTH][MAXLENGTH];
double absomega[MAXLENGTH][MAXLENGTH];
double partsum[MAXLENGTH][MAXLENGTH];
double abspartsum[MAXLENGTH][MAXLENGTH];
double len;
int nres, n_proteins=0;
char GlobalOneLetterId[MAXLENGTH];
int SequinceInteger[MAXLENGTH];



double length(double vec[3]);
double dot(double vec1[3],double vec2[3]);
void geovec(double vec2_minus_vec1[3], double vec1[3], double vec2[3]);
void unit(double uvec[3],double vec[3]);
void cross(double crossvec[3],double vec1[3],double vec2[3]);
double spangle(double vec1[3], double vec2[3], double vec3[3]);
void createunitvectors(void);
void createomega(void);
void createpartsums(void);
double mixedsums(int a, int b, int c, int d);
double absmixedsums(int a, int b, int c, int d);
double int12(void);
double inta12(void);

double int12_34(void);
double inta12_34(void);
double int12_a34(void);
double inta12_a34(void);

double int13_24(void);
double inta13_24(void);
double int13_a24(void);
double inta13_a24(void);

double int14_23(void);
double inta14_23(void);
double int14_a23(void);
double inta14_a23(void);

double int12_34_56(void);
double int12_35_46(void);
double int12_36_45(void);
double int13_24_56(void);
double int13_25_46(void);
double int13_26_45(void);
double int14_23_56(void);
double int14_25_36(void);
double int14_26_35(void);
double int15_23_46(void);
double int15_24_36(void);
double int15_26_34(void);
double int15_26_34(void);
double int16_23_45(void);
double int16_24_35(void);
double int16_25_34(void);



void GaussIntegrals(void);
void print_dipeptide_frequeces(void);
void print_peptide_sequence(void);
int Three2IntegerProteinId(char ThreeLeterAmno[3]);


void read_pdb_file(char *dir, char *fname);

void dirwalk(char *dir)
{
  char name[1024];
  struct dirent *dp;
  DIR *dfd;

  if ((dfd=opendir(dir))==NULL) 
    {
      printf("Sorry - could not open directory <%s>\n",dir);
      return; 
    }
  while ((dp=readdir(dfd))!=NULL)
  {
    if (strcmp(dp->d_name,".")==0 || strcmp(dp->d_name,"..")==0) continue;
    sprintf(name,"%s/%s",dir,dp->d_name);
    if (dp->d_name[strlen(dp->d_name)-4]=='.' && dp->d_name[strlen(dp->d_name)-3]=='p' && dp->d_name[strlen(dp->d_name)-2]=='d'
     && dp->d_name[strlen(dp->d_name)-1]=='b')
    {
      printf("After %i chains, it's %s with ",n_proteins,dp->d_name);
      read_pdb_file(dir,dp->d_name);
    }
  }
  closedir(dfd);
}

char line[1000];

int getatom(FILE *fi)
{
  char *cp;
  while (cp=fgets(line,999,fi),strstr(line,"ATOM")!=line && cp!=NULL);
  if (cp==NULL) return 0;
  return 1;
}

void read_pdb_file(char *dir, char *fname)
{
  FILE *fi;
  double x, y, z;
  int nr, no=-1500, link_count=0, start_chain,i;
  char chain_name, new_chain_name;
  char ThreeLeterAmno[3];
  char dir_fname[1024], *cp;

  strcpy(dir_fname,dir);
  strcat(dir_fname,fname); 

  fi=fopen(dir_fname,"r");
  chain_name='a';
loop:
  if (getatom(fi)==0) 
    { 
      if(link_count<MAXLENGTH && link_count>9)
	{
	  if(dipeptide==0){fprintf(farrayo,"%s %c %i %i",fname,chain_name,link_count,no-start_chain+1-link_count);}
	  if(dipeptide==1){fprintf(farrayo,"%s ",fname);}
	  nres=link_count;
	  printf("%i residues\n",link_count);
	  GaussIntegrals();
	  fclose(fi);
	  return;
	}
      fclose(fi);  
      return;
    }
  sscanf(line+22,"%4i",&nr);
  sscanf(line+21,"%1c",&new_chain_name);
  if(new_chain_name==' ') new_chain_name='0';
  
  if (nr!=no && line[13]=='C' && line[14]=='A')
    {
      if(new_chain_name!=chain_name && chain_name!='a')
	{
	  if(link_count<MAXLENGTH && link_count>9)
	    {
	      if(dipeptide==0){fprintf(farrayo,"%s %c %i %i",fname,chain_name,link_count,no-start_chain+1-link_count);}
	      if(dipeptide==1){fprintf(farrayo,"%s ",fname);}
	      nres=link_count;  
	      printf("%i residues ",link_count);
	      GaussIntegrals();
	    }  
	  link_count=0;
	}
      link_count++;
      if (link_count>MAXLENGTH-1) {
	printf(" >MAXLENGTH=%i residues\n",MAXLENGTH);
	fprintf(ferroro, "WARNING!  Chain %c of <%s> is longer than MAXLENGTH=%i. \n          No further attempt to calculate Gauss Integrals for chains of\n          <%s> has been made. Consider redefine MAXLENGTH in <GI.c>.\n", chain_name, fname, MAXLENGTH, fname);
	fclose(fi); return;
      }
      no=nr;
      if(link_count==1) start_chain=no;
      sscanf(line+30,"%lf %lf %lf",&x,&y,&z);
      ThreeLeterAmno[0]=line[17];
      ThreeLeterAmno[1]=line[18];
      ThreeLeterAmno[2]=line[19];
      SequinceInteger[link_count]=Three2IntegerProteinId(ThreeLeterAmno);
      polyl[link_count][0]=x;
      polyl[link_count][1]=y;
      polyl[link_count][2]=z;
      chain_name=new_chain_name;
      
    }
  goto loop;
}


int Three2IntegerProteinId(char ThreeLeterAmno[3])
{
  if( strncmp(ThreeLeterAmno,"ALA",3)==0) return 5;
  if( strncmp(ThreeLeterAmno,"GLU",3)==0) return 4;
  if( strncmp(ThreeLeterAmno,"GLN",3)==0) return 7;
  if( strncmp(ThreeLeterAmno,"ASP",3)==0) return 3;
  if( strncmp(ThreeLeterAmno,"ASN",3)==0) return 6;
  if( strncmp(ThreeLeterAmno,"LEU",3)==0) return 17;
  if( strncmp(ThreeLeterAmno,"GLY",3)==0) return 1;
  if( strncmp(ThreeLeterAmno,"LYS",3)==0) return 10;
  if( strncmp(ThreeLeterAmno,"SER",3)==0) return 8;
  if( strncmp(ThreeLeterAmno,"VAL",3)==0) return 13;
  if( strncmp(ThreeLeterAmno,"ARG",3)==0) return 11;
  if( strncmp(ThreeLeterAmno,"THR",3)==0) return 9;
  if( strncmp(ThreeLeterAmno,"PRO",3)==0) return 2;
  if( strncmp(ThreeLeterAmno,"ILE",3)==0) return 14;
  if( strncmp(ThreeLeterAmno,"MET",3)==0) return 15;
  if( strncmp(ThreeLeterAmno,"PHE",3)==0) return 18;
  if( strncmp(ThreeLeterAmno,"TYR",3)==0) return 19;
  if( strncmp(ThreeLeterAmno,"CYS",3)==0) return 16;
  if( strncmp(ThreeLeterAmno,"TRP",3)==0) return 20;
  if( strncmp(ThreeLeterAmno,"HIS",3)==0) return 12;
  if( strncmp(ThreeLeterAmno,"ASX",3)==0) return 3;
  if( strncmp(ThreeLeterAmno,"Glx",3)==0) return 1;
  if( strncmp(ThreeLeterAmno,"PCA",3)==0) return 8;
  fprintf(ferroro, "WARNING! Do not know the amino acid residue %c%c%c \n",ThreeLeterAmno[0],ThreeLeterAmno[1],ThreeLeterAmno[2]);
  return 21;
}





/* BASIC VECTOR CALCULOS  */


/* length: euklidean lenght in R^3*/

double length(double vec[3])
{
  return sqrt(dot(vec,vec));
}


/* dot: dot product in R^3
 */
double dot(double vec1[3],double vec2[3])
{
  return vec1[0]*vec2[0]+ vec1[1]*vec2[1]+ vec1[2]*vec2[2];
}


/* unit: gives the unit vector v/|v|, uvec, in R^3
 */
void unit(double uvec[3], double vec[3])
{
  int i;
  double len;
  len=length(vec);
  for (i = 0; i<3;i++)
     uvec[i]=vec[i]/len;
}

/* geovec: Gives the vector from vec1 to vec2
 */
void geovec(double vec2_minus_vec1[3], double vec1[3], double vec2[3])
{
  int i;
  for (i=0; i<3;i++)
  vec2_minus_vec1[i]=vec2[i]-vec1[i];
}

/* cross: gives the cross product, crossvec, of vec1 and vec2 in R^3
 */
void cross(double crossvec[3],double vec1[3],double vec2[3])
{
  crossvec[0]=vec1[1]*vec2[2]-vec2[1]*vec1[2];
  crossvec[1]=vec1[2]*vec2[0]-vec2[2]*vec1[0];
  crossvec[2]=vec1[0]*vec2[1]-vec2[0]*vec1[1];
}


/* cerateunitvectors creates the unit vectors from polyl[i] to polyl[j].
 */
void createunitvectors(void)
{
  int i,j;
  extern double polyl[MAXLENGTH][3];
  extern double unitvectors[MAXLENGTH][MAXLENGTH][3];
  extern int nres;
  double varvec[3];
  for (i=1;i<nres;i++)
    {
      for (j=i+1;j<nres+1;j++)
	{
	  geovec(varvec,polyl[i],polyl[j]);
	  unit(unitvectors[i][j],varvec);
	}
    }
}


/* SPHERICAL GEOMETRY  */

/* spangle: For  vec1, vec2, and, vec3 on the unit 2-sphere 
 *          spangle gives the exterior angle 
 *          spangle(vec1,vec2,vec3) between the geodesic segments 
 *          from vec1 to vec2 and from vec2 to vec3. 
 */

double spangle(double vec1[3], double vec2[3], double vec3[3])
{
  double vv[3];
  double ab,bc;
  cross(vv,vec1,vec2);
  ab=dot(vec1,vec2);
  bc=dot(vec2,vec3);
  return atan2(dot(vv,vec3),ab*bc- dot(vec1,vec3) );
}



/* cerateomega creates half of the mesure of area on the unit 2-sphere
 *             from which a crossing is seen between 
 *             the line segment polyl[i] to  polyl[i+1] and 
 *             the line segment polyl[j] to polyl[j+1].
 *             Cf.  P. R{\o}gen & H.G. Bohr, A new family of global protein
 *             shape descriptors,  Mathematical Biosciences Volume 182,
 *             Issue 2, April 2003, Pages 167-181.
 */


void createomega(void)
{
  int i,j;
  extern double omega[MAXLENGTH][MAXLENGTH];
  extern double absomega[MAXLENGTH][MAXLENGTH];
  extern double unitvectors[MAXLENGTH][MAXLENGTH][3];
  extern int nres;
  double anglesum, cross1[3], cross2[3], cross3[3];
  for (i=1;i<nres;i++)
    {
      omega[i][i+1]=0;
      absomega[i][i+1]=0;
    }
  for (i=1;i<nres-1;i++)
    {
      for (j=i+2;j<nres;j++)
	{
	  cross(cross1,unitvectors[i][j],unitvectors[i][j+1]);
	  cross(cross2,unitvectors[i+1][j],unitvectors[i+1][j+1]);
	  cross(cross3,cross1,cross2);
	  if (dot(cross3,cross3)<pow(10,-16)) 
	    {
	      omega[i][j]=0; 
	      absomega[i][j]=0;
	    }
	  else
	    {
	      anglesum=-spangle(unitvectors[i][j],unitvectors[i][j+1],unitvectors[i+1][j+1])-spangle(unitvectors[i][j+1],unitvectors[i+1][j+1],unitvectors[i+1][j])-spangle(unitvectors[i+1][j+1],unitvectors[i+1][j],unitvectors[i][j])-spangle(unitvectors[i+1][j],unitvectors[i][j],unitvectors[i][j+1]);
	      if ( anglesum>=0 )
		omega[i][j]=twoPI-anglesum;
	      else
		omega[i][j]=-twoPI-anglesum;
	      absomega[i][j]=fabs(omega[i][j]);
	    }
	}
    }
}




/* createpartsums creats partsum[i][j]=writhe of the segment
 *                from polyl[i] to polyl[j+1]
 */
void createpartsums(void)
{
  extern double omega[MAXLENGTH][MAXLENGTH];
  extern double absomega[MAXLENGTH][MAXLENGTH];
  extern int nres;
  extern double partsum[MAXLENGTH][MAXLENGTH]; 
  extern double abspartsum[MAXLENGTH][MAXLENGTH]; 
  int a, k;
  for (k=-1;k<2;k++)
    for (a=1;a<nres -k;a++)
      partsum[a][a+k]=0;
  for (a=1;a<nres-2;a++)
    partsum[a][a+2]=omega[a][a+2];
  for (k=3;k<nres-1;k++)
    for (a=1; a<nres-k;a++)
      partsum[a][a+k]=partsum[a][a+k-1]+partsum[a+1][a+k]+omega[a][a+k]-partsum[a+1][a+k-1];
  for (k=-1;k<2;k++)
    for (a=1;a<nres -k;a++)
      abspartsum[a][a+k]=0;
  for (a=1;a<nres-2;a++)
    abspartsum[a][a+2]=absomega[a][a+2];
  for (k=3;k<nres-1;k++)
    for (a=1; a<nres-k;a++)
      abspartsum[a][a+k]=abspartsum[a][a+k-1]+abspartsum[a+1][a+k]+absomega[a][a+k]-abspartsum[a+1][a+k-1];
}

/* Calculation of "the mixed sums" equivalent the the 
 * writhe contribution from the inteval [a;b]
 * to the interval [c;d]. 
 * NB: [a;a] is the a'th line segment 
 * - the one from polyl[a] to polyl[a+1].
 */

double mixedsum(int a,int b,int c,int d)
{
  double tempm, tempp;
  tempm=partsum[a][c-1]+partsum[b+1][d];
  tempp=partsum[a][d]+partsum[b+1][c-1];
  return tempp-tempm;
}

double absmixedsum(int a,int b,int c,int d)
{
  double tempm, tempp;
  tempm=abspartsum[a][c-1]+abspartsum[b+1][d];
  tempp=abspartsum[a][d]+abspartsum[b+1][c-1];
  return tempp-tempm;
}




/*            CALCULATION OF THE GAUSS INTEGRALS
 */

/******** FIRST ORDER GAUSS INTEGRALS **********************/

double int12(void)
{
  return partsum[1][nres-1]/(twoPI);
}

double inta12(void)
{
  return abspartsum[1][nres-1]/(twoPI);
}

/********* SECOND ORDER GAUSS INTEGRALS ***********************/


double int12_34(void)
{
  int b;
  double temp;
  temp=0;
  for(b=4;b<nres-2; b++)
    temp=temp+partsum[1][b-1]*mixedsum(b,b,b+2,nres-1);
  return temp/(twoPI*twoPI);
}

double inta12_34(void)
{
  int b;
  double temp;
  temp=0;
  for(b=4;b<nres-2; b++)
    temp=temp+abspartsum[1][b-1]*mixedsum(b,b,b+2,nres-1);
  return temp/(twoPI*twoPI);
}

double int12_a34(void)
{
  int b;
  double temp;
  temp=0;
  for(b=4;b<nres-2; b++)
    temp=temp+partsum[1][b-1]*absmixedsum(b,b,b+2,nres-1);
  return temp/(twoPI*twoPI);
}

double inta12_a34(void)
{
  int b;
  double temp;
  temp=0;
  for(b=4;b<nres-2; b++)
    temp=temp+abspartsum[1][b-1]*absmixedsum(b,b,b+2,nres-1);
  return temp/(twoPI*twoPI);
}

/*************************************************************/


double int13_24(void)
{
  int a, b;
  double temp;
  temp=0;
  for(a=2;a<nres-2; a++)
    for(b=a+2;b<nres;b++)
      temp=temp+mixedsum(1,a-1,a+1,b-1)*omega[a][b];
  return temp/(twoPI*twoPI);
}

double inta13_24(void)
{
  int a, b;
  double temp;
  temp=0;
  for(a=2;a<nres-2; a++)
    for(b=a+2;b<nres;b++)
      temp=temp+absmixedsum(1,a-1,a+1,b-1)*omega[a][b];
  return temp/(twoPI*twoPI);
}

double int13_a24(void)
{
  int a, b;
  double temp;
  temp=0;
  for(a=2;a<nres-2; a++)
    for(b=a+2;b<nres;b++)
      temp=temp+mixedsum(1,a-1,a+1,b-1)*absomega[a][b];
  return temp/(twoPI*twoPI);
}

double inta13_a24(void)
{
  int a, b;
  double temp;
  temp=0;
  for(a=2;a<nres-2; a++)
    for(b=a+2;b<nres;b++)
      temp=temp+absmixedsum(1,a-1,a+1,b-1)*absomega[a][b];
  return temp/(twoPI*twoPI);
}

/***************************************************************/

double int14_23(void)
{
  int a, b;
  double temp;
  temp=0;
  for(a=1;a<nres-4; a++)
    for(b=a+4;b<nres;b++)
      temp=temp+partsum[a+1][b-1]*omega[a][b];
  return temp/(twoPI*twoPI);
}

double inta14_23(void)
{
  int a, b;
  double temp;
  temp=0;
  for(a=1;a<nres-4; a++)
    for(b=a+4;b<nres;b++)
      temp=temp+partsum[a+1][b-1]*absomega[a][b];
  return temp/(twoPI*twoPI);
}

double int14_a23(void)
{
  int a, b;
  double temp;
  temp=0;
  for(a=1;a<nres-4; a++)
    for(b=a+4;b<nres;b++)
      temp=temp+abspartsum[a+1][b-1]*omega[a][b];
  return temp/(twoPI*twoPI);
}

double inta14_a23(void)
{
  int a, b;
  double temp;
  temp=0;
  for(a=1;a<nres-4; a++)
    for(b=a+4;b<nres;b++)
      temp=temp+abspartsum[a+1][b-1]*absomega[a][b];
  return temp/(twoPI*twoPI);
}

/**************THIRD ORDER GAUSS INTEGRALS********************/


double int12_34_56(void)
{
  int a, b;
  double temp1, temp2;
  temp1=0;
  for(a=4;a<nres-5;a++)
    {
      temp2=0;
      for(b=a+3;b<nres-2;b++)
	temp2=temp2+mixedsum(a,a,a+2,b-1)*mixedsum(b,b,b+2,nres-1);
      temp1=temp1+temp2*partsum[1][a-1];
    }
  return temp1/(twoPI*twoPI*twoPI);
}



double int12_35_46(void)
{
  int a, b;
  double temp1, temp2;
  temp1=0;
  for(a=4;a<nres-3;a++)
    {
      temp2=0;
      for(b=a+2;b<nres-1;b++)
	temp2=temp2+omega[a][b]*mixedsum(a+1,b-1,b+1,nres-1);
      temp1=temp1+temp2*partsum[1][a-1];
    }
  return temp1/(twoPI*twoPI*twoPI);
}


double int12_36_45(void)
{
  int a, b;
  double temp1, temp2;
  temp1=0;
  for(a=4;a<nres-4;a++)
    {
      temp2=0;
      for(b=a+4; b<nres; b++)
	temp2=temp2+omega[a][b]*partsum[a+1][b-1];
      temp1=temp1+temp2*partsum[1][a-1];
    }
  return temp1/(twoPI*twoPI*twoPI);
}


double int13_24_56(void)
{
  int a, b;
  double temp1, temp2;
  temp1=0;
  for(a=4;a<nres-3;a++)
    {
      temp2=0;
      for(b=2;b<a-1;b++)
	temp2=temp2+omega[b][a]*mixedsum(1,b-1,b+1,a-1);
      temp1=temp1+temp2*partsum[a+1][nres-1];
    }
  return temp1/(twoPI*twoPI*twoPI);
}


double int13_25_46(void)
{
  int a, b, c;
  double temp1, temp2;
  temp1=0;
  for(a=2;a<nres-4;a++)
    {
      for(b=a+2;b<nres-2;b++)
	{
	  temp2=0;
	  for(c=b+1;c<nres-1;c++)
	    temp2=temp2+omega[a][c]*mixedsum(b,b,c+1,nres-1);
	  temp1=temp1+temp2*mixedsum(1,a-1,a+1,b-1);
	}
    }
  return temp1/(twoPI*twoPI*twoPI);
}


double int13_26_45(void)
{
  int a, b, c;
  double temp1, temp2;
  temp1=0;
  for(a=2;a<nres-5;a++)
    {
      for(b=a+5; b<nres; b++)
	{
	  temp2=0;
	  for(c=a+1;c<b-3;c++)
	    temp2=temp2+partsum[c+1][b-1]*mixedsum(1,a-1,c,c);
	  temp1=temp1+temp2*omega[a][b];
	}
    }
  return temp1/(twoPI*twoPI*twoPI);
}

double int14_23_56(void)
{
  int a, b;
  double temp1, temp2;
  temp1=0;
  for(a=5;a<nres-3;a++)
    {
      temp2=0;
      for(b=1; b<a-3; b++)
	temp2=temp2+omega[b][a]*partsum[b+1][a-1];
      temp1=temp1+temp2*partsum[a+1][nres-1];
    }
  return temp1/(twoPI*twoPI*twoPI);
}


double int14_25_36(void)
{
  int a, b, c;
  double temp1, temp2;
  temp1=0;
  for(a=2;a<nres-4;a++)
    {
      for(b=a+3; b<nres-1; b++)
	{
	  temp2=0;
	  for(c=a+1; c<b-1; c++)
	    temp2=temp2+mixedsum(1,a-1,c+1,b-1)*mixedsum(c,c,b+1,nres-1);
	  temp1=temp1+temp2*omega[a][b];
	}
    }
  return temp1/(twoPI*twoPI*twoPI);
}


double int14_26_35(void)
{
  int a, b, c;
  double temp1, temp2;
  temp1=0;
  for(a=5;a<nres-1;a++)
    {
      for(b=3; b<a-1; b++)
	{
	  temp2=0;
	  for(c=1; c<b-1; c++)
	    temp2=temp2+mixedsum(c+1,b-1,a+1,nres-1)*mixedsum(c,c,b+1,a-1);
	  temp1=temp1+temp2*omega[b][a];
	}
    }
  return temp1/(twoPI*twoPI*twoPI);
}

double int15_23_46(void)
{
  int a, b, c;
  double temp1, temp2;
  temp1=0;
  for(a=1;a<nres-6;a++)
    {
      for(b=a+5; b<nres-1; b++)
	{
	  temp2=0;
	  for(c=a+4; c<b; c++)
	    temp2=temp2+partsum[a+1][c-1]*mixedsum(c,c,b+1,nres-1);
	  temp1=temp1+temp2*omega[a][b];
	}
    }
  return temp1/(twoPI*twoPI*twoPI);
}



double int15_24_36(void)
{
  int a, b, c;
  double temp1, temp2;
  temp1=0;
  for(a=2;a<nres-4;a++)
    {
      for(b=a+2; b<nres-2; b++)
	{
	  temp2=0;
	  for(c=b+2; c<nres; c++)
	    temp2=temp2+mixedsum(1,a-1,b+1,c-1)*mixedsum(a+1,b-1,c,c);
	  temp1=temp1+temp2*omega[a][b];
	}
    }
  return temp1/(twoPI*twoPI*twoPI);
}

double int15_26_34(void)
{
  int a, b, c;
  double temp1, temp2;
  temp1=0;
  for(a=2; a<nres-5; a++)
    {
      for(b=a+5; b<nres; b++)
	{
	  temp2=0;
	  for(c=a+4; c<b; c++)
	    temp2=temp2+partsum[a+1][c-1]*mixedsum(1,a-1,c,c);
	  temp1=temp1+temp2*omega[a][b];
	}
    }
  return temp1/(twoPI*twoPI*twoPI);
}

double int16_23_45(void)
{
  int a, b, c;
  double temp1, temp2;
  temp1=0;
  for(a=1; a<nres-7; a++)
    {
      for(b=a+7; b<nres; b++)
	{
	  temp2=0;
	  for(c=a+4; c<b-2; c++)
	    temp2=temp2+partsum[a+1][c-1]*mixedsum(c,c,c+2,b-1);
	  temp1=temp1+temp2*omega[a][b];
	}
    }
  return temp1/(twoPI*twoPI*twoPI);
}


double int16_24_35(void)
{
  int a, b, c;
  double temp1, temp2;
  temp1=0;
  for(a=3; a<nres-3; a++)
    {
      for(b=a+2; b<nres-1; b++)
	{
	  temp2=0;
	  for(c=1; c<a-1; c++)
	    temp2+=mixedsum(c,c,b+1,nres-1)*mixedsum(c+1,a-1,a+1,b-1);
	  temp1+=temp2*omega[a][b];
	}
    }
  return temp1/(twoPI*twoPI*twoPI);
}


double int16_25_34(void)
{
  int a, b;
  double temp1;
  temp1=0;
  for(a=2; a<nres-5; a++)
    for(b=a+4; b<nres-1; b++)
      temp1=temp1+mixedsum(1,a-1,b+1,nres-1)*omega[a][b]*partsum[a+1][b-1];
  return temp1/(twoPI*twoPI*twoPI);
}

/**********************PRINTING*******************/

void GaussIntegrals(void)
{
  extern int nres;
  extern double polyl[MAXLENGTH][3];
  extern double partsum[MAXLENGTH][MAXLENGTH];

  createunitvectors();
  createomega();
  createpartsums();

  if(dipeptide==1){print_dipeptide_frequeces();}
  if(peptidesequence==1){print_peptide_sequence();}

  if(CylinderTransform==0){
    fprintf(farrayo," %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e\n",int12(),inta12(),int12_34(),inta12_34(),int12_a34(),inta12_a34(),int13_24(),inta13_24(),int13_a24(),inta13_a24(),int14_23(),inta14_23(),int14_a23(),inta14_a23(),int12_34_56(),int12_35_46(),int12_36_45(),int13_24_56(),int13_25_46(),int13_26_45(),int14_23_56(),int14_25_36(),int14_26_35(),int15_23_46(),int15_24_36(),int15_26_34(),int16_23_45(),int16_24_35(),int16_25_34());}

    if(CylinderTransform==1){
    fprintf(farrayo," %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e\n",nres/70.2135/1.339875,int12()/nres/0.0401388/1.030110,inta12()/nres/0.0730344/1.222925,int12_34()/(nres*nres)/0.0014755/1.200300,inta12_34()/(nres*nres)/0.00533629/1.050380,int12_a34()/(nres*nres)/0.00513504/1.058926,inta12_a34()/(nres*nres)/0.0157794/1.174822,int13_24()/(nres*nres)/5.7905e-05/1.462006,inta13_24()/(nres*nres)/0.000460471/1.044327,int13_a24()/(nres*nres)/0.000485036/1.095845,inta13_a24()/(nres*nres)/0.00418937/1.282707,int14_23()/(nres*nres)/0.000188964/1.097088,inta14_23()/(nres*nres)/0.00269016/0.9816835,int14_a23()/(nres*nres)/0.000892068/1.153978,inta14_a23()/(nres*nres)/0.00758828/1.241605,int12_34_56()/(nres*nres*nres)/3.7779e-05/1.394955,int12_35_46()/(nres*nres*nres)/2.17081e-06/1.615541,int12_36_45()/(nres*nres*nres)/2.98437e-06/1.098159,int13_24_56()/(nres*nres*nres)/2.21604e-06/1.592069,int13_25_46()/(nres*nres*nres)/5.8425e-08/1.692098,int13_26_45()/(nres*nres*nres)/3.64209e-07/1.070658,int14_23_56()/(nres*nres*nres)/3.60442e-06/1.049276,int14_25_36()/(nres*nres*nres)/6.07681e-08/1.313312,int14_26_35()/(nres*nres*nres)/6.60069e-08/1.148898,int15_23_46()/(nres*nres*nres)/2.95267e-07/1.402297,int15_24_36()/(nres*nres*nres)/5.84497e-08/1.369921,int15_26_34()/(nres*nres*nres)/2.98807e-07/1.291011,int16_23_45()/(nres*nres*nres)/3.66585e-06/1.254290,int16_24_35()/(nres*nres*nres)/2.25853e-07/1.359014,int16_25_34()/(nres*nres*nres)/4.20434e-07/1.271251);}


  n_proteins++;

}
void print_peptide_sequence(void)
{
  int i;
  fprintf(farrayo," ");
  for(i=1; i<nres+1; i++){
    fprintf(farrayo,"%c",SequinceInteger[i]+'a');
  }
}

void print_dipeptide_frequeces(void)
{
  int i,j;
  int dipepfreq[21][21];
  
  for(i=0; i<21; i++){
    for(j=0; j<21; j++){
      dipepfreq[i][j]=0;
    }
  }

  for(i=1; i<nres;i++) {dipepfreq[SequinceInteger[i]-1][SequinceInteger[i+1]-1]++;}
  fprintf(farrayo,"%i",dipepfreq[0][0]);
  for(j=1;j<20;j++){fprintf(farrayo," %i",dipepfreq[0][j]);}
  for(i=1; i<20; i++){
    for(j=0; j<20; j++){
      fprintf(farrayo," %i",dipepfreq[i][j]);
    }
  }
}



int main(int argc, char *argv[])
{
  if (argc<4)
    {
      printf("\nGauss Integrals: Please input directory (with trailing \"/\") followed by \n                 output and error file names as arguments.\n\n");
      return 1;
    }


  ferroro=fopen(argv[3],"w");
  if (ferroro==(FILE *)NULL)
    {
      printf("Sorry - could not open file <%s> \n",argv[3]);
      return 1;
    }
  farrayo=fopen(argv[2],"w");
  if (farrayo==(FILE *)NULL)
    {
      printf("Sorry - could not open file <%s>\n",argv[2]);      
      return 1;
    }
  
  dirwalk(argv[1]);

  printf("There has been calculated 29 Gauss-integrals for %i protein chains.\n", n_proteins);
  fclose(farrayo);
  fclose(ferroro);
  return 0;
}
