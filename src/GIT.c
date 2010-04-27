//////////////////////////////////////////////////////////////////////////////
//  Gauss Integrals Tuned- A tool for describing protein backbone geometry  //
//                                                                          //
//  Author: Peter R{\o}gen                                                  //
//////////////////////////////////////////////////////////////////////////////

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

P. R{\o}gen, Evaluating protein structure descriptors and tuning Gauss integral based descriptors: Journal of Physics Condensed Matter, vol: 17, pages: 1523-1538, 2005.

*/



#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>

#define MAXLENGTH 874
#define dipeptide 0
#define CylinderTransform 0
#define SubtractAverageProtein 1
#define SmoothenBackBone 1 /* 1 for smoothing of the carbon alpha curve  */
			   /* 0 for the carbon alpha curve               */


       
 FILE *farrayo, *ferroro, *favgauss;
double twoPI=2.0*M_PI;
double epsi=0.01;
double polyl[MAXLENGTH][3];
double testpolyl[14][3]={{10,10,10},{0,0,0},{1,0,0},{1,1.1,0},{1,2,0},{0,2,0},{0,2,1},{1-0.00009,2,1},{1,1,1},{1-0.00009,1,0.5},{1-0.00009,1,0.3},{1-0.00009,1,0.1},{1-0.00009,1,-0.5},{1-0.00009,1,-1}};

double unitvectors[MAXLENGTH][MAXLENGTH][3];
double omega[MAXLENGTH][MAXLENGTH];
double absomega[MAXLENGTH][MAXLENGTH];
double partsum[MAXLENGTH][MAXLENGTH];
double abspartsum[MAXLENGTH][MAXLENGTH];
double refGauss[MAXLENGTH+1000][29];

double len;
int nres, n_proteins=0;

char GlobalOneLetterId[MAXLENGTH];
int SequinceInteger[MAXLENGTH];

double rnd(void);
double length(double vec[3]);
double dot(double vec1[3],double vec2[3]);
void geovec(double vec2_minus_vec1[3], double vec1[3], double vec2[3]);
void unit(double uvec[3],double vec[3]);
void cross(double crossvec[3],double vec1[3],double vec2[3]);
double determinant(double vec1[3],double vec2[3],double vec3[3]);
double spangle(double vec1[3], double vec2[3], double vec3[3]);
int intersectiontrianglelinesegment(double a[3],double b[3],double c[3],double p1[3],double p2[3]);
void PerturbBackBone(double ampitude);
void OrtogonalPerturbBackBone(double ampitude);
void SmoothenTheBackbone(void);
double presplinedeformationsmoothorg(double measurevalue, int measurenumber);
double presplinedeformationsmoothnew(double measurevalue, int measurenumber);
double splinedeformationsmoothorg(double measurevalue, int measurenumber);
double splinedeformationsmoothnew(double measurevalue, int measurenumber);

double presplinedeformationCalphaorg(double measurevalue, int measurenumber);
double presplinedeformationCalphanew(double measurevalue, int measurenumber);
double splinedeformationCalphaorg(double measurevalue, int measurenumber);
double splinedeformationCalphanew(double measurevalue, int measurenumber);
double splinefun(double t);


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

double distance(void);

void GaussIntegrals(void);
void print_dipeptide_frequeces(void);
int Three2IntegerProteinId(char ThreeLeterAmno[3]);


void peters_funktion(char *dir, char *fname);

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
      printf("%i It's %s\n",n_proteins,dp->d_name);
      peters_funktion(dir,dp->d_name);
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

int getatomchain(FILE *fi, char chain)/* || (line[21]==' ' && chain='0' */
{
  char *cp;
  while ((cp=fgets(line,999,fi))!=NULL && strstr(line,"ATOM")!=line && cp!=NULL && !(line[21]==chain || chain=='0'));
  if (cp==NULL) return 0;
  return 1;
}
void peters_funktion(char *dir, char *fname)
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
	  if(dipeptide==0){fprintf(farrayo,"%s %c %i ", fname,chain_name,no-start_chain+1-link_count);}
	  if(dipeptide==1){fprintf(farrayo,"%s ",fname);}
	  nres=link_count;
	  printf("  number of residues in chain %c is %i\n",chain_name,link_count);
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
	      /* printing pdb entry  fprintf(farrayo,"%c%c%c%c:%c  %i %i",line[72],line[73],line[74],line[75],chain_name,link_count,no-start_chain+1-link_count);            */
	      if(dipeptide==0){fprintf(farrayo,"%s %c %i ", fname,chain_name,no-start_chain+1-link_count);}
	      if(dipeptide==1){fprintf(farrayo,"%s ",fname);}
	     
	      nres=link_count;  
	      printf("  number of residues in chain %c is %i\n",chain_name,link_count);
	      GaussIntegrals();
	    }  
	  link_count=0;
	}
      link_count++;
      if (link_count>MAXLENGTH-1) {
	printf(" >MAXLENGTH=%i residues\n",MAXLENGTH);
	fprintf(ferroro, "WARNING!  Chain %c of <%s> is longer than MAXLENGTH=%i. \n          No further attempt to calculate Gauss Integrals for chains of\n          <%s> has been made.\n", chain_name, fname, MAXLENGTH, fname);fclose(fi); return;}
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
  fclose(fi);
  return;
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
/*  return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);*/
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

/* determinant: gives the determinant of vec1, vec2, and vec3 in R^3
 */
double determinant(double vec1[3],double vec2[3],double vec3[3])
{
  double crossvec[3];
  cross(crossvec,vec2,vec3);
  return dot(vec1,crossvec);
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


/* intersectiontrianglelinesegment returns 
 * 1 (one) if the triangle
 * with vertices in the first three points intersects the line
 * segment with vertices in the last two points.
 * 0 if not and 
 * 2 if indecific
 */

int intersectiontrianglelinesegment(double a[3],double b[3],double c[3],double p1[3],double p2[3])
{
  int i;
  double det, det1, det2, det3;
  double crossvec[3], vec1[3], vec2[3], vec3[3];
  double epsilon=pow(10,-15);

  geovec(vec1,p2,p1);
  geovec(vec2,a,c);
  geovec(vec3,b,c);
  cross(crossvec,vec2,vec3);
  det=dot(vec1,crossvec);

  if(fabs(det)<epsilon){fprintf(ferroro,"Warning - numerical problems in intersectiontrianglelinesegment det=%e  , epsilon=%e  ",fabs(det), epsilon); return 2;}
  geovec(vec1,p2,c);
  det1=dot(vec1,crossvec);
  if(det1/det>1 || det1/det<0){return 0;}

  geovec(vec1,p2,p1);
  geovec(vec2,p2,c);
  cross(crossvec,vec2,vec3);
  det2=dot(vec1,crossvec);
  if(det2/det>1 || det2/det<0){return 0;}

  geovec(vec2,a,c);
  geovec(vec3,p2,c);
  cross(crossvec,vec2,vec3);
  det3=dot(vec1,crossvec);
  if((det3/det)<=1 && det3/det>=0 && 1-(det2+det3)/det<=1 && 1-(det2+det3)/det>=0){return 1;}
  return 0;
}

void PerturbBackBone(double ampitude)
{
  int i,j;
  for(j=1;j<nres+1;j++){
    for(i=0;i<3;i++){
      polyl[j][i]+= ampitude*rnd();
    }
  }
}

/*
void geovec(double vec2_minus_vec1[3], double vec1[3], double vec2[3]);
void unit(double uvec[3],double vec[3]);
void cross(double crossvec[3],double vec1[3],double vec2[3]);
*/

void OrtogonalPerturbBackBone(double ampitude)
{
  int i,j;
  double backorthogonalpoly[MAXLENGTH][3], crossvec[3],v1[3],v2[3],unitort[3];
  double rndnbr;
  for(j=1;j<=nres;j++){
    for(i=0;i<3;i++){
      backorthogonalpoly[j][i]=polyl[j][i];
    }}
  
  
  for(j=2;j<nres;j++){
    geovec(v1,backorthogonalpoly[j-1],backorthogonalpoly[j]);
    geovec(v2,backorthogonalpoly[j],backorthogonalpoly[j+1]);
    cross(crossvec, v1, v2);
    unit(unitort, crossvec);
    rndnbr=rnd();
    for(i=0;i<3;i++){
      polyl[j][i]+= ampitude*rndnbr*unitort[i];
    }
  }
}



/* SmoothenTheBackbone smoothens the backbone if the deformation of
 * backbone can be done without selfintersections
 */
void SmoothenTheBackbone(void)
{
  int i, j, k, sum;
  extern double polyl[MAXLENGTH][3];
  extern int nres;
  double varpoly[MAXLENGTH][3], newpoint[3], newpointback[3];
  double tmax, tmin, tnew;  

  for(i=0;i<3;i++){varpoly[1][i]=polyl[1][i];}
  for(i=0;i<3;i++){varpoly[2][i]=polyl[2][i];}
  
  for(j=3;j<=nres-2;j++){
    for(i=0;i<3;i++){
      newpoint[i]=(polyl[j-2][i]+2.4*polyl[j-1][i]+2.1*polyl[j][i]+2.4*polyl[j+1][i]+polyl[j+2][i])/(2+2*2.4+2.1);
    }
    sum=0;
    for(i=1;i<j-2;i++){
      sum+=intersectiontrianglelinesegment(varpoly[j-1],polyl[j],newpoint,varpoly[i],varpoly[i+1]);
    }
    for(i=j+1;i<nres;i++){
      sum+=intersectiontrianglelinesegment(varpoly[j-1],polyl[j],newpoint,polyl[i],polyl[i+1]);
    }
    
    for(i=1;i<j-2;i++){
      sum+=intersectiontrianglelinesegment(polyl[j],polyl[j+1],newpoint,varpoly[i],varpoly[i+1]);
    }
    for(i=j+2;i<nres;i++){
      sum+=intersectiontrianglelinesegment(varpoly[j-1],polyl[j],newpoint,polyl[i],polyl[i+1]);
    }
    
    if(sum==0){
      for(i=0;i<3;i++){varpoly[j][i]=newpoint[i];}
    }
    else{
      tmax=1; tmin=0;
      for(i=0;i<3;i++){newpointback[i]=newpoint[i];}
      
      for(k=0;k<10;k++){
	tnew=(tmax+tmin)*.5;
	printf("tnew=%e\n",tnew);
	for(i=0;i<3;i++){newpoint[i]=tnew*newpointback[i]+(1-tnew)*polyl[j][i];}
	
	sum=0;
	for(i=1;i<j-2;i++){
	  sum+=intersectiontrianglelinesegment(varpoly[j-1],polyl[j],newpoint,varpoly[i],varpoly[i+1]);
	}
	for(i=j+1;i<nres;i++){
	  sum+=intersectiontrianglelinesegment(varpoly[j-1],polyl[j],newpoint,polyl[i],polyl[i+1]);
	}
	
	for(i=1;i<j-2;i++){
	  sum+=intersectiontrianglelinesegment(polyl[j],polyl[j+1],newpoint,varpoly[i],varpoly[i+1]);
	}
	for(i=j+2;i<nres;i++){
	  sum+=intersectiontrianglelinesegment(varpoly[j-1],polyl[j],newpoint,polyl[i],polyl[i+1]);
	}
	fprintf(ferroro,"sum=%i\n",sum);
	if(sum>1){tmax=tnew;}
	if(sum==0){tmin=tnew;}
      }
      tnew=tmin*.9;
      for(i=0;i<3;i++){varpoly[j][i]=tnew*newpointback[i]+(1-tnew)*polyl[j][i];}
    }
  }
  for(j=3;j<=nres-2;j++){
    for(i=0;i<3;i++){
      polyl[j][i]=varpoly[j][i];
    }
  }
}




/* SPHERICAL GEOMETRY  */

/* spangle: For  vec1, vec2, and, vec3 on the unit 2-sphere 
 *          spangle gives the exterior angle 
 *          spangle(vec1,vec2,vec3) between the geodesic segments 
 *          from vec1 to vec2 and from vec2 to vec3. 
 *          Notice, that there is no normalization
 *          before atan2 since atan2 does not need sin and cos
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
 *
 *tempm=partsum[a][c-1]+partsum[b+1][d];
      tempp=partsum[a][d]+partsum[b+1][c-1];
 *
 */
double mixedsum(int a,int b,int c,int d)
{
  double tempm, tempp;
  tempm=partsum[a][c-1]+partsum[b+1][d];
  tempp=partsum[a][d]+partsum[b+1][c-1];
  /*  tempm=partsum[a][c]+partsum[b][d];
  tempp=partsum[a][d]+partsum[b][c];
  */
  return tempp-tempm;
}

double absmixedsum(int a,int b,int c,int d)
{
  double tempm, tempp;
  tempm=abspartsum[a][c-1]+abspartsum[b+1][d];
  tempp=abspartsum[a][d]+abspartsum[b+1][c-1];
  return tempp-tempm;
}




/*            CALCULATION OF THE INVARIANTS
 */

/******** FIRST ORDER INVARIANTS **********************/

double int12(void)
{
  return partsum[1][nres-1]/(twoPI);
}

double inta12(void)
{
  return abspartsum[1][nres-1]/(twoPI);
}

/********* SECOND ORDER INVAIANTS ***********************/


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

/**************THIRD ORDER INVARIANTS********************/


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

/******************SPLINE DEFORMATIONS***************************/


double splinefun(double t)
{
  if(t<0) return 0;
  if(t<1) return t*t*t/6.0;
  if(t<2) return 2.0/3.0-2*t+2*t*t-t*t*t/2.0;
  if(t<3) return 2.0/3.0-2*(4-t)+2*(4-t)*(4-t)-(4-t)*(4-t)*(4-t)/2.0;
  if(t<4) return (4-t)*(4-t)*(4-t)/6.0;
  return 0;    
}

double presplinedeformationsmoothorg(double measurevalue, int measurenumber)
{
  double CoefSmoothOrg[30][16]={{-3.197436e+02, -2.883224e+02, -2.802069e+02, -2.744320e+02, -2.667900e+02, -2.604485e+02, -2.530480e+02, -2.449718e+02, -2.386135e+02, -2.313762e+02, -2.251291e+02, -2.185578e+02, -2.125435e+02, -2.059686e+02, -1.968521e+02, 0.000000e+00},/* 
 */{-4.323760e+02, -3.795051e+02, -3.723412e+02, -3.637384e+02, -3.554168e+02, -3.457790e+02, -3.373883e+02, -3.279643e+02, -3.182622e+02, -3.079369e+02, -2.979298e+02, -2.864793e+02, -2.763838e+02, -2.636304e+02, -2.558967e+02, 0.000000e+00},/* 
 */{-2.271238e+02, -1.716181e+02, -1.535303e+02, -1.373270e+02, -1.217554e+02, -1.056342e+02, -9.317003e+01, -7.909398e+01, -6.736479e+01, -5.558477e+01, -4.461352e+01, -3.432521e+01, -2.195757e+01, -1.358209e+01, -2.432806e+00, 0.000000e+00},/* 
 */{-4.503526e+02, -6.259510e+01, -6.063335e+01, -5.084926e+01, -4.203841e+01, -3.226741e+01, -2.660817e+01, -2.090845e+01, -1.603911e+01, -1.193387e+01, -8.198907e+00, -4.320842e+00, -1.001338e+00, 1.530018e+00, 6.845232e+00, 0.000000e+00},/* 
 */{-3.822391e+02, -3.408675e+02, -3.375999e+02, -3.320636e+02, -3.265147e+02, -3.210941e+02, -3.156815e+02, -3.089459e+02, -3.017987e+02, -2.925882e+02, -2.811764e+02, -2.726509e+02, -2.635936e+02, -2.563481e+02, -2.505648e+02, 0.000000e+00},/* 
 */{-3.548247e+02, -2.855791e+02, -2.808305e+02, -2.746463e+02, -2.693253e+02, -2.638404e+02, -2.563025e+02, -2.490772e+02, -2.412517e+02, -2.303537e+02, -2.220719e+02, -2.131713e+02, -2.029255e+02, -1.964540e+02, -1.912145e+02, 0.000000e+00},/* 
 */{-1.436239e+03, -3.720038e+02, -3.350692e+02, -3.160579e+02, -3.055444e+02, -2.931797e+02, -2.847378e+02, -2.759497e+02, -2.686665e+02, -2.623692e+02, -2.543896e+02, -2.490555e+02, -2.446986e+02, -2.396938e+02, -2.385790e+02, 0.000000e+00},/* 
 */{-3.961375e+02, -1.628548e+02, -1.616330e+02, -1.564161e+02, -1.525293e+02, -1.470919e+02, -1.421920e+02, -1.374435e+02, -1.337700e+02, -1.299775e+02, -1.259029e+02, -1.219635e+02, -1.169896e+02, -1.157990e+02, -1.097794e+02, 0.000000e+00},/* 
 */{-8.066202e+02, -6.517772e+02, -6.462230e+02, -6.372391e+02, -6.278323e+02, -6.176387e+02, -6.064545e+02, -5.928883e+02, -5.786677e+02, -5.627565e+02, -5.469211e+02, -5.337980e+02, -5.164564e+02, -5.007353e+02, -4.963110e+02, 0.000000e+00},/* 
 */{-2.460880e+02, -7.152731e+02, -6.988728e+02, -6.894818e+02, -6.790632e+02, -6.667783e+02, -6.540722e+02, -6.359570e+02, -6.183053e+02, -5.980434e+02, -5.832666e+02, -5.698090e+02, -5.545383e+02, -5.355806e+02, -5.328104e+02, 0.000000e+00},/* 
 */{-3.579592e+03, -4.408524e+02, -3.359427e+02, -3.058190e+02, -2.861726e+02, -2.731977e+02, -2.630218e+02, -2.538499e+02, -2.445866e+02, -2.360675e+02, -2.295014e+02, -2.223295e+02, -2.167535e+02, -2.106743e+02, -2.072773e+02, 0.000000e+00},/* 
 */{1.560031e+02, 8.465603e+01, 8.812625e+01, 9.469884e+01, 1.059010e+02, 1.156418e+02, 1.225411e+02, 1.287197e+02, 1.349452e+02, 1.389226e+02, 1.446053e+02, 1.484463e+02, 1.529033e+02, 1.559380e+02, 1.650088e+02, 0.000000e+00},/* 
 */{-7.866648e+01, -1.714202e+02, -1.660737e+02, -1.603419e+02, -1.551662e+02, -1.494192e+02, -1.440249e+02, -1.369109e+02, -1.286982e+02, -1.211608e+02, -1.077709e+02, -9.218865e+01, -8.189200e+01, -6.986134e+01, -5.982144e+01, 0.000000e+00},/* 
 */{-2.359462e+02, -1.131961e+03, -1.098386e+03, -1.093354e+03, -1.080364e+03, -1.069432e+03, -1.058697e+03, -1.040698e+03, -1.021851e+03, -9.957950e+02, -9.650870e+02, -9.478457e+02, -9.311873e+02, -9.134451e+02, -9.087646e+02, 0.000000e+00},/* 
 */{-1.619740e+03, -3.699248e+02, -3.147734e+02, -2.805806e+02, -2.642955e+02, -2.494710e+02, -2.384339e+02, -2.277311e+02, -2.202632e+02, -2.126506e+02, -2.064482e+02, -2.000697e+02, -1.938865e+02, -1.885999e+02, -1.837621e+02, 0.000000e+00},/* 
 */{-1.015086e+02, -2.465279e+01, -2.485352e+01, -2.181504e+01, -1.892179e+01, -1.653119e+01, -1.283194e+01, -8.428750e+00, -4.625742e-01, 1.052858e+01, 1.647150e+01, 2.004696e+01, 2.451622e+01, 2.915077e+01, 3.347932e+01, 0.000000e+00},/* 
 */{-2.486634e+02, -4.191579e+02, -4.089102e+02, -4.092159e+02, -4.028386e+02, -3.980484e+02, -3.924299e+02, -3.845550e+02, -3.742693e+02, -3.706606e+02, -3.645034e+02, -3.594337e+02, -3.563740e+02, -3.516078e+02, -3.506483e+02, 0.000000e+00},/* 
 */{7.169783e+01, -1.859575e+01, -1.356244e+01, -1.143467e+01, -8.497706e+00, -4.490052e+00, -1.081137e-01, 6.759394e+00, 1.093239e+01, 2.049624e+01, 3.023360e+01, 3.547064e+01, 4.044911e+01, 4.489441e+01, 4.886720e+01, 0.000000e+00},/* 
 */{2.139739e+02, 1.539314e+02, 1.571150e+02, 1.610361e+02, 1.643516e+02, 1.704245e+02, 1.762477e+02, 1.865355e+02, 1.914519e+02, 1.971046e+02, 2.020277e+02, 2.064675e+02, 2.093675e+02, 2.116468e+02, 2.239609e+02, 0.000000e+00},/* 
 */{-2.559139e+02, -1.508357e+02, -1.486667e+02, -1.458230e+02, -1.426569e+02, -1.391657e+02, -1.362815e+02, -1.316947e+02, -1.273580e+02, -1.237468e+02, -1.204404e+02, -1.182622e+02, -1.156993e+02, -1.121576e+02, -1.134177e+02, 0.000000e+00},/* 
 */{3.320789e+01, -1.260140e+01, -7.354663e+00, -4.334913e+00, -1.078672e+00, 3.828298e+00, 8.423274e+00, 1.690975e+01, 1.889181e+01, 2.424566e+01, 2.742190e+01, 3.141993e+01, 3.356502e+01, 3.736835e+01, 4.017256e+01, 0.000000e+00},/* 
 */{-3.368250e+02, -1.086963e+02, -1.101313e+02, -1.046192e+02, -1.018014e+02, -9.643735e+01, -9.237054e+01, -8.653599e+01, -7.951754e+01, -6.896251e+01, -6.258895e+01, -5.696923e+01, -5.307860e+01, -4.739636e+01, -4.654458e+01, 0.000000e+00},/* 
 */{3.220168e+01, 1.215428e+01, 1.283006e+01, 1.619144e+01, 1.962525e+01, 2.393409e+01, 3.082398e+01, 3.479396e+01, 3.813276e+01, 4.234134e+01, 4.560099e+01, 4.856591e+01, 5.056359e+01, 5.425393e+01, 5.548827e+01, 0.000000e+00},/* 
 */{-1.870530e+01, -2.798152e+01, -2.495826e+01, -2.264876e+01, -1.944912e+01, -1.512430e+01, -1.200502e+01, -7.226342e+00, -3.727455e+00, 3.294273e-01, 3.519920e+00, 5.965627e+00, 8.932295e+00, 1.289682e+01, 1.410634e+01, 0.000000e+00},/* 
 */{-7.351438e+01, -9.954315e+00, -7.106551e+00, -3.510217e+00, 1.598041e+00, 5.986602e+00, 1.181914e+01, 1.931187e+01, 2.317409e+01, 2.772706e+01, 3.020093e+01, 3.265570e+01, 3.597612e+01, 3.798790e+01, 4.071223e+01, 0.000000e+00},/* 
 */{-1.541581e+02, 1.614397e+01, 1.622479e+01, 1.961916e+01, 2.246212e+01, 2.570981e+01, 2.971659e+01, 3.312190e+01, 3.892171e+01, 4.527941e+01, 4.967889e+01, 5.375146e+01, 5.719296e+01, 5.935404e+01, 6.416701e+01, 0.000000e+00},/* 
 */{-2.829540e+02, -4.107829e+02, -4.056872e+02, -4.032253e+02, -3.991498e+02, -3.979315e+02, -3.933413e+02, -3.893524e+02, -3.824630e+02, -3.764142e+02, -3.679231e+02, -3.616823e+02, -3.571693e+02, -3.532505e+02, -3.511791e+02, 0.000000e+00},/* 
 */{-9.758395e+01, -2.614835e+01, -2.559130e+01, -2.017768e+01, -1.736653e+01, -1.154634e+01, -6.962113e+00, -4.259342e+00, 9.103369e+00, 2.185224e+01, 2.409766e+01, 3.206438e+01, 3.611873e+01, 4.016034e+01, 4.305715e+01, 0.000000e+00},/* 
 */{-9.107451e+01, -8.292318e+01, -7.838891e+01, -7.835102e+01, -7.359929e+01, -7.155745e+01, -6.622321e+01, -6.112757e+01, -5.528490e+01, -5.038604e+01, -4.587816e+01, -4.171027e+01, -4.082923e+01, -3.746770e+01, -3.688622e+01, 0.000000e+00},/* 
  */{-1.760486e+02, -5.053728e+01, -5.124014e+01, -4.654084e+01, -4.394787e+01, -4.059502e+01, -3.662170e+01, -3.178512e+01, -2.655638e+01, -1.818261e+01, -1.135904e+01, 1.919248e+00, 1.439862e+01, 1.914439e+01, 2.730479e+01, 0.000000e+00}} ;

  double knotsstartend[30][2]={{1.543090e+00, 3.388490e+00},/* 
 */{-1.156866e-01, 1.097668e-01},/* 
 */{-7.055847e-03, 5.420859e-02},/* 
 */{-7.931673e-04, 1.649079e-03},/* 
 */{-4.893340e-04, 3.492860e-04},/* 
 */{-4.361313e-04, 3.548503e-04},/* 
 */{-1.906141e-05, 8.196310e-05},/* 
 */{-6.992287e-04, 1.669920e-03},/* 
 */{-5.420780e-03, 4.712220e-03},/* 
 */{-3.584573e-03, 3.791193e-03},/* 
 */{-3.611285e-04, 1.480395e-03},/* 
 */{-1.965247e-03, 4.791937e-03},/* 
 */{-8.625597e-04, 5.382987e-04},/* 
 */{-1.999305e-03, 1.351821e-03},/* 
 */{-9.790955e-05, 4.140336e-04},/* 
 */{-3.514717e-06, 3.021167e-06},/* 
 */{-2.890233e-05, 3.497233e-05},/* 
 */{-1.021406e-04, 7.640340e-05},/* 
 */{-9.727883e-06, 1.213748e-05},/* 
 */{-2.635703e-05, 2.656313e-05},/* 
 */{-2.306577e-04, 2.954257e-04},/* 
 */{-2.916900e-05, 2.454150e-05},/* 
 */{-4.524907e-04, 7.236277e-04},/* 
 */{-5.537503e-04, 5.786913e-04},/* 
 */{-9.869520e-05, 1.422288e-04},/* 
 */{-1.583209e-05, 1.310987e-05},/* 
 */{-1.382874e-04, 9.903297e-05},/* 
 */{-1.148076e-05, 1.038139e-05},/* 
 */{-3.913683e-05, 4.114933e-05},/* 
 */{-2.182509e-04, 1.190886e-04}};

  double prepost[30][6]={{1.912170e+00, -2.820608e+02, 8.064372e+01, 3.019410e+00, -2.036781e+02, 1.241466e+02},/* 
 */{-7.059590e-02, -3.737121e+02, 7.195778e+02, 6.467610e-02, -2.621812e+02, 1.205403e+03},/* 
 */{5.197040e-03, -1.573586e+02, 5.654660e+03, 4.195570e-02, -1.112716e+01, 3.128396e+03},/* 
 */{-3.047180e-04, -6.087259e+01, 9.476334e+04, 1.160630e-03, 2.814403e+00, 3.252834e+04},/* 
 */{-3.216100e-04, -3.381106e+02, 1.077460e+05, 1.815620e-04, -2.549898e+02, 2.525465e+05},/* 
 */{-2.779350e-04, -2.817437e+02, 1.560857e+05, 1.966540e-04, -1.952633e+02, 2.246726e+05},/* 
 */{1.143490e-06, -3.453553e+02, 9.731883e+06, 6.175820e-05, -2.394862e+02, 1.336995e+06},/* 
 */{-2.253990e-04, -1.618459e+02, 5.727285e+04, 1.196090e-03, -1.140447e+02, 5.252556e+04},/* 
 */{-3.394180e-03, -6.472402e+02, 1.807092e+04, 2.685620e-03, -5.001585e+02, 3.352050e+04},/* 
 */{-2.109420e-03, -7.021035e+02, 1.145197e+04, 2.316040e-03, -5.357012e+02, 4.778736e+04},/* 
 */{7.176110e-06, -3.673515e+02, 1.461530e+06, 1.112090e-03, -2.099520e+02, 8.597539e+04},/* 
 */{-6.138100e-04, 8.779670e+01, 8.341635e+03, 3.440500e-03, 1.580351e+02, 1.028460e+04},/* 
 */{-5.823880e-04, -1.669793e+02, 4.757696e+04, 2.581270e-04, -6.795686e+01, 1.583164e+05},/* 
 */{-1.329080e-03, -1.106297e+03, 3.316582e+04, 6.815960e-04, -9.124181e+02, 1.527535e+05},/* 
 */{4.479090e-06, -3.291850e+02, 2.684713e+06, 3.116450e-04, -1.874320e+02, 3.176879e+05},/* 
 */{-2.207540e-06, -2.463970e+01, 6.909846e+06, 1.713990e-06, 2.993996e+01, 1.066970e+07},/* 
 */{-1.612740e-05, -4.116193e+02, 9.684580e+05, 2.219740e-05, -3.512745e+02, 2.728113e+06},/* 
 */{-6.643180e-05, -1.468994e+01, 2.337964e+05, 4.069460e-05, 4.557820e+01, 3.263996e+05},/* 
 */{-5.354810e-06, 1.566255e+02, 1.872643e+06, 7.764410e-06, 2.146509e+02, 3.803804e+06},/* 
 */{-1.577300e-05, -1.491801e+02, 1.595213e+06, 1.597910e-05, -1.126664e+02, 9.227285e+05},/* 
 */{-1.254410e-04, -8.525972e+00, 1.268746e+05, 1.902090e-04, 3.778663e+01, 8.178060e+04},/* 
 */{-1.842690e-05, -1.095576e+02, 1.902931e+06, 1.379940e-05, -4.757728e+01, 1.171395e+06},/* 
 */{-2.172670e-04, 1.295171e+01, 1.810487e+04, 4.884040e-04, 5.421559e+01, 1.431010e+04},/* 
 */{-3.272620e-04, -2.560742e+01, 4.276463e+04, 3.522030e-04, 1.288321e+01, 2.994523e+04},/* 
 */{-5.051040e-05, -7.693175e+00, 3.392335e+05, 9.404400e-05, 3.853565e+01, 1.283562e+05},/* 
 */{-1.004370e-05, 1.626442e+01, 2.984113e+06, 7.321480e-06, 6.045096e+01, 1.696862e+06},/* 
 */{-9.082330e-05, -4.067538e+02, 1.276477e+05, 5.156890e-05, -3.525284e+02, 7.721599e+05},/* 
 */{-7.108330e-06, -2.538279e+01, 3.066748e+06, 6.008960e-06, 4.058074e+01, 2.037188e+06},/* 
 */{-2.307960e-05, -7.968427e+01, 7.496635e+05, 2.509210e-05, -3.754023e+01, 5.199738e+05},/* 
											    */{-1.507830e-04, -5.080687e+01, 1.999535e+05, 5.162070e-05, 2.102056e+01, 3.581254e+05}};/*this is prepostsmoothorg*/

  double ss;
  int a;
  if (measurenumber>29){
    printf("WARNING measure number is to high in smooth original\n");
    return -99999999999.9 ;
  }
  if (measurevalue<prepost[measurenumber][0]){
    return prepost[measurenumber][2]*(measurevalue-prepost[measurenumber][0])+prepost[measurenumber][1];
  }
  if (measurevalue>prepost[measurenumber][3]){
    return prepost[measurenumber][5]*(measurevalue-prepost[measurenumber][3])+prepost[measurenumber][4];
  }
  ss=1+19*(measurevalue-knotsstartend[measurenumber][0])/(knotsstartend[measurenumber][1]-knotsstartend[measurenumber][0]);
  a=floor(ss);
  if(a>16) {
    printf("wRNING PROBLEM IN PREDEFORMATION AS a>16\n");
    return 0;
    }
  
  return splinefun(ss-a+3)*CoefSmoothOrg[measurenumber][a-3-1]+splinefun(ss-a+2)*CoefSmoothOrg[measurenumber][a-2-1]+splinefun(ss-a+1)*CoefSmoothOrg[measurenumber][a-1-1]+splinefun(ss-a)*CoefSmoothOrg[measurenumber][a-1];
}


double presplinedeformationsmoothnew(double measurevalue, int measurenumber)
{
  double CoefSmoothNew[30][16]={{-8.944632e+01, -1.099187e+02, -1.033060e+02, -9.569429e+01, -9.012278e+01, -8.301296e+01, -7.623848e+01, -6.786135e+01, -6.141938e+01, -5.487420e+01, -4.849957e+01, -4.272373e+01, -3.644640e+01, -3.069523e+01, -2.202762e+01, 0.000000e+00},/* 
 */{-4.323760e+02, -3.795051e+02, -3.723412e+02, -3.637384e+02, -3.554168e+02, -3.457790e+02, -3.373883e+02, -3.279643e+02, -3.182622e+02, -3.079369e+02, -2.979298e+02, -2.864793e+02, -2.763838e+02, -2.636304e+02, -2.558967e+02, 0.000000e+00},/* 
 */{-2.206509e+02, -4.086930e+02, -3.885715e+02, -3.694945e+02, -3.516591e+02, -3.366257e+02, -3.230932e+02, -3.116065e+02, -2.991483e+02, -2.879707e+02, -2.768270e+02, -2.666279e+02, -2.579784e+02, -2.493960e+02, -2.408323e+02, 0.000000e+00},/* 
 */{-1.627798e+03, -1.653651e+03, -1.649482e+03, -1.646591e+03, -1.642891e+03, -1.639919e+03, -1.636694e+03, -1.633282e+03, -1.628551e+03, -1.624535e+03, -1.619286e+03, -1.613368e+03, -1.605346e+03, -1.596190e+03, -1.563467e+03, 0.000000e+00},/* 
 */{-3.922498e+02, -3.522696e+02, -3.456155e+02, -3.401583e+02, -3.336302e+02, -3.269996e+02, -3.202883e+02, -3.113555e+02, -3.029190e+02, -2.907618e+02, -2.786194e+02, -2.693874e+02, -2.604276e+02, -2.505540e+02, -2.439578e+02, 0.000000e+00},/* 
 */{-1.725072e+02, -2.552874e+02, -2.493978e+02, -2.424006e+02, -2.355072e+02, -2.286147e+02, -2.199379e+02, -2.111326e+02, -2.014180e+02, -1.875711e+02, -1.783861e+02, -1.677901e+02, -1.594134e+02, -1.486498e+02, -1.420985e+02, 0.000000e+00},/* 
 */{-5.981506e+03, -2.807756e+02, -1.738971e+02, -1.479372e+02, -1.246978e+02, -1.087282e+02, -9.485573e+01, -8.212603e+01, -7.035389e+01, -6.069769e+01, -5.144344e+01, -4.218746e+01, -3.386648e+01, -2.608900e+01, -1.839906e+01, 0.000000e+00},/* 
 */{-1.763959e+03, -1.700768e+03, -1.699934e+03, -1.696725e+03, -1.693572e+03, -1.690549e+03, -1.686964e+03, -1.683650e+03, -1.679759e+03, -1.675204e+03, -1.669710e+03, -1.664418e+03, -1.656446e+03, -1.647501e+03, -1.615750e+03, 0.000000e+00},/* 
 */{-2.423587e+02, -1.787628e+02, -1.717232e+02, -1.628667e+02, -1.532275e+02, -1.434835e+02, -1.316439e+02, -1.183518e+02, -1.001996e+02, -7.680384e+01, -6.130349e+01, -4.697156e+01, -3.250732e+01, -2.227616e+01, -7.768331e+00, 0.000000e+00},/* 
 */{-2.255406e+02, -2.578660e+02, -2.487574e+02, -2.403399e+02, -2.308197e+02, -2.212080e+02, -2.092860e+02, -1.973281e+02, -1.781195e+02, -1.543145e+02, -1.396794e+02, -1.257571e+02, -1.106626e+02, -9.918504e+01, -8.721038e+01, 0.000000e+00},/* 
 */{-2.022817e+04, -3.169053e+02, -1.568568e+02, -1.205469e+02, -9.492570e+01, -7.663713e+01, -6.123187e+01, -4.754060e+01, -3.603088e+01, -2.532938e+01, -1.666340e+01, -7.449999e+00, -9.321766e-02, 8.757506e+00, 1.448935e+01, 0.000000e+00},/* 
 */{-1.691708e+03, -1.701541e+03, -1.698722e+03, -1.695657e+03, -1.692697e+03, -1.689432e+03, -1.685727e+03, -1.682531e+03, -1.678279e+03, -1.674045e+03, -1.668536e+03, -1.663085e+03, -1.655342e+03, -1.646076e+03, -1.614782e+03, 0.000000e+00},/* 
 */{1.208610e+01, -2.402713e+01, -1.594536e+01, -9.689928e+00, -2.729414e+00, 4.469296e+00, 1.330338e+01, 2.366246e+01, 3.251238e+01, 5.075630e+01, 6.467103e+01, 7.697915e+01, 8.706296e+01, 9.778898e+01, 1.075953e+02, 0.000000e+00},/* 
 */{-4.187916e+02, -3.763817e+02, -3.686820e+02, -3.595731e+02, -3.496542e+02, -3.369328e+02, -3.265093e+02, -3.091941e+02, -2.964498e+02, -2.712538e+02, -2.519308e+02, -2.318562e+02, -2.166059e+02, -1.988975e+02, -1.848086e+02, 0.000000e+00},/* 
 */{-1.617682e+04, -3.479836e+02, -1.795255e+02, -1.410600e+02, -1.134702e+02, -9.449111e+01, -7.734345e+01, -6.426249e+01, -5.120941e+01, -3.990345e+01, -3.076247e+01, -2.157896e+01, -1.321230e+01, -5.211887e+00, 2.654697e+00, 0.000000e+00},/* 
 */{-2.038454e+02, -1.605521e+02, -1.585855e+02, -1.533240e+02, -1.499669e+02, -1.442162e+02, -1.391387e+02, -1.342552e+02, -1.253144e+02, -1.183658e+02, -1.101778e+02, -1.032587e+02, -9.876589e+01, -9.222172e+01, -9.156987e+01, 0.000000e+00},/* 
 */{1.194045e+02, -4.640461e+01, -3.937670e+01, -3.665371e+01, -3.374323e+01, -2.981953e+01, -2.632391e+01, -2.190136e+01, -1.655628e+01, -1.393417e+01, -7.693308e+00, -4.414165e+00, -1.166574e+00, 2.177935e+00, 5.364475e+00, 0.000000e+00},/* 
 */{4.914698e+01, -6.925862e+01, -6.115743e+01, -6.098089e+01, -5.765389e+01, -5.392566e+01, -4.914807e+01, -4.474478e+01, -3.638328e+01, -2.925467e+01, -2.327735e+01, -1.854936e+01, -1.488345e+01, -1.066282e+01, -8.361774e+00, 0.000000e+00},/* 
 */{1.005746e+02, 5.750706e+01, 6.033379e+01, 6.281674e+01, 6.558924e+01, 6.830041e+01, 7.231943e+01, 7.540429e+01, 8.169666e+01, 8.545651e+01, 9.081339e+01, 9.404398e+01, 9.857993e+01, 9.944875e+01, 1.056589e+02, 0.000000e+00},/* 
 */{-7.210167e+01, -9.239858e+01, -8.935188e+01, -8.628471e+01, -8.396635e+01, -8.055506e+01, -7.638363e+01, -7.247419e+01, -6.838569e+01, -6.336126e+01, -5.960200e+01, -5.513016e+01, -5.195936e+01, -4.839843e+01, -4.609455e+01, 0.000000e+00},/* 
 */{7.054937e+01, 1.309963e+02, 1.329702e+02, 1.371429e+02, 1.408555e+02, 1.469020e+02, 1.502321e+02, 1.569482e+02, 1.616573e+02, 1.648894e+02, 1.696164e+02, 1.729594e+02, 1.767021e+02, 1.790866e+02, 1.862448e+02, 0.000000e+00},/* 
 */{-6.383472e+01, -4.518012e+01, -4.322435e+01, -4.079681e+01, -3.735081e+01, -3.362176e+01, -2.855168e+01, -2.375436e+01, -1.479947e+01, -9.443971e+00, -3.288144e+00, -1.063605e+00, 2.357577e+00, 5.503577e+00, 9.837346e+00, 0.000000e+00},/* 
 */{-1.319450e+02, -1.373329e+02, -1.351566e+02, -1.320714e+02, -1.297339e+02, -1.259316e+02, -1.213749e+02, -1.171764e+02, -1.134295e+02, -1.095467e+02, -1.062760e+02, -1.027397e+02, -9.999116e+01, -9.779185e+01, -9.503813e+01, 0.000000e+00},/* 
 */{-8.504188e+01, -7.249933e+01, -6.896536e+01, -6.603217e+01, -6.262269e+01, -5.930658e+01, -5.495930e+01, -5.031085e+01, -4.596648e+01, -4.046198e+01, -3.720763e+01, -3.376877e+01, -2.991338e+01, -2.732565e+01, -2.208567e+01, 0.000000e+00},/* 
 */{-7.572016e+01, -3.680928e+01, -3.193137e+01, -2.793961e+01, -2.308964e+01, -1.705898e+01, -1.097472e+01, -6.306664e+00, -1.679170e+00, 1.158472e+00, 3.907689e+00, 6.529308e+00, 8.724441e+00, 1.055803e+01, 1.386877e+01, 0.000000e+00},/* 
 */{-1.531824e+02, 1.047709e+02, 1.052261e+02, 1.084142e+02, 1.121933e+02, 1.147943e+02, 1.190171e+02, 1.221900e+02, 1.277416e+02, 1.331911e+02, 1.385639e+02, 1.425704e+02, 1.456560e+02, 1.485131e+02, 1.524278e+02, 0.000000e+00},/* 
 */{-3.078750e+02, -2.430555e+02, -2.390526e+02, -2.375404e+02, -2.327583e+02, -2.318849e+02, -2.275092e+02, -2.237270e+02, -2.176982e+02, -2.129316e+02, -2.040267e+02, -1.981699e+02, -1.943631e+02, -1.888787e+02, -1.894915e+02, 0.000000e+00},/* 
 */{-1.114843e+02, -2.497651e+02, -2.409053e+02, -2.361774e+02, -2.290416e+02, -2.194772e+02, -2.079200e+02, -2.000440e+02, -1.940763e+02, -1.891448e+02, -1.846660e+02, -1.812537e+02, -1.780991e+02, -1.748716e+02, -1.734410e+02, 0.000000e+00},/* 
 */{-3.179037e+02, -1.017927e+02, -9.525953e+01, -9.248529e+01, -8.894468e+01, -8.660980e+01, -8.159995e+01, -7.839347e+01, -7.294461e+01, -6.963943e+01, -6.580943e+01, -6.227030e+01, -5.836783e+01, -5.593004e+01, -5.522874e+01, 0.000000e+00},/* 
 */{-1.533690e+02, -1.605210e+02, -1.599891e+02, -1.550206e+02, -1.537064e+02, -1.499839e+02, -1.467191e+02, -1.419631e+02, -1.371803e+02, -1.316141e+02, -1.210307e+02, -1.166583e+02, -1.111252e+02, -1.066626e+02, -1.051770e+02, 0.000000e+00}};

  double knotsstartend[30][2]={{-1.211098e+00, 1.180243e+00},/* 
 */{-1.156866e-01, 1.097668e-01},/* 
 */{-1.145264e-01, 1.028351e-01},/* 
 */{-1.569844e-01, 3.882746e-02},/* 
 */{-3.296010e-03, 2.577440e-03},/* 
 */{-2.685857e-03, 2.372027e-03},/* 
 */{-3.413006e-02, 1.381300e-01},/* 
 */{-4.221381e+00, 1.049125e+00},/* 
 */{-8.223077e-02, 7.266607e-02},/* 
 */{-5.689723e-02, 5.124143e-02},/* 
 */{-5.554268e-01, 2.229844e+00},/* 
 */{-1.404638e+00, 3.487506e-01},/* 
 */{-2.370320e-02, 1.843030e-02},/* 
 */{-7.323573e-03, 5.515793e-03},/* 
 */{-6.256137e-01, 2.516828e+00},/* 
 */{-2.081713e-03, 1.806953e-03},/* 
 */{-3.642003e-04, 2.947163e-04},/* 
 */{-7.760810e-04, 6.564940e-04},/* 
 */{-1.881483e-05, 1.486431e-05},/* 
 */{-8.250873e-05, 6.106043e-05},/* 
 */{-3.195707e-04, 3.580527e-04},/* 
 */{-6.016903e-05, 5.669613e-05},/* 
 */{-5.989873e-04, 7.595493e-04},/* 
 */{-1.482195e-03, 1.344951e-03},/* 
 */{-9.868160e-05, 1.654604e-04},/* 
 */{-1.844706e-05, 1.436125e-05},/* 
 */{-1.199999e-04, 8.086873e-05},/* 
 */{-1.183927e-04, 1.913507e-04},/* 
 */{-3.689207e-05, 3.800077e-05},/* 
 */{-1.644483e-04, 1.118182e-04}};

  double prepost[30][6]={{-7.328300e-01, -1.045072e+02, 5.077594e+01, 7.019750e-01, -2.869502e+01, 6.357496e+01},/* 
 */{-7.059590e-02, -3.737121e+02, 7.195778e+02, 6.467610e-02, -2.621812e+02, 1.205403e+03},/* 
 */{-7.105410e-02, -3.924074e+02, 1.365723e+03, 5.936280e-02, -2.473752e+02, 1.155153e+03},/* 
 */{-1.178220e-01, -1.650384e+03, 3.065597e+02, -3.349070e-04, -1.585593e+03, 5.414025e+03},/* 
 */{-2.121320e-03, -3.470929e+02, 2.244242e+04, 1.402750e-03, -2.491980e+02, 4.008743e+04},/* 
 */{-1.674280e-03, -2.503630e+02, 1.679412e+04, 1.360450e-03, -1.475183e+02, 3.985723e+04},/* 
 */{3.219530e-04, -2.096363e+02, 2.127231e+04, 1.036780e-01, -2.454420e+01, 8.749016e+02},/* 
 */{-3.167280e+00, -1.699981e+03, 1.024202e+01, -4.975990e-03, -1.637093e+03, 2.023582e+02},/* 
 */{-5.125140e-02, -1.730515e+02, 1.073556e+03, 4.168670e-02, -1.901864e+01, 1.595162e+03},/* 
 */{-3.526950e-02, -2.505828e+02, 1.415930e+03, 2.961370e-02, -9.664737e+01, 2.340381e+03},/* 
 */{1.627420e-03, -2.257604e+02, 3.516373e+03, 1.672790e+00, 9.610768e+00, 4.314975e+01},/* 
 */{-1.053960e+00, -1.699248e+03, 2.865658e+01, -1.927050e-03, -1.635826e+03, 6.059005e+02},/* 
 */{-1.527650e-02, -1.765863e+01, 2.982316e+03, 1.000360e-02, 9.951524e+01, 3.496023e+03},/* 
 */{-4.755700e-03, -3.701480e+02, 1.308881e+04, 2.947920e-03, -1.961609e+02, 2.761578e+04},/* 
 */{2.874740e-03, -2.451903e+02, 2.660763e+03, 1.888340e+00, -3.664019e+00, 4.655001e+01},/* 
 */{-1.303980e-03, -1.587528e+02, 1.879904e+04, 1.029220e-03, -9.247294e+01, 2.128217e+04},/* 
 */{-2.324170e-04, -4.091918e+01, 6.325273e+04, 1.629330e-04, 2.790361e+00, 8.841077e+04},/* 
 */{-4.895660e-04, -6.328523e+01, 4.025483e+04, 3.699790e-04, -1.035833e+01, 4.027293e+04},/* 
 */{-1.207900e-05, 5.980030e+01, 1.014809e+06, 8.128480e-06, 1.009974e+02, 1.276975e+06},/* 
 */{-5.379490e-05, -8.992835e+01, 3.422824e+05, 3.234660e-05, -4.798653e+01, 4.740366e+05},/* 
 */{-1.840460e-04, 1.326851e+02, 1.078663e+05, 2.225280e-04, 1.806678e+02, 4.942046e+04},/* 
 */{-3.679600e-05, -4.359751e+01, 3.968127e+05, 3.332310e-05, 6.452792e+00, 5.967129e+05},/* 
 */{-3.272800e-04, -1.355042e+02, 3.238845e+04, 4.878420e-04, -9.707075e+01, 6.184445e+04},/* 
 */{-9.167660e-04, -6.973543e+01, 2.366909e+04, 7.795220e-04, -2.602887e+01, 3.177595e+04},/* 
 */{-4.585320e-05, -3.302794e+01, 3.794361e+05, 1.126320e-04, 1.132332e+01, 1.794291e+05},/* 
 */{-1.188540e-05, 1.050249e+02, 3.752511e+06, 7.799590e-06, 1.491779e+02, 2.602551e+05},/* 
 */{-7.982620e-05, -2.401468e+02, 4.182932e+05, 4.069500e-05, -1.892681e+02, 4.862248e+05},/* 
 */{-5.644400e-05, -2.428337e+02, 2.818479e+05, 1.294020e-04, -1.745095e+02, 3.340530e+05},/* 
 */{-2.191350e-05, -9.716638e+01, 2.415668e+06, 2.302220e-05, -5.586526e+01, 5.955569e+05},/* 
 */{-1.091950e-04, -1.597066e+02, 1.236513e+05, 5.656490e-05, -1.064813e+02, 3.103162e+05}};/*this is prepostsmoothorg*/

  double ss;
  int a;
  if (measurenumber>29){
    printf("WARNING measure number is to high in smooth new\n");
    return -99999999999.9 ;
  }
  if (measurevalue<prepost[measurenumber][0]){
    return prepost[measurenumber][2]*(measurevalue-prepost[measurenumber][0])+prepost[measurenumber][1];
  }
  if (measurevalue>prepost[measurenumber][3]){
    return prepost[measurenumber][5]*(measurevalue-prepost[measurenumber][3])+prepost[measurenumber][4];
  }
  ss=1+19*(measurevalue-knotsstartend[measurenumber][0])/(knotsstartend[measurenumber][1]-knotsstartend[measurenumber][0]);
  a=floor(ss);
  if(a>16) {
    printf("wRNING PROBLEM IN PREDEFORMATION AS a>16\n");
    return 0;
    }
  
  return splinefun(ss-a+3)*CoefSmoothNew[measurenumber][a-3-1]+splinefun(ss-a+2)*CoefSmoothNew[measurenumber][a-2-1]+splinefun(ss-a+1)*CoefSmoothNew[measurenumber][a-1-1]+splinefun(ss-a)*CoefSmoothNew[measurenumber][a-1];
}


double splinedeformationsmoothorg(double measurevalue, int measurenumber)
{
  return presplinedeformationsmoothorg(measurevalue, measurenumber)-presplinedeformationsmoothorg(0.0, measurenumber);
}

double splinedeformationsmoothnew(double measurevalue, int measurenumber)
{
  return presplinedeformationsmoothnew(measurevalue, measurenumber)-presplinedeformationsmoothnew(0.0, measurenumber);
}

double presplinedeformationCalphaorg(double measurevalue, int measurenumber)
{
  double CoefCalphaOrg[29][16]={{-1.022043e+01, -5.673251e+00, -5.109977e+00, -4.515572e+00, -3.920791e+00, -3.315927e+00, -2.694542e+00, -2.072482e+00, -1.459846e+00, -8.388197e-01, -2.126423e-01, 4.441569e-01, 1.154193e+00, 1.859077e+00, 2.844650e+00, 0.000000e+00},/* 
 */{-2.145919e+00, -1.298058e+00, -1.100987e+00, -8.558102e-01, -6.269785e-01, -3.821712e-01, -1.641000e-01, 7.047936e-02, 3.073944e-01, 5.574438e-01, 8.041141e-01, 1.053503e+00, 1.313411e+00, 1.564413e+00, 1.839354e+00, 0.000000e+00},/* 
 */{-3.272319e+01, -5.679388e+00, -3.027266e+00, -2.184208e+00, -1.524601e+00, -1.023559e+00, -5.868286e-01, -1.624693e-01, 2.265211e-01, 5.659684e-01, 9.532884e-01, 1.290577e+00, 1.618150e+00, 1.929699e+00, 2.376741e+00, 0.000000e+00},/* 
 */{-7.094200e+00, -7.711310e+00, -7.252447e+00, -6.629255e+00, -6.059834e+00, -5.538843e+00, -5.032154e+00, -4.544455e+00, -4.071022e+00, -3.629301e+00, -3.170317e+00, -2.731559e+00, -2.311583e+00, -1.869839e+00, -1.501668e+00, 0.000000e+00},/* 
 */{-1.193950e+01, -9.692399e+00, -9.182438e+00, -8.633668e+00, -8.090662e+00, -7.563740e+00, -7.052701e+00, -6.583849e+00, -6.129935e+00, -5.685335e+00, -5.241571e+00, -4.799894e+00, -4.392813e+00, -3.960548e+00, -3.518674e+00, 0.000000e+00},/* 
 */{-1.029792e+01, -7.280154e+00, -6.965742e+00, -6.570639e+00, -6.225599e+00, -5.932427e+00, -5.621598e+00, -5.315214e+00, -5.019830e+00, -4.739372e+00, -4.447961e+00, -4.170449e+00, -3.899608e+00, -3.629545e+00, -3.404594e+00, 0.000000e+00},/* 
 */{-1.083239e+01, -7.655407e+00, -7.314142e+00, -6.977336e+00, -6.711711e+00, -6.450447e+00, -6.181713e+00, -5.948476e+00, -5.646119e+00, -5.414567e+00, -5.151326e+00, -4.937171e+00, -4.734597e+00, -4.522417e+00, -4.150915e+00, 0.000000e+00},/* 
 */{-2.801743e+01, -1.467049e+01, -1.403551e+01, -1.331707e+01, -1.260394e+01, -1.181148e+01, -1.095126e+01, -1.009105e+01, -9.187503e+00, -8.271633e+00, -7.363734e+00, -6.392157e+00, -5.534488e+00, -4.775230e+00, -4.061104e+00, 0.000000e+00},/* 
 */{-1.637997e+01, -1.593480e+01, -1.504607e+01, -1.430765e+01, -1.344713e+01, -1.248755e+01, -1.156632e+01, -1.052623e+01, -9.459432e+00, -8.375925e+00, -7.296309e+00, -6.370489e+00, -5.483069e+00, -4.709157e+00, -3.977414e+00, 0.000000e+00},/* 
 */{-2.134964e+01, -6.187588e+00, -5.376331e+00, -4.447713e+00, -3.787271e+00, -3.205523e+00, -2.730491e+00, -2.273153e+00, -1.884532e+00, -1.524673e+00, -1.194215e+00, -8.477538e-01, -5.420715e-01, -2.672812e-01, 7.954958e-02, 0.000000e+00},/* 
 */{1.123691e+01, -4.912341e-01, 1.571220e-01, 6.084592e-01, 1.157967e+00, 1.721950e+00, 2.399714e+00, 2.927840e+00, 3.386638e+00, 3.902328e+00, 4.422849e+00, 4.921272e+00, 5.379551e+00, 5.829236e+00, 6.325006e+00, 0.000000e+00},/* 
 */{-2.258716e+00, -3.649282e+00, -3.280947e+00, -2.857969e+00, -2.489245e+00, -2.078607e+00, -1.599944e+00, -1.183162e+00, -6.581905e-01, -2.165667e-01, 2.391557e-01, 7.078326e-01, 1.185226e+00, 1.615434e+00, 2.149000e+00, 0.000000e+00},/* 
 */{-4.153297e+00, 3.719297e+00, 4.514258e+00, 5.177461e+00, 5.990068e+00, 6.797599e+00, 7.926182e+00, 9.087798e+00, 1.036058e+01, 1.164706e+01, 1.284663e+01, 1.397711e+01, 1.497106e+01, 1.590517e+01, 1.680064e+01, 0.000000e+00},/* 
 */{-1.088303e+01, -6.487401e+00, -5.009219e+00, -4.476568e+00, -3.870230e+00, -3.486995e+00, -3.019160e+00, -2.624107e+00, -2.258624e+00, -1.902028e+00, -1.591907e+00, -1.281103e+00, -9.788344e-01, -6.913639e-01, -4.364186e-01, 0.000000e+00},/* 
 */{3.618158e+00, -1.329220e+01, -8.101641e+00, -7.408401e+00, -6.714083e+00, -6.246186e+00, -5.821315e+00, -5.475716e+00, -5.136751e+00, -4.809867e+00, -4.542470e+00, -4.269489e+00, -4.039346e+00, -3.766466e+00, -3.597852e+00, 0.000000e+00},/* 
 */{6.869959e+01, -3.406791e+00, -8.466190e-01, -1.191625e-01, 3.858997e-01, 7.331306e-01, 1.068747e+00, 1.350197e+00, 1.644191e+00, 1.887744e+00, 2.151213e+00, 2.432307e+00, 2.656245e+00, 2.804493e+00, 3.183354e+00, 0.000000e+00},/* 
 */{4.410526e+00, 2.504039e+00, 2.973518e+00, 3.331052e+00, 3.797271e+00, 4.747109e+00, 5.209009e+00, 5.603954e+00, 6.007756e+00, 6.392360e+00, 6.700559e+00, 6.970694e+00, 7.297954e+00, 7.466077e+00, 7.932276e+00, 0.000000e+00},/* 
 */{1.109778e+01, -2.048909e+01, -1.800015e+01, -1.738757e+01, -1.689822e+01, -1.656360e+01, -1.625002e+01, -1.595724e+01, -1.567505e+01, -1.539908e+01, -1.510793e+01, -1.483476e+01, -1.459531e+01, -1.440614e+01, -1.432965e+01, 0.000000e+00},/* 
 */{-1.213667e+01, -1.101015e+01, -1.096272e+01, -1.084501e+01, -1.073446e+01, -1.060882e+01, -1.048110e+01, -1.034132e+01, -1.017507e+01, -1.001295e+01, -9.822326e+00, -9.695835e+00, -9.526139e+00, -9.379993e+00, -9.216112e+00, 0.000000e+00},/* 
 */{1.061639e+01, 3.125266e+00, 3.455574e+00, 3.541255e+00, 3.709320e+00, 3.861714e+00, 4.064479e+00, 4.306507e+00, 4.579029e+00, 4.722457e+00, 4.920702e+00, 5.055426e+00, 5.236318e+00, 5.355793e+00, 5.623136e+00, 0.000000e+00},/* 
 */{-1.448791e+01, -5.502174e+00, -5.204192e+00, -4.652717e+00, -4.151636e+00, -3.238154e+00, -2.852045e+00, -2.313122e+00, -1.897296e+00, -1.482405e+00, -1.244742e+00, -8.451148e-01, -5.815543e-01, -3.032072e-01, -1.189653e-02, 0.000000e+00},/* 
 */{-2.707325e+00, -2.584367e+00, -2.409549e+00, -2.284688e+00, -2.133712e+00, -1.959734e+00, -1.790931e+00, -1.567629e+00, -1.323261e+00, -1.052899e+00, -7.708424e-01, -6.108231e-01, -4.304744e-01, -2.479093e-01, -7.504138e-02, 0.000000e+00},/* 
 */{3.218393e-01, 6.194589e-01, 7.327192e-01, 8.538689e-01, 9.833153e-01, 1.122690e+00, 1.226398e+00, 1.417761e+00, 1.557927e+00, 1.798402e+00, 1.987562e+00, 2.272449e+00, 2.409217e+00, 2.623668e+00, 2.785613e+00, 0.000000e+00},/* 
 */{1.277505e+00, 1.027284e+00, 1.221877e+00, 1.338912e+00, 1.493059e+00, 1.674597e+00, 1.846230e+00, 2.075303e+00, 2.316357e+00, 2.548926e+00, 2.714932e+00, 2.884652e+00, 3.050184e+00, 3.163896e+00, 3.391667e+00, 0.000000e+00},/* 
 */{-1.207120e+01, -8.156598e+00, -8.075807e+00, -7.910472e+00, -7.745444e+00, -7.563576e+00, -7.354440e+00, -7.192403e+00, -6.984791e+00, -6.711497e+00, -6.404515e+00, -6.186597e+00, -5.961673e+00, -5.867626e+00, -5.755507e+00, 0.000000e+00},/* 
 */{-3.478837e+00, 3.369522e+00, 3.437537e+00, 3.633164e+00, 3.808192e+00, 4.035446e+00, 4.398657e+00, 4.623378e+00, 4.931857e+00, 5.251209e+00, 5.399956e+00, 5.672224e+00, 5.834579e+00, 5.993530e+00, 6.307445e+00, 0.000000e+00},/* 
 */{-1.363835e+01, -8.990702e+00, -8.610201e+00, -8.279481e+00, -7.804694e+00, -7.504485e+00, -6.471679e+00, -5.562699e+00, -5.145096e+00, -4.614474e+00, -4.228829e+00, -3.849233e+00, -3.496867e+00, -3.249384e+00, -2.905300e+00, 0.000000e+00},/* 
 */{1.011276e+00, -3.992894e+00, -3.755986e+00, -3.468003e+00, -3.286940e+00, -2.935687e+00, -2.579957e+00, -2.117530e+00, -1.744415e+00, -1.480147e+00, -1.229790e+00, -9.693668e-01, -8.101008e-01, -6.026895e-01, -4.826443e-01, 0.000000e+00},/* 
 */{-2.219402e+01, -3.265143e+00, -3.481955e+00, -3.170724e+00, -2.974825e+00, -2.660089e+00, -2.243203e+00, -1.893135e+00, -1.458174e+00, -1.209661e+00, -9.273700e-01, -8.921670e-01, -5.080333e-01, -4.366704e-01, 1.002739e-01, 0.000000e+00}};

  double knotsstartend[29][2]={{-1.710778e+00, 3.721247e+00},/* 
 */{6.197472e+00, 1.148431e+01},/* 
 */{-1.087419e+00, 3.953681e+00},/* 
 */{-1.503606e+00, 3.799374e+00},/* 
 */{-1.512775e+00, 3.843967e+00},/* 
 */{6.092954e-01, 6.079065e+00},/* 
 */{-1.047053e+00, 4.398747e+00},/* 
 */{-2.400468e+00, 3.011205e+00},/* 
 */{-2.242638e+00, 3.328439e+00},/* 
 */{-2.574953e-02, 5.628301e+00},/* 
 */{-2.144638e+00, 3.324762e+00},/* 
 */{-1.843618e+00, 3.555420e+00},/* 
 */{-2.959608e+00, 2.574180e+00},/* 
 */{5.633782e-03, 6.010784e+00},/* 
 */{-1.027839e+00, 3.892822e+00},/* 
 */{-1.087398e+00, 3.956212e+00},/* 
 */{-1.868798e+00, 3.660620e+00},/* 
 */{-1.099687e+00, 4.019378e+00},/* 
 */{-2.601223e+00, 2.429138e+00},/* 
 */{-2.647598e+00, 2.664254e+00},/* 
 */{-1.864454e+00, 3.602390e+00},/* 
 */{-2.526619e+00, 2.607906e+00},/* 
 */{-3.293609e+00, 2.049398e+00},/* 
 */{-2.624522e+00, 2.725601e+00},/* 
 */{-3.205237e+00, 2.210755e+00},/* 
 */{-2.626096e+00, 2.721633e+00},/* 
 */{-1.952496e+00, 2.662643e+00},/* 
 */{-2.372406e+00, 2.989596e+00},/* 
 */{-2.151897e+00, 2.715111e+00}};

  double prepost[29][6]={{-6.243733e-01, -5.225287e+00, 2.283746e+00, 2.634842e+00, 2.075037e+00, 2.865197e+00},/* 
 */{7.254840e+00, -1.137164e+00, 8.103390e-01, 1.042694e+01, 1.618625e+00, 8.085908e-01},/* 
 */{-7.919896e-02, -3.744586e+00, 9.652664e+00, 2.945461e+00, 2.026905e+00, 1.308636e+00},/* 
 */{-4.430099e-01, -7.328762e+00, 1.755362e+00, 2.738778e+00, -1.800972e+00, 1.484694e+00},/* 
 */{-4.414269e-01, -9.283434e+00, 1.976076e+00, 2.772619e+00, -3.867251e+00, 1.774655e+00},/* 
 */{1.703249e+00, -7.025343e+00, 1.369657e+00, 4.985111e+00, -3.584165e+00, 1.052441e+00},/* 
 */{4.210681e-02, -7.386557e+00, 1.383544e+00, 3.309587e+00, -4.429482e+00, 1.381988e+00},/* 
 */{-1.318133e+00, -1.417234e+01, 3.215761e+00, 1.928871e+00, -4.631793e+00, 2.792972e+00},/* 
 */{-1.128423e+00, -1.523605e+01, 2.836680e+00, 2.214224e+00, -4.562079e+00, 2.762995e+00},/* 
 */{1.105061e+00, -5.547702e+00, 3.816863e+00, 4.497491e+00, -1.923361e-01, 1.059374e+00},/* 
 */{-1.050758e+00, 2.714048e-02, 1.173405e+00, 2.230882e+00, 5.923228e+00, 1.197121e+00},/* 
 */{-7.638106e-01, -3.347606e+00, 1.233962e+00, 2.475612e+00, 1.727391e+00, 1.572494e+00},/* 
 */{-1.852851e+00, 4.334585e+00, 3.070713e+00, 1.467422e+00, 1.605737e+01, 1.901853e+00},/* 
 */{1.206664e+00, -5.389431e+00, 3.904192e+00, 4.809754e+00, -6.429083e-01, 8.510484e-01},/* 
 */{-4.370668e-02, -9.494055e+00, 1.277847e+01, 2.908689e+00, -3.737068e+00, 1.044720e+00},/* 
 */{-7.867628e-02, -1.415490e+00, 1.809630e+00, 2.947490e+00, 2.895195e+00, 8.808356e-01},/* 
 */{-7.629147e-01, 2.873238e+00, 1.326832e+00, 2.554736e+00, 7.573555e+00, 6.970094e-01},/* 
 */{-7.587419e-02, -1.861262e+01, 4.479866e+00, 2.995565e+00, -1.438145e+01, 1.475787e+00},/* 
 */{-1.595151e+00, -1.096765e+01, 3.456212e-01, 1.423066e+00, -9.333634e+00, 1.281372e+00},/* 
 */{-1.585227e+00, 3.379066e+00, 3.419473e-01, 1.601884e+00, 5.414026e+00, 3.656166e-01},/* 
 */{-7.710849e-01, -5.253741e+00, 1.921445e+00, 2.509022e+00, -2.442114e-01, 9.786097e-01},/* 
 */{-1.499714e+00, -2.448707e+00, 5.839086e-01, 1.581001e+00, -2.142937e-01, 6.439301e-01},/* 
 */{-2.225008e+00, 7.104946e-01, 4.248487e-01, 9.807968e-01, 2.647647e+00, 4.259992e-01},/* 
 */{-1.554497e+00, 1.176934e+00, 5.713303e-01, 1.655577e+00, 3.214357e+00, 4.221930e-01},/* 
 */{-2.122038e+00, -8.089862e+00, 6.473244e-01, 1.127557e+00, -5.836135e+00, 7.689939e-01},/* 
 */{-1.556550e+00, 3.425783e+00, 8.685352e-01, 1.652087e+00, 6.060708e+00, 4.686255e-01},/* 
 */{-1.029468e+00, -8.696239e+00, 1.852246e+00, 1.739615e+00, -3.168909e+00, 1.500176e+00},/* 
 */{-1.300006e+00, -3.792021e+00, 5.259576e-01, 1.917196e+00, -5.856522e-01, 5.501364e-01},/* 
 */{-1.178495e+00, -3.419060e+00, 1.308079e+00, 1.741709e+00, -2.904015e-01, 1.464774e+00}};      /*this is prepostCalphaorg*/

  double ss;
  int a;
  if (measurenumber>28){
    printf("WARNING measure number is to high in Calpha original\n");
    return -99999999999.9 ;
  }
  if (measurevalue<prepost[measurenumber][0]){
    return prepost[measurenumber][2]*(measurevalue-prepost[measurenumber][0])+prepost[measurenumber][1];
  }
  if (measurevalue>prepost[measurenumber][3]){
    return prepost[measurenumber][5]*(measurevalue-prepost[measurenumber][3])+prepost[measurenumber][4];
  }
  ss=1+19*(measurevalue-knotsstartend[measurenumber][0])/(knotsstartend[measurenumber][1]-knotsstartend[measurenumber][0]);
  a=floor(ss);
  if(a>16) {
    printf("wRNING PROBLEM IN PREDEFORMATION AS a>16\n");
    return 0;
    }
  
  return splinefun(ss-a+3)*CoefCalphaOrg[measurenumber][a-3-1]+splinefun(ss-a+2)*CoefCalphaOrg[measurenumber][a-2-1]+splinefun(ss-a+1)*CoefCalphaOrg[measurenumber][a-1-1]+splinefun(ss-a)*CoefCalphaOrg[measurenumber][a-1];
  
}

double presplinedeformationCalphanew(double measurevalue, int measurenumber)
{
  double CoefCalphaNew[29][16]={{-1.022043e+01, -5.673251e+00, -5.109977e+00, -4.515572e+00, -3.920791e+00, -3.315927e+00, -2.694542e+00, -2.072482e+00, -1.459846e+00, -8.388197e-01, -2.126423e-01, 4.441569e-01, 1.154193e+00, 1.859077e+00, 2.844650e+00, 0.000000e+00},/* 
 */{-6.945792e+00, -4.441694e-01, -2.081807e-01, 4.892523e-02, 3.033886e-01, 5.427508e-01, 7.881182e-01, 1.041829e+00, 1.290606e+00, 1.539794e+00, 1.795005e+00, 2.041712e+00, 2.288341e+00, 2.543794e+00, 2.812216e+00, 0.000000e+00},/* 
 */{-3.808358e+01, -1.452142e+01, -1.416696e+01, -1.367478e+01, -1.305820e+01, -1.246337e+01, -1.183278e+01, -1.135032e+01, -1.091899e+01, -1.052440e+01, -1.017863e+01, -9.747768e+00, -9.429296e+00, -8.999715e+00, -8.824979e+00, 0.000000e+00},/* 
 */{-2.667715e+00, -2.359095e+00, -2.112318e+00, -1.854287e+00, -1.554363e+00, -1.280617e+00, -9.654002e-01, -6.694114e-01, -3.721960e-01, -8.002989e-02, 2.655362e-01, 5.351787e-01, 8.459453e-01, 1.150311e+00, 1.484096e+00, 0.000000e+00},/* 
 */{-2.202387e+00, -6.054824e+00, -5.622242e+00, -5.384090e+00, -5.119058e+00, -4.852244e+00, -4.564886e+00, -4.268020e+00, -3.995904e+00, -3.693291e+00, -3.433267e+00, -3.118993e+00, -2.802019e+00, -2.479632e+00, -2.247915e+00, 0.000000e+00},/* 
 */{-1.027688e+01, -5.503177e+00, -5.083463e+00, -4.633020e+00, -4.196990e+00, -3.768040e+00, -3.317712e+00, -2.906712e+00, -2.504089e+00, -2.086427e+00, -1.681782e+00, -1.270625e+00, -8.264659e-01, -4.492907e-01, -5.198454e-02, 0.000000e+00},/* 
 */{-6.223286e+00, -6.755129e+00, -6.593793e+00, -6.431788e+00, -6.242316e+00, -6.025544e+00, -5.807117e+00, -5.581122e+00, -5.354590e+00, -5.122051e+00, -4.862186e+00, -4.623467e+00, -4.433281e+00, -4.227595e+00, -4.108058e+00, 0.000000e+00},/* 
 */{-2.199577e+01, -1.202079e+01, -1.191170e+01, -1.134749e+01, -1.093022e+01, -1.043141e+01, -9.918041e+00, -9.382728e+00, -8.777832e+00, -8.231580e+00, -7.579214e+00, -7.039768e+00, -6.366825e+00, -5.793309e+00, -5.222663e+00, 0.000000e+00},/* 
 */{-4.132907e+00, -8.980722e+00, -8.483384e+00, -8.022353e+00, -7.519826e+00, -6.965669e+00, -6.325899e+00, -5.769631e+00, -5.069582e+00, -4.481171e+00, -3.777932e+00, -3.293969e+00, -2.725822e+00, -2.285528e+00, -1.648954e+00, 0.000000e+00},/* 
 */{-4.375374e+00, -3.637166e+00, -3.260730e+00, -2.888041e+00, -2.471268e+00, -2.063043e+00, -1.671798e+00, -1.288467e+00, -8.985132e-01, -4.971453e-01, -8.729825e-02, 2.917771e-01, 6.360047e-01, 1.070199e+00, 1.431675e+00, 0.000000e+00},/* 
 */{5.587004e+00, 5.069940e+00, 5.566110e+00, 5.923896e+00, 6.310691e+00, 6.757923e+00, 7.188497e+00, 7.650106e+00, 8.149477e+00, 8.729015e+00, 9.298448e+00, 9.769806e+00, 1.025161e+01, 1.067721e+01, 1.125619e+01, 0.000000e+00},/* 
 */{-5.415667e+00, -5.001212e+00, -4.587147e+00, -4.216132e+00, -3.817904e+00, -3.451488e+00, -3.008831e+00, -2.618388e+00, -2.203575e+00, -1.806360e+00, -1.396439e+00, -1.026924e+00, -5.885159e-01, -2.406653e-01, 2.305789e-01, 0.000000e+00},/* 
 */{-1.611904e+01, -7.715713e+00, -7.379839e+00, -6.883516e+00, -6.413935e+00, -5.854990e+00, -5.371484e+00, -4.774634e+00, -4.168858e+00, -3.563709e+00, -2.944724e+00, -2.328611e+00, -1.670275e+00, -1.002568e+00, -5.702408e-01, 0.000000e+00},/* 
 */{-1.444988e+01, -7.489050e+00, -7.181202e+00, -6.813808e+00, -6.468570e+00, -6.049876e+00, -5.669575e+00, -5.285868e+00, -4.927988e+00, -4.547121e+00, -4.168992e+00, -3.785502e+00, -3.429022e+00, -3.059825e+00, -2.721411e+00, 0.000000e+00},/* 
 */{-1.073081e+01, -9.210469e+00, -8.853406e+00, -8.738926e+00, -8.582025e+00, -8.315546e+00, -8.157126e+00, -7.725160e+00, -7.347063e+00, -6.737670e+00, -6.454792e+00, -6.139445e+00, -5.924127e+00, -5.792245e+00, -5.544614e+00, 0.000000e+00},/* 
 */{-6.592760e+00, 7.031428e-01, 6.665296e-01, 9.171217e-01, 1.017995e+00, 1.180821e+00, 1.405540e+00, 1.550663e+00, 1.708191e+00, 1.823997e+00, 1.968231e+00, 2.078781e+00, 2.218312e+00, 2.335107e+00, 2.510290e+00, 0.000000e+00},/* 
 */{1.563435e+00, -5.825979e-01, -2.251505e-01, -1.288077e-02, 3.098490e-01, 5.928416e-01, 9.729463e-01, 1.323093e+00, 2.040682e+00, 2.579530e+00, 2.934737e+00, 3.272113e+00, 3.603578e+00, 3.815528e+00, 4.173744e+00, 0.000000e+00},/* 
 */{-6.113174e-01, -7.846237e-01, -7.597642e-01, -6.161209e-01, -5.233465e-01, -3.777681e-01, -2.315315e-01, -4.992162e-02, 1.312869e-01, 2.538269e-01, 3.949310e-01, 4.953781e-01, 6.256746e-01, 7.234966e-01, 8.663756e-01, 0.000000e+00},/* 
 */{-1.268815e+00, -5.214119e-01, -4.949137e-01, -3.756819e-01, -2.746838e-01, -1.375973e-01, -1.410390e-02, 1.345969e-01, 2.762125e-01, 4.223013e-01, 5.667723e-01, 7.087235e-01, 8.672988e-01, 9.566240e-01, 1.089657e+00, 0.000000e+00},/* 
 */{-2.713323e+00, -2.835140e+00, -2.672170e+00, -2.545607e+00, -2.405050e+00, -2.284063e+00, -2.118973e+00, -1.941636e+00, -1.745351e+00, -1.510684e+00, -1.285044e+00, -1.134234e+00, -9.248740e-01, -7.652693e-01, -6.376662e-01, 0.000000e+00},/* 
 */{-7.049540e+00, -1.230206e+01, -1.192029e+01, -1.161157e+01, -1.129534e+01, -1.092185e+01, -1.059090e+01, -9.642770e+00, -9.077614e+00, -8.660075e+00, -8.253944e+00, -7.921857e+00, -7.670662e+00, -7.376979e+00, -7.261573e+00, 0.000000e+00},/* 
 */{3.554987e+00, 3.860075e+00, 3.960303e+00, 4.081331e+00, 4.227219e+00, 4.351970e+00, 4.564150e+00, 4.730513e+00, 4.927617e+00, 5.230843e+00, 5.471028e+00, 5.635116e+00, 5.789697e+00, 5.946798e+00, 6.207319e+00, 0.000000e+00},/* 
 */{-5.630399e+00, -6.700020e+00, -6.598295e+00, -6.455081e+00, -6.338798e+00, -6.240210e+00, -6.062392e+00, -5.913847e+00, -5.681412e+00, -5.467231e+00, -5.293321e+00, -5.087409e+00, -4.934734e+00, -4.731242e+00, -4.697878e+00, 0.000000e+00},/* 
 */{-1.559733e+01, -4.735866e+00, -4.585365e+00, -4.469057e+00, -4.276035e+00, -4.117168e+00, -3.968017e+00, -3.799680e+00, -3.625358e+00, -3.407268e+00, -3.161770e+00, -2.907989e+00, -2.681693e+00, -2.558691e+00, -2.395861e+00, 0.000000e+00},/* 
 */{-6.617084e+00, -6.595266e+00, -6.416685e+00, -6.302219e+00, -6.136272e+00, -5.952989e+00, -5.816858e+00, -5.613805e+00, -5.396605e+00, -5.125987e+00, -4.921291e+00, -4.651057e+00, -4.536098e+00, -4.410255e+00, -4.200978e+00, 0.000000e+00},/* 
 */{-6.273604e+00, -4.171191e+00, -4.051720e+00, -3.888250e+00, -3.701981e+00, -3.441993e+00, -3.161095e+00, -2.896528e+00, -2.588697e+00, -2.309606e+00, -2.128251e+00, -1.897308e+00, -1.688176e+00, -1.522667e+00, -1.327367e+00, 0.000000e+00},/* 
 */{-1.146375e+01, -1.228110e+01, -1.221622e+01, -1.180033e+01, -1.155654e+01, -1.121121e+01, -1.082499e+01, -1.054093e+01, -9.878760e+00, -9.434193e+00, -9.041862e+00, -8.749843e+00, -8.431162e+00, -8.192864e+00, -7.977172e+00, 0.000000e+00},/* 
 */{9.868993e+00, 7.188876e+00, 7.403959e+00, 7.483963e+00, 7.655147e+00, 7.842711e+00, 8.034760e+00, 8.227436e+00, 8.500959e+00, 8.831192e+00, 9.072816e+00, 9.219329e+00, 9.464403e+00, 9.551125e+00, 9.842555e+00, 0.000000e+00},/* 
 */{-1.300258e+01, -8.143445e+00, -8.049160e+00, -7.807858e+00, -7.621204e+00, -7.381452e+00, -7.039130e+00, -6.704561e+00, -6.410524e+00, -6.024031e+00, -5.887703e+00, -5.584039e+00, -5.538617e+00, -5.292546e+00, -5.189083e+00, 0.000000e+00}};

  double knotsstartend[29][2]={{-1.710778e+00, 3.721247e+00},/* 
 */{-2.762178e+00, 2.572634e+00},/* 
 */{-2.175692e+00, 3.299203e+00},/* 
 */{-2.444458e+00, 3.188739e+00},/* 
 */{-2.499671e+00, 3.149118e+00},/* 
 */{-2.063703e+00, 3.268983e+00},/* 
 */{-1.858791e+00, 3.542937e+00},/* 
 */{-3.335414e+00, 2.165310e+00},/* 
 */{-3.136738e+00, 2.445955e+00},/* 
 */{-3.736891e+00, 1.764543e+00},/* 
 */{-3.305861e+00, 2.320500e+00},/* 
 */{-3.116074e+00, 2.046173e+00},/* 
 */{-3.966153e+00, 1.646521e+00},/* 
 */{-3.429630e+00, 1.796802e+00},/* 
 */{-2.695726e+00, 2.345148e+00},/* 
 */{-2.105666e+00, 3.240836e+00},/* 
 */{-2.736646e+00, 2.615384e+00},/* 
 */{-2.363839e+00, 2.798067e+00},/* 
 */{-2.867205e+00, 1.909386e+00},/* 
 */{-3.126988e+00, 2.294584e+00},/* 
 */{-2.535714e+00, 2.984478e+00},/* 
 */{-2.342632e+00, 2.634388e+00},/* 
 */{-3.027273e+00, 2.367289e+00},/* 
 */{-3.346432e+00, 1.934190e+00},/* 
 */{-2.948940e+00, 2.303805e+00},/* 
 */{-2.620212e+00, 2.694336e+00},/* 
 */{-2.289288e+00, 2.245806e+00},/* 
 */{-2.893821e+00, 2.440978e+00},/* 
 */{-2.172856e+00, 2.506878e+00}};

  double prepost[29][6]={{-6.243733e-01, -5.225287e+00, 2.283746e+00, 2.634842e+00, 2.075037e+00, 2.865197e+00},/* 
 */{-1.695216e+00, -2.619306e-01, 1.310846e+00, 1.505672e+00, 2.594478e+00, 7.217729e-01},/* 
 */{-1.080713e+00, -1.425704e+01, 2.993852e+00, 2.204224e+00, -8.974981e+00, 1.489806e+00},/* 
 */{-1.317819e+00, -2.160795e+00, 8.486626e-01, 2.062100e+00, 1.217154e+00, 9.714297e-01},/* 
 */{-1.369913e+00, -5.719636e+00, 9.574822e-01, 2.019360e+00, -2.438337e+00, 1.012615e+00},/* 
 */{-9.971658e-01, -5.170589e+00, 1.840708e+00, 2.202446e+00, -3.685721e-01, 1.368016e+00},/* 
 */{-7.784450e-01, -6.625079e+00, 5.194738e-01, 2.462591e+00, -4.205721e+00, 7.980071e-01},/* 
 */{-2.235269e+00, -1.190784e+01, 1.561394e+00, 1.065165e+00, -5.673222e+00, 2.295605e+00},/* 
 */{-2.020199e+00, -8.578823e+00, 1.289255e+00, 1.329417e+00, -2.140114e+00, 2.021647e+00},/* 
 */{-2.636604e+00, -3.336819e+00, 1.320924e+00, 6.642562e-01, 1.133898e+00, 1.204916e+00},/* 
 */{-2.180589e+00, 5.456418e+00, 1.457571e+00, 1.195228e+00, 1.079031e+01, 9.901023e-01},/* 
 */{-2.083624e+00, -4.673634e+00, 1.473320e+00, 1.013723e+00, -1.368226e-01, 1.537454e+00},/* 
 */{-2.843618e+00, -7.444079e+00, 1.857005e+00, 5.239864e-01, -9.360043e-01, 1.727827e+00},/* 
 */{-2.384344e+00, -7.246561e+00, 1.672131e+00, 7.515155e-01, -2.991592e+00, 1.439331e+00},/* 
 */{-1.687551e+00, -8.947070e+00, 1.140940e+00, 1.336973e+00, -5.725779e+00, 1.193065e+00},/* 
 */{-1.036365e+00, 6.885838e-01, 7.176493e-01, 2.171536e+00, 2.371545e+00, 3.652866e-01},/* 
 */{-1.666240e+00, -3.056905e-01, 9.262838e-01, 1.544978e+00, 3.893610e+00, 7.837513e-01},/* 
 */{-1.331458e+00, -7.543357e-01, 2.168256e-01, 1.765686e+00, 7.545716e-01, 3.985423e-01},/* 
 */{-1.911887e+00, -4.932613e-01, 2.807923e-01, 9.540676e-01, 9.853301e-01, 3.762646e-01},/* 
 */{-2.042673e+00, -2.707491e+00, 5.103429e-01, 1.210270e+00, -7.417994e-01, 5.188261e-01},/* 
 */{-1.431675e+00, -1.199537e+01, 8.456939e-01, 1.880440e+00, -7.359583e+00, 1.085501e+00},/* 
 */{-1.347228e+00, 3.941759e+00, 4.236768e-01, 1.638984e+00, 5.999103e+00, 3.743844e-01},/* 
 */{-1.948360e+00, -6.613538e+00, 3.225326e-01, 1.288377e+00, -4.732868e+00, 6.378293e-01},/* 
 */{-2.290308e+00, -4.632664e+00, 1.272917e+00, 8.780654e-01, -2.519749e+00, 7.007070e-01},/* 
 */{-1.898391e+00, -6.457663e+00, 5.604023e-01, 1.253256e+00, -4.355958e+00, 9.491859e-01},/* 
 */{-1.557303e+00, -4.074504e+00, 6.192396e-01, 1.631426e+00, -1.479556e+00, 7.450802e-01},/* 
 */{-1.382269e+00, -1.219807e+01, 6.684800e-01, 1.338787e+00, -8.141306e+00, 1.584302e+00},/* 
 */{-1.826861e+00, 7.353276e+00, 4.058494e-01, 1.374019e+00, 9.613367e+00, 8.278267e-02},/* 
 */{-1.236909e+00, -8.061824e+00, 9.607217e-01, 1.570931e+00, -5.277242e+00, 1.018303e+00}};      /*this is prepostCalphaoNew*/

  double ss;
  int a;
  if (measurenumber>28){
    printf("WARNING measure number is to high in Calpha New\n");
    return -99999999999.9 ;
  }
  if (measurevalue<prepost[measurenumber][0]){
    return prepost[measurenumber][2]*(measurevalue-prepost[measurenumber][0])+prepost[measurenumber][1];
  }
  if (measurevalue>prepost[measurenumber][3]){
    return prepost[measurenumber][5]*(measurevalue-prepost[measurenumber][3])+prepost[measurenumber][4];
  }
  ss=1+19*(measurevalue-knotsstartend[measurenumber][0])/(knotsstartend[measurenumber][1]-knotsstartend[measurenumber][0]);
  a=floor(ss);
  if(a>16) {
    printf("wRNING PROBLEM IN PREDEFORMATION AS a>16\n");
    return 0;
    }
  
  return splinefun(ss-a+3)*CoefCalphaNew[measurenumber][a-3-1]+splinefun(ss-a+2)*CoefCalphaNew[measurenumber][a-2-1]+splinefun(ss-a+1)*CoefCalphaNew[measurenumber][a-1-1]+splinefun(ss-a)*CoefCalphaNew[measurenumber][a-1];
}

double splinedeformationCalphaorg(double measurevalue, int measurenumber)
{
  return presplinedeformationCalphaorg(measurevalue, measurenumber)-presplinedeformationCalphaorg(0.0, measurenumber);
}
double splinedeformationCalphanew(double measurevalue, int measurenumber)
{
  return presplinedeformationCalphanew(measurevalue, measurenumber)-presplinedeformationCalphanew(0.0, measurenumber);
}








/**********************PRINTING*******************/

void GaussIntegrals(void)
{
  extern int nres;
  extern double polyl[MAXLENGTH][3];
  extern double partsum[MAXLENGTH][MAXLENGTH];
  double invar[29], invar2[29], R1,R12, R2,R22, smoothbacbonelength,smoothbacbonelength_back, v1[3];
  double gaussint[29], gaussint2[29];

  double scales2order[29]={0.,1.0059, 0.6025, 0.7274, 0.7338, 1.0524, 13.6944, 4.0666, 3.4855, 0.8277, 3.8296, 0.6918, 2.0215, 0.9243, 0.3171, 6.7371, 3.4645, 6.7988, -1.8539, 23.2397, 3.2125, -6.7321, -60.2435, 11.6135, -39.5212, -0.3665, 2.6708, 39.9532, 1.3327};

  double scales2ordersmooth[29]={ 0.000000e+00, 1.039329e+00, 9.771578e-03, -2.624606e+00, -2.761606e+00, 4.055914e+02, 1.744519e-03, -3.143907e-01, -2.198699e-01, 1.384338e+02, 4.572647e-03, -2.832035e+00, -6.652118e-01, 2.216269e+02, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00};

  double fitcoef3order[4][30]={{0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 3.300254e-01, -6.445656e-03, -3.663870e-02, -4.704765e-03, 2.035057e-05, -6.936558e-04, 2.484670e-02, -3.957696e-06, 2.436689e-05, 6.354898e-04, 3.108727e-05, -1.663249e-04, -2.461276e-03, -7.235432e-05, -4.429750e-04},{0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, -1.445703e+00, 3.358944e-01, 1.082980e-01, 3.762040e-01, 9.632395e-03, 1.817898e-02, -1.527511e-01, -2.477205e-03, 6.961623e-04, 6.263854e-02, 2.203401e-04, -2.872356e-03, -8.179391e-02, 8.723373e-03, -2.767734e-03},{0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, -4.552952e-01, -1.327581e-02, 1.842301e-01, -1.113365e-02, -1.734630e-04, -5.123921e-04, 2.695453e-01, -2.744619e-04, -8.456795e-04, -1.928650e-03, -1.131743e-03, 7.353213e-04, 2.148590e-01, 3.772707e-03, 9.536288e-03},{0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, -2.389632e-02, 3.705009e-01, 1.061140e+00, 2.680487e-01, -1.058196e-02, 4.529179e-02, -5.493136e-01, 5.274348e-03, -3.910647e-03, -5.922653e-02, -3.572589e-03, 1.368161e-02, -5.621730e-03, -1.940846e-03, 1.820422e-02}};

  double fitcoef3ordersmooth[4][31]={{0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 2.891192e-01, 6.776607e-03, 2.171263e-02, 7.035961e-03, -1.396743e-03, -3.727756e-04, 2.469942e-02, -3.696977e-04, -4.516382e-04, 1.320793e-03, 9.492686e-05, -3.673666e-03, -3.862823e-03, -2.694790e-03, -5.967507e-03},/* 
 */{0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, -4.673419e-01, 2.628907e-01, 7.915791e-02, 2.651258e-01, 5.228681e-02, 5.545261e-02, -8.421819e-04, 2.789268e-02, 1.169106e-02, 1.844081e-02, 7.532571e-03, 4.739005e-02, -1.475045e-01, 4.220526e-02, -3.022622e-02},/* 
 */{0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, -1.448568e-01, -1.981682e-02, 3.805253e-02, -1.624695e-02, -1.448439e-03, 5.549353e-03, 2.975958e-02, -2.610117e-03, 3.983089e-04, 1.113884e-02, -3.345207e-03, 3.354085e-03, 2.547726e-01, 1.331286e-02, 8.324710e-02},/* 
 */{0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, -2.036573e-05, 1.887763e-05, 2.676509e-05, -2.694326e-06, -2.417911e-05, 5.769087e-05, -1.331947e-05, 2.895806e-07, -5.701623e-06, -1.509898e-06, -1.363764e-05, -5.481478e-05, -5.027138e-05, -2.042450e-05, 4.670799e-05}};

  double scalelawsorg[29]={.92, 1.37, 1.89, 2.26, 2.26, 2.68, 0.90, 1.64, 1.76, 2.56, 1.32, 2.43, 1.80, 2.76, 2.90, 1.88, 2.65, 1.86, 1.69, 2.09, 2.81, 1.24, 1.38, 2.19, 2.14, 1.8, 2.04, 1.19, 2.24}; 
  
  double scalelawsnew[29]={0.92, 1.00, 1.01, 1.92, 2.03, 2.28, 0.98, 1.68, 1.69, 2.17, 1.48, 1.93, 1.99, 2.36, 1.53, 1.96, 2.11, 2.07, 1.69, 2.14, 1.86, 1.19, 1.23, 2.09, 2.10, 1.83, 2.16, 0.96, 2.23};
  
  double scalelawssmoothorg[30]={0.98, 0.79, 1.28, 1.68, 2.23, 2.25, 2.87, 1.29, 1.53, 1.63, 2.28, 1.36, 2.11, 1.87, 2.55, 2.74, 2.02, 1.99, 2.24, 1.76, 1.50, 2.24, 1.15, 1.15, 1.66, 1.94, 1.74, 2.42, 1.78, 1.80};

  double scalelawssmoothnew[30]={0.92, 0.79, 1.01, 1.56, 1.88, 1.92, 1.50, 1.24, 1.34, 1.48, 0.93, 1.27, 1.47, 1.76, 0.94, 1.37, 1.45, 1.51, 2.06, 1.57, 1.44, 2.01, 1.09, 0.99, 1.64, 1.91, 1.76, 1.88, 1.75, 1.72};

  double variansernew[30]={1.015899e+00, 6.046441e-02, 3.494607e-02, 3.406224e-02, 2.585240e-03, 1.441085e-03, 1.243288e-03, 4.466983e-03, 1.570046e-03, 1.685361e-03, 1.175530e-03, 2.416459e-03, 1.782516e-03, 1.133752e-03, 4.762324e-04, 4.953829e-03, 5.593157e-05, 1.956618e-04, 3.128713e-05, 2.600544e-05, 1.784686e-05, 7.642103e-04, 2.781596e-04, 2.617582e-04, 2.663494e-05, 4.370924e-06, 8.149216e-05, 1.113821e-04, 2.377183e-03, 1.769148e-05};
  
  double varianserorg[30]={1.015899e+00, 6.046441e-02, 5.884097e-03, 2.993718e-03, 1.502210e-03, 1.459503e-03, 2.992914e-04, 8.094224e-03, 2.725560e-03, 1.716925e-03, 1.982405e-04, 5.378608e-03, 2.947282e-04, 2.712258e-03, 1.039000e-04, 8.584549e-05, 4.481765e-04, 1.760766e-05, 4.933158e-04, 2.891304e-05, 2.412347e-05, 9.506914e-06, 2.199716e-04, 1.306271e-04, 1.597722e-05, 3.813277e-06, 9.418928e-05, 4.578298e-04, 1.168567e-03, 1.948569e-05};
  
  double vs[30]={2.178927e+02, 2.632709e+00, 1.896596e+01, 7.482618e+00, 5.122873e+01, 5.539920e+01, 4.408378e+02, 6.640930e-01, 8.953117e+00, 1.013272e+01, 2.048619e+02, 1.840360e+00, 3.586241e+01, 1.570184e+01, 2.857427e+02, 2.037704e+01, 1.966589e+00, 3.038134e+00, 2.048440e+00, 2.777770e-01, 5.623737e-01, 2.804170e+00, 1.716416e-01, 1.191713e-01, 4.936448e-01, 1.773537e-01, 4.841054e-01, 3.842461e+00, 3.365983e-01, 1.054532e+00};  /*variance_smooth_CATH2.4_Hclasses*/

  double backpolyl[MAXLENGTH][3], psum[3][3],crossvec[3];
  double backcenterecdpolyl[MAXLENGTH][3]; 
  double invarlength, invarlength2;
  double invarsmoothorg[30],invarsmoothnew[30];
  double invarCalphaorg[29],invarCalphanew[29];
  double distance, vcenter_of_mass[3], center_of_mass[3], v_centred[3];
  double A1[3],A2[3],A3[3],B[3],det,det1,det2,det3,rot_vec[3];
  int i,j,k,ii;
  

  if(dipeptide==1){print_dipeptide_frequeces();}
  
  
  if(SmoothenBackBone==1){
    PerturbBackBone(pow(10,-9));
  }
  for(j=1;j<=nres;j++){
    for(i=0;i<3;i++){
      backpolyl[j][i]=polyl[j][i];
    }
  }
  if(SmoothenBackBone==1){
    SmoothenTheBackbone();
    smoothbacbonelength=0;
    for(i=1;i<nres;i++){
      geovec(v1,polyl[i],polyl[i+1]);
      smoothbacbonelength+=length(v1);
    }
    smoothbacbonelength_back=smoothbacbonelength;
  }
  
  createunitvectors();
  createomega();
  createpartsums();
  
  
  gaussint[0]=int12();
  gaussint[1]=inta12();
  gaussint[2]=int12_34();
  gaussint[3]=inta12_34();
  gaussint[4]=int12_a34();
  gaussint[5]=inta12_a34();
  
  gaussint[6]=int13_24();
  gaussint[7]=inta13_24();
  gaussint[8]=int13_a24();
  gaussint[9]=inta13_a24();
  
  gaussint[10]=int14_23();
  gaussint[11]=inta14_23();
  gaussint[12]=int14_a23();
  gaussint[13]=inta14_a23();
  
  gaussint[14]=int12_34_56();
  gaussint[15]=int12_35_46();
  gaussint[16]=int12_36_45();
  gaussint[17]=int13_24_56();
  gaussint[18]=int13_25_46();
  gaussint[19]=int13_26_45();
  gaussint[20]=int14_23_56();
  gaussint[21]=int14_25_36();
  gaussint[22]=int14_26_35();
  gaussint[23]=int15_23_46();
  gaussint[24]=int15_24_36();
  gaussint[25]=int15_26_34();
  gaussint[26]=int16_23_45();
  gaussint[27]=int16_24_35();
  gaussint[28]=int16_25_34();
  
  if(SmoothenBackBone==1){
    invarlength=smoothbacbonelength-2.2539*nres;
    invar[0]=gaussint[0];
    R1=gaussint[0]/refGauss[nres][0];
    R2=gaussint[1]/refGauss[nres][1];
    
    invar[1]=gaussint[1]-(1/scales2ordersmooth[1])*refGauss[nres][1]; 
    invar[2]=gaussint[2]-(1/scales2ordersmooth[2])*refGauss[nres][2]*R1*R1;
    invar[3]=gaussint[3]-(1/scales2ordersmooth[3])*refGauss[nres][3]*R2*R1;
    invar[4]=gaussint[4]-(1/scales2ordersmooth[4])*refGauss[nres][4]*R1*R2;
    invar[5]=gaussint[5]-(1/scales2ordersmooth[5])*refGauss[nres][5]*R2*R2;
    
    invar[6]=gaussint[6]-(1/scales2ordersmooth[6])*refGauss[nres][6]*R1*R1;
    invar[7]=gaussint[7]-(1/scales2ordersmooth[7])*refGauss[nres][7]*R2*R1;
    invar[8]=gaussint[8]-(1/scales2ordersmooth[8])*refGauss[nres][8]*R1*R2;
    invar[9]=gaussint[9]-(1/scales2ordersmooth[9])*refGauss[nres][9]*R2*R2;
    
    invar[10]=gaussint[10]-(1/scales2ordersmooth[10])*refGauss[nres][10]*R1*R1;
    invar[11]=gaussint[11]-(1/scales2ordersmooth[11])*refGauss[nres][11]*R2*R1;
    invar[12]=gaussint[12]-(1/scales2ordersmooth[12])*refGauss[nres][12]*R1*R2;
    invar[13]=gaussint[13]-(1/scales2ordersmooth[13])*refGauss[nres][13]*R2*R2;
    
    for(i=14;i<29;i++)
      {
	invar[i]=gaussint[i]-fitcoef3ordersmooth[0][i+2]*gaussint[0]*gaussint[2]-fitcoef3ordersmooth[1][i+2]*gaussint[0]*gaussint[6]-fitcoef3ordersmooth[2][i+2]*gaussint[0]*gaussint[10]-fitcoef3ordersmooth[3][i+2]*refGauss[nres][i]*R1*R1*R1;
      }
    
  }
  else{
    
    invar[0]=gaussint[0];
    R1=gaussint[0]/refGauss[nres][0];
    R2=gaussint[1]/refGauss[nres][1];
    
    invar[1]=gaussint[1]-(1/scales2order[1])*refGauss[nres][1]; 
    invar[2]=gaussint[2]-(1/scales2order[2])*refGauss[nres][2]*R1*R1;
    invar[3]=gaussint[3]-(1/scales2order[3])*refGauss[nres][3]*R2*R1;
    invar[4]=gaussint[4]-(1/scales2order[4])*refGauss[nres][4]*R1*R2;
    invar[5]=gaussint[5]-(1/scales2order[5])*refGauss[nres][5]*R2*R2;
    
    invar[6]=gaussint[6]-(1/scales2order[6])*refGauss[nres][6]*R1*R1;
    invar[7]=gaussint[7]-(1/scales2order[7])*refGauss[nres][7]*R2*R1;
    invar[8]=gaussint[8]-(1/scales2order[8])*refGauss[nres][8]*R1*R2;
    invar[9]=gaussint[9]-(1/scales2order[9])*refGauss[nres][9]*R2*R2;
    
    invar[10]=gaussint[10]-(1/scales2order[10])*refGauss[nres][10]*R1*R1;
    invar[11]=gaussint[11]-(1/scales2order[11])*refGauss[nres][11]*R2*R1;
    invar[12]=gaussint[12]-(1/scales2order[12])*refGauss[nres][12]*R1*R2;
    invar[13]=gaussint[13]-(1/scales2order[13])*refGauss[nres][13]*R2*R2;
    
    for(i=14;i<29;i++)
      {
	invar[i]=gaussint[i]-fitcoef3order[0][i+1]*gaussint[0]*gaussint[2]-fitcoef3order[1][i+1]*gaussint[0]*gaussint[6]-fitcoef3order[2][i+1]*gaussint[0]*gaussint[10]-fitcoef3order[3][i+1]*refGauss[nres][i]*R1*R1*R1;
      }
  }

  
  
  if(SmoothenBackBone==1){
    invarsmoothorg[0]=splinedeformationsmoothorg(smoothbacbonelength_back/pow(nres,scalelawssmoothorg[0]),0);
    for(i=0;i<29;i++){
      invarsmoothorg[i+1]=splinedeformationsmoothorg(gaussint[i]/pow(nres,scalelawssmoothorg[i+1]),i+1);
    }
    invarsmoothnew[0]=splinedeformationsmoothnew(invarlength/pow(nres,scalelawssmoothnew[0]),0);
    for(i=0;i<29;i++){
      invarsmoothnew[i+1]=splinedeformationsmoothnew(invar[i]/pow(nres,scalelawssmoothnew[i+1]),i+1);
    } 
    
    /* PRINTING ALL INVARIANTS
    fprintf(farrayo,"%i %5.4e",nres, 19.11*pow(nres,1.0/3.0));
    for(i=0;i<30;i++){fprintf(farrayo," %5.4e",invarsmoothorg[i]);}
    for(i=0;i<30;i++){fprintf(farrayo," %5.4e",invarsmoothnew[i]);}
    fprintf(farrayo,"\n");
    */

    /* PRINTING ACCORDING TO PRIORITY  */
    fprintf(farrayo,"%i %5.4e",nres, 19.11*pow(nres,1.0/3.0));
    fprintf(farrayo," %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e\n",invarsmoothorg[13],invarsmoothorg[9],invarsmoothorg[10],invarsmoothnew[2],invarsmoothorg[14],invarsmoothorg[8],invarsmoothorg[1],invarsmoothorg[6],invarsmoothorg[12],invarsmoothnew[0],invarsmoothorg[5],invarsmoothorg[4],invarsmoothorg[11],invarsmoothnew[27],invarsmoothorg[3],invarsmoothorg[29],invarsmoothnew[15],invarsmoothorg[7],invarsmoothorg[26],invarsmoothnew[17],invarsmoothnew[21],invarsmoothorg[24],invarsmoothorg[25],invarsmoothorg[20],invarsmoothnew[18],invarsmoothorg[28],invarsmoothnew[16],invarsmoothorg[22],invarsmoothorg[23],invarsmoothorg[19]);
  }
  else
    {
      for(i=0;i<29;i++){
      invarCalphaorg[i]=splinedeformationCalphaorg(gaussint[i]/pow(nres,scalelawsorg[i])/varianserorg[i+1],i);
      }
      for(i=0;i<29;i++){
      invarCalphanew[i]=splinedeformationCalphanew(invar[i]/pow(nres,scalelawsnew[i])/variansernew[i+1],i);
      }

      /*PRINTING ALL
      fprintf(farrayo,"%i %5.4e",nres, pow(nres,1.0/3.0)/varianserorg[0]);
      for(i=0;i<29;i++){fprintf(farrayo," %5.4e",invarCalphaorg[i]);}
      for(i=0;i<29;i++){fprintf(farrayo," %5.4e",invarCalphanew[i]);}
      fprintf(farrayo,"\n");
      */
      
      /* PRINTING ACCORDING TO PRIORITY  */
      fprintf(farrayo,"%i %5.4e",nres, pow(nres,1.0/3.0)/varianserorg[0]);
      fprintf(farrayo," %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e\n",invarCalphaorg[12],invarCalphaorg[8],invarCalphaorg[7],invarCalphaorg[0],invarCalphaorg[10],invarCalphanew[5],invarCalphanew[2],invarCalphanew[9],invarCalphanew[13],invarCalphanew[11],invarCalphanew[20],invarCalphanew[16],invarCalphanew[26],invarCalphanew[3],invarCalphaorg[27],invarCalphanew[4],invarCalphaorg[6],invarCalphanew[14],invarCalphanew[1],invarCalphaorg[28],invarCalphaorg[25],invarCalphaorg[24],invarCalphaorg[21],invarCalphaorg[19],invarCalphaorg[23],invarCalphaorg[22],invarCalphanew[15],invarCalphaorg[18],invarCalphanew[17]);
    }

  /*
  if(SmoothenBackBone==1){
    if(CylinderTransform==0){
      fprintf(farrayo," %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e\n",smoothbacbonelength/vs[0], int12()/vs[1],inta12()/vs[2],int12_34()/vs[3],inta12_34()/vs[4],int12_a34()/vs[5],inta12_a34()/vs[6],int13_24()/vs[7],inta13_24()/vs[8],int13_a24()/vs[9],inta13_a24()/vs[10],int14_23()/vs[11],inta14_23()/vs[12],int14_a23()/vs[13],inta14_a23()/vs[14],int12_34_56()/vs[15],int12_35_46()/vs[16],int12_36_45()/vs[17],int13_24_56()/vs[18],int13_25_46()/vs[19],int13_26_45()/vs[20],int14_23_56()/vs[21],int14_25_36()/vs[22],int14_26_35()/vs[23],int15_23_46()/vs[24],int15_24_36()/vs[25],int15_26_34()/vs[26],int16_23_45()/vs[27],int16_24_35()/vs[28],int16_25_34())/vs[29];}
  } 
  */
  n_proteins++;
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

double rnd()
{
  return rand()/(0.5*RAND_MAX)-1.0;
}


main(int argc, char *argv[])
{
  double a[3]={0.0,0.0,0.0}, b[3]={0.000000001,1.,0.},c[3]={0.,1.,0.0};
  double p1[3]={0.,.5,1.}, p2[3]={0.,.5,-1.0};




  double inv[29];
  int i,j;
  char line[1000];

  if (argc<5)
    {
      printf("Please input directory (with trailing \"/\") followed by \n Averge-gauss_integral_file,\n output and error file names as arguments.\n");
      return;
    }


  ferroro=fopen(argv[4],"w");
  if (ferroro==(FILE *)NULL)
    {
      printf("Sorry - could not open file <%s> \n",argv[4]);
      return;
    }
  farrayo=fopen(argv[3],"w");
  if (farrayo==(FILE *)NULL)
    {
      printf("Sorry - could not open file <%s>\n",argv[3]);      
      return;
    }

  favgauss=fopen(argv[2],"r");
  if (favgauss==(FILE *)NULL)
    {
      printf("Sorry - could not open file <%s>\n",argv[2]);      
      return;
    }

  
  while(fgets(line,999,favgauss)!=NULL)
    {
      sscanf(line,"%li %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&nres,&inv[0],&inv[1],&inv[2],&inv[3],&inv[4],&inv[5],&inv[6],&inv[7],&inv[8],&inv[9],&inv[10],&inv[11],&inv[12],&inv[13],&inv[14],&inv[15],&inv[16],&inv[17],&inv[18],&inv[19],&inv[20],&inv[21],&inv[22],&inv[23],&inv[24],&inv[25],&inv[26],&inv[27],&inv[28]);
      
      for(i=0;i<29;i++) 
	{
	  if(inv[i]<=0 && nres>9)
	    {
	      printf("Warning, In table in <%s> with pre-calculated\n average Gauss integrals the term %i for proteins of length %i is zero or negative\n This is NOT POSSIBLE\n",argv[2],i,nres);
	      fprintf(ferroro,"Warning, In table in <%s> with pre-calculated\n average Gauss integrals the term %i for proteins of length %i is zero or negative\n This is NOT POSSIBLE\n",argv[2],i,nres);
	      return;
	    }
	  refGauss[nres][i]=inv[i];
	}
    }
  if(nres<MAXLENGTH)
    {
      printf("Warning, the table in <%s> with pre-calculated\n average Gauss integrals is not large enough\n",argv[2]);
      fprintf(ferroro, "Warning, the table in <%s> with pre-calculated\n average Gauss integrals is not large enough\n",argv[2]);
      return ;
    }

  /*for(i=1;i<200000;i++){rnd();}*/
  dirwalk(argv[1]);

  /*for(i=-200;i<200;i++){
    fprintf(farrayo,"%e %e\n",i/100.0, presplinedeformationCalphanew(i/100.0 ,28));
  }
  */

  printf("There has been calculated 29 Tuned Gauss Integrals for %i protein strands.\n", n_proteins);
  fclose(farrayo);
  fclose(favgauss);
  fclose(ferroro);

}




