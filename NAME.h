#ifndef __NAME__
#define __NAME__

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <time.h>
#include <iomanip>


using namespace std;
 float p=3.14159;
int AM[10000], AR[1000],R,AT[10000],NUMBER;


int amin(string AMIN){
 NUMBER=1;
if (AMIN=="PRO"){ AM[NUMBER]=1;AR[R]=1;}
if (AMIN=="ALA"){ AM[NUMBER]=2;AR[R]=2;}
if (AMIN=="LEU"){ AM[NUMBER]=3;AR[R]=3;}
if (AMIN=="TRP"){ AM[NUMBER]=4;AR[R]=4;}
if (AMIN=="GLU"){ AM[NUMBER]=5;AR[R]=5;}
if (AMIN=="LYS"){ AM[NUMBER]=6;AR[R]=6;}
if (AMIN=="ILE"){ AM[NUMBER]=7;AR[R]=7;}
if (AMIN=="ARG"){ AM[NUMBER]=8;AR[R]=8;}
if (AMIN=="ASP"){ AM[NUMBER]=9;AR[R]=9;}
if (AMIN=="VAL"){ AM[NUMBER]=10;AR[R]=10;}
if (AMIN=="MET"){ AM[NUMBER]=11;AR[R]=11;}
if (AMIN=="GLN"){ AM[NUMBER]=12;AR[R]=12;}
if (AMIN=="TYR"){ AM[NUMBER]=13;AR[R]=12;}
if (AMIN=="PHE"){ AM[NUMBER]=14;AR[R]=14;}
if (AMIN=="CYS"){ AM[NUMBER]=15;AR[R]=15;}
if (AMIN=="HSE"){ AM[NUMBER]=16;AR[R]=16;}
if (AMIN=="ASN"){ AM[NUMBER]=17;AR[R]=17;}
if (AMIN=="SER"){ AM[NUMBER]=18;AR[R]=18;}
if (AMIN=="THR"){ AM[NUMBER]=19;AR[R]=19;}
if (AMIN=="GLY"){ AM[NUMBER]=20;AR[R]=20;}
return AM[NUMBER];
}
double atom(string ATOM){
	if (ATOM=="H"){AT[NUMBER]=1;}
if (ATOM=="HA"){AT[NUMBER]=2;}
if (ATOM=="HB"){AT[NUMBER]=3;}
if (ATOM=="C"){AT[NUMBER]=4;}
if (ATOM=="CA"){AT[NUMBER]=5;}
if (ATOM=="CB"){AT[NUMBER]=6;}
if (ATOM=="N"){AT[NUMBER]=7;}
if (ATOM=="HB2"){AT[NUMBER]=8;}
if (ATOM=="HB3"){AT[NUMBER]=9;}
if (ATOM=="HG2"){AT[NUMBER]=10;}
if (ATOM=="HG22"){AT[NUMBER]=11;}
	if (ATOM=="HG21"){AT[NUMBER]=12;}
	if (ATOM=="HG23"){AT[NUMBER]=13;}

if (ATOM=="HG3"){AT[NUMBER]=14;}
if (ATOM=="HD2"){AT[NUMBER]=15;}
if (ATOM=="HD3"){AT[NUMBER]=16;}
if (ATOM=="HE"){AT[NUMBER]=17;}
if (ATOM=="HH11"){AT[NUMBER]=18;}
if (ATOM=="HH12"){AT[NUMBER]=19;}
if (ATOM=="HH21"){AT[NUMBER]=20;}
if (ATOM=="HH22"){AT[NUMBER]=21;}
if (ATOM=="CG"){AT[NUMBER]=22;}
if (ATOM=="CD"){AT[NUMBER]=23;}
if (ATOM=="CZ"){AT[NUMBER]=24;}
if (ATOM=="NE"){AT[NUMBER]=25;}
if (ATOM=="NH1"){AT[NUMBER]=26;}
if (ATOM=="NH2"){AT[NUMBER]=27;}
if (ATOM=="HD21"){AT[NUMBER]=28;}
if (ATOM=="HD22"){AT[NUMBER]=29;}
if (ATOM=="ND2"){AT[NUMBER]=30;}
if (ATOM=="HG"){AT[NUMBER]=31;}
if (ATOM=="HE2"){AT[NUMBER]=32;}
if (ATOM=="HE21"){AT[NUMBER]=33;}
if (ATOM=="HE22"){AT[NUMBER]=34;}
if (ATOM=="NE2"){AT[NUMBER]=35;}
if (ATOM=="HA2"){AT[NUMBER]=36;}
if (ATOM=="HA3"){AT[NUMBER]=37;}
if (ATOM=="HD1"){AT[NUMBER]=38;}
if (ATOM=="HE1"){AT[NUMBER]=39;}
if (ATOM=="CD2"){AT[NUMBER]=40;}
if (ATOM=="CE1"){AT[NUMBER]=41;}
if (ATOM=="ND1"){AT[NUMBER]=42;}
if (ATOM=="HG12"){AT[NUMBER]=43;}
if (ATOM=="HG13"){AT[NUMBER]=44;}
if (ATOM=="CG1"){AT[NUMBER]=45;}
if (ATOM=="CG2"){AT[NUMBER]=46;}
if (ATOM=="CD1"){AT[NUMBER]=47;}
if (ATOM=="HE3"){AT[NUMBER]=48;}
if (ATOM=="HZ"){AT[NUMBER]=49;}

if (ATOM=="CE"){AT[NUMBER]=50;}
if (ATOM=="NZ"){AT[NUMBER]=51;}
if (ATOM=="CE2"){AT[NUMBER]=52;}
if (ATOM=="HG1"){AT[NUMBER]=53;}
if (ATOM=="H2"){AT[NUMBER]=54;}
if (ATOM=="H3"){AT[NUMBER]=55;}


if (ATOM=="HZ2"){AT[NUMBER]=56;}
if (ATOM=="HZ3"){AT[NUMBER]=57;}
if (ATOM=="HH2"){AT[NUMBER]=58;}
if (ATOM=="CE3"){AT[NUMBER]=59;}
if (ATOM=="CZ2"){AT[NUMBER]=60;}
if (ATOM=="CZ3"){AT[NUMBER]=61;}
if (ATOM=="CH2"){AT[NUMBER]=62;}
if (ATOM=="NE1"){AT[NUMBER]=63;}
if (ATOM=="HH"){AT[NUMBER]=64;}
if (ATOM[0]=='O'){AT[NUMBER]=100;}
if (ATOM[0]=='S'){AT[NUMBER]=100;}
if (ATOM[0]=='D'){AT[NUMBER]=100;}

	if (ATOM=="HZ1"){AT[NUMBER]=65;}

	if (ATOM=="HG11"){AT[NUMBER]=66;}
	if (ATOM=="HG13"){AT[NUMBER]=67;}

	if (ATOM=="HB1"){AT[NUMBER]=68;}



	if (ATOM=="HD23"){AT[NUMBER]=69;}
	if (ATOM=="HD11"){AT[NUMBER]=70;}
	if (ATOM=="HD12"){AT[NUMBER]=71;}
	if (ATOM=="HD13"){AT[NUMBER]=72;}
	double er;
	er=AT[NUMBER];
	return er/100;
}
/////////////////////////////////////////////////////
float Angel(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3,float x4, float y4,float z4){
 float direct, vecx,vecy,vecz,x12,x34,x32,y12,y34,y32,z12,z34,z32,D1,D2,proj1,proj2,dot1,dot2,dot3,dot4,d32,Angell;
 x12=x1-x2;
  y12=y1-y2;
  z12=z1-z2;
  
  x34=x3-x4;
  y34=y3-y4;
  z34=z3-z4;
  
  x32=x3-x2;
  y32=y3-y2;
  z32=z3-z2;
  
  dot1=(x12*x32) + (y12*y32) + (z12*z32);
  d32=(pow(x32,2))+(pow(y32,2))+(pow(z32,2));
  proj1=dot1/d32;
  x12=x12-(proj1*x32);
  y12=y12-(proj1*y32);
  z12=z12-(proj1*z32);
  
  dot2=(x34*x32)+(y34*y32)+(z34*z32);
  proj2=dot2/d32;
  x34=x34-(proj2*x32);
  y34=y34-(proj2*y32);
  z34=z34-(proj2*z32);
  
  vecx=-(y12*z34)+(z12*y34);
  vecy=-(z12*x34)+(x12*z34);
  vecz=-(x12*y34)+(y12*x34);
  
  dot3=(vecx*x32)+(vecy*y32)+(vecz*z32);

    if (dot3<0)
    {
      direct=-1;
    }
  else
    {
      direct=1;
    }
  
  dot4=-((x12*x34)+(y12*y34)+(z12*z34));
  D1=pow(((pow(x12,2))+(pow(y12,2))+(pow(z12,2))),0.5);
  D2=pow(((pow(x34,2))+(pow(y34,2))+(pow(z34,2))),0.5);

  Angell=((acos(dot4/(D1*D2))*180))*direct/p;

  Angell= (Angell/360.0)+0.50;
 

  return Angell;

  }
  
///////////////////////////////////////////////////////////

		double  Gvector(double X1,double Y1, double Z1, double X2,double Y2, double Z2, double dc, double rs){
			double Fc,g;
			double NN=4.0;
    	double	d=pow(((pow((X1-X2),2))+(pow((Y1-Y2),2))+(pow((Z1-Z2),2))),0.5);
      
    
    		if(d<=dc)
    		{
    			Fc=(0.5*cos((p*d)/dc))+0.5;
    	
			}
			if (d>dc)
			{
				Fc=0;
			}
		
			g=(exp((-NN)*(pow((d-rs),2))))*Fc;
			
			return g;
			
		    
		    }
		    	
		    	double AMU(int x, int y, int z){
		    		
		    		double s;
		    		s=((x*10000)+(y*100)+(z));
		    		return s/212121;
				}
#endif
