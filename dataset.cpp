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
#include "NAME.h"


using namespace std;

int main(){
	// open file in read-only mode
	ifstream test_file("pppp.txt"); // get the test file contains 5 Protein File Names

	// check for file error
	if(!test_file)
	{
		cerr<<"Error Opening File"<<endl;
		exit(1);
	}
	else{
		count << test_file<< "File is sucessfully opened"<< endl;
	}

	// display data in test_file
	char data[100];
	test_file >> data;
	cout << "display " << test_file << " :" << endl;
	cout << data << endl;

	// create an output file
	ofstream myfile1; 						// ofstream: create files and write info to file
	myfile1.open ("final.txt");		// why do we need this final.txt?

	std::string PDB;
	while(test_file>>PDB) 	// reading from test_file
	{
		int sl;
		cout<<PDB<<"\n"; // priting: protein file names

		std:: string pdbbmr("C:\\Users\\mreza\\Desktop\\lg\\"+PDB+".txt"); // file path of the name
		fstream fp(pdbbmr.c_str());

		/* Example of the table
		/ ---------------------------------------------------------------------------------------------
		/ Atom_num     Res_Num     ATOM_Name	  AmincAcid_Name	 X Y Z         CHEMICAL_SHIFT
		/ n  			n		 	name1		  name1				 x,y,z		   num1
		/ n++			n++			name2		  name2				 x,y,za2	   num2
		/ -------------------------------------------------------------------------------------------- */

		// // Variables
		// string ATOM[10000];		// atom name
		// string AMIN[10000];		// amino acid name
		// // double x[10000];		// x,y,z coordinates
		// // double y[10000];
		// // double z[10000];
		// double chem[10000];		// chemical shift

		double X[3000],Y[3000],Z[3000]; // why is there another x,y,z array for

		// var for calculating
		double An[3000],An2[3000];
		double g,G,rs[10],dc,gh,gH,Gh,GH,gO,GO,gOH,GOH,gC,GC,gcd1,Gcd1,gce1,Gce1,gce2,Gce2,gz,Gz,gg,Gg,gcd2,Gcd2;


		// resgister number value
		int atomnumber,resnumber,res[3000];
		int b,m;

		// what this for???
		rs[1]=4.90; rs[2]=5.00; rs[3]=5.10; rs[4]=5.20;
		rs[5]=5.30; rs[6]=5.40; rs[7]=5.50;rs[8]=5.6;

		// hmn.....
		double direct,vecx,vecy,vecz,x12,x34,x32,y12,y34,y32,z12,z34,z32,D1,D2,proj1,proj2,dot1,dot2,dot3,dot4,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,d32;
		double directt,vecxx,vecyy,veczz,xx12,xx34,xx32,yy12,yy34,yy32,zz12,zz34,zz32,DD1,DD2,projj1,projj2,dott1,dott2,dott3,dott4,xx1,xx2,xx3,xx4,yy1,yy2,yy3,yy4,zz1,zz2,zz3,zz4,dd32;

		/* --------------------------------------------------
		/ Dihedral Angle Calculation:
		/
		/ The angle between two planes, where each plane
		/ is defined by three atoms.
		/ --------------------------------------------------*/

		// Var to calculate DHA
		// Perhaps can put this into a function and just taking the values???

		int i,j; // counters i,j
		string ATOM[10000];		// atom name
		string AMIN[10000];		// amino acid name
		double x[10000];		// x,y,z coordinates
		double y[10000];
		double z[10000];
		double chem[10000];		// chemical shift

		int R,NUMBER;
		string ATOM2,AMIN2;
		double xa2,ya2,za2;
		double chem2;

		/* --------------------------------------------------
		/ what is this for loop for?????
		/
		// Checking for N, CA, C
		// have counter for N, CA, C
		/ --------------------------------------------------*/

		while(fp>>NUMBER>>R>>AMIN2>>ATOM2>>xa2>>ya2>>za2>>chem2)
		{
			ATOM[NUMBER] = ATOM2;
			AMIN[NUMBER] = AMIN2;
			x[NUMBER] = xa2;
			y[NUMBER] = ya2;
			z[NUMBER] = za2;
			chem[NUMBER] = chem2;

			i=NUMBER;
			//cout<<NUMBER<<"\n";
			atomnumber=NUMBER;
			resnumber=R;

			if (ATOM[NUMBER]=="N")
			{
				X[R]=x[NUMBER];
				Y[R]=y[NUMBER];
				Z[R]=z[NUMBER];
				res[R]=NUMBER;
				//	cout<<chem[NUMBER]<<"\n";
			}

			if (R!=1)
			{
				if (ATOM2=="N")
				{
					//	cout<<ATOM2<<"\n";
					x2=x[i];
					y2=y[i];
					z3=z[i];
					b++;
				}
				if (ATOM[i]=="CA")
				{
					x3=x[i];
					y3=y[i];
					z3=z[i];
					b++;
				}

				if (ATOM[i]=="C")
				{
					x4=x[i];
					y4=y[i];
					z4=z[i];
					b++;
				}

				// b is a counter for.....?
				if(b>3)
				{
					An[R]=Angel(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4);
					//cout<<An[R]<<"\n";
					b=0;	// set b back to 0
				}
			} // end of if (R!=1)

			if (ATOM[i]=="C")
			{
				b++;
				x1=x[i];	y1=y[i];	z1=z[i];
			}

			if (R!=1)
			{
				if (ATOM[i]=="N")
				{
					xx4=x[i]; yy4=y[i]; zz4=z[i];
					j++;
				}

				if(j>3)
				{
					An2[R]=Angel(xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4);
					//cout<<An2[R]<<"\n";
					j=0;
				}

			}

			if (ATOM[i]=="CA")
			{
				xx2=x[i];	yy2=y[i];	zz2=z[i];
				j++;
			}
			if (ATOM[i]=="C")
			{
				xx3=x[i];	yy3=y[i];	zz3=z[i];
				j++;
			}

			if (ATOM[i]=="N")
			{
				j++;
				xx1=x[i];	yy1=y[i];	zz1=z[i];

			}
			////////////////////////////////////////////////

		} // COMPLETE calculating ????; end while loop(fp>>NUMBER>>R>>AMIN2>>ATOM2>>xa2>>ya2>>za2>>chem2)

		i=0;
		res[0]=0;
		cout << NUMBER; // ???

		for (R=1;R<resnumber+1;R++)
		{
			for(i=res[R-1]+1;i<res[R]+1;i++)
			{
				double G,g,gh,gH,Gh,GH,gO,GO,gOH,GOH,gC,GC,gcd1,Gcd1,gce1,Gce1,gce2,Gce2,gz,Gz,gg,Gg,gcd2,Gcd2;
				//cout<<i<<"\n";

				string avv,ex,ey,ez;
				avv=ATOM[i];
				ex=AMIN[res[R]];
				ey=AMIN[res[R-1]];
				ez=AMIN[res[R+1]];
				//cout<<amin(agg)<<"\n";

				myfile1<<std::fixed << std::setprecision(8) <<AMU(amin(ey),amin(ex),amin(ez))<<"  "<<atom(avv)<<"  "<<An[R-1]<<"  "<<An[R]<<"  "<<An[R+1]<<"  "<<An2[R-1]<<"  "<<An2[R]<<"  "<<An2[R+1]<<"  ";

				if (i==NUMBER)
				{
					cout<<NUMBER<<"\n";
					cout<<AMU(amin(ey),amin(ex),amin(ez))<<"  "<<atom(avv)<<"  "<<An[R-1]<<"  "<<An[R]<<"  "<<An[R+1]<<"  "<<An2[R-1]<<"  "<<An2[R]<<"  "<<An2[R+1]<<"  ";
				}

				//////////////////////////////////
				for(m=1;m<resnumber+1;m++)
				{

					if(R!=m)
					{
						g=Gvector(x[R],y[R],z[R],x[m],y[m],z[m],5.5,rs[1]);
					}

					G=G+g;
				}

				myfile1<<std::fixed << std::setprecision(8) <<G<<"  ";
				//cout<<G<<"  "<<i<<"\n";
				G=0;


				for (j=1;j<9;j++)
				{
					for(m=1;m<atomnumber+1;m++){
						if(m!=i)
						{
							if ( ATOM[m]=="N"){
								g=Gvector(x[i],y[i],z[i],x[m],y[m],z[m],6.7,rs[j]);
								G=G+g;
							}
						}
					} // end

					if (G>1){
						sl++;
					}

					myfile1<<std::fixed << std::setprecision(8) <<G<<"  ";
					//cout<<Gh<<"  "<<i<<"  "<<j<<"\n";
					G=0;
				}


				for (j=1;j<9;j++)
				{
					for(m=1;m<atomnumber+1;m++)
					{
						if(m!=i)
						{
							if ( ATOM[m]=="C")
							{

								g=Gvector(x[i],y[i],z[i],x[m],y[m],z[m],6.7,rs[j]);
								G=G+g;
							}
						}
					} // end inner for loop

					myfile1<<std::fixed << std::setprecision(8) <<G<<"  ";
					//cout<<Gh<<"  ";
					G=0;
				}

				for (j=1;j<9;j++)
				{

					for(m=1;m<atomnumber+1;m++)
					{
						if(m!=i)
						{
							if ( ATOM[m]=="H")
							{
								g=Gvector(x[i],y[i],z[i],x[m],y[m],z[m],6.7,rs[j]);
								G=G+g;
							}
						}
					}

					myfile1<<std::fixed << std::setprecision(8) <<G<<"  ";
					G=0;
				}

				////////////////////////////////////////////
				for(m=1;m<atomnumber+1;m++)
				{
					if(m!=i)
					{
						if ( ATOM[m]=="H"){
							gH=Gvector(x[i],y[i],z[i],x[m],y[m],z[m],6.0,rs[1]);
							GH=GH+gH;
						}
						if ( ATOM[m]=="O"){
							gO=Gvector(x[i],y[i],z[i],x[m],y[m],z[m],6.3,rs[1]);
							GO=GO+gO;
						}
						if ( ATOM[m]=="C"){
							gC=Gvector(x[i],y[i],z[i],x[m],y[m],z[m],6.0,rs[1]);
							GC=GC+gC;
						}
						if ( ATOM[m]=="OH"){
							gOH=Gvector(x[i],y[i],z[i],x[m],y[m],z[m],6.0,rs[1]);
							GOH=GOH+gOH;
						}
						if ( ATOM[m]=="CD1"){
							gcd1=Gvector(x[i],y[i],z[i],x[m],y[m],z[m],5.5,rs[1]);
							Gcd1=Gcd1+gcd1;
						}
						if ( ATOM[m]=="CD2"){
							gcd2=Gvector(x[i],y[i],z[i],x[m],y[m],z[m],5.5,rs[1]);
							Gcd2=Gcd2+gcd2;
						}
						if ( ATOM[m]=="CE1"){
							gce1=Gvector(x[i],y[i],z[i],x[m],y[m],z[m],5.5,rs[1]);
							Gce1=Gce1+gce1;
						}
						if ( ATOM[m]=="CE2"){
							gce2=Gvector(x[i],y[i],z[i],x[m],y[m],z[m],5.5,rs[1]);
							Gce2=Gce2+gce2;
						}
						if ( ATOM[m]=="CZ"){
							gz=Gvector(x[i],y[i],z[i],x[m],y[m],z[m],5.5,rs[1]);
							Gz=Gz+gz;
						}
						if ( ATOM[m]=="CG"){
							gg=Gvector(x[i],y[i],z[i],x[m],y[m],z[m],5.5,rs[1]);
							Gg=Gg+gg;
						}
					}
				} // end for loop

				myfile1<<std::fixed << std::setprecision(8) <<GH<<"  "<<GO<<"  "<<GOH<<"  "<<GC<<"  "<<Gcd1<<"  "<<Gcd2<<"  "<<Gce1<<"  "<<Gce2<<"  "<<Gz<<"  "<<Gg<<"  ";

				if (i==200){
					//	cout <<std::fixed << std::setprecision(8) <<GH<<"  "<<GO<<"  "<<GOH<<"  "<<GC<<"  "<<Gcd1<<"  "<<Gcd2<<"  "<<Gce1<<"  "<<Gce2<<"  "<<Gz<<"  "<<Gg<<" \n ";
				}

				GH=0; 	GO=0;  GOH=0; 	GC=0; 	Gcd1=0; 	Gcd2=0; 	Gce1=0; 	Gce2=0; 	Gz=0; 	Gg=0;
				myfile1<<std::fixed << std::setprecision(8) <<chem[i]/300<<"\n";
				//	cout<< sl<<"\n";

			}
		} // 	what is this for loop??????

	}

	myfile1.close();
	return 0;

} // end main
