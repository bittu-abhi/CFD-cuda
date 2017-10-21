#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "lid.h"

using namespace std;

void visu(point *ptr,double t,double deltat,int points,double deltax, double deltay)
{
	fstream w;
	stringstream ss;
	string filename="fi0_";
	string end=".csv";
	ss<<filename<<(t/deltat/1000)<<end;
	filename=ss.str();
	w.open(filename.c_str(),ios::out);
	w<<"X"<<","<<"Y"<<","<<"Z"<<","<<"Omega"<<","<<"psi"<<","<<"U"<<","<<"V"<<endl;
	if(w.is_open())
	{
		cout<<"Writing final values....."<<endl;
		w << fixed;
		w << setprecision(15);
		for(int i=0;i<points;i++)
		{
			for (int j = 0; j < points; j++)
			{
				w<<i*deltax<<","<<j*deltay<<","<<0<<","<<ptr[j+i*points].w<<","<<ptr[j+i*points].psi<<","<<ptr[j+i*points].delpsi[1]<<","<<-ptr[j+i*points].delpsi[0]<<endl;
			
			}		
		}
	}
	else
	{
		cout<<"Cannot open the required file"<<endl;
	}
	w.close();
}