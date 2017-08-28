#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include "ausmPlus.h"

using namespace std;

void visual(cell *domain,float t)
{
	fstream myfile1,myfile2,myfile3,myfile4;
	stringstream ss;
	string filename="finalvalues00";
	string end=".csv";
	ss<<filename<<(t/0.000001)<<end;
	filename=ss.str();
	myfile1.open(filename.c_str(),ios::out);
	myfile1<<"X"<<","<<"Y"<<","<<"Z"<<","<<"Rho"<<","<<"U"<<","<<"V"<<","<<"E"<<","<<"flag"<<","<<"nx1"<<","<<"ny1"\
	<<","<<"nx2"<<","<<"ny2"<<","<<"nx3"<<","<<"ny3"<<endl;
	if(myfile1.is_open())
	{
		cout<<"Writing final values....."<<endl;
		myfile1 << fixed;
		myfile1 << setprecision(15);
		for(int i=0;i<50266;i++)
		{
			myfile1<<(1.0/3.0)*(domain[i].nodes[0][0]+domain[i].nodes[1][0]+domain[i].nodes[2][0])<<","\
			<<(1.0/3.0)*(domain[i].nodes[0][1]+domain[i].nodes[1][1]+domain[i].nodes[2][1])<<","<<\
			0<<","<<\
			domain[i].stateVar[0]<<","<<\
			domain[i].stateVar[1]/domain[i].stateVar[0]<<","<<\
			domain[i].stateVar[2]/domain[i].stateVar[0]<<","<<\
			domain[i].stateVar[3]/domain[i].stateVar[0]<<","<<\
			domain[i].flag<<","<<\
			domain[i].norms[0][0]<<","<<domain[i].norms[0][1]<<","<<\
			domain[i].norms[1][0]<<","<<domain[i].norms[1][1]<<","<<\
			domain[i].norms[2][0]<<","<<domain[i].norms[2][1]<<","<<\
			endl;
		}
	}
	else
	{
		cout<<"Cannot open the required file"<<endl;
	}
	myfile1.close();
/*
	string filename1="convectiveflux00";
	string end1=".csv";
	stringstream ss1;
	ss1<<filename1<<(t/0.000001)<<end1;
	filename1=ss1.str();

	myfile2.open(filename1.c_str(),ios::out);
	myfile2<<"X"<<","<<"Y"<<","<<"Rho"<<","<<"Rho*U"<<","<<"Rho*V"<<","<<"Rho*E"<<","<<"flag"<<endl;
	if(myfile2.is_open())
	{
		cout<<"Writing convective fluxes....."<<endl;
		myfile2 << fixed;
		myfile2 << setprecision(15);
		for(int i=0;i<50266;i++)
		{
			myfile2<<(1.0/3.0)*(domain[i].nodes[0][0]+domain[i].nodes[1][0]+domain[i].nodes[2][0])<<","\
			<<(1.0/3.0)*(domain[i].nodes[0][1]+domain[i].nodes[1][1]+domain[i].nodes[2][1])<<","<<\
			domain[i].convflux[0][0]+domain[i].convflux[1][0]+domain[i].convflux[2][0]<<","<<\
			domain[i].convflux[0][1]+domain[i].convflux[1][1]+domain[i].convflux[2][1]<<","<<\
			domain[i].convflux[0][2]+domain[i].convflux[1][2]+domain[i].convflux[2][2]<<","<<\
			domain[i].convflux[0][3]+domain[i].convflux[1][3]+domain[i].convflux[2][3]<<","<<\
			domain[i].flag<<endl;
		}
	}
	else
	{
		cout<<"Cannot open the required file"<<endl;
	}
	myfile2.close();
	string filename2="diffusiveflux00";
	string end2=".csv";
	stringstream ss2;
	ss2<<filename2<<(t/0.000001)<<end2;
	filename2=ss2.str();

	myfile3.open(filename2.c_str(),ios::out);
	myfile3<<"X"<<","<<"Y"<<","<<"Rho"<<","<<"Rho*U"<<","<<"Rho*V"<<","<<"Rho*E"<<","<<"flag"<<endl;
	if(myfile3.is_open())
	{
		myfile3 << fixed;
		myfile3 << setprecision(15);
		cout<<"Writing diffusive fluxes....."<<endl;
		for(int i=0;i<50266;i++)
		{
			myfile3<<(1.0/3.0)*(domain[i].nodes[0][0]+domain[i].nodes[1][0]+domain[i].nodes[2][0])<<","\
			<<(1.0/3.0)*(domain[i].nodes[0][1]+domain[i].nodes[1][1]+domain[i].nodes[2][1])<<","<<\
			domain[i].diffflux[0][0]+domain[i].diffflux[1][0]+domain[i].diffflux[2][0]<<","<<\
			domain[i].diffflux[0][1]+domain[i].diffflux[1][1]+domain[i].diffflux[2][1]<<","<<\
			domain[i].diffflux[0][2]+domain[i].diffflux[1][2]+domain[i].diffflux[2][2]<<","<<\
			domain[i].diffflux[0][3]+domain[i].diffflux[1][3]+domain[i].diffflux[2][3]<<","<<\
			domain[i].flag<<endl;
		}
	}
	else
	{
		cout<<"Cannot open the required file"<<endl;
	}
	myfile3.close();
	string filename3="pressureflux00";
	string end3=".csv";
	stringstream ss3;
	ss3<<filename3<<(t/0.000001)<<end3;
	filename3=ss3.str();
	myfile4.open(filename3.c_str(),ios::out);
	myfile4<<"X"<<","<<"Y"<<","<<"Rho*U"<<","<<"Rho*V"<<","<<"flag"<<endl;
	if(myfile4.is_open())
	{
		myfile4 << fixed;
		myfile4 << setprecision(15);
		cout<<"Writing pressure fluxes....."<<endl;
		for(int i=0;i<50266;i++)
		{
			myfile4<<(1.0/3.0)*(domain[i].nodes[0][0]+domain[i].nodes[1][0]+domain[i].nodes[2][0])<<","\
			<<(1.0/3.0)*(domain[i].nodes[0][1]+domain[i].nodes[1][1]+domain[i].nodes[2][1])<<","<<\
			domain[i].presflux[0][0]+domain[i].presflux[1][0]+domain[i].presflux[2][0]<<","<<\
			domain[i].presflux[0][1]+domain[i].presflux[1][1]+domain[i].presflux[2][1]<<","<<\
			domain[i].flag<<endl;
		}
	}
	else
	{
		cout<<"Cannot open the required file"<<endl;
	}
	myfile4.close();
	*/
}
