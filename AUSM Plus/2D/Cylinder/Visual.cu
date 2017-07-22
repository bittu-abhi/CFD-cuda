#include <iostream>
#include <fstream>
#include "ausmPlus.h"

using namespace std;

void visual(cell *domain)
{
	fstream myfile1,myfile2,myfile3,myfile4;

	myfile1.open("finalvalues.csv",ios::out);
	myfile1<<"X"<<","<<"Y"<<","<<"Z"<<","<<"Rho"<<","<<"U"<<","<<"V"<<","<<"E"<<","<<"flag"<<endl;
	if(myfile1.is_open())
	{
		cout<<"Writing final values....."<<endl;
		for(int i=0;i<25000;i++)
		{
			myfile1<<0.25*(domain[i].nodes[0][0]+domain[i].nodes[1][0]+domain[i].nodes[2][0]+domain[i].nodes[3][0])<<","\
			<<0.25*(domain[i].nodes[0][1]+domain[i].nodes[1][1]+domain[i].nodes[2][1]+domain[i].nodes[3][1])<<","<<\
			0<<","<<\
			domain[i].stateVar[0]<<","<<\
			domain[i].stateVar[1]/domain[i].stateVar[0]<<","<<\
			domain[i].stateVar[2]/domain[i].stateVar[0]<<","<<\
			domain[i].stateVar[3]/domain[i].stateVar[0]<<","<<\
			domain[i].flag<<endl;
		}
	}
	else
	{
		cout<<"Cannot open the required file"<<endl;
	}
	myfile1.close();

	myfile2.open("convectiveflux.csv",ios::out);
	myfile2<<"X"<<","<<"Y"<<","<<"Rho"<<","<<"Rho*U"<<","<<"Rho*V"<<","<<"Rho*E"<<","<<"flag"<<endl;
	if(myfile2.is_open())
	{
		cout<<"Writing convective fluxes....."<<endl;
		for(int i=0;i<25000;i++)
		{
			myfile2<<0.25*(domain[i].nodes[0][0]+domain[i].nodes[1][0]+domain[i].nodes[2][0]+domain[i].nodes[3][0])<<","\
			<<0.25*(domain[i].nodes[0][1]+domain[i].nodes[1][1]+domain[i].nodes[2][1]+domain[i].nodes[3][1])<<","<<\
			domain[i].convflux[0][0]+domain[i].convflux[1][0]+domain[i].convflux[2][0]+domain[i].convflux[3][0]<<","<<\
			domain[i].convflux[0][1]+domain[i].convflux[1][1]+domain[i].convflux[2][1]+domain[i].convflux[3][1]<<","<<\
			domain[i].convflux[0][2]+domain[i].convflux[1][2]+domain[i].convflux[2][2]+domain[i].convflux[3][2]<<","<<\
			domain[i].convflux[0][3]+domain[i].convflux[1][3]+domain[i].convflux[2][3]+domain[i].convflux[3][3]<<","<<\
			domain[i].flag<<endl;
		}
	}
	else
	{
		cout<<"Cannot open the required file"<<endl;
	}
	myfile2.close();

	myfile3.open("diffusiveflux.csv",ios::out);
	myfile3<<"X"<<","<<"Y"<<","<<"Rho"<<","<<"Rho*U"<<","<<"Rho*V"<<","<<"Rho*E"<<","<<"flag"<<endl;
	if(myfile3.is_open())
	{
		cout<<"Writing diffusive fluxes....."<<endl;
		for(int i=0;i<25000;i++)
		{
			myfile3<<0.25*(domain[i].nodes[0][0]+domain[i].nodes[1][0]+domain[i].nodes[2][0]+domain[i].nodes[3][0])<<","\
			<<0.25*(domain[i].nodes[0][1]+domain[i].nodes[1][1]+domain[i].nodes[2][1]+domain[i].nodes[3][1])<<","<<\
			domain[i].diffflux[0][0]+domain[i].diffflux[1][0]+domain[i].diffflux[2][0]+domain[i].diffflux[3][0]<<","<<\
			domain[i].diffflux[0][1]+domain[i].diffflux[1][1]+domain[i].diffflux[2][1]+domain[i].diffflux[3][1]<<","<<\
			domain[i].diffflux[0][2]+domain[i].diffflux[1][2]+domain[i].diffflux[2][2]+domain[i].diffflux[3][2]<<","<<\
			domain[i].diffflux[0][3]+domain[i].diffflux[1][3]+domain[i].diffflux[2][3]+domain[i].diffflux[3][3]<<","<<\
			domain[i].flag<<endl;
		}
	}
	else
	{
		cout<<"Cannot open the required file"<<endl;
	}
	myfile3.close();
	
	myfile4.open("pressureflux.csv",ios::out);
	myfile4<<"X"<<","<<"Y"<<","<<"Rho*U"<<","<<"Rho*V"<<","<<"flag"<<endl;
	if(myfile4.is_open())
	{
		cout<<"Writing pressure fluxes....."<<endl;
		for(int i=0;i<25000;i++)
		{
			myfile4<<0.25*(domain[i].nodes[0][0]+domain[i].nodes[1][0]+domain[i].nodes[2][0]+domain[i].nodes[3][0])<<","\
			<<0.25*(domain[i].nodes[0][1]+domain[i].nodes[1][1]+domain[i].nodes[2][1]+domain[i].nodes[3][1])<<","<<\
			domain[i].presflux[0][0]+domain[i].presflux[1][0]+domain[i].presflux[2][0]+domain[i].presflux[3][0]<<","<<\
			domain[i].presflux[0][1]+domain[i].presflux[1][1]+domain[i].presflux[2][1]+domain[i].presflux[3][1]<<","<<\
			domain[i].flag<<endl;
		}
	}
	else
	{
		cout<<"Cannot open the required file"<<endl;
	}
	myfile4.close();	
}
