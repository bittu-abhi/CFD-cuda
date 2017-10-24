#ifndef sbli_H
#define sbli_H

class point
{
public:
	double stateVar[4];
	double interface[4][2];

	point();
};

void WENO(point *pt,int points,int flagXY);

#endif