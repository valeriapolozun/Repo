#ifndef __COORDINATES__
#define __COORDINATES__

struct Coordinates  // record all the information of points
{
	double x,y;//coordinates
	int quantity; //quantity ( + if pick up, - if delivery) 
	int profit;//value
};

#endif