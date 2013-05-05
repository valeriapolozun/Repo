#ifndef __TOURNEARESTNEIGHBOR__
#define __TOURNEARESTNEIGHBOR__


#include "OrienteeringProblemWithPickupsAndDeliveries.h"

#include <string>

class TourNearestNeighbor: public OrienteeringProblemWithPickupsAndDeliveries
{
	public:
	TourNearestNeighbor(std::string inputFile);
	~TourNearestNeighbor();

	private:
	void calcTourNearestNeighbor();

};



#endif