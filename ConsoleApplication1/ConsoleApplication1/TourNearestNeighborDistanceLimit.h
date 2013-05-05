#ifndef __TOURNEARESTNEIGHBORDISTANCELIMIT__
#define __TOURNEARESTNEIGHBORDISTANCELIMIT__


#include "OrienteeringProblemWithPickupsAndDeliveries.h"

#include <string>

class TourNearestNeighborDistanceLimit: public OrienteeringProblemWithPickupsAndDeliveries
{
	public:
	TourNearestNeighborDistanceLimit(std::string inputFile);
	~TourNearestNeighborDistanceLimit();

	private:
	void calcTourNearestNeighborDistanceLimit();

};


#endif

