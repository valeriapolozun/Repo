#ifndef __TOURDEFAULT__
#define __TOURDEFAULT__


#include "OrienteeringProblemWithPickupsAndDeliveries.h"

#include <string>

class TourDefault: public OrienteeringProblemWithPickupsAndDeliveries
{
	public:
	TourDefault(std::string inputFile);
	~TourDefault();

	private:
	void calcTourDefault();

};



#endif