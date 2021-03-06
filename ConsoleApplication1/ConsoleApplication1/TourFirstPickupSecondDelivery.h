#ifndef __TOURFIRSTPICKUP__
#define __TOURFIRSTPICKUP__


#include "OrienteeringProblemWithPickupsAndDeliveries.h"

#include <string>
#include <vector>

class TourFirstPickupSecondDelivery: public OrienteeringProblemWithPickupsAndDeliveries
{
	public:
	TourFirstPickupSecondDelivery(std::string inputFile);
	~TourFirstPickupSecondDelivery();

	private:
	void pickUpPointToChoose(std::vector<int> & nodes);
	void calcTourChoosePickupAndDeliveryPointPairs();
	int getPickUpDeliveryPointPairsOnePointAdded (std::vector<int> unvisitedCities, int startNode);
	std::vector <int> bestPair;

	

};


#endif

