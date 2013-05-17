#ifndef __TOURPICKUP__
#define __TOURPICKUP__


#include "OrienteeringProblemWithPickupsAndDeliveries.h"

#include <string>
#include <vector>

class TourPickupDeliveryPointPairs: public OrienteeringProblemWithPickupsAndDeliveries
{
	public:
	TourPickupDeliveryPointPairs(std::string inputFile);
	~TourPickupDeliveryPointPairs();

	private:
	void pickUpPointToChoose(std::vector<int> & nodes);
	void calcTourChoosePickupAndDeliveryPointPairs();
	int getPickUpDeliveryPointPairsOnePointAdded (std::vector<int> unvisitedCities, int startNode);
	void getPickUpDeliveryPointPairsTwoPointsAdded (std::vector<int> unvisitedCities, int startNode, std::vector <int> & bestPair);
	std::vector <int> bestPairs;
	void calcTourChoosePickupAndDeliveryPointPairs2();
	

};


#endif
