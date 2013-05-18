#ifndef __TOURPICKUPDELIVERYPARALLEL__
#define __TOURPICKUPDELIVERYPARALLEL__


#include "OrienteeringProblemWithPickupsAndDeliveries.h"

#include <string>
#include <vector>

class TourPickupDeliveryPointPairsParallel: public OrienteeringProblemWithPickupsAndDeliveries
{
	public:
	TourPickupDeliveryPointPairsParallel(std::string inputFile);
	~TourPickupDeliveryPointPairsParallel();

	private:
	void pickUpPointToChoose(std::vector<int> & nodes);
	void calcTourChoosePickupAndDeliveryPointPairs();
	void getPickUpDeliveryPointPairsTwoPointsAdded (std::vector<int> unvisitedCities, int startNode, std::vector <int> & bestPair);
	std::vector <int> bestPairs;
	void calcTourChoosePickupAndDeliveryPointPairs3();
	void getPickUpDeliveryPointPairs (std::vector<int> unvisitedCities, std::vector <int> & bestPair, int whichTour);
	

};






#endif