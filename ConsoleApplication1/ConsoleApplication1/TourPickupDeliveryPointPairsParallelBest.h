#ifndef __TOURPICKUPPARALLEL__
#define __TOURPICKUPPARALLEL__


#include "OrienteeringProblemWithPickupsAndDeliveries.h"

#include <string>
#include <vector>

class TourPickupDeliveryPointPairsParallelBest: public OrienteeringProblemWithPickupsAndDeliveries
{
	public:
	TourPickupDeliveryPointPairsParallelBest(std::string inputFile);
	~TourPickupDeliveryPointPairsParallelBest();

	private:
	void pickUpPointToChoose(std::vector<int> & nodes);
	void calcTourChoosePickupAndDeliveryPointPairs();
	int getPickUpDeliveryPointPairsOnePointAdded (std::vector<int> unvisitedCities, int startNode);
	void getPickUpDeliveryPointPairsTwoPointsAdded (std::vector<int> unvisitedCities, int startNode, std::vector <int> & bestPair, int whichTour);
	std::vector <int> bestPairs;
	void calcTourChoosePickupAndDeliveryPointPairs2();
	void getPickUpDeliveryPointPairs (std::vector<int> unvisitedCities, std::vector <int> & bestPair);
	

};






#endif