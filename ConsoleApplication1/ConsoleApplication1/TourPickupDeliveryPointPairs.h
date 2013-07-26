#ifndef __TOURPICKUP__
#define __TOURPICKUP__


#include "OrienteeringProblemWithPickupsAndDeliveries.h"

#include <string>
#include <vector>

class TourPickupDeliveryPointPairs: public OrienteeringProblemWithPickupsAndDeliveries
{
	public:
	TourPickupDeliveryPointPairs(std::string inputFile, int selectionPop);
	~TourPickupDeliveryPointPairs();

	private:
	void pickUpPointToChoose(std::vector<int> & nodes);
	void calcTourChoosePickupAndDeliveryPointPairs();
	int getPickUpDeliveryPointPairsOnePointAdded (std::vector<int> unvisitedCities, int startNode);
	void getPickUpDeliveryPointPairsTwoPointsAdded (std::vector<int> unvisitedCities, int startNode, std::vector <int> & bestPair, int whichTour);
	std::vector <int> bestPairs;
	void calcTourChoosePickupAndDeliveryPointPairs2(int whichTour, int selectionPop);
	void getPickUpDeliveryPointPairsTwoPointsAddedRandomised (vector<int> unvisitedCities, int startNode, vector <int> & bestPairs, int whichTour);
	void getPickUpDeliveryPointPairsTwoPointsAddedRandomisedBest15 (vector<int> unvisitedCities, int startNode, vector <int> & bestPairs, int whichTour, int selectionPop);

};


#endif
