#ifndef __TOURPICKUPDELIVERYPARALLEL__
#define __TOURPICKUPDELIVERYPARALLEL__


#include "OrienteeringProblemWithPickupsAndDeliveries.h"

#include <string>
#include <vector>

typedef std::pair <double, vector <int>> mypair2;

class TourPickupDeliveryPointPairsParallel: public OrienteeringProblemWithPickupsAndDeliveries
{
	public:
	TourPickupDeliveryPointPairsParallel(std::string inputFile, int selectionPop);
	~TourPickupDeliveryPointPairsParallel();

	private:
	bool twoOpt;
	void pickUpPointToChoose(std::vector<int> & nodes);
	void calcTourChoosePickupAndDeliveryPointPairs();
	void getPickUpDeliveryPointPairsTwoPointsAdded (std::vector<int> unvisitedCities, int startNode, std::vector <int> & bestPair);
	std::vector <int> bestPairs;
	void calcTourChoosePickupAndDeliveryPointPairs3(int selectionPop);
	void getPickUpDeliveryPointPairs (std::vector<int> unvisitedCities, std::vector <int> & bestPair, int whichTour);
	void getPickUpDeliveryPointPairsTwoPointsAddedRandomised (vector<int> unvisitedCities, int startNode, vector <int> & bestPairs, int whichTour, int selectionPop);
	vector <vector <mypair2>> probabilities;

};






#endif