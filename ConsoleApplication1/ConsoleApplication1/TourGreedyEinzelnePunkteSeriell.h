#ifndef __TOURFIRSTPICKUPSERIELL__
#define __TOURFIRSTPICKUPSERIELL__



#include "OrienteeringProblemWithPickupsAndDeliveries.h"

#include <string>
#include <time.h>

static clock_t clockStart;

class TourGreedyEinzelnePunkteSeriell: public OrienteeringProblemWithPickupsAndDeliveries
{
	public:
	TourGreedyEinzelnePunkteSeriell(std::string inputFile);
	~TourGreedyEinzelnePunkteSeriell();

	private:
	vector <int> unvisitedNodesForOneTour;
	void calcTourChoosePickupAndDeliveryPointPairs(int whichTour);
	void calcTourChoosePickupAndDeliveryPointPairs2();
	void pickupPoint(std::vector<int> & nodes);
	bool deliveryPointInserted;
	std::vector <int> pickupPoints;
	void deliveryPoint (std::vector<int> & nodes);
	std::vector <int> deliveryPoints;
	int pickupInsertedPosition;
	int lastPickupInserted;
	int getNextPickupPoint(std::vector <int> unvisitedPickups, std::vector <int> unvisitedDeliveries, int startNode);
	int getNextPickupPointRandomised(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode, vector<int> tour);
	int getNextDeliveryPoint(std::vector <int> unvisitedPickups, std::vector <int> unvisitedDeliveries, int startNode);
	int getNextDeliveryPointRandomised(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode, vector <int> tour);
	void putPointInBestPosition(int whichTour, int pointToInsert);
	void putDeliveryPointInBestPosition(int whichTour, int pointToInsert, int pickupInsertedPosition);
	int closestPoint( int node, vector <int> tour);
};


#endif