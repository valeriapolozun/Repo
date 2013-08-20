#ifndef __TOURFIRSTPICKUPSERIELL__
#define __TOURFIRSTPICKUPSERIELL__



#include "OrienteeringProblemWithPickupsAndDeliveries.h"

#include <string>
#include <time.h>

static clock_t clockStart;

class TourGreedyEinzelnePunkteSeriell: public OrienteeringProblemWithPickupsAndDeliveries
{
	public:
	TourGreedyEinzelnePunkteSeriell(std::string inputFile, int selectionPop);
	~TourGreedyEinzelnePunkteSeriell();

	protected:
	bool twoOpt;

	private:
	vector <int> unvisitedNodesForOneTour;
	void calcTourChoosePickupAndDeliveryPointPairs(int whichTour, int selectionPop);
	void calcTourChoosePickupAndDeliveryPointPairs2();
	void pickupPoint(std::vector<int> & nodes);
	bool deliveryPointInserted;
	bool pickupPointInserted;
	std::vector <int> pickupPoints;
	void deliveryPoint (std::vector<int> & nodes);
	std::vector <int> deliveryPoints;
	int pickupInsertedPosition;
	int lastPickupInserted;
	int getNextPickupPoint(std::vector <int> unvisitedPickups, std::vector <int> unvisitedDeliveries, int startNode, vector <int> tour);
	int getNextPickupPointRandomised(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode, vector<int> tour);
	int getNextDeliveryPoint(std::vector <int> unvisitedPickups, std::vector <int> unvisitedDeliveries, int startNode, vector <int> tour);
	int getNextDeliveryPointRandomised(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode, vector <int> tour);
	int getNextPickupPointRandomisedBest15(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode, vector <int> tour, int selectionPop);
	int getNextDeliveryPointRandomisedBest15(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode, vector <int> tour, int selectionPop);
	void putPointInBestPosition(int whichTour, int pointToInsert);
	void putDeliveryPointInBestPosition(int whichTour, int pointToInsert, int pickupInsertedPosition);
	int closestPoint( int node, vector <int> tour);
};


#endif