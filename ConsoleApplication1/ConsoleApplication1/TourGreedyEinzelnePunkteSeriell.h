#ifndef __TOURFIRSTPICKUPSERIELL__
#define __TOURFIRSTPICKUPSERIELL__


#include "OrienteeringProblemWithPickupsAndDeliveries.h"

#include <string>
#include <vector>

class TourGreedyEinzelnePunkteSeriell: public OrienteeringProblemWithPickupsAndDeliveries
{
	public:
	TourGreedyEinzelnePunkteSeriell(std::string inputFile);
	~TourGreedyEinzelnePunkteSeriell();

	private:
	void calcTourChoosePickupAndDeliveryPointPairs();
	void calcTourChoosePickupAndDeliveryPointPairs2();
	void pickupPoint(std::vector<int> & nodes);
	std::vector <int> pickupPoints;
	void deliveryPoint (std::vector<int> & nodes);
	std::vector <int> deliveryPoints;
	int getNextPickupPoint(std::vector <int> unvisitedPickups, std::vector <int> unvisitedDeliveries, int startNode);
	int getNextPickupPointRandomised(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode);
	int getNextDeliveryPoint(std::vector <int> unvisitedPickups, std::vector <int> unvisitedDeliveries, int startNode);
	int getNextDeliveryPointRandomised(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode);
	void TourGreedyEinzelnePunkteSeriell::putPointInBestPosition(int whichTour, int pointToInsert);
};


#endif