#ifndef __ORIENTEERINGPROBLEMWITHPICKUPSANDDELIVERIES__
#define __ORIENTEERINGPROBLEMWITHPICKUPSANDDELIVERIES__

#include <vector>
#include <string>
#include "InputDataProcessor.h"
#include "Coordinates.h"
#include "ProfitCalculator.h"

class OrienteeringProblemWithPickupsAndDeliveries
{
public:
	OrienteeringProblemWithPickupsAndDeliveries(std::string inputFile);
	~OrienteeringProblemWithPickupsAndDeliveries();

	//void calcOtherToursNearestNeighborDistanceLimit(vector<int> tour);
	
	bool isThereUnvisitedNodes();
	void profitsOfAllTheTours();

	void runNearestNeighborDistanceLimit();
	void runFirstPickupSecondDeliveryPoint();
	void runTwoOpt();
	void runShaking();
	void runShakingAndTwoopt();
	double getObjectiveValue(std::vector <int> tour); // Profit minus Cost
	void runStringExchangesAndTwoopt();

protected:

	InputDataProcessor inputDataProcessor;
	std::vector <std::vector <int> > solutionTours; // stores the tours
	int getNearestNeighbor(std::vector<int> unvisitedCities, int startNode);
	int problemSize; // the number of nodes
	std::vector<int> unvisitedNodes;
	bool isTotalLengthUnderLimit(std::vector<int> currentTour, int nodeToAdd);
	bool isTotalLengthUnderLimit2Nodes(std::vector<int> currentTour, std::vector <int> nodesToAdd);
	std::vector < std::vector <double> > profitPerDistanceMatrix;
	void getProfitMatrixForPickupAndDeliveryPairs(int startNode);
	int numberOfTours;
	std::vector <std::vector <int> > wrongPairs; 
	double getTourLength(std::vector <int> tour);

private:
	
	void initProfitPerDistanceMatrix();
	void calcPickupDeliveryPointPairs();
	int getHighestDemandInNeighborhood(std::vector<int> unvisitedCities, int startNode);
	void doTwoOpt(std::vector <int>  & tour);
	void changeOrderPartOfTour( std::vector <int> & tour, int from, int to);
	void shaking(std::vector <int> & tour, int neighborhoodSize);
	void makeTourFeasible(std::vector <int> & solutionTour);
	void stringExchanges(std::vector <std::vector <int> > & tours,int tourNumber, int neighborhoodSize);


	
};

#endif