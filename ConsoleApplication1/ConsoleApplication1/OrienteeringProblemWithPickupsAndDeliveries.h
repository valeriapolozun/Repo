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
	string rexe;
	string rpath;
	string filename;
	//string extension;
	
	
	void Rprintsol(string rexe, string rpath, string filename, ProfitCalculator solution, int i);

	void doTwoOpt(int whichTour);
	void erasePoints(int neighborhoodSize, double percentageToErase);
	void insertPoints(vector<Coordinates> basicData, int whichTour);
	bool isZero (int i);
	void doInsertion();
	double maxCapacity;




protected:

	InputDataProcessor inputDataProcessor;
	std::vector <std::vector <int> > solutionTours; // stores the tours
	int getNearestNeighbor(std::vector<int> unvisitedCities, int startNode);
	int problemSize; // the number of nodes
	std::vector<int> unvisitedNodes;
	bool isTotalLengthUnderLimit(std::vector<int> currentTour, int nodeToAdd);
	bool isTotalLengthUnderLimit2Nodes(std::vector<int> currentTour, std::vector <int> nodesToAdd);
	bool isTotalLengthUnderLimit (vector<int> currentTour, int nodeToAdd, int position);
	std::vector < std::vector <double> > profitPerDistanceMatrix;
	void getProfitMatrixForPickupAndDeliveryPairs(int startNode, int whichTour);
	int numberOfTours;
	std::vector <std::vector <std::vector <int> > > wrongPairs;
	std::vector <std::vector <int> > wrongPairsForOneTour;  // Vector which includes the pair which can not be added to one certain tour
	double getTourLength(std::vector <int> tour);
	void getProfitMatrixForPickupAndDeliveryPairsParallel(int whichTour);
	void getProfitMatrixForPickupAndDeliveryPairsParallelBest();
	std::vector <std::vector <int> > load;
	std::vector <std::vector <int> > goodsOnTheLorry;
	
	std::vector <std::vector <int> > bufferPlus;
	std::vector <std::vector <int> >  bufferMinus;
	std::vector <std::vector <double> >  intensity;

private:
	
	void initProfitPerDistanceMatrix();
	void calcPickupDeliveryPointPairs();
	int getHighestDemandInNeighborhood(std::vector<int> unvisitedCities, int startNode);
	
	void changeOrderPartOfTour( std::vector <int> & tour, int from, int to);
	void shaking(std::vector <int> & tour, int neighborhoodSize);
	void makeTourFeasible(std::vector <int> & solutionTour);
	void stringExchanges(std::vector <std::vector <int> > & tours,int tourNumber, int neighborhoodSize);

	


	
};

#endif