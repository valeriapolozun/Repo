#ifndef __ORIENTEERINGPROBLEMWITHPICKUPSANDDELIVERIES__
#define __ORIENTEERINGPROBLEMWITHPICKUPSANDDELIVERIES__

#include <vector>
#include <string>
#include "InputDataProcessor.h"
#include "Coordinates.h"
#include "ProfitCalculator.h"
#include "ProfitCalculatorOhneGLPK.h"
#include <time.h>  

class OrienteeringProblemWithPickupsAndDeliveries
{
public:

	clock_t clock_start;

	OrienteeringProblemWithPickupsAndDeliveries(std::string inputFile);
	~OrienteeringProblemWithPickupsAndDeliveries();

	//void calcOtherToursNearestNeighborDistanceLimit(vector<int> tour);
	

	
	bool isThereUnvisitedNodes();
	void profitsOfAllTheTours(int seedNumber, double timeStart);
	void profitsOfAllTheToursOhneGLPK();
	void printSolutions();

	void runNearestNeighborDistanceLimit();
	void runFirstPickupSecondDeliveryPoint();
	void runTwoOpt(int seedNumber, double timeStart);
	void runShaking();
	void runShakingAndTwoopt();
	double getObjectiveValue(std::vector <int> tour); // Profit minus Cost
	void runStringExchangesAndTwoopt();
	string rexe;
	string rpath;
	string filename;
	std::vector <std::vector <double>> totalFinalSolutions;
	void runExcelExporter();
	void runExcelExport(string inputFile, string heurName); // csv file
	void runExcelExportStart();
	void runExcelExportFinish();
	//string extension;
	
	
	void Rprintsol(string rexe, string rpath, string filename, ProfitCalculator solution, int i, int & count);
	void Rprintsol2(string rexe, string rpath, string filename, ProfitCalculatorOhneGLPK solution, int i, int & countSolutionRuns);

	void doTwoOpt(int whichTour);
	void erasePoints(int neighborhoodSize, double percentageToErase);
	void insertPoints(vector<Coordinates> basicData, int whichTour);
	bool isZero (int i);
	void doInsertion(int seedNumber, double timeStart);
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
	bool isTotalLengthUnderLimit2NodesDifferentPlace (vector<int> currentTour, int nodeToAdd, int position, int node2ToAdd, int position2);
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
	vector <int> quantities;
	void calcPickupDeliveryPointPairs();
	vector <double> finalSolutions;
	
	int countSolutionRuns;
	int getNextPickupPointRandomised(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode, vector <int> tour);
	int getNextDeliveryPointRandomised(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode,  vector <int> tour);
	void pickupPoint (std::vector<int> & nodes);
	void deliveryPoint (std::vector<int> & nodes);
	vector <int> pickupPoints;
	vector <int> deliveryPoints;
	int closestPoint( int node, vector <int> tour);


private:
	
	void initProfitPerDistanceMatrix();
	void initIntensityMatrix();
	int getHighestDemandInNeighborhood(std::vector<int> unvisitedCities, int startNode);
	
	void changeOrderPartOfTour( std::vector <int> & tour, int from, int to);
	void changeOrderPartOfTourDouble ( vector <double> & tour, int from, int to);
	void shaking(std::vector <int> & tour, int neighborhoodSize);
	void makeTourFeasible(std::vector <int> & solutionTour);
	void stringExchanges(std::vector <std::vector <int> > & tours,int tourNumber, int neighborhoodSize);

	


	
};

#endif