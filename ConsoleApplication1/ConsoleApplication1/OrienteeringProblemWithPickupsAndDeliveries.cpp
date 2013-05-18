#include "OrienteeringProblemWithPickupsAndDeliveries.h"
#include "ProfitCalculator.h"


#include <iostream>
#include <iterator>
#include <glpk.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <time.h>  

using namespace std;

OrienteeringProblemWithPickupsAndDeliveries::OrienteeringProblemWithPickupsAndDeliveries(string inputFile)
{
	inputDataProcessor.init(inputFile);
	problemSize = inputDataProcessor.getProblemSize();
	initProfitPerDistanceMatrix();
	calcPickupDeliveryPointPairs();
	unvisitedNodes.assign (problemSize, 1);
	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;
	numberOfTours=3;
	wrongPairs.resize(numberOfTours);
	vector <int> twotimeszero;
	twotimeszero.assign(2,0);
	for (int i = 0; i < numberOfTours; i++)
	{
		wrongPairs[i].push_back(twotimeszero);
	}
	
}

OrienteeringProblemWithPickupsAndDeliveries::~OrienteeringProblemWithPickupsAndDeliveries()
{
	
}


void OrienteeringProblemWithPickupsAndDeliveries::initProfitPerDistanceMatrix()
{
	vector <double> initToZeroVector(problemSize,0);
	profitPerDistanceMatrix.assign(problemSize, initToZeroVector);
}

void OrienteeringProblemWithPickupsAndDeliveries::calcPickupDeliveryPointPairs()
	{
	double travelDistance; // stores the calculated distances
	double profitValueMax; // profit value per unit which can be reached by visiting a certain pair of Pick up and Delivery point
	vector<Coordinates> basicData(inputDataProcessor.getBasicData());
	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			travelDistance = inputDataProcessor.getDistance(i, j);
			if(basicData[i].quantity*basicData[j].quantity >= 0 && basicData[i].quantity<=0)   
			{
				profitPerDistanceMatrix[i][j]=-DBL_MAX;
			}
			else
			{
			profitPerDistanceMatrix[i][j]= ((basicData[i].profit*basicData[min(i,j)].quantity)+(basicData[j].profit*basicData[min(i,j)].quantity))-travelDistance;
			}
		}
	}

	/*
	cout << "The profitPerDistance matrix is: "<< endl;
	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			cout<<profitPerDistanceMatrix[i][j] << "   " ;
		}
		cout<<endl;
	}
	*/
}


void OrienteeringProblemWithPickupsAndDeliveries::getProfitMatrixForPickupAndDeliveryPairs(int startNode, int whichTour)
	{
	double travelDistance; // stores the calculated distances
	double profitValueMax; // profit value per unit which can be reached by visiting a certain pair of Pick up and Delivery point
	vector<Coordinates> basicData(inputDataProcessor.getBasicData());
	

	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			if(basicData[i].quantity*basicData[j].quantity >= 0 || basicData[i].quantity<=0 || unvisitedNodes[i]==0)   
			{
				profitPerDistanceMatrix[i][j]=-DBL_MAX;
			}
			else
			{
			travelDistance = inputDataProcessor.getDistance(startNode, i)+inputDataProcessor.getDistance(i, j)+inputDataProcessor.getDistance(j, 1);
			profitPerDistanceMatrix[i][j]= ((basicData[i].profit*basicData[min(i,j)].quantity)+(basicData[j].profit*basicData[min(i,j)].quantity))-travelDistance;
			}
		}
	}
	if (!wrongPairs.empty())
	{
	for (int i = 0; i < wrongPairs[whichTour].size(); i++) // Nodes which were previously tried to be added to the tour
	{
		profitPerDistanceMatrix[wrongPairs[whichTour][i][0]] [wrongPairs[whichTour][i][1]]=-DBL_MAX;
	}
	}

	/*
	cout << "The profitPerDistance matrix is: "<< endl;
	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			cout<<profitPerDistanceMatrix[i][j] << "   " ;
		}
		cout<<endl;
	}
	*/
}


void OrienteeringProblemWithPickupsAndDeliveries::getProfitMatrixForPickupAndDeliveryPairsParallel(int whichTour)
	{
	double travelDistance; // stores the calculated distances
	double profitValueMax; // profit value per unit which can be reached by visiting a certain pair of Pick up and Delivery point
	vector<Coordinates> basicData(inputDataProcessor.getBasicData());
	

	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			if(basicData[i].quantity*basicData[j].quantity >= 0 || basicData[i].quantity<=0 || unvisitedNodes[i]==0)   
			{
				profitPerDistanceMatrix[i][j]=-DBL_MAX;
			}
			else
			{
			profitPerDistanceMatrix[i][j]= ((basicData[i].profit*basicData[min(i,j)].quantity)+(basicData[j].profit*basicData[min(i,j)].quantity));
			}
		}
	}

	for (int i = 0; i < wrongPairs[whichTour].size(); i++) // Nodes which were previously tried to be added to the tour
	{
		if ( wrongPairs[whichTour][i][0]!=0 && wrongPairs[whichTour][i][1]!=0)
		{
			profitPerDistanceMatrix[wrongPairs[whichTour][i][0]] [wrongPairs[whichTour][i][1]]=-DBL_MAX;
		}
	}

	/*
	cout << "The profitPerDistance matrix is: "<< endl;
	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			cout<<profitPerDistanceMatrix[i][j] << "   " ;
		}
		cout<<endl;
	}
	*/
}


void OrienteeringProblemWithPickupsAndDeliveries::getProfitMatrixForPickupAndDeliveryPairsParallelBest()
	{
	double travelDistance; // stores the calculated distances
	double profitValueMax; // profit value per unit which can be reached by visiting a certain pair of Pick up and Delivery point
	vector<Coordinates> basicData(inputDataProcessor.getBasicData());
	

	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			if(basicData[i].quantity*basicData[j].quantity >= 0 || basicData[i].quantity<=0 || unvisitedNodes[i]==0)   
			{
				profitPerDistanceMatrix[i][j]=-DBL_MAX;
			}
			else
			{
			profitPerDistanceMatrix[i][j]= ((basicData[i].profit*basicData[min(i,j)].quantity)+(basicData[j].profit*basicData[min(i,j)].quantity));
			}
		}
	}

	/*
	for (int i = 0; i < wrongPairs[whichTour].size(); i++) // Nodes which were previously tried to be added to the tour
	{
		if ( wrongPairs[whichTour][i][0]!=0 && wrongPairs[whichTour][i][1]!=0)
		{
			profitPerDistanceMatrix[wrongPairs[whichTour][i][0]] [wrongPairs[whichTour][i][1]]=-DBL_MAX;
		}
	}

	*/

	/*
	cout << "The profitPerDistance matrix is: "<< endl;
	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			cout<<profitPerDistanceMatrix[i][j] << "   " ;
		}
		cout<<endl;
	}
	*/
}




int OrienteeringProblemWithPickupsAndDeliveries::getNearestNeighbor (vector<int> unvisitedCities, int startNode)
{
    double min = DBL_MAX;
    int nearest = startNode;
    for(int i=0; i<unvisitedCities.size(); i++)
    {
        if(unvisitedCities[i]==0) continue; 
        if (startNode == i) continue; 
        if(inputDataProcessor.getDistance(startNode, i)<min)
        {
            min = inputDataProcessor.getDistance(startNode, i);
            nearest = i;
        }
    }
    return nearest;
}

int OrienteeringProblemWithPickupsAndDeliveries::getHighestDemandInNeighborhood(vector<int> unvisitedCities, int startNode)
{
	double min = DBL_MAX;
    int nearest = startNode;
    for(int i=0; i<unvisitedCities.size(); i++)
    {
        if(unvisitedCities[i]==0) continue; 
        if (startNode == i) continue; 
        if(inputDataProcessor.getDistance(startNode, i)<min)
        {
            min = inputDataProcessor.getDistance(startNode, i);
            nearest = i;
        }
    }
    return nearest;
}





bool OrienteeringProblemWithPickupsAndDeliveries::isTotalLengthUnderLimit(vector<int> currentTour, int nodeToAdd){
	double result = 0;
	int lastNode = currentTour[0];
	currentTour.push_back(nodeToAdd);
	
	for(int i = 1; i < currentTour.size(); i++) { // The total length of the tour will be computed
		int currentNode = currentTour[i];
		result += inputDataProcessor.getDistance(lastNode, currentNode);
		lastNode = currentNode;
	}
	result= result+ inputDataProcessor.getDistance(lastNode, 1);
	return result<=inputDataProcessor.getMaximumTourLength(); // Is my tour length smaller than dMax. TRUE or FALSE will return

}

bool OrienteeringProblemWithPickupsAndDeliveries::isTotalLengthUnderLimit2Nodes(vector<int> currentTour, vector <int> bestPair){
	double result = 0;
	int lastNode = currentTour[0];

	for(int i = 0; i < bestPair.size(); i++) 
	{
		currentTour.push_back(bestPair[i]);
	}

	for(int i = 1; i < currentTour.size(); i++) { // The total length of the tour will be computed
		int currentNode = currentTour[i];
		result += inputDataProcessor.getDistance(lastNode, currentNode);
		lastNode = currentNode;
	}
	result= result+ inputDataProcessor.getDistance(lastNode, 1);
	return result<=inputDataProcessor.getMaximumTourLength(); // Is my tour length smaller than dMax. TRUE or FALSE will return

}


bool OrienteeringProblemWithPickupsAndDeliveries::isThereUnvisitedNodes()
{
	return (find(unvisitedNodes.begin(), unvisitedNodes.end(), 1) != unvisitedNodes.end());
}

void OrienteeringProblemWithPickupsAndDeliveries::profitsOfAllTheTours()
{
	double result=0;
	for (int i=0 ; i< solutionTours.size();i++)
	{
		ProfitCalculator initialToursProfit( solutionTours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity()); 
		cout << "Profit of the " << i+1 << ". tour: " << initialToursProfit.getProfit() << endl;
		result= result+ initialToursProfit.getProfit();
	}
	cout << "The total profit of the tours all together :  " << result << endl;
}

void OrienteeringProblemWithPickupsAndDeliveries::runTwoOpt()
{
	for (int i = 0; i < solutionTours.size(); i++)
	{
		doTwoOpt( solutionTours[i]);
	}
}





/*
void OrienteeringProblemWithPickupsAndDeliveries::doTwoOpt(std::vector <int>  & tour)
{  
	int n, l;
	n= tour.size();
	l=floor((n-2)/2);
	bool criterion=true;
	while (criterion){
		criterion=false;
		for (int i = 0; i < l; i++) // length of sequence
		{
			if (criterion==true)
			{ 
				break;
			}
			for (int j=0; j<n; j++) 
			{
				if ( (inputDataProcessor.getDistance(tour[i+j], tour[n+i-1])+inputDataProcessor.getDistance(tour[i], tour[i+j+1])) <  (inputDataProcessor.getDistance(tour[i], tour[n+i-1])+inputDataProcessor.getDistance(tour[i+j], tour[i+j+1])));
				{
					std:: swap (tour[i], tour[j]);
					criterion=true;
				}

			}
		}
	}
	cout << "The tour after 2-opt: "<< endl;
		for(int j = 0; j < tour.size(); j++)
		{
			cout<<tour[j] << "   " ;
		}
		cout<<endl;
}
*/

/*void OrienteeringProblemWithPickupsAndDeliveries::doTwoOpt(std::vector <int>  & tour)
{  
	int n, l;
	n= tour.size();
	l=floor((n-2)/2);
	bool criterion=true;
	while (criterion){
		criterion=false;
		for (int i = 1; i < l; i++) // 
		{
			if (criterion==true)
			{ 
				break;
			}
	
			for (int j=i+1; j<n-2; j++) 
			{
				if ( (inputDataProcessor.getDistance(tour[i-1], tour[i+j])+inputDataProcessor.getDistance(tour[i], tour[i+j+1])) <  (inputDataProcessor.getDistance(tour[i-1], tour[i])+inputDataProcessor.getDistance(tour[i+j], tour[i+j+1])))
				{
					std:: swap (tour[i], tour[j]);
					criterion=true;
				}

			}
		}
	}
	cout << "The tour after 2-opt: "<< endl;
		for(int j = 0; j < tour.size(); j++)
		{
			cout<<tour[j] << "   " ;
		}
		cout<<endl;
}
*/
void OrienteeringProblemWithPickupsAndDeliveries::doTwoOpt(std::vector <int>  & tour)
{  
	int n;
	n= tour.size();
	bool criterion=true;
	while (criterion)
	{
		criterion=false;
		for (int i = 1; i < n-2; i++) // 
		{
			if (criterion==true)
			{ 
				break;
			}
	
			for (int j=i+1; j<n-1; j++) 
			{
				vector <int> newTourTmp(tour);
				changeOrderPartOfTour(newTourTmp, i, j);
				if (getTourLength(newTourTmp)<getTourLength(tour))
				{
						tour=newTourTmp;
						criterion=true;
						break;
				}
			}
		}
	}

	/*
	cout << "The tour after 2-opt: "<< endl;
		for(int j = 0; j < tour.size(); j++)
		{
			cout<<tour[j] << "   " ;
		}
		cout<<endl;
	*/
}

double OrienteeringProblemWithPickupsAndDeliveries::getTourLength(vector <int> tour)
{
	double tourLength=0;
	for(int i = 0; i < tour.size()-1; i++)
	{
		tourLength=tourLength+ inputDataProcessor.getDistance(tour[i],tour[i+1]);
	}
return tourLength;
}

void OrienteeringProblemWithPickupsAndDeliveries::changeOrderPartOfTour( vector <int> & tour, int from, int to)
{
	while (from<to)
	{
		int tmp = tour[from];
		tour[from]=tour [to];
		tour[to]=tmp;
		from++;
		to--;
	}
}

void  OrienteeringProblemWithPickupsAndDeliveries::shaking(vector <int> & tour, int neighborhoodSize)
{
	srand (time(NULL));
	int numberOfNodesToErase;
	int posToErase;
	int posToInsert;
	vector <int> nodesToBeDisplaced;
	numberOfNodesToErase=(rand() % neighborhoodSize)+1;
	if ((tour.size()-2)>numberOfNodesToErase)
	{
		posToErase=rand()%(tour.size()-numberOfNodesToErase-1)+1;


		for (int i = 0; i < numberOfNodesToErase; i++)
		{
			nodesToBeDisplaced.push_back(tour[posToErase+i]);
		}

		tour.erase(tour.begin()+posToErase, tour.begin()+posToErase+numberOfNodesToErase);
	
		posToInsert=rand()%(tour.size()-1)+1;
	
		for (int i=0; i< nodesToBeDisplaced.size(); i++)
		{
			tour.insert(tour.begin()+posToInsert+i, nodesToBeDisplaced[i]);
			//TO AVOID INFEASIBLE SOLUTIONS - now not used,...at the end of the improvement I run a function 
			//which will reduce the tour until the solution is feasible
			/*if (getTourLength(tour)>inputDataProcessor.getMaximumTourLength())
			{
				tour.erase(tour.begin()+posToInsert+i);
				continue;
			}
			*/
		}
	}
	/*
	cout << "The tour after the shaking process: "<< endl;
	for(int j = 0; j < tour.size(); j++)
	{
		cout<<tour[j] << "   " ;
	}
	cout<<endl;
	*/
}

void OrienteeringProblemWithPickupsAndDeliveries::runShaking()
{
	shaking(solutionTours[0],3);
}

void OrienteeringProblemWithPickupsAndDeliveries::runShakingAndTwoopt()
{
	vector < vector <int>> bestTours;
	bestTours=solutionTours;
	int maxTrials=5;
	int improvement=0;
	for (int i=0; i<solutionTours.size();i++)  // The local search with shaking will be done for each tour 
	{	
		vector <int> tour = solutionTours[i];
		for (int k=1 ; k<11;k++)  // Neighborhood size is increasing
		{

			int trials=0;
			ProfitCalculator bestToursProfit ( bestTours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity());
			while (trials<maxTrials)
			{
				shaking(tour,k);
				makeTourFeasible (tour);
				doTwoOpt(tour);
				ProfitCalculator newTour( tour, inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity());

				if (newTour.getProfit() > bestToursProfit.getProfit())
				{
					bestTours[i]=tour;
					improvement=improvement+1;
				}

				else
				{
					trials=trials+1;
				}
			}	
		}
	}
	for (int i=0 ; i< solutionTours.size();i++)
	{
		ProfitCalculator initialToursProfit( solutionTours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity()); 
		ProfitCalculator bestToursProfit( bestTours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity()); 
		cout << "Profit of the " << i+1 << ". tour: " << initialToursProfit.getProfit() << endl;
		cout << "Profit of the  " << i+1 << ".tour after shaking and 2opt: " << bestToursProfit.getProfit() << endl;
		//getTourFeasible(bestTours[i]);
	}
	cout << "The number of improvements: " << improvement << endl;
}


double OrienteeringProblemWithPickupsAndDeliveries::getObjectiveValue(vector<int> tour)
{
	ProfitCalculator profit(tour, inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity());
	double objectiveValue =profit.getProfit() - getTourLength(tour);
	return objectiveValue;
}


void OrienteeringProblemWithPickupsAndDeliveries::makeTourFeasible(vector <int> & solutionTour)
{
	vector <int> zeroIntensityIndices;
	ProfitCalculator Profit (solutionTour, inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity());
	zeroIntensityIndices=Profit.getZeroIntensityIndices();
	while (getTourLength(solutionTour) > inputDataProcessor.getMaximumTourLength())
	{
		if (zeroIntensityIndices.size()!=0)
		{
			solutionTour.erase(solutionTour.begin()+zeroIntensityIndices.back());
			zeroIntensityIndices.erase( zeroIntensityIndices.end()-1);
		}
		else
		{
			solutionTour.erase(solutionTour.begin()+ rand()%(solutionTour.size()-1)+1);
		}
	}
}



void OrienteeringProblemWithPickupsAndDeliveries::runStringExchangesAndTwoopt()
{
	vector < vector <int>> bestTours;
	bestTours=solutionTours;
	int maxTrials=1;
	int improvement=0;
	int tourNumber;
	for (int i=0; i<solutionTours.size()-1;i++)  // The local search with shaking will be done for each tour 
	{	
		tourNumber=i;
		vector <vector <int> > tours = solutionTours;
		for (int k=1 ; k<4;k++)  // Neighborhood size is increasing
		{

			int trials=0;
		ProfitCalculator bestToursProfit1 ( bestTours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity());
		ProfitCalculator bestToursProfit2 ( bestTours[i+1], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity());
			while (trials<maxTrials)
			{
				stringExchanges(tours, tourNumber, k);
				makeTourFeasible (tours[i+1]);
				doTwoOpt(tours[i+1]);
				ProfitCalculator newTour1( tours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity());
				ProfitCalculator newTour2( tours[i+1], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity());

				if (( newTour2.getProfit() - bestToursProfit2.getProfit())> 0)
				{
					bestTours[i+1]=tours [i+1];
					bestTours[i]=tours[i];
					improvement=improvement+1;
				}

				else
				{
					trials=trials+1;
				}
			}	
		}
	}

	for (int i=0 ; i< solutionTours.size();i++)
	{
		ProfitCalculator initialToursProfit( solutionTours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity()); 
		ProfitCalculator bestToursProfit( bestTours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity()); 
		cout << "Profit of the " << i+1 << ". tour: " << initialToursProfit.getProfit() << endl;
		cout << "Profit of the  " << i+1 << ".tour after shaking and 2opt: " << bestToursProfit.getProfit() << endl;
		//getTourFeasible(bestTours[i]);
	}
	cout << "The number of improvements: " << improvement << endl;
}


void  OrienteeringProblemWithPickupsAndDeliveries::stringExchanges (vector < vector <int>> & tours, int tourNumber, int neighborhoodSize)
{
	srand (time(NULL));
	int numberOfNodesToErase;
	int posToErase;
	int posToInsert;
	vector <int> nodesToBeDisplaced;
	numberOfNodesToErase=(rand() % neighborhoodSize)+1;

	//for (int tourNumber=0; tourNumber < solutionTours.size()-1; tourNumber+= 2) // tourNumber is the index of the tour from the solutionTours from which nodes will be erased
	//{
		if ((tours[tourNumber].size()-2)>numberOfNodesToErase)
		{
			//ProfitCalculator tourInwhichNodesWillBeAddedProfit ( solutionTours[tourNumber+1], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity());
			// 1. Erase
			posToErase=rand()%(tours[tourNumber].size()-numberOfNodesToErase-1)+1;
			for (int i = 0; i < numberOfNodesToErase; i++)
			{
				nodesToBeDisplaced.push_back(tours[tourNumber][posToErase+i]);
			}

			tours[tourNumber].erase(tours[tourNumber].begin()+posToErase, tours[tourNumber].begin()+posToErase+numberOfNodesToErase);
	
			
			// 2. Choosing the best place to insert in the other tour
			double minTourLength=DBL_MAX;
			double TourLengthExtension;
			for (int i = 0; i < tours[tourNumber+1].size()-1; i++)
			{
				TourLengthExtension= getTourLength (nodesToBeDisplaced) + inputDataProcessor.getDistance (tours[tourNumber+1][i], nodesToBeDisplaced[0])+ inputDataProcessor.getDistance (tours[tourNumber+1][i+1], nodesToBeDisplaced.back());
				if (TourLengthExtension < minTourLength)
				{
					minTourLength=TourLengthExtension;
					posToInsert= i+1;
				}
			
			}

			// 3. Cheapest Insertion

			for (int i=0; i< nodesToBeDisplaced.size(); i++)
			{
				tours[tourNumber+1].insert(tours[tourNumber+1].begin()+posToInsert+i, nodesToBeDisplaced[i]);
			}
			//ProfitCalculator tourInwhichNodesAddedNewProfit ( solutionTours[tourNumber+1], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity());
		}
	//}
	/*
	cout << "The tour after the shaking process: "<< endl;
	for(int j = 0; j < tour.size(); j++)
	{
		cout<<tour[j] << "   " ;
	}
	cout<<endl;
	*/
}


