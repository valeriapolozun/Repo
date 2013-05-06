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
			if (i == j) travelDistance = DBL_MAX;    // define the distance from itself is a big number, so that itself can not be the nearest neighbor
			else
			{
				travelDistance = inputDataProcessor.getDistance(i, j);
			}
				if(basicData[i].quantity*basicData[j].quantity >= 0) 
					profitPerDistanceMatrix[i][j]=-DBL_MAX;
				else
				{
				profitPerDistanceMatrix[i][j]= ((basicData[i].profit*basicData[i].quantity)+(basicData[j].profit*basicData[j].quantity))/travelDistance;
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
			/*  TO AVOID INFEASIBLE SOLUTIONS - now not used,...at the end of the improvement I run a function 
			which will reduce the tour until the solution is feasible
			if (getTourLength(tour)>inputDataProcessor.getMaximumTourLength())
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

/*
void OrienteeringProblemWithPickupsAndDeliveries::getTourFeasible(vector <int> solutionTour)
{
	while (getTourLength(solutionTour) > inputDataProcessor.getMaximumTourLength())
	{
		for (int i=1; i<intensity.size()-1;i++)
		{
			if (intensity[i]=0)
			{
				intensity.erase(intensity.begin()+i);
				solutionTour.erase(solutionTour.begin()+i);
				basicDataModified.erase(basicDataModified.begin()+i);
				break;
			}

		}
	}
}
*/


