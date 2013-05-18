#include "TourPickupDeliveryPointPairsParallel.h"
#include <iostream>


using namespace std;


TourPickupDeliveryPointPairsParallel::TourPickupDeliveryPointPairsParallel(string inputFile): OrienteeringProblemWithPickupsAndDeliveries(inputFile)
{
	bestPairs.assign(2,0);
	unvisitedNodes.assign(problemSize,1);
	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;
	calcTourChoosePickupAndDeliveryPointPairs3();
	
	
	

	profitsOfAllTheTours();
	for (int i=0; i<solutionTours.size();i++)
	{
		cout << "The tour length of the " << i+1 << ". tour is: " << getTourLength(solutionTours[i]) << endl;
	} 



}

TourPickupDeliveryPointPairsParallel::~TourPickupDeliveryPointPairsParallel()
{
}

void TourPickupDeliveryPointPairsParallel::pickUpPointToChoose(vector<int> & nodes)
{
	for (int i=0; i<nodes.size();i++)
	{
		if (inputDataProcessor.getBasicData()[i].quantity<0)
		nodes[i]=0;
	}	
}





void TourPickupDeliveryPointPairsParallel::calcTourChoosePickupAndDeliveryPointPairs3() {
	vector<int> unvisitedNodesForOneTour;

	unvisitedNodesForOneTour= unvisitedNodes;
	double additionalDistance;
	double min;
	int minpos;


	int startNode=0;
	vector<int> tour;
	tour.push_back(0);


	for (int i = 0; i < numberOfTours; i++)
	{
		solutionTours.push_back(tour);
	}


	unvisitedNodesForOneTour[0]=0;
	unvisitedNodesForOneTour[1]=0;
	
	
	for (int i = 0; i < (problemSize-2); i++)  /// TO DO : Calculate how many times it should run (problemSize-2) is wrong
	{
		getPickUpDeliveryPointPairs(unvisitedNodes, bestPairs, i%3);
		

		if (!bestPairs[0]==0)
		{
			if (isTotalLengthUnderLimit2Nodes(solutionTours[i%numberOfTours], bestPairs))
			{
				for (int j = 0; j < bestPairs.size(); j++)
				{
					solutionTours[i%numberOfTours].push_back(bestPairs[j]);
					//unvisitedNodesForOneTour[bestPairs[i]]=0;
					unvisitedNodes[bestPairs[j]]=0;
				}
			}
			else
			{
				wrongPairs[i%numberOfTours].push_back(bestPairs);
			}
			
			bestPairs[0]=0;
			bestPairs[1]=0;
			
		}
	}
	
	for (int i = 0; i < numberOfTours; i++)
	{
		solutionTours[i].push_back(1);
	}
	/*
	for (int i=0; i < tour.size(); i++)
	{
		unvisitedNodes[tour[i]]=0;
	}
	*/
	/*
	if (tour.size()==2)
	{
		unvisitedNodes.assign(problemSize, 0);
	}
	else 
	{
		solutionTours.push_back(tour);

		for (int i = 0; i < tour.size(); i++)
		{
			cout  << i+1 << ". place in the tour " << tour[i] << endl;
		}
	}
	*/

	for (int i=0; i < solutionTours.size(); i++)
	{
		for (int j=0; j < solutionTours[i].size(); j++)
		{
			cout << j+1 << ". place in the " << i+1 << " . tour is: " << solutionTours[i][j] << endl;
		}
	}


}





void TourPickupDeliveryPointPairsParallel::getPickUpDeliveryPointPairs (std::vector<int> unvisitedCities, std::vector <int> & bestPair, int whichTour) 
{
    double max = -DBL_MAX;
    //int best = startNode;
    
	for(int i=0; i<unvisitedCities.size(); i++)
    {
		if(unvisitedCities[i]==0) continue; 
		for(int j=0; j<unvisitedCities.size(); j++)
		{
			if(unvisitedCities[j]==0) continue; 
			getProfitMatrixForPickupAndDeliveryPairsParallel(whichTour);
			//cout << "the i: " << i << "the profitperdistancematrix[startnode][i]" << profitPerDistanceMatrix[startNode][i] << " a max value: " << max << endl;
			if(profitPerDistanceMatrix[i][j]>max)
			{			
            max = profitPerDistanceMatrix[i][j];
			bestPairs[0]=i;
			bestPairs[1]=j;
			}
		 }
    }

	//cout << "The next nearest point is: " << nearest << endl;
   return;
}



