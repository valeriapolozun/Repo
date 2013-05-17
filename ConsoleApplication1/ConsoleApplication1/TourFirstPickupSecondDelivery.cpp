#include "TourFirstPickupSecondDelivery.h"
#include <iostream>


using namespace std;


TourFirstPickupSecondDelivery::TourFirstPickupSecondDelivery(string inputFile): OrienteeringProblemWithPickupsAndDeliveries(inputFile)
{
	unvisitedNodes.assign(problemSize,1);
	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;
	while (isThereUnvisitedNodes())
	{
		calcTourChoosePickupAndDeliveryPointPairs();
	} 

	profitsOfAllTheTours();
}

TourFirstPickupSecondDelivery::~TourFirstPickupSecondDelivery()
{
}

void TourFirstPickupSecondDelivery::pickUpPointToChoose(vector<int> & nodes)
{
	for (int i=0; i<nodes.size();i++)
	{
		if (inputDataProcessor.getBasicData()[i].quantity<0)
		nodes[i]=0;
	}	
}



void TourFirstPickupSecondDelivery::calcTourChoosePickupAndDeliveryPointPairs(){
	vector<int> unvisitedNodesForOneTour;

	unvisitedNodesForOneTour= unvisitedNodes;
	
	int startNode=0;
	vector<int> tour;
	tour.push_back(0);
	
	pickUpPointToChoose(unvisitedNodesForOneTour);
	
	startNode= getNearestNeighbor(unvisitedNodesForOneTour, startNode);
	if (startNode!=0){	
		tour.push_back(startNode);
		
		unvisitedNodesForOneTour= unvisitedNodes;
		unvisitedNodesForOneTour[0]=0;
		unvisitedNodesForOneTour[1]=0;
		unvisitedNodesForOneTour[startNode]=0;

	for (int i = 0; i < problemSize-3; i++)
	{
		startNode=getPickUpDeliveryPointPairsOnePointAdded(unvisitedNodesForOneTour, tour.back());
		if (startNode==0){
			continue;
		}
		else
		{
			if (isTotalLengthUnderLimit(tour, startNode)){
			tour.push_back(startNode);
			}
		unvisitedNodesForOneTour[startNode]=0;
		}
	}
	}
	tour.push_back(1);

	
	for (int i=0; i < tour.size(); i++)
	{
		unvisitedNodes[tour[i]]=0;
	}

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
	
}


/*
void TourFirstPickupSecondDelivery::calcTourChoosePickupAndDeliveryPointPairs2(){
	vector<int> unvisitedNodesForOneTour;

	unvisitedNodesForOneTour= unvisitedNodes;
	
	int startNode=0;
	vector<int> tour;
	tour.push_back(0);

	unvisitedNodesForOneTour= unvisitedNodes;
	unvisitedNodesForOneTour[0]=0;
	unvisitedNodesForOneTour[1]=0;



	/*
	pickUpPointToChoose(unvisitedNodesForOneTour);
	
	startNode= getNearestNeighbor(unvisitedNodesForOneTour, startNode);
	if (startNode!=0){	
		tour.push_back(startNode);
		
		unvisitedNodesForOneTour= unvisitedNodes;
		unvisitedNodesForOneTour[0]=0;
		unvisitedNodesForOneTour[1]=0;
		unvisitedNodesForOneTour[startNode]=0;


	for (int i = 0; i < (problemSize-2); i++)  /// TO DO : Calculate how many times it should run (problemSize-2) is wrong
	{
		getPickUpDeliveryPointPairsTwoPointsAdded(unvisitedNodesForOneTour, tour.back(), bestPair);

	/*
		startNode=getPickUpDeliveryPointPairsOnePointAdded(unvisitedNodesForOneTour, tour.back());
		if (startNode==0){
			continue;
		}
		else
		{
	
		if (!bestPair.empty())
		{
			if (isTotalLengthUnderLimit2Nodes(tour, bestPair))
			{
				for (int i = 0; i < bestPair.size(); i++)
				{
					tour.push_back(bestPair[i]);
				}
			}
			unvisitedNodesForOneTour[bestPair[i]]=0;
			bestPair.erase (bestPair.begin(),bestPair.begin()+1);
		}
	}
	
	tour.push_back(1);

	
	for (int i=0; i < tour.size(); i++)
	{
		unvisitedNodes[tour[i]]=0;
	}

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
	
}

*/

int TourFirstPickupSecondDelivery::getPickUpDeliveryPointPairsOnePointAdded (vector<int> unvisitedCities, int startNode)
{
    double max = -DBL_MAX;
    int best = startNode;
    for(int i=0; i<unvisitedCities.size(); i++)
    {
        if(unvisitedCities[i]==0) continue; 
        if (startNode == i) continue; 
        //cout << "the i: " << i << "the profitperdistancematrix[startnode][i]" << profitPerDistanceMatrix[startNode][i] << " a max value: " << max << endl;
		if(profitPerDistanceMatrix[startNode][i]>max)
        {			
            max = profitPerDistanceMatrix[startNode][i];
            best = i;
        }
    }
	//cout << "The next nearest point is: " << nearest << endl;
    return best;
}


void TourFirstPickupSecondDelivery::getPickUpDeliveryPointPairsTwoPointsAdded (vector<int> unvisitedCities, int startNode, vector <int> & bestPair)
{
    double max = -DBL_MAX;
    //int best = startNode;
    
	for(int i=0; i<unvisitedCities.size(); i++)
    {
		for(int j=0; j<unvisitedCities.size(); j++)
		{
			if(unvisitedCities[i]==0) continue; 
			if (startNode == i) continue; 
			if (startNode == j) continue; 
			getProfitMatrixForPickupAndDeliveryPairs (startNode);
			//cout << "the i: " << i << "the profitperdistancematrix[startnode][i]" << profitPerDistanceMatrix[startNode][i] << " a max value: " << max << endl;
			if(profitPerDistanceMatrix[i][j]>max)
			{			
            max = profitPerDistanceMatrix[i][j];
			bestPair[0]=i;
			bestPair[1]=j;
			}
		 }
    }

	//cout << "The next nearest point is: " << nearest << endl;
   return;
}

