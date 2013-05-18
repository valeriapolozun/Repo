#include "TourPickupDeliveryPointPairs.h"
#include <iostream>


using namespace std;


TourPickupDeliveryPointPairs::TourPickupDeliveryPointPairs(string inputFile): OrienteeringProblemWithPickupsAndDeliveries(inputFile)
{
	bestPairs.assign(2,0);
	unvisitedNodes.assign(problemSize,1);
	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;
	for (int i=0; i<numberOfTours;i++)
	{
		calcTourChoosePickupAndDeliveryPointPairs2(i);
	} 

	profitsOfAllTheTours();
	for (int i=0; i<solutionTours.size();i++)
	{
		cout << "The tour length of the " << i+1 << ". tour is: " << getTourLength(solutionTours[i]) << endl;
	} 



}

TourPickupDeliveryPointPairs::~TourPickupDeliveryPointPairs()
{
}

void TourPickupDeliveryPointPairs::pickUpPointToChoose(vector<int> & nodes)
{
	for (int i=0; i<nodes.size();i++)
	{
		if (inputDataProcessor.getBasicData()[i].quantity<0)
		nodes[i]=0;
	}	
}



void TourPickupDeliveryPointPairs::calcTourChoosePickupAndDeliveryPointPairs(){
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



void TourPickupDeliveryPointPairs::calcTourChoosePickupAndDeliveryPointPairs2(int whichTour){
	vector<int> unvisitedNodesForOneTour;

	unvisitedNodesForOneTour= unvisitedNodes;
	
	int startNode=0;
	vector<int> tour;
	tour.push_back(0);

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
	*/

	for (int i = 0; i < (problemSize-2); i++)  /// TO DO : Calculate how many times it should run (problemSize-2) is wrong
	{
		getPickUpDeliveryPointPairsTwoPointsAdded(unvisitedNodesForOneTour, tour.back(), bestPairs, whichTour);

	/*
		startNode=getPickUpDeliveryPointPairsOnePointAdded(unvisitedNodesForOneTour, tour.back());
		if (startNode==0){
			continue;
		}
		else
		{
	*/
		if (!bestPairs[0]==0)
		{
			if (isTotalLengthUnderLimit2Nodes(tour, bestPairs))
			{
				for (int i = 0; i < bestPairs.size(); i++)
				{
					tour.push_back(bestPairs[i]);
					unvisitedNodesForOneTour[bestPairs[i]]=0;
					unvisitedNodes[bestPairs[i]]=0;
				}
			}
			else
			{
			wrongPairsForOneTour.push_back(bestPairs);
			}
			bestPairs[0]=0;
			bestPairs[1]=0;
		}
	}
	
	tour.push_back(1);

	/*
	for (int i=0; i < tour.size(); i++)
	{
		unvisitedNodes[tour[i]]=0;
	}
	*/

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



int TourPickupDeliveryPointPairs::getPickUpDeliveryPointPairsOnePointAdded (vector<int> unvisitedCities, int startNode)
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


void TourPickupDeliveryPointPairs::getPickUpDeliveryPointPairsTwoPointsAdded (vector<int> unvisitedCities, int startNode, vector <int> & bestPairs, int whichTour)
{
    double max = -DBL_MAX;
    //int best = startNode;
    
	for(int i=0; i<unvisitedCities.size(); i++)
    {
		if(unvisitedCities[i]==0) continue; 
		for(int j=0; j<unvisitedCities.size(); j++)
		{
			if(unvisitedCities[j]==0) continue; 
			if (startNode == i) continue; 
			if (startNode == j) continue; 
			getProfitMatrixForPickupAndDeliveryPairs (startNode, whichTour);
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

