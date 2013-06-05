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


