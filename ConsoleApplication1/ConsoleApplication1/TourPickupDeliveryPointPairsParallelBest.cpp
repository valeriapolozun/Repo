#include "TourPickupDeliveryPointPairsParallelBest.h"
#include <iostream>


using namespace std;


TourPickupDeliveryPointPairsParallelBest::TourPickupDeliveryPointPairsParallelBest(string inputFile): OrienteeringProblemWithPickupsAndDeliveries(inputFile)
{
	bestPairs.assign(2,0);
	unvisitedNodes.assign(problemSize,1);
	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;
	calcTourChoosePickupAndDeliveryPointPairs2();
	

	profitsOfAllTheTours();
	for (int i=0; i<solutionTours.size();i++)
	{
		cout << "The tour length of the " << i+1 << ". tour is: " << getTourLength(solutionTours[i]) << endl;
	} 



}

TourPickupDeliveryPointPairsParallelBest::~TourPickupDeliveryPointPairsParallelBest()
{
}

void TourPickupDeliveryPointPairsParallelBest::pickUpPointToChoose(vector<int> & nodes)
{
	for (int i=0; i<nodes.size();i++)
	{
		if (inputDataProcessor.getBasicData()[i].quantity<0)
		nodes[i]=0;
	}	
}



void TourPickupDeliveryPointPairsParallelBest::calcTourChoosePickupAndDeliveryPointPairs(){
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



void TourPickupDeliveryPointPairsParallelBest::calcTourChoosePickupAndDeliveryPointPairs2(){
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
		//getPickUpDeliveryPointPairsTwoPointsAdded(unvisitedNodes, tour.back(), bestPairs);

	/*
		startNode=getPickUpDeliveryPointPairsOnePointAdded(unvisitedNodesForOneTour, tour.back());
		if (startNode==0){
			continue;
		}
		else
		{
	*/
		getPickUpDeliveryPointPairs(unvisitedNodes, bestPairs);

		if (!bestPairs[0]==0)
		{
			min=DBL_MAX;
			for (int i = 0; i < numberOfTours; i++)
			{
				
				additionalDistance=inputDataProcessor.getDistance( solutionTours[i].back(), bestPairs[0]) + inputDataProcessor.getDistance( bestPairs[1], 1);
				if (additionalDistance <min)
				{
					min=additionalDistance;
					minpos= i;
				}
			}

			if (isTotalLengthUnderLimit2Nodes(solutionTours[minpos], bestPairs))
			{
				for (int i = 0; i < bestPairs.size(); i++)
				{
					solutionTours[minpos].push_back(bestPairs[i]);
					//unvisitedNodesForOneTour[bestPairs[i]]=0;
					unvisitedNodes[bestPairs[i]]=0;
				}
			}
			else
			{
			wrongPairs[minpos].push_back(bestPairs);
			}
			bestPairs[0]=0;
			bestPairs[1]=0;
			;
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
}



int TourPickupDeliveryPointPairsParallelBest::getPickUpDeliveryPointPairsOnePointAdded (vector<int> unvisitedCities, int startNode)
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


void TourPickupDeliveryPointPairsParallelBest::getPickUpDeliveryPointPairsTwoPointsAdded (vector<int> unvisitedCities, int startNode, vector <int> & bestPairs, int whichTour)
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


void TourPickupDeliveryPointPairsParallelBest::getPickUpDeliveryPointPairs (std::vector<int> unvisitedCities, std::vector <int> & bestPair) 
{
    double max = -DBL_MAX;
    //int best = startNode;
    
	for(int i=0; i<unvisitedCities.size(); i++)
    {
		if(unvisitedCities[i]==0) continue; 
		for(int j=0; j<unvisitedCities.size(); j++)
		{
			if(unvisitedCities[j]==0) continue; 
			getProfitMatrixForPickupAndDeliveryPairsParallelBest();
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
