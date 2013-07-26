#include "TourNearestNeighborDistanceLimit.h"
#include <iostream>


using namespace std;


TourNearestNeighborDistanceLimit::TourNearestNeighborDistanceLimit(string inputFile): OrienteeringProblemWithPickupsAndDeliveries(inputFile)
{
	unvisitedNodes.assign(problemSize,1);
	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;
	while (isThereUnvisitedNodes())
	{
		calcTourNearestNeighborDistanceLimit();
	} 

	profitsOfAllTheTours(1,0);
}

TourNearestNeighborDistanceLimit::~TourNearestNeighborDistanceLimit()
{
}


void TourNearestNeighborDistanceLimit::calcTourNearestNeighborDistanceLimit(){
	vector<int> unvisitedNodesForOneTour;

	unvisitedNodesForOneTour= unvisitedNodes;
	
	int startNode=0;
	vector<int> tour;
	tour.push_back(0);
	for (int i = 0; i < problemSize-2; i++)
	{
		startNode=getNearestNeighbor(unvisitedNodesForOneTour, startNode);
		if (isTotalLengthUnderLimit(tour, startNode))
		{
			tour.push_back(startNode);
		}
		unvisitedNodesForOneTour[startNode]=0;
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

		cout << "The " << solutionTours.size() << ". tour: "<< endl;
		for (int i = 0; i < tour.size(); i++)
		{
			cout  << i+1 << ". place in the tour " << tour[i] << endl;
		}
	}
}
