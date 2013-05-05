#include "TourNearestNeighbor.h"
#include <iostream>


using namespace std;


TourNearestNeighbor::TourNearestNeighbor(string inputFile): OrienteeringProblemWithPickupsAndDeliveries(inputFile)
{
	calcTourNearestNeighbor();
	profitsOfAllTheTours();
}

TourNearestNeighbor::~TourNearestNeighbor()
{
}


void TourNearestNeighbor::calcTourNearestNeighbor(){
	vector<int> unvisitedNodes;

	unvisitedNodes.assign (problemSize, 1);
	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;

	int startNode=0;
	vector<int> tour;
	tour.push_back(0);
	for (int i = 0; i < problemSize-2; i++)
	{
		startNode=getNearestNeighbor(unvisitedNodes, startNode);
		tour.push_back(startNode);
		unvisitedNodes[startNode]=0;
	}
	tour.push_back(1);

	unvisitedNodes.assign (problemSize, 1);
	for (int i=0; i < tour.size(); i++)
	{
		unvisitedNodes[tour[i]]=0;
	}

	for (int i = 0; i < tour.size(); i++)
	{
		cout  << i+1 << ". place in the tour " << tour[i] << endl;
	}
	
}

