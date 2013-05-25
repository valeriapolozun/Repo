#include "TourGreedyEinzelnePunkteSeriell.h"
#include <iostream>


using namespace std;


TourGreedyEinzelnePunkteSeriell::TourGreedyEinzelnePunkteSeriell(string inputFile): OrienteeringProblemWithPickupsAndDeliveries(inputFile)
{
	unvisitedNodes.assign(problemSize,1);
	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;

	for (int i=0; i<numberOfTours;i++)
	{
		calcTourChoosePickupAndDeliveryPointPairs();
	} 

	profitsOfAllTheTours();
		for (int i=0; i<solutionTours.size();i++)
	{
		cout << "The tour length of the " << i+1 << ". tour is: " << getTourLength(solutionTours[i]) << endl;
	} 

}

TourGreedyEinzelnePunkteSeriell::~TourGreedyEinzelnePunkteSeriell()
{
}




void TourGreedyEinzelnePunkteSeriell::calcTourChoosePickupAndDeliveryPointPairs(){
	vector<int> unvisitedNodesForOneTour;

	unvisitedNodesForOneTour= unvisitedNodes;
	unvisitedNodesForOneTour[0]=0;
	unvisitedNodesForOneTour[1]=0;
	
	int startNode=0;
	vector<int> tour;
	tour.push_back(0);

	pickupPoint (unvisitedNodesForOneTour);
	deliveryPoint (unvisitedNodesForOneTour);
	int maxRun= max( pickupPoints.size(), deliveryPoints.size());

	for (int i = 0; i < maxRun; i++)
	{
		startNode=getNextPickupPoint( pickupPoints, deliveryPoints, startNode);
		if (isTotalLengthUnderLimit(tour, startNode))
		{
			tour.push_back(startNode);
		}

		//pickupPoints.erase(pickupPoints.begin() + startNode);

		for (int j = 0; j < pickupPoints.size(); j++)
		{
			if (pickupPoints[j]==startNode)
			{
				pickupPoints.erase(pickupPoints.begin()+j);
			}
		}


		startNode=getNextDeliveryPoint(pickupPoints, deliveryPoints, startNode);
	
		if (isTotalLengthUnderLimit(tour, startNode))
		{
			tour.push_back(startNode);
		}

		for (int j = 0; j < deliveryPoints.size(); j++)
		{
			if (startNode== deliveryPoints[j])
			{
				deliveryPoints.erase(deliveryPoints.begin() + j);
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



void TourGreedyEinzelnePunkteSeriell::pickupPoint (std::vector<int> & nodes)
{
	for (int i=0; i<nodes.size();i++)
	{
		if (nodes[i]==1)
		{
			if (inputDataProcessor.getBasicData()[i].quantity<0)
			{
			pickupPoints.push_back(i);
			}
		}
	}
}


void TourGreedyEinzelnePunkteSeriell::deliveryPoint (std::vector<int> & nodes)
{
	for (int i=0; i<nodes.size();i++)
	{
		if (nodes[i]==1)
		{
			if (inputDataProcessor.getBasicData()[i].quantity>0)
			{
			deliveryPoints.push_back(i);
			}
		}
	}
}


int TourGreedyEinzelnePunkteSeriell::getNextPickupPoint(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode)
{
    double min = DBL_MAX;
    int bestToChoose = startNode;
	double distancesToDeliveryPoints=0;
    for(int i=0; i< unvisitedPickups.size(); i++)
    {
         for(int j=0; j< unvisitedDeliveries.size(); j++)
		{
			distancesToDeliveryPoints= distancesToDeliveryPoints + inputDataProcessor.getDistance(unvisitedPickups[i], unvisitedDeliveries[j]);	
		}
	if(distancesToDeliveryPoints<min)
		{
		min =distancesToDeliveryPoints;
		bestToChoose = unvisitedPickups[i];
		}
	distancesToDeliveryPoints=0;
    }
    return bestToChoose;
}


int TourGreedyEinzelnePunkteSeriell::getNextDeliveryPoint(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode)
{
    double min = DBL_MAX;
    int bestToChoose = startNode;
	double distancesToPickupPoints=0;
    for(int i=0; i< unvisitedDeliveries.size(); i++)
    {
         for(int j=0; j< unvisitedPickups.size(); j++)
		{
			distancesToPickupPoints= distancesToPickupPoints + inputDataProcessor.getDistance(unvisitedDeliveries[i], unvisitedPickups[j]);	
		}
	if(distancesToPickupPoints<min)
		{
		min = distancesToPickupPoints;
		bestToChoose = unvisitedDeliveries[i];
		}
	distancesToPickupPoints=0;
    }
    return bestToChoose;
}


