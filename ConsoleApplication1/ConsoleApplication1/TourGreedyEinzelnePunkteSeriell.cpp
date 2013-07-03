#include "TourGreedyEinzelnePunkteSeriell.h"
#include "InputDataProcessor.h"
#include "LoadCalculator.h"
#include <iostream>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <algorithm>


using namespace std;


TourGreedyEinzelnePunkteSeriell::TourGreedyEinzelnePunkteSeriell(string inputFile): OrienteeringProblemWithPickupsAndDeliveries(inputFile)
{
	unvisitedNodes.assign(problemSize,1);
	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;
	srand(time(NULL));

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
		//startNode=getNextPickupPoint( pickupPoints, deliveryPoints, startNode);
		startNode=getNextPickupPointRandomised( pickupPoints, deliveryPoints, startNode);
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


		//startNode=getNextDeliveryPoint(pickupPoints, deliveryPoints, startNode);
		startNode=getNextDeliveryPointRandomised(pickupPoints, deliveryPoints, startNode);

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
			if (inputDataProcessor.getBasicData()[i].quantity>0)
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
			if (inputDataProcessor.getBasicData()[i].quantity<0)
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

int TourGreedyEinzelnePunkteSeriell::getNextPickupPointRandomised(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode)
{
    double min = DBL_MAX;
    int bestToChoose = startNode;
	double distancesToDeliveryPoints=0;
	double profitValue;
	vector <double> cumulatedDistances;
	cumulatedDistances.push_back(0);
	double randNumber;
    for(int i=0; i< unvisitedPickups.size(); i++)
    {
         for(int j=0; j< unvisitedDeliveries.size(); j++)
		{
			profitValue=profitPerDistanceMatrix [unvisitedPickups[i]][unvisitedDeliveries[j]] /inputDataProcessor.getDistance(startNode, unvisitedPickups[i]);	
			
			if (profitValue>0)
			{
			distancesToDeliveryPoints= distancesToDeliveryPoints + profitValue;
			}
		}
		if (distancesToDeliveryPoints>0)
		{
			cumulatedDistances.push_back(cumulatedDistances.back()+distancesToDeliveryPoints);	
		}
		else
		{
			cumulatedDistances.push_back(cumulatedDistances.back());
		}
		
	
		distancesToDeliveryPoints=0;
    }
   
	double total=cumulatedDistances.back();
	//srand(time(NULL));
	//int proba=rand();
	randNumber= (double) rand() / (double) (RAND_MAX + 1)* total  ;
	for(int k=1; k<cumulatedDistances.size(); k++)
	{
		if (randNumber<=cumulatedDistances[k] && randNumber>cumulatedDistances[k-1])
		{
			bestToChoose = unvisitedPickups[k-1];
			break;
		}
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

int TourGreedyEinzelnePunkteSeriell::getNextDeliveryPointRandomised(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode)
{
   
    int bestToChoose = startNode;
	double profitsToGet=0;
	double quantityToTransfer;
	vector <double> cumulatedDistances;
	double travelDistance;
	cumulatedDistances.push_back(0);
	double randNumber;
    for(int i=0; i< unvisitedDeliveries.size(); i++)
    {
        /*for(int j=0; j< unvisitedPickups.size(); j++)
		{*/
		quantityToTransfer= min(-(inputDataProcessor.basicData[unvisitedDeliveries[i]].quantity), inputDataProcessor.basicData[startNode].quantity);
		profitsToGet= -(inputDataProcessor.basicData[unvisitedDeliveries[i]].profit) * quantityToTransfer;
		travelDistance=inputDataProcessor.getDistance(startNode, unvisitedDeliveries[i]);

		cumulatedDistances.push_back(cumulatedDistances.back()+profitsToGet-travelDistance);
	
		profitsToGet=0;
    }
   
	double total=cumulatedDistances.back();
	randNumber= (double) rand() / (double) (RAND_MAX + 1) * total;
	for(int k=1; k<cumulatedDistances.size(); k++)
	{
		if (randNumber<=cumulatedDistances[k] && randNumber>cumulatedDistances[k-1])
		{
			bestToChoose = unvisitedDeliveries[k-1];
			break;
		}
	}

	return bestToChoose;
}





/*
void TourGreedyEinzelnePunkteSeriell::putPointInBestPosition(int whichTour, int pointToInsert)
{
	double minTourLengthExtension=DBL_MAX;
	double TourLengthExtension;
	int posToInsert;
		
	for (int i = 0; i < solutionTours[whichTour].size()-1; i++)
	{
		if (isTotalLengthUnderLimit( solutionTours[whichTour], pointToInsert, i+1) && (inputDataProcessor.getQuantity(pointToInsert)) < bufferPlus[whichTour][i+1] )
		{
			TourLengthExtension= inputDataProcessor.getDistance (solutionTours[whichTour][i], pointToInsert)+ inputDataProcessor.getDistance (pointToInsert, solutionTours[whichTour][i+1]);
			if (TourLengthExtension < minTourLengthExtension) 
			{
				minTourLengthExtension=TourLengthExtension;
				posToInsert= i+1;
			}
		}
	}
	if (minTourLengthExtension!=DBL_MAX)
	{
		solutionTours[whichTour].insert(solutionTours[whichTour].begin()+posToInsert, pointToInsert);
		intensity[whichTour].insert(intensity[whichTour].begin()+posToInsert, 1);
		LoadCalculator newLoad( load[whichTour], goodsOnTheLorry[whichTour], bufferPlus[whichTour], bufferMinus[whichTour], intensity[whichTour]);
		unvisitedNodes[pointToInsert]=0;
	}
						
	unvisitedNodesForOneTour[pointToInsert]=0;


}


*/