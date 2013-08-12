#include "TourGreedyEinzelnePunkteSeriell.h"
#include "OrienteeringProblemWithPickupsAndDeliveries.h"
#include "InputDataProcessor.h"
#include "LoadCalculator.h"
#include <iostream>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <algorithm>
#include <fstream>


using namespace std;

typedef std::pair <double, int> mypair;
bool comparator (const mypair& l, const mypair& r )
{ return l.first> r.first;}


TourGreedyEinzelnePunkteSeriell::TourGreedyEinzelnePunkteSeriell(string inputFile, int selectionPop): OrienteeringProblemWithPickupsAndDeliveries(inputFile)
{

	
	for (int seedNumber=0; seedNumber<100; seedNumber++) // 100 seed run
	{
	//int seedNumber=44;
		solutionTours.clear();
		unvisitedNodes.assign(problemSize,1);
		unvisitedNodes[0]=0;
		unvisitedNodes[1]=0;
		srand(seedNumber);
		deliveryPointInserted=false;
		double timeStart=clock();
		//double timeStart=0;
		getProfitMatrixForPickupAndDeliveryPairsParallel(0);
		for (int i=0; i<numberOfTours;i++)
		{
			calcTourChoosePickupAndDeliveryPointPairs(i, selectionPop);
		} 
		profitsOfAllTheTours(seedNumber, timeStart);
		
	runExcelExport(inputFile, "heurSeriellPoint" + std::to_string(selectionPop));
	runTwoOpt(seedNumber, timeStart);
	runExcelExport(inputFile, "heurSeriellPoint" + std::to_string(selectionPop) + "+2opt");
	doInsertion(seedNumber, timeStart);
	runExcelExport(inputFile, "heurSeriellPoint" + std::to_string(selectionPop) + "+2opt+insertion");

		for (int i=0; i<solutionTours.size();i++)
		{
			cout << "The tour length of the " << i+1 << ". tour is: " << getTourLength(solutionTours[i]) << endl;
		} 
		if (selectionPop==1)
		{
		break;
		}


	}
	

}

TourGreedyEinzelnePunkteSeriell::~TourGreedyEinzelnePunkteSeriell()
{
}




void TourGreedyEinzelnePunkteSeriell::calcTourChoosePickupAndDeliveryPointPairs(int whichTour, int selectionPop){
	//vector<int> unvisitedNodesForOneTour;

	unvisitedNodesForOneTour= unvisitedNodes;
	unvisitedNodesForOneTour[0]=0;
	unvisitedNodesForOneTour[1]=0;
	
	int startNode=0;
	vector<int> tour;
	tour.push_back(0);
	tour.push_back(1);
	solutionTours.push_back(tour);

	pickupPoint (unvisitedNodesForOneTour);
	deliveryPoint (unvisitedNodesForOneTour);
	int maxRun= max( pickupPoints.size(), deliveryPoints.size());

	for (int i = 0; i < maxRun; i++)
	{

		pickupPointInserted=false;
		while (pickupPointInserted!=true && pickupPoints.size()!=0)//(std::find(unvisitedNodesForOneTour.begin(), unvisitedNodesForOneTour.end(), 1)!=unvisitedNodesForOneTour.end()))
		{
			startNode=0;
			if (selectionPop==1)
			{
			startNode=getNextPickupPoint(pickupPoints, deliveryPoints, startNode, solutionTours[whichTour]);
			}
			else
			{
				if (selectionPop==3 || selectionPop==15)
				{
				startNode=getNextPickupPointRandomisedBest15(pickupPoints, deliveryPoints, lastPickupInserted, solutionTours[whichTour], selectionPop);
				}
				else
				{
				startNode=getNextPickupPointRandomised(pickupPoints, deliveryPoints, lastPickupInserted, solutionTours[whichTour]);
				}
			}

			/*
			//startNode=getNextPickupPoint( pickupPoints, deliveryPoints, startNode, solutionTours[whichTour]);
			//startNode=getNextPickupPointRandomised( pickupPoints, deliveryPoints, *(solutionTours[whichTour].end()-2), solutionTours[whichTour]);
			startNode=getNextPickupPointRandomisedBest15(pickupPoints, deliveryPoints, *(solutionTours[whichTour].end()-2), solutionTours[whichTour], 15);
			*/

			if (startNode!=0)
			{
				putPointInBestPosition( whichTour, startNode);

				unvisitedNodesForOneTour[startNode]=0;

				if (startNode==1 || startNode==0)
				{
				cout << "error";
				}

				/*
				if (isTotalLengthUnderLimit(tour, startNode))
				{
					tour.push_back(startNode);
				}
				*/
				//pickupPoints.erase(pickupPoints.begin() + startNode);

				for (int j = 0; j < pickupPoints.size(); j++)
				{
					if (pickupPoints[j]==startNode)
					{
						pickupPoints.erase(pickupPoints.begin()+j);
					}
				}
			}
			else
			{
				pickupPoints.clear();
			}
		
		}

		deliveryPointInserted=false;
		while (deliveryPointInserted!=true && deliveryPoints.size()!=0)//(std::find(unvisitedNodesForOneTour.begin(), unvisitedNodesForOneTour.end(), 1)!=unvisitedNodesForOneTour.end()))
		{
			startNode=0;


			if (selectionPop==1)
			{
			startNode=getNextDeliveryPoint(pickupPoints, deliveryPoints, startNode, solutionTours[whichTour]);
			}
			else
			{
				if (selectionPop==3 || selectionPop==15)
				{
				startNode=getNextDeliveryPointRandomisedBest15(pickupPoints, deliveryPoints, lastPickupInserted, solutionTours[whichTour], selectionPop);
				}
				else
				{
				startNode=getNextDeliveryPointRandomised(pickupPoints, deliveryPoints, lastPickupInserted, solutionTours[whichTour]);
				}
			}


			//putPointInBestPosition( whichTour, startNode);

			if (startNode==1 || startNode==0)
			{
			cout << "error";
			}

			if (startNode!=0)
			{
				putDeliveryPointInBestPosition(whichTour, startNode, pickupInsertedPosition);
				unvisitedNodesForOneTour[startNode]=0;

				for (int j = 0; j < deliveryPoints.size(); j++)
				{
					if (startNode== deliveryPoints[j])
					{
						deliveryPoints.erase(deliveryPoints.begin() + j);
					}
				}
			}
			else // no more delivery can be inserted!
			{
			deliveryPoints.clear();
			}

		//unvisitedNodesForOneTour[startNode]=0;
		}
		/*
		if (isTotalLengthUnderLimit(tour, startNode))
		{
			tour.push_back(startNode);
		}
		*/


		
		
	}

	//solutionTours[whichTour].push_back(1);

	for (int i=0; i < solutionTours[whichTour].size(); i++)
	{
		unvisitedNodes[solutionTours[whichTour][i]]=0;
	}

	if (solutionTours[whichTour].size()==2)
	{
		unvisitedNodes.assign(problemSize, 0);
	}
	/*
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



void TourGreedyEinzelnePunkteSeriell::pickupPoint (std::vector<int> & nodes)
{
	pickupPoints.clear();
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
	deliveryPoints.clear();
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


int TourGreedyEinzelnePunkteSeriell::getNextPickupPoint(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode, vector <int> tour)
{
    double min = DBL_MAX;
    int bestToChoose = startNode;
	double distancesToDeliveryPoints=0;
	double max=0;
	double profitValue;

	for(int i=0; i< unvisitedPickups.size(); i++)
	{
		
        for(int j=0; j< unvisitedDeliveries.size(); j++)
		{
			//if (inputDataProcessor.getDistance(closestPoint(unvisitedDeliveries[j],tour), unvisitedDeliveries[j])< inputDataProcessor.getDistance(closestPoint(unvisitedPickups[i],tour), unvisitedPickups[i]))
			//{
			profitValue=profitPerDistanceMatrix[unvisitedPickups[i]][unvisitedDeliveries[j]] / (inputDataProcessor.getDistance(closestPoint(unvisitedPickups[i],tour), unvisitedPickups[i]));
			//}
			if (profitValue>max)
			{
				//max=profitPerDistanceMatrix[unvisitedPickups[i]][unvisitedDeliveries[j]]/(inputDataProcessor.getDistance(closestPoint(unvisitedPickups[i],tour), unvisitedPickups[i]));
				max=profitValue;
				bestToChoose=unvisitedPickups[i];
			}
		}
	}
	return bestToChoose;

	/* old code
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
	old code end */

}

int TourGreedyEinzelnePunkteSeriell::getNextPickupPointRandomised(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode, vector <int> tour)
{
    double min = DBL_MAX;
    int bestToChoose = 0;
	double profitsWithOtherDeliveryPoints=0;
	double profitValue=0;
	vector <double> cumulatedProfits;
	cumulatedProfits.push_back(0);
	double randNumber;
	double max;
    for(int i=0; i< unvisitedPickups.size(); i++)
    {
		max=0;
        for(int j=0; j< unvisitedDeliveries.size(); j++)
		{
			//if (inputDataProcessor.getDistance(closestPoint(unvisitedDeliveries[j],tour), unvisitedDeliveries[j])< inputDataProcessor.getDistance(closestPoint(unvisitedPickups[i],tour), unvisitedPickups[i]))
			//{
			profitValue=profitPerDistanceMatrix[unvisitedPickups[i]][unvisitedDeliveries[j]] / (inputDataProcessor.getDistance(closestPoint(unvisitedPickups[i],tour), unvisitedPickups[i]));
			//}
			if (profitValue>max)
			{
				//max=profitPerDistanceMatrix[unvisitedPickups[i]][unvisitedDeliveries[j]]/(inputDataProcessor.getDistance(closestPoint(unvisitedPickups[i],tour), unvisitedPickups[i]));
				max=profitValue;
			}
		}
		
			//profitValue=profitPerDistanceMatrix [unvisitedPickups[i]][unvisitedDeliveries[j]] /((inputDataProcessor.getDistance(startNode, unvisitedPickups[i])*10));	
			
		/*
			if (max>0)
			{
			profitsWithOtherDeliveryPoints= profitsWithOtherDeliveryPoints + max;
			}
		*/
		//
		/*
		if (profitsWithOtherDeliveryPoints>0)
		{
		*/
		cumulatedProfits.push_back(cumulatedProfits.back() + max);//+profitsWithOtherDeliveryPoints);	
	/*	
	}
		else
		{
			cumulatedProfits.push_back(cumulatedProfits.back());
		}
		

		profitsWithOtherDeliveryPoints=0;
    */
	}
	
   
	double total=cumulatedProfits.back();
	//srand(time(NULL));
	//int proba=rand();
	randNumber= (double) rand() / (double) (RAND_MAX + 1)* total  ;
	for(int k=1; k<cumulatedProfits.size(); k++)
	{
		if (randNumber<=cumulatedProfits[k] && randNumber>cumulatedProfits[k-1])
		{
			bestToChoose = unvisitedPickups[k-1];
			break;
		}
	}

	return bestToChoose;
}


int TourGreedyEinzelnePunkteSeriell::getNextPickupPointRandomisedBest15(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode, vector <int> tour, int selectionPop)
{
    double min = DBL_MAX;
    int bestToChoose = 0;
	double profitsWithOtherDeliveryPoints=0;
	double profitValue=0;
	mypair best15;
	vector <mypair> best15Pairs;
	double randNumber;
	double max;
    for(int i=0; i< unvisitedPickups.size(); i++)
    {
		max=0;
        for(int j=0; j< unvisitedDeliveries.size(); j++)
		{
			//if (inputDataProcessor.getDistance(closestPoint(unvisitedDeliveries[j],tour), unvisitedDeliveries[j])< inputDataProcessor.getDistance(closestPoint(unvisitedPickups[i],tour), unvisitedPickups[i]))
			//{
			profitValue=profitPerDistanceMatrix[unvisitedPickups[i]][unvisitedDeliveries[j]] / (inputDataProcessor.getDistance(closestPoint(unvisitedPickups[i],tour), unvisitedPickups[i]));
			//}
			if (profitValue>max)
			{
				//max=profitPerDistanceMatrix[unvisitedPickups[i]][unvisitedDeliveries[j]]/(inputDataProcessor.getDistance(closestPoint(unvisitedPickups[i],tour), unvisitedPickups[i]));
				max=profitValue;
			}

		}
		
		best15Pairs.push_back(std::make_pair( max, i));
	}
	std::sort(best15Pairs.begin(), best15Pairs.end(), comparator);
	for(int k=1; k<best15Pairs.size(); k++)
	{
		best15Pairs[k].first=best15Pairs[k-1].first+best15Pairs[k].first;
	}

	if (best15Pairs.size()>selectionPop)
	{
	best15Pairs.erase(best15Pairs.begin()+selectionPop, best15Pairs.end());
	}
	best15Pairs.insert(best15Pairs.begin(),std::make_pair(0,0));
   
	double total=best15Pairs.back().first;
	randNumber= (double) rand() / (double) (RAND_MAX + 1)* total  ;
	for(int k=1; k<best15Pairs.size(); k++)
	{
		if (randNumber<=best15Pairs[k].first && randNumber>best15Pairs[k-1].first)
		{
			bestToChoose = unvisitedPickups[best15Pairs[k-1].second];
			break;
		}
	}

	return bestToChoose;
}





int TourGreedyEinzelnePunkteSeriell::getNextDeliveryPoint(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode, vector <int> tour)
{
    double min = DBL_MAX;
    int bestToChoose = startNode;
	double distancesToPickupPoints=0;
	double max=0;
	double profitValue=0;

	for(int i=0; i< unvisitedDeliveries.size(); i++)
    {
		for(int j=0; j< unvisitedPickups.size(); j++)
		{
			//if (inputDataProcessor.getDistance(closestPoint(unvisitedPickups[j],tour), unvisitedPickups[j])< inputDataProcessor.getDistance(closestPoint(unvisitedDeliveries[i],tour), unvisitedDeliveries[i]))
			//{
			profitValue=profitPerDistanceMatrix[unvisitedPickups[j]][unvisitedDeliveries[i]] / (inputDataProcessor.getDistance(closestPoint(unvisitedDeliveries[i],tour), unvisitedDeliveries[i]));
			//}
			if (profitValue>max)
			{
				//max=profitPerDistanceMatrix[unvisitedPickups[i]][unvisitedDeliveries[j]]/(inputDataProcessor.getDistance(closestPoint(unvisitedPickups[i],tour), unvisitedPickups[i]));
				max=profitValue;
				bestToChoose=unvisitedDeliveries[i];
			}
		}
		
	}
	return bestToChoose;



	/*old code
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
	*/
}

int TourGreedyEinzelnePunkteSeriell::getNextDeliveryPointRandomised(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode,vector <int> tour)
{
   
    int bestToChoose = 0;
	double profitsToGet=0;
	double profitValue=0;
	double quantityToTransfer;
	vector <double> cumulatedDistances;
	double travelDistance;
	cumulatedDistances.push_back(0);
	double randNumber;
	double max;
    for(int i=0; i< unvisitedDeliveries.size(); i++)
    {
		max=0;
		 for(int j=0; j< unvisitedPickups.size(); j++)
		{
			//if (inputDataProcessor.getDistance(closestPoint(unvisitedPickups[j],tour), unvisitedPickups[j])< inputDataProcessor.getDistance(closestPoint(unvisitedDeliveries[i],tour), unvisitedDeliveries[i]))
			//{
			profitValue=profitPerDistanceMatrix[unvisitedPickups[j]][unvisitedDeliveries[i]] / (inputDataProcessor.getDistance(closestPoint(unvisitedDeliveries[i],tour), unvisitedDeliveries[i]));
			//}
			if (profitValue>max)
			{
				//max=profitPerDistanceMatrix[unvisitedPickups[i]][unvisitedDeliveries[j]]/(inputDataProcessor.getDistance(closestPoint(unvisitedPickups[i],tour), unvisitedPickups[i]));
				max=profitValue;
			}
		}
		
        /*for(int j=0; j< unvisitedPickups.size(); j++)
		{*/
		//quantityToTransfer= min(-(inputDataProcessor.basicData[unvisitedDeliveries[i]].quantity), inputDataProcessor.basicData[startNode].quantity);
		/*
		 profitsToGet= -(inputDataProcessor.basicData[unvisitedDeliveries[i]].profit); // * quantityToTransfer;
		travelDistance=inputDataProcessor.getDistance(closestPoint(unvisitedDeliveries[i],tour), unvisitedDeliveries[i]);
		profitValue=profitsToGet/travelDistance;
		*/
		cumulatedDistances.push_back(cumulatedDistances.back()+max);
		//profitsToGet=0;
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

int TourGreedyEinzelnePunkteSeriell::getNextDeliveryPointRandomisedBest15(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode, vector <int> tour, int selectionPop)
{
    double min = DBL_MAX;
    int bestToChoose = 0;
	double profitsWithOtherDeliveryPoints=0;
	double profitValue=0;
	mypair best15;
	vector <mypair> best15Pairs;
	double randNumber;
	double max;

	  for(int i=0; i< unvisitedDeliveries.size(); i++)
    {
		max=0;
		 for(int j=0; j< unvisitedPickups.size(); j++)
		{
			//if (inputDataProcessor.getDistance(closestPoint(unvisitedPickups[j],tour), unvisitedPickups[j])< inputDataProcessor.getDistance(closestPoint(unvisitedDeliveries[i],tour), unvisitedDeliveries[i]))
			//{
			profitValue=profitPerDistanceMatrix[unvisitedPickups[j]][unvisitedDeliveries[i]] / (inputDataProcessor.getDistance(closestPoint(unvisitedDeliveries[i],tour), unvisitedDeliveries[i]));
			//}
			if (profitValue>max)
			{
				//max=profitPerDistanceMatrix[unvisitedPickups[i]][unvisitedDeliveries[j]]/(inputDataProcessor.getDistance(closestPoint(unvisitedPickups[i],tour), unvisitedPickups[i]));
				max=profitValue;
			}
		}

		best15Pairs.push_back(std::make_pair( max, i));
	}
	std::sort(best15Pairs.begin(), best15Pairs.end(), comparator);
	for(int k=1; k<best15Pairs.size(); k++)
	{
		best15Pairs[k].first=best15Pairs[k-1].first+best15Pairs[k].first;
	}
	if (best15Pairs.size()>selectionPop)
	{
	best15Pairs.erase(best15Pairs.begin()+selectionPop, best15Pairs.end());
	}

	best15Pairs.insert(best15Pairs.begin(),std::make_pair(0,0));
	
   
	double total=best15Pairs.back().first;
	randNumber= (double) rand() / (double) (RAND_MAX + 1)* total  ;
	for(int k=1; k<best15Pairs.size(); k++)
	{
		if (randNumber<=best15Pairs[k].first && randNumber>best15Pairs[k-1].first)
		{
			bestToChoose = unvisitedDeliveries[best15Pairs[k-1].second];
			break;
		}
	}

	return bestToChoose;
}





void TourGreedyEinzelnePunkteSeriell::putPointInBestPosition(int whichTour, int pointToInsert)
{
	double minTourLengthExtension=DBL_MAX;
	double TourLengthExtension;
	int posToInsert;

	if (solutionTours[whichTour].size()!=1)
	{
		for (int i = 0; i < solutionTours[whichTour].size()-1; i++)
		{
			if (isTotalLengthUnderLimit( solutionTours[whichTour], pointToInsert, i+1) /*&& (inputDataProcessor.getQuantity(pointToInsert)) < bufferPlus[whichTour][i+1]*/ )
			{
				TourLengthExtension= inputDataProcessor.getDistance (solutionTours[whichTour][i], pointToInsert)+ inputDataProcessor.getDistance (pointToInsert, solutionTours[whichTour][i+1])- inputDataProcessor.getDistance(solutionTours[whichTour][i], solutionTours[whichTour][i+1]);
				if (TourLengthExtension < minTourLengthExtension) 
				{
					minTourLengthExtension=TourLengthExtension;
					posToInsert= i+1;
				}
			}
		}
	}
	else
	{
	solutionTours[whichTour].insert(solutionTours[whichTour].begin()+1, pointToInsert);
	pickupInsertedPosition=1;
	}
	
	if (minTourLengthExtension!=DBL_MAX)
	{
		solutionTours[whichTour].insert(solutionTours[whichTour].begin()+posToInsert, pointToInsert);
		pickupPointInserted=true;
		lastPickupInserted=pointToInsert;
		pickupInsertedPosition=posToInsert;
		//intensity[whichTour].insert(intensity[whichTour].begin()+posToInsert, 1);
		//LoadCalculator newLoad( load[whichTour], goodsOnTheLorry[whichTour], bufferPlus[whichTour], bufferMinus[whichTour], intensity[whichTour]);
		unvisitedNodes[pointToInsert]=0;
	}
	else
	{
		pickupInsertedPosition=false;
	}
						
	//unvisitedNodesForOneTour[pointToInsert]=0;
	
}

void TourGreedyEinzelnePunkteSeriell::putDeliveryPointInBestPosition(int whichTour, int pointToInsert, int pickupInsertedPosition)
{
	deliveryPointInserted=false;
	vector <int> wrongDeliveries;
	double minTourLengthExtension=DBL_MAX;
	double TourLengthExtension;
	int posToInsert;

	if (solutionTours[whichTour].size()!=1)
	{
		for (int i = pickupInsertedPosition; i < solutionTours[whichTour].size()-1; i++)
		{
			if (isTotalLengthUnderLimit( solutionTours[whichTour], pointToInsert, i+1) /*&& (inputDataProcessor.getQuantity(pointToInsert)) < bufferPlus[whichTour][i+1]*/ )
			{
				TourLengthExtension= inputDataProcessor.getDistance (solutionTours[whichTour][i], pointToInsert)+ inputDataProcessor.getDistance (pointToInsert, solutionTours[whichTour][i+1])- inputDataProcessor.getDistance (solutionTours[whichTour][i], solutionTours[whichTour][i+1]);
				if (TourLengthExtension < minTourLengthExtension) 
				{
					minTourLengthExtension=TourLengthExtension;
					posToInsert= i+1;
				}
			}
		}
	}

	else
	{
	solutionTours[whichTour].insert(solutionTours[whichTour].begin()+1, pointToInsert);
	}
	
	
	if (minTourLengthExtension!=DBL_MAX)
	{
		solutionTours[whichTour].insert(solutionTours[whichTour].begin()+posToInsert, pointToInsert);
		//intensity[whichTour].insert(intensity[whichTour].begin()+posToInsert, 1);
		//LoadCalculator newLoad( load[whichTour], goodsOnTheLorry[whichTour], bufferPlus[whichTour], bufferMinus[whichTour], intensity[whichTour]);
		unvisitedNodes[pointToInsert]=0;
		deliveryPointInserted=true;

		for (int i=0; i< wrongDeliveries.size(); i++)
		{
			unvisitedNodes[wrongDeliveries[i]]=1;
			deliveryPoints.push_back(wrongDeliveries[i]);
		}
		wrongDeliveries.erase(wrongDeliveries.begin(), wrongDeliveries.end());
	}
	else
	{
		unvisitedNodesForOneTour[pointToInsert]=0;
		wrongDeliveries.push_back(pointToInsert);
		deliveryPointInserted=false;
	}


	//unvisitedNodesForOneTour[pointToInsert]=0;a
	
}

int TourGreedyEinzelnePunkteSeriell::closestPoint( int node, vector <int> tour)
{
	double minLength=DBL_MAX;
	double min=tour.back();
	for (int i=0; i<tour.size();i++)
	{
		if (inputDataProcessor.getDistance(tour[i], node)<minLength)
		{
			minLength=inputDataProcessor.getDistance(tour[i], node);
			min=i;
		}
	}
return min;
}