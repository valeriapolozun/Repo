#include "TourPickupDeliveryPointPairsParallelBest.h"
#include <iostream>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace std;


TourPickupDeliveryPointPairsParallelBest::TourPickupDeliveryPointPairsParallelBest(string inputFile): OrienteeringProblemWithPickupsAndDeliveries(inputFile)
{
	
	for (int s=0; s<1; s++)
	{
		solutionTours.clear();
		bestPairs.assign(2,0);
		unvisitedNodes.assign(problemSize,1);
		unvisitedNodes[0]=0;
		unvisitedNodes[1]=0;
		srand(time(NULL));
		calcTourChoosePickupAndDeliveryPointPairs2();
		profitsOfAllTheTours();
		profitsOfAllTheToursOhneGLPK();
		for (int i=0; i<solutionTours.size();i++)
		{
			cout << "The tour length of the " << i+1 << ". tour is: " << getTourLength(solutionTours[i]) << endl;
		}
		
	}
	printSolutions();
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





void TourPickupDeliveryPointPairsParallelBest::calcTourChoosePickupAndDeliveryPointPairs2()
{
	vector<int> unvisitedNodesForOneTour;

	unvisitedNodesForOneTour= unvisitedNodes;
	double additionalDistance;
	double min;
	int minpos;


	int startNode=0;
	vector<int> tour;
	tour.push_back(0);
	tour.push_back(1);


	for (int i = 0; i < numberOfTours; i++)
	{
		solutionTours.push_back(tour);
	}


	unvisitedNodesForOneTour[0]=0;
	unvisitedNodesForOneTour[1]=0;

	getProfitMatrixForPickupAndDeliveryPairsParallelBest();

	for (int i = 0; i < (problemSize-2); i++)  /// TO DO : Calculate how many times it should run (problemSize-2) is wrong
	{

		getPickUpDeliveryPointPairs(unvisitedNodes, bestPairs);

		if (!bestPairs[0]==0)
		{
			
			double minTourLengthExtension=DBL_MAX;
			double TourLengthExtension;
			int posToInsert;
			int posToInsert2;
			int tourToInsert;
		
			for (int m=0; m<numberOfTours;m++)
			{
				for (int j = 0; j < solutionTours[m].size()-1; j++)
				{
					for (int k = j; k < solutionTours[m].size()-1; k++)
					{

						if (isTotalLengthUnderLimit2NodesDifferentPlace( solutionTours[m], bestPairs[0], j+1, bestPairs[1], k+2))
						{
							int node2pos;
							if (k==j)
							{
								node2pos=bestPairs[0];
							}
							else
							{
								node2pos=solutionTours[m][k];
							}
							TourLengthExtension= inputDataProcessor.getDistance (solutionTours[m][j], bestPairs[0])+ inputDataProcessor.getDistance (bestPairs[0], solutionTours[m][j+1])-inputDataProcessor.getDistance(solutionTours[m][j], solutionTours[m] [j+1])+inputDataProcessor.getDistance (node2pos, bestPairs[1])+ inputDataProcessor.getDistance (bestPairs[1], solutionTours[m][k+1])-inputDataProcessor.getDistance(node2pos,solutionTours[m][k+1]);
							if (TourLengthExtension < minTourLengthExtension) 
							{
								minTourLengthExtension=TourLengthExtension;
								posToInsert= j+1;
								posToInsert2=k+2;
								tourToInsert=m;
							}
						}
					}
				}
			}

			if (minTourLengthExtension!=DBL_MAX)
			{
				solutionTours[tourToInsert].insert(solutionTours[tourToInsert].begin()+posToInsert, bestPairs[0]);
				solutionTours[tourToInsert].insert(solutionTours[tourToInsert].begin()+posToInsert2, bestPairs[1]);
				unvisitedNodesForOneTour[bestPairs[0]]=0;
				unvisitedNodes[bestPairs[0]]=0;
				unvisitedNodesForOneTour[bestPairs[1]]=0;
				unvisitedNodes[bestPairs[1]]=0;
			}
			else
			{
				wrongPairs[0].push_back(bestPairs);
				profitPerDistanceMatrix[bestPairs[0]] [bestPairs[1]]=-DBL_MAX;			
			}
			bestPairs[0]=0;
			bestPairs[1]=0;
		}
	}
		
	/*old code	
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
	*/
	
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
	vector <double> probabilities;
	vector <vector <int>> positions;
	int randNumber;
	probabilities.push_back(0);
    //int best = startNode;
    
	for(int i=0; i<unvisitedCities.size(); i++)
    {
		if(unvisitedCities[i]==0) continue; 
		for(int j=0; j<unvisitedCities.size(); j++)
		{
			if(unvisitedCities[j]==0) continue; 
			//getProfitMatrixForPickupAndDeliveryPairsParallelBest();
			//cout << "the i: " << i << "the profitperdistancematrix[startnode][i]" << profitPerDistanceMatrix[startNode][i] << " a max value: " << max << endl;
			
			/* old code - here the absolute best pair will be chosen
			if(profitPerDistanceMatrix[i][j]>max)
			{			
            max = profitPerDistanceMatrix[i][j];
			bestPairs[0]=i;
			bestPairs[1]=j;
			}
			*/

			if (profitPerDistanceMatrix[i][j]>0)
			{
				probabilities.push_back(probabilities.back()+profitPerDistanceMatrix[i][j]);
				std::vector <int> pair;
				pair.push_back(i);
				pair.push_back(j);
				positions.push_back(pair);
			}	
		}
	}
		
	if (probabilities.back()>0)
		
	{
		
		int total=floor(probabilities.back());
	
		randNumber= (double) rand() / (RAND_MAX + 1) * total;
		for(int k=1; k<probabilities.size(); k++)
		{
			if (randNumber<=probabilities[k] && randNumber>probabilities[k-1])
			{
				bestPairs[0]=positions[k-1][0];
				bestPairs[1]=positions[k-1][1];
				break;
			}
		}
	}

	//cout << "The next nearest point is: " << nearest << endl;
   return;
}
