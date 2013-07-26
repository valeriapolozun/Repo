#include "TourPickupDeliveryPointPairsParallel.h"
#include <iostream>
#include <algorithm> 
#include <math.h>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */



using namespace std;


TourPickupDeliveryPointPairsParallel::TourPickupDeliveryPointPairsParallel(string inputFile): OrienteeringProblemWithPickupsAndDeliveries(inputFile)
{
	bestPairs.assign(2,0);
	unvisitedNodes.assign(problemSize,1);
	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;
	calcTourChoosePickupAndDeliveryPointPairs3();
	
	profitsOfAllTheTours(1,0);
	/*
	for (int i=0; i<solutionTours.size();i++)
	{
		cout << "The tour length of the " << i+1 << ". tour is: " << getTourLength(solutionTours[i]) << endl;
	} 
	*/


}

TourPickupDeliveryPointPairsParallel::~TourPickupDeliveryPointPairsParallel()
{
}

void TourPickupDeliveryPointPairsParallel::pickUpPointToChoose(vector<int> & nodes)
{
	for (int i=0; i<nodes.size();i++)
	{
		if (inputDataProcessor.getBasicData()[i].quantity<0)
		nodes[i]=0;
	}	
}





void TourPickupDeliveryPointPairsParallel::calcTourChoosePickupAndDeliveryPointPairs3() {
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
		getProfitMatrixForPickupAndDeliveryPairsParallel(numberOfTours);
	}

	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;
	/*
	unvisitedNodesForOneTour[0]=0;
	unvisitedNodesForOneTour[1]=0;
	*/
	
	for (int i = 0; i < (problemSize-2); i++)  /// TO DO : Calculate how many times it should run (problemSize-2) is wrong
	{
		//getPickUpDeliveryPointPairs(unvisitedNodes, bestPairs, i%3);
		getPickUpDeliveryPointPairsTwoPointsAddedRandomised(unvisitedNodes, tour.back(),bestPairs, i%3);
		

		if (!bestPairs[0]==0)
		{
			double minTourLengthExtension=DBL_MAX;
			double TourLengthExtension;
			int posToInsert;
			int posToInsert2;
			for (int j = 0; j < solutionTours[i%numberOfTours].size()-1; j++)
			{
				for (int k = j; k < solutionTours[i%numberOfTours].size()-1; k++)
				{

					if (isTotalLengthUnderLimit2NodesDifferentPlace( solutionTours[i%numberOfTours], bestPairs[0], j+1, bestPairs[1], k+2))
					{
						int node2pos;
						if (k==j)
						{
							node2pos=bestPairs[0];
						}
						else
						{
							node2pos=solutionTours[i%numberOfTours][k];
						}
						TourLengthExtension= inputDataProcessor.getDistance (solutionTours[i%numberOfTours][j], bestPairs[0])+ inputDataProcessor.getDistance (bestPairs[0], solutionTours[i%numberOfTours][j+1])+inputDataProcessor.getDistance (node2pos, bestPairs[1])+ inputDataProcessor.getDistance (bestPairs[1], solutionTours[i%numberOfTours][k+1]);
						if (TourLengthExtension < minTourLengthExtension) 
						{
							minTourLengthExtension=TourLengthExtension;
							posToInsert= j+1;
							posToInsert2=k+2;
						}

					}
				}
			}

			if (minTourLengthExtension!=DBL_MAX)
			{
				solutionTours[i%numberOfTours].insert(solutionTours[i%numberOfTours].begin()+posToInsert, bestPairs[0]);
				solutionTours[i%numberOfTours].insert(solutionTours[i%numberOfTours].begin()+posToInsert2, bestPairs[1]);
				unvisitedNodesForOneTour[bestPairs[0]]=0;
				unvisitedNodes[bestPairs[0]]=0;
				unvisitedNodesForOneTour[bestPairs[1]]=0;
				unvisitedNodes[bestPairs[1]]=0;
			}
			else
			{
				wrongPairs[i%numberOfTours].push_back(bestPairs);
				profitPerDistanceMatrix[bestPairs[0]] [bestPairs[1]]=-DBL_MAX;
			
					
			}

			bestPairs[0]=0;
			bestPairs[1]=0;
		}
		
	}

	/* old code
			if (isTotalLengthUnderLimit2Nodes(solutionTours[i%numberOfTours], bestPairs))
			{
				for (int j = 0; j < bestPairs.size(); j++)
				{
					solutionTours[i%numberOfTours].push_back(bestPairs[j]);
					//unvisitedNodesForOneTour[bestPairs[i]]=0;
					unvisitedNodes[bestPairs[j]]=0;
				}
			}
			else
			{
				wrongPairs[i%numberOfTours].push_back(bestPairs);
			}
			
			bestPairs[0]=0;
			bestPairs[1]=0;
			
		}
	}
	*/
	/*
	for (int i = 0; i < numberOfTours; i++)
	{
		solutionTours[i].push_back(1);
	}
	*/

	for (int i=0; i < tour.size(); i++)
	{
		unvisitedNodes[tour[i]]=0;
	}
	
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

	for (int i=0; i < solutionTours.size(); i++)
	{
		for (int j=0; j < solutionTours[i].size(); j++)
		{
			cout << j+1 << ". place in the " << i+1 << " . tour is: " << solutionTours[i][j] << endl;
		}
	}


}





void TourPickupDeliveryPointPairsParallel::getPickUpDeliveryPointPairs (std::vector<int> unvisitedCities, std::vector <int> & bestPair, int whichTour) 
{
    double max = -DBL_MAX;
    //int best = startNode;
	bool wrongPairsFound=false;
    
	for(int i=0; i<unvisitedCities.size(); i++)
    {
		if(unvisitedCities[i]==0) continue; 
		for(int j=0; j<unvisitedCities.size(); j++)
		{
			if(unvisitedCities[j]==0) continue; 

			//getProfitMatrixForPickupAndDeliveryPairsParallel(whichTour);
			//cout << "the i: " << i << "the profitperdistancematrix[startnode][i]" << profitPerDistanceMatrix[startNode][i] << " a max value: " << max << endl;
			for(int k=0; k<wrongPairs[whichTour].size(); k++)
			{
				if ( wrongPairs[whichTour][k][0]!=0 && wrongPairs[whichTour][k][1]!=0)
				{
					wrongPairsFound=true;
				}
			}
			
			
			//bool wrongPairsfound=(std::find(wrongPairs[whichTour].begin(), wrongPairs[whichTour].end(), (i,j)) != wrongPairs[whichTour].end());
			if ((profitPerDistanceMatrix[i][j]>max) && (wrongPairsFound==false))
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

void TourPickupDeliveryPointPairsParallel::getPickUpDeliveryPointPairsTwoPointsAddedRandomised (vector<int> unvisitedCities, int startNode, vector <int> & bestPairs, int whichTour)
{
    double max = -DBL_MAX;
	bool wrongPairsFound=false;
	vector <double> probabilities;
	vector <vector <int>> positions;
	int randNumber;
	probabilities.push_back(0);
    //int best = startNode;
    
	for(int i=0; i<unvisitedCities.size(); i++)
    {
		if(unvisitedCities[i]==0) continue; 
		for(int j=1; j<unvisitedCities.size(); j++)
		{
			if(unvisitedCities[j]==0) continue; 
			//if (startNode == i) continue; 
			//if (startNode == j) continue; 
			//getProfitMatrixForPickupAndDeliveryPairs (startNode, whichTour);
			//cout << "the i: " << i << "the profitperdistancematrix[startnode][i]" << profitPerDistanceMatrix[startNode][i] << " a max value: " << max << endl;
			for(int k=0; k<wrongPairs[whichTour].size(); k++)
			{
				if ( wrongPairs[whichTour][k][0]!=0 && wrongPairs[whichTour][k][1]!=0)
				{
					wrongPairsFound=true;
				}
			}
			
			
			
			
			if (profitPerDistanceMatrix[i][j]>0 && (wrongPairsFound==false))
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
		srand(time(NULL));
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
	//cout << "The best pair is: " << bestPairs[0] << " and " << bestPairs[1] << endl;
   return;
}




