#include "TourPickupDeliveryPointPairsParallel.h"
#include <iostream>
#include <algorithm> 
#include <math.h>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */



using namespace std;

bool comparator2 (const mypair2& l, const mypair2& r )
{ return l.first> r.first;}

TourPickupDeliveryPointPairsParallel::TourPickupDeliveryPointPairsParallel(string inputFile, int selectionPop): OrienteeringProblemWithPickupsAndDeliveries(inputFile)
{
	
	for (int seedNumber=0; seedNumber<100; seedNumber++) // 100 seed run
	{
			
		solutionTours.clear();
		srand(seedNumber);
		double timeStart=clock();
		bestPairs.assign(2,0);
		unvisitedNodes.assign(problemSize,1);
		unvisitedNodes[0]=0;
		unvisitedNodes[1]=0;
		vector <mypair2> initMatrix;
		probabilities.assign(3, initMatrix);
		calcTourChoosePickupAndDeliveryPointPairs3(selectionPop);
	
		profitsOfAllTheTours( seedNumber,timeStart);
		/*
		for (int i=0; i<solutionTours.size();i++)
		{
			cout << "The tour length of the " << i+1 << ". tour is: " << getTourLength(solutionTours[i]) << endl;
		} 
		*/
		runTwoOpt(seedNumber, timeStart);
		//profitsOfAllTheTours( seedNumber,timeStart);

		if (selectionPop==1)
		{
		break;
		}
	}
	runExcelExport(inputFile, "heurParallelPairs" + std::to_string(selectionPop));


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





void TourPickupDeliveryPointPairsParallel::calcTourChoosePickupAndDeliveryPointPairs3(int selectionPop) {
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

		if (selectionPop==1)
			{
			getPickUpDeliveryPointPairsTwoPointsAddedRandomised(unvisitedNodes, tour.back(),bestPairs, i%3, selectionPop);
			}
			else
			{
				if (selectionPop==3 || selectionPop==15)
				{
				
				getPickUpDeliveryPointPairsTwoPointsAddedRandomised(unvisitedNodes, tour.back(),bestPairs, i%3, selectionPop);
		
				}
				else
				{
				getPickUpDeliveryPointPairsTwoPointsAddedRandomised(unvisitedNodesForOneTour, tour.back(), bestPairs, i%3,probabilities[i%3].size());
				}
			}
		




		//getPickUpDeliveryPointPairs(unvisitedNodes, bestPairs, i%3);
		//getPickUpDeliveryPointPairsTwoPointsAddedRandomised(unvisitedNodes, tour.back(),bestPairs, i%3, 3);
		

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
				//wrongPairs[i%numberOfTours].push_back(bestPairs);
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

void TourPickupDeliveryPointPairsParallel::getPickUpDeliveryPointPairsTwoPointsAddedRandomised (vector<int> unvisitedCities, int startNode, vector <int> & bestPairs, int whichTour, int selectionPop)
{
    
	double max = -DBL_MAX;
	bool wrongPairsFound=false;
	//vector <double> probabilities;
	vector <vector <int>> positions;
	int randNumber;
	std::vector <int> a;
	a.push_back(0);
	a.push_back(0);
	//probabilities.push_back(0);
    //int best = startNode;
	vector <double> probabilitiesCumulated;
	probabilitiesCumulated.assign(selectionPop+1,0);
    
	if (probabilities[whichTour].size()==0)
	{
		//for(int n=0; n<numberOfTours; n++)
		//{
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
				
					/*
					for(int k=0; k<wrongPairs[whichTour].size(); k++)
					{
						if ( wrongPairs[whichTour][k][0]!=0 && wrongPairs[whichTour][k][1]!=0)
						{
							wrongPairsFound=true;
						}
					}
					*/
			
					if (profitPerDistanceMatrix[i][j]>0)// && (wrongPairsFound==false))
					{
					std::vector <int> pair;
					pair.push_back(i);
					pair.push_back(j);
					probabilities[whichTour].push_back(std::make_pair(profitPerDistanceMatrix[i][j], pair));
		
					//positions.push_back(pair);
					}	
				}
			}
			std::sort (probabilities[whichTour].begin(), probabilities[whichTour].end(), comparator2);
		
			probabilities[whichTour].insert(probabilities[whichTour].begin(),std::make_pair(0,a));
		
	}

	for(int k=1; k<probabilities[whichTour].size(); k++)
	{
		if ((unvisitedCities[probabilities[whichTour][k].second[0]]==0) || (unvisitedCities[probabilities[whichTour][k].second[1]]==0))
		{
			probabilities[whichTour].erase(probabilities[whichTour].begin()+k);
		}	
	}


	if (selectionPop+1<probabilities[whichTour].size())
	{
		for(int k=1; k<selectionPop+1 ; k++)
		{
			probabilitiesCumulated[k]=probabilitiesCumulated[k-1]+probabilities[whichTour][k].first;
		}
		int total=floor(probabilitiesCumulated[selectionPop]);
		
	
	}
	else 
	{
		for(int k=1; k<probabilities[whichTour].size() ; k++)
		{
			probabilitiesCumulated[k]=probabilitiesCumulated[k-1]+probabilities[whichTour][k].first;
		}
		int total=floor(probabilitiesCumulated[probabilities.size()-1]);
		for(int m=probabilities.size(); m<selectionPop ; m++)
		{
			probabilitiesCumulated.erase(probabilitiesCumulated.begin()+probabilities.size(), probabilitiesCumulated.end());
		}
	}


		if (probabilitiesCumulated.back()>0)
		{
		
		int total=floor(probabilitiesCumulated.back());
		//srand(time(NULL));
		randNumber= (double) rand() / (RAND_MAX + 1) * total;
			for(int k=1; k<probabilitiesCumulated.size(); k++)
			{
				if (randNumber<=probabilitiesCumulated[k] && randNumber>probabilitiesCumulated[k-1])
				{
					bestPairs[0]=probabilities[whichTour][k].second[0];
					bestPairs[1]=probabilities[whichTour][k].second[1];
					probabilities[whichTour].erase(probabilities[whichTour].begin()+k);
					break;
				}
			}
		}
	//cout << "The best pair is: " << bestPairs[0] << " and " << bestPairs[1] << endl;
   return;
}




