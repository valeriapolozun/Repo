#include "TourPickupDeliveryPointPairs.h"
#include "LoadCalculator.h"
#include <iostream>
#include <math.h>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <algorithm> 

using namespace std;

typedef std::pair <double, vector <int>> mypair2;
bool comparator (const mypair2& l, const mypair2& r )
{ return l.first> r.first;}

TourPickupDeliveryPointPairs::TourPickupDeliveryPointPairs(string inputFile, int selectionPop): OrienteeringProblemWithPickupsAndDeliveries(inputFile)
{
	for (int seedNumber=0; seedNumber<5; seedNumber++) // 100 seed run
	{
	//int seedNumber=44;
	solutionTours.clear();


	bestPairs.assign(2,0);
	unvisitedNodes.assign(problemSize,1);
	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;
	srand(seedNumber);
	double timeStart=clock();
	for (int i=0; i<numberOfTours;i++)
	{
		calcTourChoosePickupAndDeliveryPointPairs2(i, selectionPop);
	} 

	profitsOfAllTheTours(seedNumber,timeStart);
	for (int i=0; i<solutionTours.size();i++)
	{
		cout << "The tour length of the " << i+1 << ". tour is: " << getTourLength(solutionTours[i]) << endl;
	}
	}

	runExcelExport(inputFile, "heurSeriellPairs" + std::to_string(selectionPop));


}

TourPickupDeliveryPointPairs::~TourPickupDeliveryPointPairs()
{
}

void TourPickupDeliveryPointPairs::pickUpPointToChoose(vector<int> & nodes)
{
	for (int i=0; i<nodes.size();i++)
	{
		if (inputDataProcessor.getBasicData()[i].quantity<0)
		nodes[i]=0;
	}	
}




void TourPickupDeliveryPointPairs::calcTourChoosePickupAndDeliveryPointPairs2(int whichTour, int selectionPop){
	
	
	vector<int> unvisitedNodesForOneTour;
	probabilities.clear();

	unvisitedNodesForOneTour= unvisitedNodes;
	
	int startNode=0;
	vector<int> tour;
	tour.push_back(0);
	
	unvisitedNodesForOneTour[0]=0;
	unvisitedNodesForOneTour[1]=0;
	bool pickupInserted=false;
	tour.push_back(1);
	
	getProfitMatrixForPickupAndDeliveryPairs (startNode, whichTour);
	
	bool criterion=true;
	//for (int i = 0; i < (problemSize-2); i++)  /// TO DO : Calculate how many times it should run (problemSize-2) is wrong
	while (criterion)
	{
			if (selectionPop==1)
			{
			getPickUpDeliveryPointPairsTwoPointsAdded(unvisitedNodesForOneTour, tour.back(), bestPairs, whichTour);
			}
			else
			{
				if (selectionPop==3 || selectionPop==15)
				{
				getPickUpDeliveryPointPairsTwoPointsAddedRandomisedBest15(unvisitedNodesForOneTour, tour.back(), bestPairs, whichTour, selectionPop);
				}
				else
				{
				getPickUpDeliveryPointPairsTwoPointsAddedRandomised(unvisitedNodesForOneTour, tour.back(), bestPairs, whichTour);
				}
			}
		
		
		
		
		if (bestPairs[0]==0)
		{
			criterion=false;
		}
		
		if (!bestPairs[0]==0)
		{
		
			double minTourLengthExtension=DBL_MAX;
			double TourLengthExtension;
			int posToInsert;
			int posToInsert2;
			for (int j = 0; j < tour.size()-1; j++)
			{
				for (int k = j; k < tour.size()-1; k++)
				{

					if (isTotalLengthUnderLimit2NodesDifferentPlace( tour, bestPairs[0], j+1, bestPairs[1], k+2))
					{
						int node2pos;
						if (k==j)
						{
						node2pos=bestPairs[0];
						}
						else
						{
						node2pos=tour[k];
						}
						TourLengthExtension= inputDataProcessor.getDistance (tour[j], bestPairs[0])+ inputDataProcessor.getDistance (bestPairs[0], tour[j+1])-inputDataProcessor.getDistance (tour[j], tour[j+1])+inputDataProcessor.getDistance (node2pos, bestPairs[1])+ inputDataProcessor.getDistance (bestPairs[1], tour[k+1])-+ inputDataProcessor.getDistance (node2pos, tour[k+1]);;
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
				tour.insert(tour.begin()+posToInsert, bestPairs[0]);
				tour.insert(tour.begin()+posToInsert2, bestPairs[1]);
				unvisitedNodesForOneTour[bestPairs[0]]=0;
				unvisitedNodes[bestPairs[0]]=0;
				unvisitedNodesForOneTour[bestPairs[1]]=0;
				unvisitedNodes[bestPairs[1]]=0;
			}
			else
			{
				wrongPairs[whichTour].push_back(bestPairs);
				profitPerDistanceMatrix[bestPairs[0]] [bestPairs[1]]=-DBL_MAX;
					
			}

			bestPairs[0]=0;
			bestPairs[1]=0;
		}
		
	}
	
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




int TourPickupDeliveryPointPairs::getPickUpDeliveryPointPairsOnePointAdded (vector<int> unvisitedCities, int startNode)
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


void TourPickupDeliveryPointPairs::getPickUpDeliveryPointPairsTwoPointsAdded (vector<int> unvisitedCities, int startNode, vector <int> & bestPairs, int whichTour)
{
    double max = -DBL_MAX;
    //int best = startNode;
    
	for(int i=0; i<unvisitedCities.size(); i++)
    {
		if(unvisitedCities[i]==0) continue; 
		for(int j=0; j<unvisitedCities.size(); j++)
		{
			if(unvisitedCities[j]==0) continue; 
			//if (startNode == i) continue; 
			//if (startNode == j) continue; 
			//getProfitMatrixForPickupAndDeliveryPairs (startNode, whichTour);
			//cout << "the i: " << i << "the profitperdistancematrix[startnode][i]" << profitPerDistanceMatrix[startNode][i] << " a max value: " << max << endl;
			if(profitPerDistanceMatrix[i][j]>max)
			{			
            max = profitPerDistanceMatrix[i][j];
			bestPairs[0]=i;
			bestPairs[1]=j;
			}


		 }
    }

	cout << "The best pair is: " << bestPairs[0] << " and " << bestPairs[1] << endl;
   return;
}


void TourPickupDeliveryPointPairs::getPickUpDeliveryPointPairsTwoPointsAddedRandomised (vector<int> unvisitedCities, int startNode, vector <int> & bestPairs, int whichTour)
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
		for(int j=1; j<unvisitedCities.size(); j++)
		{
			if(unvisitedCities[j]==0) continue; 
			//if (startNode == i) continue; 
			//if (startNode == j) continue; 
			//getProfitMatrixForPickupAndDeliveryPairs (startNode, whichTour);
			//cout << "the i: " << i << "the profitperdistancematrix[startnode][i]" << profitPerDistanceMatrix[startNode][i] << " a max value: " << max << endl;
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
	cout << "The best pair is: " << bestPairs[0] << " and " << bestPairs[1] << endl;
   return;
}


void TourPickupDeliveryPointPairs::getPickUpDeliveryPointPairsTwoPointsAddedRandomisedBest15 (vector<int> unvisitedCities, int startNode, vector <int> & bestPairs, int whichTour, int selectionPop)
{
    double max = -DBL_MAX;
	//vector <mypair2> probabilities;
	vector <double> probabilitiesCumulated;
	probabilitiesCumulated.assign(selectionPop+1,0);
	vector <vector <int>> positions;
	std::vector <int> a;
	a.push_back(0);
	a.push_back(0);
	int randNumber;
	//probabilities.push_back(0);
    //int best = startNode;
    

	if( probabilities.size()==0)
	{

		for(int i=0; i<unvisitedCities.size(); i++)
		{
			if(unvisitedCities[i]==0) continue; 
			for(int j=0; j<unvisitedCities.size(); j++)
			{
				if(unvisitedCities[j]==0) continue; 
				//if (startNode == i) continue; 
				//if (startNode == j) continue; 
				//getProfitMatrixForPickupAndDeliveryPairs (startNode, whichTour);
				//cout << "the i: " << i << "the profitperdistancematrix[startnode][i]" << profitPerDistanceMatrix[startNode][i] << " a max value: " << max << endl;
				if (profitPerDistanceMatrix[i][j]>0)
				{
				std::vector <int> pair;
				pair.push_back(i);
				pair.push_back(j);
				positions.push_back(pair);
				probabilities.push_back(std::make_pair(profitPerDistanceMatrix[i][j], pair));
				}	
			}
		}
	//best15Pairs.push_back(std::make_pair( max, i));
	std::sort (probabilities.begin(), probabilities.end(), comparator);

	/*
	if (probabilities.size()>selectionPop)
	{
	probabilities.erase(probabilities.begin()+selectionPop, probabilities.end());
	}
	*/

	probabilities.insert(probabilities.begin(),std::make_pair(0,a));

		/*
		for(int k=1; k<probabilities.size(); k++)
		{
			probabilities[k].first=probabilities[k-1].first+probabilities[k].first;
		}
		*/

	}

	for(int k=1; k<probabilities.size(); k++)
	{
		if ((unvisitedCities[probabilities[k].second[0]]==0) || (unvisitedCities[probabilities[k].second[1]]==0))
		{
			probabilities.erase(probabilities.begin()+k);
		}	
	}


	if (selectionPop+1<probabilities.size())
	{
		for(int k=1; k<selectionPop+1 ; k++)
		{
			probabilitiesCumulated[k]=probabilitiesCumulated[k-1]+probabilities[k].first;
		}
		int total=floor(probabilitiesCumulated[selectionPop]);
		
	
	}
	else 
	{
		for(int k=1; k<probabilities.size() ; k++)
		{
			probabilitiesCumulated[k]=probabilitiesCumulated[k-1]+probabilities[k].first;
		}
		int total=floor(probabilitiesCumulated[probabilities.size()-1]);
		for(int m=probabilities.size(); m<selectionPop ; m++)
		{
			probabilitiesCumulated.erase(probabilitiesCumulated.begin()+probabilities.size(), probabilitiesCumulated.end());
		}
	}

	//probabilitiesCumulated.insert(probabilitiesCumulated.begin(),0);

		
	if (probabilitiesCumulated.back()>0)
		{
		
		int total=floor(probabilitiesCumulated.back());
		//srand(time(NULL));
		randNumber= (double) rand() / (RAND_MAX + 1) * total;
			for(int k=1; k<probabilitiesCumulated.size(); k++)
			{
				if (randNumber<=probabilitiesCumulated[k] && randNumber>probabilitiesCumulated[k-1])
				{
					bestPairs[0]=probabilities[k].second[0];
					bestPairs[1]=probabilities[k].second[1];
					probabilities.erase(probabilities.begin()+k);
					break;
				}
			}
		}
	cout << "The best pair is: " << bestPairs[0] << " and " << bestPairs[1] << endl;
   return;
}

