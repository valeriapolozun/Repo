#include "OrienteeringProblemWithPickupsAndDeliveries.h"
#include "ProfitCalculator.h"
#include "ProfitCalculatorOhneGLPK.h"
#include "LoadCalculator.h"
#include "ExcelExporter.h"


#include <iostream>
#include <fstream>
#include <iterator>
#include <glpk.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <time.h>  
#include <windows.h>
#include <algorithm>  

using namespace std;

//static clock_t clock_start;
static clock_t clockStartThisSolution;


OrienteeringProblemWithPickupsAndDeliveries::OrienteeringProblemWithPickupsAndDeliveries(string inputFile)
{

	countSolutionRuns=0;
	clock_start=clock();
	inputDataProcessor.init(inputFile);
	problemSize = inputDataProcessor.getProblemSize();
	maxCapacity= inputDataProcessor.getMaximumLoadCapacity();
	initProfitPerDistanceMatrix();
	calcPickupDeliveryPointPairs();
	unvisitedNodes.assign (problemSize, 1);
	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;
	numberOfTours=3;
	wrongPairs.resize(numberOfTours);
	load.resize(numberOfTours);
	bufferPlus.resize(numberOfTours);
	bufferMinus.resize(numberOfTours);
	goodsOnTheLorry.resize(numberOfTours);
	vector <int> twotimeszero;
	twotimeszero.assign(2,0);
	for (int i = 0; i < numberOfTours; i++)
	{
		wrongPairs[i].push_back(twotimeszero);
	}
	rexe="\"C:/Program Files/R/R-3.0.1/bin/R\"";
	rpath="C:/Users/User/Documents/Rfiles/";
	filename="empty";
	system ("del \"C:\\Users\\User\\Documents\\Rfiles\\tour*\""); 

	
}

OrienteeringProblemWithPickupsAndDeliveries::~OrienteeringProblemWithPickupsAndDeliveries()
{
	
}

bool OrienteeringProblemWithPickupsAndDeliveries::isZero (int i) {
  return (i==0);
}



void OrienteeringProblemWithPickupsAndDeliveries::initProfitPerDistanceMatrix()
{
	vector <double> initToZeroVector(problemSize,0);
	profitPerDistanceMatrix.assign(problemSize, initToZeroVector);
}

void OrienteeringProblemWithPickupsAndDeliveries::initIntensityMatrix()
{
	vector <double> initToZeroVector(solutionTours[0].size(),0);
	intensity.assign(solutionTours.size(), initToZeroVector);
}



void OrienteeringProblemWithPickupsAndDeliveries::calcPickupDeliveryPointPairs()
	{
	double travelDistance; // stores the calculated distances
	double profitValueMax; // profit value per unit which can be reached by visiting a certain pair of Pick up and Delivery point
	vector<Coordinates> basicData(inputDataProcessor.getBasicData());
	for (int i = 0; i < problemSize; i++)
	{
		for(int j = i+1; j < problemSize; j++)
		{
			if(basicData[i].quantity*basicData[j].quantity >= 0 || basicData[i].quantity<=0 )   
			{
				profitPerDistanceMatrix[i][j]=-DBL_MAX;
			}
			else
			{
				//travelDistance = inputDataProcessor.getDistance(i, j);
				if ( basicData[i].quantity>0)
				{
					profitPerDistanceMatrix[i][j]= (basicData[i].profit-basicData[j].profit)*min(basicData[i].quantity,-(basicData[j].quantity));
				}
				else
				{
					profitPerDistanceMatrix[i][j]=-DBL_MAX;
					//profitPerDistanceMatrix[i][j]= (-basicData[i].profit+basicData[j].profit)*min(-( basicData[i].quantity),basicData[j].quantity)/travelDistance;
				}

			}
		}
	}


			/*
			//travelDistance = inputDataProcessor.getDistance(i, j);
			if(basicData[i].quantity*basicData[j].quantity >= 0 && basicData[i].quantity<=0)   
			{
				profitPerDistanceMatrix[i][j]=-DBL_MAX;
			}
			else
			{
			//profitPerDistanceMatrix[i][j]= ((basicData[i].profit*basicData[min(i,j)].quantity)+(basicData[j].profit*basicData[min(i,j)].quantity))-travelDistance;
			profitPerDistanceMatrix[i][j]= ((basicData[i].profit*basicData[min(i,j)].quantity)+(basicData[j].profit*basicData[min(i,j)].quantity));
			}
		}
	}

	/*
	cout << "The profitPerDistance matrix is: "<< endl;
	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			cout<<[i][j] << "   " ;
		}
		cout<<endl;
	}
*/	
}



void OrienteeringProblemWithPickupsAndDeliveries::getProfitMatrixForPickupAndDeliveryPairs(int startNode, int whichTour)
	{
	double travelDistance; // stores the calculated distances
	double profitValueMax; // profit value per unit which can be reached by visiting a certain pair of Pick up and Delivery point
	vector<Coordinates> basicData(inputDataProcessor.getBasicData());
	

	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			if(basicData[i].quantity*basicData[j].quantity >= 0 || basicData[i].quantity<=0 || unvisitedNodes[i]==0)   
			{
				profitPerDistanceMatrix[i][j]=-DBL_MAX;
			}
			else
			{
			//travelDistance = inputDataProcessor.getDistance(startNode, i)+inputDataProcessor.getDistance(i, j)+inputDataProcessor.getDistance(j, 1);
			if ( basicData[i].quantity>0)
			{
			profitPerDistanceMatrix[i][j]= (basicData[i].profit-basicData[j].profit)*min(basicData[i].quantity,-(basicData[j].quantity));
			}
			else
			{
			profitPerDistanceMatrix[i][j]=-DBL_MAX;
			//profitPerDistanceMatrix[i][j]= (-basicData[i].profit+basicData[j].profit)*min(-( basicData[i].quantity),basicData[j].quantity)/travelDistance;
			}


			}
		}
	}
	/*
	if (!wrongPairs.empty())
	{
	for (int i = 0; i < wrongPairs[whichTour].size(); i++) // Nodes which were previously tried to be added to the tour
	{
		profitPerDistanceMatrix[wrongPairs[whichTour][i][0]] [wrongPairs[whichTour][i][1]]=-DBL_MAX;
	}
	
	}
	*/

	/*
	cout << "The profitPerDistance matrix is: "<< endl;
	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			cout<<profitPerDistanceMatrix[i][j] << "   " ;
		}
		cout<<endl;
	}
	*/	
}




void OrienteeringProblemWithPickupsAndDeliveries::getProfitMatrixForPickupAndDeliveryPairsParallel(int whichTour)
	{
	double travelDistance; // stores the calculated distances
	double profitValueMax; // profit value per unit which can be reached by visiting a certain pair of Pick up and Delivery point
	vector<Coordinates> basicData(inputDataProcessor.getBasicData());
	//double travelCostFactor;

	for (int i = 0; i < problemSize; i++)
	{
		for(int j = i+1; j < problemSize; j++)
		{
			if(basicData[i].quantity*basicData[j].quantity >= 0 || basicData[i].quantity<=0 || unvisitedNodes[i]==0)   
			{
				profitPerDistanceMatrix[i][j]=-DBL_MAX;
			}
			else
			{
			//travelCostFactor=inputDataProcessor.getDistance(0,i)+inputDataProcessor.getDistance(0,j)+inputDataProcessor.getDistance(i,1)+inputDataProcessor.getDistance(j,1);
			profitPerDistanceMatrix[i][j]= (basicData[i].profit-basicData[j].profit)* min(basicData[i].quantity,-(basicData[j].quantity));
			}
		}
	}
	/*
	for (int i = 0; i < wrongPairs[whichTour].size(); i++) // Nodes which were previously tried to be added to the tour
	{
		if ( wrongPairs[whichTour][i][0]!=0 && wrongPairs[whichTour][i][1]!=0)
		{
			profitPerDistanceMatrix[wrongPairs[whichTour][i][0]] [wrongPairs[whichTour][i][1]]=-DBL_MAX;
		}
	}
	*/
	/*
	cout << "The profitPerDistance matrix is: "<< endl;
	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			cout<<profitPerDistanceMatrix[i][j] << "   " ;
		}
		cout<<endl;
	}
	*/
}


void OrienteeringProblemWithPickupsAndDeliveries::getProfitMatrixForPickupAndDeliveryPairsParallelBest()
	{
	double travelCostFactor; // stores the calculated distances
	double profitValueMax; // profit value per unit which can be reached by visiting a certain pair of Pick up and Delivery point
	vector<Coordinates> basicData(inputDataProcessor.getBasicData());
	

	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			if(basicData[i].quantity*basicData[j].quantity >= 0 || basicData[i].quantity<=0 || unvisitedNodes[i]==0)   
			{
				profitPerDistanceMatrix[i][j]=-DBL_MAX;
			}
			else
			{
			travelCostFactor=inputDataProcessor.getDistance(i,j);
			profitPerDistanceMatrix[i][j]= (basicData[i].profit-basicData[j].profit)*min(basicData[i].quantity,-(basicData[j].quantity))/(travelCostFactor*2);
			//profitPerDistanceMatrix[i][j]= ((basicData[i].profit*basicData[min(i,j)].quantity)+(basicData[j].profit*basicData[min(i,j)].quantity));
			}
		}
	}

	/*
	for (int i = 0; i < wrongPairs[whichTour].size(); i++) // Nodes which were previously tried to be added to the tour
	{
		if ( wrongPairs[whichTour][i][0]!=0 && wrongPairs[whichTour][i][1]!=0)
		{
			profitPerDistanceMatrix[wrongPairs[whichTour][i][0]] [wrongPairs[whichTour][i][1]]=-DBL_MAX;
		}
	}

	*/

	/*
	cout << "The profitPerDistance matrix is: "<< endl;
	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			cout<<profitPerDistanceMatrix[i][j] << "   " ;
		}
		cout<<endl;
	}
	*/
}




int OrienteeringProblemWithPickupsAndDeliveries::getNearestNeighbor (vector<int> unvisitedCities, int startNode)
{
    double min = DBL_MAX;
    int nearest = startNode;
    for(int i=0; i<unvisitedCities.size(); i++)
    {
        if(unvisitedCities[i]==0) continue; 
        if (startNode == i) continue; 
        if(inputDataProcessor.getDistance(startNode, i)<min)
        {
            min = inputDataProcessor.getDistance(startNode, i);
            nearest = i;
        }
    }
    return nearest;
}

int OrienteeringProblemWithPickupsAndDeliveries::getHighestDemandInNeighborhood(vector<int> unvisitedCities, int startNode)
{
	double min = DBL_MAX;
    int nearest = startNode;
    for(int i=0; i<unvisitedCities.size(); i++)
    {
        if(unvisitedCities[i]==0) continue; 
        if (startNode == i) continue; 
        if(inputDataProcessor.getDistance(startNode, i)<min)
        {
            min = inputDataProcessor.getDistance(startNode, i);
            nearest = i;
        }
    }
    return nearest;
}





bool OrienteeringProblemWithPickupsAndDeliveries::isTotalLengthUnderLimit(vector<int> currentTour, int nodeToAdd){
	double result = 0;
	int lastNode = currentTour[0];
	currentTour.push_back(nodeToAdd);
	
	for(int i = 1; i < currentTour.size(); i++) { // The total length of the tour will be computed
		int currentNode = currentTour[i];
		result += inputDataProcessor.getDistance(lastNode, currentNode);
		lastNode = currentNode;
	}
	result= result+ inputDataProcessor.getDistance(lastNode, 1);
	return result<=inputDataProcessor.getMaximumTourLength(); // Is my tour length smaller than dMax. TRUE or FALSE will return

}


bool OrienteeringProblemWithPickupsAndDeliveries::isTotalLengthUnderLimit (vector<int> currentTour, int nodeToAdd, int position){
	double result = 0;
	int lastNode = currentTour[0];
	currentTour.insert(currentTour.begin()+position, nodeToAdd);
	
	for(int i = 1; i < currentTour.size(); i++) { // The total length of the tour will be computed
		int currentNode = currentTour[i];
		result += inputDataProcessor.getDistance(lastNode, currentNode);
		lastNode = currentNode;
	}
	//result= result+ inputDataProcessor.getDistance(lastNode, 1);
	return result<=inputDataProcessor.getMaximumTourLength(); // Is my tour length smaller than dMax. TRUE or FALSE will return

}

bool OrienteeringProblemWithPickupsAndDeliveries::isTotalLengthUnderLimit2NodesDifferentPlace (vector<int> currentTour, int nodeToAdd, int position, int node2ToAdd, int position2){
	double result = 0;
	int lastNode = currentTour[0];
	currentTour.insert(currentTour.begin()+position, nodeToAdd);
	currentTour.insert(currentTour.begin()+position2, node2ToAdd);

	for(int i = 1; i < currentTour.size(); i++) { // The total length of the tour will be computed
		int currentNode = currentTour[i];
		result += inputDataProcessor.getDistance(lastNode, currentNode);
		lastNode = currentNode;
	}
	//result= result+ inputDataProcessor.getDistance(lastNode, 1);
	return result<=inputDataProcessor.getMaximumTourLength(); // Is my tour length smaller than dMax. TRUE or FALSE will return

}








bool OrienteeringProblemWithPickupsAndDeliveries::isTotalLengthUnderLimit2Nodes(vector<int> currentTour, vector <int> bestPair){
	double result = 0;
	int lastNode = currentTour[0];

	for(int i = 0; i < bestPair.size(); i++) 
	{
		currentTour.push_back(bestPair[i]);
	}

	for(int i = 1; i < currentTour.size(); i++) { // The total length of the tour will be computed
		int currentNode = currentTour[i];
		result += inputDataProcessor.getDistance(lastNode, currentNode);
		lastNode = currentNode;
	}
	result= result+ inputDataProcessor.getDistance(lastNode, 1);
	return result<=inputDataProcessor.getMaximumTourLength(); // Is my tour length smaller than dMax. TRUE or FALSE will return

}


bool OrienteeringProblemWithPickupsAndDeliveries::isThereUnvisitedNodes()
{
	return (find(unvisitedNodes.begin(), unvisitedNodes.end(), 1) != unvisitedNodes.end());
}


void OrienteeringProblemWithPickupsAndDeliveries::profitsOfAllTheToursOhneGLPK()
{
	finalSolutions.clear();
	double result=0;
	double time=0;
	for (int i=0 ; i< solutionTours.size();i++)
	{
		if (intensity.size()>0)
		{
		ProfitCalculatorOhneGLPK profitohneGLPK( solutionTours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity(), getTourLength(solutionTours[i]), intensity[i], clock_start, clockStartThisSolution); 
		finalSolutions.push_back( profitohneGLPK.getProfit());
		finalSolutions.push_back (getTourLength(solutionTours[i]));
		cout << "Profit of the " << i+1 << ". tour ohne GLPK: " << profitohneGLPK.getProfit() << endl;
		result= result+ profitohneGLPK.getProfit();
		Rprintsol2(rexe, rpath, filename, profitohneGLPK, i, countSolutionRuns);
		time=profitohneGLPK.time_c;
		}

	}
	finalSolutions.push_back(result);
	finalSolutions.push_back(time);
	totalFinalSolutions.push_back(finalSolutions);
	cout << "The total profit of the tours all together :  " << result << endl;
	cout << "The tours:    "<< endl;
	for(int i = 0; i < solutionTours.size(); i++)
	{
		cout << " The tour nr. " << i+1 << "  " ;
		for(int j = 0; j < solutionTours[i].size(); j++)
		{
			cout<< solutionTours[i][j] << "   " ;
		}
		cout<<endl;

	}
}


void OrienteeringProblemWithPickupsAndDeliveries::profitsOfAllTheTours(int seedNumber, double timeStart)
{
	totalFinalSolutions.clear();
	intensity.clear();
	double result=0;
	double time=0;
	for (int i=0 ; i< solutionTours.size();i++)
	{
		finalSolutions.push_back(seedNumber);
		ProfitCalculator initialToursProfit( solutionTours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity(), getTourLength(solutionTours[i]), clock_start, clockStartThisSolution); 
		finalSolutions.push_back(i+1); // Tour
		finalSolutions.push_back( initialToursProfit.getProfit()); // Profit
		finalSolutions.push_back (getTourLength(solutionTours[i])); // Tourlength
		vector <double> intensityToInsert=  initialToursProfit.getIntensity();
		intensity.push_back(intensityToInsert);
		
		cout << "Profit of the " << i+1 << ". tour: " << initialToursProfit.getProfit() << endl;
		//Rprintsol(rexe, rpath, filename, initialToursProfit, i, countSolutionRuns);
		result= result+ initialToursProfit.getProfit();
		//time=initialToursProfit.time_c- timeStart;
		time=(clock()-timeStart)/CLOCKS_PER_SEC;
		finalSolutions.push_back(time);
		totalFinalSolutions.push_back(finalSolutions);
		finalSolutions.clear();

	}
	//RprintsolNEW(rexe, rpath, filename,countSolutionRuns, result);
	

	//finalSolutions.push_back(result);
	//finalSolutions.push_back(time);

	//totalFinalSolutions.push_back(finalSolutions);

	cout << "The total profit of the tours all together :  " << result << endl;
	cout << "The tours:    "<< endl;
	for(int i = 0; i < solutionTours.size(); i++)
	{
		cout << " The tour nr. " << i+1 << "  " ;
		for(int j = 0; j < solutionTours[i].size(); j++)
		{
			cout<< solutionTours[i][j] << "   " ;
		}
		cout<<endl;

	}





}

void OrienteeringProblemWithPickupsAndDeliveries::runTwoOpt(int seedNumber, double timeStart)
{
	for (int i = 0; i < solutionTours.size(); i++)
	{
		doTwoOpt(i);
	}
	//profitsOfAllTheToursOhneGLPK();
	profitsOfAllTheTours(seedNumber, timeStart);
}





/*
void OrienteeringProblemWithPickupsAndDeliveries::doTwoOpt(std::vector <int>  & tour)
{  
	int n, l;
	n= tour.size();
	l=floor((n-2)/2);
	bool criterion=true;
	while (criterion){
		criterion=false;
		for (int i = 0; i < l; i++) // length of sequence
		{
			if (criterion==true)
			{ 
				break;
			}
			for (int j=0; j<n; j++) 
			{
				if ( (inputDataProcessor.getDistance(tour[i+j], tour[n+i-1])+inputDataProcessor.getDistance(tour[i], tour[i+j+1])) <  (inputDataProcessor.getDistance(tour[i], tour[n+i-1])+inputDataProcessor.getDistance(tour[i+j], tour[i+j+1])));
				{
					std:: swap (tour[i], tour[j]);
					criterion=true;
				}

			}
		}
	}
	cout << "The tour after 2-opt: "<< endl;
		for(int j = 0; j < tour.size(); j++)
		{
			cout<<tour[j] << "   " ;
		}
		cout<<endl;
}
*/

/*void OrienteeringProblemWithPickupsAndDeliveries::doTwoOpt(std::vector <int>  & tour)
{  
	int n, l;
	n= tour.size();
	l=floor((n-2)/2);
	bool criterion=true;
	while (criterion){
		criterion=false;
		for (int i = 1; i < l; i++) // 
		{
			if (criterion==true)
			{ 
				break;
			}
	
			for (int j=i+1; j<n-2; j++) 
			{
				if ( (inputDataProcessor.getDistance(tour[i-1], tour[i+j])+inputDataProcessor.getDistance(tour[i], tour[i+j+1])) <  (inputDataProcessor.getDistance(tour[i-1], tour[i])+inputDataProcessor.getDistance(tour[i+j], tour[i+j+1])))
				{
					std:: swap (tour[i], tour[j]);
					criterion=true;
				}

			}
		}
	}
	cout << "The tour after 2-opt: "<< endl;
		for(int j = 0; j < tour.size(); j++)
		{
			cout<<tour[j] << "   " ;
		}
		cout<<endl;
}
*/


void OrienteeringProblemWithPickupsAndDeliveries::doTwoOpt(int whichTour)
{  
	cout << "The tour before 2-opt: "<< endl;
		for(int j = 0; j < solutionTours[whichTour].size(); j++)
		{
			cout<<solutionTours[whichTour][j] << "   " ;
		}
	cout<<endl;
	
	
	int n;
	n= solutionTours[whichTour].size();
	bool criterion=true;
	double temp;
	

	LoadCalculator loadCalculation(load[whichTour], goodsOnTheLorry[whichTour], bufferPlus[whichTour], bufferMinus[whichTour], intensity[whichTour], inputDataProcessor.getTourQuantities(solutionTours[whichTour]));

	while (criterion)
	{
		criterion=false;
		for (int i = 1; i < n-2; i++) // 
		{
			if (criterion==true)
			{ 
				break;
			}
	
			for (int j=i+1; j<n-1; j++) 
			{
				if ( (inputDataProcessor.getDistance(solutionTours[whichTour][i-1], solutionTours[whichTour][j])+inputDataProcessor.getDistance(solutionTours[whichTour][i], solutionTours[whichTour][j+1])) <  (inputDataProcessor.getDistance(solutionTours[whichTour][i-1], solutionTours[whichTour][i])+inputDataProcessor.getDistance(solutionTours[whichTour][j], solutionTours[whichTour][j+1])))
				{
					
					//LoadCalculator loadCalculation(load[whichTour], goodsOnTheLorry[whichTour], bufferPlus[whichTour], bufferMinus[whichTour], intensity[whichTour], inputDataProcessor.getTourQuantities(solutionTours[whichTour]));
					
					changeOrderPartOfTour(load[whichTour], i, j);
					vector <int> goodsOnLoad(goodsOnTheLorry[whichTour]);
					int feasible=1;
					goodsOnLoad[0]= 0;
					for (int k=i; k < goodsOnLoad.size(); k++) // for (int k=1; k < load[whichTour].size(); k++)
					{
						goodsOnLoad[k]=goodsOnLoad[k-1]+load[whichTour][k];
					}
					

					/*
					for (int k=i; k < j+1; k++)
					{
						goodsOnTheLorry[whichTour][k]=goodsOnTheLorry[whichTour][k-1] + load[whichTour][k];
					}
					*/

					for (int m=0; m < goodsOnLoad.size(); m++)
					{
						if ((goodsOnLoad[m]<0) || goodsOnLoad[m]>inputDataProcessor.getMaximumLoadCapacity())
						{
						changeOrderPartOfTour(load[whichTour], i, j);
						feasible=0; // Der Tausch ist nicht zulaessig
						break;
						}	
					}
		
					if (feasible==1) // Der Tausch ist zulässig
					{
						changeOrderPartOfTour(solutionTours[whichTour], i, j);
						changeOrderPartOfTourDouble(intensity[whichTour], i, j);
						for (int k=i; k < j+1; k++)
						{
							goodsOnTheLorry[whichTour][k]=goodsOnTheLorry[whichTour][k-1] + load[whichTour][k];
						}
						ProfitCalculatorOhneGLPK newProfit ( solutionTours[whichTour], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity(), getTourLength(solutionTours[whichTour]), intensity[whichTour], clock_start, clockStartThisSolution);
						//cout << "The new profit is:   " << newProfit.getProfit() << endl;
						
						criterion=true;
						break;
					}

				}

			}
			
		}
	}

	cout << "The tour after 2-opt: "<< endl;
		for(int j = 0; j < solutionTours[whichTour].size(); j++)
		{
			cout<<solutionTours[whichTour][j] << "   " ;
		}
		cout<<endl;
}

double OrienteeringProblemWithPickupsAndDeliveries::getTourLength(vector <int> tour)
{
	double tourLength=0;
	for(int i = 0; i < tour.size()-1; i++)
	{
		tourLength=tourLength+ inputDataProcessor.getDistance(tour[i],tour[i+1]);
	}
return tourLength;
}

void OrienteeringProblemWithPickupsAndDeliveries::changeOrderPartOfTour( vector <int> & tour, int from, int to)
{
	while (from<to)
	{
		int tmp = tour[from];
		tour[from]=tour [to];
		tour[to]=tmp;
		from++;
		to--;
	}
}

void OrienteeringProblemWithPickupsAndDeliveries::changeOrderPartOfTourDouble ( vector <double> & tour, int from, int to)
{
	while (from<to)
	{
		double tmp = tour[from];
		tour[from]=tour [to];
		tour[to]=tmp;
		from++;
		to--;
	}
}
void  OrienteeringProblemWithPickupsAndDeliveries::shaking(vector <int> & tour, int neighborhoodSize)
{
	srand (time(NULL));
	int numberOfNodesToErase;
	int posToErase;
	int posToInsert;
	vector <int> nodesToBeDisplaced;
	numberOfNodesToErase=(rand() % neighborhoodSize)+1;
	if ((tour.size()-2)>numberOfNodesToErase)
	{
		posToErase=rand()%(tour.size()-numberOfNodesToErase-1)+1;


		for (int i = 0; i < numberOfNodesToErase; i++)
		{
			nodesToBeDisplaced.push_back(tour[posToErase+i]);
		}

		tour.erase(tour.begin()+posToErase, tour.begin()+posToErase+numberOfNodesToErase);
	
		posToInsert=rand()%(tour.size()-1)+1;
	
		for (int i=0; i< nodesToBeDisplaced.size(); i++)
		{
			tour.insert(tour.begin()+posToInsert+i, nodesToBeDisplaced[i]);
			//TO AVOID INFEASIBLE SOLUTIONS - now not used,...at the end of the improvement I run a function 
			//which will reduce the tour until the solution is feasible
			/*if (getTourLength(tour)>inputDataProcessor.getMaximumTourLength())
			{
				tour.erase(tour.begin()+posToInsert+i);
				continue;
			}
			*/
		}
	}
	/*
	cout << "The tour after the shaking process: "<< endl;
	for(int j = 0; j < tour.size(); j++)
	{
		cout<<tour[j] << "   " ;
	}
	cout<<endl;
	*/
}

void OrienteeringProblemWithPickupsAndDeliveries::runShaking()
{
	shaking(solutionTours[0],3);
}

void OrienteeringProblemWithPickupsAndDeliveries::erasePoints(int neighborhoodSize, double percentageToErase)
{
		cout << "The tours before erasing points  "<< endl;
	for(int i = 0; i < solutionTours.size(); i++)
	{
		cout << " The tour nr. " << i+1 << "  " ;
		for(int j = 0; j < solutionTours[i].size(); j++)
		{
			cout<< solutionTours[i][j] << "   " ;
		}
		cout<<endl;

	}
	
		
	int posToErase;
	for (int i=0; i<solutionTours.size();i++)  // The local search with shaking will be done for each tour 
	{	
		LoadCalculator currentLoad(load[i],goodsOnTheLorry[i], bufferPlus[i], bufferMinus[i], intensity[i], inputDataProcessor.getQuantities());
		//vector <int> tour = solutionTours[i];
		for (int k=1 ; k< (max(floor(percentageToErase*solutionTours[i].size()/neighborhoodSize), 2)) ; k++)  // Neighborhood size is increasing
		{
			/* old code
			posToErase=rand()%(solutionTours[i].size()-neighborhoodSize-1)+1;
			*/
			int minLoad=INT_MAX;
			for (int l=1; l<solutionTours[i].size()-1;l++)
			{
				if (abs(load[i][l+1]) < minLoad)
				{
					minLoad=load [i] [l+1];
					posToErase=l;
				}
			}



			for (int m = 0; m < neighborhoodSize; m++)
			{
				unvisitedNodes[solutionTours[i][posToErase+m]]=1;
			}

			solutionTours[i].erase(solutionTours[i].begin()+posToErase, solutionTours[i].begin()+posToErase+neighborhoodSize);
			LoadCalculator newLoad(load[i],goodsOnTheLorry[i], bufferPlus[i], bufferMinus[i], intensity[i], inputDataProcessor.getQuantities());
		}
	}

	cout << "The tours after erasing points  "<< endl;
	for(int i = 0; i < solutionTours.size(); i++)
	{
		cout << " The tour nr. " << i+1 << "  " ;
		for(int j = 0; j < solutionTours[i].size(); j++)
		{
			cout<< solutionTours[i][j] << "   " ;
		}
		cout<<endl;

	}

}



void OrienteeringProblemWithPickupsAndDeliveries::insertPoints(vector <Coordinates> basicData, int whichTour)
{
	vector <int> tour=solutionTours[whichTour];
	//int best=DBL_MIN;
	int countTrials=0;
	vector <int> unvisitedNodesForOneTour = unvisitedNodes;
	int maxTrials=count(unvisitedNodesForOneTour.begin(), unvisitedNodesForOneTour.end(), 1);
	bool insertionEnd=false;
	unvisitedNodesForOneTour[0]=0;
	pickupPoint(unvisitedNodes);
	deliveryPoint(unvisitedNodes);
	//LoadCalculator loadCalculation(load[whichTour], goodsOnTheLorry[whichTour], bufferPlus[whichTour], bufferMinus[whichTour], intensity[whichTour], inputDataProcessor.getTourQuantities(solutionTours[whichTour]));
	getProfitMatrixForPickupAndDeliveryPairs(0 , whichTour);

	while (insertionEnd!=true)
	{
		int best=DBL_MIN;
		/*
		double minCost=DBL_MAX;
		for (int i=0; i< unvisitedNodesForOneTour.size(); i++)
		{
			if (unvisitedNodesForOneTour[i]==1)
			{
				if (basicData[i].quantity >0) // if it is a pick up point
				{
					if (abs(basicData[i].profit) < minCost)   // The cheaper it costs, the better it is
					{
					best= i;
					minCost=abs(basicData[i].profit);
					}
				}
			}
		}
		*/
		
		if(countTrials%2==1)
		{
		//pickupPoint(unvisitedNodesForOneTour);
		best=getNextPickupPointRandomised(pickupPoints, deliveryPoints, best, solutionTours[whichTour] );
		}
		//
		else
		{
		//deliveryPoint(unvisitedNodes);
		best=getNextDeliveryPointRandomised(pickupPoints, deliveryPoints, best,solutionTours[whichTour] );
		}

			// 2. Choosing the best place to insert in the other tour
			double minTourLengthExtension=DBL_MAX;
			double TourLengthExtension;
			int posToInsert;
			bool criterion=false;
			for (int i = 0; i < tour.size()-1; i++)
			{

				// pickup point to be inserted
				if((countTrials%2==1) && isTotalLengthUnderLimit( tour, best, i+1) && 0 < bufferPlus[whichTour][i+1] && best!=0  )
			
				{
					TourLengthExtension= inputDataProcessor.getDistance (tour[i], best)+ inputDataProcessor.getDistance (best, tour[i+1]);
					if (TourLengthExtension < minTourLengthExtension) 
					{
						minTourLengthExtension=TourLengthExtension;
						posToInsert= i+1;
					}
				}
				// delivery point to be inserted
				else
				{
					if((countTrials%2==0) && isTotalLengthUnderLimit( tour, best, i+1) &&  bufferMinus[whichTour][i+1]>0 && best!=0  )
			
					{
						TourLengthExtension= inputDataProcessor.getDistance (tour[i], best)+ inputDataProcessor.getDistance (best, tour[i+1]);
						if (TourLengthExtension < minTourLengthExtension) 
						{
							minTourLengthExtension=TourLengthExtension;
							posToInsert= i+1;
						}
					}
				
				
				
				}
			}
			if (minTourLengthExtension!=DBL_MAX && best!=0)
			{
				solutionTours[whichTour].insert(solutionTours[whichTour].begin()+posToInsert, best);
				intensity[whichTour].insert(intensity[whichTour].begin()+posToInsert, 1);
				LoadCalculator newLoad( load[whichTour], goodsOnTheLorry[whichTour], bufferPlus[whichTour], bufferMinus[whichTour], intensity[whichTour], inputDataProcessor.getQuantities());
				unvisitedNodes[best]=0;
				//insertionEnd=true;
				//countTrials=countTrials+1;	
				tour.insert(tour.begin()+posToInsert, best);
			}
						
			unvisitedNodesForOneTour[best]=0;
			countTrials=countTrials+1;		
		

		if (countTrials==maxTrials || maxTrials==0)
		{
		insertionEnd=true;
		}
	}
	/*
	cout << "The tours after inserting points  "<< endl;
	for(int i = 0; i < solutionTours.size(); i++)
	{
		cout << " The tour nr. " << i+1 << "  " ;
		for(int j = 0; j < solutionTours[i].size(); j++)
		{
			cout<< solutionTours[i][j] << "   " ;
		}
		cout<<endl;
	}
	*/
	

}

void OrienteeringProblemWithPickupsAndDeliveries::doInsertion(int seedNumber, double timeStart)
{
	for (int i=0; i < solutionTours.size(); i++)
	{
	insertPoints( inputDataProcessor.getBasicData(), i);
	}
	profitsOfAllTheTours(seedNumber, timeStart);
}




void OrienteeringProblemWithPickupsAndDeliveries::runShakingAndTwoopt()
{
	vector < vector <int>> bestTours;
	bestTours=solutionTours;
	int maxTrials=5;
	int improvement=0;
	for (int i=0; i<solutionTours.size();i++)  // The local search with shaking will be done for each tour 
	{	
		vector <int> tour = solutionTours[i];
		for (int k=1 ; k<11;k++)  // Neighborhood size is increasing
		{

			int trials=0;
			ProfitCalculator bestToursProfit ( bestTours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity(), getTourLength(bestTours[i]), clock_start, clockStartThisSolution);
			while (trials<maxTrials)
			{
				shaking(tour,k);
				makeTourFeasible (tour);
				ProfitCalculator currentTour ( bestTours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity(), getTourLength(bestTours[i]), clock_start, clockStartThisSolution);
				doTwoOpt(0);
				ProfitCalculator newTour( tour, inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity(), getTourLength(tour), clock_start, clockStartThisSolution);

				if (newTour.getProfit() > bestToursProfit.getProfit())
				{
					bestTours[i]=tour;
					improvement=improvement+1;
				}

				else
				{
					trials=trials+1;
				}
			}	
		}
	}
	for (int i=0 ; i< solutionTours.size();i++)
	{
		ProfitCalculator initialToursProfit( solutionTours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity(), getTourLength(solutionTours[i]), clock_start, clockStartThisSolution); 
		ProfitCalculator bestToursProfit( bestTours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity(), getTourLength(bestTours[i]), clock_start, clockStartThisSolution); 
		cout << "Profit of the " << i+1 << ". tour: " << initialToursProfit.getProfit() << endl;
		cout << "Profit of the  " << i+1 << ".tour after shaking and 2opt: " << bestToursProfit.getProfit() << endl;
		//getTourFeasible(bestTours[i]);
	}
	cout << "The number of improvements: " << improvement << endl;
}


double OrienteeringProblemWithPickupsAndDeliveries::getObjectiveValue(vector<int> tour)
{
	ProfitCalculator profit(tour, inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity(), getTourLength(tour), clock_start, clockStartThisSolution);
	double objectiveValue =profit.getProfit() - getTourLength(tour);
	return objectiveValue;
}


void OrienteeringProblemWithPickupsAndDeliveries::makeTourFeasible(vector <int> & solutionTour)
{
	vector <int> zeroIntensityIndices;
	ProfitCalculator Profit (solutionTour, inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity(), getTourLength(solutionTour), clock_start, clockStartThisSolution);
	zeroIntensityIndices=Profit.getZeroIntensityIndices();
	while (getTourLength(solutionTour) > inputDataProcessor.getMaximumTourLength())
	{
		if (zeroIntensityIndices.size()!=0)
		{
			solutionTour.erase(solutionTour.begin()+zeroIntensityIndices.back());
			zeroIntensityIndices.erase( zeroIntensityIndices.end()-1);
		}
		else
		{
			solutionTour.erase(solutionTour.begin()+ rand()%(solutionTour.size()-1)+1);
		}
	}
}



void OrienteeringProblemWithPickupsAndDeliveries::runStringExchangesAndTwoopt()
{
	vector < vector <int>> bestTours;
	bestTours=solutionTours;
	int maxTrials=1;
	int improvement=0;
	int tourNumber;
	for (int i=0; i<solutionTours.size()-1;i++)  // The local search with shaking will be done for each tour 
	{	
		tourNumber=i;
		vector <vector <int> > tours = solutionTours;
		for (int k=1 ; k<4;k++)  // Neighborhood size is increasing
		{

			int trials=0;
		ProfitCalculator bestToursProfit1 ( bestTours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity(), getTourLength(tours[i]), clock_start, clockStartThisSolution);
		ProfitCalculator bestToursProfit2 ( bestTours[i+1], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity(), getTourLength(tours[i+1]), clock_start, clockStartThisSolution);
			while (trials<maxTrials)
			{
				stringExchanges(tours, tourNumber, k);
				makeTourFeasible (tours[i+1]);
				doTwoOpt(1);
				ProfitCalculator newTour1( tours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity(), getTourLength(tours[i]), clock_start, clockStartThisSolution);
				ProfitCalculator newTour2( tours[i+1], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity(), getTourLength(tours[i+1]), clock_start, clockStartThisSolution);

				if (( newTour2.getProfit() - bestToursProfit2.getProfit())> 0)
				{
					bestTours[i+1]=tours [i+1];
					bestTours[i]=tours[i];
					improvement=improvement+1;
				}

				else
				{
					trials=trials+1;
				}
			}	
		}
	}

	for (int i=0 ; i< solutionTours.size();i++)
	{
		ProfitCalculator initialToursProfit( solutionTours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity(), getTourLength(solutionTours[i]), clock_start, clockStartThisSolution); 
		ProfitCalculator bestToursProfit( bestTours[i], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity(), getTourLength(bestTours[i]), clock_start, clockStartThisSolution); 
		cout << "Profit of the " << i+1 << ". tour: " << initialToursProfit.getProfit() << endl;
		cout << "Profit of the  " << i+1 << ".tour after shaking and 2opt: " << bestToursProfit.getProfit() << endl;
		//getTourFeasible(bestTours[i]);
	}
	cout << "The number of improvements: " << improvement << endl;
}


void  OrienteeringProblemWithPickupsAndDeliveries::stringExchanges (vector < vector <int>> & tours, int tourNumber, int neighborhoodSize)
{
	srand (time(NULL));
	int numberOfNodesToErase;
	int posToErase;
	int posToInsert;
	vector <int> nodesToBeDisplaced;
	numberOfNodesToErase=(rand() % neighborhoodSize)+1;

	//for (int tourNumber=0; tourNumber < solutionTours.size()-1; tourNumber+= 2) // tourNumber is the index of the tour from the solutionTours from which nodes will be erased
	//{
		if ((tours[tourNumber].size()-2)>numberOfNodesToErase)
		{
			
			// 1. Erase
			posToErase=rand()%(tours[tourNumber].size()-numberOfNodesToErase-1)+1;
			for (int i = 0; i < numberOfNodesToErase; i++)
			{
				nodesToBeDisplaced.push_back(tours[tourNumber][posToErase+i]);
			}

			tours[tourNumber].erase(tours[tourNumber].begin()+posToErase, tours[tourNumber].begin()+posToErase+numberOfNodesToErase);
	
			
			// 2. Choosing the best place to insert in the other tour
			double minTourLength=DBL_MAX;
			double TourLengthExtension;
			for (int i = 0; i < tours[tourNumber+1].size()-1; i++)
			{
				TourLengthExtension= getTourLength (nodesToBeDisplaced) + inputDataProcessor.getDistance (tours[tourNumber+1][i], nodesToBeDisplaced[0])+ inputDataProcessor.getDistance (tours[tourNumber+1][i+1], nodesToBeDisplaced.back());
				if (TourLengthExtension < minTourLength)
				{
					minTourLength=TourLengthExtension;
					posToInsert= i+1;
				}
			
			}

			// 3. Cheapest Insertion

			for (int i=0; i< nodesToBeDisplaced.size(); i++)
			{
				tours[tourNumber+1].insert(tours[tourNumber+1].begin()+posToInsert+i, nodesToBeDisplaced[i]);
			}
			//ProfitCalculator tourInwhichNodesAddedNewProfit ( solutionTours[tourNumber+1], inputDataProcessor.getBasicData(), inputDataProcessor.getMaximumLoadCapacity());
		}
	//}
	/*
	cout << "The tour after the shaking process: "<< endl;
	for(int j = 0; j < tour.size(); j++)
	{
		cout<<tour[j] << "   " ;
	}
	cout<<endl;
	*/
}


void  OrienteeringProblemWithPickupsAndDeliveries::Rprintsol(string rexe, string rpath, string filename, ProfitCalculator solution, int i, int & countSolutionRuns){
	
	countSolutionRuns=countSolutionRuns+1;
	solution.savesol( rpath + "example1.txt");
	inputDataProcessor.Rsavesinstance(rpath + "example2.txt");		
	string syscommand("\"C:\\Program Files\\R\\R-3.0.1\\bin\\R\" -f C:\\Users\\User\\Documents\\Rfiles\\start.r");
	//replace(syscommand.begin(),syscommand.end(),'/','\\');
	system(syscommand.c_str());

	string rename ("rename \"C:\\Users\\User\\Documents\\Rfiles\\test123.pdf\" \"tour" + std::to_string(countSolutionRuns) + "_" + std::to_string(i+1) + ".pdf\"");
	system( rename.c_str());
	
	
double time_c=(double)(clock()- clock_start)/CLOCKS_PER_SEC;
cout << time_c;

}


void  OrienteeringProblemWithPickupsAndDeliveries::Rprintsol2(string rexe, string rpath, string filename, ProfitCalculatorOhneGLPK solution, int i, int & countSolutionRuns){
	
	countSolutionRuns=countSolutionRuns+1;
	solution.savesol( rpath + "example1.txt");
	inputDataProcessor.Rsavesinstance(rpath + "example2.txt");		
	string syscommand("\"C:\\Program Files\\R\\R-3.0.1\\bin\\R\" -f C:\\Users\\User\\Documents\\Rfiles\\start.r");
	//replace(syscommand.begin(),syscommand.end(),'/','\\');
	system(syscommand.c_str());

	string rename ("rename \"C:\\Users\\User\\Documents\\Rfiles\\test123.pdf\" \"tour" + std::to_string(countSolutionRuns) + "_" + std::to_string(i+1) + ".pdf\"");
	system( rename.c_str());
	
	
double time_c=(double)(clock()- clock_start)/CLOCKS_PER_SEC;
cout << time_c;

}




void  OrienteeringProblemWithPickupsAndDeliveries::printSolutions()
{
	for (int i=0; i<totalFinalSolutions.size();i++)
	{
		cout << "The profits of the " << i+1 << ". tour is: " << endl;
		for (int j=0; j<totalFinalSolutions[i].size();j++)
		{
			cout << totalFinalSolutions[i][j] << " " ;
		}
		cout << endl;
	}
}


int OrienteeringProblemWithPickupsAndDeliveries::getNextPickupPointRandomised(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode, vector <int> tour)
{
    double min = DBL_MAX;
    int bestToChoose = startNode;
	double profitsWithOtherDeliveryPoints=0;
	double profitValue;
	vector <double> cumulatedProfits;
	cumulatedProfits.push_back(0);
	double randNumber;
	vector<Coordinates> basicData(inputDataProcessor.getBasicData());
    for(int i=0; i< unvisitedPickups.size(); i++)
    {
         //for(int j=0; j< unvisitedDeliveries.size(); j++)
		//{
			profitValue=-1 / double(basicData[unvisitedPickups[i]]. profit); // inputDataProcessor. getDistance (unvisitedPickups[i], tour[closestPoint(unvisitedPickups[i],tour)]); // /inputDataProcessor.getDistance(startNode, unvisitedPickups[i]);	
			
			if (profitValue>0)
			{
			profitsWithOtherDeliveryPoints= profitsWithOtherDeliveryPoints + profitValue;
			}
		//}
		if (profitsWithOtherDeliveryPoints>0)
		{
			cumulatedProfits.push_back(cumulatedProfits.back()+profitsWithOtherDeliveryPoints);	
		}
		else
		{
			cumulatedProfits.push_back(cumulatedProfits.back());
		}
		
	
		profitsWithOtherDeliveryPoints=0;
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


int OrienteeringProblemWithPickupsAndDeliveries::getNextDeliveryPointRandomised(vector <int> unvisitedPickups, vector <int> unvisitedDeliveries, int startNode, vector <int> tour)
{
   
    int bestToChoose = 0;
	double profitsToGet=0;
	double profitValue=0;
	double quantityToTransfer;
	vector <double> cumulatedDistances;
	double travelDistance;
	cumulatedDistances.push_back(0);
	double randNumber;
	vector<Coordinates> basicData(inputDataProcessor.getBasicData());
	double max;
    for(int i=0; i< unvisitedDeliveries.size(); i++)
    {
		max=0;
		 //for(int j=0; j< unvisitedPickups.size(); j++)
		//{
			//if (inputDataProcessor.getDistance(closestPoint(unvisitedPickups[j],tour), unvisitedPickups[j])< inputDataProcessor.getDistance(closestPoint(unvisitedDeliveries[i],tour), unvisitedDeliveries[i]))
			//{
			profitValue=-1* double(basicData[unvisitedDeliveries[i]].profit); // inputDataProcessor.getDistance(closestPoint(unvisitedDeliveries[i], tour),unvisitedDeliveries[i]);  // / (inputDataProcessor.getDistance(closestPoint(unvisitedDeliveries[i],tour), unvisitedDeliveries[i]));
			//}
			if (profitValue>max)
			{
				//max=profitPerDistanceMatrix[unvisitedPickups[i]][unvisitedDeliveries[j]]/(inputDataProcessor.getDistance(closestPoint(unvisitedPickups[i],tour), unvisitedPickups[i]));
				max=profitValue;
			}
		//}
		
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



void OrienteeringProblemWithPickupsAndDeliveries::pickupPoint (std::vector<int> & nodes)
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


void OrienteeringProblemWithPickupsAndDeliveries::deliveryPoint (std::vector<int> & nodes)
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


int OrienteeringProblemWithPickupsAndDeliveries::closestPoint( int node, vector <int> tour)
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

void OrienteeringProblemWithPickupsAndDeliveries::runExcelExporter()
{
ExcelExporter excelExporter(totalFinalSolutions);
};


void OrienteeringProblemWithPickupsAndDeliveries::runExcelExportStart()
{
	ofstream MyExcelFile;
	MyExcelFile.open("C:\\Users\\User\\Documents\\Rfiles\\example.csv", ios::out);
	MyExcelFile << "Instance"<< ";";
	MyExcelFile << "Seed" << ";";
	MyExcelFile << "Tour" << ";";

	//for(int i=0; i<numberOfTours; i++)
	//{ 
	MyExcelFile << "Profit"<< ";";
	MyExcelFile << "Tourlength" << ";";
	//}

	//MyExcelFile << "Total Profit" << ";";
	MyExcelFile << "Run time" << ";";
	MyExcelFile << "Method" << ";" << endl;
}

void OrienteeringProblemWithPickupsAndDeliveries::runExcelExportFinish()
{
	ofstream MyExcelFile;
	MyExcelFile.close();
}

void OrienteeringProblemWithPickupsAndDeliveries::runExcelExport(string inputFile, string heurName)
{
	
	ofstream MyExcelFile;
	MyExcelFile.open("C:\\Users\\User\\Documents\\Rfiles\\example.csv", ios::out | ios::app);

	int rowCount=totalFinalSolutions.size();
	int colCount=totalFinalSolutions[0].size();
	/*
	MyExcelFile << "Instance"<< ";";
	MyExcelFile << "Seed" << ";";
	MyExcelFile << "Tour" << ";";

	//for(int i=0; i<numberOfTours; i++)
	//{ 
	MyExcelFile << "Profit"<< ";";
	MyExcelFile << "Tourlength" << ";";
	//}

	//MyExcelFile << "Total Profit" << ";";
	MyExcelFile << "Run time" << ";" << endl;
	*/
	 for(int i=0; i<rowCount; i++)
	 { 
		MyExcelFile << "=\"" << inputFile << "\"" << ";";
		//MyExcelFile << "=\"" << seedNumber << "\"" << ";";
		for(int j=0; j<colCount; j++) 
		{ 
			MyExcelFile << "=\"" << totalFinalSolutions[i][j] << "\"" << ";";
		}
		MyExcelFile << "=\"" << heurName << "\"" << ";";
		MyExcelFile	<< endl;
	 }
	 totalFinalSolutions.clear();
	//MyExcelFile.close();

};



void OrienteeringProblemWithPickupsAndDeliveries::saveSolution(string fname, double result)
{
		vector<Coordinates> basicData(inputDataProcessor.getBasicData());
		ofstream myfile;
		myfile.open (fname);
		myfile << "result = list() " << endl;


		myfile << "result$profit =  " << result << endl;
		myfile << "result$bound = " << -1 << endl;
		myfile << "result$length = list()" << endl;
		myfile << "result$length[[1]] = " << inputDataProcessor.getMaximumTourLength() << endl;
		myfile << "result$length[[2]] = " << inputDataProcessor.getMaximumTourLength() << endl;
		myfile << "result$length[[3]] = " << inputDataProcessor.getMaximumTourLength() << endl;
		myfile << "result$maxload = list()" << endl;
		myfile << "result$maxload[[1]] = " << inputDataProcessor.getMaximumLoadCapacity() << endl;
		myfile << "result$maxload[[2]] = " << inputDataProcessor.getMaximumLoadCapacity() << endl;
		myfile << "result$maxload[[3]] = " << inputDataProcessor.getMaximumLoadCapacity() << endl;
		myfile << "result$cputime = " << 0 << endl;//CPUtime
		myfile << "result$cputimetotal = " << 0 << endl; //CPUtimetotal
		myfile << "result$iterations = " << 0 << endl;
		myfile << "result$tour = list()" << endl;
		myfile << "result$intensity = list()" << endl;
		myfile << "result$quantity = list()" << endl;
		for (int m=0; m<numberOfTours;m++)
		{
			myfile << "result$tour[[" << m+1 << "]]= c(";
			for(int i = 0; i < solutionTours[m].size()-1; i++){
				myfile << 1 + solutionTours[m][i] << ", ";
			}
			myfile << 1 + solutionTours[m][solutionTours[m].size()-1] << ")" << endl;
			myfile << "result$intensity[[" << m+1 << "]]= c(";
			for(int i = 0; i < intensity[m].size()-1; i++){
				myfile << intensity[m][i] << ", ";
			}
			myfile << intensity[m][intensity[m].size()-1] << ")" << endl;
			myfile << "result$quantity[[" << m+1 << "]]= c(";
			for(int i = 0; i < solutionTours[m].size()-1; i++){
				myfile << intensity[m][i]* basicData[solutionTours[m][i]].quantity << ", ";
			}
			myfile << intensity[m][intensity[m].size()-1]* basicData[solutionTours[m][solutionTours[m].size()-1]].quantity << ")" << endl;
		}
		myfile.close();
	};

void  OrienteeringProblemWithPickupsAndDeliveries::RprintsolNEW(string rexe, string rpath, string filename, int & countSolutionRuns, double result){
	
	countSolutionRuns=countSolutionRuns+1;
	saveSolution( rpath + "example1.txt", result);
	inputDataProcessor.Rsavesinstance(rpath + "example2.txt");		
	string syscommand("\"C:\\Program Files\\R\\R-3.0.1\\bin\\R\" -f C:\\Users\\User\\Documents\\Rfiles\\start.r");
	//replace(syscommand.begin(),syscommand.end(),'/','\\');
	system(syscommand.c_str());

	string rename ("rename \"C:\\Users\\User\\Documents\\Rfiles\\test123.pdf\" \"tour" + std::to_string(countSolutionRuns) + ".pdf\"");
	system( rename.c_str());
	
	
double time_c=(double)(clock()- clock_start)/CLOCKS_PER_SEC;
cout << time_c;

}





