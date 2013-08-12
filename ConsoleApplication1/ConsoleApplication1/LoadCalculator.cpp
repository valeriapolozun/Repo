
#include "LoadCalculator.h"
#include "InputDataProcessor.h"
#include "OrienteeringProblemWithPickupsAndDeliveries.h"

using namespace std;

LoadCalculator:: LoadCalculator (vector <int> & load, vector <int> & goodsOnTheLorry, vector <int> & bufferPlus, vector <int> & bufferMinus, vector<double> intensity, vector <int> quantities)
{
	load.clear();
	goodsOnTheLorry.clear();
	bufferPlus.clear();
	bufferMinus.clear();

	//maxCapacity = & OrienteeringProblemWithPickupsAndDeliveries::maxCapacity;
	maxCapacity= 300;
	loadCalculation (load, goodsOnTheLorry, intensity, quantities);
	bufferPlusCalculation( goodsOnTheLorry, bufferPlus);
	bufferMinusCalculation (  goodsOnTheLorry, bufferMinus);
}



void LoadCalculator::loadCalculation(vector <int> & load, vector <int> & goodsOnTheLorry, vector <double> intensity, vector <int> quantities)
{
	load.clear();
	load.push_back(0);
	goodsOnTheLorry.push_back(0);
	for (int i=1; i < intensity.size(); i++)
	{
		 load.push_back(intensity[i]*quantities[i]);
		 goodsOnTheLorry.push_back( goodsOnTheLorry[i-1]+ load[i]);
	}
}


void LoadCalculator::bufferPlusCalculation(vector <int> goodsOnTheLorry, vector <int> & bufferPlus)
{
	double min;
	for (int i=0; i < goodsOnTheLorry.size(); i++)
	{
		min= maxCapacity;
		for (int j=i; j < goodsOnTheLorry.size(); j++)
		{
			if (min> maxCapacity-goodsOnTheLorry[j])
			{
				min=maxCapacity-goodsOnTheLorry[j];
			}
		}	
		bufferPlus.push_back(min);
		
	}
}

void LoadCalculator::bufferMinusCalculation ( vector <int> goodsOnTheLorry, vector <int> & bufferMinus)
{
	double max;
	for (int i=0; i < goodsOnTheLorry.size(); i++)
	{
		max=0;
		for (int j=i; j < goodsOnTheLorry.size(); j++)
		{
			if (max < goodsOnTheLorry[j])
			{
				max=goodsOnTheLorry[j];
			}
		}
		
		 bufferMinus.push_back(max);
	}	
}

