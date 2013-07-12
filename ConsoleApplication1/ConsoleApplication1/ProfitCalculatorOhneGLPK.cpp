#include "ProfitCalculatorOhneGLPK.h"
#include "InputDataProcessor.h"
#include "OrienteeringProblemWithPickupsAndDeliveries.h"
#include <iostream>
#include <fstream>
#include <time.h>

using namespace std;


ProfitCalculatorOhneGLPK:: ProfitCalculatorOhneGLPK (std::vector <int> tourInput,std::vector <Coordinates> basicDataInput, double maxCapacityInput, double tourLengthInput, vector <double> intensityInput,clock_t clock_start, clock_t clockStartThisSolution ) 
{
	tour=tourInput;
	tourLength=tourLengthInput;
	intensity=intensityInput;
	maxCapacity= maxCapacityInput;
	for (int i = 0; i < tour.size(); i++)
	{
		basicData.push_back(basicDataInput[tour[i]]);
	}
	result=0;
	calculateProfit2(tour, intensity);
	
	result=result-tourLength ;
	time_c= ((clock()- clockStartThisSolution)*1000)/CLOCKS_PER_SEC/1000;
	CPUtimetotal=(double)((clock()- clock_start)*1000)/CLOCKS_PER_SEC/1000;
}


void ProfitCalculatorOhneGLPK:: calculateProfit()
{
	
	for (int i = 0; i < basicData.size(); i++)
	{
		vectorMaxQuantities.push_back(basicData[i].quantity);
	}


	for (int i = 0; i < basicData.size(); i++)
	{
		maxProfits.push_back(basicData[i].quantity*basicData[i].profit);
	}


	
	if (basicData.size()>2)
	{
		for(int i = 0; i < basicData.size(); i++) 
		{
			result += basicData[i].quantity*basicData[i].profit*intensity[i];
		}
	}

}


void ProfitCalculatorOhneGLPK:: calculateProfit2( vector<int> tour, vector <double> intensity)
{
	
	for (int i = 0; i < tour.size(); i++)
	{
		vectorMaxQuantities.push_back( basicData[i].quantity);
	}


	for (int i = 0; i < tour.size(); i++)
	{
		maxProfits.push_back(basicData[i].quantity*basicData[i].profit);
	}


	
	if (basicData.size()>2)
	{
		for(int i = 0; i < basicData.size(); i++) 
		{
			result += basicData[i].quantity*basicData[i].profit*intensity[i];
		}
	}

}




void ProfitCalculatorOhneGLPK::savesol(string fname)
{
		ofstream myfile;
		myfile.open (fname);
		myfile << "result = list() " << endl;


		myfile << "result$profit =  " << result << endl;
		myfile << "result$bound = " << upperbound << endl;
		myfile << "result$length = " << 50 << endl;
		myfile << "result$maxload = " << 300 << endl;
		myfile << "result$cputime = " << time_c << endl; // CPU time of this solution
		myfile << "result$cputimetotal = " << CPUtimetotal << endl; // CPU time total, TODO:  correct total time
		myfile << "result$iterations = " << Iteratations << endl;
		myfile << "result$tour = c(";
		for(int i = 0; i < tour.size()-1; i++){
			myfile << 1 + tour[i] << ", ";
		}
		myfile << 1 + tour[tour.size()-1] << ")" << endl;
		myfile << "result$intensity = c(";
		for(int i = 0; i < intensity.size()-1; i++){
			myfile << intensity[i] << ", ";
		}
		myfile << intensity[intensity.size()-1] << ")" << endl;
		myfile << "result$quantity = c(";
		for(int i = 0; i < tour.size()-1; i++){
			myfile << intensity[i]* basicData[i].quantity << ", ";
		}
		myfile << intensity[tour.size()-1]*basicData[tour.size()-1].quantity << ")" << endl;
		myfile.close();
	};


double ProfitCalculatorOhneGLPK::getProfit()
{
	return result;
}