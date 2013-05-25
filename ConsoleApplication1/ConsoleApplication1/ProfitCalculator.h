#ifndef __PROFITCALCULATOR__
#define __PROFITCALCULATOR__


#include "Coordinates.h"
#include <vector>
#include <string>
using namespace std;

class ProfitCalculator
{
	public:
	ProfitCalculator(std::vector <int> tourInput,std::vector <Coordinates> basicDataInput, double maxCapacityInput);
	double getProfit();
	std:: vector <int> getZeroIntensityIndices();
	void savesol(string fname);

	private:
	void calculateProfit();
	double intensityCalculation();

	std:: vector<int> tour;
	double maxCapacity;
	std:: vector<Coordinates> basicData;
	std::vector <double> vectorMaxQuantities;
	std::vector <double> maxProfits;
	std::vector<double> intensity;
	double result;
	int upperbound;
	double CPUtime;
	int Iteratations;
	double CPUtimetotal;



};


#endif