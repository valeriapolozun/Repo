#ifndef __PROFITCALCULATOR__
#define __PROFITCALCULATOR__


#include "Coordinates.h"
#include <vector>

class ProfitCalculator
{
	public:
	ProfitCalculator(std::vector <int> tourInput,std::vector <Coordinates> basicDataInput, double maxCapacityInput);
	double getProfit();
	std:: vector <int> getZeroIntensityIndices();

	private:
	void calculateProfit();
	double intensityCalculation();

	std:: vector<int> tour;
	double maxCapacity;
	std:: vector<Coordinates> basicData;
	std::vector <double> initialTour; // TODO change name
	std::vector <double> maxProfits;
	std::vector<double> intensity;
	double result;
};



#endif