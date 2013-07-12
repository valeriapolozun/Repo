#ifndef __PROFITCALCULATOR__
#define __PROFITCALCULATOR__


#include "Coordinates.h"
#include "time.h"
#include <vector>
#include <string>
using namespace std;

class ProfitCalculator
{
	public:
	double time_c;
	ProfitCalculator(std::vector <int> tourInput,std::vector <Coordinates> basicDataInput, double maxCapacityInput, double tourLength, clock_t clock_start, clock_t clockStartThisSolution );
	double getProfit();
	std:: vector <int> getZeroIntensityIndices();
	void savesol(string fname);

	void loadCalculation( vector <int> & load, vector <int> & goodsOnTheLorry);
	void bufferPlusCalculation(vector <int> goodsOnTheLorry, vector <int> & bufferPlus);
	void bufferMinusCalculation(vector <int> goodsOnTheLorry, vector <int> & bufferMinus);
	vector <double> getIntensity();

	private:
	void calculateProfit();
	double intensityCalculation();	

	std:: vector<int> tour;
	double tourLength;
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