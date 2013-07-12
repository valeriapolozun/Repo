#ifndef __PROFITCALCULATOROHNEGLPK__
#define __PROFITCALCULATOROHNEGLPK__


#include "Coordinates.h"
#include "time.h"
#include <vector>
#include <string>
using namespace std;



class ProfitCalculatorOhneGLPK
{
	public:
		ProfitCalculatorOhneGLPK(std::vector <int> tourInput,std::vector <Coordinates> basicDataInput, double maxCapacityInput, double tourLength, vector <double> intensityInput, clock_t clock_start, clock_t clockStartThisSolution);
		double getProfit();
		void calculateProfit();
		void calculateProfit2(vector<int> tour, vector<double> intensity);
		void savesol(string fname);
		double time_c;

private:
	
	std:: vector<int> tour;
	double tourLength;
	std:: vector<double> intensity;
	double maxCapacity;
	std:: vector<Coordinates> basicData;
	std::vector <double> vectorMaxQuantities;
	std::vector <double> maxProfits;
	double result;
	int upperbound;
	double CPUtime;
	int Iteratations;
	double CPUtimetotal;
	
};


#endif


