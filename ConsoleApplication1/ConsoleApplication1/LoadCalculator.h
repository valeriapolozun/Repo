
#ifndef __LOADCALCULATOR__
#define __LOADCALCULATOR__


#include "Coordinates.h"
#include <vector>
#include <string>

using namespace std;


class LoadCalculator
{
	public:
	LoadCalculator (std::vector <int> & load, std::vector <int> & goodsOnTheLorry, std::vector <int> & bufferPlus, std:: vector <int> & bufferMinus, std:: vector <double> intensity, vector <int> quantities);
	void loadCalculation (vector <int> & load, vector <int> & goodsOnTheLorry, vector <double> intensity, vector <int> quantities);
	vector <int> load;
	void bufferPlusCalculation(vector <int> goodsOnTheLorry, vector <int> & bufferPlus);
	void bufferMinusCalculation ( vector <int> goodsOnTheLorry, vector <int> & bufferMinus);




	private:
	std:: vector<Coordinates> basicData;
	vector <int> tour;
	vector <double> intensity;
	vector <double> maxQuantities;
	double maxCapacity;

};

#endif

