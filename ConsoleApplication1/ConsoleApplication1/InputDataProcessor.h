#ifndef __INPUTDATAPROCESSOR__
#define __INPUTDATAPROCESSOR__

#include <vector>
#include <string>
#include "Coordinates.h"
using namespace std;


class InputDataProcessor
{
public:
	InputDataProcessor();
	void init(const std::string& inputFile);
	double getMaximumTourLength();
	double getMaximumLoadCapacity();
	const std::vector<Coordinates>& getBasicData();
	int getProblemSize();
	double getDistance(int from, int to);
	void Rsavesinstance(std::string fname);
	std::vector<Coordinates> basicData;
	vector <int> getQuantities();
	vector <int> getTourQuantities(vector <int> tour);
	int getQuantity(int node);


protected:
	double dMax; // maximum tour length
	vector <int> quantities;
		vector <int> tourquantities;

private:
	void readData(const std::string& inputFile);
	void initDistanceMatrix();
	void calcDistanceMatrix();

	
	double lMax; // maximum load capacity
	
	int problemSize;
	std::vector < std::vector <double> > distanceMatrix;
};

#endif