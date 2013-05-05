#ifndef __INPUTDATAPROCESSOR__
#define __INPUTDATAPROCESSOR__

#include <vector>
#include <string>
#include "Coordinates.h"


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

protected:
	double dMax; // maximum tour length

private:
	void readData(const std::string& inputFile);
	void initDistanceMatrix();
	void calcDistanceMatrix();

	
	double lMax; // maximum load capacity
	std::vector<Coordinates> basicData;
	int problemSize;
	std::vector < std::vector <double> > distanceMatrix;
};

#endif