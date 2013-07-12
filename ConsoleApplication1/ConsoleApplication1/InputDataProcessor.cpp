#include "InputDataProcessor.h"
#include <iostream>
#include <fstream>
using namespace std;

int InputDataProcessor::getProblemSize()
{
	return problemSize;
}

double InputDataProcessor::getDistance(int from, int to)
{
	return distanceMatrix[from][to];
}

const vector<Coordinates>& InputDataProcessor::getBasicData()
{
	return basicData;
}

double InputDataProcessor::getMaximumTourLength()
{
	return dMax;
}

double InputDataProcessor::getMaximumLoadCapacity()
{
	return lMax;
}

vector <int> InputDataProcessor::getQuantities()
{
	for (int i=0; i<basicData.size();i++)
	{
	quantities.push_back(basicData[i].quantity);
	}
	return quantities;
}

vector <int> InputDataProcessor::getTourQuantities(vector <int> tour)
{
	for (int i=0; i<tour.size();i++)
	{
	tourquantities.push_back(basicData[tour[i]].quantity);
	}
	return tourquantities;
}





int InputDataProcessor::getQuantity(int node)
{
	return basicData[node].quantity;
}





InputDataProcessor::InputDataProcessor()
{
}

void InputDataProcessor::init(const string& inputFile)
{
		readData(inputFile);
		problemSize=basicData.size();
		initDistanceMatrix();
	    calcDistanceMatrix ();
}





void InputDataProcessor::readData(const string& inputFile)
{
	Coordinates c;
	char line [80];

	ifstream file (inputFile); //open file

	if (!file)
	{ 
		cout << "Error opening file!"<<endl;
	}

	file >> dMax; // read in dMax
	file >> lMax; // read in load capacity
	file.getline(line, 80);

	while ( ! file.fail())
	{
		file >> c.x; // read in x coordinates
		file >> c.y; // read in y coordinates
		file >> c.quantity; // read in quantity
		file >> c.profit; //read in profit
		file.getline(line,80); //escape the recent line
		//note: the first and second coordinates are the staring and ending points
			basicData.push_back(c);
	}
	
	file.close();
	for(int i=0; i<basicData.size();i++)
	{
		cout<<basicData[i].x<<"  "<<basicData[i].y<<"  "<<basicData[i].quantity<< basicData[i].profit<<endl;
	}

}

		


void InputDataProcessor::initDistanceMatrix()
{
	vector <double> initToZeroVector(problemSize,0);
	distanceMatrix.assign(problemSize, initToZeroVector);
}


void InputDataProcessor::calcDistanceMatrix()
{
	double dx;     // xi.cood - xj.cood
	double dy;     // yi.cood - yj.cood
	double travelDistance;  // stores the calculated time

	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			if (i == j) travelDistance = DBL_MAX;    // define the distance from itself is a big number, so that itself can not be the nearest neighbor
			else
			{
				dx=basicData[i].x-basicData[j].x;   
				dy=basicData[i].y-basicData[j].y;
				travelDistance = sqrt (pow(dx, 2)+pow(dy,2));  //calculating the distance		    
			}
			distanceMatrix[i][j] = travelDistance;
		}
	}
	
	/*
	cout << "The distance matrix is: "<< endl;
	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			cout<<distanceMatrix[i][j] << "   " ;
		}
		cout<<endl;
	}
	*/
	
}


void InputDataProcessor::Rsavesinstance(string fname){
		ofstream myfile;
		myfile.open (fname);
		myfile << "tmax = " << getMaximumTourLength()<< endl;
		myfile << "lmax = " << getMaximumLoadCapacity() << endl;
		myfile << "instance1 = rbind(" << endl;
		int n=basicData.size();
		for(int i = 0; i < n-1; i++){
			myfile << "c(";
			myfile << basicData[i].x << ", " << basicData[i].y<< ", " << basicData[i].quantity<< ", "  << basicData[i].profit ;
			myfile << "), " << endl;
		}
		myfile << "c(" << basicData[n-1].x << ", " << basicData[n-1].y<< ", " << basicData[n-1].quantity << ", "  << basicData[n-1].profit  << endl;
		myfile << ") ) " << endl;
		
		myfile << "tcost = rbind(" << endl;
		for(int i = 0; i < n-1; i++){
			myfile << "c(";
			for(int j = 0; j < n-1; j++){
				myfile << distanceMatrix[i][j] << ", ";
			}
			myfile << distanceMatrix[i][n-1] << " ), " << endl;
		}
		myfile << "c(";
		for(int j = 0; j < n-1; j++){
			myfile << distanceMatrix[n-1][j] << ", ";
		}
		myfile << distanceMatrix[n-1][n-1] << " ) " << endl;
		myfile << ") " << endl;

		myfile.close();
};



