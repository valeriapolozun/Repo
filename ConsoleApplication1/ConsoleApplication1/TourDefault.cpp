#include "TourDefault.h"
#include <iostream>


using namespace std;


TourDefault::TourDefault(string inputFile): OrienteeringProblemWithPickupsAndDeliveries(inputFile)
{
	calcTourDefault();
	profitsOfAllTheTours();
}

TourDefault::~TourDefault()
{
}


void TourDefault::calcTourDefault(){
	
	cout << "The tour, in which we visite all the nodes in order" << endl;
	vector<Coordinates> basicData(inputDataProcessor.getBasicData());

	vector<int> tourDefault;
	tourDefault.push_back(0);
	for (int i = 2; i < basicData.size(); i++)
	{
		tourDefault.push_back(i);
	}
	tourDefault.push_back(1);

	for (int i = 0; i < tourDefault.size(); i++)
	{
		cout  << i+1 << ". place in the tour " << tourDefault[i] << endl;
	}

	solutionTours.push_back(tourDefault);	
}

