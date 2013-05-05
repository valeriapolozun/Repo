#include <iostream>

#include "OrienteeringProblemWithPickupsAndDeliveries.h"
#include "TourDefault.h"
#include "TourNearestNeighbor.h"
#include "TourNearestNeighborDistanceLimit.h"
#include "TourFirstPickupSecondDelivery.h"

using namespace std;

int main(int argc, char **argv) {
	

	//OrienteeringProblemWithPickupsAndDeliveries problemWithSmallInput("kicsiteszt.txt");
	
	//TourDefault td("kicsiteszt.txt");

	//td.runTwoOpt();
	
	//TourDefault td("set_64_1_50_300.txt");
	// TourNearestNeighbor tnn("kicsiteszt.txt");



	TourNearestNeighborDistanceLimit tnndl("set_64_1_50_300.txt");
	
	tnndl.runShakingAndTwoopt();
	


	

	

	
	//TourNearestNeighbor tnn("kicsiteszt.txt");
	//TourFirstPickupSecondDelivery tfpsd("kicsiteszt.txt");
	//problemWithSmallInput.runTourDefault();

	/*problemWithSmallInput.runTwoOpt();
	
	problemWithSmallInput.runNearestNeighborDistanceLimit();

	problemWithSmallInput.runFirstPickupSecondDeliveryPoint();*/

	int k;
	cin >> k;
	return 0;
}
