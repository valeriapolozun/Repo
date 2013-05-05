#include <iostream>

#include "OrienteeringProblemWithPickupsAndDeliveries.h"
#include "TourDefault.h"
#include "TourNearestNeighbor.h"
#include "TourNearestNeighborDistanceLimit.h"
#include "TourFirstPickupSecondDelivery.h"

using namespace std;

int main(int argc, char **argv) {
	
	//TourDefault td("kicsitesztdl3.txt"); // Expected tour: {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}

	//TourNearestNeighbor tnd("kicsitesztdl3.txt"); // Expected tour : { 0,2,3,4, 6, 5, 7, 8, 9, 10}

	//TourFirstPickupSecondDelivery tfirstpickupseconddelivery("testFirstPickupSecondDeliverydl1000.txt"); // Expected tour: { 0,2,4,5,6,1}

	//TourNearestNeighborDistanceLimit tnndl3("kicsitesztdl3.txt"); // Expected tours : no tours
	//TourNearestNeighborDistanceLimit tnndl5("kicsitesztdl5.txt"); // Expected tours : {0,2,1} 
	//TourNearestNeighborDistanceLimit tnndl8("kicsitesztdl8.txt"); // Expected tours : {0,2,3,4,1} , {0,5,1}
	
	TourNearestNeighborDistanceLimit tnndlbig("set_64_1_50_300.txt"); // 
	tnndlbig.runShakingAndTwoopt();

	int k;
	cin >> k;
	return 0;
}
