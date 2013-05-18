#include <iostream>

#include "OrienteeringProblemWithPickupsAndDeliveries.h"
#include "TourDefault.h"
#include "TourNearestNeighbor.h"
#include "TourNearestNeighborDistanceLimit.h"
#include "TourFirstPickupSecondDelivery.h"
#include "TourPickupDeliveryPointPairs.h"
#include "TourPickupDeliveryPointPairsParallelBest.h"
#include "TourPickupDeliveryPointPairsParallel.h"

using namespace std;

int main(int argc, char **argv) {
	

	// The tours are generated in different algorithms based on the input files.

	//TourDefault td("smalltest_dl3.txt"); // Expected tour: {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}

	//TourNearestNeighbor tnd("smalltest_dl3.txt"); // Expected tour : { 0,2,3,4, 6, 5, 7, 8, 9, 10}

	//TourFirstPickupSecondDelivery tfirstpickupseconddelivery("testFirstPickupSecondDeliverydl1000.txt"); // Expected tour: { 0,2,4,5,6,1}

	//TourNearestNeighborDistanceLimit tnndl3("smalltest_dl3.txt"); // Expected tours : no tours
	//TourNearestNeighborDistanceLimit tnndl5("smalltest_dl5.txt"); // Expected tours : {0,2,1} 
	
	//TourNearestNeighborDistanceLimit tnndl8("smalltest_dl8.txt"); // Expected tours : {0,2,3,4,1} , {0,5,1}
	//tnndl8.runShakingAndTwoopt();

	//TourNearestNeighborDistanceLimit tnndlbig("set_64_1_50_300kicsi.txt"); // 
	
	//tnndlbig.runShakingAndTwoopt();
	//tnndlbig.runStringExchangesAndTwoopt();




	//TourPickupDeliveryPointPairs pickupdeliverypairs("set_64_1_50_300small.txt");

	//TourPickupDeliveryPointPairsParallel pickupdeliverypairsparallel("set_64_1_50_300small.txt");

	TourPickupDeliveryPointPairsParallelBest pickupdeliverypairsparallelbestofalltours("set_64_1_50_300small.txt");

	int k;
	cin >> k;
	return 0;
}
