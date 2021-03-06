#include <iostream>

#include "OrienteeringProblemWithPickupsAndDeliveries.h"
#include "TourDefault.h"
#include "TourNearestNeighbor.h"
#include "TourNearestNeighborDistanceLimit.h"
#include "TourFirstPickupSecondDelivery.h"
#include "TourPickupDeliveryPointPairs.h"
#include "TourPickupDeliveryPointPairsParallelBest.h"
#include "TourPickupDeliveryPointPairsParallel.h"
#include "TourGreedyEinzelnePunkteSeriell.h"
#include "ExcelExporter.h"

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



	//TourNearestNeighborDistanceLimit tnndlbig("smalltest_dl8test2opt.txt"); // 
	
	//tnndlbig.runShakingAndTwoopt();
	/*
	std::vector <int> probatour;
	probatour.push_back(0);
	probatour.push_back(2);
	probatour.push_back(3);
	probatour.push_back(4);
	probatour.push_back(5);
	probatour.push_back(1);

	tnndlbig.doTwoOpt(probatour);
	*/


	
	//tnndlbig.runStringExchangesAndTwoopt();

	//TourFirstPickupSecondDelivery tfirstpickupseconddelivery("set_64_1_50_300.txt");
	
	

	//TourPickupDeliveryPointPairsParallel pickupdeliverypairsparallel("set_64_1_50_300.txt");

	//TourPickupDeliveryPointPairsParallelBest pickupdeliverypairsparallelbestofalltours("set_64_1_50_300.txt");
	//pickupdeliverypairsparallelbestofalltours.runTwoOpt();

	/*
	for (int i=0; i<10; i++)
	{
	pickupdeliverypairsparallelbestofalltours.erasePoints(1, 0.1);
	pickupdeliverypairsparallelbestofalltours.doInsertion();
	pickupdeliverypairsparallelbestofalltours.runTwoOpt();
	}
	*/
	
	int start=1;
	for (int i=10; i<86; i+=5)
	{
	OrienteeringProblemWithPickupsAndDeliveries proba("tsiligirides_problem_1_budget_" + std::to_string(i) + "_300.txt");
	//OrienteeringProblemWithPickupsAndDeliveries proba ("tsiligirides_problem_1_budget_30_300.txt");
	if (start==1)
	{
	proba.runExcelExportStart();
	start=start+1;
	}
	
		TourGreedyEinzelnePunkteSeriell tourgreedyeinzeln("tsiligirides_problem_1_budget_" + std::to_string(i) + "_300.txt", 1);
		
		TourGreedyEinzelnePunkteSeriell tourgreedyeinzeln2("tsiligirides_problem_1_budget_" + std::to_string(i) + "_300.txt", 3);
		
		
		TourGreedyEinzelnePunkteSeriell tourgreedyeinzeln3("tsiligirides_problem_1_budget_" + std::to_string(i) + "_300.txt", 15);
		
		TourGreedyEinzelnePunkteSeriell tourgreedyeinzeln4("tsiligirides_problem_1_budget_" + std::to_string(i) + "_300.txt", 20);
		
		TourPickupDeliveryPointPairs pickupdeliverypairs("tsiligirides_problem_1_budget_" + std::to_string(i) + "_300.txt", 1);
		
		TourPickupDeliveryPointPairs pickupdeliverypairs2("tsiligirides_problem_1_budget_" + std::to_string(i) + "_300.txt", 3);
	
		TourPickupDeliveryPointPairs pickupdeliverypairs3("tsiligirides_problem_1_budget_" + std::to_string(i) + "_300.txt", 15);
		TourPickupDeliveryPointPairs pickupdeliverypairs4("tsiligirides_problem_1_budget_" + std::to_string(i) + "_300.txt", 20);
		
		TourPickupDeliveryPointPairsParallel pickupdeliverypairsparallel("tsiligirides_problem_1_budget_" + std::to_string(i) + "_300.txt", 1);
		TourPickupDeliveryPointPairsParallel pickupdeliverypairsparallel2("tsiligirides_problem_1_budget_" + std::to_string(i) + "_300.txt", 3);
	
		TourPickupDeliveryPointPairsParallel pickupdeliverypairsparallel3("tsiligirides_problem_1_budget_" + std::to_string(i) + "_300.txt", 15);
		
		TourPickupDeliveryPointPairsParallel pickupdeliverypairsparallel4("tsiligirides_problem_1_budget_" + std::to_string(i) + "_300.txt", 20);
		
	}
	
	//proba.runExcelExportFinish();

	//tourgreedyeinzeln.runTwoOpt();
	//tourgreedyeinzeln.runExcelExport();

	/*
	tourgreedyeinzeln.erasePoints(1, 0.1);
	tourgreedyeinzeln.erasePoints(1, 0.1);
	tourgreedyeinzeln.erasePoints(1, 0.1);


	tourgreedyeinzeln.doInsertion();
	*/


	cout << "Finished! Press something to exit";
	int m;
	cin >> m;
	return 0;
}
