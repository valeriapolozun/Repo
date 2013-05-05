#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>

#include <glpk.h>

using namespace std;

double dMax;

template <typename valueType>
bool SameSign(typename valueType x, typename valueType y)
{
    return (x >= 0) ^ (y < 0);
}


struct Coordinates  // record all the information of points
{
	double x,y;//coordinates
	int quantity; //quantity ( + if pick up, - if delivery) 
	int profit;//value
};

class OrienteeringProblemWithPickupsAndDeliveries
{
public:
	OrienteeringProblemWithPickupsAndDeliveries(string inputFile);
	~OrienteeringProblemWithPickupsAndDeliveries();

	void calcTourDefault();
	void calcTourNearestNeighbor();
	void calcTourNearestNeighborDistanceLimit();
	void calcTourChoosePickupAndDeliveryPointPairs();
	//void calcOtherToursNearestNeighborDistanceLimit(vector<int> tour);
	
	bool isThereUnvisitedNodes();
	int getTotalProfit();
	void profitsOfAllTheTours();

private:
	void readData(string inputFile);
	void initDistanceMatrix();
	void calcDistanceMatrix ();
	int getNearestNeighbor(vector<int> unvisitedCities, int startNode);
	void getInitialTour();
	void getInitialTourMaxProfit();
	double intensityCalculation();
	bool isTotalLengthUnderLimit(vector<int> currentTour, int nodeToAdd);
	void calculateTourProfitAndIntensity();
	void initProfitPerDistanceMatrix();
	void calcPickupDeliveryPointPairs();

	int getPickUpDeliveryPointPairs(vector<int> unvisitedCities, int startNode);
	void pickUpPointToChoose(vector<int> & nodes, vector<Coordinates> basicData);
	int getHighestDemandInNeighborhood(vector<int> unvisitedCities, int startNode);
	


	void initVariablesBetweenCalculations();

	double dMax; // maximum tour length
	double lMax; // maximum load capacity
	vector<Coordinates> basicData; 
	vector<Coordinates> basicDataModified;
	vector <vector<Coordinates> > allBasicDataModified;
	vector<double> maxProfits;
	vector<double> intensity;
	vector<double> initialTour;
	vector < vector <double> > distanceMatrix;
	int problemSize; // the number of nodes
	vector < vector <double> > profitPerDistanceMatrix;
	vector<int> unvisitedNodes;
	vector <vector <int> > solutionTours; // stores the tours
	vector <int> totalProfits;
	
};

OrienteeringProblemWithPickupsAndDeliveries::~OrienteeringProblemWithPickupsAndDeliveries()
{
	
}

OrienteeringProblemWithPickupsAndDeliveries::OrienteeringProblemWithPickupsAndDeliveries(string inputFile)
{
	readData(inputFile);
	problemSize=basicData.size();
	initDistanceMatrix();
	calcDistanceMatrix ();
	initProfitPerDistanceMatrix();
	calcPickupDeliveryPointPairs();
	unvisitedNodes.assign (problemSize, 1);
	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;
	
}

void OrienteeringProblemWithPickupsAndDeliveries::readData(string inputFile)
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

void OrienteeringProblemWithPickupsAndDeliveries::pickUpPointToChoose(vector<int> & nodes, vector<Coordinates> basicData)
{
	for (int i=0; i<nodes.size();i++)
	{
		if (basicData[i].quantity<0)
		nodes[i]=0;
	}	
return;
}		


void OrienteeringProblemWithPickupsAndDeliveries::initDistanceMatrix()
{
	vector <double> initToZeroVector(problemSize,0);
	distanceMatrix.assign(problemSize, initToZeroVector);
}

void OrienteeringProblemWithPickupsAndDeliveries::initProfitPerDistanceMatrix()
{
	vector <double> initToZeroVector(problemSize,0);
	profitPerDistanceMatrix.assign(problemSize, initToZeroVector);
}

void OrienteeringProblemWithPickupsAndDeliveries::calcDistanceMatrix()
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
	
	/*cout << "The distance matrix is: "<< endl;
	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			cout<<distanceMatrix[i][j] << "   " ;
		}
		cout<<endl;
	}*/
	
}

void OrienteeringProblemWithPickupsAndDeliveries::calcPickupDeliveryPointPairs()
	{
	double dx;     // xi.cood - xj.cood
	double dy;     // yi.cood - yj.cood
	double travelDistance; // stores the calculated distances
	double profitValueMax; // profit value per unit which can be reached by visiting a certain pair of Pick up and Delivery point

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
				if (SameSign(basicData[i].quantity, basicData[j].quantity))
					profitPerDistanceMatrix[i][j]=-DBL_MAX;
				else
				{
				profitPerDistanceMatrix[i][j]= ((basicData[i].profit*basicData[i].quantity)+(basicData[j].profit*basicData[j].quantity))/travelDistance;
				}
		}
	}

	/*
	cout << "The profitPerDistance matrix is: "<< endl;
	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			cout<<profitPerDistanceMatrix[i][j] << "   " ;
		}
		cout<<endl;
	}
	*/
}




void OrienteeringProblemWithPickupsAndDeliveries::calcTourNearestNeighbor(){
	initVariablesBetweenCalculations();
	vector<int> unvisitedNodes;

	unvisitedNodes.assign (problemSize, 1);
	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;

	int startNode=0;
	vector<int> tour;
	tour.push_back(0);
	for (int i = 0; i < problemSize-2; i++)
	{
		startNode=getNearestNeighbor(unvisitedNodes, startNode);
		tour.push_back(startNode);
		unvisitedNodes[startNode]=0;
	}
	tour.push_back(1);

	unvisitedNodes.assign (problemSize, 1);
	for (int i=0; i < tour.size(); i++)
	{
		unvisitedNodes[tour[i]]=0;
	}

	for (int i = 0; i < tour.size(); i++)
	{
		//cout  << i+1 << ". place in the tour " << tour[i] << endl;
		basicDataModified.push_back(basicData[tour[i]]);
	}
	
	calculateTourProfitAndIntensity();
}

void OrienteeringProblemWithPickupsAndDeliveries::calcTourNearestNeighborDistanceLimit(){
	initVariablesBetweenCalculations();
	vector<int> unvisitedNodesForOneTour;

	unvisitedNodesForOneTour= unvisitedNodes;
	
	int startNode=0;
	vector<int> tour;
	tour.push_back(0);
	for (int i = 0; i < problemSize-2; i++)
	{
		startNode=getNearestNeighbor(unvisitedNodesForOneTour, startNode);
		if (isTotalLengthUnderLimit(tour, startNode)){
		tour.push_back(startNode);
		}
		unvisitedNodesForOneTour[startNode]=0;
	}
	tour.push_back(1);

	for (int i=0; i < tour.size(); i++)
	{
		unvisitedNodes[tour[i]]=0;
	}

	solutionTours.push_back(tour);	
	
	
	if (tour.size()==2)
	{
		unvisitedNodes.assign(problemSize, 0);
	}
	else 
	{
		solutionTours.push_back(tour);

		for (int i = 0; i < tour.size(); i++)
		{
			cout  << i+1 << ". place in the tour " << tour[i] << endl;
			basicDataModified.push_back(basicData[tour[i]]);
		}
	
		allBasicDataModified.push_back(basicDataModified);

		calculateTourProfitAndIntensity();
	}
}


void OrienteeringProblemWithPickupsAndDeliveries::calcTourChoosePickupAndDeliveryPointPairs(){
	initVariablesBetweenCalculations();
	
	vector<int> unvisitedNodesForOneTour;

	unvisitedNodesForOneTour= unvisitedNodes;
	
	int startNode=0;
	vector<int> tour;
	tour.push_back(0);
	
	pickUpPointToChoose(unvisitedNodesForOneTour, basicData);
	
	startNode= getNearestNeighbor(unvisitedNodesForOneTour, startNode);
	tour.push_back(startNode);


	for (int i = 0; i < problemSize-3; i++)
	{
		startNode=getPickUpDeliveryPointPairs(unvisitedNodesForOneTour, startNode);
		if (isTotalLengthUnderLimit(tour, startNode)){
		tour.push_back(startNode);
		}
		unvisitedNodesForOneTour[startNode]=0;
	}
	tour.push_back(1);

	
	for (int i=0; i < tour.size(); i++)
	{
		unvisitedNodes[tour[i]]=0;
	}

	solutionTours.push_back(tour);	
	
	
	for (int i = 0; i < tour.size(); i++)
	{
		cout  << i+1 << ". place in the tour " << tour[i] << endl;
		basicDataModified.push_back(basicData[tour[i]]);
	}
	
	allBasicDataModified.push_back(basicDataModified);

	calculateTourProfitAndIntensity();

}


int OrienteeringProblemWithPickupsAndDeliveries::getTotalProfit()
{
	int result = 0;
	for(int i = 0; i < basicDataModified.size(); i++) {// this: Solution whose total score is just being computed.
		cout << "basicDataModified[i].quantity: " << basicDataModified[i].quantity 
			<<" basicDataModified.profit[i]: " << basicDataModified[i].profit 
			<< "  intensity[i]:" << intensity[i] 
		    << "  basicDataModified[i].profit*intensity[i]:" << basicDataModified[i].profit*intensity[i]<<endl;
		result += basicDataModified[i].quantity*basicDataModified[i].profit*intensity[i];
	}
	totalProfits.push_back(result);
	return totalProfits.back();
}

int OrienteeringProblemWithPickupsAndDeliveries::getNearestNeighbor (vector<int> unvisitedCities, int startNode)
{
    double min = DBL_MAX;
    int nearest = startNode;
    for(int i=0; i<unvisitedCities.size(); i++)
    {
        if(unvisitedCities[i]==0) continue; 
        if (startNode == i) continue; 
        if(distanceMatrix[startNode][i]<min)
        {
            min = distanceMatrix[startNode][i];
            nearest = i;
        }
    }
    return nearest;
}

int OrienteeringProblemWithPickupsAndDeliveries::getHighestDemandInNeighborhood(vector<int> unvisitedCities, int startNode)
{
	double min = DBL_MAX;
    int nearest = startNode;
    for(int i=0; i<unvisitedCities.size(); i++)
    {
        if(unvisitedCities[i]==0) continue; 
        if (startNode == i) continue; 
        if(distanceMatrix[startNode][i]<min)
        {
            min = distanceMatrix[startNode][i];
            nearest = i;
        }
    }
    return nearest;
}

int OrienteeringProblemWithPickupsAndDeliveries::getPickUpDeliveryPointPairs (vector<int> unvisitedCities, int startNode)
{
    double max = -DBL_MAX;
    int nearest = startNode;
    for(int i=0; i<unvisitedCities.size(); i++)
    {
        if(unvisitedCities[i]==0) continue; 
        if (startNode == i) continue; 
        //cout << "the i: " << i << "the profitperdistancematrix[startnode][i]" << profitPerDistanceMatrix[startNode][i] << " a max value: " << max << endl;
		if(profitPerDistanceMatrix[startNode][i]>max)
        {			
            max = profitPerDistanceMatrix[startNode][i];
            nearest = i;
        }
    }
	//cout << "The next nearest point is: " << nearest << endl;
    return nearest;
}


void OrienteeringProblemWithPickupsAndDeliveries::getInitialTour()
{
	for (int i = 0; i < basicDataModified.size(); i++)
	{
		initialTour.push_back(basicDataModified[i].quantity);
	}

	/*for (int i = 0; i < initialTour.size(); i++)
	{
		cout << "Initial tour quantity " << i << ": " << initialTour[i] << endl;
	}*/
	
}

void OrienteeringProblemWithPickupsAndDeliveries::getInitialTourMaxProfit()
{
	for (int i = 0; i < basicDataModified.size(); i++)
	{
		maxProfits.push_back(basicDataModified[i].quantity*basicDataModified[i].profit);
	}


	/*for (int i = 0; i < maxProfits.size(); i++)
	{
		cout << "Initial tour max profits " << i << ": " << initialTourMaxProfits[i] << endl;
	}*/
	
}

double OrienteeringProblemWithPickupsAndDeliveries::intensityCalculation (){
	glp_prob *lp = glp_create_prob();
	glp_term_out(1);
	double Lmax = 300;

	const int N = initialTour.size();
	const int K = N*3 + 1;

	// int *[K], ja[K];
	int* ia = new int[K+1]; 
	int* ja = new int[K+1]; 
	double* ar= new double[K+1];
	double z, x1, x2;
	
	lp = glp_create_prob();
	glp_set_prob_name(lp, "MFPsub");
	glp_set_obj_dir(lp, GLP_MAX);
	glp_add_rows(lp, N+1);
	
	
	for(int i = 0; i <= N; ++i){
		char strbuffer[50];
		sprintf_s (strbuffer, "LB%d",i);
		glp_set_row_name(lp, i+1, strbuffer);
		glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);
	}

	glp_add_cols(lp, N+1+N);

	for(int i = 0; i < N+1; ++i){
		char strbuffer[50];
		sprintf_s (strbuffer, "L%d",i);
		glp_set_col_name(lp, i+1, strbuffer);
		glp_set_col_bnds(lp, i+1, GLP_DB, 0.0, Lmax);
		glp_set_obj_coef(lp, i+1, 0.0);
	}
	
	for(int i = 0; i < N; ++i){
		char strbuffer[50];
		sprintf_s (strbuffer, "y%d",i+1);
		glp_set_col_name(lp, N+1+i+1, strbuffer);
		
		
		glp_set_col_bnds(lp, N+1+i+1, GLP_DB, 0.0, 1);
		
		glp_set_obj_coef(lp, N+1+i+1, maxProfits[i]);
	}
	

	int k = 1;
	ia[k] = 1;
	ja[k] = 1;
	ar[k] = 1;
	++k;
	for(int i = 2; i <= N+1; ++i){ 
		ia[k] = i;
		ja[k] = i-1;
		ar[k] = -1;
		++k;
		ia[k] = i;
		ja[k] = i;
		ar[k] = 1;
		++k;
		ia[k] = i;
		ja[k] = N+i;
		ar[k] = -initialTour[i-2];

		++k;
	}

	/*for (int i = 0; i <=  K;  ++i) {
		cout  << i << "\t" << ia[i]<< "\t"  << ja[i] << "\t" << ar[i] << endl;
	}*/

	glp_load_matrix(lp, K, ia, ja, ar);
	 glp_write_lp(lp, NULL,"test.lp");
	glp_simplex(lp, NULL );

    for (int i = 0; i <  N;  ++i) {
		cout << " Goods on the lorry in the timeperiod " << i << ": " << glp_get_col_prim(lp, i+1) << endl;
	}

	
	 for (int i = N+1; i <  2*N+1;  ++i) {
		 intensity.push_back(glp_get_col_prim(lp, i+1));
	}


	delete ia; 
	delete ja; 
	delete ar;

	return 0;
}

void OrienteeringProblemWithPickupsAndDeliveries::calcTourDefault(){
	
	vector<Coordinates> basicData(this->basicData);
	initVariablesBetweenCalculations();
	Coordinates lastNode=basicData[1];
	basicData.erase(basicData.begin()+1);
	basicData.push_back(lastNode);

	vector<int> tourDefault;
	for (int i = 0; i < basicData.size(); i++)
	{
		tourDefault.push_back(i);
	}
	
	for (int i = 0; i < tourDefault.size(); i++)
	{
		basicDataModified.push_back(basicData[tourDefault[i]]);
	}
	
	calculateTourProfitAndIntensity();
}

bool OrienteeringProblemWithPickupsAndDeliveries::isTotalLengthUnderLimit(vector<int> currentTour, int nodeToAdd){
	double result = 0;
	int lastNode = currentTour[0];
	currentTour.push_back(nodeToAdd);
	
	for(int i = 1; i < currentTour.size(); i++) { // The total length of the tour will be computed
		int currentNode = currentTour[i];
		result += distanceMatrix[lastNode][currentNode];
		lastNode = currentNode;
	}
	result= result+ distanceMatrix[lastNode][1];
	return result<=dMax; // Is my tour length smaller than dMax. TRUE or FALSE will return

}

void OrienteeringProblemWithPickupsAndDeliveries::calculateTourProfitAndIntensity()
{
	getInitialTour();
	getInitialTourMaxProfit();
	intensityCalculation();

}

void OrienteeringProblemWithPickupsAndDeliveries::initVariablesBetweenCalculations()
{
	basicDataModified.clear();
	maxProfits.clear();
	intensity.clear();
	initialTour.clear();
}

bool OrienteeringProblemWithPickupsAndDeliveries::isThereUnvisitedNodes()
{
	return (find(unvisitedNodes.begin(), unvisitedNodes.end(), 1) != unvisitedNodes.end());
}

void OrienteeringProblemWithPickupsAndDeliveries::profitsOfAllTheTours()
{
	for (int i = 0; i < totalProfits.size(); i++)
	{
		cout << "The " << i+1 << ". tour has the profit of " << totalProfits[i] << endl;
	}
}


int main(int argc, char **argv) {

	
	OrienteeringProblemWithPickupsAndDeliveries problemWithSmallInput("set_64_1_50_300kicsi.txt");

	// Tours will be generated until all the points are visited. Each tour has the maximum distance limit.
	while (problemWithSmallInput.isThereUnvisitedNodes())
	{
	problemWithSmallInput.calcTourNearestNeighborDistanceLimit();
	problemWithSmallInput.getTotalProfit();
	} 

	problemWithSmallInput.profitsOfAllTheTours();


	//problemWithSmallInput.initProfitPerDistanceMatrix();
	//problemWithSmallInput.calcPickupDeliveryPointPairs();
	/*
	do  // Tours will be generated until all the points are visited. Each tour has the maximum distance limit.
	{
	problemWithSmallInput.calcTourChoosePickupAndDeliveryPointPairs();
	} while (problemWithSmallInput.isThereUnvisitedNodes());

	problemWithSmallInput.profitsOfAllTheTours();
	*/
//	problemWithSmallInput.calcTourChoosePickupAndDeliveryPointPairs();
	//cout << "Tour with Pick up and Delivery Point Pairs max profit: " << problemWithSmallInput.getTotalProfit() << endl;
	

	int k;
	cin >> k;
	return 0;
}
