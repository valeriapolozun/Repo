#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>

#include <glpk.h>

using namespace std;

double dMax;

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

	int getTotalProfit();

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

	void initVariablesBetweenCalculations();

	double dMax;
	double lMax;
	vector<Coordinates> basicData;
	vector<Coordinates> basicDataModified;
	vector<double> maxProfits;
	vector<double> intensity;
	vector<double> initialTour;
	double** distanceMatrix;
	int problemSize;
};

OrienteeringProblemWithPickupsAndDeliveries::~OrienteeringProblemWithPickupsAndDeliveries()
{
	for( int i = 0 ; i < problemSize ; i++ )
	{
		delete [] distanceMatrix[i] ;
	}

	delete [] distanceMatrix ;
}

OrienteeringProblemWithPickupsAndDeliveries::OrienteeringProblemWithPickupsAndDeliveries(string inputFile)
{
	readData(inputFile);
	problemSize=basicData.size();
	initDistanceMatrix();
	calcDistanceMatrix ();
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
}

void OrienteeringProblemWithPickupsAndDeliveries::initDistanceMatrix()
{
	distanceMatrix = new double*[problemSize];
	for(int i =0; i < problemSize; i++)
	{
		distanceMatrix[i] = new double[problemSize];
	}
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



	for (int i = 0; i < tour.size(); i++)
	{
		//cout  << i+1 << ". place in the tour " << tour[i] << endl;
		basicDataModified.push_back(basicData[tour[i]]);
	}
	
	calculateTourProfitAndIntensity();
}

void OrienteeringProblemWithPickupsAndDeliveries::calcTourNearestNeighborDistanceLimit(){
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
		if (isTotalLengthUnderLimit(tour, startNode)){
		tour.push_back(startNode);
		}
		unvisitedNodes[startNode]=0;
	}
	tour.push_back(1);



	for (int i = 0; i < tour.size(); i++)
	{
		//cout  << i+1 << ". place in the tour " << tour[i] << endl;
		basicDataModified.push_back(basicData[tour[i]]);
	}
	
	calculateTourProfitAndIntensity();

}

int OrienteeringProblemWithPickupsAndDeliveries::getTotalProfit()
{
	int result = 0;
	for(int i = 0; i < basicDataModified.size(); i++) {// this: Solution whose total score is just being computed.
		/*cout << "basicDataModified[i].quantity: " << basicDataModified[i].quantity 
			<<" basicDataModified.profit[i]: " << basicDataModified[i].profit 
			<< "  intensity[i]:" << intensity[i] 
		    << "  basicDataModified[i].profit*intensity[i]:" << basicDataModified[i].profit*intensity[i]<<endl;*/
		result += basicDataModified[i].quantity*basicDataModified[i].profit*intensity[i];
	}
	return result;
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

void OrienteeringProblemWithPickupsAndDeliveries::getInitialTour()
{
	for (int i = 0; i < basicDataModified.size(); i++)
	{
		initialTour.push_back(basicDataModified[i].quantity);
	}

	for (int i = 0; i < initialTour.size(); i++)
	{
		cout << "Initial tour quantity " << i << ": " << initialTour[i] << endl;
	}
	
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

	for (int i = 0; i <=  K;  ++i) {
		cout  << i << "\t" << ia[i]<< "\t"  << ja[i] << "\t" << ar[i] << endl;
	}

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

double intensityCalculation (vector<double> balancetour,vector<double> profittour, vector<double>& intensity){
	glp_prob *lp = glp_create_prob();
	glp_term_out(1);
	double Lmax = 300;

	const int N = balancetour.size();
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
		
		glp_set_obj_coef(lp, N+1+i+1, profittour[i]);
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
		ar[k] = -balancetour[i-2];

		++k;
	}

	for (int i = 0; i <=  K;  ++i) {
		cout  << i << "\t" << ia[i]<< "\t"  << ja[i] << "\t" << ar[i] << endl;
	}

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



struct Solutions 
{
	vector<int> tour;
	//int getTotalProfit(const vector<Coordinates> &basicData, const vector<double>& intensity); // record the total score
	//double getTotalDistance(double** distanceMatrix); // record the distance
};


double getTotalProfit(const vector<Coordinates> &basicData, const vector<double>& intensity) {
	int result = 0;
	for(int i = 0; i < basicData.size(); i++) {// this: Solution whose total score is just being computed.
		cout << "basicData[i].quantity: " << basicData[i].quantity 
			<<" basicData.profit[i]: " << basicData[i].profit 
			<< "  intensity[i]:" << intensity[i] 
		    << "  basicData[i].profit*intensity[i]:" << basicData[i].profit*intensity[i]<<endl;
		result += basicData[i].quantity*basicData[i].profit*intensity[i];
	}
	return result;
}


bool isTotalLengthUnderLimit(double dMax, vector<int> currentTour, int nodeToAdd, double** distanceMatrix){
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


void readData (vector<Coordinates> & basicData,double &LMax, double &dMax)
{
	Coordinates c;
	char line [80];

	ifstream file ("set_64_1_50_300kicsi.txt"); //open file

	if (!file)
	{ 
		cout << "Error opening file!"<<endl;
	}

	file >> dMax; // read in dMax
	file >> LMax; // read in load capacity
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
}


void calcDistanceMatrix (vector<Coordinates> basicData,double** distanceMatrix)
{
	double dx;     // xi.cood - xj.cood
	double dy;     // yi.cood - yj.cood
	double travelDistance;  // stores the calculated time
	int problemSize=basicData.size();
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
	cout << "The distance matrix is: "<< endl;
	for (int i = 0; i < problemSize; i++)
	{
		for(int j = 0; j < problemSize; j++)
		{
			cout<<distanceMatrix[i][j] << "   " ;
		}
		cout<<endl;
	}
	
}

void getInitialTour(vector<double>& initialTourQuantities, const vector<Coordinates>& basicData)
{
	for (int i = 0; i < basicData.size(); i++)
	{
		initialTourQuantities.push_back(basicData[i].quantity);
	}

	for (int i = 0; i < initialTourQuantities.size(); i++)
	{
		cout << "Initial tour quantity " << i << ": " << initialTourQuantities[i] << endl;
	}
	
}


void getInitialTourMaxProfit(vector<double>& initialTourMaxProfits, const vector<Coordinates>& basicData)
{
	for (int i = 0; i < basicData.size(); i++)
	{
		initialTourMaxProfits.push_back(basicData[i].quantity*basicData[i].profit);
	}


	for (int i = 0; i < initialTourMaxProfits.size(); i++)
	{
		cout << "Initial tour max profits " << i << ": " << initialTourMaxProfits[i] << endl;
	}
	
}

void calcTourDefault(vector<Coordinates> basicData, vector<Coordinates>& basicDataModified){
	
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
	
}

int getNearestNeighbor (double** distanceMatrix,vector<int> unvisitedCities, int startNode)
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


void calcTourNearestNeighbor(double** distanceMatrix,vector<Coordinates> basicData, vector<Coordinates>& basicDataModified){
	
	vector<int> unvisitedNodes;
	int problemSize=basicData.size();
	unvisitedNodes.assign (problemSize, 1);
	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;

	int startNode=0;
	vector<int> tour;
	tour.push_back(0);
	for (int i = 0; i < problemSize-2; i++)
	{
		startNode=getNearestNeighbor(distanceMatrix, unvisitedNodes, startNode);
		tour.push_back(startNode);
		unvisitedNodes[startNode]=0;
	}
	tour.push_back(1);



	for (int i = 0; i < tour.size(); i++)
	{
		cout  << i+1 << ". place in the tour " << tour[i] << endl;
		basicDataModified.push_back(basicData[tour[i]]);
	}
	
}

void calcTourNearestNeighborDistanceLimit(double** distanceMatrix,vector<Coordinates> basicData, vector<Coordinates>& basicDataModified){
	
	vector<int> unvisitedNodes;
	int problemSize=basicData.size();
	unvisitedNodes.assign (problemSize, 1);
	unvisitedNodes[0]=0;
	unvisitedNodes[1]=0;

	int startNode=0;
	vector<int> tour;
	tour.push_back(0);
	for (int i = 0; i < problemSize-2; i++)
	{
		startNode=getNearestNeighbor(distanceMatrix, unvisitedNodes, startNode);
		if (isTotalLengthUnderLimit(dMax, tour, startNode, distanceMatrix)){
		tour.push_back(startNode);
		}
		unvisitedNodes[startNode]=0;
	}
	tour.push_back(1);



	for (int i = 0; i < tour.size(); i++)
	{
		cout  << i+1 << ". place in the tour " << tour[i] << endl;
		basicDataModified.push_back(basicData[tour[i]]);
	}
	
}

	

int main(int argc, char **argv) {

	OrienteeringProblemWithPickupsAndDeliveries problemWithSmallInput("set_64_1_50_300kicsi.txt");
	problemWithSmallInput.calcTourDefault();
	cout << problemWithSmallInput.getTotalProfit() << endl;
	
	problemWithSmallInput.calcTourNearestNeighbor();
	cout << problemWithSmallInput.getTotalProfit() << endl;
	
	problemWithSmallInput.calcTourNearestNeighborDistanceLimit();
	cout << problemWithSmallInput.getTotalProfit() << endl;
	
	problemWithSmallInput.calcTourDefault();
	cout << problemWithSmallInput.getTotalProfit() << endl;


	/*double LMax; // load max
	vector<Coordinates> basicData;

	readData (basicData,LMax, dMax);
	
	for(int i=0; i<basicData.size();i++)
	{
		cout<<basicData[i].x<<"  "<<basicData[i].y<<"  "<<basicData[i].quantity<< basicData[i].profit<<endl;
	}

	const int problemSize=basicData.size();
	double** distanceMatrix = new double*[problemSize];
	for(int i =0; i < problemSize; i++)
	{
		distanceMatrix[i] = new double[problemSize];
	}
	//double distanceMatrix[problemSize][problemSize];

	calcDistanceMatrix (basicData,distanceMatrix);


	vector<Coordinates> basicDataModified;

	// calcTourDefault(basicData, basicDataModified);

	//calcTourNearestNeighbor(distanceMatrix, basicData, basicDataModified);

	calcTourNearestNeighborDistanceLimit(distanceMatrix, basicData, basicDataModified);

	vector<double> initialTour;
	getInitialTour(initialTour, basicDataModified);

	vector<double> maxProfits;
	getInitialTourMaxProfit(maxProfits, basicDataModified);
	


	vector<double> intensity;
	intensityCalculation(initialTour, maxProfits, intensity);
	for(int i=0; i<intensity.size();i++)
	{
		cout<< " the intensity " << i+1 << ":  " << intensity [i]<<endl;
	}
    
	cout << getTotalProfit(basicDataModified, intensity)<< endl;
	

	for( int i = 0 ; i < problemSize ; i++ )
	{
		delete [] distanceMatrix[i] ;
	}

	delete [] distanceMatrix ;
	*/
	int k;
	cin >> k;
	return 0;
}
