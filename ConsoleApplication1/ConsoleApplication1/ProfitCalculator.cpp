#include "ProfitCalculator.h"
#include <glpk.h>
#include <iostream>

using namespace std;

ProfitCalculator:: ProfitCalculator(std::vector <int> tourInput,std::vector <Coordinates> basicDataInput, double maxCapacityInput) 
{
	tour=tourInput;
	maxCapacity= maxCapacityInput;
	for (int i = 0; i < tour.size(); i++)
	{
		basicData.push_back(basicDataInput[tour[i]]);
	}
	result=0;
	calculateProfit();
}

void ProfitCalculator:: calculateProfit()
{
	
	for (int i = 0; i < basicData.size(); i++)
	{
		initialTour.push_back(basicData[i].quantity);
	}


	for (int i = 0; i < basicData.size(); i++)
	{
		maxProfits.push_back(basicData[i].quantity*basicData[i].profit);
	}

	intensityCalculation();

	
	if (basicData.size()>2)
	{
		for(int i = 0; i < basicData.size(); i++) 
		{
			result += basicData[i].quantity*basicData[i].profit*intensity[i];
		}
	}


}


double ProfitCalculator::intensityCalculation (){
	glp_prob *lp = glp_create_prob();
	glp_term_out(1);

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
		glp_set_col_bnds(lp, i+1, GLP_DB, 0.0, maxCapacity); // Lmax
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
	// glp_write_lp(lp, NULL,"test.lp");
	glp_simplex(lp, NULL );

	/*
    for (int i = 0; i <  N;  ++i) {
		cout << " Goods on the lorry in the timeperiod " << i << ": " << glp_get_col_prim(lp, i+1) << endl;
	}
	*/
	
	 for (int i = N+1; i <  2*N+1;  ++i) {
		 intensity.push_back(glp_get_col_prim(lp, i+1));
	}


	delete ia; 
	delete ja; 
	delete ar;

	return 0;
}

double ProfitCalculator::getProfit()
{
	return result;
}