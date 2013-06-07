#include "ProfitCalculator.h"
#include "InputDataProcessor.h"
#include <glpk.h>
#include <iostream>
#include <fstream>

using namespace std;

ProfitCalculator:: ProfitCalculator(std::vector <int> tourInput,std::vector <Coordinates> basicDataInput, double maxCapacityInput) 
{
	upperbound = - 1;
	CPUtime = 0;
	Iteratations = 0;
	CPUtimetotal = 0;

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
		vectorMaxQuantities.push_back(basicData[i].quantity);
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
	glp_term_out(GLP_OFF);

	const int N = vectorMaxQuantities.size();
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
		ar[k] = -vectorMaxQuantities[i-2];

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

vector <int> ProfitCalculator::getZeroIntensityIndices()
{
	vector <int> zeroIntensityIndices;
	for (int i=1; i < intensity.size()-1; i++)
	{
		if (intensity[i]==0)
		{
			zeroIntensityIndices.push_back(i);
		}
	}
	return zeroIntensityIndices;
}

void ProfitCalculator::loadCalculation(vector <int> & load, vector <int> & goodsOnTheLorry)
{
	load.push_back(0);
	goodsOnTheLorry.push_back(0);
	for (int i=1; i < intensity.size(); i++)
	{
		 load.push_back(intensity[i]*basicData[i].quantity);
		 goodsOnTheLorry.push_back( goodsOnTheLorry[i-1]+ load[i]);
	}
}


void ProfitCalculator::bufferPlusCalculation(vector <int> goodsOnTheLorry, vector <int> & bufferPlus)
{
	double min;
	for (int i=0; i < goodsOnTheLorry.size(); i++)
	{
		min=maxCapacity;
		for (int j=i; j < goodsOnTheLorry.size(); j++)
		{
			if (min> maxCapacity-goodsOnTheLorry[j])
			{
				min=maxCapacity-goodsOnTheLorry[j];
			}
		}	
		bufferPlus.push_back(min);
		
	}
}

void ProfitCalculator::bufferMinusCalculation ( vector <int> goodsOnTheLorry, vector <int> & bufferMinus)
{
	double max;
	for (int i=0; i < goodsOnTheLorry.size(); i++)
	{
		max=0;
		for (int j=i; j < goodsOnTheLorry.size(); j++)
		{
			if (max < goodsOnTheLorry[j])
			{
				max=goodsOnTheLorry[j];
			}
		}
		
		 bufferMinus.push_back(max);
	}	
}









void ProfitCalculator::savesol(string fname)
{
		ofstream myfile;
		myfile.open (fname);
		myfile << "result = list() " << endl;


		myfile << "result$profit =  " << result << endl;
		myfile << "result$bound = " << upperbound << endl;
		myfile << "result$length = " << 50 << endl;
		myfile << "result$maxload = " << 300 << endl;
		myfile << "result$cputime = " << CPUtime << endl;
		myfile << "result$cputimetotal = " << CPUtimetotal << endl;
		myfile << "result$iterations = " << Iteratations << endl;
		myfile << "result$tour = c(";
		for(int i = 0; i < tour.size()-1; i++){
			myfile << 1 + tour[i] << ", ";
		}
		myfile << 1 + tour[tour.size()-1] << ")" << endl;
		myfile << "result$intensity = c(";
		for(int i = 0; i < intensity.size()-1; i++){
			myfile << intensity[i] << ", ";
		}
		myfile << intensity[intensity.size()-1] << ")" << endl;
		myfile << "result$quantity = c(";
		for(int i = 0; i < tour.size()-1; i++){
			myfile << intensity[i]* basicData[i].quantity << ", ";
		}
		myfile << intensity[tour.size()-1]*basicData[tour.size()-1].quantity << ")" << endl;
		myfile.close();
	};




