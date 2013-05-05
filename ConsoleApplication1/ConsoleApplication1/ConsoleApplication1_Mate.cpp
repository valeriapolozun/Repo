// -------------------------------------------------------------- -*- C++ -*-
// File: PADOP
// Martin Romauch
// Version 1.0  
// --------------------------------------------------------------------------

//#undef DEBUG
//#define NDEBUG

//#define CPLEXAVAILABLE
#define GLPKAVAILABLE
#undef CPLEXAVAILABLE
//#undef GLPKAVAILABLE

#ifdef CPLEXAVAILABLE
#include <ilcplex/ilocplex.h>
#endif

#ifndef CPLEXAVAILABLE
#include <ctime>
#endif


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <functional>
#include <iterator>
#include <numeric>
#include <assert.h>     /* assert */

#ifdef GLPKAVAILABLE
#include <glpk.h>
#endif

using namespace std;


#ifdef CPLEXAVAILABLE
ILOSTLBEGIN
typedef IloArray<IloNumArray>    FloatMatrix;
typedef IloArray<IloNumVarArray> NumVarMatrix;
#endif

typedef double myfloattype;
typedef int myinttype;

#define TOL 1e-10 

// template fetch subset of indices from a vector 
template <class T> void select(vector<myinttype> indices, vector<T> values, vector<T> & result) {
  int m = indices.size();
  for(int i = 0; i<m; i++){
	  result.push_back(values[indices[i]]);
  }
}

template <class T> T GetMax (T a, T b) {
  T result;
  result = (a>b)? a : b;
  return (result);
}

// TOPD will contain the instance and methods for 
class TOPD;

struct point{
	myfloattype x;
	myfloattype y;
};

class tabu_list{
public:
	vector<myinttype> list;
	TOPD *topd;
};

// solution 
class sol {
public:
	vector<myinttype> tour;
	vector<myfloattype> intensity;
	vector<myfloattype> quantity;
	vector<myfloattype> load;
	vector<myfloattype> bufferplus; // max quantity that is possible to add
	vector<myfloattype> bufferminus; // max quantity that is possible to remove
	myfloattype length;
	myfloattype maxload;
	myfloattype profit;
	myfloattype upperbound;
	TOPD *topd;
	sol();
	void savesol(string fname);
};

sol::sol(){
	upperbound = - 1;
}

void sol::savesol(string fname){
		ofstream myfile;
		myfile.open (fname);
		myfile << "result = list() " << endl;
		myfile << "result$profit =  " << profit << endl;
		myfile << "result$bound = " << upperbound << endl;
		myfile << "result$length = " << length << endl;
		myfile << "result$maxload = " << maxload << endl;
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
		for(int i = 0; i < quantity.size()-1; i++){
			myfile << quantity[i] << ", ";
		}
		myfile << quantity[quantity.size()-1] << ")" << endl;
		myfile.close();
	};

struct infomove{
	int status;				// =1 if o.k. 
	int type;				// =1 if pickup, =-1 if delivery
	myinttype node;			// node
	myfloattype heurval;	// heur information
	myinttype nodepos;		// original position in vector
	myinttype tourpos;		// position in tour
};


// used for initializing TOPD
struct TOPDparameters {
	TOPDparameters();
	myfloattype TourCostFactor; // objective: influence of transport costs 
	myfloattype  ILOMaxCompTime; // CPX limit for comp. time 
	string setILOParHeur; // for CPX "heu"
	string path;  //= "C:/Users/Martin/Desktop/projects/pickupanddeliveryteamorienteering/orienteering instances/dat2/";
	string rpath; // = "C:/Users/martin/Desktop/projects/pickupanddeliveryteamorienteering/R/";
	string rexe; //= "\"C:/Program Files/R/R-2.15.0/bin/R\"";
	string filename; //= "tsiligirides_problem_2_budget_32_300.txt";
    string filename2;
	string CPXpar; //= "heu";
	string type; // = "CPLEX";
};

TOPDparameters::TOPDparameters() { 
	TourCostFactor = 1; // objective: influence of transport costs 
	ILOMaxCompTime = 10; // CPX limit for comp. time 
	setILOParHeur = "heu"; // for CPX "heu"
	path  = "C:/Users/Martin/Desktop/projects/pickupanddeliveryteamorienteering/orienteering instances/dat2/";
	rpath  = "C:/Users/martin/Desktop/projects/pickupanddeliveryteamorienteering/R/";
	rexe = "\"C:/Program Files/R/R-2.15.0/bin/R\"";
	filename = "tsiligirides_problem_2_budget_32_300.txt";
    filename2;
	CPXpar = "heu";
	type = "CPLEX";
}

// heurinfo is a struct which conatins heurinfo1 heurinfo2 
struct heurinfo2 {
	set <myinttype> index; 
	vector <myfloattype> value;
	multimap<myfloattype, myinttype,std::greater<myfloattype> >  mm;
	myfloattype sum;
};

struct heurinfo1 {
	heurinfo2 selected;
	heurinfo2 available;
};

struct heurinfo {
	heurinfo1 pickup;
	heurinfo1 delivery;// represents heursol
	vector<myinttype> type; // -2,-1,0,1,2 ... delivery.sel,delivery.av,start&end,pickup.av,pickup.sel
	vector<myfloattype> heurval; // contains the heurvalues for the node
};

// vector of elitist solutions
struct ElitLst {
	vector <sol> member;
	myfloattype lb; //lower bound
	myfloattype ub; //upper bound
};

struct detourretvalue {myinttype index; myfloattype val;};

class TOPD {
public:
	// instance
	void setfilename(string fname1,string fname2 , string path1);
	void load();
	void displayfile();
	void setTourCostFactor(myfloattype number);
	
	// general
	
	myfloattype calcprofit();
	myfloattype calcprofit(sol & solution);
	myfloattype calclength(vector< myinttype> tour);

	#ifdef  CPLEXAVAILABLE
	// algortihm CPLEX
	IloEnv env;
	IloModel model; 
	IloCplex cplex;
	void convert2ILO();
	void ILOcreateandsolve();
	sol ILOsol;
	void setComptime(myfloattype number);
	void setILOParHeur(string par);
	~TOPD();
	myfloattype calclength();
    #endif

	
	TOPD(TOPDparameters parameters);

	// R output
	void Rsavesinstance(string fname);
	void Rprintsol(string rexe, string rpath, string filename, string extension,sol solution);
	void setRpar(string str1,string str2);

	// heuristic 
	sol heursol; 
	ElitLst elitists;
	void initheursol(); 
	void initsol(sol & solution);
	void initializheuristicinfo(int type = 1);
	void updateheuristicinfo(int type=1, myinttype par1=0);
	int repair(sol & solution);
	myfloattype optimizequantities(sol & solution);
	myfloattype optimizequantities2(sol & solution);

	// 
	infomove addpickup(sol & solution );
	infomove adddelivery(sol & solution );
	int checkcalcbufferplus();
	int checksolution();
	int checksolution(sol & solution);

	int insertpickup(sol & solution,detourretvalue ret, myinttype newpickup, myfloattype intensity = 1);
	int insertdelivery(sol & solution,detourretvalue ret, myinttype newpickup, myfloattype intensity = 1);
	int deletenode(sol & solution, myinttype index);
	int deletenode(infomove retval);
	int deletenode(myinttype tourpos);
	int destroy(sol & solution, myinttype type = 1);

	int emptysol(); // start-end
	int greedyconstruct(sol & solution ); // uses heuristicinfo
	int randconstruct(sol & solution ); // uses heuristicinfo as prob
	void pertubate(myfloattype epsilon); // pertubates heuristicinfo

private:
	string filename;
	string filename2;
	string path;

	string rexe;
	string rpath;

	myinttype n;
	myinttype start; 
	myinttype end;
	myfloattype Tmax;
	myfloattype Lmax;

	myfloattype TourCostFactor;

	vector< point > location;
	vector< myfloattype > profit;
	vector< myfloattype > balance;
	vector< vector<myfloattype> >  distancematrix;
	
	vector< myinttype > pickup;
	vector< myinttype > delivery;
	
	myfloattype calcdist(point a, point b);
	void calcdistmat();
	myfloattype calcmaxload();
	myfloattype calcmaxload(vector< myinttype> tour, vector<myfloattype> intensity);	

	#ifdef  CPLEXAVAILABLE
	// CPLEX
	IloInt	nbRetailer;
	IloNum  ILOTourCostFactor;
	IloNum  ILOParHeur;
	IloNum  ILOMaxCompTime;
	IloNumArray ILOprofit;
	IloNumArray ILObalance;
    IloArray<IloNumArray> ILOdistancematrix;
	#endif
	
	// HEURISTIC
	myfloattype L_toleranceminus;
	myfloattype L_toleranceplus;
	myfloattype T_tolerance;
	heurinfo heuristicinfo; // e.g. profit // multimap structure
	detourretvalue detourlength(myinttype cand);
	detourretvalue detourlength(myinttype cand, int type);
	detourretvalue detourlengthold(myinttype cand, int type);
	int undo_addpickup(infomove retval);
	int undo_adddelivery(infomove retval);
	void updatesol();
	void updatesol(sol & solution);
	friend class TWOOPT;
	friend class TWOEXCHANGE;
};

int op_multt(int i) { return ++i; }

void TOPD::initializheuristicinfo(int type){
	if(type == 1){
		heuristicinfo.type = vector<myinttype>(n,0);
		heuristicinfo.heurval = vector<myfloattype>(n,0);
		heuristicinfo.pickup.available.sum = 0;
		heuristicinfo.delivery.available.sum = 0;
		heuristicinfo.pickup.selected.sum = 0;
		heuristicinfo.delivery.selected.sum = 0;
		for(int i = 0 ; i < profit.size(); i++){	
			if(balance[i]>0){
				heuristicinfo.pickup.available.value.push_back(-profit[i]);
				heuristicinfo.pickup.available.index.insert(i);
				heuristicinfo.type[i] = 1;
				heuristicinfo.heurval[i] = -profit[i];
				heuristicinfo.pickup.available.mm.insert(pair<myfloattype, myinttype>(-profit[i],i));
				heuristicinfo.pickup.available.sum += -profit[i];
			} 
			if(balance[i]<0){
				heuristicinfo.delivery.available.value.push_back(-profit[i]);
				heuristicinfo.delivery.available.index.insert(i);
				heuristicinfo.type[i] = -1;
				heuristicinfo.heurval[i] = -profit[i];
				heuristicinfo.delivery.available.mm.insert(pair<myfloattype, myinttype>(-profit[i],i));
				heuristicinfo.delivery.available.sum += -profit[i];
			}
		}
	}
}


// work in progress
void TOPD::updateheuristicinfo(int type, myinttype par1){
	if(type == 1){
		// myfloattype Lbuffer = Lmax - heursol.maxload;
		// myfloattype Tbuffer = Tmax - heursol.length;
		// myfloattype d2 = numeric_limits<myfloattype>::infinity();

		heuristicinfo.pickup.available.mm.clear();
		
		
		for (set <myinttype>::iterator iit = heuristicinfo.pickup.available.index.begin(); iit!= heuristicinfo.pickup.available.index.end(); ++iit){
		//for(int i = 0 ; i < heuristicinfo.pickup.available.index.size(); ++i){
			pair<myfloattype, myinttype> entry; 
			myinttype pickupnode = *iit;
			entry.second = pickupnode;
			myfloattype d1 = detourlength(pickupnode).val; 
				
				
			set <myinttype>::iterator jit = heuristicinfo.delivery.available.index.begin();
			myinttype deliverynode = *(jit);
				

			myfloattype prftmx = - TourCostFactor * distancematrix[pickupnode][deliverynode] + min(balance[pickupnode],-balance[deliverynode])*(profit[pickupnode]-profit[deliverynode]);
			entry.first = - 2 * TourCostFactor * d1 + prftmx;
			myfloattype prft;
				
			++jit;
			for (; jit!= heuristicinfo.delivery.available.index.end(); ++jit){
			//for(int j = 1; j <heuristicinfo.delivery.available.index.size(); ++j){
				deliverynode = *jit;
				prft = - TourCostFactor * distancematrix[pickupnode][deliverynode] + min(balance[pickupnode],-balance[deliverynode])*(profit[pickupnode]-profit[deliverynode]);
				if(prftmx < prft){
					prftmx = prft;
				}
			}	

			entry.first = - 2 * TourCostFactor * d1 + prftmx;
			heuristicinfo.pickup.available.mm.insert(entry);
			heuristicinfo.heurval[pickupnode] = entry.first;
		}
	}
	if(type == 2){
		const myinttype center = par1; //center of deletion
		// myfloattype Lbuffer = Lmax - heursol.maxload;
		// myfloattype Tbuffer = Tmax - heursol.length;
		// myfloattype d2 = numeric_limits<myfloattype>::infinity();

		heuristicinfo.pickup.selected.mm.clear();
		heuristicinfo.pickup.selected.sum = 0;
		for (set <myinttype>::iterator iit = heuristicinfo.pickup.selected.index.begin(); iit!= heuristicinfo.pickup.selected.index.end(); ++iit){
		//for(int i = 0 ; i < heuristicinfo.pickup.available.index.size(); ++i){
			pair<myfloattype, myinttype> entry; 
			const myinttype pickupnode = *iit;
			entry.second = pickupnode;
			entry.first = TOPD::distancematrix[center][pickupnode];
			heuristicinfo.pickup.selected.mm.insert(entry);
			heuristicinfo.heurval[pickupnode] = entry.first;
			heuristicinfo.pickup.selected.sum += entry.first;
		}
	}
}
	

detourretvalue TOPD::detourlength(myinttype cand){
	detourretvalue ret;
	myfloattype d;
	ret.val = distancematrix[heursol.tour[0]][cand] + distancematrix[cand][heursol.tour[1]] - distancematrix[heursol.tour[0]][heursol.tour[1]];
	ret.index = 0;
	for(int i = 1; i < (heursol.tour.size()-1); i++){
		d = distancematrix[heursol.tour[i]][cand] + distancematrix[cand][heursol.tour[i+1]] - distancematrix[heursol.tour[i]][heursol.tour[i+1]];
		if(d < ret.val){
			ret.val = d;
			ret.index = i;
		}
	}
	return(ret);
}

detourretvalue TOPD::detourlengthold(myinttype cand, int type){
	detourretvalue ret;
	if (type == 1){ // delivery
		myfloattype d;
		ret.val = distancematrix[heursol.tour[0]][cand] + distancematrix[cand][heursol.tour[1]] - distancematrix[heursol.tour[0]][heursol.tour[1]];
		ret.index = 0;
		if( heursol.bufferminus[ret.index] <= 0){
			ret.index = - 1;
		}
		for(int i = 1; i < (heursol.tour.size()-1); i++){
			d = distancematrix[heursol.tour[i]][cand] + distancematrix[cand][heursol.tour[i+1]] - distancematrix[heursol.tour[i]][heursol.tour[i+1]];
			if((ret.index = - 1 && heursol.bufferminus[ret.index] > 0) || (d < ret.val && heursol.bufferminus[ret.index] > 0)){
				ret.val = d;
				ret.index = i;
			}
		}
		return(ret);
	}
	if (type == 2){ // pickup
		myfloattype d;
		ret.val = distancematrix[heursol.tour[0]][cand] + distancematrix[cand][heursol.tour[1]] - distancematrix[heursol.tour[0]][heursol.tour[1]];
		ret.index = 0;
		if( heursol.bufferplus[ret.index] <= 0){
			ret.index = - 1;
		}
		for(int i = 1; i < (heursol.tour.size()-1); i++){
			d = distancematrix[heursol.tour[i]][cand] + distancematrix[cand][heursol.tour[i+1]] - distancematrix[heursol.tour[i]][heursol.tour[i+1]];
			if((ret.index = - 1 && heursol.bufferplus[ret.index] > 0) || (d < ret.val && heursol.bufferplus[ret.index] > 0 )){
				ret.val = d;
				ret.index = i;
			}
		}
		return(ret);
	}
}

detourretvalue TOPD::detourlength(myinttype cand, int type){
	detourretvalue ret;
	if (type == 1){ // delivery
		myfloattype d;
		ret.val = distancematrix[heursol.tour[0]][cand] + distancematrix[cand][heursol.tour[1]] - distancematrix[heursol.tour[0]][heursol.tour[1]];
		ret.index = 0;
		if( heursol.bufferminus[ret.index] < L_toleranceminus + TOL){
			ret.index = - 1;
		}
		if(  heursol.length + ret.val > Tmax + T_tolerance + TOL ){
			ret.index = - 2;
		}
		for(int i = 1; i < (heursol.tour.size()-1); i++){
			d = distancematrix[heursol.tour[i]][cand] + distancematrix[cand][heursol.tour[i+1]] - distancematrix[heursol.tour[i]][heursol.tour[i+1]];
			if( heursol.bufferminus[i] >=  L_toleranceminus + TOL &&   heursol.length + d <= Tmax + T_tolerance  + TOL){
				if(ret.index < 0) {
					ret.val = d;
					ret.index = i;
					continue;
				}
				if(ret.index >=0 && d < ret.val ){
					ret.val = d;
					ret.index = i;
				}
			}
			
		}
		return(ret);
	}
	if (type == 2){ // pickup
		myfloattype d;
		ret.val = distancematrix[heursol.tour[0]][cand] + distancematrix[cand][heursol.tour[1]] - distancematrix[heursol.tour[0]][heursol.tour[1]];
		ret.index = 0;
		if( heursol.bufferplus[ret.index] <  L_toleranceplus + TOL){
			ret.index = - 1;
		} 
		if(  heursol.length + ret.val  > Tmax + T_tolerance ){
			ret.index = - 2;
		}
		for(int i = 1; i < heursol.tour.size()-1; ++i) { // find the earliest tolerable position in the tour
			d = distancematrix[heursol.tour[i]][cand] + distancematrix[cand][heursol.tour[i+1]] - distancematrix[heursol.tour[i]][heursol.tour[i+1]];
			if(ret.index < 0  &&  heursol.bufferplus[i] >=  L_toleranceplus &&   heursol.length + d <= Tmax + T_tolerance){ 
				ret.val = d;
				ret.index = i;
				break;
			}
			// if(ret.index >= 0 && d < ret.val && heursol.bufferplus[ret.index] > 0 )){
			// 	 ret.val = d;
			// 	 ret.index = i;
			//	 continue;
			// }
		}
		return(ret);
	}
}


/*
detourretvalue TOPD::detourlength(myinttype cand, int type){
	detourretvalue ret;
	if (type == 1){ // delivery
		myfloattype d;
		ret.val = distancematrix[heursol.tour[0]][cand] + distancematrix[cand][heursol.tour[1]] - distancematrix[heursol.tour[0]][heursol.tour[1]];
		ret.index = 0;
		if( heursol.bufferminus[ret.index] <= 0){
			ret.index = - 1;
		}
		for(int i = 1; i < (heursol.tour.size()-1); i++){
			d = distancematrix[heursol.tour[i]][cand] + distancematrix[cand][heursol.tour[i+1]] - distancematrix[heursol.tour[i]][heursol.tour[i+1]];
			if((ret.index = - 1 && heursol.load[ret.index] > 0) || (d < ret.val && heursol.load[ret.index] > 0)){
				ret.val = d;
				ret.index = i;
			}
		}
		return(ret);
	}
	if (type == 2){ // pickup
		myfloattype d;
		ret.val = distancematrix[heursol.tour[0]][cand] + distancematrix[cand][heursol.tour[1]] - distancematrix[heursol.tour[0]][heursol.tour[1]];
		ret.index = 0;
		if( heursol.bufferplus[ret.index] <= 0){
			ret.index = - 1;
		}
		for(int i = 1; i < (heursol.tour.size()-1); i++){
			d = distancematrix[heursol.tour[i]][cand] + distancematrix[cand][heursol.tour[i+1]] - distancematrix[heursol.tour[i]][heursol.tour[i+1]];
			if((ret.index = - 1 && heursol.bufferplus[ret.index] > 0) || (d < ret.val && heursol.bufferplus[ret.index] > 0 )){
				ret.val = d;
				ret.index = i;
			}
		}
		return(ret);
	}
}
*/

//static bool compare(myfloattype a, myfloattype b)
//{
//   return (a< b);
//}

template<class _Ty>
	struct mymin
		: public binary_function<_Ty, _Ty, _Ty>
	{	// functor for operator+
	_Ty operator()(const _Ty& _Left, const _Ty& _Right) const
		{	// apply minus
		return (min(_Left,_Right));
		}
	};

int TOPD::insertpickup(sol & solution,detourretvalue ret, myinttype newpickup, myfloattype intensity){
	//if(intensity<1){
		solution.length = solution.length + ret.val;
		//std::vector<myfloattype>::iterator result;
		//result = max_element(solution.load.begin()+ret.index+1, solution.load.end());
		//solution.maxload = max((*result) + intensity * balance[newpickup], solution.maxload);
		myfloattype bmvold = 0.0 + *(solution.bufferminus.begin()+ret.index+1);
		myfloattype bmvnew = intensity * balance[newpickup] + *(solution.bufferminus.begin()+ret.index);
		
		transform(solution.load.begin()+ret.index+1,solution.load.end(), solution.load.begin()+ret.index+1, bind1st(std::plus<myfloattype>(),intensity * balance[newpickup]));  
		solution.load.insert(solution.load.begin()+ret.index+1,*(solution.load.begin()+ret.index)+intensity * balance[newpickup]);
		solution.tour.insert(solution.tour.begin()+ret.index+1,newpickup);
		solution.quantity.insert(solution.quantity.begin()+ret.index+1,intensity * balance[newpickup]);
		solution.intensity.insert(solution.intensity.begin()+ret.index+1,intensity);
		transform (solution.bufferplus.begin()+ret.index,solution.bufferplus.end(), solution.bufferplus.begin()+ret.index, bind1st(std::plus<myfloattype>(),-intensity * balance[newpickup]));
		solution.bufferplus.insert(solution.bufferplus.begin()+ret.index+1,0+*(solution.bufferplus.begin()+ret.index)); 
		transform (solution.bufferplus.begin(),solution.bufferplus.begin()+ret.index+1, solution.bufferplus.begin(), bind1st(mymin<myfloattype>(), *(solution.bufferplus.begin()+ret.index+1))); //
		
		transform (solution.bufferminus.begin()+ret.index+1,solution.bufferminus.end(), solution.bufferminus.begin()+ret.index+1, bind1st(std::plus<myfloattype>(),intensity * balance[newpickup]));
		//transform (solution.bufferminus.begin(),solution.bufferminus.begin()+ret.index+1, solution.bufferminus.begin(), bind1st(mymin<myfloattype>(), *(solution.bufferminus.begin()+ret.index)+intensity * balance[newpickup])); //
		solution.bufferminus.insert(solution.bufferminus.begin()+ret.index+1,bmvnew);
		for(int i = ret.index+1; i >= 0; i--){
			if( solution.load[i] > bmvold ){
				solution.bufferminus[i] = min(solution.load[i],solution.bufferminus[i+1]);
			} else {
				break;
			}
		}

		solution.maxload = Lmax - solution.bufferplus[0];
		solution.profit += intensity * profit[newpickup] * balance[newpickup] - TourCostFactor*ret.val;
	return 1;
}

int TOPD::insertdelivery(sol & solution,detourretvalue ret, myinttype newdelivery, myfloattype intensity){
	solution.length = solution.length + ret.val;
	
	myfloattype bpvold = 0.0 + *(solution.bufferplus.begin()+ret.index+1);
	myfloattype bpvnew = - intensity * balance[newdelivery] + *(solution.bufferplus.begin()+ret.index);
	
	transform(solution.load.begin()+ret.index+1,solution.load.end(), solution.load.begin()+ret.index+1, bind1st(std::plus<myfloattype>(),intensity * balance[newdelivery]));  
	solution.load.insert(solution.load.begin()+ret.index+1,*(solution.load.begin()+ret.index)+intensity * balance[newdelivery]);
	solution.tour.insert(solution.tour.begin()+ret.index+1,newdelivery);
	solution.quantity.insert(solution.quantity.begin()+ret.index+1,intensity * balance[newdelivery]);
	solution.intensity.insert(solution.intensity.begin()+ret.index+1,intensity);
	
	transform (solution.bufferminus.begin()+ret.index+1,solution.bufferminus.end(), solution.bufferminus.begin()+ret.index+1, bind1st(std::plus<myfloattype>(),intensity * balance[newdelivery]));
	solution.bufferminus.insert(solution.bufferminus.begin()+ret.index+1,min(solution.load[ret.index+1], solution.bufferminus[ret.index+1])); 
	transform (solution.bufferminus.begin(),solution.bufferminus.begin()+ret.index+1, solution.bufferminus.begin(), bind1st(mymin<myfloattype>(), *(solution.bufferminus.begin()+ret.index+1))); //
		
	transform (solution.bufferplus.begin()+ret.index+1,solution.bufferplus.end(), solution.bufferplus.begin()+ret.index+1, bind1st(std::plus<myfloattype>(),-intensity * balance[newdelivery]));
	//transform (solution.bufferminus.begin(),solution.bufferminus.begin()+ret.index+1, solution.bufferminus.begin(), bind1st(mymin<myfloattype>(), *(solution.bufferminus.begin()+ret.index)+intensity * balance[newpickup])); //
	solution.bufferplus.insert(solution.bufferplus.begin()+ret.index+1,bpvnew);
	for(int i = ret.index; i >= 0; --i){
		if( Lmax - solution.load[i] > bpvold ){
			solution.bufferplus[i] = min(Lmax - solution.load[i],solution.bufferplus[i+1]);
		} else {
			break;
		}
	}

	//transform (solution.bufferplus.begin()+ret.index+1,solution.bufferplus.end(), solution.bufferplus.begin()+ret.index+1, bind1st(std::plus<myfloattype>(),-intensity * balance[newpickup]));
	//solution.bufferplus.insert(solution.bufferplus.begin()+ret.index+1,*(solution.bufferplus.begin()+ret.index)-intensity * balance[newdelivery]);
	//transform (solution.bufferminus.begin()+ret.index+1,solution.bufferminus.end(), solution.bufferminus.begin()+ret.index+1, bind1st(std::plus<myfloattype>(),intensity * balance[newpickup]));
	//transform (solution.bufferminus.begin(),solution.bufferminus.begin()+ret.index+1, solution.bufferminus.begin(), bind1st(mymin<myfloattype>(), *(solution.bufferminus.begin()+ret.index)+intensity * balance[newpickup])); //
	//solution.bufferminus.insert(solution.bufferminus.begin()+ret.index+1,*(solution.bufferminus.begin()+ret.index)+intensity * balance[newpickup]);
	
	//result = max_element(solution.load.begin()+ret.index+1, solution.load.end());
	//result2 = max_element(solution.load.begin(), solution.load.begin()+ret.index+1);
	//solution.maxload = max((*result) + intensity * balance[newdelivery], (*result2));
	
	solution.maxload = Lmax - solution.bufferplus[0];
	solution.profit += intensity * profit[newdelivery] * balance[newdelivery] - TourCostFactor*ret.val;
	return 1;
}

//typename vector<T>::iterator it

template <class T> void shift(vector<T> &c1,  typename vector<T>::iterator  it1, vector<T> &c2, typename vector<T>::iterator it2)
{
	c2.insert(it2,*it1);
    c1.erase(it1);
}

bool myfunction (myfloattype i,myfloattype j) { return (i>j); }

infomove TOPD::addpickup(sol & solution){
	infomove retval;
	retval.type = 1;
	if(heuristicinfo.pickup.available.index.size() == 0){	 
		retval.status = -1;
		return retval;
	} 

	const set <myinttype>::iterator ite = heuristicinfo.pickup.available.index.end();
	set <myinttype>::iterator it = heuristicinfo.pickup.available.index.begin();

	int N = heuristicinfo.pickup.available.index.size();

	const multimap<myfloattype, myinttype,std::greater<myfloattype>>::iterator iiend = heuristicinfo.pickup.available.mm.end();
	multimap<myfloattype, myinttype,std::greater<myfloattype>>::iterator ii=heuristicinfo.pickup.available.mm.begin();
	int mmn = heuristicinfo.pickup.available.mm.size(); 
	vector<int> skip(mmn,0); //skip(mmn,1);	skip[2]=0;
	int nskip = 0;

	myfloattype vsum = heuristicinfo.pickup.available.sum;
	int test = 1; 

	// for(int k=0;k<98;++k){
	//   cout << (float) rand()/RAND_MAX << endl;
	// }

	myinttype node;
	detourretvalue ret;
	int i;

	while(test==1 && nskip<mmn){
		float rnd = ((float)rand()/RAND_MAX);
		myfloattype	rnd2 =  rnd*vsum;
		myfloattype psum = 0;
		i = 0;
		int iend = mmn;
		int enter = mmn+1;
		test = -1;
		for( ; ii!=iiend && i < iend; ){
			psum += (*ii).first;
			if(rnd2 <= psum){
				if(enter == mmn+1){
					enter = i; //remember the entering point
				}
				if(skip[i]==0){
					//cout << "  [" << (*ii).first << ", " << (*ii).second << "]" << endl;
					node = (*ii).second;
					ret = detourlength(node,2);
					if(solution.length + ret.val <= Tmax && ret.index != -1){
						test=0; // selected node is o.k. - exit for and while loop
					} else {
						skip[i]=1;
						++nskip;
						test = 1; // exit for loop and select a new rand candidate 
					}
					break;
				}
				if(i == mmn-1){
					i = 0; 
					ii=heuristicinfo.pickup.available.mm.begin();
					iend = enter;
					test = 1;
					continue; // restart for loop at the beginning and search unil entering point is reached
				}
			}
			++ii;
			++i;
		}
	}

	if(test == 0){	
		retval.node = node;
		retval.nodepos = heuristicinfo.pickup.selected.index.size();   //cout << ret.val << endl; ret.index << endl;
		retval.tourpos = ret.index;
		//retval.heurval = heuristicinfo.pickup.available.value[i];
		retval.heurval = heuristicinfo.heurval[node];
		heuristicinfo.type[node] = 2;

		myfloattype qty = min( heursol.bufferplus[ret.index], balance[node])/balance[node]; 
		insertpickup(solution,ret,node,qty);

		// heuristicinfo.pickup.available.mm.insert(
		heuristicinfo.pickup.selected.mm.insert(*ii);
		heuristicinfo.pickup.selected.sum += retval.heurval; 
        heuristicinfo.pickup.available.mm.erase(ii);
		heuristicinfo.pickup.available.sum -= retval.heurval; 
		
		heuristicinfo.pickup.selected.index.insert(node);
		heuristicinfo.pickup.available.index.erase(node);
		//shift(heuristicinfo.pickup.available.index, heuristicinfo.pickup.available.index.begin() + i ,heuristicinfo.pickup.selected.index, heuristicinfo.pickup.selected.index.end()); //begin()+insertindex);
		// shift(heuristicinfo.pickup.available.value, heuristicinfo.pickup.available.value.begin() + i ,heuristicinfo.pickup.selected.value, heuristicinfo.pickup.selected.value.end()); //+insertindex);
		
		retval.status = 1;
		// undo_addpickup(retval);
		// deletenode(solution, ret.index);
		return retval;
	} else {
		retval.status = 0;
		return retval;
	}
}

/* {
	while( i < N){
		myinttype node = heuristicinfo.pickup.available.index[i];

		detourretvalue ret = detourlength(node,2);
		if( solution.length + ret.val > Tmax || ret.index == -1){
			i++;
			continue;
		};
		retval.node = node;
		retval.nodepos = i;   //cout << ret.val << endl; ret.index << endl;
		retval.tourpos = ret.index;
		retval.heurval = heuristicinfo.pickup.available.value[i];

		myfloattype qty = min( heursol.bufferplus[ret.index], balance[node])/balance[node]; 
		insertpickup(solution,ret,node,qty);

		// heuristicinfo.pickup.available.mm.insert(
		// c2.insert(it2,*it1);
        // c1.erase(it1);
		
		shift(heuristicinfo.pickup.available.index, heuristicinfo.pickup.available.index.begin() + i ,heuristicinfo.pickup.selected.index, heuristicinfo.pickup.selected.index.end()); //begin()+insertindex);
		shift(heuristicinfo.pickup.available.value, heuristicinfo.pickup.available.value.begin() + i ,heuristicinfo.pickup.selected.value, heuristicinfo.pickup.selected.value.end()); //+insertindex);
		
		

		retval.status = 1;
		// undo_addpickup(retval);
		// deletenode(solution, ret.index);
		return retval;
	}
	retval.status = 0;
	return retval;
} */


int TOPD::undo_addpickup(infomove retval)
{
   //std::vector<myfloattype>::iterator low=std::lower_bound(heuristicinfo.pickup.available.value.begin(), heuristicinfo.pickup.available.value.end(), retval.heurval,myfunction); //
   //int index = (low- heuristicinfo.pickup.available.value.begin());
   //heuristicinfo.pickup.available.value.insert( heuristicinfo.pickup.available.value.begin() + index ,retval.heurval);
   //heuristicinfo.pickup.available.index.insert( heuristicinfo.pickup.available.index.begin() + index ,retval.node);
   //heuristicinfo.pickup.selected.index.erase(heuristicinfo.pickup.selected.index.begin()+retval.nodepos);
   //heuristicinfo.pickup.selected.value.erase(heuristicinfo.pickup.selected.value.begin()+retval.nodepos);

   heuristicinfo.pickup.available.mm.insert(pair<myfloattype, myinttype>(retval.heurval,retval.node));
   const pair < multimap<myfloattype, myinttype,std::greater<myfloattype>>::iterator,multimap<myfloattype, myinttype,std::greater<myfloattype>>::iterator > mmrange = heuristicinfo.pickup.selected.mm.equal_range(retval.heurval);
   multimap<myfloattype, myinttype,std::greater<myfloattype>>::iterator it = mmrange.first;
   for (; it!=mmrange.second; ++it){
	   if((*it).second == retval.node){
		   heuristicinfo.pickup.selected.mm.erase(it);
		   break;
	   }
   }
   heuristicinfo.pickup.available.sum += retval.heurval; 
   heuristicinfo.pickup.selected.sum -= retval.heurval; 
   heuristicinfo.type[retval.node] = 1;
   heuristicinfo.pickup.available.index.insert(retval.node);
   heuristicinfo.pickup.selected.index.erase(retval.node);
   return 1;
}

int TOPD::undo_adddelivery(infomove retval)
{
  // std::vector<myfloattype>::iterator low=std::lower_bound(heuristicinfo.delivery.available.value.begin(), heuristicinfo.delivery.available.value.end(), retval.heurval,myfunction); //
  // int index = (low- heuristicinfo.delivery.available.value.begin());
  // heuristicinfo.delivery.available.value.insert(heuristicinfo.delivery.available.value.begin() + index ,retval.heurval);
  // heuristicinfo.delivery.available.index.insert(heuristicinfo.delivery.available.index.begin() + index ,retval.node);

  // heuristicinfo.delivery.selected.index.erase(heuristicinfo.delivery.selected.index.begin()+retval.nodepos);
  // heuristicinfo.delivery.selected.value.erase(heuristicinfo.delivery.selected.value.begin()+retval.nodepos);

   heuristicinfo.delivery.available.mm.insert(pair<myfloattype, myinttype>(retval.heurval,retval.node));
   const pair < multimap<myfloattype, myinttype,std::greater<myfloattype>>::iterator,multimap<myfloattype, myinttype,std::greater<myfloattype>>::iterator > mmrange = heuristicinfo.delivery.selected.mm.equal_range(retval.heurval);
   multimap<myfloattype, myinttype,std::greater<myfloattype>>::iterator it = mmrange.first;
   for (; it!=mmrange.second; ++it){
	   if((*it).second == retval.node){
		   heuristicinfo.delivery.selected.mm.erase(it);
		   break;
	   }
   }
   heuristicinfo.delivery.available.sum += retval.heurval; 
   heuristicinfo.delivery.selected.sum -= retval.heurval; 
   heuristicinfo.type[retval.node] = -1;
   heuristicinfo.delivery.available.index.insert(retval.node);
   heuristicinfo.delivery.selected.index.erase(retval.node);
   return 1;
}

// deletes index + 1 from tour plus corresponding changes on sol struct
int TOPD::deletenode(sol & solution, myinttype index)
{
	myinttype node = solution.tour[index+1];
	myfloattype qty = solution.intensity[index+1] * balance[node];
	myfloattype intensity = - solution.intensity[index+1];
	detourretvalue ret;
	ret.index = index;
	ret.val = distancematrix[solution.tour[index]][solution.tour[index+2]]- distancematrix[solution.tour[index]][node] - distancematrix[node][solution.tour[index+2]];
	solution.length = solution.length + ret.val;
	
	if(qty>0){	
	    myinttype newdelivery = node;
		myfloattype bpvold = 0.0 + *(solution.bufferplus.begin()+ret.index+1);
		myfloattype bpvnew = - intensity * balance[newdelivery] + *(solution.bufferplus.begin()+ret.index+1);
	
		transform(solution.load.begin()+ret.index+1,solution.load.end(), solution.load.begin()+ret.index+1, bind1st(std::plus<myfloattype>(),intensity * balance[newdelivery]));  
		transform(solution.bufferminus.begin()+ret.index+1,solution.bufferminus.end(), solution.bufferminus.begin()+ret.index+1, bind1st(std::plus<myfloattype>(),intensity * balance[newdelivery]));

		transform(solution.bufferminus.begin(),solution.bufferminus.begin()+ret.index+1, solution.bufferminus.begin(), bind1st(mymin<myfloattype>(), *(solution.bufferminus.begin()+ret.index+1)));		
		transform(solution.bufferplus.begin()+ret.index+1,solution.bufferplus.end(), solution.bufferplus.begin()+ret.index+1, bind1st(std::plus<myfloattype>(),-intensity * balance[newdelivery]));

		for(int i = ret.index+1; i >= 0; i--){
			if( Lmax - solution.load[i] > bpvold ){
				solution.bufferplus[i] = min(Lmax - solution.load[i],solution.bufferplus[i+1]);
			} else {
				break;
			}
		}
		
		solution.load.erase(solution.load.begin()+ret.index+1);
		solution.tour.erase(solution.tour.begin()+ret.index+1);
		solution.quantity.erase(solution.quantity.begin()+ret.index+1);
		solution.intensity.erase(solution.intensity.begin()+ret.index+1);	
		solution.bufferminus.erase(solution.bufferminus.begin()+ret.index+1); 
		solution.bufferplus.erase(solution.bufferplus.begin()+ret.index+1);

		solution.maxload = Lmax - solution.bufferplus[0];
		solution.profit += intensity * profit[newdelivery] * balance[newdelivery] - TourCostFactor*ret.val;
		
		//insertdelivery(solution,ret,node,-solution.intensity[index+1]);
		return 1;
	} else {
        myinttype newpickup = node;
       
		myfloattype bmvold = 0.0 + *(solution.bufferminus.begin()+ret.index+1);
		myfloattype bmvnew = intensity * balance[newpickup] + *(solution.bufferminus.begin()+ret.index+1);
		
		transform(solution.load.begin()+ret.index+1,solution.load.end(), solution.load.begin()+ret.index+1, bind1st(std::plus<myfloattype>(),intensity * balance[newpickup]));  
			
		transform (solution.bufferplus.begin()+ret.index+1,solution.bufferplus.end(), solution.bufferplus.begin()+ret.index+1, bind1st(std::plus<myfloattype>(),-intensity * balance[newpickup]));
		transform (solution.bufferplus.begin(),solution.bufferplus.begin()+ret.index+1, solution.bufferplus.begin(), bind1st(mymin<myfloattype>(), *(solution.bufferplus.begin()+ret.index+1))); //
		
		transform (solution.bufferminus.begin()+ret.index+1,solution.bufferminus.end(), solution.bufferminus.begin()+ret.index+1, bind1st(std::plus<myfloattype>(),intensity * balance[newpickup]));

		for(int i = ret.index+1; i >= 0; i--){
			if( solution.load[i] > bmvold ){
				solution.bufferminus[i] = min(solution.load[i],solution.bufferminus[i+1]);
			} else {
				break;
			}
		}
		
		solution.load.erase(solution.load.begin()+ret.index+1);
		solution.quantity.erase(solution.quantity.begin()+ret.index+1);
		solution.tour.erase(solution.tour.begin()+ret.index+1);
		solution.intensity.erase(solution.intensity.begin()+ret.index+1);
		solution.bufferplus.erase(solution.bufferplus.begin()+ret.index+1);
		solution.bufferminus.erase(solution.bufferminus.begin()+ret.index+1);

		solution.maxload = Lmax - solution.bufferplus[0];
		solution.profit += intensity * profit[newpickup] * balance[newpickup] - TourCostFactor*ret.val;

		// insertpickup(solution,ret,node,-solution.intensity[index+1]);
		return 1;

	}
   
}

infomove TOPD::adddelivery(sol & solution){

	infomove retval;
	retval.type = -1;
	if(heuristicinfo.delivery.available.index.size() == 0){	 
		retval.status = -1;
		return retval;
	} 

	const set <myinttype>::iterator ite = heuristicinfo.delivery.available.index.end();
	set <myinttype>::iterator it = heuristicinfo.delivery.available.index.begin();

	int N = heuristicinfo.delivery.available.index.size();

	const multimap<myfloattype, myinttype,std::greater<myfloattype>>::iterator iiend = heuristicinfo.delivery.available.mm.end();
	multimap<myfloattype, myinttype,std::greater<myfloattype>>::iterator ii=heuristicinfo.delivery.available.mm.begin();
	int mmn = heuristicinfo.delivery.available.mm.size(); 
	vector<int> skip(mmn,0); //skip(mmn,1);	skip[2]=0;
	int nskip = 0;

	myfloattype vsum = heuristicinfo.delivery.available.sum;
	int test = 1; 

	myinttype node;
	detourretvalue ret;
	int i;

	while(test==1 && nskip<mmn){
		float rnd = ((float)rand()/RAND_MAX);
		myfloattype	rnd2 =  rnd*vsum;
		myfloattype psum = 0;
		i = 0;
		int iend = mmn;
		int enter = mmn+1;
		test = -1;
		for( ; ii!=iiend && i < iend; ){
			psum += (*ii).first;
			if(rnd2 <= psum){
				if(enter == mmn+1){
					enter = i; //remember the entering point
				}
				if(skip[i]==0){
					//cout << "  [" << (*ii).first << ", " << (*ii).second << "]" << endl;
					node = (*ii).second;
					ret = detourlength(node,1);
					if(solution.length + ret.val <= Tmax && ret.index != -1){
						test=0; // selected node is o.k. - exit for and while loop
					} else {
						skip[i]=1;
						++nskip;
						test = 1; // exit for loop and select a new rand candidate 
					}
					break;
				}
				if(i == mmn-1){
					i = 0; 
					ii=heuristicinfo.delivery.available.mm.begin();
					iend = enter;
					test = 1;
					continue; // restart for loop at the beginning and search unil entering point is reached
				}
			}
			++ii;
			++i;
		}
	}

	if(test == 0){	
		retval.node = node;
		retval.nodepos = heuristicinfo.delivery.selected.index.size();  //cout << ret.val << endl; ret.index << endl;
		retval.tourpos = ret.index;
		//retval.heurval = heuristicinfo.delivery.available.value[i];
		retval.heurval = heuristicinfo.heurval[node];
		heuristicinfo.type[node] = -2;

		myfloattype qty = -min(heursol.bufferminus[ret.index], -balance[node])/balance[node]; 
		insertdelivery(solution,ret,node,qty);

		// heuristicinfo.delivery.available.mm.insert(
		heuristicinfo.delivery.selected.mm.insert(*ii);
		heuristicinfo.delivery.selected.sum += retval.heurval; 
        heuristicinfo.delivery.available.mm.erase(ii);
		heuristicinfo.delivery.available.sum -= retval.heurval; 
		
	    heuristicinfo.delivery.selected.index.insert(node);
		heuristicinfo.delivery.available.index.erase(node);
		// shift(heuristicinfo.delivery.available.index, heuristicinfo.delivery.available.index.begin() + i ,heuristicinfo.delivery.selected.index, heuristicinfo.delivery.selected.index.end()); //begin()+insertindex);
		// shift(heuristicinfo.delivery.available.value, heuristicinfo.delivery.available.value.begin() + i ,heuristicinfo.delivery.selected.value, heuristicinfo.delivery.selected.value.end()); //+insertindex);
		
		retval.status = 1;
		// undo_adddelivery(retval);
		// deletenode(solution, ret.index);
		return retval;
	} else {
		retval.status = 0;
		return retval;
	}

}

int TOPD::deletenode(infomove retval){
	if(balance[retval.node] > 0){
		undo_addpickup(retval); 
	} else {
		undo_adddelivery(retval); 
	}
	deletenode(heursol,retval.tourpos);
	return 1;
}

int TOPD::deletenode(myinttype tourpos){
	infomove retval;
	retval.node = heursol.tour[tourpos+1]; 
 	retval.heurval = heuristicinfo.heurval[heursol.tour[tourpos+1]];
	if(balance[retval.node] > 0){
		undo_addpickup(retval); 
	} else {
		undo_adddelivery(retval); 
	}
	deletenode(heursol,tourpos);
	return 1;
}


TOPD::TOPD(TOPDparameters parameters){
	setfilename(parameters.filename,parameters.filename2,parameters.path);
	setRpar(parameters.rexe,parameters.rpath);
	load();
	setTourCostFactor(parameters.TourCostFactor);
	#ifdef  CPLEXAVAILABLE
	setComptime(parameters.ILOMaxCompTime);
	setILOParHeur(parameters.setILOParHeur);
	#endif
	initheursol();
}

void TOPD::initheursol(){
	L_toleranceplus = 0;
	L_toleranceminus = 0;
	T_tolerance = 0;
	elitists.lb = 0;
	elitists.ub = 0;
	heursol.intensity.push_back(1);
	heursol.intensity.push_back(1);
	heursol.tour.push_back(0);
	heursol.tour.push_back(1);
	
	heursol.length = calclength(heursol.tour);

	heursol.profit = -TourCostFactor * heursol.length;
	heursol.maxload = 0;
	heursol.quantity.push_back(0); 
	heursol.quantity.push_back(0);
	heursol.load.push_back(0); 
	heursol.load.push_back(0);
	heursol.bufferplus.push_back(Lmax); 
	heursol.bufferplus.push_back(Lmax);
	heursol.bufferminus.push_back(0); 
	heursol.bufferminus.push_back(0);
	heursol.topd = this;
	initializheuristicinfo();
}

void TOPD::initsol(sol & solution){
	solution.intensity.push_back(1);
	solution.intensity.push_back(1);
	solution.tour.push_back(0);
	solution.tour.push_back(1);
	
	#ifdef CPLEXAVAILABLE
	solution.length = calclength(heursol.tour);
	#endif

	solution.profit = -TourCostFactor * heursol.length;
	solution.maxload = 0;
	solution.quantity.push_back(0); 
	solution.quantity.push_back(0);
	solution.load.push_back(0); 
	solution.load.push_back(0);
	solution.bufferplus.push_back(Lmax); 
	solution.bufferplus.push_back(Lmax);
	solution.bufferminus.push_back(0); 
	solution.bufferminus.push_back(0);
	solution.topd = this;
}

#ifdef CPLEXAVAILABLE
void TOPD::setComptime(myfloattype number){
	ILOMaxCompTime = number;
}



void TOPD::setILOParHeur(string par){
	if(par.compare("heu")){
		ILOParHeur = 1;
	} else {
		ILOParHeur = 0;
	}
}
#endif

void TOPD::setTourCostFactor(myfloattype number){
	TourCostFactor = number;
	#ifdef CPLEXAVAILABLE
	ILOTourCostFactor = number;
	#endif
}

int TOPD::checkcalcbufferplus(){
	myfloattype val = heursol.load[heursol.load.size()-1];
	int i;
	for(i =  heursol.bufferplus.size()-1; i > 0; ){
		if(abs( (Lmax-heursol.bufferplus[i]) - val) > TOL){
			return 0;
		}
		--i;
		val = max(heursol.load[i],val);
	}
	if(abs((Lmax-heursol.bufferplus[i]) - val) > TOL){
		return 0;
	}
	return 1;
}

#ifdef CPLEXAVAILABLE
myfloattype TOPD::calcmaxload(){
	ILOsol.length = 0;
	vector<myfloattype> load(ILOsol.tour.size());
	partial_sum(ILOsol.quantity.begin(), ILOsol.quantity.end(),load.begin());
	ILOsol.maxload = *(max_element(load.begin(), load.end()));
	vector<myfloattype> profits;
	select(ILOsol.tour,profit,profits);
	return ILOsol.maxload;
}

myfloattype TOPD::calclength(){
	myfloattype length = 0;
	for(int i = 1; i < ILOsol.tour.size(); i++){
		length += distancematrix[ILOsol.tour[i-1]][ILOsol.tour[i]];;
	}
	return length;
}
#endif

myfloattype TOPD::calcprofit(){
	myfloattype vprofit = 0;
	for(int i = 1; i < (heursol.tour.size()-1); i++){
		vprofit += balance[heursol.tour[i]]*profit[heursol.tour[i]]*heursol.intensity[i];
	}
	return vprofit -  TourCostFactor * calclength(heursol.tour);
}

myfloattype TOPD::calcprofit(sol & solution){
	myfloattype vprofit = 0;
	for(int i = 1; i < (solution.tour.size()-1); i++){
		vprofit += (*solution.topd).balance[solution.tour[i]]*profit[solution.tour[i]]*solution.intensity[i];
	}
	return vprofit -  TourCostFactor * calclength(solution.tour);
}


myfloattype TOPD::calclength(vector<myinttype> tour){
	myfloattype length = 0;
	for(int i = 1; i < tour.size(); i++){
		length += distancematrix[tour[i-1]][tour[i]];
	}
	return length;
}

#ifdef  CPLEXAVAILABLE
void TOPD::ILOcreateandsolve(){
	try {
		convert2ILO();
		IloModel model1(env);
		model = model1;
		IloInt i, j;

		IloInt	nbRetailer;
	 
		nbRetailer   = balance.size(); 
	 
		IloBool consistentData = (ILOdistancematrix.getSize() == nbRetailer);
		for(i = 0; consistentData && (i < nbRetailer); i++)
			consistentData = (ILOdistancematrix[i].getSize() == nbRetailer);
		if (!consistentData) {
			cerr << "ERROR: data file '" 
				<< filename << "' contains inconsistent data" << endl;
			throw(-1);
		}

		IloNumVarArray visit(env, nbRetailer, 0, 1, ILOINT);
		IloNumVarArray y(env, nbRetailer, 0, 1);
		NumVarMatrix x(env, nbRetailer);
      
		for(i = 0; i < nbRetailer; i++){
			x[i] = IloNumVarArray(env, nbRetailer, 0, 1, ILOINT);
		}

		IloNumVarArray start(env, nbRetailer, 0, 1000); 
		IloNumVarArray load(env, nbRetailer, 0, 300); //maximum capacity of the vehicle
		IloNumVar tourdesc(env, 0, 1, ILOINT);

		// for(i = 0; i < nbRetailer; i++)
		//      model.add(IloSum(x[i]) - x[i][i] == visit[i]);
	 
		// visited nodes have an out-degree 1 
		for(j = 0; j < nbRetailer; j++) {
			IloExpr v(env);
			for(i = 0; i < nbRetailer; i++)
			if(i != j) {
				v += x[j][i];
			}
			model.add(v == visit[j]);
			v.end();
		}

		/*
		for(j = 6; j < 17; j++) {
			model.add( visit[j] == 0);
		}
		*/

		// you may load or unload parts at visited nodes 
		for(i = 0; i < nbRetailer; i++) {
			model.add( y[i] <= visit[i] );
		}

		// visited nodes have an in-degree 1 
		for(j = 0; j < nbRetailer; j++) {
			IloExpr v(env);
			for(i = 0; i < nbRetailer; i++)
			if(i != j) {
				v += x[i][j];
			}
			model.add(v == visit[j]);
			v.end();
		}

		// sub tour elimination 
		for(j = 0; j < nbRetailer; j++) {
			for(i = 0; i < nbRetailer; i++)
			{
				if(i != j && j!= 0) { // the predecessors of 0 cannot have a lower rank 
					//model.add( start[i] + ILOdistance[i][j] <= start[j] + 1000*(1-x[i][j]) );
					model.add( start[i] + 1 <= start[j] + 1000*(1-x[i][j]) );
				}
			}    
		}

		//pick up and delivery 
		for(j = 0; j < nbRetailer; j++) {
			for(i = 0; i < nbRetailer; i++)
			{
				if(i != j && j != 0) { // the load equation need not be satisfied for the return
					model.add( load[j] <= load[i] + ILObalance[j]*y[j] + 1000*(1-x[i][j]) );
					model.add( load[i] + ILObalance[j]*y[j] <=  load[j] + 1000*(1-x[i][j]) );
				}
			}    
		}

		// initial load is zero 
		model.add( load[0] == 0 );

		// start and end 
		model.add( visit[0] == tourdesc);
		model.add( visit[1] == tourdesc);
		model.add( x[1][0] == tourdesc);
	  
		// maximum tour length 
		IloExpr v(env);
		for(i = 0; i < nbRetailer; i++) {
			v += IloScalProd(ILOdistancematrix[i], x[i]);
		}
		model.add(v <= Tmax);
		v.end();


		IloNumArray balancevalue(env);
		for(i = 0; i < nbRetailer; i++) {
			balancevalue.add( ILOprofit[i] * ILObalance[i]);
		}

		// objective 
		IloExpr obj = IloScalProd(balancevalue, y);
	 
		for(i = 0; i < nbRetailer; i++) {
			obj -= ILOTourCostFactor*IloScalProd(ILOdistancematrix[i], x[i]);
		}
		obj += ILOTourCostFactor*ILOdistancematrix[1][0]*tourdesc; // compensation: there is no need to return

		model.add(IloMaximize(env, obj));
		obj.end();

		IloCplex cplex1(env);
		cplex = cplex1;
		if (ILOParHeur > 0) {
			cplex.setParam( IloCplex::FPHeur,2); 
			cplex.setParam( IloCplex::PolishAfterIntSol,1);
		}

		cplex.extract(model);

		// cplex.setParam( IloCplex::TiLim,1200 ); 
		// FPHeur = 2

		cplex.setParam( IloCplex::TiLim,ILOMaxCompTime);
		cplex.solve();
		nbRetailer = balance.size();

		IloNum tolerance = cplex.getParam(IloCplex::EpInt);
		vector<int> vstart;
		vector<int> vend;
		
		cplex.out() << "Optimal value: " << cplex.getObjValue() << endl;

		// write solution

		ILOsol.profit =  cplex.getObjValue();
		ILOsol.upperbound = cplex.getBestObjValue();

		int k = 0;	
		
		for(j=0; j < nbRetailer; j++) {
			if (cplex.getValue(visit[j]) >= 1 - tolerance) {
				break;
			}
		}
		IloInt tourstart = j; 
		ILOsol.tour.push_back(j);
		ILOsol.intensity.push_back(cplex.getValue(y[j]));
		ILOsol.quantity.push_back(cplex.getValue(y[j])*balance[j]);
		k++;

		int crit = 1;
		while(crit == 1){
			crit = 0;
			for (i = 0; i < nbRetailer; i++) {
				if (i != j && cplex.getValue(x[j][i]) >= 1 - tolerance){
					if(tourstart != i){
						vstart.push_back(j);
						vend.push_back(i);
						j = i;
						ILOsol.tour.push_back(j);
						ILOsol.intensity.push_back(cplex.getValue(y[j]));
						ILOsol.quantity.push_back(cplex.getValue(y[j])*balance[j]);
						crit = 1;
						k++;
					}
					break;
				}
			}
		}
		cout << ILOsol.tour[1];

		calclength();
		calcmaxload();

		//ILOsol.intensity;
		//ILOsol.tour;
		//ILOsol.length;
		//ILOsol.maxload;
	}
	catch(IloException& e) {
		cerr  << " ERROR: " << e << endl;   
	}
	catch(...) {
		cerr  << " ERROR" << endl;   
	}
}
#endif

#ifdef  CPLEXAVAILABLE
TOPD::~TOPD(){
	env.end();
}
#endif

myfloattype TOPD::calcdist(point a, point b){
	return pow((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y),0.5); 
}

void TOPD::setfilename(string fname1,string fname2 , string path1){
	filename =	fname1;
	filename2 =	fname2;
	path =	path1;
	n = 0;
}

void TOPD::setRpar(string str1,string str2){
	rexe = str1;
	rpath = str2;
}

void TOPD::load(){
	string line;

	if(n == 0){ //count number of customers 
		ifstream myfile (filename);
		if (myfile.is_open())
		{
			while ( myfile.good() )
			{
				getline (myfile,line);
				if(!line.empty()){
					n++;
				}
			}
			myfile.close();
			n = n - 2;
		}
		else cout << "Unable to open file";
	}

	ifstream myfile (filename);
	start = 0; //index of start node
	end = 1; //index of end node 
	if (myfile.is_open())
	{
		if (myfile.good() )
		{
			getline (myfile,line);
			//cout << line << endl;
			istringstream iss;
			iss = istringstream(line);
			myfloattype s;
			if (!(iss).eof())
			{
				iss >> s;
				Tmax = s;
			}
			//delete iss;
			getline (myfile,line);
			iss = istringstream(line);
			//cout << line << endl;

			if (!iss.eof())
			{
				iss >> s;
				Lmax = s;
			}

			myinttype index = 0;
			while (!myfile.eof()){
				point mypoint;
				myfloattype xfloat;
				
				getline (myfile,line);

				if(!line.empty()){
					iss = istringstream(line);
					iss >> xfloat; 
					mypoint.x = xfloat;
					iss >> xfloat; 
					mypoint.y = xfloat;

					location.push_back(mypoint);
				
					iss >> xfloat;
					if(xfloat > 0){
						pickup.push_back(index);
					} else
					{
						delivery.push_back(index);
					}
					index ++;
					balance.push_back(xfloat);
				
					iss >> xfloat; 
					profit.push_back(xfloat);
				}
			}
		}
		myfile.close();
		calcdistmat();
	}
	else cout << "Unable to open file";
	cout << "end reading matrix from file" << endl;
}

void TOPD::calcdistmat(){
	myinttype i,j;
	for(i = 0; i < n; i++){
		vector<myfloattype> v(n);
		distancematrix.push_back(v);
	}	
	for(i = 0; i < n; i++){
		distancematrix[i][i]=0;
		for(j = i+1; j < n; j++){
			myfloattype d = calcdist(location[i],location[j]);
			distancematrix[i][j]=d;
			distancematrix[j][i]=d;
		}
	}	
}

void TOPD::displayfile(){
  string line;
  ifstream myfile (filename);
  if (myfile.is_open())
  {
    while ( myfile.good() )
    {
      getline (myfile,line);
      cout << line << endl;
    }
    myfile.close();
  }
  else cout << "Unable to open file";
}

void TOPD::Rprintsol(string rexe, string rpath, string filename, string extension,sol solution){
	
	solution.savesol( rpath + "example1.txt");
	Rsavesinstance(rpath + "example2.txt");		

	int lastchar = filename.find_first_of(".")-1;
	string filetrunc = string(filename.substr(0,lastchar));
	string fnamesolcplex =  filetrunc + "_" + extension + ".pdf";

	string syscommand(rexe);
	syscommand = syscommand + " -f " + rpath + "start.r";
	replace(syscommand.begin(),syscommand.end(),'/','\\');
	system(syscommand.c_str());
	syscommand = std::string("echo ") + syscommand;
	system(syscommand.c_str());
	// system("\"C:/Program Files/R/R-2.15.0/bin/R\" -f C:/Users/martin/Desktop/projects/pickupanddeliveryteamorienteering/R/start.r");
		
	syscommand = string("del ");
	syscommand = syscommand + rpath + fnamesolcplex;
	replace(syscommand.begin(),syscommand.end(),'/','\\');
	system(syscommand.c_str());
	syscommand = std::string("echo ") + syscommand;
	system(syscommand.c_str());
	
	syscommand = string("ren ");
	syscommand = syscommand + rpath + "test123.pdf " + fnamesolcplex;
	replace(syscommand.begin(),syscommand.end(),'/','\\');
	system(syscommand.c_str()); 
	syscommand = std::string("echo ") + syscommand;
	system(syscommand.c_str());
}

int loaddata(string filename) {
  string line;
  ifstream myfile (filename);
  if (myfile.is_open())
  {
    while ( myfile.good() )
    {
      getline (myfile,line);
      //cout << line << endl;
    }
    myfile.close();
  }
  else cout << "Unable to open file"; 
  return 0;
}
    
double to_double ( const char *p )
{
	std::stringstream ss ( p );
	double result = 0;
	ss>> result;
	return result;
}

template<class T> class sorter {
    const vector<T> &values;
public:
    sorter(const std::vector<T> &v) : values(v) {}
    bool operator()(int a, int b) { return values[a] < values[b]; }
};

template<class T> std::vector<int> order(const std::vector<T> &values)
{
    std::vector<int> rv(values.size());
    int idx = 0;
    for (std::vector<int>::iterator i = rv.begin(); i != rv.end(); i++)
        *i = idx++;
    std::sort(rv.begin(), rv.end(), sorter<T>(values));
    return rv;
}

#ifdef CPLEXAVAILABLE
void TOPD::convert2ILO(){
	ILObalance = IloNumArray(env, balance.size());
	ILOprofit = IloNumArray(env, profit.size());
	ILOdistancematrix = IloArray <IloNumArray> (env, profit.size());
	for(int i = 0; i<balance.size(); i++){
		ILObalance[i] = balance[i];
		ILOprofit[i] = profit[i];
	}
	
	for(int i = 0; i<balance.size(); i++){
		ILOdistancematrix[i] = IloNumArray(env, profit.size());
		for(int j = 0; j<balance.size(); j++){
			ILOdistancematrix[i][j] = distancematrix[i][j];
		}
	}
}
#endif

void TOPD::Rsavesinstance(string fname){
		ofstream myfile;
		myfile.open (fname);
		myfile << "tmax = " << Tmax << endl;
		myfile << "lmax = " << Lmax << endl;
		myfile << "instance1 = rbind(" << endl;
		for(int i = 0; i < n-1; i++){
			myfile << "c(";
			myfile << location[i].x << ", " << location[i].y<< ", " << balance[i]<< ", "  << profit[i] ;
			myfile << "), " << endl;
		}
		myfile << "c(" << location[n-1].x << ", " << location[n-1].y<< ", " << balance[n-1]<< ", "  << profit[n-1]  << endl;
		myfile << ") ) " << endl;
		
		myfile << "tcost = rbind(" << endl;
		for(int i = 0; i < n-1; i++){
			myfile << "c(";
			for(int j = 0; j < n-1; j++){
				myfile << distancematrix[i][j] << ", ";
			}
			myfile << distancematrix[i][n-1] << " ), " << endl;
		}
		myfile << "c(";
		for(int j = 0; j < n-1; j++){
			myfile << distancematrix[n-1][j] << ", ";
		}
		myfile << distancematrix[n-1][n-1] << " ) " << endl;
		myfile << ") " << endl;

		myfile.close();
};


// To shuffle an array a of n elements (indices 0..n-1):

inline void permutation1(vector<int> & a,const int n){
  for(int i=0;i<n;i++){
     int j = i + rand() % (n-i);
     int v = a[i];
     a[i] = a[j];
	 a[j] = v;
  }
}

inline void permutation2(vector<int> & a,const int n){
  for(int i=n;i>0;){
     int j = rand() % i; 
	 i--;
	 int v = a[i]; 
     a[i] = a[j];
	 a[j] = v;
  }
}

inline void permutation3(vector<int> & a,const int n){
  for(int i=n;i>0;){
     int j = rand() % i; 
	 i--;
	 if(i != j){
		 int v = a[i]; 
		 a[i] = a[j];
		 a[j] = v;
	 }
  }
}

void testshuffle(const int n = 200,const int N = 50000){
	
	vector<vector<int>> h(n);
	for(int i = 0; i < n; i++){
	   h[i] = vector<int>(n,0);
	}
	vector<vector<int>> f(n);
	for(int i = 0; i < n; i++){
	   f[i] = vector<int>(n,0);
	}
	vector<vector<int>> g(n);
	for(int i = 0; i < n; i++){
	   g[i] = vector<int>(n,0);
	}

	vector<int> a;
	for(int i = 0; i < n; i++){
		a.push_back(i+1);
	}
	double start = std::clock();
	for(int i = 1; i < N; i++){
		permutation2(a,n);
		for(int i = 0; i < n; i++){
			f[i][a[i]-1]++;
		}
	}
    double end = std::clock();
    std::cout << "Elapsed time: " << (end - start) / CLOCKS_PER_SEC << "\n";
	
	start = std::clock();
	for(int i = 1; i < N; i++){
		permutation3(a,n);
		for(int i = 0; i < n; i++){
			h[i][a[i]-1]++;
		}
	}
    end = std::clock();
    std::cout << "Elapsed time: " << (end - start) / CLOCKS_PER_SEC << "\n";


	start = std::clock();
	for(int i = 1; i < N; i++){
		random_shuffle(a.begin(),a.end());
		for(int i = 0; i < n; i++){
			g[i][a[i]-1]++;
		}
	}
	    
	end = std::clock();
    std::cout << "Elapsed time: " << (end - start) / CLOCKS_PER_SEC << "\n";
	
	

	int myints[] = {1,2,3};

  std::sort (myints,myints+3);

  std::cout << "The 3! possible permutations with 3 elements:\n";
  do {
    std::cout << myints[0] << ' ' << myints[1] << ' ' << myints[2] << '\n';
  } while ( std::next_permutation(myints,myints+3) );

  std::cout << "After loop: " << myints[0] << ' ' << myints[1] << ' ' << myints[2] << '\n';

}

void testmultimap()
	{
  // Compare (<) function not required since it is built into string class.
  // else declaration would comparison function in multimap definition.
  // i.e. multimap<string, int, compare> m;

  multimap<string, int> m;

  m.insert(pair<string, int>("a", 1));
  m.insert(pair<string, int>("c", 2));
  m.insert(pair<string, int>("b", 3));
  m.insert(pair<string, int>("b", 4));
  m.insert(pair<string, int>("a", 5));
  m.insert(pair<string, int>("b", 6));

  cout << "Number of elements with key a: " << m.count("a") << endl;
  cout << "Number of elements with key b: " << m.count("b") << endl;
  cout << "Number of elements with key c: " << m.count("c") << endl;

  cout << "Elements in m: " << endl;
  for (multimap<string, int>::iterator it = m.begin();
       it != m.end();
       ++it)
  {
     cout << "  [" << (*it).first << ", " << (*it).second << "]" << endl;
  }

   pair<multimap<string, int>::iterator, multimap<string, int>::iterator> ppp;

   // equal_range(b) returns pair<iterator,iterator> representing the range
   // of element with key b
   ppp = m.equal_range("b");

   // Loop through range of maps of key "b"
   cout << endl << "Range of \"b\" elements:" << endl;
   for (multimap<string, int>::iterator it2 = ppp.first;
       it2 != ppp.second;
       ++it2)
   {
       cout << "  [" << (*it2).first << ", " << (*it2).second << "]" << endl;
   }

   //cout << m[2] << endl;
   multimap<string, int>::iterator it2 = m.begin();
   advance(it2,2);
   cout << (*it2).first << endl;

// Can't do this (??)
//   cout << ppp.first << endl;
//   cout << ppp.second << endl;

   m.clear();

	int aa[5] =  {5,4,3,2,1};
	vector<int> a(aa,aa+sizeof(aa)/sizeof(*aa));
	vector<int>::iterator b = min_element(a.begin()+1,a.begin()+3);
	b = max_element(a.begin()+1,a.begin()+3);
	std::reverse(a.begin()+1,a.begin()+3);
	std::reverse(a.begin()+1,a.begin()+3);

}

int TOPD::checksolution()
{
	    int retval = 0;

		myfloattype d = calclength(heursol.tour);
		if(abs (d - heursol.length) > TOL){
			retval += 1;
			cerr << "tlength " << d << " " << abs (d - heursol.length) << endl;
		}

		d = calclength(heursol.tour);
		if(heursol.length > Tmax + TOL){
			retval += 10;
			cerr << "Tmax " << Tmax << " " << heursol.length << endl;
		}

		myfloattype p = calcprofit();
		if(abs (p - heursol.profit) > TOL){
			retval += 100;
			cerr << "profit " << p << " " << abs (p - heursol.profit) << endl;
		}
		
		std::vector<myfloattype>::iterator mmm = max_element(heursol.load.begin(),heursol.load.end());
		if(abs (*mmm - heursol.maxload) > TOL){
			retval += 1000;
			cerr << "Lmax " << Lmax << " " << abs (*mmm - heursol.maxload) << endl;
		}

		if (checkcalcbufferplus()==0){
			retval += 10000;
			cerr << "bufferplus inconsistent" << endl;
		}

		myfloattype val = heursol.load[heursol.load.size()-1];
		int i;
		for(i =  heursol.bufferminus.size()-1; i > 0; ){
			if(abs( heursol.bufferminus[i] - val) > TOL){
				retval += 100000;
				cerr << "bufferminus inconsistent" << endl;
				break;
			}
			--i;
			val = min(heursol.load[i],val);
		}
		if(abs( heursol.bufferminus[i]  - val) > TOL){
			retval += 100000;
			cerr << "bufferminus inconsistent" << endl;
		}

		val = heursol.quantity[0];
		for(i = 0; i < heursol.load.size();++i ){
			val += heursol.quantity[i];
			if(abs( heursol.load[i] - val) > TOL){
				retval += 1000000;
				cerr << "load inconsistent with quantity" << endl;
				break;
			}
		}

		val = heursol.quantity[0];
		for(i = 0; i < heursol.load.size(); ++i){
			if( heursol.load[i] < - TOL){
				retval += 10000000;
				cerr << "load is negative" << endl;
				break;
			}
		}

		set <myinttype> selected;
		for(i = 0; i < heursol.tour.size(); ++i){
			selected.insert(heursol.tour[i]);
			if( selected.size() < i+1){
				retval += 100000000;
				cerr << "heursol.tour contains redundant entries" << endl;
				break;
			}
			if(balance[heursol.tour[i]]>0){
				if(heuristicinfo.type[heursol.tour[i]]!=2){
					retval += 1000000000;
					cerr << "heuristicinfo.type is inconsitent with balance or tour (!=2)" << endl;
					break;
				}
			}
			if(balance[heursol.tour[i]]<0){
				if(heuristicinfo.type[heursol.tour[i]]!=-2){
					retval += 1000000000;
					cerr << "heuristicinfo.type is inconsitent with balance or tour (!=-2)" << endl;
					break;
				}
			}
		}

		val = 0;
		for(i = 0; i < heursol.tour.size(); ++i){
			if( balance[heursol.tour[i]] > 0){
				val += heuristicinfo.heurval[heursol.tour[i]];
			}
		}
		if(abs(val - heuristicinfo.pickup.selected.sum)> 100*TOL){
			retval += 10000000000;
			cerr << "heuristicinfo.pickup.selected.sum is inconsitent heuristicinfo.heurval" << endl;
		}

		for(i = 0; i < heursol.quantity.size(); ++i){
			if( abs(heursol.quantity[i]) > abs(balance[heursol.tour[i]]) + TOL){
				retval += 1000000;
				cerr << "quanitiy too large" << endl;
				break;
			}
		}

		return retval;
}


int TOPD::checksolution(sol & solution)
{
	    int retval = 0;

		myfloattype d = calclength(solution.tour);
		if(abs (d - solution.length) > TOL){
			retval += 1;
			cerr << "tlength " << d << " " << abs (d - solution.length) << endl;
		}

		d = calclength(solution.tour);
		if(solution.length > Tmax + TOL){
			retval += 10;
			cerr << "Tmax " << Tmax << " " << solution.length << endl;
		}

		myfloattype p = calcprofit(solution);
		if(abs (p - solution.profit) > TOL){
			retval += 100;
			cerr << "profit " << p << " " << abs (p - solution.profit) << endl;
		}
		
		std::vector<myfloattype>::iterator mmm = max_element(solution.load.begin(),solution.load.end());
		if(abs (*mmm - solution.maxload) > TOL){
			retval += 1000;
			cerr << "Lmax " << Lmax << " " << abs (*mmm - solution.maxload) << endl;
		}

		if (checkcalcbufferplus()==0){
			retval += 10000;
			cerr << "bufferplus inconsistent" << endl;
		}

		myfloattype val = solution.load[solution.load.size()-1];
		int i;
		for(i =  solution.bufferminus.size()-1; i > 0; ){
			if(abs( solution.bufferminus[i] - val) > TOL){
				retval += 100000;
				cerr << "bufferminus inconsistent" << endl;
				break;
			}
			--i;
			val = min(solution.load[i],val);
		}
		if(abs( solution.bufferminus[i]  - val) > TOL){
			retval += 100000;
			cerr << "bufferminus inconsistent" << endl;
		}

		val = solution.quantity[0];
		for(i = 0; i < solution.load.size(); ++i ){
			val += solution.quantity[i];
			if(abs( solution.load[i] - val) > TOL){
				retval += 1000000;
				cerr << "load inconsistent with quantity" << endl;
				break;
			}
		}

		val = solution.quantity[0];
		for(i = 0; i < solution.load.size(); ++i){
			if( solution.load[i] < - TOL){
				retval += 10000000;
				cerr << "load is negative" << endl;
				break;
			}
		}

		set <myinttype> selected;
		for(i = 0; i < solution.tour.size(); ++i){
			selected.insert(solution.tour[i]);
			if( selected.size() < i+1){
				retval += 100000000;
				cerr << "solution.tour contains redundant entries" << endl;
				break;
			}
			if(balance[solution.tour[i]]>0){
				if(heuristicinfo.type[solution.tour[i]]!=2){
					retval += 1000000000;
					cerr << "heuristicinfo.type is inconsitent with balance or tour (!=2)" << endl;
					break;
				}
			}
			if(balance[solution.tour[i]]<0){
				if(heuristicinfo.type[solution.tour[i]]!=-2){
					retval += 1000000000;
					cerr << "heuristicinfo.type is inconsitent with balance or tour (!=-2)" << endl;
					break;
				}
			}
		}

		val = 0;
		for(i = 0; i < solution.tour.size(); ++i){
			if( balance[solution.tour[i]] > 0){
				val += heuristicinfo.heurval[solution.tour[i]];
			}
		}
		if(abs(val - heuristicinfo.pickup.selected.sum)> 100*TOL){
			retval += 10000000000;
			cerr << "heuristicinfo.pickup.selected.sum is inconsitent heuristicinfo.heurval" << endl;
		}

		for(i = 0; i < solution.quantity.size(); ++i){
			if( abs(solution.quantity[i]) > abs(balance[solution.tour[i]]) + TOL){
				retval += 100000000000;
				cerr << "quanitiy too large" << endl;
				break;
			}
		}

		return retval;
}

void TOPD::updatesol()
{
		// intensity, quantity and tour are known
		// this update changes the other compentens of the class sol accordingly. 
		const int N = heursol.load.size();
		heursol.maxload = heursol.load[N-1];
		heursol.bufferplus[N-1] = Lmax - heursol.load[N-1];
		heursol.bufferminus[N-1] = heursol.load[N-1];
		for(int i = N-2; i >= 0; --i){
			heursol.maxload = max(heursol.maxload, heursol.load[i]);
			heursol.bufferplus[i] = min(heursol.bufferplus[i+1], Lmax - heursol.load[i]);
			heursol.bufferminus[i] = min(heursol.bufferminus[i+1], heursol.load[i]);
		}
}

void TOPD::updatesol(sol & solution)
{
		// quantity and tour are known
		// this update changes the other compentens of the class sol accordingly. 

		const int N = solution.load.size();
		solution.load[0] = solution.quantity[0];
		solution.profit = - solution.length * TourCostFactor;
		for(int i = 1; i < N; ++i){
			solution.load[i] = solution.load[i-1]  + solution.quantity[i];
			if(abs(balance[solution.tour[i]])>TOL){
				solution.intensity[i] = solution.quantity[i]/balance[solution.tour[i]];
			} else {
				solution.intensity[i] = 0;
			}
			solution.profit += solution.quantity[i]*profit[solution.tour[i]];
		}

		solution.maxload = solution.load[N-1];
		solution.bufferplus[N-1] = Lmax - solution.load[N-1];
		solution.bufferminus[N-1] = solution.load[N-1];
		for(int i = N-2; i >= 0; --i){
			solution.maxload = max(solution.maxload, solution.load[i]);
			solution.bufferplus[i] = min(solution.bufferplus[i+1], Lmax - solution.load[i]);
			solution.bufferminus[i] = min(solution.bufferminus[i+1], solution.load[i]);
		}
}



int TOPD::repair(sol & solution)
{
	int retval = 1;
	if(solution.maxload >= Lmax ){
		// times Lmax/solution.maxload
		transform(solution.load.begin(),solution.load.end(), solution.load.begin(), bind1st(std::multiplies<myfloattype>(),Lmax/solution.maxload)); 
	}	

	return retval;
}

#ifdef  CPLEXAVAILABLE
myfloattype TOPD::optimizequantities(sol & solution)
{
	try {

		IloEnv envNFP;
		IloModel model(envNFP);

		IloInt	Tlength;
	 
		Tlength = solution.tour.size(); 

		IloNumVarArray intensity(envNFP, Tlength, 0, 1);
		IloNumVarArray load(envNFP, Tlength, 0, Lmax);

		//pick up and delivery 
		for(int i = 1; i < Tlength; i++){
			model.add( load[i] == load[i-1] + ILObalance[solution.tour[i]]*intensity[i] );
		}
		
		// initial load is zero 
		model.add( load[0] == 0 );
	  
		IloNumArray profit(envNFP);
		for(int i = 0; i < Tlength; i++) {
			profit.add( ILOprofit[solution.tour[i]] * ILObalance[solution.tour[i]]);
		}

		// objective 
		IloExpr obj = IloScalProd(profit, intensity);
	 
		model.add(IloMaximize(envNFP, obj));
		obj.end();

		IloCplex cplex(envNFP);
		cplex.setOut(envNFP.getNullStream());

		// cplex.setParam( IloCplex::FPHeur,2); 
		
		cplex.extract(model);

		// cplex.setParam( IloCplex::TiLim,1200 ); 
		// FPHeur = 2
		cplex.setParam( IloCplex::RootAlg, IloCplex::Network);

		cplex.setParam( IloCplex::TiLim,ILOMaxCompTime);
		cplex.solve();

		IloNum tolerance = cplex.getParam(IloCplex::EpInt);
		
		cplex.out() << "Optimal value: " << cplex.getObjValue() - TourCostFactor*solution.length  << endl;

		//cout << "Optimal value: \t" << cplex.getObjValue() - TourCostFactor*solution.length << endl;

		// write model
		// cplex.exportModel ("lpex1.lp");

		//ILOsol.profit =  cplex.getObjValue();
		//ILOsol.upperbound = cplex.getBestObjValue();

		sol mcfsol = heursol;
		vector<myfloattype> ld;
		for(int i=0; i < Tlength; ++i) {
			//IloNum val = cplex.getValue(intensity[i]);
			//if (val >= 1 - tolerance) {
				ld.push_back(cplex.getValue(load[i]));
				mcfsol.load[i] = cplex.getValue(load[i]);
			//}
		}


		vector<myfloattype> inten;
		inten.push_back(0);
		for(int i=1; i < Tlength-1; ++i) {
			inten.push_back(cplex.getValue(intensity[i]));
			IloNum x1 = cplex.getValue(intensity[i]);
			mcfsol.intensity[i] = cplex.getValue(intensity[i]);
			mcfsol.quantity[i] = cplex.getValue(intensity[i])*balance[mcfsol.tour[i]];

			// cout << cplex.getValue(intensity[i]);
		}
		inten.push_back(0);

		myfloattype a = cplex.getObjValue();
		
		updatesol(mcfsol);
		
		for(myinttype i = 1; i < mcfsol.tour.size()-1;){
			 if(mcfsol.intensity[i] < TOL){
			    deletenode(mcfsol,i-1);
			 } else{
				 ++i;
			 }
		}
		 
		updatesol(mcfsol);

		char strbuffer[50];
		sprintf_s (strbuffer, "newsolcplex-%f-",mcfsol.profit);
		Rprintsol(rexe,rpath,filename2,strbuffer,mcfsol);
		
		envNFP.end();
		return a;

	}
	catch(IloException& e) {
		cerr  << " ERROR: " << e << endl;   
		return 0;
	}
	catch(...) {
		cerr  << " ERROR" << endl;   
		return 0;
	}
}
#endif

myfloattype TOPD::optimizequantities2(sol & solution)
{
	
	myinttype m = solution.tour.size(); 
	vector<myfloattype> profitplus(m);
	vector<myinttype> profitplusindex(m);
	vector<myfloattype> profitminus(m);
	vector<myinttype> profitminusindex(m);
	
	vector<myfloattype> bufferminus(m);
	vector<myfloattype> bufferplus(m);

	profitminus[0] = 1e10;
	profitminusindex[0] = 0;
	for(myinttype i = 1; i < m; ++i){
		if(balance[solution.tour[i]]>0 && -profit[solution.tour[i]] < profitminus[i-1]){
			profitminus[i] = -profit[solution.tour[i]];
			profitminusindex[i] = i; 
		} else {
			profitminus[i] = profitminus[i-1];
			profitminusindex[i] = profitminusindex[i-1];
		}
	}

	profitplus[m-1] = 0;
	profitplusindex[m-1] = m-1;
	for(myinttype i = m-2; i >= 0; --i){
		if(balance[solution.tour[i]]<0 && -profit[solution.tour[i]] + TOL >= profitplus[i+1]){
			profitplus[i] = -profit[solution.tour[i]];
			profitplusindex[i] = i;
		} else {
			profitplus[i] = profitplus[i+1];
			profitplusindex[i] = profitplusindex[i+1];
		}
	}

	vector<myfloattype> result(m);
	std::transform(profitplus.begin(), profitplus.end(), profitminus.begin(), result.begin(),minus<myfloattype>());

	sol newsolution;
	newsolution = heursol;
	newsolution.intensity = vector<myfloattype>(m,0);
	newsolution.quantity = vector<myfloattype>(m,0);
	newsolution.profit = 0;
	updatesol(newsolution);

	vector<myfloattype>::iterator minel;
	minel = max_element(result.begin(),result.end());
	myfloattype index = distance(result.begin(),minel);

	if(*minel>TOL){
		newsolution.quantity[profitminusindex[index]] = min(balance[solution.tour[profitminusindex[index]]],-balance[solution.tour[profitplusindex[index]]]);
		newsolution.quantity[profitplusindex[index]] = -newsolution.quantity[profitminusindex[index]];
		updatesol(newsolution);
		assert(checksolution(newsolution)==0);
	} else {
		return 0;
	}

	// add most profitable pairs (+,-)
	while(1==1){
		profitminus[0] = 1e10;
		profitminusindex[0] = 0;
		for(myinttype i = 1; i < m; ++i){
			if(solution.bufferplus[i] > TOL){
				if(balance[solution.tour[i]]>0 && -profit[solution.tour[i]] < profitminus[i-1] && newsolution.quantity[i] + TOL < balance[solution.tour[i]]   ){
					profitminus[i] = -profit[solution.tour[i]];
					profitminusindex[i] = i; 
				} else {
					profitminus[i] = profitminus[i-1];
					profitminusindex[i] = profitminusindex[i-1];
				}
			} else {
				profitminus[i] = 1e10;
				profitminusindex[i] = i; 
			}
		}

		profitplus[m-1] = 0;
		profitplusindex[m-1] = m-1;
		for(myinttype i = m-2; i >= 0; --i){
			if(solution.bufferplus[i] > TOL){
				if(balance[solution.tour[i]]<0 && -profit[solution.tour[i]] + TOL >= profitplus[i+1] && -newsolution.quantity[i] + TOL < -balance[solution.tour[i]] ){
					profitplus[i] = -profit[solution.tour[i]];
					profitplusindex[i] = i;
				} else {
					profitplus[i] = profitplus[i+1];
					profitplusindex[i] = profitplusindex[i+1];
				}
			} else {
				profitplus[i] = 0;
				profitplusindex[i] = i; 
			}
		}

		std::transform(profitplus.begin(), profitplus.end(), profitminus.begin(), result.begin(),minus<myfloattype>());
		//vector<myfloattype>::iterator minel;
		minel = max_element(result.begin(),result.end());
		index = distance(result.begin(),minel);

		if(*minel>TOL){
			const myfloattype delta = min(balance[solution.tour[profitminusindex[index]]]-newsolution.quantity[profitminusindex[index]],-balance[solution.tour[profitplusindex[index]]]+newsolution.quantity[profitplusindex[index]]);
			newsolution.quantity[profitminusindex[index]] += delta;
			newsolution.quantity[profitplusindex[index]] -= delta;
			updatesol(newsolution);
			assert(checksolution(newsolution)==0);
		} else {
			break;
		}
	}

	// check for inverse pairs (- +) intensity < 1 forward
	myinttype candi = 1;
	while(candi<m-1){
		if(balance[newsolution.tour[candi]]<0 && newsolution.intensity[candi] + TOL < 1){
			int test = 0;
			myfloattype minload = newsolution.load[candi];
			myfloattype maxprft = 0;
			myinttype candj = 0;
			myfloattype candminload = 0;
			for(myinttype j = candi+1; j < m-1; ++j){
				minload = min(minload,newsolution.load[j]);
				if( balance[newsolution.tour[j]]>0 && newsolution.intensity[j] + TOL < 1 && minload > TOL){
					myfloattype prft = - profit[newsolution.tour[candi]] + profit[newsolution.tour[j]];
					if( prft + TOL >= maxprft ){
						maxprft = prft;
						candj = j;
						candminload = minload;
						test = 1;
					}
				}
			}
			if(test == 1){
				myfloattype delta = min(balance[newsolution.tour[candj]]-newsolution.quantity[candj],-balance[newsolution.tour[candi]]+newsolution.quantity[candi]);
				delta = min(candminload,delta);
				newsolution.quantity[candj] += delta;
				newsolution.quantity[candi] -= delta;
				updatesol(newsolution);
				assert(checksolution(newsolution)==0);
			} else {
				++candi;
			}
		} else {
			++candi;
		}
	}

	// check for inverse pairs (- +) intensity < 1 backwards
	candi = m-2;
	while(candi>0){
		if(balance[newsolution.tour[candi]]>0 && newsolution.intensity[candi] + TOL < 1){
			int test = 0;
			myfloattype maxload = 0;
			myfloattype maxprft = 0;
			myinttype candj = 0;
			myfloattype candminload = 0;
			
			for(myinttype j = candi+1; j < m-1; ++j){
				maxload = max(maxload,newsolution.load[j]);
				if(maxload + TOL < Lmax) {
					break;
				}
				if( balance[newsolution.tour[j]] < 0){
					myfloattype prft = - profit[newsolution.tour[candi]] + profit[newsolution.tour[j]];
					if( prft + TOL >= maxprft ){
						maxprft = prft;
						candj = j;
						candminload = Lmax - maxload;
						test = 1;
					}
				}
			}

			myfloattype minload = Lmax;
			for(myinttype j = candi -1 ; j > 0; --j){
				minload = min(minload,newsolution.load[j]);
				if(minload  < TOL) {
					break;
				}
				if( balance[newsolution.tour[j]]<0 ){
					myfloattype prft = - profit[newsolution.tour[candi]] + profit[newsolution.tour[j]];
					if( prft + TOL >= maxprft ){
						maxprft = prft;
						candj = j;
						candminload = minload;
						test = 1;
					}
				}
			}

			if(test == 1){
				myfloattype delta = min(newsolution.quantity[candj],newsolution.quantity[candi]);
				delta = min(candminload,delta);
				newsolution.quantity[candj] -= delta;
				newsolution.quantity[candi] += delta;
				updatesol(newsolution);
				assert(checksolution(newsolution)==0);
			} else {
				--candi;
			}
		} else {
			--candi;
		}
	}

	// check if there is a better delivery customer (-,-)
	candi = 1;
	while(candi<m-1){
		if(balance[newsolution.tour[candi]] < 0 && newsolution.intensity[candi] + TOL < 1){
			int test = 0;
			myfloattype minload = newsolution.load[candi];
			myfloattype maxprft = 0;
			myinttype candj = 0;
			myfloattype candminload = 0;
			for(myinttype j = candi+1; j < m-1; ++j){
				minload = min(minload,newsolution.load[j-1]);
				if(minload  < TOL) {
					break;
				}
				if( newsolution.quantity[j]< -TOL ){
					myfloattype prft = - profit[newsolution.tour[candi]] + profit[newsolution.tour[j]];
					if( prft + TOL >= maxprft ){
						maxprft = prft;
						candj = j;
						candminload = minload;
						test = 1;
					}
				}
			}

			myfloattype maxload = 0;
			
			for(myinttype j = candi-1; j > 0; --j){
				maxload = max(maxload,newsolution.load[j]);
				if(maxload + TOL > Lmax ) {
					break;
				}
				if(  newsolution.quantity[j] < -TOL ){
					myfloattype prft = - profit[newsolution.tour[candi]] + profit[newsolution.tour[j]];
					if( prft + TOL >= maxprft ){
						maxprft = prft;
						candj = j;
						candminload = Lmax - maxload;
						test = 1;
					}
				}
			}

			if(test == 1){
				myfloattype delta = min(-newsolution.quantity[candj],-balance[newsolution.tour[candi]]+newsolution.quantity[candi]);
				delta = min(candminload,delta);
				newsolution.quantity[candj] += delta;
				newsolution.quantity[candi] -= delta;
				updatesol(newsolution);
				assert(checksolution(newsolution)==0);
			} else {
				++candi;
			}
		} else {
			++candi;
		}
	}

	// check if there is a better pickup customer (+,+)
	candi = 1;
	while(candi<m-1){
		if(balance[newsolution.tour[candi]] > 0 && newsolution.intensity[candi] + TOL < 1){
			int test = 0;
			myfloattype maxload = newsolution.load[candi];
			myfloattype maxprft = 0;
			myinttype candj = 0;
			myfloattype candminload = 0;
			for(myinttype j = candi+1; j < m-1; ++j){
				maxload = max( maxload,newsolution.load[j]);
				if(maxload + TOL > Lmax ) {
					break;
				}
				if( newsolution.quantity[j]>TOL ){
					myfloattype prft = profit[newsolution.tour[candi]] - profit[newsolution.tour[j]];
					if( prft + TOL >= maxprft ){
						maxprft = prft;
						candj = j;
						candminload = Lmax - maxload;
						test = 1;
					}
				}
			}

			myfloattype minload = newsolution.load[candi];;
	
			for(myinttype j = candi-1; j > 0; --j){
				minload = min(minload,newsolution.load[j]);
				if(minload < TOL) {
					break;
				}
				if( newsolution.quantity[j]>TOL ){
					myfloattype prft = profit[newsolution.tour[candi]] - profit[newsolution.tour[j]];
					if( prft + TOL >= maxprft ){
						maxprft = prft;
						candj = j;
						candminload = minload;
						test = 1;
					}
				}
			}

			if(test == 1){
				myfloattype delta = min(newsolution.quantity[candj],balance[newsolution.tour[candi]]-newsolution.quantity[candi]);
				delta = min(candminload,delta);
				newsolution.quantity[candj] -= delta;
				newsolution.quantity[candi] += delta;
				updatesol(newsolution);
				assert(checksolution(newsolution)==0);
			} else {
				++candi;
			}
		} else {
			++candi;
		}
	}



	//updatesol(newsolution);
	char strbuffer [50];
	sprintf_s(strbuffer,"heu1-%f-",newsolution.profit + newsolution.length);
	//Rprintsol(rexe,rpath,Filename,strbuffer,newsolution);

	if(newsolution.profit > heursol.profit){
		heursol = newsolution;
		for(myinttype i = 1; i < heursol.tour.size()-1;){
			if(heursol.intensity[i] < TOL){
			deletenode(i-1);
			} else {
				++i;
			}
		}
		char strbuffer2 [50];
		sprintf_s(strbuffer2,"heunew-%f-",newsolution.profit+ newsolution.length);
		//Rprintsol(rexe,rpath,filename2,strbuffer2,newsolution);
		//Rprintsol(rexe,rpath,filename2,strbuffer2,newsolution);
	}

	return newsolution.profit + newsolution.length;
}


 /*   		fmat = matrix(0,3*(length(tournew)-2)+1,2*(length(tournew)-2)+1)
			fmat[1,1]=1
			for(i in 2:(length(tournew)-1)) {fmat[i,i-1]=1;fmat[i,i]=-1;fmat[i,length(tournew)-2+i]=if(instance1[tournew[i],3]>0){1}else{-1} }
			for(i in 2:(2*(length(tournew)-2)+1)) {fmat[length(tournew)-1+i-1,i]=1 }

			fdir = c();
			fdir[1:(length(tournew)-1)] = "=="
			fdir[length(tournew):(3*(length(tournew)-2)+1)] = "<="

			fobj=c(); 
			fobj[(length(tournew)):(2*(length(tournew)-2)+1)] = 
			    instance1[tournew[2:(length(tournew)-1)],4]*(1*(instance1[tournew[2:(length(tournew)-1)],3]>0))-
				instance1[tournew[2:(length(tournew)-1)],4]*(1*(instance1[tournew[2:(length(tournew)-1)],3]<0))
			fobj[is.na(fobj)]=0;
			fobj = fobj * (lmaxold/tmax);
			
			fb=c();
			fb[1:((length(tournew)-2)+1)]=0
			fb[length(tournew):(2*(length(tournew)-2)+1)]=lmax
			fb[(2*(length(tournew)-2)+2):(3*(length(tournew)-2)+1)]=abs(instance1[tournew[2:(length(tournew)-1)],3]);
			
			
			##debug##      print("????????????")
			##debug##      print(fobj)
			##debug##      print(fmat)
			##debug##      print(fdir)
			##debug##      print(fb)
			##debug##      print("????????????")
			
			resopt = Rglpk_solve_LP(fobj, fmat, fdir, rhs=fb, max = TRUE);
			##debug##      cat("/././././././././././././././././\n")
			##debug##      print(resopt)
			##debug##      cat("/././././././././././././././././\n")
			reszeros = 0;
			if(resopt$status == 0){
			   reszeros =  sum(abs(resopt$solution[((length(tournew)-2)+2):((2*(length(tournew)-2)+1))])<1e-10)
			}
			if(resopt$status == 0 & ((profit  < resopt$optimum - ctourfactor*tlength) | (profit == resopt$optimum - ctourfactor*tlength &  reszeros>0) )){
			   ##debug##       cat("/////////////////\n")
			   ##debug##       cat("espilon3",epsilon3,"\n");
			   ##debug##       print(resopt)
			   ##debug##       cat("solution \t",resopt$solution[((length(tournew)-2)+2):((2*(length(tournew)-2)+1))],"\n");
			   ##debug##       cat("max abs bal \t",abs(instance1[tournew[2:(length(tournew)-1)],3]),"\n");
			   ##debug##       cat("bal (old) \t",instance1[tournew,3]*intensity,"\n");
			   ##debug##       cat("intensity (old) \t",intensity,"\n");
			   ##debug##       cat("profitsnew*int(old) \t",cumsum(profitsnew*intensity)- ctourfactor*tlength,"\n");
			   ##debug##       cat("profitsnew*int(old) bd \t",cumsum(instance1[tournew,3]*instance1[tournew,4]*intensity*(lmaxold/tmax))- ctourfactor*tlength,"\n");
			   ##debug##       cat("profitsnew \t",profitsnew,"\n");
			   ##debug##       cat("profitsnew bd \t",instance1[tournew,3]*instance1[tournew,4]*(lmaxold/tmax),"\n");
			   
			   intensity[2:(length(tournew)-1)] = resopt$solution[((length(tournew)-2)+2):((2*(length(tournew)-2)+1))]/abs(instance1[tournew[2:(length(tournew)-1)],3]);
			   
			   ##debug##       cat("intensity \t",intensity,"\n");
			   ##debug##       cat("profit \t",profit,"\n");
			   ##debug##       cat("profitsnew*int \t",cumsum(profitsnew*intensity),"\n")
			   ##debug##       cat("/////////////////\n")
			   
			   profitscumnew = cumsum(profitsnew*intensity);
			   profit = sum(profitsnew*intensity)- ctourfactor*tlength;
			   
			   ##debug##       cat("profit \t",profit,"\n");
			   ##debug##       cat("/////////////////\n")
			   
			   loadingnew = cumsum(instance1[tournew,3]*intensity);
			   maxload = max(loadingnew);*/




// this routine destroys parts of the solution
// return = 1 returns a proven feasible solution -> sould not be followed by a repair
int TOPD::destroy(sol & solution, myinttype type){
	if(type < 5){ // delete randomly parts of the solution (random backwards)
		int delstart = 1 + (rand() % (solution.tour.size()));
		for(int k = solution.tour.size()-2; k >= delstart;--k){ 
			(*solution.topd).deletenode(k-1);
			assert((*solution.topd).checksolution()==0);
		}
		return 1;
	}
	if(type >= 5 && type < 80){ // delete randomly parts of the solution (random 10%)
		myfloattype q = 0.1; // delete approx. 10% 
		myinttype k = 0; // number of deleted pickup customers
		myinttype i = 0;
		//while ( i < (*solution.topd).heuristicinfo.pickup.selected.index.size()){
		while ( i < solution.tour.size()){
			const myinttype node = solution.tour[i];
			if((*solution.topd).heuristicinfo.type[node]!=2){
				++i;
				continue;
			}
			const myfloattype rnd = ((float)rand()/RAND_MAX);
			if( rnd <= q  * (1 - (*solution.topd).heuristicinfo.heurval[node]/(*solution.topd).heuristicinfo.pickup.selected.sum) ){
				(*solution.topd).deletenode(i-1);
				//assert((*solution.topd).checksolution()==0);
				++k;	
				if(k >= q*solution.tour.size()){
					break;
				}
			} else {
				++i;
			}
		}

		i = 0;
		while ( i < solution.tour.size()){
			const myinttype node = solution.tour[i];
			if((*solution.topd).heuristicinfo.type[node]!=-2){
				++i;
				continue;
			}
			if( solution.bufferminus[i] + TOL <= 0 ){
				if( abs(solution.bufferminus[i]-(*solution.topd).balance[node]) <= TOL){
					(*solution.topd).deletenode(i-1);
					// assert((*solution.topd).checksolution()==0);
				} else {
					(*solution.topd).deletenode(i-1);
					// assert((*solution.topd).checksolution()==0);
					// WIP - hier könnte man auch "partiell" rückeinfügen.
				}
				++k;
			} else {
				++i;
			}
		}
	}
	if(type >= 80 && type < 100){ // delete randomly parts of the solution (random 10%)
		myfloattype q = 0.1; // delete approx. 10% 
		myinttype k = 0; // number of deleted pickup customers
		myinttype i = 0;
		//while ( i < (*solution.topd).heuristicinfo.pickup.selected.index.size()){
		while ( i < solution.tour.size()){
			const myinttype node = solution.tour[i];
			if((*solution.topd).heuristicinfo.type[node]!=2){
				++i;
				continue;
			}
			const myfloattype rnd = ((float)rand()/RAND_MAX);
			if( rnd <= q  ){
				(*solution.topd).deletenode(i-1);
				//assert((*solution.topd).checksolution()==0);
				++k;	
				if(k >= q*solution.tour.size()){
					break;
				}
			} else {
				++i;
			}
		}

		i = 0;
		while ( i < solution.tour.size()){
			const myinttype node = solution.tour[i];
			if((*solution.topd).heuristicinfo.type[node]!=-2){
				++i;
				continue;
			}
			if( solution.bufferminus[i] + TOL <= 0 ){
				if( abs(solution.bufferminus[i]-(*solution.topd).balance[node]) <= TOL){
					(*solution.topd).deletenode(i-1);
					//assert((*solution.topd).checksolution()==0);
				} else {
					(*solution.topd).deletenode(i-1);
					//assert((*solution.topd).checksolution()==0);
					// WIP - hier könnte man auch "partiell" rückeinfügen.
				}
				++k;
			} else {
				++i;
			}
		}
	}
	if(type >= 100 && type < 200){ // delete randomly parts of the solution (random 10%)
		myfloattype q = 0.1; // delete approx. 10% 
		myinttype k = 0; // number of deleted pickup customers
		myinttype i = 0;
		//while ( i < (*solution.topd).heuristicinfo.pickup.selected.index.size()){
		while ( i < solution.tour.size()){
			const myinttype node = solution.tour[i];
			if((*solution.topd).heuristicinfo.type[node]!=2){
				++i;
				continue;
			}
			const myfloattype rnd = ((float)rand()/RAND_MAX);
			if( rnd <= q  * (1 - (*solution.topd).heuristicinfo.heurval[node]/(*solution.topd).heuristicinfo.pickup.selected.sum) ){
				const myinttype node = solution.tour[i];
				(*solution.topd).deletenode(i-1);
				//assert((*solution.topd).checksolution()==0);
				++k;	
				if(k == 1){
					updateheuristicinfo(2,node);
				}
			} else {
				++i;
			}
		}

		i = 0;
		while ( i < solution.tour.size()){
			const myinttype node = solution.tour[i];
			if((*solution.topd).heuristicinfo.type[node]!=-2){
				++i;
				continue;
			}
			if( solution.bufferminus[i] + TOL <= 0 ){
				if( abs(solution.bufferminus[i]-(*solution.topd).balance[node]) <= TOL){
					(*solution.topd).deletenode(i-1);
					//assert((*solution.topd).checksolution()==0);
				} else {
					(*solution.topd).deletenode(i-1);
					//assert((*solution.topd).checksolution()==0);
					// WIP - hier könnte man auch "partiell" rückeinfügen.
				}
				++k;
			} else {
				++i;
			}
		}

	}
}

////////////////////////////
////////////////////////////

class TWOEXCHANGE {
public:
	TOPD *topd;
	sol *solution;
	int type;
	myinttype N; // #selected	
	myinttype M; // #selected	
	set<myinttype>::iterator M1; // #available
	set<myinttype>::iterator M2; // #available
	myfloattype delta;
	pair <myinttype,set<myinttype>::iterator > actualmove; // first:selected - second:available
	int nextmove();
	int isfeasible();
	int loop();
	int LS();
	TWOEXCHANGE(TOPD *xxxx);
	friend class TOPD;
};

TWOEXCHANGE::TWOEXCHANGE(TOPD *xxxx){
	delta = 0;
	topd = xxxx;
	solution = & ((*xxxx).heursol);
	N = (*xxxx).heursol.tour.size()-2;
	M1 = (*topd).heuristicinfo.pickup.available.index.end();
	M2 = (*topd).heuristicinfo.delivery.available.index.end();
	M1--;
	M2--;
	type = (*xxxx).heuristicinfo.type[(*xxxx).heursol.tour[1]];
	actualmove = pair <myinttype,set<myinttype>::iterator>(1,(*xxxx).heuristicinfo.pickup.available.index.begin());
}

int TWOEXCHANGE::nextmove(){
	delta = 0;
	if(type == 2){
		if( actualmove.second != M1){
			++actualmove.second;
		} else {
			if ( actualmove.first != N) {
				++actualmove.first;
				type = (*topd).heuristicinfo.type[(*topd).heursol.tour[actualmove.first]];
				if(type == 2){
					actualmove.second = (*topd).heuristicinfo.pickup.available.index.begin();
				} else {
					actualmove.second = (*topd).heuristicinfo.delivery.available.index.begin();
				}
			} else {
				return 0;
			}
		} 
	} else { // type ==-2
		if( actualmove.second != M2){
			++actualmove.second;
		} else {
			if ( actualmove.first != N) {
				++actualmove.first;
				type = (*topd).heuristicinfo.type[(*topd).heursol.tour[actualmove.first]];
				if(type == 2){
					actualmove.second = (*topd).heuristicinfo.pickup.available.index.begin();
				} else {
					actualmove.second = (*topd).heuristicinfo.delivery.available.index.begin();
				}
			} else {
				return 0;
			}
		} 
	}
	return 1;
}

/*
int TWOEXCHANGE::isfeasible(){
	const int i = actualmove.first;
	const int j = actualmove.second;
	myfloattype deltaload = (*topd).heursol.quantity[i]-(*topd).heursol.quantity[j-1];
	
	std::vector<myfloattype>::iterator result;
	result =  max_element((*topd).heursol.load.begin()+i, (*topd).heursol.load.begin()+j);
	const myfloattype lmin  = (*topd).heursol.load[i-1] + (*topd).heursol.load[j-1] - *result;
	result =  min_element((*topd).heursol.load.begin()+i, (*topd).heursol.load.begin()+j);
	const myfloattype lmax  = (*topd).heursol.load[i-1] + (*topd).heursol.load[j-1] - *result;

	if (  lmin + TOL >= 0 ){					// bufferminus feasible
		if (    lmax  <= (*topd).Lmax + TOL){   // bufferplus feasible
			return 1;
		} else {
			return 0;
		}
	} else {
		return 0;
	}
}
*/

int TWOEXCHANGE::loop(){
	cout << "----------"<< endl; 
	if (N==3) return 0; 
	int start = 1;
	while(start == 1 || nextmove()==1){
		if(type==2){
			cout << "P: " <<(*topd).heursol.tour[1 + actualmove.first ] << " vs "<<  *(actualmove.second) << endl;
				//pickup.available.index[actualmove.second] << endl; // (*topd).heursol.tour[1 + actualmove.first ] << " vs "<< (*topd).heuristicinfo.pickup.available.index[actualmove.second] << endl;
		} else {
			cout << "D: " <<(*topd).heursol.tour[1 + actualmove.first ] << " vs "<<  *(actualmove.second) << endl;
			//(*topd).heursol.tour[1 +actualmove.first] << " vs "<< (*topd).heuristicinfo.delivery.available.index[actualmove.second] << endl;
		}
		start = 0;
	}
	return 1;
}

/////////////////////
/////////////////////

class TWOOPT {
public:
	TOPD *topd;
	sol *solution;
	myinttype N;
	myfloattype delta;
	pair <myinttype,myinttype> actualmove; // first < second
	int nextmove();
	int isfeasible();
	int loop();
	int LS();
	TWOOPT(TOPD *xxxx);
	friend class TOPD;
};

TWOOPT::TWOOPT(TOPD *xxxx){
	delta = 0;
	topd = xxxx;
	solution = & ((*xxxx).heursol);
	N = (*solution).tour.size();
	actualmove = pair <myinttype,myinttype>(2,N);
}

int TWOOPT::nextmove(){
	delta = 0;
	if(actualmove.first!=1){
		--actualmove.first;
		--actualmove.second;
	} else {
		actualmove.first = N-actualmove.second+1;
		actualmove.second = N-1;
		if(actualmove.first==N-2){
			actualmove = pair <myinttype,myinttype>(2,N);
			return 0;
		}
	}
	return 1;
}

int TWOOPT::isfeasible(){
	const int i = actualmove.first;
	const int j = actualmove.second;
	myfloattype deltaload = (*topd).heursol.quantity[i]-(*topd).heursol.quantity[j-1];
	
	std::vector<myfloattype>::iterator result;
	result =  max_element((*topd).heursol.load.begin()+i, (*topd).heursol.load.begin()+j);
	const myfloattype lmin  = (*topd).heursol.load[i-1] + (*topd).heursol.load[j-1] - *result;
	result =  min_element((*topd).heursol.load.begin()+i, (*topd).heursol.load.begin()+j);
	const myfloattype lmax  = (*topd).heursol.load[i-1] + (*topd).heursol.load[j-1] - *result;

	if (  lmin + TOL >= 0 ){					// bufferminus feasible
		if (    lmax  <= (*topd).Lmax + TOL){   // bufferplus feasible
			return 1;
		} else {
			return 0;
		}
	} else {
		return 0;
	}
}

int TWOOPT::loop(){
	cout << "----------"<< endl; 
	if (N==3) return 0; 
	while(nextmove()==1){
		cout << "(" << actualmove.first << "," << actualmove.second << ")" << endl;
	}
	return 1;
}

int TWOOPT::LS(){
	//cout << "----------"<< endl; 
	if (N==3) return 0; 
	while(nextmove()==1){
		//cout << "(" << actualmove.first << "," << actualmove.second << ")" << endl;
		//cout << "(" << actualmove.first << "," << actualmove.second << ")" << endl; 
		const int i = actualmove.first;
		const int j = actualmove.second;
		myfloattype deltadist = - (*topd).distancematrix[(*solution).tour[i-1]][(*solution).tour[i]] - (*topd).distancematrix[(*solution).tour[j-1]][(*solution).tour[j]];
		deltadist += (*topd).distancematrix[(*solution).tour[i-1]][(*solution).tour[j-1]] + (*topd).distancematrix[(*solution).tour[i]][(*solution).tour[j]];
		if(deltadist <= -TOL){
			 //cout << "(" << actualmove.first << "," << actualmove.second << ")" << " "<< deltadist << endl;
			if ( isfeasible()==1 ) {
				// 
				//    I N      P R O G R E S S :::: load anpassen un dann umdrehen 
				// 
				//
				// std::vector<myfloattype>::iterator result;
				// result = max_element((*topd).heursol.load.begin()+i, (*topd).heursol.load.begin()+j);
				// if(*result + TOL >= (*topd).heursol.maxload ){ // redundanz bei update
				//	(*topd).heursol.maxload = *result;
				// }


				// myfloattype deltaload = (*topd).heursol.quantity[i]-(*topd).heursol.quantity[j-1];
				
				reverse((*topd).heursol.quantity.begin()+i,(*topd).heursol.quantity.begin()+j);
				reverse((*topd).heursol.tour.begin()+i,(*topd).heursol.tour.begin()+j);
				reverse((*topd).heursol.intensity.begin()+i,(*topd).heursol.intensity.begin()+j);				
				for(int k = i ; k < j; k++){
					 (*topd).heursol.load[k] =  (*topd).heursol.load[k-1] +  (*topd).heursol.quantity[k];
				}
				
				//transform((*topd).heursol.load.begin()+i,(*topd).heursol.load.begin()+j-1,(*topd).heursol.load.begin()+i,bind1st(std::plus<myfloattype>(),-(*topd).heursol.load[j-1]));
				//transform((*topd).heursol.load.begin()+i,(*topd).heursol.load.begin()+j-1,(*topd).heursol.load.begin()+i,bind1st(std::multiplies<myfloattype>(),-1));
				//transform((*topd).heursol.load.begin()+i,(*topd).heursol.load.begin()+j-1,(*topd).heursol.load.begin()+i,bind1st(std::plus<myfloattype>(),+(*topd).heursol.load[i]));
				//transform((*topd).heursol.load.begin()+i,(*topd).heursol.load.begin()+j-1,(*topd).heursol.load.begin()+i,bind1st(std::plus<myfloattype>(),-deltaload));
				//reverse((*topd).heursol.load.begin()+i,(*topd).heursol.load.begin()+j-1);


				(*topd).updatesol(); //  hier wird auch das max bereuchnet.
				(*topd).heursol.profit -=  (*topd).TourCostFactor * deltadist;
				(*topd).heursol.length +=  deltadist;
				assert((*topd).checksolution()==0);
			}
		}
	}
	return 1;
}

void reverseheu(TOPD & topd1,int itlim){
	vector < infomove >  vmoves;
	for(int i = 1; i <= itlim; i++){
		// cout << i << ": " << topd1.heursol.tour[0] << "-" << topd1.heursol.tour[1] << "-" << topd1.heursol.tour[2] << endl;
		topd1.updateheuristicinfo();
		topd1.checksolution();
		infomove move_p = topd1.addpickup(topd1.heursol);
		topd1.checksolution();
		// TWOOPT twoopt(&topd1);
		// twoopt.loop();
		//twoopt.LS();
		//topd1.checksolution();
		  // d = topd1.calclength(topd1.heursol.tour);
		if ( move_p.status == 1){
			vmoves.push_back(move_p);
			int kstep = 1;
			infomove move_d = topd1.adddelivery(topd1.heursol);
			topd1.checksolution();
			 // d = topd1.calclength(topd1.heursol.tour);
			int test = 0;
			while( move_d.status == 1){
				++kstep ;
				test = 1;
				vmoves.push_back(move_d);
				if(topd1.elitists.ub < topd1.heursol.profit){
					topd1.elitists.member.push_back(topd1.heursol);
					topd1.elitists.ub = topd1.heursol.profit;
					topd1.checksolution();
					cout << i <<  "\t" << topd1.heursol.profit << endl;
					char buffer [50];
				    sprintf_s (buffer, "heu-%f-",topd1.heursol.profit);
					//topd1.Rprintsol(rexe,rpath,filename,buffer,topd1.heursol);
				}
				move_d = topd1.adddelivery(topd1.heursol);
				topd1.checksolution();
			    // d = topd1.calclength(topd1.heursol.tour);	 
			}
			if (test == 1){
				continue;
			}
			//infomove move = vmoves[vmoves.size()-1];
			//topd1.deletenode(move);
			//vmoves.erase(vmoves.begin() + vmoves.size() - 1);
			//d = topd1.calclength(topd1.heursol.tour);

			int delstart = vmoves.size() - ( rand() % kstep );
			for(int k = vmoves.size()-1; k >= delstart;--k){
				infomove move = vmoves[k];
				topd1.deletenode(move);
				topd1.checksolution();
				//d = topd1.calclength(topd1.heursol.tour);
			}
			vmoves.erase(vmoves.begin() + delstart, vmoves.end());
			continue;
		}
		if(vmoves.size() == 0){
			continue;
		}
		int delstart = rand() % vmoves.size();
		for(int k = vmoves.size()-1; k >= delstart;--k){
			infomove move = vmoves[k];
			topd1.deletenode(move);
			topd1.checksolution();
			//d = topd1.calclength(topd1.heursol.tour);
		}
		vmoves.erase(vmoves.begin() + delstart, vmoves.end());
	}
}

void reverseheu2(TOPD & topd1,int itlim){
	vector < infomove >  vmoves;

	for(int i = 1; i <= itlim; i++){
		//cout << i << endl;
		topd1.updateheuristicinfo();
		infomove move_p = topd1.addpickup(topd1.heursol);
		topd1.checksolution(); 

		if ( move_p.status == 1){
			vmoves.push_back(move_p);
			int kstep = 1;
			infomove move_d = topd1.adddelivery(topd1.heursol);
			topd1.checksolution();
			int test = 0;
			while( move_d.status == 1){
				++kstep ;
				test = 1;
				vmoves.push_back(move_d);
				if(topd1.elitists.ub < topd1.heursol.profit){
					topd1.elitists.member.push_back(topd1.heursol);
					topd1.elitists.ub = topd1.heursol.profit;
					topd1.checksolution();
					cout << i <<  "\t" << topd1.heursol.profit << endl;
					char buffer [50];
				    sprintf_s (buffer, "heu-%f-",topd1.heursol.profit);
					//topd1.Rprintsol(rexe,rpath,filename,buffer,topd1.heursol);
				}
				move_d = topd1.adddelivery(topd1.heursol);
				topd1.checksolution();
			}
			if (test == 1){
				continue;
			}
			//infomove move = vmoves[vmoves.size()-1];
			//topd1.deletenode(move);
			//vmoves.erase(vmoves.begin() + vmoves.size() - 1);

			int delstart = vmoves.size() - ( rand() % kstep );
			for(int k = vmoves.size()-1; k >= delstart;--k){
				infomove move = vmoves[k];
				topd1.deletenode(move.tourpos);
				topd1.checksolution();
			}
			vmoves.erase(vmoves.begin() + delstart, vmoves.end());
			continue;
		}
		if(vmoves.size() == 0){
			continue;
		}
		int delstart = rand() % vmoves.size();
		for(int k = vmoves.size()-1; k >= delstart;--k){
			infomove move = vmoves[k];
			topd1.deletenode(move.tourpos);
			topd1.checksolution();
		}
		vmoves.erase(vmoves.begin() + delstart, vmoves.end());
	}
}


void simpleLNS(TOPD & topd1,int itlim){
	myfloattype rat1 = 0;
	myfloattype rat2 = 0;
	vector < infomove >  vmoves;
	
	for(int i = 1; i <= itlim; i++){
		//cout << i << endl;
		//topd1.updateheuristicinfo();
		infomove move_p = topd1.addpickup(topd1.heursol);
		assert(topd1.checksolution()==0);
		
		if( topd1.heursol.tour.size()>=4){
			TWOOPT twoopt(&topd1);
			//twoopt.loop();
			twoopt.LS();
			assert(topd1.checksolution()==0);

			//myfloattype a = topd1.optimizequantities( topd1.heursol);
			//myfloattype b = topd1.optimizequantities2( topd1.heursol);


			/*
			if(abs(a-b)>TOL && b>TOL){
			//if(a - TOL >  topd1.heursol.profit){
			//if(a + TOL < b && a>0){
			 cout << a << "," << b << endl << "a ";
			 cin >> a;
			//	++rat1;
			//} else {
			//	++rat2;
			}
			*/

			//cout << rat1/(rat1+rat2) << endl;
			
			//assert(topd1.checksolution()==0);
		}
		
		if ( move_p.status == 1){
			vmoves.push_back(move_p);
			int kstep = 1;
			infomove move_d = topd1.adddelivery(topd1.heursol);
			assert(topd1.checksolution()==0);
			int test = 0;
			while( move_d.status == 1){
				++kstep ;
				test = 1;

				//myfloattype b = topd1.optimizequantities2( topd1.heursol);

				TWOOPT twoopt(&topd1);
				//twoopt.LS();
				//cout << i << endl;
				myfloattype b = topd1.optimizequantities2(topd1.heursol);
				myfloattype a = b;
				//cout << -i << endl;
				//a = topd1.optimizequantities( topd1.heursol);			
				if(abs(a-b)>TOL){
				   cout << "i \t" << i << " " << a << "," << b << endl;
				}
				TWOOPT twoopt2(&topd1);
				twoopt2.LS();
				if(topd1.elitists.ub + TOL < topd1.heursol.profit){
					topd1.elitists.member.push_back(topd1.heursol);
					topd1.elitists.ub = topd1.heursol.profit;
					assert(topd1.checksolution()==0);
					cout << i <<  "\t" << topd1.heursol.profit << endl;
					char buffer [50];
				    sprintf_s (buffer,"heu-%f-",topd1.heursol.profit);
					//topd1.Rprintsol(rexe,rpath,filename,buffer,topd1.heursol);
				}
				move_d = topd1.adddelivery(topd1.heursol);
				assert(topd1.checksolution()==0);
			}
			if (test == 1){
				continue;
			}

			topd1.destroy(topd1.heursol,rand()%140);
			assert(topd1.checksolution()==0);
			continue;
		}
		if(vmoves.size() == 0){
			continue;
		}
		//int delstart = 1 + (rand() % (topd1.heursol.tour.size()-2));
		//for(int k = topd1.heursol.tour.size()-2; k >= delstart;--k){
		//	topd1.deletenode(k-1);
		//	topd1.checksolution();
		//}
		topd1.destroy(topd1.heursol,rand()%140);
		assert(topd1.checksolution()==0);
	}
}
	
//static int foobar(void *info, char *buf) { return 1; }

/* disable output from glpk api routines */
//_glp_lib_print_hook(foobar, NULL);


/* enable output from glpk api routines */
//_glp_lib_print_hook(NULL, NULL);

#ifdef GLPKAVAILABLE
int testglpk(){
	cout.precision(6); //Make the output precision to 6 significant digits
	glp_prob *lp = glp_create_prob();
	glp_term_out(0);
	int ia[1+1000], ja[1+1000];
	double ar[1+1000], z, x1, x2;
	lp = glp_create_prob();
	glp_set_prob_name(lp, "Carpinteria");
	glp_set_obj_dir(lp, GLP_MAX);
	glp_add_rows(lp, 3);
	glp_set_row_name(lp, 1, "f");
	glp_set_row_bnds(lp, 1, GLP_UP, 0.0, 100.0);
	glp_set_row_name(lp, 2, "c");
	glp_set_row_bnds(lp, 2, GLP_UP, 0.0, 80.0);
	glp_set_row_name(lp, 3, "r");
	glp_set_row_bnds(lp, 3, GLP_UP, 0.0, 40.0);
	glp_add_cols(lp, 2);
	glp_set_col_name(lp, 1, "x1");
	glp_set_col_bnds(lp, 1, GLP_LO, 0.0, 0.0);
	glp_set_obj_coef(lp, 1, 3.0);
	glp_set_col_name(lp, 2, "x2");
	glp_set_col_bnds(lp, 2, GLP_LO, 0.0, 0.0);
	glp_set_obj_coef(lp, 2, 2.0);

	ia[1] = 1, ja[1] = 1, ar[1] = 2.0; /* a[1,1] = 2 */
	ia[2] = 1, ja[2] = 2, ar[2] = 1.0; /* a[1,2] = 1 */
	ia[3] = 2, ja[3] = 1, ar[3] = 1.0; /* a[2,1] = 1 */
	ia[4] = 2, ja[4] = 2, ar[4] = 1.0; /* a[2,2] = 1 */
	ia[5] = 3, ja[5] = 1, ar[5] = 1.0; /* a[3,1] = 1 */
	ia[6] = 3, ja[6] = 2, ar[6] = 0.0; /* a[3,2] = 0 */

	glp_load_matrix(lp, 6, ia, ja, ar);

	glp_simplex(lp, NULL );
	z = glp_get_obj_val(lp);
	x1 = glp_get_col_prim(lp, 1);
	x2 = glp_get_col_prim(lp, 2);

	printf("\nz = %g; x1 = %g; x2 = %g; \n",
	z, x1, x2);
	glp_delete_prob(lp);
	return 0;


}


int testglpk2(vector<myfloattype> balancetour,vector<myfloattype> profittour){ 
	//cout.precision(6); //Make the output precision to 6 significant digits
	glp_prob *lp = glp_create_prob();
	glp_term_out(1);
	myfloattype Lmax = 300;

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
		sprintf_s (strbuffer, "q%d",i+1);
		glp_set_col_name(lp, N+1+i+1, strbuffer);
		if( balancetour[i]>0 ){
			glp_set_col_bnds(lp, N+1+i+1, GLP_DB, 0.0, balancetour[i]);
		}
		else {
			glp_set_col_bnds(lp, N+1+i+1, GLP_DB, 0.0, -balancetour[i]);
		}
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
		if( balancetour[i-2]>0 ){
			ar[k] = -1;
		}
		else {
			ar[k] = 1;
		};
		++k;
	}

	for (int i = 0; i <  K;  ++i) {
		cout  << i << "\t" << ia[i]<< "\t"  << ja[i] << "\t" << ar[i] << endl;
	}

	glp_load_matrix(lp, K, ia, ja, ar);
	 glp_write_lp(lp, NULL,"test.lp");
	glp_simplex(lp, NULL );

    for (int i = 0; i <  2*N+1;  ++i) {
		cout << glp_get_col_prim(lp, i+1) << endl;
	}

	delete ia; 
	delete ja; 
	delete ar;

	return 0;
}
#endif

int main(int argc, char **argv) {
	// -cf 1 -type CPLEX -CPXtime 15 -CPXpar heu -path "C:/Users/Martin/Desktop/projects/pickupanddeliveryteamorienteering/orienteering instances/dat2/" -file tsiligirides_problem_2_budget_32_300.txt -rexe "C:\Program Files\R\R-2.14.2\bin\R"
	// faster implementation of shuffle
	// testshuffle();
	
    #ifdef GLPKAVAILABLE
	printf("GLPK Version: %s\n", glp_version());
	testglpk();

	vector<double> a;
	a.push_back(rand()%100-50);
	a.push_back(rand()%100-50);
	a.push_back(rand()%100-50);
	a.push_back(rand()%100-50);
	a.push_back(rand()%100-50);
	a.push_back(rand()%100-50);

	vector<double> b;
	b.push_back(-(rand()%100)*a[0]);
	b.push_back(-(rand()%100)*a[1]);
	b.push_back(-(rand()%100)*a[2]);
	b.push_back(-(rand()%100)*a[3]);
	b.push_back(-(rand()%100)*a[4]);
	b.push_back(-(rand()%100)*a[5]);

	testglpk2(a,b);
	#endif
    
	testmultimap();

	//filename2 = path.append(filename);	


	TOPDparameters parameters;

	double cpxtime = 20;
	if( argc > 0){	
		string line = "";
		int i = 1; 
		while( i+1 < argc){
			line = argv[i];
			if(line.compare("-type") == 0){i++; parameters.type=argv[i];}; 
			if(line.compare("-CPXtime") == 0){i++; parameters.ILOMaxCompTime =to_double(argv[i]);}; 
			if(line.compare("-path") == 0){i++; parameters.path=argv[i];}; 
			if(line.compare("-rpath") == 0){i++; parameters.rpath=argv[i];}; 
			if(line.compare("-rexe") == 0){i++;	parameters.rexe = std::string("\"") + argv[i] + std::string("\"");}; 
			if(line.compare("-file") == 0){i++; parameters.filename2 = argv[i]; parameters.filename = parameters.path.append(argv[i]);}; 
			if(line.compare("-CPXpar") == 0){i++; parameters.CPXpar=argv[i];}; 
			if(line.compare("-cf") == 0){i++; parameters.TourCostFactor =to_double(argv[i]);};
			i++;
		}
	}

	srand(1000);
	TOPD topd1(parameters);
	reverseheu(topd1,1000);
	srand(1000);
	TOPD topd2(parameters);
	reverseheu2(topd2,1000);
	srand(2345);
	TOPD topd3(parameters);
	// 2 zeros were deleted
	simpleLNS(topd3,100000);

	int h = 0;
	cin >> h;
	
	#ifdef CPLEXAVAILABLE
	for(int i = 1; i < 1+3; i++){
		TOPD topd(parameters);
		topd.ILOcreateandsolve();
		topd.Rprintsol(parameters.rexe,parameters.rpath,parameters.filename2,"cplex",topd.ILOsol);
	}
	return 0;
	#endif
}