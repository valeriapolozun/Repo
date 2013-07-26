#ifndef __EXCELCALCULATOR__
#define __EXCELCALCULATOR__


#include "Coordinates.h"
#include <vector>
#include <string>
using namespace std;

class ExcelExporter
{
public: 
	ExcelExporter(vector <std::vector <double>> totalFinalSolutions);
	void runExcelExporting(vector <std::vector <double>> totalFinalSolutions);

private:
//	HRESULT AutoWrap(int autoType, VARIANT *pvResult, IDispatch *pDisp, LPOLESTR ptName, int cArgs...);

};


#endif