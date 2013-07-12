
#include "ExcelExporter.h"
#include "OrienteeringProblemWithPickupsAndDeliveries.h"
#include <iostream>
#include <fstream>
#include <ole2.h> // OLE2 Definitions

using namespace std;


// AutoWrap() - Automation helper function...
HRESULT AutoWrap(int autoType, VARIANT *pvResult, IDispatch *pDisp, LPOLESTR ptName, int cArgs...)
{
    // Begin variable-argument list...
    va_list marker;
    va_start(marker, cArgs);

    if(!pDisp) {
        cout << "NULL IDispatch passed to AutoWrap()";
        _exit(0);
    }

    // Variables used...
    DISPPARAMS dp = { NULL, NULL, 0, 0 };
    DISPID dispidNamed = DISPID_PROPERTYPUT;
    DISPID dispID;
    HRESULT hr;
    char buf[200];
    char szName[200];

    
    // Convert down to ANSI
    WideCharToMultiByte(CP_ACP, 0, ptName, -1, szName, 256, NULL, NULL);
    
    // Get DISPID for name passed...
    hr = pDisp->GetIDsOfNames(IID_NULL, &ptName, 1, LOCALE_USER_DEFAULT, &dispID);
    if(FAILED(hr)) {
        sprintf_s(buf, "IDispatch::GetIDsOfNames(\"%s\") failed w/err 0x%08lx", szName, hr);
        cout << buf;
        _exit(0);
        return hr;
    }
    
    // Allocate memory for arguments...
    VARIANT *pArgs = new VARIANT[cArgs+1];
    // Extract arguments...
    for(int i=0; i<cArgs; i++) {
        pArgs[i] = va_arg(marker, VARIANT);
    }
    
    // Build DISPPARAMS
    dp.cArgs = cArgs;
    dp.rgvarg = pArgs;
    
    // Handle special-case for property-puts!
    if(autoType & DISPATCH_PROPERTYPUT) {
        dp.cNamedArgs = 1;
        dp.rgdispidNamedArgs = &dispidNamed;
    }
    
    // Make the call!
    hr = pDisp->Invoke(dispID, IID_NULL, LOCALE_SYSTEM_DEFAULT, autoType, &dp, pvResult, NULL, NULL);
    if(FAILED(hr)) {
        sprintf_s(buf, "IDispatch::Invoke(\"%s\"=%08lx) failed w/err 0x%08lx", szName, dispID, hr);
        cout << buf;
        _exit(0);
        return hr;
    }
    // End variable-argument section...
    va_end(marker);
    
    delete [] pArgs;
    
    return hr;
}



ExcelExporter::ExcelExporter(vector <std::vector <double>> totalFinalSolutions)
{

	// Initialize COM for this thread...
   CoInitialize(NULL);

   // Get CLSID for our server...
   CLSID clsid;
   HRESULT hr = CLSIDFromProgID(L"Excel.Application", &clsid);

   if(FAILED(hr)) {

      cout << "CLSIDFromProgID() failed";
      return;
   }

   // Start server and get IDispatch...
   IDispatch *pXlApp;
   hr = CoCreateInstance(clsid, NULL, CLSCTX_LOCAL_SERVER, IID_IDispatch, (void **)&pXlApp);
   if(FAILED(hr)) {
      cout << "Excel not registered properly";
      return;
   }

   // Make it visible (i.e. app.visible = 1)
   {

      VARIANT x;
      x.vt = VT_I4;
      x.lVal = 1;
      AutoWrap(DISPATCH_PROPERTYPUT, NULL, pXlApp, L"Visible", 1, x);
   }

   // Get Workbooks collection
   IDispatch *pXlBooks;
   {
      VARIANT result;
      VariantInit(&result);
      AutoWrap(DISPATCH_PROPERTYGET, &result, pXlApp, L"Workbooks", 0);
      pXlBooks = result.pdispVal;
   }

   // Call Workbooks.Add() to get a new workbook...
   IDispatch *pXlBook;
   {
      VARIANT result;
      VariantInit(&result);
      AutoWrap(DISPATCH_PROPERTYGET, &result, pXlBooks, L"Add", 0);
      pXlBook = result.pdispVal;
   }

   // Create a 15x15 safearray of variants...
   // TODO: Estimate how large the safearray has to be (how large is vector totalFinalSolutions)
   VARIANT arr;
   ULONG rowCount;
   ULONG colCount;
   arr.vt = VT_ARRAY | VT_VARIANT;
   {
		

		rowCount = totalFinalSolutions.size();
		colCount= totalFinalSolutions[0].size(); // +1 is the plus column for the time

      SAFEARRAYBOUND sab[2];
      sab[0].lLbound = 1; sab[0].cElements = rowCount; // statt 15 rowcount
      sab[1].lLbound = 1; sab[1].cElements = colCount; // statt 15 colcount
      arr.parray = SafeArrayCreate(VT_VARIANT, 2, sab);
   }

   // Fill safearray with some values...
   // TODO: Copy data from totalFinalSolutions into our safearray
   for(int i=0; i<rowCount; i++) { //rowcount
      for(int j=0; j<colCount; j++) { // colcount
         // Create entry value for (i,j)
         VARIANT tmp;
         tmp.vt = VT_I4;
         tmp.lVal = totalFinalSolutions[i][j]; // TODO: Apply a value here // totalFinalSolutions [i] [j]
         // Add to safearray...
         long indices[] = {i,j};
         SafeArrayPutElement(arr.parray, indices, (void *)&tmp);
      }
   }

   // Get ActiveSheet object
   IDispatch *pXlSheet;
   {
      VARIANT result;
      VariantInit(&result);
      AutoWrap(DISPATCH_PROPERTYGET, &result, pXlApp, L"ActiveSheet", 0);
      pXlSheet = result.pdispVal;
   }

   // Get Range object for the Range A1:O15...
   IDispatch *pXlRange;
   {
      VARIANT parm;
      parm.vt = VT_BSTR;





      parm.bstrVal = ::SysAllocString(L"A1:O15"); /// ????

      VARIANT result;
      VariantInit(&result);
      AutoWrap(DISPATCH_PROPERTYGET, &result, pXlSheet, L"Range", 1, parm);
      VariantClear(&parm);

      pXlRange = result.pdispVal;
   }

   // Set range with our safearray...
   AutoWrap(DISPATCH_PROPERTYPUT, NULL, pXlRange, L"Value", 1, arr);

   // Wait for user...
   cout << "All done.";

   // Set .Saved property of workbook to TRUE so we aren't prompted
   // to save when we tell Excel to quit...
   /*{
      VARIANT x;
      x.vt = VT_I4;
      x.lVal = 1;
      AutoWrap(DISPATCH_PROPERTYPUT, NULL, pXlBook, L"Saved", 1, x);
   }*/

   // TODO: Instead of setting saved property please set name property for a filename and then call save!

   // Tell Excel to quit (i.e. App.Quit)
//   AutoWrap(DISPATCH_METHOD, NULL, pXlApp, L"Quit", 0);

   // Release references...
   pXlRange->Release();
   pXlSheet->Release();
   pXlBook->Release();
   pXlBooks->Release();
   pXlApp->Release();
   VariantClear(&arr);

   // Uninitialize COM for this thread...
   CoUninitialize();




}


