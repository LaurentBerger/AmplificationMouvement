#ifndef __IIRFILTER__
#define IIRFILTER__

#include <vector>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>

class IIRFilter 
{
public :
    std::vector <double> a; // denominator  
    std::vector <double> b; // Numerator 
    int n;                  // Filter order y= bx-ay Y=B/A

    double fLow,fHigh;

    IIRFilter(std::string filterType,int order,double fs,std::vector<double> frequency);
    std::vector<int> CoefNumeratorBP();
    std::vector<int> CoefNumeratorLP();
    std::vector<int> CoefNumeratorHP();
    std::vector<double> CoefDenominatorBP();
    std::vector<double> CoefDenominatorLP();
    std::vector<double> CoefDenominatorHP();
private :
    std::vector<double> MultPoly3( std::vector<double> &b, std::vector<double> &c);
    std::vector<double> MultPoly2( std::vector<double> &b);

    double NormalisationLP();

    double NormalisationHP();

    double NormalisationBP();

};


#endif
