#include "IIRFilter.hpp"
#include<iostream>
using namespace std;

int main(int argc, char **argv)
{
    std::vector<double> pb={5,10};
    IIRFilter f("butterworth",4,30,pb);
    cout <<  "Numerator = ";
    for (int i = 0; i<f.b.size();i++)
        cout << f.b[i] << "\t";
    cout <<  "\n                  -----------------------\nDenominator = ";
    for (int i = 0; i<f.a.size();i++)
        cout << f.a[i] << "\t";
    cout <<  "\n  ";
}