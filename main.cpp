#include "IIRFilter.hpp"
#include<iostream>
using namespace std;

int main(int argc, char **argv)
{
    std::vector<double> pb={10,12};
    IIRFilter f("butterworth",2,30,pb);
    for (int i = 0; i<f.a.size();i++)
        cout << f.a[i] << "\t";
        cout <<  "\n";
    for (int i = 0; i<f.b.size();i++)
        cout << f.b[i] << "\t";
        cout <<  "\n";
}