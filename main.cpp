#include "opencv2/opencv.hpp"
#include "IIRFilter.hpp"
#include<iostream>

using namespace std;
using namespace cv;


class PyramideGaussienne {
vector<Mat> pyramide;
public :
    PyramideGaussienne(Mat );
    vector<Mat> &get(){return pyramide;};

};

class PyramideLaplacienne {
vector<Mat> pyramide;
public :
    PyramideLaplacienne(vector<Mat> &);
    vector<Mat> &get(){return pyramide;};

};


PyramideGaussienne::PyramideGaussienne(Mat m)
{
    Mat x=m;
    Mat y;
    pyramide.push_back(x);
    while (x.rows >= 4 && x.cols > 4)
    {
        pyrDown(x,y);
        pyramide.push_back(y);
        x=y;
    }

}

PyramideLaplacienne::PyramideLaplacienne(vector<Mat> &m)
{
    Mat y;
    for (int i = 1; i<m.size();i++)
    {   
        Mat tmp1,tmp2;
        pyrUp(m[i],tmp1);
        subtract(m[i-1],tmp1,tmp2,noArray(),CV_32F);
        pyramide.push_back(tmp2);
    }

}


int main(int argc, char **argv)
{
    Mat m = imread("f:/lib/opencv/samples/data/lena.jpg");

    PyramideGaussienne pg(m);
    PyramideLaplacienne pl(pg.get());


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