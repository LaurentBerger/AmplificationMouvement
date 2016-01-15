#include "opencv2/opencv.hpp"
#include "IIRFilter.hpp"
#include<iostream>

using namespace std;
using namespace cv;


class PyramideGaussienne {
vector<Mat> pyr;
public :
    PyramideGaussienne(Mat );
    vector<Mat> &get(){return pyr;};

};

class PyramideLaplacienne {
vector<Mat> pyr;
public :
    PyramideLaplacienne(vector<Mat> &);
    vector<Mat> &get(){return pyr;};

};

class PyramideRiesz {
vector<Mat> xPyr;
vector<Mat> yPyr;
public :
    PyramideRiesz(vector<Mat> &);
    vector<Mat> &get(){return xPyr;};
    vector<Mat> &getx(){return xPyr;};
    vector<Mat> &gety(){return yPyr;};

};


PyramideGaussienne::PyramideGaussienne(Mat m)
{
    Mat x=m;
    Mat y;
    pyr.push_back(x);
    while (x.rows >= 32 && x.cols > 32)
    {
        pyrDown(x,y);
        pyr.push_back(y);
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
        pyr.push_back(tmp2);
    }

}

PyramideRiesz::PyramideRiesz(vector<Mat> &m)
{
    Mat xKernel=(Mat_<float>(3,3) << 0, 0, 0, -0.5, 0, 0.5, 0, 0, 0);
    Mat yKernel=(Mat_<float>(3,3) << 0, -0.5, 0, 0, 0, 0, 0, 0.5, 0);
    for (int i = 0; i<m.size()-1;i++)
    {   
        Mat tmp;
        filter2D(m[i],tmp,CV_32F,xKernel);
        xPyr.push_back(tmp);
        filter2D(m[i],tmp,CV_32F,yKernel);
        yPyr.push_back(tmp);
    }

}


int main(int argc, char **argv)
{
    VideoCapture vid;


    vid.open("C:\\Users\\laurent\\Documents\\Visual Studio 2015\\AmplificationMouvement\\baby_mp4.mp4");
    if (!vid.isOpened())
    {
        cout << "Video not opened!\n";
        exit(0);
    }
    Mat m;

    vid.read(m);
    PyramideGaussienne pgPre(m);
    PyramideLaplacienne plPre(pgPre.get());
    PyramideRiesz prPre(plPre.get());

	while (vid.read(m))
	{
		PyramideGaussienne pg(m);
		PyramideLaplacienne pl(pg.get());
		PyramideRiesz pr(pl.get());

	}
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