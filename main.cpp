#include "opencv2/opencv.hpp"
#include "IIRFilter.hpp"
#include<iostream>

using namespace std;
using namespace cv;

class Pyramide {
protected:
    vector<Mat> pyr;
public :
    Pyramide(){};
	Pyramide(Pyramide &p, bool zero);
    vector<Mat> &get(){return pyr;};

    void push_back(Mat m){ pyr.push_back(m); };
	int size() { return pyr.size(); };
	void swap(Pyramide x)
	{
		pyr = x.pyr;
	}

    Mat & operator [](int i) {return pyr[i];}
	Pyramide& operator=(Pyramide x)
	{
		swap(x);
		return *this;
	}


	Pyramide& operator+=(Pyramide& a)
	{
		if (pyr.size() != a.size())
		{
			cv::Exception e;
			e.code = -2;
			e.msg = "Pyramide size must be equal";
			throw e;
		}
		if (pyr[0].size() != a.get()[0].size())
		{
			cv::Exception e;
			e.code = -3;
			e.msg = "level 0 size  size must be equal";
			throw e;
		}
		Pyramide p(a, false);
		for (int i = 0; i < pyr.size(); i++)
		{
			p.get()[i] = pyr[i] + a.get()[i];
		}
		return p;
	}
	friend Pyramide operator+(Pyramide a, const Pyramide& b)
	{

	}
};

Pyramide::Pyramide(Pyramide &p, bool zero)
{
	pyr.resize(p.size());
	for (int i = 0; i < p.size(); i++)
	{
		Mat m;
		if (zero)
			m = Mat::zeros(p.get()[i].size(), p.get()[i].type());
		else
			m = p.get()[i].clone();
		pyr[i] = m;
	}
}

class PyramideGaussienne:public Pyramide {
public :
    PyramideGaussienne(Mat );

};

class PyramideLaplacienne:public Pyramide {
public :
    PyramideLaplacienne(vector<Mat> &);
    vector<Mat> &get(){return pyr;};

};

class PyramideRiesz {
Pyramide xPyr;
Pyramide yPyr;
public :
    PyramideRiesz(vector<Mat> &);
    Pyramide &get(){return xPyr;};
    Pyramide &getx(){return xPyr;};
    Pyramide &gety(){return yPyr;};

};


PyramideGaussienne::PyramideGaussienne(Mat m):Pyramide()
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


Mat ArcCos(Mat m)
{
    if (m.depth() != CV_32F )
    {
        cv::Exception e;
        e.code = -1;
        e.msg = "Mat must be real with one channel for ArcCos function";
        throw e;
    }
    vector<double> minMat(m.channels()),maxMat(m.channels());
    minMaxLoc(m,minMat.data(),maxMat.data());
    for (int i = 0; i<m.channels();i++)
        if (abs(minMat[i])>1 || abs(maxMat[i])> 1)
        {
            cv::Exception e;
            e.code = -1;
            e.msg = "mat values must be in range -1 to 1 for ArcCos function";
            throw e;
        }


    Mat t(m.size(),CV_32FC(m.channels()));

    for (int i = 0; i < m.rows; i++)
    {
        float *ptr1 = (float*)m.ptr(i);
        float *ptr2 = (float*)t.ptr(i);
        for (int j=0;j<m.cols*m.channels();j++,ptr1++,ptr2++)
            *ptr2 = acos(*ptr1);
    }

    return t;
}

vector<Pyramide> DifferencePhaseAmplitude(PyramideLaplacienne &plAct, PyramideRiesz &prAct, PyramideLaplacienne &plPre, PyramideRiesz &prPre )
{
    int level=prAct.get().size();
    vector<Pyramide> v(3);

    for (int i = 0; i < level; i++)
    {
        Mat qReal= plAct.get()[i].mul(plPre.get()[i]) + prAct.getx()[i].mul(prPre.getx()[i]) + prAct.gety()[i].mul(prPre.gety()[i]);
        Mat qX= -plAct.get()[i].mul(prPre.getx()[i])+plPre.get()[i].mul(prAct.getx()[i]);
        Mat qY= -plAct.get()[i].mul(prPre.gety()[i])+plPre.get()[i].mul(prAct.gety()[i]);

        Mat num=qX.mul(qX)+qY.mul(qY);
        Mat qAmplitude;
        sqrt(qReal.mul(qReal)+num,qAmplitude);
        Mat diffPhase = ArcCos(qReal/qAmplitude);
        Mat cosAngle=qX/num;
        Mat sinAngle=qY/num;
        Mat diffPhaseCos=diffPhase.mul(cosAngle);
        Mat diffPhaseSin=diffPhase.mul(sinAngle);
       
        Mat amplitude;
        sqrt(qAmplitude,amplitude);
        v[0].push_back(diffPhaseCos);
        v[1].push_back(diffPhaseSin);
        v[2].push_back(amplitude);
    }
return v;
}

Mat IIRtemporalFilter(IIRFilter &f, Mat phase, Mat *r0, Mat *r1)
{
    Mat tmp;
    tmp = f.b[0] * phase + (*r0);
    *r0 = f.b[1] * phase + (*r1) - f.a[1]*tmp;
    *r1 = f.b[2] * phase  - f.a[2]*tmp;
    return tmp;
}


Mat AmplitudeWeightedblur(Mat img,Mat weighted,Mat kernel)
{
Mat num,den;
Mat m;


sepFilter2D(weighted,den,CV_32F,kernel,kernel);
multiply(img,weighted,m);
sepFilter2D(m,num,CV_32F,kernel,kernel);
divide(num,den,m);
return m;
}


Mat PhaseShiftCoefficientRealPart(Mat laplevel, Mat rieszXLevel, Mat rieszYLevel, Mat phasemagnifiedCos, Mat phaseMaginifiedSin)
{
    Mat r;
    return r;
}

int main(int argc, char **argv)
{
     std::vector<double> pb={5,10};
    IIRFilter f("butterworth",1,30,pb);
    /* H(z) =\frac{Y(z)}{X(z)}= \frac{\sum_{k=0}^{N}b_k z^{-k}}{\sum_{k=0}^{M}a_k z^{-k}} */
    /* y(n)=b_0x(n)+...b_N x(n-N)-a_1 y(n-1)-...-a_M y(n-M) */
    VideoCapture vid;
    double amplificationfactor=0.1;

    vid.open("C:\\Users\\Laurent.PC-LAURENT-VISI\\Documents\\Visual Studio 2013\\AmplificationMouvement\\baby_mp4.mp4");
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
	Pyramide phaseCos( prPre.getx(), true);
	Pyramide phaseSin(prPre.getx(),true);
	Pyramide r0Cos( prPre.getx(), true);
	Pyramide r1Cos(prPre.getx(),true);
	Pyramide r0Sin( prPre.getx(), true);
	Pyramide r1Sin(prPre.getx(),true);
	Pyramide motionMagnified(prPre.getx());
    Mat kernel;

    kernel = getGaussianKernel(2,-1);

	while (vid.read(m))
	{
		PyramideGaussienne pgAct(m);
		PyramideLaplacienne plAct(pgAct.get());
		PyramideRiesz prAct(plAct.get());
        Mat amplitude;
		// Valeur de retour 2 pyramides : DiffPaseCos DiffPhaseSin et amplitude
        vector<Pyramide> w=DifferencePhaseAmplitude(plAct,prAct,plPre,prPre);
		for (int i = 0; i < phaseCos.size(); i++)
		{
			phaseCos[i] += w[0][i];
			phaseSin[i] += w[1][i];
            Mat phaseFilterdCos=IIRtemporalFilter(f,phaseCos[i],&r0Cos[i],&r1Cos[i]);
            Mat phaseFilterdSin=IIRtemporalFilter(f,phaseSin[i],&r0Sin[i],&r1Sin[i]);

            phaseFilterdCos = AmplitudeWeightedblur(phaseFilterdCos,w[0][i],kernel);
            phaseFilterdSin = AmplitudeWeightedblur(phaseFilterdSin,w[0][i],kernel);
            Mat phaseMagnifiedFilteredCos;
            Mat phaseMagnifiedFilteredSin;

            phaseMagnifiedFilteredCos = amplificationfactor*phaseFilterdCos;
            phaseMagnifiedFilteredSin = amplificationfactor*phaseFilterdSin;


		}

	}
    /* H(z) =\frac{Y(z)}{X(z)}= \frac{\sum_{k=0}^{N}b_k z^{-k}}{\sum_{k=0}^{M}a_k z^{-k}} */
    /* y(n)=b_0x(n)+...b_N x(n-N)-a_1 y(n-1)-...-a_M y(n-M) */
    cout <<  "Numerator = ";
    for (int i = 0; i<f.b.size();i++)
        cout << f.b[i] << "\t";
    cout <<  "\n                  -----------------------\nDenominator = ";
    for (int i = 0; i<f.a.size();i++)
        cout << f.a[i] << "\t";
    cout <<  "\n  ";
}