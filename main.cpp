#include "opencv2/opencv.hpp"
#include "IIRFilter.hpp"
#include<iostream>

using namespace std;
using namespace cv;

#define MIN_ROW_COL_PYRAMID 32

void DisplayImage(Mat x,string s)
{
	return;
	vector<Mat> sx;
	split(x, sx);
	vector<double> minVal(3), maxVal(3);
	for (int i = 0; i < sx.size(); i++)
	{
		minMaxLoc(sx[i], &minVal[i], &maxVal[i]);

	}
	maxVal[0] = *max_element(maxVal.begin(), maxVal.end());
	minVal[0] = *min_element(minVal.begin(), minVal.end());
	Mat uc;
	x.convertTo(uc, CV_8U,255/(maxVal[0]-minVal[0]),-255*minVal[0]/(maxVal[0]-minVal[0]));
	imshow(s, uc);
	waitKey(); 

}

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
static Mat lowPassFilter;
static Mat highPassFilter;

    void InitFilters();

public :
    PyramideLaplacienne(vector<Mat> &);
    PyramideLaplacienne(Mat &); // construct Laplacian pyramid using http://people.csail.mit.edu/nwadhwa/riesz-pyramid/RieszPyrSupplemental.zip
    Mat Collpase(Pyramide &gauss);
    Mat Collpase();
    vector<Mat> &get(){return pyr;};

};

class PyramideRiesz {

    Pyramide xPyr;
    Pyramide yPyr;



    public :
        PyramideRiesz(vector<Mat> &); // construct Riesz pyramid using laplacian pyramid
        Pyramide &get(){return xPyr;};
        Pyramide &getx(){return xPyr;};
        Pyramide &gety(){return yPyr;};

};


PyramideGaussienne::PyramideGaussienne(Mat m):Pyramide()
{
    Mat x=m;
    Mat y;
    pyr.push_back(x);
    while (x.rows >= MIN_ROW_COL_PYRAMID && x.cols > MIN_ROW_COL_PYRAMID)
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
    pyr.push_back(m[m.size()-1]);

}


PyramideLaplacienne::PyramideLaplacienne(Mat &m)
{
    if (lowPassFilter.empty())
        InitFilters();

    
    Mat tmpLow,tmpHigh;
    Mat tmp=m;
    while (tmp.rows >= MIN_ROW_COL_PYRAMID && tmp.cols > MIN_ROW_COL_PYRAMID)

    {   
        filter2D(tmp,tmpHigh,CV_32F,highPassFilter);
        pyr.push_back(tmpHigh);
        filter2D(tmp,tmpLow,CV_32F,lowPassFilter);
/*        Mat x=tmpLow+tmpHigh;
		vector<Mat> sx;
		split(x, sx);
		vector<double> minVal(3), maxVal(3);
		for (int i = 0; i < sx.size(); i++)
		{
			minMaxLoc(sx[i], &minVal[i], &maxVal[i]);

		}
		maxVal[0] = *max_element(maxVal.begin(), maxVal.end());
		minVal[0] = *min_element(minVal.begin(), minVal.end());
		Mat uc;
		x.convertTo(uc,CV_8U);
		imshow("lap",uc);
		waitKey();*/
		resize(tmpLow,tmp,Size(tmpLow.cols/2,tmpLow.rows/2));
//		pyrDown(tmpLow,tmp);
    }
    pyr.push_back(tmp);
}

Mat PyramideLaplacienne::Collpase(Pyramide &gauss)
{
    Mat x,y;

    gauss[ gauss.size()-1].convertTo(y,CV_32F);
    for (int i = pyr.size()-1; i>=0;i--)
    {   
        pyrUp(y,y);
        add(y,pyr[i],y);
    }
    return y;
}

Mat PyramideLaplacienne::Collpase()
{
    Mat x,y;

    y = pyr[pyr.size()-1];
    for (int i = pyr.size()-2; i>=0;i--)
    {   
        pyrUp(y,y);
        add(y,pyr[i],y);
    }
    return y;
}



void PyramideLaplacienne::InitFilters()
{
    // Supplemental for Riesz Pyramid for Fast Phase-Based Video Magnification 
    lowPassFilter = (Mat_<float>(9,9)<< -0.0001,   -0.0007,  -0.0023,  -0.0046,  -0.0057,  -0.0046,  -0.0023,  -0.0007,  -0.0001,
	                                    -0.0007,   -0.0030,  -0.0047,  -0.0025,  -0.0003,  -0.0025,  -0.0047,  -0.0030,  -0.0007,
	                                    -0.0023,   -0.0047,   0.0054,   0.0272,   0.0387,   0.0272,   0.0054,  -0.0047,  -0.0023,
	                                    -0.0046,   -0.0025,   0.0272,   0.0706,   0.0910,   0.0706,   0.0272,  -0.0025,  -0.0046,
	                                    -0.0057,   -0.0003,   0.0387,   0.0910,   0.1138,   0.0910,   0.0387,  -0.0003,  -0.0057,
	                                    -0.0046,   -0.0025,   0.0272,   0.0706,   0.0910,   0.0706,   0.0272,  -0.0025,  -0.0046,
	                                    -0.0023,   -0.0047,   0.0054,   0.0272,   0.0387,   0.0272,   0.0054,  -0.0047,  -0.0023,
	                                    -0.0007,   -0.0030,  -0.0047,  -0.0025,  -0.0003,  -0.0025,  -0.0047,  -0.0030,  -0.0007,
	                                    -0.0001,   -0.0007,  -0.0023,  -0.0046,  -0.0057,  -0.0046,  -0.0023,  -0.0007,  -0.0001);
    highPassFilter=(Mat_<float>(9,9)<<   0.0000,   0.0003,   0.0011,   0.0022,   0.0027,   0.0022,   0.0011,   0.0003,   0.0000,
	                                     0.0003,   0.0020,   0.0059,   0.0103,   0.0123,   0.0103,   0.0059,   0.0020,   0.0003,
	                                     0.0011,   0.0059,   0.0151,   0.0249,   0.0292,   0.0249,   0.0151,   0.0059,   0.0011,
	                                     0.0022,   0.0103,   0.0249,   0.0402,   0.0469,   0.0402,   0.0249,   0.0103,   0.0022,
	                                     0.0027,   0.0123,   0.0292,   0.0469,  -0.9455,   0.0469,   0.0292,   0.0123,   0.0027,
	                                     0.0022,   0.0103,   0.0249,   0.0402,   0.0469,   0.0402,   0.0249,   0.0103,   0.0022,
	                                     0.0011,   0.0059,   0.0151,   0.0249,   0.0292,   0.0249,   0.0151,   0.0059,   0.0011,
	                                     0.0003,   0.0020,   0.0059,   0.0103,   0.0123,   0.0103,   0.0059,   0.0020,   0.0003,
	                                     0.0000,   0.0003,   0.0011,   0.0022,   0.0027,   0.0022,   0.0011,   0.0003,   0.0000);
}




PyramideRiesz::PyramideRiesz(vector<Mat> &m)
{

	//    Mat xKernel=(Mat_<float>(3,3) << 0, 0, 0, 0.5, 0, -0.5, 0, 0, 0);
	//    Mat yKernel=(Mat_<float>(3,3) << 0, .5, 0, 0, 0, 0, 0, -0.5, 0);
	    Mat yKernel=(Mat_<float>(3,3) << -0.12, -0.34, -0.12, 0,0,0, 0.12, 0.34, 0.12);
	    Mat xKernel=yKernel.t();

		cout << xKernel << yKernel << endl;
	for (int i = 0; i<m.size()-1;i++)
    {   
        Mat tmp1,tmp2;
        filter2D(m[i],tmp1,CV_32F,xKernel);
        xPyr.push_back(tmp1);
		filter2D(m[i],tmp2,CV_32F,yKernel);
        yPyr.push_back(tmp2);
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

Mat Cos(Mat m)
{
    if (m.depth() != CV_32F )
    {
        cv::Exception e;
        e.code = -1;
        e.msg = "Mat must be real with one channel for Cos function";
        throw e;
    }


    Mat t(m.size(),CV_32FC(m.channels()));

    for (int i = 0; i < m.rows; i++)
    {
        float *ptr1 = (float*)m.ptr(i);
        float *ptr2 = (float*)t.ptr(i);
        for (int j=0;j<m.cols*m.channels();j++,ptr1++,ptr2++)
            *ptr2 = cos(*ptr1);
    }

    return t;
}

Mat Sin(Mat m)
{
    if (m.depth() != CV_32F )
    {
        cv::Exception e;
        e.code = -1;
        e.msg = "Mat must be real with one channel for Cos function";
        throw e;
    }


    Mat t(m.size(),CV_32FC(m.channels()));

    for (int i = 0; i < m.rows; i++)
    {
        float *ptr1 = (float*)m.ptr(i);
        float *ptr2 = (float*)t.ptr(i);
        for (int j=0;j<m.cols*m.channels();j++,ptr1++,ptr2++)
            *ptr2 = sin(*ptr1);
    }

    return t;
}



vector<Mat> DifferencePhaseAmplitude(Mat &c_real, Mat &cRzX, Mat &cRzY, Mat &p_real, Mat &pRzX, Mat &pRzY )
{
    vector<Mat> v(3);

    {
        Mat qReal=c_real.mul(p_real) + cRzX.mul(pRzX) + cRzY.mul(pRzY);
        Mat qX= -c_real.mul(pRzX)+p_real.mul(cRzX);
        Mat qY= -c_real.mul(pRzY)+p_real.mul(cRzY);

        Mat num=qX.mul(qX)+qY.mul(qY);
        Mat qAmplitude;
        sqrt(qReal.mul(qReal)+num,qAmplitude);
		Mat ratio;
		divide(qReal, qAmplitude, ratio);
        Mat diffPhase = ArcCos(ratio);
        Mat cosAngle;
        Mat sinAngle;
        divide(qX,num,cosAngle);
        divide(qY,num,sinAngle);
        Mat diffPhaseCos=diffPhase.mul(cosAngle);
        Mat diffPhaseSin=diffPhase.mul(sinAngle);
       
        Mat amplitude;
        sqrt(qAmplitude,amplitude);
        v[0]=diffPhaseCos;
        v[1]=diffPhaseSin;
        v[2]=amplitude;
		DisplayImage(c_real, "phasecos");
		DisplayImage(cRzX, "phaseSion");
		DisplayImage(cRzY, "p");
	}
return v;
}

Mat IIRtemporalFilter(IIRFilter &f, Mat phase, Mat *r0, Mat *r1)
{
    Mat tmp;
    tmp = f.b[0] * phase + (*r0);
    *r0 = f.b[1] * phase + (*r1) - f.a[1]*tmp;
    if (f.b.size()==3)
        *r1 = f.b[2] * phase  - f.a[2]*tmp;
    else
        *r1 = 0;
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


Mat PhaseShiftCoefficientRealPart(Mat laplevel, Mat rieszXLevel, Mat rieszYLevel, Mat phaseMagnifiedCos, Mat phaseMagnifiedSin)
{
    Mat r;
    Mat pm,expPhaseReal,expPhaseImag;
    Mat expphix,expphiy,tmp;

    sqrt(phaseMagnifiedCos.mul(phaseMagnifiedCos)+phaseMagnifiedSin.mul(phaseMagnifiedSin),pm);
    expPhaseReal = Cos(pm);
    expPhaseImag = Sin(pm);
    divide(expPhaseImag,pm,tmp);
    expphix = phaseMagnifiedCos.mul(tmp);
    expphiy = phaseMagnifiedSin.mul(tmp);

    r = expPhaseReal.mul(laplevel)-expphix.mul(rieszXLevel)-expphiy.mul(rieszYLevel);

    return r;
}

Mat PyramideLaplacienne::lowPassFilter;
Mat PyramideLaplacienne::highPassFilter;



int main(int argc, char **argv)
{
	if (0==1)
	{
		Mat m = imread("c:/lib/opencv/samples/data/lena.jpg",CV_LOAD_IMAGE_GRAYSCALE);// (512, 512, CV_32FC1);
/*    float pi = acos(-1);
    for (int i = 0; i < m.rows; i++)
    {
        float *ptf = (float*)m.ptr(i);
        for (int j = 0; j < m.cols; j++,ptf++)
        {
            double r = (i - 256)*(i - 256) + (j - 256)*(j-256);
            r = sqrt(r);
            if (r>=100 && r<120)
                *ptf=255;
            else
                *ptf=128;
        }
    }*/
    Mat uc;
    m.convertTo(uc,CV_8U,1);
    imshow("wave",uc);
    PyramideGaussienne pgPre(m);
    PyramideLaplacienne plPre(m);
    PyramideRiesz prPre(plPre.get());
    imshow("riesz",prPre.getx()[0]);
    waitKey();
    }


    // std::vector<double> pb={0,0.6};
     std::vector<double> pb={0,100};
    IIRFilter f("butterworth",1,300,pb);
    /* H(z) =\frac{Y(z)}{X(z)}= \frac{\sum_{k=0}^{N}b_k z^{-k}}{\sum_{k=0}^{M}a_k z^{-k}} */
    /* y(n)=b_0x(n)+...b_N x(n-N)-a_1 y(n-1)-...-a_M y(n-M) */
    VideoCapture vid;
    VideoWriter vidWrite;
    double amplificationfactor=1;

//	vid.open("C:\\Users\\Laurent.PC-LAURENT-VISI\\Documents\\Visual Studio 2013\\AmplificationMouvement\\camera.avi");
	vid.open("C:\\Users\\Laurent\\Documents\\Visual Studio 2015\\AmplificationMouvement\\baby_mp4.mp4");
	if (!vid.isOpened())
    {
        cout << "Video not opened!\n";
        exit(0);
    }
    Mat m;

    vid.read(m);
    vidWrite.open("write.avi",CV_FOURCC('M','J','P','G'),30,m.size());
    //PyramideGaussienne pgPre(m);
    PyramideLaplacienne plPre(m);
	if (0==1)
    {
        Mat x = plPre.Collpase();
        vector<Mat> sx;
        split(x,sx);
        imshow("video",m);
        vector<double> minVal(3), maxVal(3);
        minMaxLoc(sx[0],&minVal[0],&maxVal[0]);
        minMaxLoc(sx[1],&minVal[1],&maxVal[1]);
        minMaxLoc(sx[2],&minVal[2],&maxVal[2]);
        maxVal[0] = *max_element(maxVal.begin(),maxVal.end());
        minVal[0] = *min_element(minVal.begin(),minVal.end());
        Mat uc;
        x.convertTo(uc,CV_8U,255/(maxVal[0]-minVal[0]),-255*minVal[0]/(maxVal[0]-minVal[0]));
        imshow("Laplacian Motion",uc);
        waitKey();
    }

    PyramideRiesz prPre(plPre.get());
	Pyramide phaseCos( prPre.getx(), true);
	Pyramide phaseSin(prPre.getx(),true);
	Pyramide r0Cos( prPre.getx(), true);
	Pyramide r1Cos(prPre.getx(),true);
	Pyramide r0Sin( prPre.getx(), true);
	Pyramide r1Sin(prPre.getx(),true);
	Pyramide motionMagnified(prPre.getx());
    Mat kernel;

	kernel = Mat::ones(1, 1, CV_32F);// getGaussianKernel(3, -1);
    int numLevels = plPre.size()-1;
	while (vid.read(m))
	{
		PyramideLaplacienne plAct(m);
		PyramideRiesz prAct(plAct.get());
		PyramideLaplacienne prMotionMagnifiedLap(plAct);
		for (int i = 0; i < numLevels; i++)
		{
            vector<Mat> w=DifferencePhaseAmplitude(plAct[i],prAct.getx()[i],prAct.gety()[i],plPre[i],prPre.getx()[i],prPre.gety()[i]);
			phaseCos[i] += w[0];
			phaseSin[i] += w[1];
			Mat phaseFilterdCos=IIRtemporalFilter(f,phaseCos[i],&r0Cos[i],&r1Cos[i]);
			Mat phaseFilterdSin=IIRtemporalFilter(f,phaseSin[i],&r0Sin[i],&r1Sin[i]);
			//Mat phaseFilterdCos=phaseCos[i];
			//Mat phaseFilterdSin=phaseSin[i];

            phaseFilterdCos = AmplitudeWeightedblur(phaseFilterdCos,w[2],kernel);
            phaseFilterdSin = AmplitudeWeightedblur(phaseFilterdSin,w[2],kernel);
            Mat phaseMagnifiedFilteredCos;
            Mat phaseMagnifiedFilteredSin;

			phaseMagnifiedFilteredCos = amplificationfactor*phaseFilterdCos;
			phaseMagnifiedFilteredSin = amplificationfactor*phaseFilterdSin;
			prMotionMagnifiedLap[i]=PhaseShiftCoefficientRealPart(plAct[i], prAct.getx()[i], prAct.gety()[i], phaseMagnifiedFilteredCos, phaseMagnifiedFilteredSin);

		}
        prMotionMagnifiedLap[numLevels]=plAct[numLevels];
        Mat x = prMotionMagnifiedLap.Collpase();
        vector<Mat> sx;
        split(x,sx);
        imshow("video",m);
        vector<double> minVal(3), maxVal(3);
        minMaxLoc(sx[0],&minVal[0],&maxVal[0]);
        minMaxLoc(sx[1],&minVal[1],&maxVal[1]);
        minMaxLoc(sx[2],&minVal[2],&maxVal[2]);
        maxVal[0] = *max_element(maxVal.begin(),maxVal.end());
        minVal[0] = *min_element(minVal.begin(),minVal.end());
        Mat uc;
        x.convertTo(uc,CV_8U,255/(maxVal[0]-minVal[0]),-255*minVal[0]/(maxVal[0]-minVal[0]));
        imshow("Laplacian Motion",uc);
        vidWrite << uc;
        waitKey(5);
        cout << vid.get(CAP_PROP_POS_MSEC)<<endl;
        prMotionMagnifiedLap[phaseCos.size()] = plAct[phaseCos.size()];
        plPre=plAct;
        prPre=prAct;


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