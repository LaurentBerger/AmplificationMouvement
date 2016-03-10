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
	x.convertTo(uc, CV_8U,255/(maxVal[0]-minVal[0]),-0*minVal[0]/(maxVal[0]-minVal[0]));
	imshow(s, uc);
	waitKey(); 

}


class Pyramid {
protected:
	int nbBand;
    vector<vector<Mat > > pyr;
public :
	Pyramid() {};
	Pyramid(Mat m,int level=-1);
	Pyramid(Pyramid &p, bool zero);
    vector <vector<Mat> > &get(){return pyr;};

    void push_back(Mat m){ pyr.push_back(m); };
	int size() { return pyr.size(); };
	int NbBand() { if (pyr.size() == 0) return 0; return pyr[0].size(); };
	void swap(Pyramid x)
	{
		pyr = x.pyr;
	}

    vector<Mat> & operator [](int i) {return pyr[i];}
	Pyramid& operator=(Pyramid x)
	{
		swap(x);
		return *this;
	}


	Pyramid& operator+=(Pyramid& a)
	{
		if (pyr.size() != a.size() )
		{
			cv::Exception e;
			e.code = -2;
			e.msg = "Pyramid size must be equal";
			throw e;
		}
		if (pyr[0].size() != a.get()[0].size())
		{
			cv::Exception e;
			e.code = -3;
			e.msg = "level 0 size  size must be equal";
			throw e;
		}
		Pyramid p(a, false);
		for (int i = 0; i < pyr.size(); i++)
		{
			for (int j=0; j<a.get()[i].size();j++)
				p.get()[i][j] = pyr[i][j] + a.get()[i][j];
		}
		return p;
	}
	friend Pyramid operator+(Pyramid a, const Pyramid& b)
	{

	}
	virtual void collapse(int nbLevel=-1);
	virtual void reduce(int nbLevel = -1);
};

Pyramid::Pyramid(Pyramid &p, bool zero)
{
	pyr.resize(p.size());
	for (int i = 0; i < p.size(); i++)
	{
		for (int j = 0; j < p.get()[i].size(); j++)
		{
			Mat m;
			if (zero)
				m = Mat::zeros(p.get()[i][j].size(), p.get()[i][j].type());
			else
				m = p.get()[i][j].clone();
			pyr[i] = m;
		}
	}
}

class GaussianPyramid:public Pyramid {
public :
	GaussianPyramid(Mat );

};

class LaplacianPyramid :public Pyramid {
public:
	LaplacianPyramid(Mat);
	Mat Collapse();

};



class SteerableLaplacianPyramid:public Pyramid {
static Mat lowPassFilter;
static Mat highPassFilter;

    void InitFilters();

public :
	SteerableLaplacianPyramid(vector<Mat> &);
	SteerableLaplacianPyramid(Mat &); // construct Laplacian pyramid using http://people.csail.mit.edu/nwadhwa/riesz-pyramid/RieszPyrSupplemental.zip
//    Mat Collapse(Pyramid &gauss);
    Mat Collapse();
 
};

class PyramidRiesz:public Pyramid {




    public :
        PyramidRiesz(vector<Mat> &); // construct Riesz pyramid using laplacian pyramid
		Mat Collapse();

};


GaussianPyramid::GaussianPyramid(Mat m):Pyramid()
{
    Mat x=m;
    Mat y;
    pyr.push_back(x);
    while (x.rows >= MIN_ROW_COL_PYRAMID && x.cols > MIN_ROW_COL_PYRAMID)
    {
		vector<Mat> v;
		pyrDown(x,y);
		v.push_back(y);
        pyr.push_back(v);
        x=y;
    }

}

SteerableLaplacianPyramid::SteerableLaplacianPyramid(vector<Mat> &m)
{
    Mat y;
    for (int i = 1; i<m.size();i++)
    {   
        Mat tmp1,tmp2;
        pyrUp(m[i],tmp1);
        subtract(m[i-1],tmp1,tmp2,noArray(),CV_32F);
        pyr.push_back(tmp2);
    }
    Mat x;
    m[m.size() - 1].convertTo(x,CV_32F);
    pyr.push_back(x);

}


SteerableLaplacianPyramid::SteerableLaplacianPyramid(Mat &m)
{
    if (lowPassFilter.empty())
        InitFilters();

    
    Mat tmp=m;
    while (tmp.rows >= MIN_ROW_COL_PYRAMID && tmp.cols > MIN_ROW_COL_PYRAMID)
    {   
	    Mat tmpLow,tmpHigh;
        filter2D(tmp,tmpHigh,CV_32F,highPassFilter);
        pyr.push_back(tmpHigh);
        filter2D(tmp,tmpLow,CV_32F,lowPassFilter);
//		resize(tmpLow,tmp,Size(),0.5,0.5);
		pyrDown(tmpLow,tmp);
    }
    pyr.push_back(tmp);
}

Mat LaplacianPyramid::Collapse()
{
	Mat x, y;

	y = pyr[pyr.size() - 1][0];
	for (int i = pyr.size() - 2; i >= 0; i--)
	{
		pyrUp(y, x);
		add(x, pyr[i], y);
	}
	return y;

}



Mat SteerableLaplacianPyramid::Collapse()
{
    Mat x,y;
    y = pyr[pyr.size()-1][0];
    for (int i = pyr.size()-2; i>=0;i--)
    {   
        Mat tmp1,tmp2;
        resize(y,x,pyr[i][0].size(),-1,-1,CV_INTER_NN);//pyrUp(y,x);
        filter2D(x,tmp1,CV_32F,lowPassFilter);
        filter2D(pyr[i],tmp2,CV_32F,highPassFilter);
        add(tmp1,tmp2,y);
    }
    return y;
}



void SteerableLaplacianPyramid::InitFilters()
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




PyramidRiesz::PyramidRiesz(vector<Mat> &m)
{

	    Mat xKernel=(Mat_<float>(3,3) << 0, 0, 0, 0.5, 0, -0.5, 0, 0, 0);
	    Mat yKernel=(Mat_<float>(3,3) << 0, .5, 0, 0, 0, 0, 0, -0.5, 0);
	//    Mat yKernel=(Mat_<float>(3,3) << -0.12, -0.34, -0.12, 0,0,0, 0.12, 0.34, 0.12);
	//    Mat xKernel=yKernel.t();

	for (int i = 0; i<m.size()-1;i++)
    {   
		vector<Mat> v;
        Mat tmp1,tmp2;
        filter2D(m[i],tmp1,CV_32F,xKernel);
        v.push_back(tmp1);
		filter2D(m[i],tmp2,CV_32F,yKernel);
        v.push_back(tmp2);
		pyr.push_back(v);
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
        Mat qconjReal=c_real.mul(p_real) + cRzX.mul(pRzX) + cRzY.mul(pRzY);
        Mat qconjX= -c_real.mul(pRzX)+p_real.mul(cRzX);
        Mat qconjY= -c_real.mul(pRzY)+p_real.mul(cRzY);

        Mat num=qconjX.mul(qconjX)+qconjY.mul(qconjY);
        Mat qconjAmplitude;
        sqrt(qconjReal.mul(qconjReal)+num,qconjAmplitude);
		Mat ratio;
		divide(qconjReal, qconjAmplitude, ratio);
        Mat diffPhase = ArcCos(ratio);
        Mat cosAngle;
        Mat sinAngle;
        divide(qconjX,num,cosAngle);
        divide(qconjY,num,sinAngle);
        Mat diffPhaseCos=diffPhase.mul(cosAngle);
        Mat diffPhaseSin=diffPhase.mul(sinAngle);
       
        Mat amplitude;
        sqrt(qconjAmplitude,amplitude);
        v[0]=diffPhaseCos;
        v[1]=diffPhaseSin;
        v[2]=amplitude;
		DisplayImage(c_real, "phasecos");
		DisplayImage(cRzX, "phaseSion");
		DisplayImage(cRzY, "p");
	}
return v;
}

Mat IIRtemporalFilter(IIRFilter &f, Mat phase, vector<Mat> ri)
{
    Mat tmp;
    tmp = f.b[0] * phase + (ri[0]);
    for (int i = 0; i<ri.size()-2;i++)
        ri[i] = f.b[i+1] * phase + (ri[i+1]) + f.a[i+1]*tmp;
    ri[ri.size()-2] = f.b[ri.size()-1] * phase  + f.a[ri.size()-1]*tmp;
    return tmp;
}


Mat AmplitudeWeightedblur(Mat img,Mat weighted,Mat kernelx,Mat kernely)
{
Mat num,den;
Mat m;


sepFilter2D(weighted,den,CV_32F,kernelx,kernely);
multiply(img,weighted,m);
sepFilter2D(m,num,CV_32F,kernelx,kernely);
divide(num,den,m);
//divide(m,weighted,m);
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

Mat SteerableLaplacianPyramid::lowPassFilter;
Mat SteerableLaplacianPyramid::highPassFilter;



int main(int argc, char **argv)
{

	if (0==1)
	{
        std::vector<double> pb={36,62};
        IIRFilter f("butterworth",4,300,pb);
        Mat x=(Mat_<float>(1,1) << 1);
        Mat x0=(Mat_<float>(1,1) << 0);
        Mat x1=(Mat_<float>(1,1) << 0);
        Mat x2=(Mat_<float>(1,1) << 0);
        Mat x3=(Mat_<float>(1,1) << 0);
        Mat x4=(Mat_<float>(1,1) << 0);
        Mat x5=(Mat_<float>(1,1) << 0);
        Mat x6=(Mat_<float>(1,1) << 0);
        Mat x7=(Mat_<float>(1,1) << 0);
        Mat x8=(Mat_<float>(1,1) << 0);
        Mat x9=(Mat_<float>(1,1) << 0);
        Mat x10=(Mat_<float>(1,1) << 0);
        vector<Mat> ri = {x0,x1,x2,x3,x4,x5,x6,x7,x8};

        Mat y;
        for (int i = 2; i < 20; i++)
        {
            y=IIRtemporalFilter(f,x,ri);
            cout<<"x ="<<x<<" y="<<y<<endl;
            x=(Mat_<float>(1,1) << float(i));

        }

    }
    vector<string> listeVideo;
    vector<vector<double> > frequencyBand;
    std::vector<double> pb(3);

    listeVideo.push_back("../Video/baby.avi");pb[0]=0.2;pb[1]=5;pb[2]=30; frequencyBand.push_back(pb);
    listeVideo.push_back("../Video/baby_mp4.mp4");pb[0]=0.2;pb[1]=5;pb[2]=30; frequencyBand.push_back(pb);
    listeVideo.push_back("../Video/baby_blanket.avi");pb[0]=0.2;pb[1]=5;pb[2]=30; frequencyBand.push_back(pb);
    listeVideo.push_back("../Video/camera.avi");pb[0]=6;pb[1]=30;pb[2]=300; frequencyBand.push_back(pb);
    listeVideo.push_back("../Video/balance.avi");pb[0]=1;pb[1]=7;pb[2]=300; frequencyBand.push_back(pb);
    listeVideo.push_back("../Video/drum.avi");pb[0]=74;pb[1]=78;pb[2]=1900; frequencyBand.push_back(pb);
    listeVideo.push_back("../Video/guitar.avi");pb[0]=78;pb[1]=86;pb[2]=600; frequencyBand.push_back(pb);
    listeVideo.push_back("../Video/smoke.avi");pb[0]=9;pb[1]=15;pb[2]=200; frequencyBand.push_back(pb);
    for (int i = 0; i<listeVideo.size();i++)
        cout << i << " for " << listeVideo[i] << " with butterworth band pass filter ["<<frequencyBand[i][0]<<","<<frequencyBand[i][1]<<"]\n";

    string nomFichier,codeFichier;
    int fileIndx=-1; 
    VideoCapture vid;
    do
    {
        cout << "Video file or number : ";
        cin>>codeFichier;
        if (codeFichier.length() == 1)
        {
            fileIndx = stoi(codeFichier);
            if (fileIndx>=0 && fileIndx<listeVideo.size())
            {
                nomFichier = listeVideo[fileIndx];
                pb = frequencyBand[fileIndx];
            }
            else
                fileIndx=-1;
        }
        else
        {
            nomFichier=codeFichier;
            pb[0]=5;
            pb[1]=0;
            pb[2]=30;
        }
        if (nomFichier.length())
            vid.open(nomFichier);
        if(!vid.isOpened())
            cout << "File not found "<<nomFichier<<"\n";
    }
    while (!vid.isOpened());
    cout << "butterworth band pass filter\n";
    cout << "Example Low pass at 3Hz give 3 for f1 and 0 for f2\n";
    cout << "Example band pass at [3Hz,5hz] give 3 for f1 and 5 for f2\n";
    cout << "Example high pass at [3Hz,15 give 3 for f1 and 15 for f2 (fs=30hz=30images per second)\n";
    cout << "Example high pass at [3Hz,15 give 3 for f1 and 15 for f2 (fs=30hz=30images per second)\n";
    if (fileIndx!=-1)
        cout << "give -1 for f1 for default values\n";
    double x;
    cout << "f1 = ";
    cin>>x;
    if (x!=-1)
        pb[0]=x;
    cout << "f2 = ";
    cin>>x;
    if (x!=-1)
        pb[1]=x;
    cout << "Sampling frequency = ";
    cin>>x;
    if (x!=-1)
        pb[2]=x;
    int ordre;
    cout << "Order filter (more than 3  unstable) = ";
    cin>>ordre;
    if (x!=-1)
        pb[2]=x;
    double fs = pb[2];
    pb.resize(2);
    double amplificationfactor=50;
    double kernelSize, std;
    cout << "Amplification factor = ";
    cin>>amplificationfactor;
    cout << "Kernel size = ";
    cin>>kernelSize;
    cout << "Std value = ";
    cin>>std;


    IIRFilter f("butterworth",ordre,fs,pb);

    VideoWriter vidWrite;

    Mat m,m1;

    vid.read(m1);
    Mat mc;
    cvtColor(m1,mc,COLOR_BGR2HSV);
    vector<Mat> sx;
    split(mc,sx);
    m = sx[2];
    vidWrite.open("write.avi",CV_FOURCC('P','I','M','1'),30,m.size());
    PyramidLaplacienne plPre(m);

    PyramidRiesz prPre(plPre.get());

	Pyramid phaseCos( prPre.getx(), true);
	Pyramid phaseSin(prPre.getx(),true);
    vector<Pyramid> rCos(std::max(f.a.size(),f.b.size()));
    vector<Pyramid> rSin(std::max(f.a.size(),f.b.size()));
    for (int i = 0; i < std::max(f.a.size(), f.b.size()); i++)
    {
        Pyramid r0( prPre.getx(), true);
        Pyramid r1( prPre.getx(), true);
        rCos[i] = r0;
        rSin[i] = r1;
    }

	Pyramid motionMagnified(prPre.getx());
    Mat kernelx,kernely;
    vector<Mat> riCos(std::max(f.a.size(), f.b.size()));
    vector<Mat> riSin(std::max(f.a.size(), f.b.size()));
	kernelx =  getGaussianKernel(kernelSize, std);
    kernely = kernelx.t();
    int numLevels = plPre.size()-1;
	while (vid.read(m1))
	{
        
        cvtColor(m1,mc,COLOR_BGR2HSV);
        vector<Mat> sx;
        split(mc,sx);
        m = sx[2];
        imshow("video",m);
        PyramidLaplacienne plAct(m);
		PyramidRiesz prAct(plAct.get());
		PyramidLaplacienne prMotionMagnifiedLap(plAct);
		for (int i = 0; i < numLevels; i++)
		{
            vector<Mat> w=DifferencePhaseAmplitude(plAct[i],prAct.getx()[i],prAct.gety()[i],plPre[i],prPre.getx()[i],prPre.gety()[i]);
			phaseCos[i] += w[0];
			phaseSin[i] += w[1];
            for (int j = 0;j < riCos.size(); j++)
            {
                riCos[j] = rCos[j][i];
                riSin[j] = rSin[j][i];
            }
			Mat phaseFilteredCos=IIRtemporalFilter(f,phaseCos[i],riCos);
			Mat phaseFilteredSin=IIRtemporalFilter(f,phaseSin[i],riSin);


            phaseFilteredCos = AmplitudeWeightedblur(phaseFilteredCos,w[2],kernelx,kernely);
            phaseFilteredSin = AmplitudeWeightedblur(phaseFilteredSin,w[2],kernelx,kernely);
            Mat phaseMagnifiedFilteredCos;
            Mat phaseMagnifiedFilteredSin;

			phaseMagnifiedFilteredCos = amplificationfactor*phaseFilteredCos;
			phaseMagnifiedFilteredSin = amplificationfactor*phaseFilteredSin;
			prMotionMagnifiedLap[i]=PhaseShiftCoefficientRealPart(plAct[i], prAct.getx()[i], prAct.gety()[i], phaseMagnifiedFilteredCos, phaseMagnifiedFilteredSin);

		}
        prMotionMagnifiedLap[numLevels]=plAct[numLevels];
        Mat x = prMotionMagnifiedLap.Collpase();
	    double minVal, maxVal;
		minMaxLoc(x, &minVal, &maxVal);

	    Mat uc;
	    x.convertTo(uc, CV_8U,255/(maxVal-minVal),-minVal/(maxVal-minVal));
	    x.convertTo(uc,CV_8U);
        sx[2]=uc;
        merge(sx,mc);
        Mat mr;
        cvtColor(mc,mr,COLOR_HSV2BGR);
        imshow("Laplacian Motion",uc);
        vidWrite << mr;
        waitKey(1);
        cout << vid.get(CAP_PROP_POS_MSEC)<<endl;
        plPre=plAct;
        prPre=prAct;


	}
}