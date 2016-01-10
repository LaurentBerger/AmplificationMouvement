#include "IIRFilter.hpp"
#include <opencv2/core/core.hpp>

using namespace std;

// Some filtering code is based on C/C++ code from Exstrom Laboratories, LLC.
// The code from Extrom Laboratories LLC is licensed under the GNU Public License, version 2:
// http://www.gnu.org/copyleft/gpl.html.  




IIRFilter::IIRFilter(string ftype, int order, double fs, vector<double> frequency)
{
    if (ftype != "butterworth")
    {
        cv::Exception e;
        e.code = -2;
        e.msg = "only butterworth filter is implemented";
    }
    // Frequency Bande passante où le gain du filtre est de 1 : 
    // 0 f1 passebas
    // f1 f2 passebande
    // f1 fs/2 passe haut
    if (frequency.size() != 2)
    {
        cv::Exception e;
        e.code = -1;
        e.msg = "you must give fLow and fHigh";

        throw e;

    }
    n =order;
    if (frequency[0] == 0)
    {
        // Low pass
        fLow = 2*frequency[1]/fs;
        a = CoefDenominatorLP();
        vector<int> bb = CoefNumeratorLP();
        double s=NormalisationLP();
        b.resize(bb.size());
        for (int i = 0; i<b.size();i++)
            b[i] = bb[i]*s;
    }
    else if (frequency[1] == fs / 2)
    {
        // High pass
        fHigh = 2/fs*frequency[0];
        a = CoefDenominatorHP();
        vector<int> bb = CoefNumeratorHP();
        b.resize(bb.size());
        double s=NormalisationHP();
        for (int i = 0; i<b.size();i++)
            b[i] = bb[i]*s;
    }
    else 
    {
        // High pass
        fLow= 2/fs*frequency[0];
        fHigh = 2/fs*frequency[1];
        a = CoefDenominatorBP();
        vector<int> bb = CoefNumeratorBP();
        b.resize(bb.size());
        double s=NormalisationBP();
        for (int i = 0; i<b.size();i++)
            b[i] = bb[i]*s;
    }

}

/**********************************************************************
    dcof_bwbp - calculates the d coefficients for a butterworth bandpass 
    filter. The coefficients are returned as an array of doubles.
*/

vector<double> IIRFilter::CoefDenominatorBP()
{
    int k;            // loop variables
    double theta;     // M_PI * (f2f - f1f) / 2.0
    double cp;        // cosine of phi
    double st;        // sine of theta
    double ct;        // cosine of theta
    double s2t;       // sine of 2*theta
    double c2t;       // cosine 0f 2*theta
    vector<double> rcof(2*n);     // z^-2 coefficients
    vector<double> tcof(2*n);     // z^-1 coefficients
    vector<double> dcof;     // dk coefficients
    double parg;      // pole angle
    double sparg;     // sine of pole angle
    double cparg;     // cosine of pole angle
    double a;         // workspace variables

    cp = cos(M_PI * (fLow + fHigh) / 2.0);
    theta = M_PI * (fHigh - fLow) / 2.0;
    st = sin(theta);
    ct = cos(theta);
    s2t = (2.0 * st * ct);        // sine of 2*theta
    c2t = (2.0 * ct * ct - 1.0);  // cosine of 2*theta


    for (k = 0; k < n; ++k)
    {
        parg = M_PI * (2 * k + 1.) / (2 * n);
        sparg = sin(parg);
        cparg = cos(parg);
        a = (1.0 + s2t * sparg);
        rcof[2 * k] = c2t / a;
        rcof[2 * k + 1] = s2t * cparg / a;
        tcof[2 * k] =  (- 2.0 * cp * (ct + st * sparg) / a);
        tcof[2 * k + 1] = (-2.0 * cp * st * cparg / a);
    }

    dcof = MultPoly3( tcof, rcof);

    dcof[1] = dcof[0];
    dcof[0] = 1.0;
    for (k = 3; k <= 2 * n; ++k)
        dcof[k] = dcof[2 * k - 2];
    return dcof;
}


/**********************************************************************
  dcof_bwlp - calculates the d coefficients for a butterworth lowpass 
  filter. The coefficients are returned as an array of doubles.
*/

vector<double> IIRFilter::CoefDenominatorLP(  )
{
    int k;            // loop variables
    double theta;     // M_PI * fcf / 2.0
    double st;        // sine of theta
    double ct;        // cosine of theta
    double parg;      // pole angle
    double sparg;     // sine of the pole angle
    double cparg;     // cosine of the pole angle
    double a;         // workspace variable
    vector<double> rcof(2 * n);     // binomial coefficients
    vector<double> dcof;     


    theta = M_PI * fLow;
    st = sin(theta);
    ct = cos(theta);

    for( k = 0; k < n; ++k )
    {
	parg = M_PI * (double)(2*k+1)/(double)(2*n);
	sparg = sin(parg);
	cparg = cos(parg);
	a = 1.0 + st*sparg;
	rcof[2*k] = -ct/a;
	rcof[2*k+1] = -st*cparg/a;
    }

    dcof = MultPoly2(  rcof );

    dcof[1] = dcof[0];
    dcof[0] = 1.0;
    for( k = 3; k <= n; ++k )
        dcof[k] = dcof[2*k-2];
    return  dcof ;
}

/**********************************************************************
  dcof_bwhp - calculates the d coefficients for a butterworth highpass 
  filter. The coefficients are returned as an array of doubles.
*/

vector<double> IIRFilter::CoefDenominatorHP(  )
{
    double tmp=fLow;
    fLow=fHigh;
    return  CoefDenominatorLP();
}



vector<double> IIRFilter::MultPoly3( vector<double> &b, vector<double> &c)
        {
            int i, j;
            double a;
           vector<double> aa(4 * n);

            aa[2] = c[0];
            aa[3] = c[1];
            aa[0] = b[0];
            aa[1] = b[1];

            for (i = 1; i < n; ++i)
            {
                aa[2 * (2 * i + 1)] += c[2 * i] * aa[2 * (2 * i - 1)] - c[2 * i + 1] * aa[2 * (2 * i - 1) + 1];
                aa[2 * (2 * i + 1) + 1] += c[2 * i] * aa[2 * (2 * i - 1) + 1] + c[2 * i + 1] * aa[2 * (2 * i - 1)];

                for (j = 2 * i; j > 1; --j)
                {
                    aa[2 * j] += b[2 * i] * aa[2 * (j - 1)] - b[2 * i + 1] * aa[2 * (j - 1) + 1] +
                    c[2 * i] * aa[2 * (j - 2)] - c[2 * i + 1] * aa[2 * (j - 2) + 1];
                    aa[2 * j + 1] += b[2 * i] * aa[2 * (j - 1) + 1] + b[2 * i + 1] * aa[2 * (j - 1)] +
                    c[2 * i] * aa[2 * (j - 2) + 1] + c[2 * i + 1] * aa[2 * (j - 2)];
                }

                aa[2] += b[2 * i] * aa[0] - b[2 * i + 1] * aa[1] + c[2 * i];
                aa[3] += b[2 * i] * aa[1] + b[2 * i + 1] * aa[0] + c[2 * i + 1];
                aa[0] += b[2 * i];
                aa[1] += b[2 * i + 1];
            }
            return (aa);
        }


vector<double> IIRFilter::MultPoly2( vector<double> &p )
{
    int i, j;
    vector<double> aa(2*n);
    double a;


    for( i = 0; i < n; ++i )
    {
	for( j = i; j > 0; --j )
	{
	    aa[2*j] += p[2*i] * aa[2*(j-1)] - p[2*i+1] * aa[2*(j-1)+1];
	    aa[2*j+1] += p[2*i] * aa[2*(j-1)+1] + p[2*i+1] * aa[2*(j-1)];
	}
	aa[0] += p[2*i];
	aa[1] += p[2*i+1];
    }
    return aa ;
}

/**********************************************************************
  ccof_bwlp - calculates the c coefficients for a butterworth lowpass 
  filter. The coefficients are returned as an array of integers.
*/

vector<int> IIRFilter::CoefNumeratorLP()
{
    vector<int> coef(n+1);
    int m;
    int i;


    coef[0] = 1;
    coef[1] = n;
    m = n/2;
    for( i=2; i <= m; ++i)
    {
        coef[i] = (n-i+1)*coef[i-1]/i;
        coef[n-i]= coef[i];
    }
    coef[n-1] = n;
    coef[n] = 1;

    return coef ;
}

/**********************************************************************
  ccof_bwhp - calculates the c coefficients for a butterworth highpass 
  filter. The coefficients are returned as an array of integers.
*/

vector<int> IIRFilter::CoefNumeratorHP(  )
{
    vector<int> coef;
    int i;

    coef = CoefNumeratorLP(  );

    for( i = 0; i <= n; ++i)
        if( i % 2 ) coef[i] = -coef[i];

    return coef ;
}

/**********************************************************************
  ccof_bwbp - calculates the c coefficients for a butterworth bandpass 
  filter. The coefficients are returned as an array of integers.
*/

vector<int> IIRFilter::CoefNumeratorBP( )
{
    vector<int> tcof;
    vector<int> coef(2*n+1);
    int i;


    tcof = CoefNumeratorHP();

    for( i = 0; i < n; ++i)
    {
        coef[2*i] = tcof[i];
        coef[2*i+1] = 0.0;
    }
    coef[2*n] = tcof[n];

    return coef ;
}




/**********************************************************************
  sf_bwlp - calculates the scaling factor for a butterworth lowpass filter.
  The scaling factor is what the c coefficients must be multiplied by so
  that the filter response has a maximum value of 1.
*/

double IIRFilter::NormalisationLP(  )
{
    int m, k;         // loop variables
    double omega;     // M_PI * fcf
    double fomega;    // function of omega
    double parg0;     // zeroth pole angle
    double sf;        // scaling factor

    omega = M_PI * fLow;
    fomega = sin(omega);
    parg0 = M_PI / (double)(2*n);

    m = n / 2;
    sf = 1.0;
    for( k = 0; k < n/2; ++k )
        sf *= 1.0 + fomega * sin((double)(2*k+1)*parg0);

    fomega = sin(omega / 2.0);

    if( n % 2 ) sf *= fomega + cos(omega / 2.0);
    sf = pow( fomega, n ) / sf;

    return(sf);
}

/**********************************************************************
  sf_bwhp - calculates the scaling factor for a butterworth highpass filter.
  The scaling factor is what the c coefficients must be multiplied by so
  that the filter response has a maximum value of 1.
*/

double IIRFilter::NormalisationHP( )
{
    int m, k;         // loop variables
    double omega;     // M_PI * fcf
    double fomega;    // function of omega
    double parg0;     // zeroth pole angle
    double sf;        // scaling factor

    omega = M_PI * fHigh;
    fomega = sin(omega);
    parg0 = M_PI / (double)(2*n);

    m = n / 2;
    sf = 1.0;
    for( k = 0; k < n/2; ++k )
        sf *= 1.0 + fomega * sin((double)(2*k+1)*parg0);

    fomega = cos(omega / 2.0);

    if( n % 2 ) sf *= fomega + sin(omega / 2.0);
    sf = pow( fomega, n ) / sf;

    return(sf);
}

/**********************************************************************
  sf_bwbp - calculates the scaling factor for a butterworth bandpass filter.
  The scaling factor is what the c coefficients must be multiplied by so
  that the filter response has a maximum value of 1.
*/

double IIRFilter::NormalisationBP(  )
{
    int k;            // loop variables
    double ctt;       // cotangent of theta
    double sfr, sfi;  // real and imaginary parts of the scaling factor
    double parg;      // pole angle
    double sparg;     // sine of pole angle
    double cparg;     // cosine of pole angle
    double a, b, c;   // workspace variables

    ctt = 1.0 / tan(M_PI * (fHigh - fLow) / 2.0);
    sfr = 1.0;
    sfi = 0.0;

    for( k = 0; k < n; ++k )
    {
	parg = M_PI * (double)(2*k+1)/(double)(2*n);
	sparg = ctt + sin(parg);
	cparg = cos(parg);
	a = (sfr + sfi)*(sparg - cparg);
	b = sfr * sparg;
	c = -sfi * cparg;
	sfr = b - c;
	sfi = a - b - c;
    }

    return( 1.0 / sfr );
}
