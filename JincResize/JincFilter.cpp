

#include <math.h>
#include "JincFilter.h"

/************************************************************************/
/* Jinc function implementations                                        */
/************************************************************************/

const double JINC_PI = 3.141592653589793238462643;
const double JINC_EPSILON = 1e-15;

/*
  See Pratt "Digital Image Processing" p.97 for Jinc/Bessel functions.
  http://mathworld.wolfram.com/JincFunction.html and page 11 of
  http://www.ph.ed.ac.uk/%7ewjh/teaching/mo/slides/lens/lens.pdf

  The original "zoom" program by Paul Heckbert called this "Bessel".  But
  really it is more accurately named "Jinc".

  Note: Most of code in this file are copied from ImageMagick Library
  which is under Apache 2.0 license
*/

// These may be defined in some version of math.h
#undef I0
#undef J1
#undef P1
#undef Q1

static double I0(double x)
{
  /*
    Zeroth order Bessel function of the first kind.
  */
  double sum = 1.0;
  double y   = x*x/4.0;
  double t   = y;

  for (register int i = 2; t > JINC_EPSILON; i++) {
    sum += t;
    t   *= y / ((double) i*i);
  }

  return sum;
}

static double J1(double x)
{
  static const double
    Pone[] = {
       0.581199354001606143928050809e+21,
      -0.6672106568924916298020941484e+20,
       0.2316433580634002297931815435e+19,
      -0.3588817569910106050743641413e+17,
       0.2908795263834775409737601689e+15,
      -0.1322983480332126453125473247e+13,
       0.3413234182301700539091292655e+10,
      -0.4695753530642995859767162166e+7,
       0.270112271089232341485679099e+4
    },
    Qone[] = {
      0.11623987080032122878585294e+22,
      0.1185770712190320999837113348e+20,
      0.6092061398917521746105196863e+17,
      0.2081661221307607351240184229e+15,
      0.5243710262167649715406728642e+12,
      0.1013863514358673989967045588e+10,
      0.1501793594998585505921097578e+7,
      0.1606931573481487801970916749e+4,
      0.1e+1
    };

  double p = Pone[8];
  double q = Qone[8];
  
  for (register int i = 7; i >= 0; i--){
    p = p*x*x + Pone[i];
    q = q*x*x + Qone[i];
  }

  return p/q;
}

static double P1(double x)
{
  static const double
    Pone[] = {
      0.352246649133679798341724373e+5,
      0.62758845247161281269005675e+5,
      0.313539631109159574238669888e+5,
      0.49854832060594338434500455e+4,
      0.2111529182853962382105718e+3,
      0.12571716929145341558495e+1
    },
    Qone[] = {
      0.352246649133679798068390431e+5,
      0.626943469593560511888833731e+5,
      0.312404063819041039923015703e+5,
      0.4930396490181088979386097e+4,
      0.2030775189134759322293574e+3,
      0.1e+1
    };

  double p = Pone[5];
  double q = Qone[5];

  for (register int i = 4; i >= 0; i--) {
    p = p*(8.0/x)*(8.0/x) + Pone[i];
    q = q*(8.0/x)*(8.0/x) + Qone[i];
  }

  return p/q;
}

static double Q1(double x) {
  static const double
    Pone[] = {
      0.3511751914303552822533318e+3,
      0.7210391804904475039280863e+3,
      0.4259873011654442389886993e+3,
      0.831898957673850827325226e+2,
      0.45681716295512267064405e+1,
      0.3532840052740123642735e-1
    },
    Qone[] = {
      0.74917374171809127714519505e+4,
      0.154141773392650970499848051e+5,
      0.91522317015169922705904727e+4,
      0.18111867005523513506724158e+4,
      0.1038187585462133728776636e+3,
      0.1e+1
    };

  double p = Pone[5];
  double q = Qone[5];

  for (register int i = 4; i >= 0; i--) {
    p = p*(8.0/x)*(8.0/x) + Pone[i];
    q = q*(8.0/x)*(8.0/x) + Qone[i];
  }

  return p/q;
}

static double BesselOrderOne(double x)
{
  if (x < 1e-8)
    return 0.0;
  
  double p = x;

  if (x < 0.0)
    x = -x;

  if (x < 8.0)
    return p * J1(x);

  double q = (
      sqrt((double)(2.0/(JINC_PI*x)))
    * (P1(x)*(1.0/sqrt(2.0)*(sin((double)x)-cos((double)x)))-8.0/x*Q1(x)*(-1.0/sqrt(2.0)*(sin((double) x)+cos((double) x))))
  );
  
  if (p < 0.0)
    q = -q;

  return q;
}

static double Jinc(double x)
{
  if (x < 1e-8)
    return 0.5*JINC_PI;
  return BesselOrderOne(JINC_PI*x)/x;
}

JincFilter::JincFilter(int taps) :
  taps(taps)
{
  // EWA Lanczos function => jinc(x)*jinc(x/tap)
  // The tap size get converted to real value by following table
  // Window function (jinc(x/tap)) is actually jinc(x*jinc_zero[0]/jinc_zero[tap-1])

  static double jinc_zeros[16] = {
    1.2196698912665045,
    2.2331305943815286,
    3.2383154841662362,
    4.2410628637960699,
    5.2427643768701817,
    6.2439216898644877,
    7.2447598687199570,
    8.2453949139520427,
    9.2458926849494673,
    10.246293348754916,
    11.246622794877883,
    12.246898461138105,
    13.247132522181061,
    14.247333735806849,
    15.247508563037300,
    16.247661874700962
  };

  if (taps > 16) {
    this->taps = taps = 16;
  }

  support = jinc_zeros[taps-1];
  window_factor = jinc_zeros[0]/support;
}

float JincFilter::factor(float dist)
{
  dist = sqrt(dist);
  if (dist > support)
    return 0;
  return Jinc(dist) * Jinc(dist * window_factor);
}
