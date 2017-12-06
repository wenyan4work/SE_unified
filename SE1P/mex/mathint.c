#include "mathint.h"

#if (defined __AVX__ || defined __SSE4_2__)
#include "math_x86.h"
#endif


double GL_value[] = {
  0.0000875560264341,
  0.0004612700113121,
  0.0011333756872430,
  0.0021036207325094,
  0.0033714435498935,
  0.0049360907541328,
  0.0067966286377069,
  0.0089519457821408,
  0.0114007542680463,
  0.0141415906264317,
  0.0171728167840174,
  0.0204926210731500,
  0.0240990193293678,
  0.0279898560848899,
  0.0321628058610418,
  0.0366153745605260,
  0.0413449009595198,
  0.0463485582991216,
  0.0516233559754209,
  0.0571661413273014,
  0.0629736015209841,
  0.0690422655302257,
  0.0753685062110155,
  0.0819485424695466,
  0.0887784415221781,
  0.0958541212460431,
  0.1031713526189034,
  0.1107257622467940,
  0.1185128349779526,
  0.1265279166014690,
  0.1347662166290456,
  0.1432228111582063,
  0.1518926458152429,
  0.1607705387761404,
  0.1698511838636770,
  0.1791291537188462,
  0.1885989030447076,
  0.1982547719207257,
  0.2080909891856185,
  0.2181016758866909,
  0.2282808487935948,
  0.2386224239744122,
  0.2491202204319278,
  0.2597679637979140,
  0.2705592900832239,
  0.2814877494814479,
  0.2925468102238625,
  0.3037298624833663,
  0.3150302223250705,
  0.3264411357011823,
  0.3379557824877933,
  0.3495672805611614,
  0.3612686899110478,
  0.3730530167886529,
  0.3849132178866700,
  0.3968422045489604,
  0.4088328470073314,
  0.4208779786428876,
  0.4329704002694061,
  0.4451028844361781,
  0.4572681797477423,
  0.4694590151979302,
  0.4816681045156332,
  0.4938881505196921,
  0.5061118494803079,
  0.5183318954843668,
  0.5305409848020698,
  0.5427318202522577,
  0.5548971155638218,
  0.5670295997305939,
  0.5791220213571124,
  0.5911671529926685,
  0.6031577954510396,
  0.6150867821133300,
  0.6269469832113471,
  0.6387313100889522,
  0.6504327194388386,
  0.6620442175122067,
  0.6735588642988177,
  0.6849697776749295,
  0.6962701375166337,
  0.7074531897761375,
  0.7185122505185521,
  0.7294407099167761,
  0.7402320362020860,
  0.7508797795680722,
  0.7613775760255878,
  0.7717191512064052,
  0.7818983241133091,
  0.7919090108143816,
  0.8017452280792743,
  0.8114010969552925,
  0.8208708462811538,
  0.8301488161363231,
  0.8392294612238596,
  0.8481073541847571,
  0.8567771888417937,
  0.8652337833709545,
  0.8734720833985310,
  0.8814871650220474,
  0.8892742377532059,
  0.8968286473810967,
  0.9041458787539569,
  0.9112215584778219,
  0.9180514575304535,
  0.9246314937889845,
  0.9309577344697743,
  0.9370263984790159,
  0.9428338586726985,
  0.9483766440245791,
  0.9536514417008783,
  0.9586550990404803,
  0.9633846254394740,
  0.9678371941389582,
  0.9720101439151101,
  0.9759009806706322,
  0.9795073789268500,
  0.9828271832159826,
  0.9858584093735683,
  0.9885992457319537,
  0.9910480542178592,
  0.9932033713622931,
  0.9950639092458672,
  0.9966285564501065,
  0.9978963792674906,
  0.9988666243127571,
  0.9995387299886880,
  0.9999124439735659
};

double GL_weight[] = {
  0.0002246904801461,
  0.0005229063396701,
  0.0008212515093345,
  0.0011191442154813,
  0.0014163757357290,
  0.0017127630204551,
  0.0020081274918693,
  0.0023022921283514,
  0.0025950809163382,
  0.0028863187714329,
  0.0031758315808536,
  0.0034634462834494,
  0.0037489909628174,
  0.0040322949452430,
  0.0043131888993084,
  0.0045915049358304,
  0.0048670767075034,
  0.0051397395079161,
  0.0054093303697515,
  0.0056756881620402,
  0.0059386536863701,
  0.0061980697719754,
  0.0064537813696337,
  0.0067056356443082,
  0.0069534820664760,
  0.0071971725020834,
  0.0074365613010736,
  0.0076715053844326,
  0.0079018643296997,
  0.0081275004548926,
  0.0083482789007946,
  0.0085640677115557,
  0.0087747379135588,
  0.0089801635925043,
  0.0091802219686657,
  0.0093747934702723,
  0.0095637618049755,
  0.0097470140293533,
  0.0099244406164154,
  0.0100959355210650,
  0.0102613962434801,
  0.0104207238903756,
  0.0105738232341107,
  0.0107206027696042,
  0.0108609747690261,
  0.0109948553342302,
  0.0111221644468999,
  0.0112428260163725,
  0.0113567679251182,
  0.0114639220718434,
  0.0115642244121935,
  0.0116576149970314,
  0.0117440380082680,
  0.0118234417922238,
  0.0118957788905017,
  0.0119610060683517,
  0.0120190843405120,
  0.0120699789945096,
  0.0121136596114076,
  0.0121501000839859,
  0.0121792786323453,
  0.0122011778169248,
  0.0122157845489250,
  0.0122230900981313,
  0.0122230900981313,
  0.0122157845489250,
  0.0122011778169248,
  0.0121792786323453,
  0.0121501000839859,
  0.0121136596114076,
  0.0120699789945096,
  0.0120190843405120,
  0.0119610060683517,
  0.0118957788905017,
  0.0118234417922238,
  0.0117440380082680,
  0.0116576149970314,
  0.0115642244121935,
  0.0114639220718434,
  0.0113567679251182,
  0.0112428260163725,
  0.0111221644468999,
  0.0109948553342302,
  0.0108609747690261,
  0.0107206027696042,
  0.0105738232341107,
  0.0104207238903756,
  0.0102613962434801,
  0.0100959355210650,
  0.0099244406164154,
  0.0097470140293533,
  0.0095637618049755,
  0.0093747934702723,
  0.0091802219686657,
  0.0089801635925043,
  0.0087747379135588,
  0.0085640677115557,
  0.0083482789007946,
  0.0081275004548926,
  0.0079018643296997,
  0.0076715053844326,
  0.0074365613010736,
  0.0071971725020834,
  0.0069534820664760,
  0.0067056356443082,
  0.0064537813696337,
  0.0061980697719754,
  0.0059386536863701,
  0.0056756881620402,
  0.0054093303697515,
  0.0051397395079161,
  0.0048670767075034,
  0.0045915049358304,
  0.0043131888993084,
  0.0040322949452430,
  0.0037489909628174,
  0.0034634462834494,
  0.0031758315808536,
  0.0028863187714329,
  0.0025950809163382,
  0.0023022921283514,
  0.0020081274918693,
  0.0017127630204551,
  0.0014163757357290,
  0.0011191442154813,
  0.0008212515093345,
  0.0005229063396701,
  0.0002246904801461
};

/* Incomplete Bessel with Adaptive Simpson rule*/
double f(double x, double bes_a, double bes_b, int der)
{
  if(x==0)
    return 0;
  if(der==1)
     return exp(-bes_a/x-bes_b*x);
  else if(der==0)
     return exp(-bes_a/x-bes_b*x)/x;
  return 0;
}

double
simpson(double a, double b, double bes_a, double bes_b, int der)
{
  double c = (a+b)/2.;
  return (b-a)/6.*(f(a, bes_a, bes_b, der)+4.*f(c, bes_a, bes_b, der)+f(b, bes_a, bes_b, der));
}

double
errest(double a, double b, double *val, double bes_a, double bes_b, int der)
{
  double c = (a+b)/2.;
  double I = simpson(a,b, bes_a, bes_b, der);
  *val      = simpson(a,c, bes_a, bes_b, der) + simpson(c,b, bes_a, bes_b, der);
  return fabs(I- (*val) );
}

double
IncompBesselK0_Simpson(double tol, int *cnt, double bes_a, double bes_b, int der)
{

  int HEAD, n;
  double res = 0.0, val;
  double a, b;
  /* FIXME: size of stack should be dynamic. One possibility is 
   * is to use vectors. Otherwise one can keep only the last value 
   * of the interval and if OK reset it to c2 again! (expensive)
   */
  double *ST = (double*) malloc(STACK_SIZE*sizeof(double));

  // First value in the stack is upper bound of the integral
  ST[0] = 1.;
  
  if( 15.*errest(0.,1.,&val, bes_a, bes_b,der)<tol)
    return val;
  a = 0.;
  b = 1.;
  
  HEAD  = 0;
  
  n=0;  // to count the number of splits
  for(;;)
    {
      if(15.*errest(a,b,&val, bes_a, bes_b,der)>tol)
        {
          if ( HEAD>=(STACK_SIZE-1) )
            {
              printf("Increase the stack size!\n");
              exit(EXIT_FAILURE);
            }
          ST[++HEAD] = (a+b)/2.;
          b          = ST[HEAD];
          n++;
        }
      else
        {
          res += val;
          // if a==b or HEAD is pointing the the first element we break
          if(HEAD == 0)
            break;
          a    = ST[HEAD];
          b    = ST[--HEAD];
        }
    }
  free(ST);
  *cnt = n;
  return res;
}

/* kernels for computing the modified bessel function of the 
 * second kind and its derivative */
double bessel_f(double x, void * p)
{
  gsl_params * params = (gsl_params *)p;
  double a = (params->a);
  double b = (params->b);
  return exp(-a/x-b*x)/x;
}

double bessel_f_der(double x, void * p)
{
  gsl_params * params = (gsl_params *)p;
  double a = (params->a);
  double b = (params->b);
  return exp(-a/x-b*x);
}

double call_gsl_bessel_integrator(double a, double b, 
				  gsl_integration_workspace *w,
				  int der)
{
  double result, error;
  gsl_function F;
  if(der==0)
    F.function = &bessel_f;
  else if(der==1)
    F.function = &bessel_f_der;

  if(a==0){
    result = 1.e308;
    return result;
  }

  if(a<b && 0){
    gsl_params  params = {b,a};
    F.params = &params;
    gsl_integration_qags (&F, 0, 1, 0, 1e-12, 10000,
     			  w, &result, &error);
    double z = 2.*gsl_sf_bessel_K0(2.*sqrt(a*b));
    result = z-result;
  }
  else{
    gsl_params  params = {a,b};
    F.params = &params;
    gsl_integration_qags (&F, 0, 1, 0, 1e-12, 10000, 
			  w, &result, &error); 
  }
  
  return result;   
}


  
double func(double x, double a, double b)
{
  if(x==0)
    return 0;
  else
    return exp(-a/x-b*x)/x;
}

double computeK0(double a, double b) {
    // integration over [0,1] with "panels" panels
  double s=0; 
  int panels = 4;
  double x0;
  double h = 1./panels;
  for (int p=0;p<panels;p++) {
    x0 = 0+p*h;
#ifdef _OPENMP
#pragma omp parallel for  reduction(+:s)
#endif
    for (int i=0; i<128;i++) {
      double x = x0 + h*GL_value[i];
      s += GL_weight[i]*func(x,a,b)*h;
    }
  }
  return s;
}

double computeINCBK0(double a, double b, int der) {
    // integration over [0,1] with "panels" panels
  double s=0; 
  int panels = 4;
  double x0;
  double h = 1./panels;
  for (int p=0;p<panels;p++) {
    x0 = 0+p*h;
/* #ifdef _OPENMP */
/* #pragma omp parallel for  reduction(+:s) */
/* #endif */
    for (int i=0; i<128;i++) {
      double x = x0 + h*GL_value[i];
      s += GL_weight[i]*f(x,a,b,der)*h;
    }
  }
  return s;
}
