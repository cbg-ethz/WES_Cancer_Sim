#ifndef __STATS_TESTS_H__
#define __STATS_TESTS_H__
#include <math.h>
#include <assert.h>
#include <string>
    using std::string;
#include <string.h>
#include <stdio.h>
#include <vector>
    using std::vector;
#include <utility> // pair
	using std::pair; 
#include <stdlib.h>  
#include <algorithm> // sort
#include <iostream>
using std::cout;


/* Log gamma function
 * \log{\Gamma(z)}
 * AS245, 2nd algorithm, http://lib.stat.cmu.edu/apstat/245
 */
double kf_lgamma(double z)
{
	double x = 0;
	x += 0.1659470187408462e-06 / (z+7);
	x += 0.9934937113930748e-05 / (z+6);
	x -= 0.1385710331296526     / (z+5);
	x += 12.50734324009056      / (z+4);
	x -= 176.6150291498386      / (z+3);
	x += 771.3234287757674      / (z+2);
	x -= 1259.139216722289      / (z+1);
	x += 676.5203681218835      / z;
	x += 0.9999999999995183;
	return log(x) - 5.58106146679532777 - z + (z-0.5) * log(z+6.5);
}

/* complementary error function
 * \frac{2}{\sqrt{\pi}} \int_x^{\infty} e^{-t^2} dt
 * AS66, 2nd algorithm, http://lib.stat.cmu.edu/apstat/66
 */
double kf_erfc(double x)
{
	const double p0 = 220.2068679123761;
	const double p1 = 221.2135961699311;
	const double p2 = 112.0792914978709;
	const double p3 = 33.912866078383;
	const double p4 = 6.37396220353165;
	const double p5 = .7003830644436881;
	const double p6 = .03526249659989109;
	const double q0 = 440.4137358247522;
	const double q1 = 793.8265125199484;
	const double q2 = 637.3336333788311;
	const double q3 = 296.5642487796737;
	const double q4 = 86.78073220294608;
	const double q5 = 16.06417757920695;
	const double q6 = 1.755667163182642;
	const double q7 = .08838834764831844;
	double expntl, z, p;
	z = fabs(x) * M_SQRT2;
	if (z > 37.) return x > 0.? 0. : 2.;
	expntl = exp(z * z * - .5);
	if (z < 10. / M_SQRT2) // for small z
	    p = expntl * ((((((p6 * z + p5) * z + p4) * z + p3) * z + p2) * z + p1) * z + p0)
			/ (((((((q7 * z + q6) * z + q5) * z + q4) * z + q3) * z + q2) * z + q1) * z + q0);
	else p = expntl / 2.506628274631001 / (z + 1. / (z + 2. / (z + 3. / (z + 4. / (z + .65)))));
	return x > 0.? 2. * p : 2. * (1. - p);
}

/* The following computes regularized incomplete gamma functions.
 * Formulas are taken from Wiki, with additional input from Numerical
 * Recipes in C (for modified Lentz's algorithm) and AS245
 * (http://lib.stat.cmu.edu/apstat/245).
 *
 * A good online calculator is available at:
 *
 *   http://www.danielsoper.com/statcalc/calc23.aspx
 *
 * It calculates upper incomplete gamma function, which equals
 * kf_gammaq(s,z)*tgamma(s).
 */

#define KF_GAMMA_EPS 1e-14
#define KF_TINY 1e-290

// regularized lower incomplete gamma function, by series expansion
static double _kf_gammap(double s, double z)
{
	double sum, x;
	int k;
	for (k = 1, sum = x = 1.; k < 100; ++k) {
		sum += (x *= z / (s + k));
		if (x / sum < KF_GAMMA_EPS) break;
	}
	return exp(s * log(z) - z - kf_lgamma(s + 1.) + log(sum));
}

static double chi2_cdf(float k, float x)
{
	// use regularized incomplete gamma function
	return _kf_gammap(k/2, x/2); 
}

static double fisher_pval_method(vector<float>* pvals)
{
	double ret = 0.0;
	for (unsigned int i=0; i<pvals->size(); i++)
	{
		ret+=-2*log(pvals->at(i)); 
	}
	//printf("ret: %f\n", ret); 
	return 1-chi2_cdf(2*pvals->size(), ret); 
}

// regularized upper incomplete gamma function, by continued fraction
static double _kf_gammaq(double s, double z)
{
	int j;
	double C, D, f;
	f = 1. + z - s; C = f; D = 0.;
	// Modified Lentz's algorithm for computing continued fraction
	// See Numerical Recipes in C, 2nd edition, section 5.2
	for (j = 1; j < 100; ++j) {
		double a = j * (s - j), b = (j<<1) + 1 + z - s, d;
		D = b + a * D;
		if (D < KF_TINY) D = KF_TINY;
		C = b + a / C;
		if (C < KF_TINY) C = KF_TINY;
		D = 1. / D;
		d = C * D;
		f *= d;
		if (fabs(d - 1.) < KF_GAMMA_EPS) break;
	}
	return exp(s * log(z) - z - kf_lgamma(s) - log(f));
}

double kf_gammap(double s, double z)
{
	return z <= 1. || z < s? _kf_gammap(s, z) : 1. - _kf_gammaq(s, z);
}

double kf_gammaq(double s, double z)
{
	return z <= 1. || z < s? 1. - _kf_gammap(s, z) : _kf_gammaq(s, z);
}

/* Regularized incomplete beta function. The method is taken from
 * Numerical Recipe in C, 2nd edition, section 6.4. The following web
 * page calculates the incomplete beta function, which equals
 * kf_betai(a,b,x) * gamma(a) * gamma(b) / gamma(a+b):
 *
 *   http://www.danielsoper.com/statcalc/calc36.aspx
 */
static double kf_betai_aux(double a, double b, double x)
{
	double C, D, f;
	int j;
	if (x == 0.) return 0.;
	if (x == 1.) return 1.;
	f = 1.; C = f; D = 0.;
	// Modified Lentz's algorithm for computing continued fraction
	for (j = 1; j < 200; ++j) {
		double aa, d;
		int m = j>>1;
		aa = (j&1)? -(a + m) * (a + b + m) * x / ((a + 2*m) * (a + 2*m + 1))
			: m * (b - m) * x / ((a + 2*m - 1) * (a + 2*m));
		D = 1. + aa * D;
		if (D < KF_TINY) D = KF_TINY;
		C = 1. + aa / C;
		if (C < KF_TINY) C = KF_TINY;
		D = 1. / D;
		d = C * D;
		f *= d;
		if (fabs(d - 1.) < KF_GAMMA_EPS) break;
	}
	return exp(kf_lgamma(a+b) - kf_lgamma(a) - kf_lgamma(b) + a * log(x) + b * log(1.-x)) / a / f;
}
double kf_betai(double a, double b, double x)
{
	return x < (a + 1.) / (a + b + 2.)? kf_betai_aux(a, b, x) : 1. - kf_betai_aux(b, a, 1. - x);
}

template<typename T>
void print_vec(vector<T>* vec, unsigned int num, const char* format, const char* tag)
{
	printf("%s: ", tag); 
	for (unsigned int i=0; i<vec->size() && i<num ; i++)
		printf(format, vec->at(i)); 
	printf("\n"); 
}

bool my_compare_second(const pair<int,float>& p1, const pair<int,float>& p2)
{
	return p1.second<p2.second; 
}

void benjamini_hochberg(vector<float>* p_vals)
{
	int m=0; 
	vector<pair<int, float> > sort_map; 
	for (unsigned int i=0; i<p_vals->size(); i++)
	{
		if (isnan(p_vals->at(i)))
			continue; 
		m++; 
		sort_map.push_back(pair<int, float>(i, p_vals->at(i))); 	
	}
	sort(sort_map.begin(), sort_map.end(), my_compare_second); 

	printf("%s first:%f, last:%f\n", __func__, sort_map.front().second, sort_map.back().second); 

	int k=1; 
	for (unsigned int i=0; i<sort_map.size(); i++)
	{
		int idx = sort_map[i].first; 
		assert(idx<p_vals->size()); 
		assert(p_vals->at(idx) == sort_map[i].second); 

		//if (k<10)
		//	printf("benjamini_hochberg: %f -> %f (k:%i m:%i) \n", p_vals->at(idx), p_vals->at(idx)*m/k, k, m); 
		p_vals->at(idx) *= m/k; 
		k++; 
	}
}

float mean(vector<float>* set)
{
	if (set->size()==0)
	{
		printf("Error: mean of empty set\n"); 
		exit(1); 
	}
		
	double sum = 0; 
	for (unsigned int i=0; i<set->size(); i++)
		sum += set->at(i); 
	return sum/set->size(); 
}

double corrected_std(vector<float>* set)
{
	assert(set->size()>=2); 

	float mu = mean(set); 

	double sum = 0; 
	for (unsigned int i=0; i<set->size(); i++)
		sum += (mu-set->at(i))*(mu-set->at(i)); 
	
	return sqrt(sum/(set->size()-1)); 
}

float ttest(vector<float>* set1, vector<float>* set2)
{
	// assuming unpaired samples, potentially unequal sample size and variance

	int n1 = set1->size(); 
	int n2 = set2->size(); 
	if (n1<2 || n2<2)
	{
		printf("[%s] one set too small (%i<2 || %i<2), return 1.0)\n", __func__, n1, n2); 
		return 1.0; 
	}

	float mu1 = mean(set1); 
	float mu2 = mean(set2); 

	float std1 = corrected_std(set1); 
	float std2 = corrected_std(set2); 

	float var1 = std1*std1; 
	float var2 = std2*std2; 

	//printf("s1(%.3f %.3f %.3f mean: %.2f std:%.2f)\ns2(%.3f %.3f %.3f mean:%.2f std:%.2f)\n", set1->at(0), set1->at(1), set1->at(2), mu1, std1, set2->at(0), set2->at(1), set2->at(2), mu2, std2); 


	double s = sqrt(var1/n1 + var2/n2); 
	double t = abs((mu1 - mu2)/s); 


	float df = pow(var1/n1 + var2/n2, 2) / (pow(var1/n1,2)/(n1-1) + pow(var2/n2,2)/(n2-1)); 


	float x = df/(pow(t, 2)+df); 
	float res = 1-0.5*kf_betai(df/2, 0.5, x); 

	// this is for the special case with df==2
	//float res2 = 0.5+t/(2*sqrt(2+t*t)); 
	
	//if (res>0.9 )
	
	//printf("t: %.5f (s: %.3f) df: %.3f res: %.4f \n", t, s, df, res); 

	if (res>1)
	{
		printf("Error: res (%f) is larger than 1.0, this is an odd pvalue\n", res); 
		res = 1.0; 
	}
	//assert(res<=1); 
	if (res<1e-30)
	{
		printf("Error: res (%f) is smaller than 1e-30, this is an odd pvalue\n", res); 
		res = 1e-30;
	}
	//assert(res>=0); 

	if (res>0.5)
		res = 1-res; 

	return res; 
}

double area_under_curve(double* xy, int len, bool reversed)
{   
	assert(len>0 && xy);

	double area = 0.0;

	if (!reversed)
	{   
		for (int i=1; i<len; i++)
			area += 0.5*(xy[2*i]-xy[2*(i-1)])*(xy[2*i+1]+xy[2*(i-1)+1]);
	}
	else
	{   
		for (int i=1; i<len; i++)
			area += 0.5*(xy[2*i+1]-xy[2*(i-1)+1])*(xy[2*i]+xy[2*(i-1)]);
	}

	return area;
}

float compute_auPRC(vector<pair<float, int> >* PR, int total_ones)
{
	float tp = 0; 
	double* curve = new double[2*PR->size()]; 
	for (unsigned int i=0; i<PR->size(); i++)
	{
		tp += PR->at(i).second==1;
		// precision 
		curve[2*i] = tp/(i+1); 
		// recall
		curve[2*i+1] = tp/total_ones; 
	}
	bool reversed = true; 
	double auPRC = area_under_curve(curve, PR->size(), reversed); 
	delete[] curve; 
	return auPRC; 
}

double my_gammaln(double x)
{
    double y,tmp,ser;
    static double cof[6] = {76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5 };

    y=x;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser = 1.000000000190015;
    for (int j=0; j<=5; j++)
        ser+= cof[j]/++y;

    return -tmp+log(2.5066282746310005*ser/x);
}

double my_factln(float h)
{
	return my_gammaln(h+1.0); 
}

double fisher_exact(int a, int b, int c, int d)
{
	// table: 
	//
	//	a  |  b 
	//  --------
	//  c  |  d  
	//

	int n = a+b+c+d; 

	double ret = 0.0; 
	
	ret += my_factln(a+b); 
	ret += my_factln(c+d); 
	ret += my_factln(a+c); 
	ret += my_factln(b+d); 
	ret -= my_factln(a); 
	ret -= my_factln(b); 
	ret -= my_factln(c); 
	ret -= my_factln(d); 
	ret -= my_factln(n); 

	return exp(ret); 
}

double binom_pdf(int k, int n, float p)
{   
	if (p==0.0)
		return 1.0 *(k==0); 

    double ln_binco = my_factln(n)-my_factln(k)-my_factln(n-k);
    double ln_ret = ln_binco + k*log(p) + (n-k)*log(1-p);
	//printf("[%s] k:%i n:%i p:%.3f log(p):%f %f %f %f\n", __func__, k, n, p, log(p), ln_binco, ln_ret, exp(ln_ret)); 
    return exp(ln_ret);
}

float binom_cdf(int cnt, int num, float p, string greater_or_less)
{
	if (cnt>=num)
		return 1.0; 

    // compute cdf
    float ret = 0;
	if(greater_or_less.compare("less")==0){
	    for (int i=0; i<=cnt; i++)
	        ret+=binom_pdf(i, num, p);
	} else if (greater_or_less.compare("greater")==0) {
		if(cnt==num){
			cout << "cnt == sum!!\n";
			ret=0.0;
		} else {
			for (int i=cnt; i<=num; i++)
				ret+=binom_pdf(i, num, p);
		}
	} else {
		std::cerr << "error: neither greater nor less\n";
		return -1.0;
	}

	assert(ret<=1+1e-3);
	if (ret>1.0)
		ret = 1.0; 
	
    return ret;
}


#endif
