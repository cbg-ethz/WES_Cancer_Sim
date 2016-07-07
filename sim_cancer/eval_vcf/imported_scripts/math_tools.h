#ifndef _MATH_TOOLS_H__
#define _MATH_TOOLS_H__
#include <vector>
#include <assert.h>
#include <math.h>


double gammaln(double x);

//double factln(int n);

double factln(float n);

double beta_pdf(double x, double a, double b); 

double beta_pdf_deriv(double x, double a, double b); 

double beta_neg_log_lik_deriv(double x, double a, double b); 
double beta_sample(double a, double b); 

double gamma_sample(int k, double theta); 

void beta_mom(double mu, double var, double* a, double* b); 

double bino_pdf(int k, int n, float p);

float bino_cfd(int cnt, int num, float p);
//float bino_cdf(int cnt, int num, float p){ return bino_cfd(cnt, num, p);};

double norm_cdf(double x); 

double my_urand();

//template <typename T>
//T my_prctile(std::vector<T>* vec, float pos);

//template <typename T>
//std::vector<T> prctile(std::vector<T>* vec, std::vector<float> pos);

//template <typename T>
//T pearson(std::vector<T>* vec1, std::vector<T>* vec2);
// move code of template here, because otherwise I would need to define which 
// instanciations of the tempalte should go into math_tools.o
template <typename T>
T pearson(std::vector<T>* vec1, std::vector<T>* vec2)
{
	double s1 = 0; 
	double s2 = 0; 
	assert(vec1->size()==vec2->size()); 
	for (int i=0; i<vec1->size(); i++) 
		s1 += vec1->at(i); 
	for (int i=0; i<vec2->size(); i++) 
		s2 += vec2->at(i); 

	s1 /= vec1->size(); 
	s2 /= vec2->size(); 

	T num = 0.0; 
	T var1 = 0.0; 
	T var2 = 0.0; 
	for (int i=0; i<vec1->size(); i++)
	{
		num += (vec1->at(i)-s1)*(vec2->at(i)-s2); 
		var1 += (vec1->at(i)-s1)*(vec1->at(i)-s1);
		var2 += (vec2->at(i)-s2)*(vec2->at(i)-s2);
	}
	return num/(sqrt(var1)*sqrt(var2)); 
};

template <typename T>
T my_prctile(std::vector<T>* vec, float pos)
{

	int max_pos = ((int) vec->size())-1;
	int idx1 = floor(max_pos*pos);
	int idx2 = ceil(max_pos*pos);

	assert(idx1>=0);
	assert(idx2>=0);
	assert(idx1<=max_pos);
	assert(idx2<=max_pos);
	
	return ((vec->at(idx1)+vec->at(idx2))/2);
};

template <typename T>
std::vector<T> prctile(std::vector<T>* vec, std::vector<float> pos)
{
	sort(vec->begin(), vec->end());

	std::vector<T> ret;
	for (unsigned int i=0; i<pos.size(); i++)
	{
		ret.push_back(my_prctile(vec, pos[i])); 
	}
	return ret;
};


float mean(std::vector<float>* set);

#endif
