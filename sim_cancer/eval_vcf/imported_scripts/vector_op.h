#ifndef _VECTOR_OP_H__
#define _VECTOR_OP_H__

#include <assert.h>

template<typename T>
T max(vector<T>* vec)
{
	assert(vec->size()>=1);
	typename vector<T>::iterator it;
	it = max_element (vec->begin(), vec->end());
	return *it;
}
template<typename T>
T min(vector<T>* vec)
{
	assert(vec->size()>=1);
	typename vector<T>::iterator it;
	it = min_element (vec->begin(), vec->end());
	return *it;
}

template<typename T>
void mult(vector<T>* vec, T val)
{
	for (int i=0; i<vec->size(); i++)
	{
		(vec->at(i))*=val;
	}
}

vector<int> range(int lb, int ub)
{
	vector<int> ret;
	for (int i=lb; i<=ub; i++)
		ret.push_back(i);
	
	return ret;
}
template<typename T>
vector<T> logspace(T lb, T ub, int num)
{
	vector<T> ret;
	T step = (log(ub)-log(lb))/num;
	assert(step>0);

	for (T f=log(lb); f<=log(ub); f+=step)
		ret.push_back(exp(f));
	
	return ret;
}

template<typename T> 
void print_vec(vector<T>* vec, const char* format_str)
{
	for (int i=0; i<vec->size(); i++)
	{
		printf(format_str, vec->at(i));
	}
	printf("\n");
}
template<typename T> 
void print_mat(vector<vector<T> >* mat, const char* format_str)
{
	for (int i=0; i<mat->size(); i++)
	{
		for (int j=0; j<mat->at(i).size(); j++)
		{
			printf(format_str, mat->at(i)[j]);
		}
		printf("\n");
	}
}
#endif
