#ifndef __VECTOR_MATH_H__
#define __VECTOR_MATH_H__

	
template <typename T>
double emp_var(std::vector<T>* vec)
{
	if (vec->size()<=1)
		return 0.0; 

	double mean = 0.0; 
	for (size_t i=0; i<vec->size(); i++)
		mean += vec->at(i)/vec->size(); 
	
	double var = 0.0; 
	for (size_t i=0; i<vec->size(); i++)
		var += ((vec->at(i)-mean)*(vec->at(i)-mean))/(vec->size()-1); 

	return var; 
}

//template <typename T>
//T my_prctile(std::vector<T>* vec, float pos, bool is_sorted)
//{
//	if (!is_sorted)
//		sort(vec->begin(), vec->end());
//
//	int max_pos = ((int) vec->size())-1;
//	int idx1 = floor(max_pos*pos);
//	int idx2 = ceil(max_pos*pos);
//
//	assert(idx1>=0);
//	assert(idx2>=0);
//	assert(idx1<=max_pos);
//	assert(idx2<=max_pos);
//	
//	return ((vec->at(idx1)+vec->at(idx2))/2);
//}


template <typename T>                                                                                                                             
T my_prctile(std::vector<T>* vec, float pos, bool is_sorted)                                                                 
{                                                                                                                       
    if (!is_sorted)                                                                                                                                               
        sort(vec->begin(), vec->end());                                                                                 
                                                                                                                                               
    int max_pos = ((int) vec->size())-1;                                                                                                                                    
                                                                                                                                                                  
    // assert that if pos is larger than one, then only due to minor numerical issues                                    
    assert(pos<1.0+1e-3);                                                                                                                                   
    if (pos>=1.0)                                                                                                                                             
    {                                                                                                                                                                       
        return vec->at(max_pos);                                                                                                                         
    }                                                                                                                                             
                                                                                                                        
    int idx1 = floor(max_pos*pos);                                                                                      
    int idx2 = ceil(max_pos*pos);                                                                                                                           
                                                                                                                                          
    assert(idx1>=0);                                                                                                    
    assert(idx2>=0);                                                                                                            
    assert(idx1<=max_pos);                                                                                                                           
    assert(idx2<=max_pos);                                                                                                                                              
                                                                                                                                                                  
    return ((vec->at(idx1)+vec->at(idx2))/2);                                                                           
}               


template <typename T>
T my_sum(std::vector<T>* vec)
{
	T sum = 0;

	for (uint64_t i=0; i<vec->size(); i++)
		sum+=vec->at(i); 
	
	return sum;
}


template <typename T>
T my_prctile(std::vector<T>* vec, float pos)
{
	return my_prctile(vec, pos, true); 
}

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
}

template <typename T>
std::vector<T> prctile(std::vector<T>* vec, int num_steps)
{
	vector<float> pos; 
	for (int i=0; i<num_steps; i++)
	{
		pos.push_back(((float) i)/(num_steps-1)); 
	}

	sort(vec->begin(), vec->end());

	std::vector<T> ret;
	for (unsigned int i=0; i<pos.size(); i++)
	{
		ret.push_back(my_prctile(vec, pos[i])); 
	}
	return ret;
}


#endif
