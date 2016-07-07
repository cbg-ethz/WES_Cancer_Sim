#ifndef __GET_INT_H__
#define __GET_INT_H__

#include <mex.h>
#include <stdlib.h>

#include <stdio.h>

int get_int(const mxArray *prhs)
{
	if (!prhs)
		mexErrMsgTxt("get_int: called with NULL pointer arg");
	if (!mxIsDouble(prhs))
		mexErrMsgTxt("get_int: input is not a double");
	if (mxGetN(prhs)!=1 || mxGetM(prhs)!=1)
		 mexErrMsgTxt( "get_int: argument should be a scalar\n");

	double* arg = mxGetPr(prhs);	
	return (int) *arg;
};
double get_double(const mxArray *prhs)
{
	if (!prhs)
		mexErrMsgTxt("get_double: called with NULL pointer arg");
	if (!mxIsDouble(prhs))
		mexErrMsgTxt("get_double: input is not a double");
	if (mxGetN(prhs)!=1 || mxGetM(prhs)!=1)
		 mexErrMsgTxt( "get_double: argument should be a scalar\n");

	double* arg = mxGetPr(prhs);	
	return *arg;
};
char *get_string(const mxArray *prhs)
{
	char *buf;
	int buflen;

	if (!prhs)
		mexErrMsgTxt("get_string called with NULL pointer arg");
	if (!mxIsChar(prhs))
		mexErrMsgTxt("input is not a string");
	if (mxGetM(prhs) != 1)
		mexErrMsgTxt("input is not a row vector");
	buflen = mxGetN(prhs) + 1;
	buf = new char[buflen];
	/* copy the string from prhs into buf and add terminating NULL char */
	if (mxGetString(prhs, buf, buflen))
		mexErrMsgTxt("not enough space");
	return buf;
};
#endif
