#ifndef _SURROUND_
#define _SURROUND_

#include "common.h"

class surroundtype
{
	public:
	vector<vtype> r;
	vector<long double> R;
	vector<vector<vector<vector<int> > > > box;
	vector<vector<int> > izmax;
	void init_grid();
};

#endif
