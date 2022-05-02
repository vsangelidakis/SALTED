#ifndef _SURROUND_
#define _SURROUND_

#include "common.h"

class surroundtype	// FIXME: Are there two surroundtype classes, here and in cluster?
{
	public:
	vector<vtype> r;
	vector<long double> R;
	vector<vector<vector<vector<int> > > > box; // 3D matrix with size [nx*ny*nz]
	vector<vector<int> > izmax; // 2D matrix with size [nx*ny]
	void init_grid();
	vector<int> id;
};

#endif
