#ifndef _FUNCTION1_
#define _FUNCTION1_
#include "common.h"
#include "surround.h"
#include "cluster.h"

/* Type		function		description */
// void		initrand:		initalises random number generator (apart from seed)
// void		load_parameters:	reads file with input parameters
// vtype	cmass:			calculates the center of mass of the packing
// long double	mass:			calculates the mass of the packing assuming that sphere with radius 30 has unit mass
// int		checkinters:		calculates the number of overlaps in the sediment
// bool		checkinters:		calculates the number of overlaps between the cluster spheres and the sediment
// void		find_cluster_contacts:	finds number of contacts that the cluster has with its surroundings
// void		unload:			adds the cluster c to the fixed particles
// void		unloadsphere:		adds the cluster c to the fixed particles
// void		output:			outputs the fixed particles to a file
// void		output:			outputs a set of clusters to a single file
// void		writedetails:		outputs to a file the number of clusters sedimented (input variable size) and the number of intersections in the sediment (also input parameter)
// void		write2cm:		outputs to a file the filling height of the sediment

void initrand();
void load_parameters(const char* inputfile);

vtype cmass(const vector<vtype>& r, const vector<long double>& R);
long double mass(const vector<vtype>& r, const vector<long double>& R);

int checkinters(surroundtype& su);
bool checkinters(const clustertype& c,const surroundtype& s);

void unload(clustertype& c,surroundtype& su);
void unloadsphere(const vtype& v,int ii,long double R,surroundtype& su);

void find_cluster_contacts(long double eps_contact,const clustertype& c,const surroundtype& su,vector<vtype>& rps,vector<vtype>& rds,vector<int>& p,vector<int>& d);

void output(surroundtype& su,string fname);
void outputVTK(surroundtype& su,string fname);
void output(vector<clustertype>& cd,string fname);
void writedetails(int size,int intersect,string name);
void write2cm(long double z,string name);

#endif
