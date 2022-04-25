#ifndef _FUNCTION1_
#define _FUNCTION1_
#include "common.h"
#include "surround.h"
#include "cluster.h"

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
void output(vector<clustertype>& cd,string fname);
void writedetails(int size,int intersect,string name);
void write2cm(long double z,string name);

#endif
