#include "surround.h"
void surroundtype::init_grid()
{
	box.erase(box.begin(),box.end());
	box.resize(nx);
	izmax.resize(nx);
	for(int i=0;i<nx;i++) {box[i].resize(ny);izmax[i].resize(ny);}
	for(int i=0;i<nx;i++) 
	for(int j=0;j<ny;j++)
	{
		box[i][j].resize(nz);
		izmax[i][j]=0;
	}
}

