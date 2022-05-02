#include "surround.h"
void surroundtype::init_grid()
{
	box.erase(box.begin(),box.end());
	box.resize(nx); // Make box with size [nx*1]
	izmax.resize(nx);
	for(int i=0; i<nx; i++) {box[i].resize(ny); izmax[i].resize(ny);} // Resize box to [nx*ny]
	for(int i=0; i<nx; i++)
	for(int j=0; j<ny; j++)
	{
		box[i][j].resize(nz); // Resize box to [nx*ny*nz]
		izmax[i][j]=0;
	}
}

