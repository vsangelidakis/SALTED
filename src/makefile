INCLUDES = -I$(GSL_INCLUDE_PATH)
LIBRARIES = -L$(GSL_LIBRARY_PATH)

deposit: deposit.o cluster.o common.o functions1.o  numerics.o surround.o mod.o matrix.o rotinf.o 
	g++  -O3  -mtune=native $(INCLUDES) $(LIBRARIES) -Wl,--no-as-needed -lgsl -lgslcblas -lm deposit.o cluster.o common.o functions1.o  numerics.o surround.o mod.o matrix.o rotinf.o  -o deposit  -lrt
deposit.o: deposit.cpp common.h surround.h functions1.h cluster.h numerics.h
	g++   $(INCLUDES) -O3  -mtune=native -c deposit.cpp  -lrt
cluster.o: cluster.cpp cluster.h functions1.h 
	g++   $(INCLUDES) -O3  -mtune=native -c cluster.cpp  -lrt
common.o: common.cpp common.h
	g++  $(INCLUDES) -O3  -mtune=native -c common.cpp  -lrt
functions1.o: functions1.cpp functions1.h common.h 
	g++  $(INCLUDES) -O3  -mtune=native -c functions1.cpp  -lrt
numerics.o: numerics.cpp numerics.h
	g++   $(INCLUDES) -O3  -mtune=native -c numerics.cpp  -lrt
surround.o: surround.cpp surround.h
	g++   $(INCLUDES) -O3  -mtune=native -c surround.cpp  -lrt
mod.o: mod.cpp mod.h
	g++  $(INCLUDES) -O3  -mtune=native -c mod.cpp -lrt
matrix.o: matrix.cpp matrix.h
	g++  $(INCLUDES) -O3  -mtune=native -c matrix.cpp -lrt
rotinf.o: rotinf.cpp rotinf.h
	g++   $(INCLUDES) -O3  -mtune=native -c rotinf.cpp -lrt


