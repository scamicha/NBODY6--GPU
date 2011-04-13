FC = gfortran
#FC = gfortran44
FFLAGS = -O3 -Wall
#FFLAGS -pg
CXX = g++ -O3 -Wall -march=core2
#CXX += pg
CUDA_PATH = /sw/keeneland/cuda/3.1/linux_binary
SDK_PATH = /nics/c/home/scamicha/NVIDIA_GPU_Computing_SDK/C
NVCC = $(CUDA_PATH)/bin/nvcc -Xcompiler "-Wall -O3 -fopenmp" -I $(SDK_PATH)/common/inc/
NVCC += -arch sm_21 -Xptxas -dlcm=cg
RUNDIR = ./run

OBJECTS = $(SOURCE:.f=.o)

all: gpu sse

gpu: gpunb.gpu.o gpupot.gpu.o intgrt.o gpucor.o nbintp.o cnbint.o nbint.o start.o fpoly2a.o adjust.o energy2.o phicor.o cmfirr.o cmfirr2.o checkl2.o kspert.o swap.o scale.o wtime.o cxvpred.o gpuirr.o jpred.o fpcorr2.o repair.o
	cp -f gpunb.gpu.o gpunb.o
	cp -f gpupot.gpu.o gpupot.o
	cp -f *.o ../Ncode
	make gpu -C ../Ncode -f ../GPU2/Makefile_gpu
	mv -f ../Ncode/nbody6 $(RUNDIR)/nbody6.gpu

sse: gpunb.sse.o gpupot.sse.o intgrt.o gpucor.o nbintp.o cnbint.o nbint.o start.o fpoly2a.o adjust.o energy2.o phicor.o cmfirr.o cmfirr2.o checkl2.o kspert.o swap.o scale.o wtime.o cxvpred.o gpuirr.o jpred.o fpcorr2.o repair.o
	cp -f gpunb.sse.o gpunb.o
	cp -f gpupot.sse.o gpupot.o
	cp -f *.o ../Ncode
	make sse -C ../Ncode -f ../GPU2/Makefile_gpu -j8
	mv -f ../Ncode/nbody6 $(RUNDIR)/nbody6.sse

cnbint.o: ./lib/cnbint.cpp
	$(CXX) -c $^

wtime.o: ./lib/wtime.cpp
	$(CXX) -c $^

cxvpred.o: ./lib/cxvpred.cpp
	$(CXX) -fopenmp -c $^

gpuirr.o: ./irrlib/gpuirr.sse.cpp
	$(CXX) -fopenmp -c $^ -o $@

gpuirr.s: ./irrlib/gpuirr.sse.cpp
	$(CXX) -fopenmp -S $^ -o $@

intgrt.o: intgrt.omp.f
	$(FC) $(FFLAGS) -fopenmp $^ -c -o $@ 

intgrt.s: intgrt.omp.f
	$(FC) $(FFLAGS) -fopenmp $^ -S -o $@

start.o: start.omp.F
	$(FC) $(FFLAGS) -fopenmp $^ -c -o $@

gpunb.gpu.o: lib/gpunb.reduce.cu
	$(NVCC) $^ -c -o $@

gpupot.gpu.o: lib/gpupot.gpu.cu
	$(NVCC) $^ -c -o $@

gpunb.sse.o: lib/gpunb.sse.cpp
	$(CXX) -fopenmp $^ -c -o $@

gpupot.sse.o: lib/gpupot.sse.cpp
	$(CXX) -fopenmp $^ -c -o $@

clean:
	rm -f *.o

archive:
	tar -czvf ../gpu2.tar.gz ../GPU2	
