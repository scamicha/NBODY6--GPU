// #include <iostream>
#include <cstdio>
// #include <vector>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <cutil.h>
#include <omp.h>
#include "cuda_pointer.h"

#define NTHREAD 64 // 64, 96, 128 or 192
// #define NJBLOCK 16 // 8800GTS/512 has 16
#define NJBLOCK 14 // for GTX 470
#define NIBLOCK 16 // 16 or 32 
#define NIMAX (NTHREAD * NIBLOCK) // 1024

#define NBMAX 64 // NNB per block, must be power of 2

#define MAX_CPU 8
#define MAX_GPU 4

template <class T>
struct myvector{
	int num;
	T *val;
	myvector(){
		num = 0;
		val = NULL;
	}
	~myvector(){
		delete [] val;
	}
	void clear(){
		num = 0;
	}
	void reserve(size_t count){
		val = new T[count];
	}
	void free(){
		delete [] val;
	}
	void push_back(const T &t){
		val[num++] = t;
	}
	size_t size(){
		return num;
	}
	void resize(size_t size){
		num = size;
	}
	T &operator[](int i){
		return val[i];
	}
};

#define PROFILE
#ifdef PROFILE
#include <sys/time.h>
static double get_wtime(){
#if 1
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.e-6 * tv.tv_usec;
#else
	struct timespec tv;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tv);
	return tv.tv_sec + 1.e-9 * tv.tv_nsec;
#endif
}
#else
static double get_wtime(){
	return 0.0;
}
#endif

static double time_send, time_grav, time_reduce;
static long long numInter;

struct Jparticle{
	float3 pos;
	float  mass;
	float3 vel;
	float  pad;
	Jparticle() {}
	Jparticle(double mj, double xj[3], double vj[3]){
		pos.x = xj[0];
		pos.y = xj[1];
		pos.z = xj[2];
		mass  = mj;
		vel.x = vj[0];
		vel.y = vj[1];
		vel.z = vj[2];
	}
};
struct Iparticle{
	float3 pos;
	float  h2;
	float3 vel;
	float  pad;
	Iparticle() {}
	Iparticle(double h2i, double xi[3], double vi[3]){
		pos.x = xi[0];
		pos.y = xi[1];
		pos.z = xi[2];
		h2    = h2i;
		vel.x = vi[0];
		vel.y = vi[1];
		vel.z = vi[2];
	}
};
struct Force{
	float3 acc;
	float  pot;
	float3 jrk;
	int    nnb;          //  8 words
	unsigned short  neib[NBMAX]; // 24 words
	__device__  Force(){
		acc.x = acc.y = acc.z = 0.f;
		jrk.x = jrk.y = jrk.z = 0.f;
		pot = 0.f;
		nnb = 0;
	}
};

__device__ void h4_kernel(
		const int j,
		const Iparticle &ip, 
		const Jparticle &jp, 
		Force &fo,
		float3 &acc,
		float3 &jrk,
		int &nnb){
	float dx = jp.pos.x - ip.pos.x;
	float dy = jp.pos.y - ip.pos.y;
	float dz = jp.pos.z - ip.pos.z;
	float dvx = jp.vel.x - ip.vel.x;
	float dvy = jp.vel.y - ip.vel.y;
	float dvz = jp.vel.z - ip.vel.z;

	float r2 = dx*dx + dy*dy + dz*dz;
	float rv = dx*dvx + dy*dvy + dz*dvz;
	float rinv1 = rsqrtf(r2);
	if(r2 < ip.h2){
		// fo.neib[fo.nnb++ % NBMAX] = j;
		fo.neib[nnb & (NBMAX-1)] = (unsigned)j;
		nnb++;
		rinv1 = 0.f;
	}
	float rinv2 = rinv1 * rinv1;
	float mrinv1 = jp.mass * rinv1;
	float mrinv3 = mrinv1 * rinv2;
	rv *= -3.f * rinv2;
	
#ifdef POTENTIAL
	fo.pot += mrinv1;
#endif
	acc.x += mrinv3 * dx;
	acc.y += mrinv3 * dy;
	acc.z += mrinv3 * dz;
	// fo.acc.z += 1.0;
	jrk.x += mrinv3 * (dvx + rv * dx);
	jrk.y += mrinv3 * (dvy + rv * dy);
	jrk.z += mrinv3 * (dvz + rv * dz);
}

#if 0
__device__ void h4_grav_kernel(
		const int j,
		const Iparticle &ip, 
		const Jparticle &jp, 
		float3 &acc,
		float3 &jrk){
	float dx = jp.pos.x - ip.pos.x;
	float dy = jp.pos.y - ip.pos.y;
	float dz = jp.pos.z - ip.pos.z;
	float dvx = jp.vel.x - ip.vel.x;
	float dvy = jp.vel.y - ip.vel.y;
	float dvz = jp.vel.z - ip.vel.z;

	float r2 = dx*dx + dy*dy + dz*dz;
	float rv = dx*dvx + dy*dvy + dz*dvz;
	float rinv1 = rsqrtf(r2);
	if(r2 < ip.h2){
		rinv1 = 0.f;
	}
	float rinv2 = rinv1 * rinv1;
	float mrinv1 = jp.mass * rinv1;
	float mrinv3 = mrinv1 * rinv2;
	rv *= -3.f * rinv2;
	
#ifdef POTENTIAL
	fo.pot += mrinv1;
#endif
	acc.x += mrinv3 * dx;
	acc.y += mrinv3 * dy;
	acc.z += mrinv3 * dz;
	// fo.acc.z += 1.0;
	jrk.x += mrinv3 * (dvx + rv * dx);
	jrk.y += mrinv3 * (dvy + rv * dy);
	jrk.z += mrinv3 * (dvz + rv * dz);
}

__device__ void h4_neib_kernel(
		const int j,
		const Iparticle &ip, 
		const Jparticle &jp, 
		Force &fo,
		int &nnb){
	float dx = jp.pos.x - ip.pos.x;
	float dy = jp.pos.y - ip.pos.y;
	float dz = jp.pos.z - ip.pos.z;

	float r2 = dx*dx + dy*dy + dz*dz;
	if(r2 < ip.h2){
		// fo.neib[fo.nnb++ % NBMAX] = j;
		fo.neib[nnb & (NBMAX-1)] = (unsigned)j;
		nnb++;
	}
}
#endif

__global__ void h4_gravity(
		int nbody,
		Iparticle ipbuf[],
		Jparticle jpbuf[],
		Force fobuf[][NJBLOCK]){
	int ibid = blockIdx.x;
	int jbid = blockIdx.y;
	int tid = threadIdx.x;
	int iaddr = tid + NTHREAD * ibid;
	int jstart = (nbody * (jbid  )) / NJBLOCK;
	int jend   = (nbody * (jbid+1)) / NJBLOCK;

	Iparticle ip = ipbuf[iaddr];
	// Force fo;
	Force &fo = fobuf[iaddr][jbid];
	float3 acc = make_float3(0.0f, 0.0f, 0.0f);
	float3 jrk = make_float3(0.0f, 0.0f, 0.0f);
	int nnb = 0;
	for(int j=jstart; j<jend; j+=NTHREAD){
		__shared__ Jparticle jpshare[NTHREAD];
		__syncthreads();
#if 0
		jpshare[tid] = jpbuf[j+tid];
#else
		float4 *src = (float4 *)&jpbuf[j];
		float4 *dst = (float4 *)jpshare;
		dst[        tid] = src[        tid];
		dst[NTHREAD+tid] = src[NTHREAD+tid];
#endif
		__syncthreads();

#if 1
		if(jend-j < NTHREAD){
#pragma unroll 4
			for(int jj=0; jj<jend-j; jj++){
				Jparticle jp = jpshare[jj];
				h4_kernel(j+jj, ip, jp, fo, acc, jrk, nnb);
			}
		}else{
#pragma unroll 4
			for(int jj=0; jj<NTHREAD; jj++){
				Jparticle jp = jpshare[jj];
				h4_kernel(j+jj, ip, jp, fo, acc, jrk, nnb);
			}
		}
#else
		if(jend-j < NTHREAD){
#pragma unroll 4
			for(int jj=0; jj<jend-j; jj++){
				Jparticle jp = jpshare[jj];
				h4_grav_kernel(j+jj, ip, jp, acc, jrk);
			}
#pragma unroll 4
			for(int jj=0; jj<jend-j; jj++){
				Jparticle jp = jpshare[jj];
				h4_neib_kernel(j+jj, ip, jp, fo, nnb);
			}
		}else{
#pragma unroll 4
			for(int jj=0; jj<NTHREAD; jj++){
				Jparticle jp = jpshare[jj];
				h4_grav_kernel(j+jj, ip, jp, acc, jrk);
			}
#pragma unroll 4
			for(int jj=0; jj<NTHREAD; jj++){
				Jparticle jp = jpshare[jj];
				h4_neib_kernel(j+jj, ip, jp, fo, nnb);
			}
		}
#endif
	}
	fo.acc = acc;
	fo.jrk = jrk;
	fo.nnb = nnb;
	// fobuf[iaddr][jbid] = fo;
}

static cudaPointer <Jparticle> jpbuf[MAX_GPU];
static cudaPointer <Iparticle> ipbuf[MAX_GPU];
static cudaPointer <Force[NJBLOCK]> fobuf[MAX_GPU];
static int numCPU, numGPU;
static int joff[MAX_GPU + 1];
static myvector<int> nblist[MAX_CPU];
static int nbody, nbodymax;
static int device_id[MAX_GPU];
// static int *nblist;
static bool is_open = false;
static bool devinit = false;

void GPUNB_devinit(){
	if(devinit) return;

	cudaGetDeviceCount(&numGPU);
	assert(numGPU <= MAX_GPU);
	char *gpu_list = getenv("GPU_LIST");
	if(gpu_list){
		// get GPU list from environment variable
		numGPU = 0;
		char *p = strtok(gpu_list, " ");
		while(p){
			device_id[numGPU++] = atoi(p);
			p = strtok(NULL, " ");
			assert(numGPU <= MAX_GPU);
		}
	}else{
		// use all GPUs
		for(int i=0; i<numGPU; i++){
			device_id[i] = i;
		}
	}
	
	// numGPU = 1;
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		if(tid == 0) numCPU = omp_get_num_threads();
	}
	assert(numCPU <= MAX_CPU);
	assert(numGPU <= numCPU);
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		if(tid < numGPU){
			cudaSetDevice(device_id[tid]);
		}
	}
#ifdef PROFILE
	fprintf(stderr, "***********************\n");
	fprintf(stderr, "Initializing NBODY6/GPU library\n");
	fprintf(stderr, "#CPU %d, #GPU %d\n", numCPU, numGPU);
	fprintf(stderr, " device:");
	for(int i=0; i<numGPU; i++){
		fprintf(stderr, " %d", device_id[i]);
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "***********************\n");
#endif
	devinit = true;
}

void GPUNB_open(int nbmax){
	time_send = time_grav = time_reduce = 0.0;
	numInter = 0;
	nbodymax = nbmax;

	if(is_open){
		fprintf(stderr, "gpunb: it is already open\n");
		return;
	}
	is_open = true;

#if 0
	cudaGetDeviceCount(&numGPU);
	assert(numGPU <= MAX_GPU);
	char *gpu_list = getenv("GPU_LIST");
	if(gpu_list){
		// get GPU list from environment variable
		numGPU = 0;
		char *p = strtok(gpu_list, " ");
		while(p){
			device_id[numGPU++] = atoi(p);
			p = strtok(NULL, " ");
			assert(numGPU <= MAX_GPU);
		}
	}else{
		// use all GPUs
		for(int i=0; i<numGPU; i++){
			device_id[i] = i;
		}
	}
	
	// numGPU = 1;
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		if(tid == 0) numCPU = omp_get_num_threads();
	}
	assert(numCPU <= MAX_CPU);
	assert(numGPU <= numCPU);
#else
	GPUNB_devinit();
#endif

	for(int id=0; id<numGPU + 1; id++){
		joff[id] = (id * nbmax) / numGPU;
	}

	// omp_set_num_threads(numGPU);
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		if(tid < numGPU){
			cudaSetDevice(device_id[tid]);
			int nj = joff[tid+1] - joff[tid];
			jpbuf[tid].allocate(nj + NTHREAD);
			ipbuf[tid].allocate(NIMAX);
			fobuf[tid].allocate(NIMAX);
		}
	}
	// omp_set_num_threads(numCPU);
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		nblist[tid].reserve(nbmax);
	}
#ifdef PROFILE
	fprintf(stderr, "***********************\n");
	fprintf(stderr, "Opened NBODY6/GPU library\n");
	fprintf(stderr, "#CPU %d, #GPU %d\n", numCPU, numGPU);
	fprintf(stderr, " device:");
	for(int i=0; i<numGPU; i++){
		fprintf(stderr, " %d", device_id[i]);
	}
	fprintf(stderr, "\n");
	for(int i=0; i<numGPU+1; i++){
		fprintf(stderr, " %d", joff[i]);
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "nbmax = %d\n", nbmax);
	fprintf(stderr, "***********************\n");
#endif
}

void GPUNB_close(){
	if(!is_open){
		fprintf(stderr, "gpunb: it is already close\n");
		return;
	}
	is_open = false;
	// omp_set_num_threads(numGPU);
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		if(tid < numGPU){
			jpbuf[tid].free();
			ipbuf[tid].free();
			fobuf[tid].free();
		}
	}
	// omp_set_num_threads(numCPU);
	nbodymax = 0;

#ifdef PROFILE
# if 0
	std::cerr << "***********************" << std::endl;
	std::cerr << "time send : " << time_send << " sec " << std::endl;
	std::cerr << "time grav : " << time_grav << " sec " << std::endl;
	std::cerr << 60.e-9 * numInter / time_grav << " Gflops (gravity part only)" << std::endl;
	std::cerr << "***********************" << std::endl;
# else
	fprintf(stderr, "***********************\n");
	fprintf(stderr, "time send   : %f sec\n", time_send);
	fprintf(stderr, "time grav   : %f sec\n", time_grav);
	fprintf(stderr, "time reduce : %f sec\n", time_reduce);
	fprintf(stderr, "%f Gflops (gravity part only)\n", 60.e-9 * numInter / time_grav);
	fprintf(stderr, "***********************\n");
# endif
#endif
}

void GPUNB_send(
		int _nbody,
		double mj[],
		double xj[][3],
		double vj[][3]){
	nbody = _nbody;
	assert(nbody <= nbodymax);
	time_send -= get_wtime();
	for(int id=0; id<numGPU + 1; id++){
		joff[id] = (id * nbody) / numGPU;
	}
	// omp_set_num_threads(numGPU);
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		if(tid < numGPU){
			int nj = joff[tid+1] - joff[tid];
			// fprintf(stderr, "%d : %d\n", tid, nj);
			for(int j=0; j<nj; j++){
				int jj = j + joff[tid];
				jpbuf[tid][j] = Jparticle(mj[jj], xj[jj], vj[jj]);
			}
			jpbuf[tid].htod(nj);
		}
	}
	// size_t jpsize = nj * sizeof(Jparticle);
	// cudaMemcpy(jp_dev, jp_host, jpsize, cudaMemcpyHostToDevice);
	time_send += get_wtime();
	// omp_set_num_threads(numCPU);
}

static void handle_overflow(
		myvector<int> &list,
		cudaPointer <Iparticle> &ip,
		cudaPointer <Jparticle> &jp,
		int nj,
		int i,
		int jb,
		int nnb,
		int joff
		){
	fprintf(stderr, "gpunb overflow: %d %d %d\n", i, jb, nnb);
	int jstart = (nj * (jb  )) / NJBLOCK;
	int jend   = (nj * (jb+1)) / NJBLOCK;

	float xi = ip[i].pos.x;
	float yi = ip[i].pos.y;
	float zi = ip[i].pos.z;
	float h2i = ip[i].h2;
	int n = 0;
	for(int j=jstart; j<jend; j++){
		float dx = jp[j].pos.x - xi;
		float dy = jp[j].pos.y - yi;
		float dz = jp[j].pos.z - zi;
		float r2 = dx*dx + dy*dy + dz*dz;
		if(r2 < h2i){
			list.push_back(j + joff);
			n++;
		}
	}
	// assert(n == nnb);
	if(n != nnb){
		fprintf(stderr, "warning, NNB_GPU != NNB_CPU\n");
		// list.resize(list.size() - n);
		list.num -= n;
		if(n > nnb){
			ip[i].h2 *= .99999988079071044922f;
		}
		if(n < nnb){
			ip[i].h2 *= 1.00000011920928955078f;
		}
		handle_overflow(list, ip, jp, nj, i, jb, nnb, joff);
	}
}

void GPUNB_regf(
		int ni,
		double h2[],
		double xi[][3],
		double vi[][3],
		double acc[][3],
		double jrk[][3],
		double pot[],
		int lmax,
		int nbmax,
		int *listbase){
	time_grav -= get_wtime();
	numInter += ni * nbody;
	assert(0 < ni && ni <= NIMAX);

	// omp_set_num_threads(numGPU);
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		if(tid < numGPU){
			// cudaSetDevice(device_id[tid]);
			int dev;
			cudaGetDevice(&dev);
			assert(dev == device_id[tid]);
			for(int i=0; i<ni; i++){
				ipbuf[tid][i] = Iparticle(h2[i], xi[i], vi[i]);
			}
			// set i-particles
			ipbuf[tid].htod(ni);

			// gravity kernel
			int niblock = 1 + (ni-1) / NTHREAD;
			dim3 grid(niblock, NJBLOCK, 1);
			dim3 threads(NTHREAD, 1, 1);
			int nj = joff[tid+1] - joff[tid];
			h4_gravity <<< grid, threads >>> 
				(nj, ipbuf[tid], jpbuf[tid], fobuf[tid]);

			// recieve force
			// fprintf(stderr, "DBG %d %d\n", tid, device_id[tid]);
			fobuf[tid].dtoh(ni);
		}
	}

	const double wt = get_wtime();
	time_grav   += wt;
	time_reduce -= wt;

	// reduction phase
	// omp_set_num_threads(numCPU);
#pragma omp parallel for
	for(int i=0; i<ni; i++){
		int tid = omp_get_thread_num();
		double ax=0, ay=0, az=0;
		double jx=0, jy=0, jz=0;
#ifdef POTENTIAL
		double poti=0;
#endif
		for(int id=0; id<numGPU; id++){
			for(int jb=0; jb<NJBLOCK; jb++){
				// Force &fo = fo_host[i][jb];
				Force &fo = fobuf[id][i][jb];
				ax += fo.acc.x;
				ay += fo.acc.y;
				az += fo.acc.z;
				jx += fo.jrk.x;
				jy += fo.jrk.y;
				jz += fo.jrk.z;
#ifdef POTENTIAL
				poti += fo.pot;
#endif
			}
		}
		acc[i][0] = ax;
		acc[i][1] = ay;
		acc[i][2] = az;
		jrk[i][0] = jx;
		jrk[i][1] = jy;
		jrk[i][2] = jz;
		// fprintf(stderr, "%f %f %f %f %f %f\n", ax, ay, az, jx, jy, jz);
		// exit(0);
#ifdef POTENTIAL
		pot[i] = poti;
#endif
		bool overflow = false;
		nblist[tid].clear();
		for(int id=0; id<numGPU; id++){
			for(int jb=0; jb<NJBLOCK; jb++){
				// Force &fo = fo_host[i][jb];
				Force &fo = fobuf[id][i][jb];
				int nj = joff[id+1] - joff[id];
				// int jstart = (nbody * jb) / NJBLOCK;
				int jstart = (nj * jb) / NJBLOCK;
				if(fo.nnb <= NBMAX){
					for(int k=0; k<fo.nnb; k++){
						int nb = fo.neib[k];
						while(nb < jstart) nb += (1<<16);
						nb += joff[id];
						nblist[tid].push_back(nb);
						// nblist.push_back(fo.neib[k]);
					}
				}else{
					// overflow = true;
					handle_overflow(nblist[tid], ipbuf[id], jpbuf[id], nj, i, jb, fo.nnb, joff[id]);
				}
			}
		}
		int *nnbp = listbase + lmax * i;
		int *nblistp = nnbp + 1;
		int nnb = nblist[tid].size();
		if(nnb > nbmax) overflow = true;
		// assert(!overflow);
		if(overflow){
			*nnbp = -1;
		}else{
			*nnbp = nnb;
			for(int k=0; k<nnb; k++){
				nblistp[k] = nblist[tid][k];
			}
		}
	}
#if 0
	if(ni > 0){
		FILE *fp = fopen("Force.gpu", "w");
		assert(fp);
		for(int i=0; i<ni; i++){
			int nnb =  listbase[i*lmax];
			fprintf(fp, "%d %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %d\n",
					i, acc[i][0], acc[i][1], acc[i][2], 
					   jrk[i][0], jrk[i][1], jrk[i][2], nnb);
		}
		fprintf(fp, "\n");
		fclose(fp);
		exit(1);
	}
#endif
	// time_grav += get_wtime();
	time_reduce += get_wtime();
}

extern "C" {
	void gpunb_devinit_(){
		GPUNB_devinit();
	}
	void gpunb_open_(int *nbmax){
		GPUNB_open(*nbmax);
	}
	void gpunb_close_(){
		GPUNB_close();
	}
	void gpunb_send_(
			int *nj,
			double mj[],
			double xj[][3],
			double vj[][3]){
		GPUNB_send(*nj, mj, xj, vj);
	}
	void gpunb_regf_(
			int *ni,
			double h2[],
			double xi[][3],
			double vi[][3],
			double acc[][3],
			double jrk[][3],
			double pot[],
			int *lmax,
			int *nbmax,
			int *list){ // list[][lmax]
		GPUNB_regf(*ni, h2, xi, vi, acc, jrk, pot, *lmax, *nbmax, list);
	}
}
