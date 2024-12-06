/***********************************************************************[grid.cuh]
Copyright(c) 2020, Muhammad Osama - Anton Wijs,
Copyright(c) 2022-present, Muhammad Osama.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**********************************************************************************/

#ifndef __GRID_INFO_
#define __GRID_INFO_

#include "definitions.cuh"
#include "datatypes.hpp"

namespace ParaFROST {

	// x
	typedef uint32 grid_t;
	struct GlobalPerfCount{
		float previous_runtime;
		float runtime;
		grid_t blocks;
		int op_type;
		GlobalPerfCount *next;
	};

	extern GlobalPerfCount *perf_count_head = nullptr;
	void addKernelPerfCount(float runtime, grid_t blocks, int op_type);
	void cleanupPerfCount();
	grid_t adjustBlocksBasedOnHistory(grid_t currentBlocks, int op_type);
	#define OP_TYPE_ERE 0
	#define OP_TYPE_ELIM 1
	#define global_bx		(grid_t)(blockDim.x * blockIdx.x)
	#define global_bx_off	(grid_t)((blockDim.x << 1) * blockIdx.x)
	#define global_tx		(grid_t)(global_bx + threadIdx.x)
	#define global_tx_off	(grid_t)(global_bx_off + threadIdx.x)
	#define stride_x        (grid_t)(blockDim.x * gridDim.x)
	#define stride_x_off	(grid_t)((blockDim.x << 1) * gridDim.x)
	// y
	#define global_by		(grid_t)(blockDim.y * blockIdx.y)
	#define global_ty		(grid_t)(global_by + threadIdx.y)
	#define stride_y		(grid_t)(blockDim.y * gridDim.y)

 #define ADAPTIVE_INIT(NTHREADS, MINTHREADS, DEVICEPROPS, FUNCATTR)          \
        {                                                                      \
			int deviceCount;													\
			cudaGetDeviceCount(&deviceCount);									\
			cudaDeviceProp deviceProp;											\
			if (deviceCount) {													\
				cudaGetDeviceProperties(&deviceProp, 0);						\
			}																	\
            NTHREADS = deviceProp.maxThreadsPerBlock;                         \
            MINTHREADS = deviceProp.warpSize;                                 \
                                                                               \
            if (FUNCATTR.sharedSizeBytes > deviceProp.sharedMemPerBlock / 2)  \
                NTHREADS /= 2;                                                 \
                                                                               \
            if (FUNCATTR.numRegs > deviceProp.regsPerBlock / NTHREADS)        \
                NTHREADS /= 2;                                                 \
                                                                               \
            NTHREADS = max(MINTHREADS,                                         \
                           (NTHREADS / DEVICEPROPS.warpSize) * DEVICEPROPS.warpSize); \
        }

    // Macros for profiling-based tuning
    #define PROFILE_BASED_TUNING(DATALEN, DEVICEPROPS, FUNCATTR, GRID, BLOCK) \
        {                                                                     \
            int minGridSize;                                                  \
            cudaOccupancyMaxPotentialBlockSize(&minGridSize, &BLOCK,          \
                                               kernelFunc, 0,                 \
                                               DEVICEPROPS.maxThreadsPerBlock); \
                                                                               \
            BLOCK = max(BLOCK, DEVICEPROPS.warpSize);                         \
            BLOCK = min(BLOCK, DEVICEPROPS.maxThreadsPerBlock);               \
            GRID = (DATALEN + BLOCK - 1) / BLOCK;                             \
                                                                               \
            GRID = min(GRID, DEVICEPROPS.multiProcessorCount * DEVICEPROPS.maxBlocksPerMultiProcessor); \
        }

    // Hybrid tuning for iterative workloads
    #define HYBRID_TUNING(DATALEN, GRID, BLOCK, ITERATION)                     \
        {                                                                     \
            if (ITERATION > 0 && DATALEN < GRID * BLOCK) {                    \
                GRID = (DATALEN + BLOCK - 1) / BLOCK;                         \
            }                                                                 \
            grid_t minBlocks = GRID * 0.5;                                    \
            GRID = max(GRID, minBlocks);                                      \
        }

    // Dynamic feedback integration
    #define DYNAMIC_FEEDBACK(KERNEL, BLOCK, GRID)                              \
        {                                                                     \
            cudaDeviceSynchronize();                                          \
            float achievedOccupancy;                                         \
            cudaOccupancy(&achievedOccupancy, KERNEL, BLOCK, GRID);           \
                                                                               \
            if (achievedOccupancy < 0.5) BLOCK /= 2;                          \
        }

    // Bandwidth-aware tuning
    #define BANDWIDTH_AWARE_TUNING(BLOCK, DEVICEPROPS, DATASIZE)               \
        {                                                                     \
            int threadsPerSM = DEVICEPROPS.maxThreadsPerMultiProcessor / BLOCK; \
            int memoryTransactions = threadsPerSM * DATASIZE;                 \
                                                                               \
            if (memoryTransactions > DEVICEPROPS.memoryBusWidth / 2) {        \
                BLOCK /= 2;                                                   \
            }                                                                 \
        }

	// macros for blocks calculation
	#define ROUNDUPBLOCKS(DATALEN, NTHREADS)							     \
			(((DATALEN) + (NTHREADS) - 1) / (NTHREADS))

	#define OPTIMIZEBLOCKS(DATALEN, NTHREADS)                                \
			assert(DATALEN);                                                 \
			assert(NTHREADS);                                                \
			assert(maxGPUThreads);                                           \
			const grid_t REALBLOCKS = ROUNDUPBLOCKS(DATALEN, NTHREADS);    \
			const grid_t MAXBLOCKS = maxGPUThreads / NTHREADS;             \
			const grid_t nBlocks = MIN(REALBLOCKS, MAXBLOCKS);             \

	#define OPTIMIZEBLOCKSELIM(NVARS, MAXTHREADS, MINOPTS)                   \
			assert(NVARS);                                                   \
			assert(MAXTHREADS);                                              \
			assert(maxGPUThreads);                                           \
			const grid_t MINTHREADS = MINOPTS ## _min_threads;		     \
			grid_t nThreads = MAXTHREADS;								     \
			grid_t realBlocks = ROUNDUPBLOCKS(NVARS, nThreads);            \
			const grid_t MAXBLOCKS = maxGPUThreads / MAXTHREADS;           \
			const grid_t MINBLOCKS =										 \
				  grid_t(MAXBLOCKS * (MINOPTS ## _min_blocks));		     \
			while (nThreads > MINTHREADS && realBlocks <= MINBLOCKS) {	     \
				nThreads >>= 1;											     \
				realBlocks = ROUNDUPBLOCKS(NVARS, nThreads);			     \
			}															     \
			const grid_t nBlocks = MIN(realBlocks, MAXBLOCKS);             \

	#define OPTIMIZEBLOCKSERE(NVARS, BLOCK2D, MINOPTS)                       \
			assert(NVARS);                                                   \
			assert(BLOCK2D.x);                                               \
			assert(BLOCK2D.y);                                               \
			assert(maxGPUThreads);                                           \
			const grid_t MINTHREADS = MINOPTS ## _min_threads;		     \
			grid_t realBlocks = ROUNDUPBLOCKS(NVARS, BLOCK2D.y);		     \
			const grid_t MAXBLOCKS = maxGPUThreads /					     \
										(BLOCK2D.x * BLOCK2D.y);             \
			const grid_t MINBLOCKS =										 \
				  grid_t(MAXBLOCKS * (MINOPTS ## _min_blocks));		     \
			while (BLOCK2D.y > MINTHREADS && realBlocks <= MINBLOCKS) {	     \
				BLOCK2D.y >>= 1;											 \
				realBlocks = ROUNDUPBLOCKS(NVARS, BLOCK2D.y);			     \
			}																 \
			const grid_t nBlocks = MIN(realBlocks, MAXBLOCKS);             \

	#define OPTIMIZEBLOCKS2(DATALEN, NTHREADS)                               \
			assert(DATALEN);                                                 \
			assert(NTHREADS);                                                \
			assert(maxGPUThreads);                                           \
			const grid_t REALTHREADS = (NTHREADS) << 1;					 \
			const grid_t REALBLOCKS = ROUNDUPBLOCKS(DATALEN, REALTHREADS); \
			const grid_t MAXBLOCKS = maxGPUThreads / REALTHREADS;          \
			const grid_t nBlocks = MIN(REALBLOCKS, MAXBLOCKS);             \

	// macros for shared memory calculation
    #define OPTIMIZESHARED(NTHREADS, MINCAP)                 \
            assert(MINCAP);                                  \
            assert(NTHREADS);                                \
            assert(maxGPUSharedMem);                         \
            const size_t smemSize = (NTHREADS) * (MINCAP);   \
            assert(maxGPUSharedMem >= smemSize);             \

}

#endif