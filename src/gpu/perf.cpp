#include "grid.cuh"

namespace ParaFROST {
    void addKernelPerfCount(float runtime, grid_t blocks, int op_type) {
		float previous_runtime = 0.0;
		GlobalPerfCount *r = (GlobalPerfCount *)malloc(sizeof(GlobalPerfCount));
		r->runtime = runtime;
		r->blocks = blocks;
		r->previous_runtime = previous_runtime;
		r->next = nullptr;

		if (perf_count_head == nullptr) {
			perf_count_head = r;
		}
		else {
			GlobalPerfCount *iterator = perf_count_head;
			while (iterator->next != nullptr) {
				iterator = iterator->next;
				if (iterator->op_type == op_type)
					previous_runtime = iterator->runtime;
			}
			iterator->next = r;
			r->previous_runtime = previous_runtime;
		}
	}

	void cleanupPerfCount() {
		GlobalPerfCount* temp = perf_count_head;
		while (temp != nullptr) {
			GlobalPerfCount* toDelete = temp;
			temp = temp->next;
			free(toDelete);
		}
		perf_count_head = nullptr;
	}

	grid_t adjustBlocksBasedOnHistory(grid_t currentBlocks, int op_type) {

		GlobalPerfCount* temp = perf_count_head;
		float totalRuntime = 0.0f;
		int count = 0;

		while (temp != nullptr) {
			if (temp->op_type == op_type) {
				totalRuntime += temp->runtime;
				count++;
			}
			temp = temp->next;
		}

		float avgRuntime = count > 0 ? totalRuntime / count : 0.0f;


		if (count > 1) {
			GlobalPerfCount *lastRun = nullptr;
			temp = perf_count_head;
			while (temp->next != nullptr) {
				if (temp->op_type == op_type)
					lastRun = temp;
				temp = temp->next; 
			}

			if (lastRun->runtime > avgRuntime * 1.1) {
				return (grid_t)std::max(1, (int)(currentBlocks * 0.9)); 
			}

			if (lastRun->runtime < avgRuntime * 0.9) {
				return currentBlocks + 1;
			}
		}

		return currentBlocks;
	}
}