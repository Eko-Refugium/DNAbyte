# gpu_dna_levenshtein.py

import numpy as np
import pycuda.autoinit
import pycuda.driver as cuda
from pycuda.compiler import SourceModule



def init_gpu_simulation(length):
    """One-time GPU init"""
    global simulate_kernel, mod
    mod = SourceModule(SIMULATION_KERNEL)
    simulate_kernel = mod.get_function("simulate_errors")
    self._gpu_simulate_ready = True

def gpu_simulate_errors(data, error_rate, batch_size=1024):
    """GPU-accelerated error simulation - COMPILES ON RTX 30xx"""
    if not data:
        return [], {'error_counter': 0}
    
    num_sequences = len(data)
    seq_len = len(data[0])
    
    # Compile and run
    mod = SourceModule(SIMULATION_KERNEL)
    simulate_kernel = mod.get_function("simulate_errors")
    
    # Memory allocation
    seq_bytes = strings_to_bytes_fixed(data, seq_len)
    d_input = cuda.mem_alloc(seq_bytes.nbytes)
    d_output = cuda.mem_alloc(seq_bytes.nbytes)
    d_errors = cuda.mem_alloc(num_sequences * 4)
    
    cuda.memcpy_htod(d_input, seq_bytes.ravel())
    
    # Launch kernel
    block_size = (256, 1, 1)
    grid_size = ((num_sequences + 255) // 256, 1)
    
    simulate_kernel(
        d_input, np.int32(num_sequences), np.int32(seq_len),
        np.float32(error_rate), d_output, d_errors,
        block=block_size, grid=grid_size
    )
    
    # Results
    output_bytes = np.empty_like(seq_bytes)
    cuda.memcpy_dtoh(output_bytes, d_output)
    error_counts = np.empty(num_sequences, dtype=np.int32)
    cuda.memcpy_dtoh(error_counts, d_errors)
    
    # Convert to DNA strings
    result_sequences = [''.join(chr(b) for b in output_bytes[i] if 65 <= b <= 84) 
                       for i in range(num_sequences)]
    
    total_errors = int(error_counts.sum())
    
    # Cleanup
    d_input.free()
    d_output.free()
    d_errors.free()
    
    return result_sequences, {'error_counter': total_errors}
KERNEL_CODE = """
__global__ void levenshtein_dna(
    const unsigned char* queries, int fixed_len, 
    const unsigned char* candidates, int num_candidates, 
    int* distances
)
{
    int query_idx = blockIdx.y * blockDim.y + threadIdx.y;
    int cand_idx  = blockIdx.x * blockDim.x + threadIdx.x;

    if (query_idx >= 64 || cand_idx >= num_candidates) return;

    int q_offset = query_idx * fixed_len;
    int c_offset = cand_idx * fixed_len;
    const unsigned char* query = queries + q_offset;
    const unsigned char* cand  = candidates + c_offset;

    int dp[2][257];
    int* prev = dp[0];
    int* curr = dp[1];

    prev[0] = 0;
    for (int j = 1; j <= fixed_len; j++) {
        prev[j] = j;
    }

    for (int i = 1; i <= fixed_len; i++) {
        curr[0] = i;
        for (int j = 1; j <= fixed_len; j++) {
            int cost = (query[i-1] == cand[j-1]) ? 0 : 1;
            int del = prev[j] + 1;
            int ins = curr[j-1] + 1;
            int sub = prev[j-1] + cost;
            curr[j] = min(del, min(ins, sub));
        }
        int* temp = prev; prev = curr; curr = temp;
    }

    int result_idx = query_idx * num_candidates + cand_idx;
    distances[result_idx] = prev[fixed_len];
}
"""

# Global cache
candidates_gpu = None
num_candidates_global = 0
kernel = None
fixed_len = 0


class MultiGPUMatcher:
    def __init__(self, candidates, length):
        self.length = length
        self.num_gpus = cuda.Device.count()
        print(f"Found {self.num_gpus} GPUs")
        
        # Create one matcher per GPU
        self.matchers = []
        for gpu_id in range(self.num_gpus):
            matcher = GPUCache(candidates, length, gpu_id)
            self.matchers.append(matcher)
    
    def find_closest(self, queries):
        """Split queries across all GPUs, return list of (idx, dist)"""
        chunk_size = len(queries) // self.num_gpus
        chunks = [queries[i:i+chunk_size] for i in range(0, len(queries), chunk_size)]
        if len(chunks) < self.num_gpus:
            chunks.extend([[] for _ in range(self.num_gpus - len(chunks))])
        
        # Run each GPU in parallel
        with Pool(self.num_gpus) as pool:
            results = pool.starmap(self._process_gpu_chunk, 
                                 enumerate(chunks))
        
        # Flatten results
        all_results = []
        for chunk_results in results:
            all_results.extend(chunk_results)
        return all_results
    
    def _process_gpu_chunk(self, gpu_id, chunk):
        """Process one chunk on one GPU"""
        if not chunk:
            return []
        return self.matchers[gpu_id].find_closest(chunk)

def gpu_levenshtein_dna_batch_multi_gpu(queries, batch_size=64):
    """Distribute across ALL available GPUs"""
    global candidates_gpu_list, num_candidates_per_gpu, kernels, fixed_len
    
    num_gpus = len(candidates_gpu_list)
    all_results = [[] for _ in range(num_gpus)]
    
    # Each GPU processes its candidates in parallel
    for gpu_id in range(num_gpus):
        candidates_gpu = candidates_gpu_list[gpu_id]
        num_cands = num_candidates_per_gpu[gpu_id]
        kernel = kernels[gpu_id]
        
        # Same batch processing per GPU
        for i in range(0, len(queries), batch_size):
            batch = queries[i:i+batch_size]
            if not batch: break
            
            batch_bytes = strings_to_bytes_fixed(batch, fixed_len)
            d_queries = cuda.mem_alloc(batch_bytes.nbytes)
            d_distances = cuda.mem_alloc(len(batch) * num_cands * 4)
            
            cuda.memcpy_htod(d_queries, batch_bytes.ravel())
            
            # Optimized block/grid (32x16)
            block_size = (32, 16, 1)
            grid_size = (
                (num_cands + 31) // 32,
                (len(batch) + 15) // 16
            )
            
            # Set active GPU
            cuda.Context.pop()  # Clear previous context
            cuda.init()         # Re-init
            cuda.Device(gpu_id).make_context()
            
            kernel(
                d_queries, np.int32(fixed_len),
                candidates_gpu, np.int32(num_cands),
                d_distances,
                block=block_size, grid=grid_size
            )
            
            batch_results = np.empty((len(batch), num_cands), dtype=np.int32)
            cuda.memcpy_dtoh(batch_results, d_distances)
            
            # Offset indices by this GPU's candidates
            for j in range(len(batch)):
                local_best_idx = np.argmin(batch_results[j])
                global_best_idx = num_candidates_per_gpu[:gpu_id].sum() + local_best_idx
                best_dist = batch_results[j][local_best_idx]
                all_results[gpu_id].append((global_best_idx, best_dist))
            
            d_queries.free()
            d_distances.free()
    
    # Merge: take global minimum across all GPUs
    final_results = []
    for query_i in range(len(queries)):
        best_across_gpus = min(all_results[gpu_id][query_i % len(all_results[gpu_id])] 
                              for gpu_id in range(num_gpus))
        final_results.append(best_across_gpus)
    
    return final_results

def init_gpu_cache_multi_gpu(candidates, length=100):
    """Split candidates across ALL GPUs"""
    global candidates_gpu_list, num_candidates_per_gpu, kernels, fixed_len
    
    fixed_len = length
    num_gpus = cuda.Device.count()
    print(f"Found {num_gpus} GPUs - splitting {len(candidates)} candidates")
    
    candidates_gpu_list = []
    num_candidates_per_gpu = []
    kernels = []
    
    # Split candidates roughly equally
    candidates_per_gpu = len(candidates) // num_gpus
    remainder = len(candidates) % num_gpus
    
    start_idx = 0
    for gpu_id in range(num_gpus):
        end_idx = start_idx + candidates_per_gpu + (1 if gpu_id < remainder else 0)
        gpu_candidates = candidates[start_idx:end_idx]
        num_cands = len(gpu_candidates)
        
        print(f"  GPU {gpu_id}: {num_cands} candidates")
        
        # Upload to this GPU
        cuda.Context.pop() if cuda.Context.get_current_context() else None
        cuda.Device(gpu_id).make_context()
        
        gpu_bytes = strings_to_bytes_fixed(gpu_candidates, length)
        gpu_mem = cuda.mem_alloc(gpu_bytes.nbytes)
        cuda.memcpy_htod(gpu_mem, gpu_bytes.ravel())
        
        candidates_gpu_list.append(gpu_mem)
        num_candidates_per_gpu.append(num_cands)
        
        # Compile kernel for this GPU
        mod = SourceModule(KERNEL_CODE)
        kernels.append(mod.get_function("levenshtein_dna"))
        
        start_idx = end_idx
    
    print("Multi-GPU ready!")


def strings_to_bytes_fixed(strings, length):
    arr = np.zeros((len(strings), length), dtype=np.uint8)
    for i, s in enumerate(strings):
        s_bytes = s.encode("ascii")
        arr[i, :len(s_bytes)] = np.frombuffer(s_bytes, dtype=np.uint8)
    return arr


def init_gpu_cache(candidates, length):
    """Upload library strings to GPU once."""
    global candidates_gpu, num_candidates_global, kernel, fixed_len

    fixed_len = length
    cand_bytes = strings_to_bytes_fixed(candidates, length)

    candidates_gpu = cuda.mem_alloc(cand_bytes.nbytes)
    cuda.memcpy_htod(candidates_gpu, cand_bytes.ravel())
    num_candidates_global = len(candidates)

    mod = SourceModule(KERNEL_CODE)
    kernel = mod.get_function("levenshtein_dna")


def gpu_levenshtein_dna_batch(queries, batch_size=64):
    """Return (best_idx, best_dist) for each query."""
    global candidates_gpu, num_candidates_global, kernel, fixed_len

    all_results = []

    for i in range(0, len(queries), batch_size):
        batch = queries[i:i + batch_size]
        if not batch:
            break

        batch_bytes = strings_to_bytes_fixed(batch, fixed_len)

        d_queries = cuda.mem_alloc(batch_bytes.nbytes)
        d_distances = cuda.mem_alloc(len(batch) * num_candidates_global * 4)

        cuda.memcpy_htod(d_queries, batch_bytes.ravel())

        threads_x = 32  # candidates (warp size = coalesced)
        threads_y = 8  # queries (double your current 8)
        block_size = (threads_x, threads_y, 1)  # 512 threads vs your 256
        grid_size = (
            (num_candidates_global + threads_x - 1) // threads_x,
            (len(batch) + threads_y - 1) // threads_y
        )
        kernel(
            d_queries, np.int32(fixed_len),
            candidates_gpu, np.int32(num_candidates_global),
            d_distances,
            block=block_size, grid=grid_size
        )

        batch_results = np.empty(
            (len(batch), num_candidates_global), dtype=np.int32
        )
        cuda.memcpy_dtoh(batch_results, d_distances)

        for j in range(len(batch)):
            best_idx = int(np.argmin(batch_results[j]))
            best_dist = int(batch_results[j][best_idx])
            all_results.append((best_idx, best_dist))

        d_queries.free()
        d_distances.free()

    return all_results
