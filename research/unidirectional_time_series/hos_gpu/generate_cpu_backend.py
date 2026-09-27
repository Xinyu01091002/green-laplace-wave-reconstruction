"""Mechanical CUDA -> OpenMP execution mapping; no physics formulas changed."""
import pathlib
import re
import sys
source=pathlib.Path(__file__).resolve().parent
target=pathlib.Path(sys.argv[1]);target.mkdir(parents=True,exist_ok=True)
for name in ['device_fft.cuh','second_order_rhs.cuh','hos_rhs.cuh','rk4.cuh']:
    text=(source/name).read_text()
    text=text.replace('#include <cuda_runtime.h>','#include "cpu_cuda.hpp"').replace('#include <cufft.h>','')
    text=text.replace('__global__ ','').replace('<<<256,256>>>','')
    # Only parallelise the outer per-cell kernel loops, not HOS degree loops.
    text=re.sub(r'(^[ \t]*)(for\s*\(size_t i\s*=)',r'\1#pragma omp parallel for schedule(static)\n\1\2',text,flags=re.M)
    if '<<<' in text or '__global__' in text:raise RuntimeError('Untranslated kernel launch')
    (target/name).write_text(text)
