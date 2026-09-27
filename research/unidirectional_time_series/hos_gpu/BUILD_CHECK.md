# Build checkpoint, 27 September 2026 (Europe/London)

Status: CUDA compilation and linking passed; GPU numerical tests NOT RUN.

Source snapshot `/home/lxy/hos-gpu-buildcheck-20260927-v2` on the existing
60093 SSH endpoint was compiled with GNU 13.3.0 and CUDA 12.9.86, targeting
sm_86. Both FP32 and FP64 template instantiations compiled and linked.
CTest marked both device tests **Skipped**, because no usable CUDA device
was visible. The CTest headline is not evidence of numerical test passes.

The earlier `v1` snapshot is retained: its CUDA template/conditional parser
errors were fixed locally before uploading a fresh `v2` snapshot.

SHA-256, verified equal locally and remotely before the successful build:

| File | SHA-256 |
|---|---|
| CMakeLists.txt | aab476a791cb67758f81327466e1d134f8295c242dc8afe99e4bef1d85cdebb8 |
| device_fft.cuh | b93301b79d9c97cade2e866abcea276aa7c2c534c6e55aade815ce07fb8a9d56 |
| fft_check.cu | 4a31ac179d96b196a1023ada576e9cb3e582af105add7edd9361cc0a142f6af8 |
| README.md | 8b807539622c9158f43c56662842efc9030205dd09309cda3179c705949f9f88 |

Separate hardware observation on `lxy@60.188.112.99:22`: hostname jfm-99;
PCI graphics device VMware SVGA II [15ad:0405]; DRM card0/renderD128 present;
no NVIDIA PCI entry, nvidia-smi, /dev/nvidia devices or NVIDIA GPU proc entries
were found by the read-only checks. This establishes no visible CUDA GPU in
that guest, not absence of physical GPUs in the underlying virtualization host.

No complete HOS RHS, time integration, production replacement, speedup or
whole-solver VRAM fit has been validated at this checkpoint.
