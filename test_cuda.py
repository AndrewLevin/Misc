from cuda import cuda, nvrtc
import numpy as np

def ASSERT_DRV(err):
    if isinstance(err,cuda.CUresult):
        if err != cuda.CUresult.CUDA_SUCCESS:
            raise RuntimeError("Cuda Error: {}".format(err))
    elif isinstance(err,nvrtc.nvrtcResult):
        if err != nvrtc.nvrtcResult.NVRTC_SUCCESS:
            raise RuntimeError("Nvrtc Error: {}".format(err))
    else:
        raise RuntimeError("Unknown error type: {}".format(err))

saxpy = """\
extern "C" __global__
void saxpy(float a, float *x, float *y, float *out, size_t n)
{
 size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
 if (tid < n) {
   out[tid] = a * x[tid] + y[tid];
 }
}
"""

err, prog = nvrtc.nvrtcCreateProgram(str.encode(saxpy), b"saxpy.cu", 0, [], [])

ASSERT_DRV(err)

opts = [b"--fmad=false", b"--gpu-architecture=compute_75"]

err, = nvrtc.nvrtcCompileProgram(prog, 2, opts)

ASSERT_DRV(err)

err, ptxSize = nvrtc.nvrtcGetPTXSize(prog)

ASSERT_DRV(err)

ptx = b" " * ptxSize

err, = nvrtc.nvrtcGetPTX(prog, ptx)

ASSERT_DRV(err)

err, = cuda.cuInit(0)

ASSERT_DRV(err)

err, cuDevice = cuda.cuDeviceGet(0)

ASSERT_DRV(err)

err, context = cuda.cuCtxCreate(0, cuDevice)

ASSERT_DRV(err)

ptx = np.char.array(ptx)

err, module = cuda.cuModuleLoadData(ptx.ctypes.data)

ASSERT_DRV(err)

err, kernel = cuda.cuModuleGetFunction(module, b"saxpy")

ASSERT_DRV(err)

NUM_THREADS = 512

NUM_BLOCKS = 32768

a = np.array([2.0], dtype=np.float32)

n = np.array(NUM_THREADS * NUM_BLOCKS, dtype=np.uint32)

bufferSize = n * a.itemsize

hX = np.random.rand(n).astype(dtype=np.float32)

hY = np.random.rand(n).astype(dtype=np.float32)

hOut = np.zeros(n).astype(dtype=np.float32)

err, dXclass = cuda.cuMemAlloc(bufferSize)

ASSERT_DRV(err)

err, dYclass = cuda.cuMemAlloc(bufferSize)

ASSERT_DRV(err)

err, dOutclass = cuda.cuMemAlloc(bufferSize)

ASSERT_DRV(err)

err, stream = cuda.cuStreamCreate(0)

ASSERT_DRV(err)

err, = cuda.cuMemcpyHtoDAsync(
    dXclass, hX.ctypes.data, bufferSize, stream
)

ASSERT_DRV(err)

err, = cuda.cuMemcpyHtoDAsync(
    dYclass, hY.ctypes.data, bufferSize, stream
)

ASSERT_DRV(err)

dX = np.array([int(dXclass)], dtype=np.uint64)

dY = np.array([int(dYclass)], dtype=np.uint64)

dOut = np.array([int(dOutclass)], dtype=np.uint64)

args = [a, dX, dY, dOut, n]

args = np.array([arg.ctypes.data for arg in args], dtype=np.uint64)

err, = cuda.cuLaunchKernel(
    kernel,
    NUM_BLOCKS, # grid x dim
    1, # grid y dim
    1, # grid z dim
    NUM_THREADS, # block x dim
    1, #block y dim
    1, #block z dim
    0, # dynamic shared memory
    stream, # stream
    args.ctypes.data, # kernel arguments
    0, # extra (ignore)
    )

ASSERT_DRV(err)

err, = cuda.cuMemcpyDtoHAsync(
    hOut.ctypes.data, dOutclass, bufferSize, stream
    )

ASSERT_DRV(err)

err, = cuda.cuStreamSynchronize(stream)

ASSERT_DRV(err)

hZ = a * hX + hY

if not np.allclose(hOut, hZ):
    raise ValueError("Error outside tolerance for host-device vectors")

err, = cuda.cuStreamDestroy(stream)

ASSERT_DRV(err)

err, = cuda.cuMemFree(dXclass)

ASSERT_DRV(err)

err, = cuda.cuMemFree(dYclass)

ASSERT_DRV(err)

err, = cuda.cuMemFree(dOutclass)

ASSERT_DRV(err)

err, = cuda.cuModuleUnload(module)

ASSERT_DRV(err)

err, = cuda.cuCtxDestroy(context)

ASSERT_DRV(err)
