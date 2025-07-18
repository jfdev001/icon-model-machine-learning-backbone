// ICON
//
// ---------------------------------------------------------------
// Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
// Contact information: icon-model.org
//
// See AUTHORS.TXT for a list of authors
// See LICENSES/ for license information
// SPDX-License-Identifier: BSD-3-Clause
// ---------------------------------------------------------------

#include <cuda.h>

#include <cub/device/device_select.cuh>
#include <cub/iterator/counting_input_iterator.cuh>
#include <memory>
#include <unordered_map>

#include "index_list.h"

// Use stream-ordered allocs if we're capturing a graph
namespace {
bool isStreamCapturing(gpuStream_t stream) {
  cudaStreamCaptureStatus captureStatus;
  cudaStreamIsCapturing(stream, &captureStatus);
  return captureStatus != cudaStreamCaptureStatusNone;
}

class Storage {
 public:
  virtual void requestSize(size_t requestedSize) = 0;
  int* getNvalidPtr() { return reinterpret_cast<int*>(data); }
  char* getScratchPtr() { return data + alignment; }
  virtual ~Storage() = default;

 protected:
  char* data                 = nullptr;
  static const int alignment = 512;
};

class AsyncStorage : public Storage {
 public:
  AsyncStorage(gpuStream_t stream) : stream(stream) {}

  void requestSize(size_t requestedSize) override final {
    if (data != nullptr) {
      cudaFreeAsync(data, stream);
    }
    cudaMallocAsync(&data, alignment + requestedSize, stream);
  }
  ~AsyncStorage() override { cudaFreeAsync(data, stream); }

 private:
  gpuStream_t stream;
};

class SyncStorage : public Storage {
 public:
  void requestSize(size_t requestedSize) override final {
    if (curSize < requestedSize + alignment) {
      cudaFree(data);
      cudaMalloc(&data, requestedSize + alignment);
      curSize = requestedSize + alignment;
    }
  }
  ~SyncStorage() override { cudaFree(data); }

 private:
  size_t curSize = 0;
};

std::unordered_map<gpuStream_t, std::shared_ptr<SyncStorage>> syncStorageMap;

// Use async storage in case we're capturing a graph
// otherwise the sync storage per-stream
std::shared_ptr<Storage> getStorage(gpuStream_t stream) {
  if (isStreamCapturing(stream)) {
    return std::make_shared<AsyncStorage>(stream);
  } else {
    if (syncStorageMap.find(stream) == syncStorageMap.end()) {
      syncStorageMap[stream] = std::make_shared<SyncStorage>();
    }
    return syncStorageMap[stream];
  }
}
}  // namespace

template <typename T>
static void c_generate_index_list_gpu_generic_device(const T* dev_conditions, const int startid, const int endid,
                                                     int* dev_indices, int* dev_nvalid, Storage* storage,
                                                     gpuStream_t stream) {
  const int n = endid - startid + 1;

  // Argument is the offset of the first element
  cub::CountingInputIterator<int> iterator(startid);

  // Determine temporary device storage requirements
  size_t storageRequirement;
  cub::DeviceSelect::Flagged(nullptr, storageRequirement, iterator, dev_conditions + startid - 1, dev_indices,
                             dev_nvalid, n, stream);

  // Allocate temporary storage
  storage->requestSize(storageRequirement);
  if (dev_nvalid == nullptr) {
    dev_nvalid = storage->getNvalidPtr();
  }

  cub::DeviceSelect::Flagged(storage->getScratchPtr(), storageRequirement, iterator, dev_conditions + startid - 1,
                             dev_indices, dev_nvalid, n, stream);
}

template <typename T>
static void c_generate_index_list_gpu_batched_generic(const int batch_size, const T* dev_conditions,
                                                      const int cond_stride, const int startid, const int endid,
                                                      int* dev_indices, const int idx_stride, int* dev_nvalid,
                                                      gpuStream_t stream) {
  auto storage = getStorage(stream);

  for (int i = 0; i < batch_size; i++)
    c_generate_index_list_gpu_generic_device(dev_conditions + cond_stride * i, startid, endid,
                                             dev_indices + idx_stride * i, dev_nvalid + i, storage.get(), stream);
}

template <typename T>
static void c_generate_index_list_gpu_generic(const T* dev_conditions, const int startid, const int endid,
                                              int* dev_indices, int* ptr_nvalid, bool copy_to_host,
                                              gpuStream_t stream) {
  auto storage = getStorage(stream);

  c_generate_index_list_gpu_generic_device(dev_conditions, startid, endid, dev_indices,
                                           copy_to_host ? storage->getNvalidPtr() : ptr_nvalid, storage.get(), stream);

  if (copy_to_host) {
    cudaMemcpyAsync(ptr_nvalid, storage->getNvalidPtr(), sizeof(int), cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);
  }
}

///
/// Exposed functions
///
/// Non-batched first
///
void c_generate_index_list_gpu_single(const void* dev_conditions, const int startid, const int endid, int* dev_indices,
                                      int* nvalid, int data_size, bool copy_to_host, gpuStream_t stream) {
  switch (data_size) {
    case 1:
      c_generate_index_list_gpu_generic(static_cast<const char*>(dev_conditions), startid, endid, dev_indices, nvalid,
                                        copy_to_host, stream);
      break;
    case 4:
      c_generate_index_list_gpu_generic(static_cast<const int*>(dev_conditions), startid, endid, dev_indices, nvalid,
                                        copy_to_host, stream);
      break;
  }
}

///
/// And now batched
///
void c_generate_index_list_gpu_batched(const int batch_size, const void* dev_conditions, const int cond_stride,
                                       const int startid, const int endid, int* dev_indices, const int idx_stride,
                                       int* dev_nvalid, int data_size, gpuStream_t stream) {
  switch (data_size) {
    case 1:
      c_generate_index_list_gpu_batched_generic(batch_size, static_cast<const char*>(dev_conditions), cond_stride,
                                                startid, endid, dev_indices, idx_stride, dev_nvalid, stream);
      break;

    case 4:
      c_generate_index_list_gpu_batched_generic(batch_size, static_cast<const int*>(dev_conditions), cond_stride,
                                                startid, endid, dev_indices, idx_stride, dev_nvalid, stream);
      break;
  }
}
