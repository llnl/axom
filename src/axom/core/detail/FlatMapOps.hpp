// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef Axom_Core_Detail_FlatMapOps_Hpp
#define Axom_Core_Detail_FlatMapOps_Hpp

#include "axom/core/detail/FlatTable.hpp"

namespace axom
{
namespace detail
{
namespace flat_map
{

inline void setSentinel(axom::ArrayView<GroupBucket> metadata)
{
#if defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_CUDA)
  // Note: HIP can access device memory from the host and does not need special
  // handling - we just defer to the host path in all cases.
  MemorySpace space = getAllocatorSpace(metadata.getAllocatorID());
  if(space == MemorySpace::Device)
  {
    GroupBucket last_bucket;
    axom::copy(&last_bucket, metadata.data() + metadata.size() - 1, sizeof(GroupBucket));
    last_bucket.setSentinel();
    axom::copy(metadata.data() + metadata.size() - 1, &last_bucket, sizeof(GroupBucket));
    return;
  }
#endif
  metadata[metadata.size() - 1].setSentinel();
}

template <typename KVPair, typename LookupPolicy, typename StoragePair = TypeErasedStorage<KVPair>>
inline void destroyBuckets(axom::ArrayView<GroupBucket> metadata, axom::ArrayView<StoragePair> buckets)
{
  if(std::is_trivially_destructible<KVPair>::value)
  {
    // Nothing to do.
    return;
  }

#if defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_CUDA)
  // Note: HIP can access device memory from the host and does not need special
  // handling - we just defer to the host path in all cases.
  MemorySpace space = getAllocatorSpace(metadata.getAllocatorID());
  // CUDA-only: buckets located in device-only memory and non-trivially
  // destructible. We'll need to "relocate" the objects to the host to
  // destroy.
  //
  // This presumes the objects are "trivially relocatable." (See a similar
  // requirement in axom::Array for device allocations on CUDA)
  axom::Array<StoragePair> buckets_host;
  axom::Array<GroupBucket> metadata_host;
  if(space == MemorySpace::Device)
  {
    int host_allocator_id = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
    metadata_host = axom::Array<GroupBucket>(metadata, host_allocator_id);
    buckets_host = axom::Array<StoragePair>(buckets, host_allocator_id);
    metadata = metadata_host.view();
    buckets = buckets_host.view();
  }
#endif
  IndexType index = LookupPolicy {}.nextValidIndex(metadata, LookupPolicy::NO_MATCH);
  while(index < buckets.size())
  {
    buckets[index].get().~KVPair();
    index = LookupPolicy {}.nextValidIndex(metadata, index);
  }
}

template <typename KVPair, typename LookupPolicy, typename StoragePair = TypeErasedStorage<KVPair>>
inline void copyBuckets(axom::ArrayView<const GroupBucket> metadata,
                        axom::ArrayView<const StoragePair> from_buckets,
                        axom::ArrayView<StoragePair> to_buckets)
{
  if(std::is_trivially_copyable<KVPair>::value)
  {
    axom::copy(to_buckets.data(), from_buckets.data(), sizeof(StoragePair) * from_buckets.size());
    return;
  }

  axom::ArrayView<StoragePair> to_buckets_stage = to_buckets;
#if defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_CUDA)
  // Note: HIP can access device memory from the host and does not need special
  // handling - we just defer to the host path in all cases.

  // Non-trivially copyable:
  // "Relocate" to the host to call host-based copy constructor.
  MemorySpace meta_space = getAllocatorSpace(metadata.getAllocatorID());
  axom::Array<GroupBucket> metadata_host;
  if(meta_space == MemorySpace::Device)
  {
    int host_allocator_id = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
    metadata_host = axom::Array<GroupBucket>(metadata, host_allocator_id);
    metadata = metadata_host.view();
  }
  MemorySpace from_space = getAllocatorSpace(from_buckets.getAllocatorID());
  axom::Array<StoragePair> from_buckets_host;
  if(from_space == MemorySpace::Device)
  {
    int host_allocator_id = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
    from_buckets_host = axom::Array<StoragePair>(from_buckets, host_allocator_id);
    from_buckets = from_buckets_host.view();
  }

  MemorySpace to_space = getAllocatorSpace(to_buckets.getAllocatorID());
  axom::Array<StoragePair> to_buckets_host;
  if(to_space == MemorySpace::Device)
  {
    int host_allocator_id = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
    to_buckets_host = axom::Array<StoragePair>(to_buckets, host_allocator_id);
    to_buckets_stage = to_buckets_host.view();
  }
#endif
  // Copy all elements.
  IndexType index = LookupPolicy {}.nextValidIndex(metadata, LookupPolicy::NO_MATCH);
  while(index < from_buckets.size())
  {
    new(&to_buckets_stage[index].data) KVPair(from_buckets[index].get());
    index = LookupPolicy {}.nextValidIndex(metadata, index);
  }

#if defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_CUDA)
  if(to_space == MemorySpace::Device)
  {
    // We created new elements in a host staging buffer, relocate them to
    // our desired buffer.
    axom::copy(to_buckets.data(), to_buckets_stage.data(), sizeof(StoragePair) * to_buckets.size());
  }
#endif
}

}  // namespace flat_map
}  // namespace detail
}  // namespace axom

#endif
