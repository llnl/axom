// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/quest/Discretize.hpp"
#include "axom/quest/SorClipper.hpp"
#include "axom/quest/GeometryClipper.hpp"
#include "axom/spin/BVH.hpp"
#include "axom/fmt.hpp"
#include "axom/core/WhereMacro.hpp"

#include <limits>

namespace axom
{
namespace quest
{

SorClipper::SorClipper(const klee::Geometry& kGeom, const std::string& name)
  : GeometryClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("Sor") : name)
  , m_maxRadius(0.0)
  , m_minRadius(std::numeric_limits<double>::max())
{
  extractClipperInfo();

  for(auto& pt : m_sorCurve)
  {
    m_maxRadius = fmax(m_maxRadius, pt[1]);
    m_minRadius = fmin(m_minRadius, pt[1]);
  }
  SLIC_ERROR_IF(m_minRadius < 0.0,
                axom::fmt::format("SorClipper '{}' has a negative radius", m_name));

  for(const auto& pt : m_sorCurve)
  {
    m_curveBb.addPoint(pt);
  }

  FSorClipper::combineRadialSegments(m_sorCurve);

  axom::Array<axom::Array<Point2DType>> sections;
  splitIntoMonotonicSections(m_sorCurve.view(), sections);
  for(int i = 0; i < sections.size(); ++i)
  {
    axom::ArrayView<const Point2DType> section = sections[i].view();
    std::string sectionName = axom::fmt::format("{}.{}", m_name, i);
    m_fsorStrategies.push_back(
      std::make_shared<FSorClipper>(kGeom,
                                    sectionName,
                                    section,
                                    m_sorOrigin,
                                    m_sorDirection,
                                    m_levelOfRefinement));
  }
}

bool SorClipper::specializedClip(quest::ShapeeMesh& shapeeMesh,
                                 axom::ArrayView<double> ovlap)
{
  AXOM_ANNOTATE_SCOPE("SorClipper::specializedClip");
  /*
    The SOR curve has been split into SOR functions that do not double
    back on itself.  We compute the overlaps for each section and
    accumulate them with the correct sign.  (Functions going backward
    remove stuff from the functions above them, so they contribute a
    negative volume.)

    By convention, backward curves should generate negative volume,
    but for some reason, the cone discretization functionality always
    generates positive volumes.  We correct this by manually applying
    the correct sign.
  */
  const axom::IndexType cellCount = ovlap.size();
  axom::Array<double> tmpOvlap(cellCount, cellCount, ovlap.getAllocatorID());
  for(auto& fsorStrategy : m_fsorStrategies)
  {
    tmpOvlap.fill(0.0);
    GeometryClipper clipper(shapeeMesh, fsorStrategy);
    clipper.setVerbose(true);
    clipper.clip(tmpOvlap);
    auto sorCurve = fsorStrategy->getSorCurve();
    int sign = axom::utilities::sign_of(sorCurve[sorCurve.size()-1][0] - sorCurve[0][0], 0.0);
    accumulateData( ovlap, tmpOvlap.view(), double(sign),
                    shapeeMesh.getRuntimePolicy() );
  }
  return true;
}

/*
  Split curve into sections that goes monotonically up or down the
  axis of symmetry.  If x changes directions at a radial segment,
  split at the end with greater y value.  A radial segment is one with
  constant x (or z), pointing in the y (or radial) direction.

  This method assumes there are no consecutive radial segments
  (combineRadialSegments has been called on pts).
*/
void SorClipper::splitIntoMonotonicSections(axom::ArrayView<const Point2DType> pts,
                                            axom::Array<axom::Array<Point2DType>>& sections)
{
  AXOM_ANNOTATE_SCOPE("SorClipper::splitIntoMonotonicSections");
  axom::Array<axom::IndexType> splitIdx =
    FSorClipper::findZSwitchbacks(pts);

  const axom::IndexType sectionCount = splitIdx.size() - 1;
  sections.clear();
  sections.resize(sectionCount);
  for(axom::IndexType i = 0; i < sectionCount; ++i)
  {
    axom::IndexType firstInSection = splitIdx[i];
    axom::IndexType lastInSection = splitIdx[i + 1];
    auto& curSection = sections[i];
    curSection.reserve(lastInSection - firstInSection + 1);
    for(axom::IndexType j = firstInSection; j <= lastInSection; ++j)
    {
      curSection.push_back(pts[j]);
    }
  }
}

// Compute a += b.
template <typename ExecSpace>
void accumulateDataImpl(axom::ArrayView<double> a,
                        axom::ArrayView<const double> b,
                        double scale)
{
  SLIC_ASSERT(a.size() == b.size());
  axom::for_all<ExecSpace>(
    a.size(),
    AXOM_LAMBDA(axom::IndexType i) { a[i] += scale * b[i]; });
}

void SorClipper::accumulateData(axom::ArrayView<double> a,
                                axom::ArrayView<const double> b,
                                double scale,
                                axom::runtime_policy::Policy runtimePolicy)
{
  using RuntimePolicy = axom::runtime_policy::Policy;
  if(runtimePolicy == RuntimePolicy::seq)
  {
    accumulateDataImpl<axom::SEQ_EXEC>(a, b, scale);
  }
#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
  else if(runtimePolicy == RuntimePolicy::omp)
  {
    accumulateDataImpl<axom::OMP_EXEC>(a, b, scale);
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
  else if(runtimePolicy == RuntimePolicy::cuda)
  {
    accumulateDataImpl<axom::CUDA_EXEC<256>>(a, b, scale);
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_HIP
  else if(runtimePolicy == RuntimePolicy::hip)
  {
    accumulateDataImpl<axom::HIP_EXEC<256>>(a, b, scale);
  }
#endif
  else
  {
    SLIC_ERROR(
      axom::fmt::format("Unrecognized runtime policy {}", runtimePolicy));
  }
  return;
}

void SorClipper::extractClipperInfo()
{
  auto sorOriginArray = m_info.fetch_existing("sorOrigin").as_double_array();
  auto sorDirectionArray = m_info.fetch_existing("sorDirection").as_double_array();
  for(int d = 0; d < 3; ++d)
  {
    m_sorOrigin[d] = sorOriginArray[d];
    m_sorDirection[d] = sorDirectionArray[d];
  }

  auto discreteFunctionArray = m_info.fetch_existing("discreteFunction").as_double_array();
  auto n = discreteFunctionArray.number_of_elements();

  SLIC_ERROR_IF(
    n % 2 != 0,
    axom::fmt::format(
      "***SorClipper: Discrete function must have an even number of values.  It has {}.",
      n));

  m_sorCurve.resize(axom::ArrayOptions::Uninitialized(), n / 2);
  for(int i = 0; i < n/2; ++i)
  {
    m_sorCurve[i] = Point2DType{discreteFunctionArray[i*2], discreteFunctionArray[i*2 + 1]};
  }

  m_levelOfRefinement = m_info.fetch_existing("levelOfRefinement").to_double();
}

}  // end namespace quest
}  // end namespace axom
