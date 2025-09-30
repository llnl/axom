// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"      // for axom compile-time definitions
#include "axom/core/Types.hpp"  // for nullptr

#include "axom/mint/config.hpp"  // for mint compile-time type definitions
#include "axom/mint/mesh/FieldData.hpp"

#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre.hpp"
#endif

namespace axom
{
namespace mint
{
//------------------------------------------------------------------------------
// INTERNAL HELPER METHODS
//------------------------------------------------------------------------------
namespace internal
{
#ifdef AXOM_MINT_USE_SIDRE
mint::Field* getFieldFromView(const std::string& name, sidre::View* view)
{
  SLIC_ASSERT(view != nullptr);
  SLIC_ASSERT(!view->isEmpty());

  using int32 = std::int32_t;
  using int64 = std::int64_t;

  mint::Field* f = nullptr;

  switch(view->getTypeID())
  {
  case sidre::INT32_ID:
    f = new mint::FieldVariable<int32>(name, view);
    break;
  case sidre::INT64_ID:
    f = new mint::FieldVariable<int64>(name, view);
    break;
  case sidre::FLOAT64_ID:
    f = new mint::FieldVariable<double>(name, view);
    break;
  case sidre::FLOAT32_ID:
    f = new mint::FieldVariable<float>(name, view);
    break;
  default:
    SLIC_ERROR("Encountered unsupported type [" << view->getTypeID() << "]");
  }  // END switch

  SLIC_ERROR_IF(f == nullptr, "null field!");
  return (f);
}

//------------------------------------------------------------------------------
void removeFromSidre(sidre::Group* grp, const std::string& name)
{
  SLIC_ASSERT(grp != nullptr);
  SLIC_ASSERT(grp->hasChildGroup(name));

  grp->destroyGroup(name);
}

#endif

}  // namespace internal

//------------------------------------------------------------------------------
//  FIELDDATA IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
FieldData::FieldData(int association)
  : m_association(association)
  , m_resize_ratio(axom::deprecated::MCArray<double>::DEFAULT_RESIZE_RATIO)
  , m_fields()
{
#ifdef AXOM_MINT_USE_SIDRE
  m_fields_group = nullptr;
#endif

  SLIC_ERROR_IF((m_association < 0) || (m_association >= NUM_FIELD_ASSOCIATIONS),
                "Invalid field association!");
}

//------------------------------------------------------------------------------
#ifdef AXOM_MINT_USE_SIDRE
FieldData::FieldData(int association, sidre::Group* fields_group, const std::string& topo)
  : m_association(association)
  , m_resize_ratio(axom::deprecated::MCArray<double>::DEFAULT_RESIZE_RATIO)
  , m_fields()
  , m_fields_group(fields_group)
  , m_topology(topo)
{
  SLIC_ERROR_IF((m_association < 0) || (m_association >= NUM_FIELD_ASSOCIATIONS),
                "Invalid field association!");

  SLIC_ERROR_IF(m_fields_group == nullptr, "NULL sidre group!");

  size_t numGroups = m_fields_group->getNumGroups();
  for(size_t i = 0; i < numGroups; ++i)
  {
    sidre::Group* gp = m_fields_group->getGroup(i);
    SLIC_ERROR_IF(gp == nullptr, "Encountered a NULL group");

    SLIC_ERROR_IF(!gp->hasChildView("topology"),
                  "field [" << gp->getName() << "] does not conform to blueprint!"
                            << " Missing 'topology' view");
    SLIC_ERROR_IF(!gp->getView("topology")->isString(), "topology view needs to hold a string.");
    if(gp->getView("topology")->getString() != m_topology)
    {
      continue;
    }

    SLIC_ERROR_IF(!gp->hasChildView("association"),
                  "field [" << gp->getName() << "] does not conform to blueprint!"
                            << " Missing 'association' view");
    SLIC_ERROR_IF(!gp->getView("association")->isString(),
                  "association view needs to hold a string.");

    SLIC_ERROR_IF(!gp->hasChildView("volume_dependent"),
                  "field [" << gp->getName() << "] does not conform to blueprint!"
                            << " Missing 'volume_dependent' view");
    SLIC_ERROR_IF(!gp->getView("volume_dependent")->isString(),
                  "volume_dependent view needs to hold a string.");

    SLIC_ERROR_IF(!gp->hasChildView("values"),
                  "field [" << gp->getName() << "] does not conform to blueprint!"
                            << " Missing 'values' view");

    // NOTE: currently the blue-print supports
    const char* assoc = gp->getView("association")->getString();
    const bool isVertex = (strcmp(assoc, "vertex") == 0);
    const bool isElement = (strcmp(assoc, "element") == 0);
    const bool isFace = (strcmp(assoc, "face") == 0);
    const bool isEdge = (strcmp(assoc, "edge") == 0);
    SLIC_ERROR_IF(!isVertex && !isElement && !isFace && !isEdge,
                  "field [" << gp->getName() << "] has invalid association!"
                            << " => association= " << assoc);

    int centering;
    if(isVertex)
    {
      centering = NODE_CENTERED;
    }
    else if(isElement)
    {
      centering = CELL_CENTERED;
    }
    else if(isFace)
    {
      centering = FACE_CENTERED;
    }
    else
    {
      centering = EDGE_CENTERED;
    }

    IndexType num_tuples = -1;
    if(centering == m_association)
    {
      const std::string name = gp->getName();
      SLIC_ERROR_IF(hasField(name), "Encountered a duplicate field!");

      sidre::View* view = gp->getView("values");
      Field* field = internal::getFieldFromView(name, view);
      if(num_tuples == -1)
      {
        num_tuples = field->getNumTuples();
      }

      SLIC_ERROR_IF(field->getNumTuples() != num_tuples, "Inconsistent number of tuples");

      m_fields[name] = field;
    }  // END if centering
  }  // END for all fields
}
#endif

//------------------------------------------------------------------------------
void FieldData::clear()
{
  const int numFields = getNumFields();
  for(int i = 0; i < numFields; ++i)
  {
    delete getField(i);
  }

  m_fields.clear();
}

//------------------------------------------------------------------------------
void FieldData::resize(IndexType newNumTuples)
{
  const IndexType numFields = getNumFields();
  for(int i = 0; i < numFields; ++i)
  {
    getField(i)->resize(newNumTuples);
  };
}

//------------------------------------------------------------------------------
void FieldData::emplace(IndexType pos, IndexType num_tuples)
{
  const IndexType numFields = getNumFields();
  for(int i = 0; i < numFields; ++i)
  {
    getField(i)->emplace(pos, num_tuples);
  };
}

//------------------------------------------------------------------------------
void FieldData::reserve(IndexType newCapacity)
{
  const IndexType numFields = getNumFields();
  for(int i = 0; i < numFields; ++i)
  {
    getField(i)->reserve(newCapacity);
  };
}

//------------------------------------------------------------------------------
void FieldData::shrink()
{
  const IndexType numFields = getNumFields();
  for(int i = 0; i < numFields; ++i)
  {
    getField(i)->shrink();
  };
}

//------------------------------------------------------------------------------
void FieldData::setResizeRatio(double ratio)
{
  m_resize_ratio = ratio;
  const IndexType numFields = getNumFields();
  for(int i = 0; i < numFields; ++i)
  {
    getField(i)->setResizeRatio(ratio);
  };
}

//------------------------------------------------------------------------------
bool FieldData::checkConsistency(IndexType num_tuples, IndexType capacity) const
{
  const int numFields = getNumFields();
  if(numFields == 0)
  {
    return true;
  }

  bool tuple_status = true;
  bool capacity_status = true;
  bool resize_status = true;
  for(int i = 1; i < numFields; ++i)
  {
    const Field* f = getField(i);
    tuple_status &= f->getNumTuples() == num_tuples;
    capacity_status &= f->getCapacity() >= f->getNumTuples();
    if(!f->isExternal())
    {
      capacity_status &= f->getCapacity() == capacity;
      resize_status &= f->getResizeRatio() == m_resize_ratio;
    }
  }

  SLIC_WARNING_IF(!tuple_status, "Inconsistent number of tuples.");
  SLIC_WARNING_IF(!capacity, "Inconsistent capacity.");
  SLIC_WARNING_IF(!resize_status, "Inconsistent resize ratio.");
  return tuple_status && capacity_status && resize_status;
}

//------------------------------------------------------------------------------
void FieldData::removeField(const std::string& name)
{
  mint::Field* f = getField(name);
  SLIC_ERROR_IF(f == nullptr, "field [" << name << "] does not exist!");

  m_fields.erase(name);
  delete f;

#ifdef AXOM_MINT_USE_SIDRE
  if(hasSidreGroup() && m_fields_group->hasChildGroup(name))
  {
    internal::removeFromSidre(m_fields_group, name);
    SLIC_ASSERT(!m_fields_group->hasChildGroup(name));
  }
#endif
}

//------------------------------------------------------------------------------
void FieldData::removeField(int i)
{
  mint::Field* f = getField(i);
  removeField(f->getName());
  f = nullptr;
}

} /* namespace mint */
} /* namespace axom */
