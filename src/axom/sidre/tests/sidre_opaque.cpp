// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include "axom/sidre.hpp"

using axom::sidre::DataStore;
using axom::sidre::DataType;
using axom::sidre::Group;
using axom::sidre::View;

//------------------------------------------------------------------------------
// Some simple types and functions used in tests
// (included in namespace to prevent clashes)
//------------------------------------------------------------------------------
namespace sidreopaquetest
{
enum Centering
{
  _Zone_,
  _Node_,
  _UnknownCentering_
};

enum DType
{
  _Double_,
  _Int_,
  _UnknownType_
};

class Extent
{
public:
  Extent(int lo, int hi) : m_ilo(lo), m_ihi(hi) { ; }

  int getNumPts(Centering cent) const
  {
    int retval = 0;

    switch(cent)
    {
    case _Zone_:
    {
      retval = (m_ihi - m_ilo + 1);
      break;
    }

    case _Node_:
    {
      retval = (m_ihi - m_ilo + 2);
      break;
    }

    default:
    {
      retval = -1;  // I know magic numbers are bad. So sue me.
    }

    }  // switch on centering

    return retval;
  }

  int m_ilo;
  int m_ihi;
};

class MeshVar
{
public:
  MeshVar(Centering cent, DType type, int depth = 1) : m_cent(cent), m_type(type), m_depth(depth)
  {
    ;
  }

  int getNumVals(const Extent* ext) const { return (ext->getNumPts(m_cent) * m_depth); }

  Centering m_cent;
  DType m_type;
  int m_depth;
};

}  // namespace sidreopaquetest

//------------------------------------------------------------------------------
//
// Simple test that adds an opaque data object, retrieves it and checks if
// the retrieved object is in the expected state.
//
TEST(sidre_opaque, basic_inout)
{
  using namespace sidreopaquetest;

  const int ihi_val = 9;

  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* problem_gp = root->createGroup("problem");

  Extent* ext = new Extent(0, ihi_val);

  View* ext_view = problem_gp->createView("ext", ext);

  EXPECT_TRUE(ext_view->isExternal());
  EXPECT_FALSE(ext_view->isApplied());
  EXPECT_TRUE(ext_view->isOpaque());
  EXPECT_EQ(ext, ext_view->getVoidPtr());

  Extent* test_extent = static_cast<Extent*>(ext_view->getVoidPtr());
  EXPECT_EQ(ext, test_extent);

  int test_ihi = test_extent->m_ihi;
  EXPECT_EQ(test_ihi, ihi_val);

  // Similar test with different view methods

  Extent* ext2 = new Extent(0, 2 * ihi_val);

  View* ext2_view = problem_gp->createView("ext2")->setExternalDataPtr(ext2);

  EXPECT_TRUE(ext2_view->isOpaque());
  EXPECT_EQ(ext2, ext2_view->getVoidPtr());

  Extent* test_extent2 = static_cast<Extent*>(ext2_view->getVoidPtr());
  EXPECT_EQ(test_extent2, ext2_view->getVoidPtr());

  int test_ihi2 = test_extent2->m_ihi;
  EXPECT_EQ(test_ihi2, 2 * ihi_val);

  // clean up...
  delete ext;
  delete ext2;
  delete ds;
}

//------------------------------------------------------------------------------
//
// Test that adds "MeshVars" as opaque data objects, creates views for their
// data on each of two domains, allocates their data (based on centering,
// domain size, and depth), and then checks to if the allocated data
// lengths match the expected values.
//
TEST(sidre_opaque, meshvar)
{
  using namespace sidreopaquetest;

  const int NUM_DOMAINS = 2;
  const int ilo_val[] = {0, 10};
  const int ihi_val[] = {9, 21};
  const std::string dom_name[] = {std::string("domain0"), std::string("domain1")};

  const int zone_var_depth = 1;
  const int node_var_depth = 2;

  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* problem_gp = root->createGroup("problem");

  // Add two different mesh vars to mesh var group
  Group* meshvar_gp = problem_gp->createGroup("mesh_var");
  MeshVar* zone_mv = new MeshVar(_Zone_, _Int_, zone_var_depth);
  View* zone_mv_view = meshvar_gp->createView("zone_mv", zone_mv);
  EXPECT_EQ(zone_mv, zone_mv_view->getVoidPtr());

  MeshVar* node_mv = new MeshVar(_Node_, _Double_, node_var_depth);
  View* node_mv_view = meshvar_gp->createView("node_mv", node_mv);
  EXPECT_EQ(node_mv, node_mv_view->getVoidPtr());

  //
  // Create domain groups, add extents
  // Create data views for mesh var data on domains and allocate
  //
  for(int idom = 0; idom < NUM_DOMAINS; ++idom)
  {
    Group* dom_gp = problem_gp->createGroup(dom_name[idom]);
    Extent* dom_ext = new Extent(ilo_val[idom], ihi_val[idom]);
    dom_gp->createView("ext", dom_ext);
    EXPECT_EQ(dom_ext, dom_gp->getView("ext")->getVoidPtr());

    MeshVar* zonemv = static_cast<MeshVar*>(zone_mv_view->getVoidPtr());
    View* dom_zone_view = dom_gp->createView("zone_data");
    dom_zone_view->allocate(DataType::c_int(zonemv->getNumVals(dom_ext)));

    MeshVar* nodemv = static_cast<MeshVar*>(node_mv_view->getVoidPtr());
    View* dom_node_view = dom_gp->createView("node_data");
    dom_node_view->allocate(DataType::c_double(nodemv->getNumVals(dom_ext)));
  }

  //
  //  Print datastore contents to see what's going on.
  //
  //  ds->print();

  //
  // Check data lengths
  //
  for(int idom = 0; idom < 2; ++idom)
  {
    Group* dom_gp = problem_gp->getGroup(dom_name[idom]);
    Extent* dom_ext = static_cast<Extent*>(dom_gp->getView("ext")->getVoidPtr());

    MeshVar* zonemv = static_cast<MeshVar*>(zone_mv_view->getVoidPtr());
    MeshVar* nodemv = static_cast<MeshVar*>(node_mv_view->getVoidPtr());

    int num_zone_vals = zonemv->getNumVals(dom_ext);
    int test_num_zone_vals = dom_gp->getView("zone_data")->getNumElements();
    EXPECT_EQ(num_zone_vals, test_num_zone_vals);

    int num_node_vals = nodemv->getNumVals(dom_ext);
    int test_num_node_vals = dom_gp->getView("node_data")->getNumElements();
    EXPECT_EQ(num_node_vals, test_num_node_vals);
  }

  // clean up...
  delete zone_mv;
  delete node_mv;
  for(int idom = 0; idom < 2; ++idom)
  {
    delete static_cast<Extent*>(problem_gp->getGroup(dom_name[idom])->getView("ext")->getVoidPtr());
  }
  delete ds;
}
