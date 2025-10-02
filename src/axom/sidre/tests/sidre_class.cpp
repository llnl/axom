// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include "axom/sidre.hpp"
#include <vector>

using axom::sidre::Buffer;
using axom::sidre::DataStore;
using axom::sidre::DataType;
using axom::sidre::Group;
using axom::sidre::View;

namespace classtest
{
class Class1
{
public:
  Class1() { }

  explicit Class1(size_t len)
  {
    m_idata = std::vector<int>(len);

    for(size_t ii = 0; ii < len; ++ii)
    {
      m_idata[ii] = 3 * ii;
    }
  }

  std::vector<int>& getIData() { return m_idata; }

  void copyToGroup(Group* gp)
  {
    gp->createView("idata", &m_idata[0])->apply(DataType::c_int(m_idata.size()));
  }

  void copyFromGroup(Group* gp)
  {
    View* iview = gp->getView("idata");
    size_t ilen = iview->getNumElements();
    m_idata = std::vector<int>(ilen);

    int* g_idata = iview->getData();
    for(size_t ii = 0; ii < ilen; ++ii)
    {
      m_idata[ii] = g_idata[ii];
    }
  }

  void checkState(const Class1& class1) { checkState(&(class1.m_idata[0])); }

  void checkState(Group* gp)
  {
    int* idata_chk = gp->getView("idata")->getData();
    checkState(idata_chk);
  }

private:
  void checkState(const int* tidata)
  {
    for(size_t ii = 0; ii < m_idata.size(); ++ii)
    {
      EXPECT_EQ(m_idata[ii], tidata[ii]);
    }
  }

  std::vector<int> m_idata;
};

class Class2
{
public:
  Class2() { }

  explicit Class2(size_t len)
  {
    m_idata = std::vector<int>(len);
    m_ddata = std::vector<double>(len);
    m_class1 = Class1(len);

    for(size_t ii = 0; ii < len; ++ii)
    {
      m_idata[ii] = ii;
      m_ddata[ii] = 2.0 * m_idata[ii];
    }
  }

  std::vector<int>& getIData() { return m_idata; }
  std::vector<double>& getDData() { return m_ddata; }
  Class1& getClass1() { return m_class1; }

  void copyToGroup(Group* gp)
  {
    gp->createView("idata", &m_idata[0])
      ->
      //     apply(DataType::c_int(m_idata.size()));
      apply(axom::sidre::INT_ID, m_idata.size());
    gp->createView("ddata", &m_ddata[0])
      ->
      //     apply(DataType::c_double(m_ddata.size()));
      apply(axom::sidre::DOUBLE_ID, m_ddata.size());

    Group* gp1 = gp->createGroup("myclass1");

    m_class1.copyToGroup(gp1);
  }

  void copyFromGroup(Group* gp)
  {
    View* iview = gp->getView("idata");
    size_t ilen = iview->getBuffer()->getTotalBytes() / sizeof(CONDUIT_NATIVE_INT);
    m_idata = std::vector<int>(ilen);

    int* g_idata = iview->getData();
    for(size_t ii = 0; ii < ilen; ++ii)
    {
      m_idata[ii] = g_idata[ii];
    }

    View* dview = gp->getView("ddata");
    size_t dlen = dview->getNumElements();
    m_ddata = std::vector<double>(dlen);

    double* g_ddata = dview->getData();
    for(size_t ii = 0; ii < dlen; ++ii)
    {
      m_ddata[ii] = g_ddata[ii];
    }

    Group* gp1 = gp->getGroup("myclass1");
    m_class1.copyFromGroup(gp1);
  }

  void checkState(const Class2& class2)
  {
    checkState(&(class2.m_idata[0]), &(class2.m_ddata[0]));
    m_class1.checkState(class2.m_class1);
  }

  void checkState(Group* gp)
  {
    int* idata_chk = gp->getView("idata")->getData();
    double* ddata_chk = gp->getView("ddata")->getData();
    checkState(idata_chk, ddata_chk);

    Group* gp1 = gp->getGroup("myclass1");
    m_class1.checkState(gp1);
  }

private:
  void checkState(const int* tidata, const double* tddata)
  {
    for(size_t ii = 0; ii < m_idata.size(); ++ii)
    {
      EXPECT_EQ(m_idata[ii], tidata[ii]);
    }
    for(size_t ii = 0; ii < m_ddata.size(); ++ii)
    {
      EXPECT_EQ(m_ddata[ii], tddata[ii]);
    }
  }

  std::vector<int> m_idata;
  std::vector<double> m_ddata;
  Class1 m_class1;
};

}  // end namespace classtest

//------------------------------------------------------------------------------
// Test copying class data to group hierarchy
//------------------------------------------------------------------------------
TEST(sidre_class, class_to_group)
{
  using namespace classtest;

  const size_t len = 13;

  Class2 myclass2(len);

  EXPECT_EQ(myclass2.getIData().size(), len);
  EXPECT_EQ(myclass2.getDData().size(), len);
  EXPECT_EQ(myclass2.getClass1().getIData().size(), len);

  myclass2.checkState(myclass2);

  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  Group* c2group = root->createGroup("myclass2");
  myclass2.copyToGroup(c2group);

  ds->print();

  EXPECT_EQ(c2group->getNumViews(), 2u);
  EXPECT_EQ(c2group->getNumGroups(), 1u);

  myclass2.checkState(c2group);

  delete ds;
}

// TODO - need to fix save/load
#if 0
//------------------------------------------------------------------------------
// Test save/load class data using group hierarchy
//------------------------------------------------------------------------------
TEST(sidre_class, save_load_class_to_group)
{
  using namespace classtest;

  const size_t len = 21;

  Class2 myclass2(len);

  DataStore* ds   = new DataStore();
  Group* root = ds->getRoot();

  Group* c2group = root->createGroup("myclass2");
  myclass2.copyToGroup(c2group);

  EXPECT_EQ(c2group->getNumViews(), 2u);
  EXPECT_EQ(c2group->getNumGroups(), 1u);

  myclass2.checkState(c2group);

  ds->print();

  ds->getRoot()->save("out_save_load_class_to_group", "conduit");


  DataStore* ds2 = new DataStore();
  ds2->getRoot()->load("out_save_load_class_to_group","conduit");

  ds2->print();

  Group* load_myclass2 = ds2->getRoot()->getGroup("myclass2");
  Class2 load_class2;

  load_class2.copyFromGroup(load_myclass2);

  myclass2.checkState(load_class2);

  delete ds;
  delete ds2;
}
#endif
