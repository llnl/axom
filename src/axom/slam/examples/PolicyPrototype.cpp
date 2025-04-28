// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//#include<string>
//#include<iostream>
#include <cstdlib>

namespace slamTemplateEx
{
struct Set
{
  //    std::string foo() { return "Base foo."; }
};

///  Offsetting policies

struct NoTranslation
{
  NoTranslation(int sz) : m_sz(sz) { }
  NoTranslation(int /*ignore*/, int sz) : m_sz(sz) { }
  inline int offset() const { return 0; }
  inline int size() const { return m_sz; }

  void setSize(int sz) { m_sz = sz; }

private:
  int m_sz;
};

struct HasTranslation
{
  HasTranslation(int lo, int hi) : m_lo(lo), m_hi(hi) { }

  inline int offset() const { return m_lo; }
  inline int size() const { return m_hi - m_lo; }

private:
  int m_lo, m_hi;
};

//// Indirection policies

struct NoIndirection
{
  inline int indirection(int pos) const { return pos; }
};

struct HasIndirection
{
  inline int indirection(int pos) const { return m_data[pos]; }
  int*& data() { return m_data; }

private:
  int* m_data;
};

/// Striding policies

struct NoStride
{
  inline int stride() const { return 1; }
};

template <int STRIDE>
struct FixedStride
{
  inline int stride() const { return STRIDE; }
};

struct RuntimeStride
{
  RuntimeStride(int stride) : m_stride(stride) { }
  inline int stride() const { return m_stride; }

private:
  int m_stride;
};

template <typename IndexType = int,
          typename TranslationPolicy = NoTranslation,
          typename IndirectionPolicy = NoIndirection,
          typename StridePolicy = NoStride>
struct OrderedSet : public Set, TranslationPolicy, IndirectionPolicy, StridePolicy
{
  using MyTranslationPolicy = TranslationPolicy;

  OrderedSet(int sz) : TranslationPolicy(0, sz) { }
  OrderedSet(int lo, int hi) : TranslationPolicy(lo, hi) { }

  inline int at(int pos) const { return indirection(pos * stride() + offset()); }

  inline int size() const { return TranslationPolicy::size(); }
  inline int offset() const { return TranslationPolicy::offset(); }
  inline int stride() const { return StridePolicy::stride(); }
  inline int indirection(int pos) const { return IndirectionPolicy::indirection(pos); }
};

////  Define concrete Set types as combinations of the above

struct PositionSet : OrderedSet<int>
{
  using ParentType = OrderedSet<int>;
  using MyTranslationPolicy = ParentType::MyTranslationPolicy;

  PositionSet(int n) : ParentType(n) { }
};

struct RangeSet : OrderedSet<int, HasTranslation>
{
  using ParentType = OrderedSet<int, HasTranslation>;
  using MyTranslationPolicy = ParentType::MyTranslationPolicy;

  RangeSet(int lo, int hi) : ParentType(lo, hi) { }
};

struct IndirectionSet : OrderedSet<int, NoTranslation, HasIndirection>
{
  using ParentType = OrderedSet<int, NoTranslation, HasIndirection>;
  using MyTranslationPolicy = ParentType::MyTranslationPolicy;

  IndirectionSet(int n) : ParentType(n) { }
};

}  // end namespace slamTemplateEx

//using namespace TemplateEx;
/*
   void printArray(int* arr, int sz)
   {
    int* off = arr;

    std::cout <<"\nPrinting position set..\n";
    for(int i=0; i< 10; ++i)
        std::cout<<"\ti: " << *off++ << std::endl;

    std::cout <<"\nPrinting range set..\n";
    for(int i=0; i< 10; ++i)
        std::cout<<"\ti: " << *off++ << std::endl;

    std::cout <<"\nPrinting indirection set..\n";
    for(int i=0; i< 10; ++i)
        std::cout<<"\ti: " << *off++ << std::endl;

   }
 */

using ResType = int;  // = long long int;

template <typename SetType>
inline ResType sumSet(const SetType& set)
{
  ResType sum = 0;

  for(int i = 0; i < set.size(); ++i)
  {
    sum += set.at(i);
  }

  return sum;
}

template <typename SetType>
inline void copySet(const SetType& set, int* buf)
{
  for(int i = 0; i < set.size(); ++i)
  {
    *buf++ = set.at(i);
  }
}

int main(int argc, char* argv[])
{
  // Process command line arguments
  // first argument is the size of the set
  // second is the begin index of the range set.
  int numElts = 10;
  int rangeBeginElt = 42;

  if(argc == 2)
  {
    numElts = std::atoi(argv[0]);
    rangeBeginElt = std::atoi(argv[1]);
  }

  // allocate and initialize the indirection array elements
  int* pVal = new int[numElts];
  for(int i = 0; i < numElts; ++i)
  {
    pVal[i] = i * i;
  }
  if(numElts > 0)
  {
    pVal[numElts - 1] = 12345;
  }

  slamTemplateEx::PositionSet pSet(numElts);
  slamTemplateEx::RangeSet rSet(rangeBeginElt, numElts + rangeBeginElt);
  slamTemplateEx::IndirectionSet iSet(numElts);
  iSet.data() = pVal;

  // Test 1 -- iterate through each set and find the sum
  ResType sumP = sumSet(pSet);
  ResType sumR = sumSet(rSet);
  ResType sumI = sumSet(iSet);

  bool test1 = (sumP < sumR && sumR < sumI);

  // Test 2 -- copy all three sets into a buffer
  int* buf = new int[3 * numElts];
  copySet(pSet, buf + numElts * 0);
  copySet(rSet, buf + numElts * 1);
  copySet(iSet, buf + numElts * 2);

  // We will now compare the three individual sums
  slamTemplateEx::IndirectionSet totIndSet(3 * numElts);
  totIndSet.data() = buf;
  bool test2 = (sumSet(totIndSet) == (sumP + sumR + sumI));

  //    std::cout <<"\nSum P: " << sumP
  //              <<"\nSum R: " << sumR
  //              <<"\nSum I: " << sumI
  //              <<"\nsumSum: " << sumP + sumR + sumI
  //              <<"\ntotal: " << sumSet(totIndSet)
  //              <<"\nTest1: " << ( test1 ? "true" : "false" )
  //              <<"\nTest2: " << ( test2 ? "true" : "false" )
  //              << std::endl;

  delete[] pVal;
  delete[] buf;

  return (test1 && test2) ? 0 : 1;
}
