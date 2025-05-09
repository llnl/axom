// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/// Tests axom's ItemCollection hierachy,
/// including MapCollection, ListCollection and IndexedCollection

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core.hpp"

#include "axom/fmt.hpp"

struct NamedItem
{
  NamedItem(const std::string& name) : m_name(name) { }
  const std::string& getName() const { return m_name; }

  friend bool operator==(const NamedItem& lhs, const NamedItem& rhs)
  {
    return lhs.getName() == rhs.getName();
  }

  std::string m_name;
};

template <typename CollectionType>
class ItemCollectionTest : public ::testing::Test
{
public:
  using ValueType = typename CollectionType::value_type;

  static constexpr bool IsNameBased =
    std::is_same<axom::MapCollection<ValueType>, CollectionType>::value;

protected:
  void SetUp() override { m_coll = new CollectionType; }

  void TearDown() override
  {
    // Remove and deallocate items
    for(auto idx = m_coll->getFirstValidIndex(); idx != axom::InvalidIndex;
        idx = m_coll->getNextValidIndex(idx))
    {
      auto* val = m_coll->removeItem(idx);
      delete val;
    }

    delete m_coll;
  }

  // Adds several initial items to the collection
  std::map<std::string, axom::IndexType> addItems(const std::vector<std::string>& strs)
  {
    std::map<std::string, axom::IndexType> insertionMap;

    for(auto& str : strs)
    {
      auto* val = this->template create_item<ValueType>(str);
      // elide warning about non-empty string in ListCollection
      std::string key = IsNameBased ? str : "";

      insertionMap[str] = m_coll->insertItem(val, key);
    }

    return insertionMap;
  }

protected:
  /// SFINAE function to create a double from a string
  /// Note: caller is responsible for deallocating the associated memory
  template <typename T>
  typename std::enable_if<std::is_same<double, T>::value, double*>::type create_item(
    const std::string& str) const
  {
    // creates arbitrary double from the input string's length
    return new double(str.size() * 1.111);
  }

  /// SFINAE function to create a NewItem from a string
  /// Note: caller is responsible for deallocating the associated memory
  template <typename T>
  typename std::enable_if<std::is_same<NamedItem, T>::value, NamedItem*>::type create_item(
    const std::string& str) const
  {
    return new NamedItem(str);
  }

protected:
  axom::ItemCollection<ValueType>* m_coll {nullptr};
};

using MyTypes =
  ::testing::Types<axom::IndexedCollection<double>,
                   axom::ListCollection<double>,
                   //axom::MapCollection<double>, // note: invalid since double doesn't have required getName()
                   axom::IndexedCollection<NamedItem>,
                   axom::ListCollection<NamedItem>,
                   axom::MapCollection<NamedItem>>;
TYPED_TEST_SUITE(ItemCollectionTest, MyTypes);

TYPED_TEST(ItemCollectionTest, testEmpty)
{
  auto* coll = this->m_coll;

  EXPECT_EQ(0, coll->getNumItems());
}

TYPED_TEST(ItemCollectionTest, insertItem)
{
  using ValueType = typename TestFixture::ValueType;

  auto* coll = this->m_coll;

  auto* val = this->template create_item<ValueType>("foo");
  std::string key = TestFixture::IsNameBased ? "foo" : "";

  coll->insertItem(val, key);
  EXPECT_EQ(1, coll->getNumItems());
}

TYPED_TEST(ItemCollectionTest, insertSeveralItems)
{
  using ValueType = typename TestFixture::ValueType;

  auto* coll = this->m_coll;

  for(auto& str : {"a", "b", "c", "aa", "bb"})
  {
    auto* val = this->template create_item<ValueType>(str);
    std::string key = TestFixture::IsNameBased ? str : "";

    coll->insertItem(val, key);
  }
  EXPECT_EQ(5, coll->getNumItems());
}

TYPED_TEST(ItemCollectionTest, insertManyItems)
{
  using ValueType = typename TestFixture::ValueType;
  const int SZ = 1000;
  auto* coll = this->m_coll;

  int num_added = 0;
  int num_removed = 0;

  for(char c = 'a'; c <= 'z'; ++c)
  {
    std::vector<axom::IndexType> indices;

    // insert SZ items
    for(int i = 0; i < SZ; ++i)
    {
      auto str = axom::fmt::format("{}_{:08}", c, i);
      auto* val = this->template create_item<ValueType>(str);
      std::string key = TestFixture::IsNameBased ? str : "";

      auto idx = coll->insertItem(val, key);
      indices.push_back(idx);
      ++num_added;
    }
    EXPECT_EQ(num_added - num_removed, coll->getNumItems());

    // remove a third of the new items
    for(std::size_t i = 0; i < indices.size(); ++i)
    {
      if(i % 3 == 0)
      {
        auto* val = coll->removeItem(indices[i]);
        delete val;
        ++num_removed;
      }
    }
    EXPECT_EQ(num_added - num_removed, coll->getNumItems());
  }
}

TYPED_TEST(ItemCollectionTest, removeNonExistent)
{
  using ValueType = typename TestFixture::ValueType;
  const int SZ = 100;
  auto* coll = this->m_coll;

  int num_added = 0;
  int num_removed = 0;

  std::vector<axom::IndexType> indices;

  // insert SZ items
  for(int i = 0; i < SZ; ++i)
  {
    auto str = axom::fmt::format("a_{:08}", i);
    auto* val = this->template create_item<ValueType>(str);
    std::string key = TestFixture::IsNameBased ? str : "";

    auto idx = coll->insertItem(val, key);
    indices.push_back(idx);
    ++num_added;
  }
  EXPECT_EQ(num_added - num_removed, coll->getNumItems());

  // remove a third of the items
  for(std::size_t i = 0; i < indices.size(); ++i)
  {
    EXPECT_TRUE(coll->hasItem(indices[i]));
    if(i % 3 == 0)
    {
      auto* val = coll->removeItem(indices[i]);
      delete val;
      ++num_removed;
    }
  }
  EXPECT_EQ(num_added - num_removed, coll->getNumItems());

  // attempt to remove same items
  for(std::size_t i = 0; i < indices.size(); ++i)
  {
    if(i % 3 == 0)
    {
      EXPECT_FALSE(coll->hasItem(indices[i]));
      auto* val = coll->removeItem(indices[i]);
      EXPECT_EQ(nullptr, val);
    }
    else
    {
      EXPECT_TRUE(coll->hasItem(indices[i]));
    }
  }
  EXPECT_EQ(num_added - num_removed, coll->getNumItems());
}

TYPED_TEST(ItemCollectionTest, hasItems)
{
  auto* coll = this->m_coll;

  std::vector<std::string> names {"a", "b", "c", "aa", "bb", "cc"};

  auto map = this->addItems(names);
  EXPECT_EQ(names.size(), coll->getNumItems());

  for(auto kv : map)
  {
    axom::IndexType idx = kv.second;
    EXPECT_TRUE(coll->hasItem(idx));
  }
}

TYPED_TEST(ItemCollectionTest, removeItems)
{
  auto* coll = this->m_coll;

  std::vector<std::string> names {"a", "b", "c", "aa", "bb", "cc"};

  auto map = this->addItems(names);
  EXPECT_EQ(names.size(), coll->getNumItems());

  for(auto kv : map)
  {
    axom::IndexType idx = kv.second;
    EXPECT_TRUE(coll->hasItem(idx));
    auto* val = coll->removeItem(idx);
    delete val;
    EXPECT_FALSE(coll->hasItem(idx));
  }
  EXPECT_EQ(0, coll->getNumItems());
}

TYPED_TEST(ItemCollectionTest, removeAndAddItems)
{
  auto* coll = this->m_coll;

  std::vector<std::string>
    names {"a", "b", "c", "aa", "bb", "cc", "aaa", "bbb", "ccc", "aaaa", "bbbb", "cccc"};
  std::vector<std::string> remove_names {"b", "bb", "bbb", "bbbb"};
  std::vector<std::string> keep_names {"a", "c", "aa", "cc", "aaa", "ccc", "aaaa", "cccc"};
  std::vector<std::string> add_names {"d", "dd", "ddd", "dddd"};

  auto map = this->addItems(names);
  EXPECT_EQ(names.size(), coll->getNumItems());

  // Remove items from 'remove' list
  for(auto n : remove_names)
  {
    if(map.find(n) != map.end())
    {
      axom::IndexType idx = map[n];
      auto* val = coll->getItem(idx);
      coll->removeItem(idx);
      delete val;
    }
  }

  // check that items from 'keep' list are still present and 'remove' are missing
  for(auto kv : map)
  {
    if(std::find(remove_names.begin(), remove_names.end(), kv.first) != remove_names.end())
    {
      EXPECT_FALSE(coll->hasItem(kv.second));
    }

    if(std::find(keep_names.begin(), keep_names.end(), kv.first) != keep_names.end())
    {
      EXPECT_TRUE(coll->hasItem(kv.second));
    }
  }

  // Add items from 'add' list
  auto secondMap = this->addItems(add_names);
  map.insert(secondMap.begin(), secondMap.end());

  // check that items from 'keep' and 'add' lists are still present
  // can no longer check remove since indices will be reused
  for(auto kv : map)
  {
    if(std::find(keep_names.begin(), keep_names.end(), kv.first) != keep_names.end())
    {
      EXPECT_TRUE(coll->hasItem(kv.second));
    }

    if(std::find(add_names.begin(), add_names.end(), kv.first) != add_names.end())
    {
      EXPECT_TRUE(coll->hasItem(kv.second));
    }
  }
}

TYPED_TEST(ItemCollectionTest, getItemsAndRemove)
{
  auto* coll = this->m_coll;

  std::vector<std::string> names {"a", "b", "c", "aa", "bb", "cc"};

  auto map = this->addItems(names);
  EXPECT_EQ(names.size(), coll->getNumItems());

  // Tests gettting items -- then remove it
  for(auto kv : map)
  {
    axom::IndexType idx = kv.second;
    auto* val = coll->getItem(idx);
    coll->removeItem(idx);
    delete val;
  }
  EXPECT_EQ(0, coll->getNumItems());
}

TYPED_TEST(ItemCollectionTest, getItemsAndRemoveAll)
{
  auto* coll = this->m_coll;

  std::vector<std::string> names {"a", "b", "c", "aa", "bb", "cc"};

  auto map = this->addItems(names);
  EXPECT_EQ(names.size(), coll->getNumItems());

  // Tests gettting items -- then remove it
  for(auto kv : map)
  {
    axom::IndexType idx = kv.second;
    auto* val = coll->getItem(idx);
    delete val;
  }
  EXPECT_EQ(names.size(), coll->getNumItems());

  coll->removeAllItems();
  EXPECT_EQ(0, coll->getNumItems());
}

TYPED_TEST(ItemCollectionTest, iterators)
{
  auto* coll = this->m_coll;

  std::vector<std::string> names {"a", "b", "c", "aa", "bb", "cc", "aaa", "bbb", "ccc"};

  auto map = this->addItems(names);
  EXPECT_EQ(names.size(), coll->getNumItems());

  // remove some items
  for(auto str : {"dddd", "ccc", "bb", "a"})
  {
    if(map.find(str) != map.end())
    {
      auto* val = coll->removeItem(map[str]);
      delete val;
    }
  }

  // add some items
  this->addItems({"dddd", "ddd", "dd", "d"});
  EXPECT_EQ(names.size() - 3 + 4, coll->getNumItems());

  // iterators
  for(auto it = coll->begin(), itEnd = coll->end(); it != itEnd; ++it)
  {
    axom::IndexType idx = it.index();
    EXPECT_TRUE(coll->hasItem(idx));
    EXPECT_EQ(*coll->getItem(idx), *it);
  }

  // const-iterators
  for(auto it = coll->cbegin(), itEnd = coll->cend(); it != itEnd; ++it)
  {
    axom::IndexType idx = it.index();
    EXPECT_TRUE(coll->hasItem(idx));
    EXPECT_EQ(*coll->getItem(idx), *it);
  }

  // iterators on const collection
  const auto* c_coll = coll;
  for(auto it = c_coll->begin(), itEnd = c_coll->cend(); it != itEnd; ++it)
  {
    axom::IndexType idx = it.index();
    EXPECT_TRUE(coll->hasItem(idx));
    EXPECT_EQ(*coll->getItem(idx), *it);
  }
}

// ----------------------------------------------------------------------------
// Adds tests specifically for MapCollection
// ----------------------------------------------------------------------------

template <typename TheValueType>
class MapCollectionTest : public ItemCollectionTest<axom::MapCollection<TheValueType>>
{
public:
  using ValueType = TheValueType;
  using MapCollectionType = axom::MapCollection<TheValueType>;
  using ItemCollectionBase = ItemCollectionTest<MapCollectionType>;

protected:
  void SetUp() override { ItemCollectionBase::SetUp(); }

  void TearDown() override { ItemCollectionBase::TearDown(); }

  MapCollectionType* getCollection() { return static_cast<MapCollectionType*>(this->m_coll); }
};

using MCollTypes = ::testing::Types<NamedItem>;
TYPED_TEST_SUITE(MapCollectionTest, MCollTypes);

TYPED_TEST(MapCollectionTest, testMapCollection)
{
  auto* map_coll = this->getCollection();
  if(map_coll != nullptr)
  {
    std::vector<std::string> names {"a", "b", "c", "aa", "bb", "cc", "aaa", "bbb", "ccc"};

    // add some items and check their properties
    this->addItems(names);
    EXPECT_EQ(names.size(), map_coll->getNumItems());

    for(auto& str : names)
    {
      EXPECT_TRUE(map_coll->hasItem(str));
      EXPECT_EQ(str, map_coll->getItem(str)->getName());

      auto idx = map_coll->getItemIndex(str);
      EXPECT_EQ(str, map_coll->getItemName(idx));
      EXPECT_EQ(*map_coll->getItem(idx), *map_coll->getItem(str));
    }

    // remove some items
    std::vector<std::string> removed_names {"bbb", "bb", "b"};
    for(auto& str : removed_names)
    {
      EXPECT_TRUE(map_coll->hasItem(str));
      auto* val = map_coll->removeItem(str);
      EXPECT_FALSE(map_coll->hasItem(str));
      delete val;
    }
    EXPECT_EQ(names.size() - removed_names.size(), map_coll->getNumItems());

    // iterate through remaining items
    for(auto& val : *map_coll)
    {
      const auto& name = val.getName();
      EXPECT_TRUE(map_coll->hasItem(name));

      const auto& idx = map_coll->getItemIndex(name);
      EXPECT_TRUE(map_coll->hasItem(idx));

      // item is in names, but not in removed_names
      EXPECT_NE(names.end(), std::find(names.begin(), names.end(), name));
      EXPECT_EQ(removed_names.end(), std::find(removed_names.begin(), removed_names.end(), name));
    }

    // add some items and iterate through collection
    {
      std::vector<std::string> added_names {"dddd", "ddd", "dd", "d"};
      this->addItems(added_names);
      auto amended_names = names;
      amended_names.insert(amended_names.end(), added_names.begin(), added_names.end());
      EXPECT_EQ(amended_names.size() - removed_names.size(), map_coll->getNumItems());

      // iterate through current items
      for(auto& val : *map_coll)
      {
        const auto& name = val.getName();
        EXPECT_TRUE(map_coll->hasItem(name));

        const auto& idx = map_coll->getItemIndex(name);
        EXPECT_TRUE(map_coll->hasItem(idx));

        // item is in (ammended) names, but not in removed_names
        EXPECT_NE(amended_names.end(), std::find(amended_names.begin(), amended_names.end(), name));
        EXPECT_EQ(removed_names.end(), std::find(removed_names.begin(), removed_names.end(), name));
      }
    }
  }
}

TYPED_TEST(MapCollectionTest, removeNonExistent)
{
  using ValueType = typename TestFixture::ValueType;
  const int SZ = 100;
  auto* map_coll = this->getCollection();

  if(map_coll != nullptr)
  {
    int num_added = 0;
    int num_removed = 0;

    std::vector<std::string> names_to_remove;

    // insert SZ items
    for(int i = 0; i < SZ; ++i)
    {
      auto str = axom::fmt::format("a_{:08}", i);
      auto* val = this->template create_item<ValueType>(str);
      map_coll->insertItem(val, str);
      ++num_added;

      if(i % 3 == 0)
      {
        names_to_remove.push_back(str);
      }
    }
    EXPECT_EQ(num_added - num_removed, map_coll->getNumItems());

    // remove a third of the items
    for(const auto& str : names_to_remove)
    {
      EXPECT_TRUE(map_coll->hasItem(str));
      auto* val = map_coll->removeItem(str);
      delete val;
      ++num_removed;
    }
    EXPECT_EQ(num_added - num_removed, map_coll->getNumItems());

    // attempt to remove same items
    for(const auto& str : names_to_remove)
    {
      EXPECT_FALSE(map_coll->hasItem(str));
      auto* val = map_coll->removeItem(str);
      EXPECT_EQ(nullptr, val);
    }
    EXPECT_EQ(num_added - num_removed, map_coll->getNumItems());
  }
}

TYPED_TEST(MapCollectionTest, insertAlreadyPresent)
{
  using ValueType = typename TestFixture::ValueType;
  const int SZ = 100;
  auto* map_coll = this->getCollection();

  if(map_coll != nullptr)
  {
    int num_added = 0;
    int num_removed = 0;

    std::vector<std::string> names_to_remove;
    std::vector<std::string> names_to_add;

    // insert SZ items
    for(int i = 0; i < SZ; ++i)
    {
      auto str = axom::fmt::format("a_{:08}", i);
      auto* val = this->template create_item<ValueType>(str);
      map_coll->insertItem(val, str);
      ++num_added;

      // create a list of items to remove and to add back
      switch(i % 3)
      {
      case 0:
        names_to_remove.push_back(str);
        names_to_add.push_back(str);
        break;
      case 1:
        names_to_add.push_back(str);
        break;
      case 2:
        //no-op
        break;
      }
    }
    EXPECT_EQ(num_added - num_removed, map_coll->getNumItems());

    // remove a third of the items
    for(const auto& str : names_to_remove)
    {
      EXPECT_TRUE(map_coll->hasItem(str));
      auto* val = map_coll->removeItem(str);
      delete val;
      ++num_removed;
      EXPECT_EQ(num_added - num_removed, map_coll->getNumItems());
    }
    EXPECT_EQ(num_added - num_removed, map_coll->getNumItems());

    // add back a bunch of items; some aleady present, others new
    for(const auto& str : names_to_add)
    {
      const bool hasItem = map_coll->hasItem(str);
      const bool wasRemoved =
        std::find(names_to_remove.begin(), names_to_remove.end(), str) != names_to_remove.end();
      EXPECT_NE(hasItem, wasRemoved);

      auto* val = this->template create_item<ValueType>(str);
      auto idx = map_coll->insertItem(val, str);

      if(hasItem)
      {
        // new item was not removed, must deallocate
        EXPECT_EQ(axom::InvalidIndex, idx);
        delete val;
        val = nullptr;
      }
      else
      {
        // new item was added
        EXPECT_NE(axom::InvalidIndex, idx);
        ++num_added;
      }
      EXPECT_TRUE(map_coll->hasItem(str));
      EXPECT_EQ(num_added - num_removed, map_coll->getNumItems());
    }
    EXPECT_EQ(num_added - num_removed, map_coll->getNumItems());
  }
}

// ----------------------------------------------------------------------------
// Adds tests specifically for IndexedCollection
// ----------------------------------------------------------------------------

template <typename TheValueType>
class IndexedCollectionTest : public ItemCollectionTest<axom::IndexedCollection<TheValueType>>
{
public:
  using ValueType = TheValueType;
  using IndexedCollectionType = axom::IndexedCollection<TheValueType>;
  using ItemCollectionBase = ItemCollectionTest<IndexedCollectionType>;

protected:
  void SetUp() override { ItemCollectionBase::SetUp(); }

  void TearDown() override { ItemCollectionBase::TearDown(); }

  IndexedCollectionType* getCollection()
  {
    return static_cast<IndexedCollectionType*>(this->m_coll);
  }
};

using ICollTypes = ::testing::Types<double, NamedItem>;
TYPED_TEST_SUITE(IndexedCollectionTest, ICollTypes);

TYPED_TEST(IndexedCollectionTest, testIndexedCollection)
{
  using ValueType = typename TestFixture::ValueType;

  auto* idx_coll = this->getCollection();
  if(idx_coll != nullptr)
  {
    std::vector<std::string> names {"a", "b", "c", "aa", "bb", "cc", "aaa", "bbb", "ccc"};

    // add some items and check their properties
    auto map = this->addItems(names);
    EXPECT_EQ(names.size(), idx_coll->getNumItems());

    auto idx_b = map["b"];
    auto idx_bb = map["bb"];
    auto idx_bbb = map["bbb"];

    // remove items by index
    for(auto rem_idx : {idx_b, idx_bb, idx_bbb})
    {
      EXPECT_TRUE(idx_coll->hasItem(rem_idx));
      auto* val = idx_coll->removeItem(rem_idx);
      delete val;
    }

    // attempt to remove items by index again
    for(auto rem_idx : {idx_b, idx_bb, idx_bbb})
    {
      EXPECT_FALSE(idx_coll->hasItem(rem_idx));
      EXPECT_EQ(nullptr, idx_coll->removeItem(rem_idx));
    }

    // get new indices without inserting anything; size should stay the same
    auto sz = idx_coll->getNumItems();
    for(int i = 0; i < 10; ++i)
    {
      idx_coll->getValidEmptyIndex();
      EXPECT_EQ(sz, idx_coll->getNumItems());
    }

    // add some new items
    for(auto str : {"ddddd", "dddd", "ddd", "dd", "d"})
    {
      auto idx = idx_coll->getValidEmptyIndex();
      EXPECT_FALSE(idx_coll->hasItem(idx));

      auto* val = this->template create_item<ValueType>(str);
      idx_coll->insertItem(val, idx);
      EXPECT_TRUE(idx_coll->hasItem(idx));
    }
    EXPECT_EQ(sz + 5, idx_coll->getNumItems());
  }
}

TYPED_TEST(IndexedCollectionTest, outOfOrderInsert)
{
  using ValueType = typename TestFixture::ValueType;

  auto* idx_coll = this->getCollection();
  if(idx_coll != nullptr)
  {
    std::vector<std::string> names {"a", "b", "c", "aa", "bb", "cc", "aaa", "bbb", "ccc"};

    // add some items and check their properties
    auto map = this->addItems(names);
    EXPECT_EQ(names.size(), idx_coll->getNumItems());

    const axom::IndexType idx_end = idx_coll->getLastAvailableEmptyIndex();
    EXPECT_EQ(idx_end, idx_coll->getNumItems());

    // remove some items by index
    const axom::IndexType idx_b = map["b"];
    const axom::IndexType idx_bb = map["bb"];
    const axom::IndexType idx_bbb = map["bbb"];
    for(auto rem_idx : {idx_b, idx_bb, idx_bbb})
    {
      EXPECT_TRUE(idx_coll->hasItem(rem_idx));
      auto* val = idx_coll->removeItem(rem_idx);
      delete val;
    }

    // Insert items into locations known to be empty but not at the top of the stack
    // Ensure that the size is correct
    {
      auto sz = idx_coll->getNumItems();
      int numAdded = 0;
      for(auto& pr : {std::make_pair(idx_b, this->template create_item<ValueType>("d")),
                      std::make_pair(idx_bb, this->template create_item<ValueType>("dd")),
                      std::make_pair(idx_end, this->template create_item<ValueType>("dddd")),
                      std::make_pair(idx_bbb, this->template create_item<ValueType>("ddd"))})
      {
        const auto idx = pr.first;
        auto* val = pr.second;
        EXPECT_FALSE(idx_coll->hasItem(idx));

        // Note: The following overload is unique to IndexedCollection
        idx_coll->insertItem(val, idx);
        EXPECT_TRUE(idx_coll->hasItem(idx));
        ++numAdded;
        EXPECT_EQ(sz + numAdded, idx_coll->getNumItems());
      }
      EXPECT_EQ(sz + 4, idx_coll->getNumItems());
    }
  }
}

TYPED_TEST(IndexedCollectionTest, insertAlreadyPresent)
{
  using ValueType = typename TestFixture::ValueType;
  const int SZ = 100;
  auto* indexed_coll = this->getCollection();

  if(indexed_coll != nullptr)
  {
    int num_added = 0;
    int num_removed = 0;

    std::vector<axom::IndexType> inds_to_remove;
    std::vector<axom::IndexType> inds_to_add;

    // insert SZ items
    for(int i = 0; i < SZ; ++i)
    {
      auto str = axom::fmt::format("a_{:08}", i);
      auto* val = this->template create_item<ValueType>(str);
      auto idx = indexed_coll->insertItem(val);
      ++num_added;

      // create a list of items to remove and to add back
      switch(i % 3)
      {
      case 0:
        inds_to_remove.push_back(idx);
        inds_to_add.push_back(idx);
        break;
      case 1:
        inds_to_add.push_back(idx);
        break;
      case 2:
        //no-op
        break;
      }
    }
    EXPECT_EQ(num_added - num_removed, indexed_coll->getNumItems());

    // remove a third of the items
    for(const auto& idx : inds_to_remove)
    {
      EXPECT_TRUE(indexed_coll->hasItem(idx));
      auto* val = indexed_coll->removeItem(idx);
      delete val;
      ++num_removed;
      EXPECT_EQ(num_added - num_removed, indexed_coll->getNumItems());
    }
    EXPECT_EQ(num_added - num_removed, indexed_coll->getNumItems());

    // add back a bunch of items; some aleady present, others new
    for(const auto& idx : inds_to_add)
    {
      const bool hasItem = indexed_coll->hasItem(idx);
      const bool wasRemoved =
        std::find(inds_to_remove.begin(), inds_to_remove.end(), idx) != inds_to_remove.end();
      EXPECT_NE(hasItem, wasRemoved);

      auto str = axom::fmt::format("a_{:08}", idx);
      auto* val = this->template create_item<ValueType>(str);
      auto newIndex = indexed_coll->insertItem(val, idx);

      if(hasItem)
      {
        // new item was not removed, must deallocate
        EXPECT_EQ(axom::InvalidIndex, newIndex);
        delete val;
        val = nullptr;
      }
      else
      {
        // new item was added
        EXPECT_NE(axom::InvalidIndex, newIndex);
        ++num_added;
      }
      EXPECT_TRUE(indexed_coll->hasItem(idx));
      EXPECT_EQ(num_added - num_removed, indexed_coll->getNumItems());
    }
    EXPECT_EQ(num_added - num_removed, indexed_coll->getNumItems());
  }
}

TYPED_TEST(IndexedCollectionTest, insertArbitraryIdx)
{
  using ValueType = typename TestFixture::ValueType;
  auto* indexed_coll = this->getCollection();

  if(indexed_coll != nullptr)
  {
    int num_added = 0;

    std::vector<axom::IndexType> indices {1, 10, 100, 1000, 500, 50, 5};
    for(auto idx : indices)
    {
      auto str = axom::fmt::format("a_{:08}", idx);
      auto* val = this->template create_item<ValueType>(str);
      auto insertedIdx = indexed_coll->insertItem(val, idx);
      ++num_added;

      EXPECT_EQ(insertedIdx, idx);
      EXPECT_EQ(num_added, indexed_coll->getNumItems());
    }

    EXPECT_EQ(indices.size(), indexed_coll->getNumItems());
  }
}

TYPED_TEST(IndexedCollectionTest, insertBadIdx)
{
  using ValueType = typename TestFixture::ValueType;
  auto* indexed_coll = this->getCollection();

  if(indexed_coll != nullptr)
  {
    int num_added = 0;

    std::vector<axom::IndexType> indices {1, -10, 100, -1000, 500, -50, 5};
    for(auto idx : indices)
    {
      auto str = axom::fmt::format("a_{:08}", idx);
      auto* val = this->template create_item<ValueType>(str);
      auto insertedIdx = indexed_coll->insertItem(val, idx);
      if(idx < 0)
      {
        EXPECT_EQ(axom::InvalidIndex, insertedIdx);
        delete val;
      }
      else
      {
        EXPECT_EQ(idx, insertedIdx);
        ++num_added;
      }
      EXPECT_EQ(num_added, indexed_coll->getNumItems());
    }

    EXPECT_EQ(num_added, indexed_coll->getNumItems());
  }
}
