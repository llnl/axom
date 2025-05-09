// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/Path.hpp"

TEST(core_Path, trivial_construct)
{
  const std::string path_str = "foo";
  axom::Path path(path_str);
  EXPECT_EQ(static_cast<std::string>(path), path_str);
}

TEST(core_Path, lead_delim_construct)
{
  const std::string path_str = "/foo";
  axom::Path path(path_str);
  EXPECT_EQ(static_cast<std::string>(path), path_str);
}

TEST(core_Path, two_component_construct)
{
  const std::string path_str = "foo/bar";
  axom::Path path(path_str);
  EXPECT_EQ(static_cast<std::string>(path), path_str);
}

TEST(core_Path, two_component_construct_alt_delim)
{
  const std::string path_str = "foo:bar";
  axom::Path path(path_str, ':');
  EXPECT_EQ(static_cast<std::string>(path), path_str);
}

TEST(core_Path, two_component_construct_lead_delim)
{
  const std::string path_str = "/foo/bar";
  axom::Path path(path_str);
  EXPECT_EQ(static_cast<std::string>(path), path_str);
}

TEST(core_Path, two_component_construct_alt_lead_delim)
{
  const std::string path_str = ":foo:bar";
  axom::Path path(path_str, ':');
  EXPECT_EQ(static_cast<std::string>(path), path_str);
}

TEST(core_Path, no_component_construct_root_delim)
{
  const std::string path_str = "/";
  axom::Path path(path_str);
  EXPECT_EQ(static_cast<std::string>(path), path_str);
}

TEST(core_Path, join_trivial)
{
  axom::Path first("foo");
  axom::Path second("bar");
  auto joined = axom::Path::join({first, second});
  EXPECT_EQ(static_cast<std::string>(joined), "foo/bar");
}

TEST(core_Path, join_trivial_alt_delim)
{
  axom::Path first("foo");
  axom::Path second("bar");
  auto joined = axom::Path::join({first, second}, ':');
  EXPECT_EQ(static_cast<std::string>(joined), "foo:bar");
}

TEST(core_Path, join_trivial_lead_delim)
{
  axom::Path first("/foo");
  axom::Path second("/bar");
  auto joined = axom::Path::join({first, second});
  EXPECT_EQ(static_cast<std::string>(joined), "/foo/bar");
}

TEST(core_Path, join_trivial_alt_lead_delim)
{
  axom::Path first(":foo", ':');
  axom::Path second(":bar", ':');
  auto joined = axom::Path::join({first, second}, ':');
  EXPECT_EQ(static_cast<std::string>(joined), ":foo:bar");
}

TEST(core_Path, join_multi_component)
{
  axom::Path first("foo/bar");
  axom::Path second("baz/quux");
  auto joined = axom::Path::join({first, second});
  EXPECT_EQ(static_cast<std::string>(joined), "foo/bar/baz/quux");
}

TEST(core_Path, join_multi_component_alt_delim)
{
  axom::Path first("foo:bar", ':');
  axom::Path second("baz:quux", ':');
  auto joined = axom::Path::join({first, second}, ':');
  EXPECT_EQ(static_cast<std::string>(joined), "foo:bar:baz:quux");
}

TEST(core_Path, join_multi_component_lead_delim)
{
  axom::Path first("/foo/bar");
  axom::Path second("/baz/quux");
  auto joined = axom::Path::join({first, second});
  EXPECT_EQ(static_cast<std::string>(joined), "/foo/bar/baz/quux");
}

TEST(core_Path, join_multi_component_alt_lead_delim)
{
  axom::Path first(":foo:bar", ':');
  axom::Path second(":baz:quux", ':');
  auto joined = axom::Path::join({first, second}, ':');
  EXPECT_EQ(static_cast<std::string>(joined), ":foo:bar:baz:quux");
}

TEST(core_Path, parent_basic)
{
  axom::Path path("foo/bar");
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), "foo");
}

TEST(core_Path, parent_basic_alt_delim)
{
  axom::Path path("foo;bar", ';');
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), "foo");
}

TEST(core_Path, parent_basic_lead_delim)
{
  axom::Path path("/foo/bar");
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), "/foo");
}

TEST(core_Path, parent_basic_alt_lead_delim)
{
  axom::Path path(":foo:bar", ':');
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), ":foo");
}

TEST(core_Path, parent_basic_root_lead_delim)
{
  axom::Path path("/foo");
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), "/");
}

TEST(core_Path, parent_multi_component)
{
  axom::Path path("foo/bar/baz");
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), "foo/bar");
}

TEST(core_Path, parent_multi_component_alt_delim)
{
  axom::Path path("foo;bar;baz", ';');
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), "foo;bar");
}

TEST(core_Path, parent_multi_component_lead_delim)
{
  axom::Path path("/foo/bar/baz");
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), "/foo/bar");
}

TEST(core_Path, parent_multi_component_alt_lead_delim)
{
  axom::Path path("|foo|bar|baz", '|');
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), "|foo|bar");
}

TEST(core_Path, dirname_basename_basic)
{
  axom::Path path("foo/bar");
  EXPECT_EQ(path.dirName(), "foo");
  EXPECT_EQ(path.baseName(), "bar");
}

TEST(core_Path, dirname_basename_basic_alt_delim)
{
  axom::Path path("foo;bar", ';');
  EXPECT_EQ(path.dirName(), "foo");
  EXPECT_EQ(path.baseName(), "bar");
}

TEST(core_Path, dirname_basename_basic_lead_delim)
{
  axom::Path path("/foo/bar");
  EXPECT_EQ(path.dirName(), "/foo");
  EXPECT_EQ(path.baseName(), "bar");
}

TEST(core_Path, dirname_basename_basic_alt_lead_delim)
{
  axom::Path path("+foo+bar", '+');
  EXPECT_EQ(path.dirName(), "+foo");
  EXPECT_EQ(path.baseName(), "bar");
}

TEST(core_Path, dirname_basename_root_lead_delim)
{
  axom::Path path("/foo");
  EXPECT_EQ(path.dirName(), "/");
  EXPECT_EQ(path.baseName(), "foo");
}

TEST(core_Path, dirname_basename_multi_component)
{
  axom::Path path("foo/bar/baz");
  EXPECT_EQ(path.dirName(), "foo/bar");
  EXPECT_EQ(path.baseName(), "baz");
}

TEST(core_Path, dirname_basename_multi_component_alt_delim)
{
  axom::Path path("foo;bar;baz", ';');
  EXPECT_EQ(path.dirName(), "foo;bar");
  EXPECT_EQ(path.baseName(), "baz");
}

TEST(core_Path, dirname_basename_multi_component_lead_delim)
{
  axom::Path path("/foo/bar/baz");
  EXPECT_EQ(path.dirName(), "/foo/bar");
  EXPECT_EQ(path.baseName(), "baz");
}

TEST(core_Path, dirname_basename_multi_component_alt_lead_delim)
{
  axom::Path path("|foo|bar|baz", '|');
  EXPECT_EQ(path.dirName(), "|foo|bar");
  EXPECT_EQ(path.baseName(), "baz");
}

TEST(core_Path, basic_equals)
{
  axom::Path first("foo/bar/baz");
  axom::Path second("foo/bar/baz");
  axom::Path third("foo:bar:baz", ':');
  axom::Path fourth("/foo/bar/baz");
  axom::Path fifth(":foo:bar:baz", ':');
  EXPECT_TRUE(first == second);
  EXPECT_FALSE(first != second);
  EXPECT_FALSE(first == third);
  EXPECT_TRUE(first != third);
  EXPECT_FALSE(first == fourth);
  EXPECT_TRUE(first != fourth);
  EXPECT_FALSE(first == fifth);
  EXPECT_FALSE(third == fifth);
  EXPECT_TRUE(third != fifth);
}

TEST(core_Path, parts_empty)
{
  axom::Path first("");
  axom::Path second("/");
  axom::Path third(":", ':');

  EXPECT_EQ(first.parts().size(), 0);
  EXPECT_EQ(second.parts().size(), 0);
  EXPECT_EQ(third.parts().size(), 0);
}

TEST(core_Path, parts_single)
{
  axom::Path first("foo");
  axom::Path second("/foo");
  axom::Path third("foo:", ':');

  EXPECT_EQ(first.parts().size(), 1);
  EXPECT_EQ(first.parts()[0], "foo");
  EXPECT_EQ(second.parts().size(), 1);
  EXPECT_EQ(second.parts()[0], "foo");
  EXPECT_EQ(third.parts().size(), 1);
  EXPECT_EQ(third.parts()[0], "foo");
}

TEST(core_Path, parts_multipart)
{
  axom::Path path("foo/bar/baz");
  auto parts = path.parts();

  EXPECT_EQ(parts.size(), 3);
  EXPECT_EQ(parts[0], "foo");
  EXPECT_EQ(parts[1], "bar");
  EXPECT_EQ(parts[2], "baz");
}

TEST(core_Path, parts_multipart_empty)
{
  axom::Path path("foo//bar/baz");
  auto parts = path.parts();

  EXPECT_EQ(parts.size(), 3);
  EXPECT_EQ(parts[0], "foo");
  EXPECT_EQ(parts[1], "bar");
  EXPECT_EQ(parts[2], "baz");
}

TEST(core_Path, basic_split)
{
  axom::Path path("foo/bar/baz");
  auto pair = path.split();

  EXPECT_EQ(pair.first, "foo/bar");
  EXPECT_EQ(pair.second, "baz");
}

TEST(core_Path, basic_split_lead_delim)
{
  axom::Path path("/foo/bar/baz");
  auto pair = path.split();

  EXPECT_EQ(pair.first, "/foo/bar");
  EXPECT_EQ(pair.second, "baz");
}

TEST(core_Path, basic_split_dirname_basename)
{
  axom::Path path("foo/bar/baz");
  auto pair = path.split();

  EXPECT_EQ(pair.first, path.dirName());
  EXPECT_EQ(pair.second, path.baseName());
}

TEST(core_Path, basic_split_dirname_basename_lead_delim)
{
  axom::Path path("/foo/bar/baz");
  auto pair = path.split();

  EXPECT_EQ(pair.first, path.dirName());
  EXPECT_EQ(pair.second, path.baseName());
}

TEST(core_Path, no_split)
{
  axom::Path path("foo");
  auto pair = path.split();

  EXPECT_EQ(pair.first, "");
  EXPECT_EQ(pair.second, "foo");
}

TEST(core_Path, split_empty_basename)
{
  axom::Path path("foo/");
  auto pair = path.split();

  EXPECT_EQ(pair.first, "");
  EXPECT_EQ(pair.second, "foo");
}

TEST(core_Path, split_root_dirname)
{
  axom::Path path("/foo");
  auto pair = path.split();

  EXPECT_EQ(pair.first, "/");
  EXPECT_EQ(pair.second, "foo");
}

TEST(core_Path, split_empty)
{
  axom::Path path("");
  auto pair = path.split();

  EXPECT_EQ(pair.first, "");
  EXPECT_EQ(pair.second, "");
}
