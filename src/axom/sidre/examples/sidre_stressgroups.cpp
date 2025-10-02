// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**************************************************************************
 *************************************************************************/

#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/utilities/Timer.hpp"
#include "axom/sidre.hpp"

#include <algorithm>
#include <random>
#include <set>

using axom::sidre::DataStore;
using axom::sidre::Group;
using axom::sidre::IndexType;

/**************************************************************************
 * Subroutine:  main
 * Purpose   :
 *************************************************************************/

int main(int argc, char* argv[])
{
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  size_t num_groups = 0;
  if(argc > 1)
  {
    num_groups = static_cast<IndexType>(atoi(argv[1]));
  }
  else
  {
    delete ds;
    return 0;
  }

  static const char alphanum[] =
    "0123456789"
    "!@#$%^&*"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz";

  int num_chars = sizeof(alphanum) - 1;

  std::set<std::string> name_set;

  while(name_set.size() < num_groups)
  {
    std::string new_string;
    for(int c = 0; c < 8; ++c)
    {
      new_string += alphanum[rand() % num_chars];
    }
    name_set.insert(new_string);
  }

  std::vector<std::string> names(num_groups);
  {
    int i = 0;
    for(std::set<std::string>::const_iterator itr = name_set.begin(); itr != name_set.end(); ++itr)
    {
      names[i] = *itr;
      ++i;
    }
  }

  /*
    for (int i = 0; i < num_groups; ++i) {
       std::stringstream sstr;
       sstr << i;

       names[i] = sstr.str();;
    }
 */
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(names.begin(), names.end(), g);

  axom::utilities::Timer create_timer(true);
  for(size_t i = 0; i < num_groups; ++i)
  {
    root->createGroup(names[i]);
  }
  create_timer.stop();
  std::cout << "Create time " << create_timer.elapsed() << ".\n";

  std::shuffle(names.begin(), names.end(), g);

  axom::utilities::Timer query_timer(true);
  for(size_t i = 0; i < num_groups; ++i)
  {
    if(!root->hasGroup(names[i]))
    {
      break;
    }
  }
  query_timer.stop();
  std::cout << "Query time " << query_timer.elapsed() << ".\n";

  std::shuffle(names.begin(), names.end(), g);

  axom::utilities::Timer destroy_timer(true);
  for(size_t i = 0; i < num_groups; ++i)
  {
    root->destroyGroup(names[i]);
  }
  destroy_timer.stop();
  std::cout << "Destroy time " << destroy_timer.elapsed() << ".\n";

  axom::utilities::Timer mixed_timer(true);
  if(num_groups > 100)
  {
    mixed_timer.stop();
    int num_shuffles = num_groups / 100;
    for(int j = 0; j < num_shuffles; ++j)
    {
      std::shuffle(names.begin(), names.end(), g);

      mixed_timer.start();
      for(int i = 0; i < 100; ++i)
      {
        if(root->hasGroup(names[i]))
        {
          root->destroyGroup(names[i]);
        }
        else
        {
          root->createGroup(names[i]);
        }
      }
      mixed_timer.stop();
    }
    mixed_timer.start();
  }
  std::cout << "Mixed time " << mixed_timer.elapsed() << ".\n";

  delete ds;

  return 0;
}
