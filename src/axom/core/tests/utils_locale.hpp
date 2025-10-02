// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
//
// file: utils_locale.hpp
// Checks locales and prints out related info
//
//-----------------------------------------------------------------------------

#include "axom/core/utilities/StringUtilities.hpp"
#include "axom/core/utilities/System.hpp"

#include "axom/fmt.hpp"

#include <locale>
#include <string>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <array>
#include <vector>
#include <set>
#include <algorithm>

#include "gtest/gtest.h"

#ifdef _WIN32
  #include <Windows.h>
#endif

namespace
{
#ifdef __linux__
/**
 * \brief Utility function to capture output of system call
 * \note Adapted from https://stackoverflow.com/a/478960
 */
std::string execute_command(const std::string& cmd)
{
  std::array<char, 128> buffer;
  std::string result;

  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);
  if(!pipe)
  {
    throw std::runtime_error(axom::fmt::format("popen() failed for command '{}'", cmd));
  }

  while(fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
  {
    result += buffer.data();
  }

  return result;
}
#endif  // __linux__

#ifdef _WIN32
/// callback function for EnumSystemLocaleEx call in enumerate_locales_windows test
/// note: We're expecting the LPARAM to be a pointer to a vector of wide strings
BOOL CALLBACK MyFuncLocaleEx(LPWSTR pStr, DWORD dwFlags, LPARAM strvec_ptr)
{
  using StrVec = std::vector<std::wstring>;
  reinterpret_cast<StrVec*>(strvec_ptr)->push_back(pStr);
  return TRUE;
}
#endif

}  // end anonymous namespace

//-----------------------------------------------------------------------------
TEST(utils_locale, default_locale)
{
  // This test is informative -- it will print the default and user locales
  std::locale loc(std::locale(), new std::ctype<char>);
  std::cout << "The default locale is " << std::locale().name() << '\n'
            << "The user's locale is " << std::locale("").name() << '\n'
            << "Nameless locale is " << loc.name() << '\n';
}

//-----------------------------------------------------------------------------
TEST(utils_locale, check_en_US_locale)
{
  // This test is informative; not having en_US should not fail the test
  try
  {
    std::locale loc("en_US.UTF-8");
    std::cout << "Name after setting locale to 'en_US.UTF-8': " << loc.name() << std::endl;
  }
  catch(std::runtime_error&)
  {
    std::cout << "Could not set locale to 'en_US.UTF-8'" << std::endl;
  }
}

TEST(utils_locale, try_catch)
{
  std::locale loc;  // initialized to locale::classic()
  try
  {
    loc = std::locale("en_US.UTF-8");
  }
  catch(std::runtime_error&)
  {
    loc = std::locale(loc, "", std::locale::ctype);
  }

  std::cout << "The selected locale after attempting to use 'en_US.UTF-8' "
               "with try/catch is: "
            << loc.name() << std::endl;
}

TEST(utils_locale, axom_locale_default)
{
  try
  {
    const auto loc = axom::utilities::locale();
    std::cout << axom::fmt::format("Axom's default locale is: '{}'\n", loc.name());
  }
  catch(std::runtime_error&)
  {
    FAIL() << "Could not initialize a valid default locale\n";
  }
}

TEST(utils_locale, axom_locale_variations)
{
  const std::int32_t ii = 10000000;
  const std::int64_t ll = 9999999999999ll;
  const double dd = 12345.6789;

  // check several different variations, including some invalid ones
  for(auto& loc_str : {"<DEFAULT>",
                       "",
                       "C",
                       "C.utf8",
                       "en_US",
                       "en_US.utf8",
                       "en_US.UTF8",
                       "en_US.UTF-8",
                       "NON-EXISTENT"})
  {
    try
    {
      const auto loc = (loc_str == std::string("<DEFAULT>")) ? axom::utilities::locale()
                                                             : axom::utilities::locale(loc_str);

      std::cout << axom::fmt::format("Locale for '{}' is '{}'\n", loc_str, loc.name());

      std::cout << axom::fmt::format(loc, "Formatting an int32 using the locale {:L}\n", ii);
      std::cout << axom::fmt::format(loc, "Formatting an int64 using the locale {:L}\n", ll);
      std::cout << axom::fmt::format(loc, "Formatting a double using the locale {:L}\n", dd);
    }
    catch(std::runtime_error&)
    {
      FAIL() << axom::fmt::format("Could not initialize a valid locale for '{}'\n", loc_str);
    }
  }
}

#ifdef __linux__
TEST(utils_locale, enumerate_locales_linux)
{
  std::vector<std::string> locale_list;
  std::vector<std::string> available_locales;

  // get list of locales by running "locale -a" system call
  try
  {
    std::string output = ::execute_command("locale -a");
    locale_list = axom::utilities::string::split(output, '\n');
  }
  catch(const std::runtime_error& e)
  {
    std::cerr << e.what() << '\n';
  }

  // remove non-ascii strings since they're causing problems on blueos
  auto is_non_ascii = [](char c) { return static_cast<unsigned char>(c) > 127; };
  locale_list.erase(std::remove_if(locale_list.begin(),
                                   locale_list.end(),
                                   [=](const std::string& s) {
                                     return std::any_of(s.begin(), s.end(), is_non_ascii);
                                   }),
                    locale_list.end());

  // sort and unique-ify the list
  std::sort(locale_list.begin(), locale_list.end());
  locale_list.erase(std::unique(locale_list.begin(), locale_list.end()), locale_list.end());

  // check if each item can successfully be used to create a locale
  for(const auto& loc_str : locale_list)
  {
    try
    {
      auto loc = std::locale(loc_str);
      available_locales.push_back(loc.name());
    }
    catch(std::runtime_error&)
    {
      // could not instantiate this locale
    }
  }

  std::cout << axom::fmt::format("{} available locales on this system: '{}'\n",
                                 available_locales.size(),
                                 axom::fmt::join(available_locales, "', '"));
}
#endif  // __linux__

#ifdef _WIN32
TEST(utils_locale, enumerate_locales_windows)
{
  std::vector<std::wstring> locale_list;
  std::vector<std::string> available_locales;

  // get vector of locales on this system
  // note: we have to cast the locale_list pointer to an LPARAM to satisfy the Windows API
  LPARAM lparam = reinterpret_cast<LPARAM>(&locale_list);
  EnumSystemLocalesEx(MyFuncLocaleEx, LOCALE_ALL, lparam, nullptr);

  // sort and unique-ify the list
  std::sort(locale_list.begin(), locale_list.end());
  locale_list.erase(std::unique(locale_list.begin(), locale_list.end()), locale_list.end());

  // check if each item can successfully be used to create a locale
  for(const auto& loc_str : locale_list)
  {
    try
    {
      auto loc = std::locale(std::string(loc_str.begin(), loc_str.end()));
      available_locales.push_back(loc.name());
    }
    catch(std::runtime_error&)
    {
      // could not instantiate this locale
    }
  }

  std::cout << axom::fmt::format("{} available locales on this system: '{}'\n",
                                 available_locales.size(),
                                 axom::fmt::join(available_locales, "', '"));
}
#endif  // _WIN32
