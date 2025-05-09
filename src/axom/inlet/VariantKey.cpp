// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/VariantKey.hpp"

namespace axom
{
namespace inlet
{
VariantKey& VariantKey::operator=(const int key)
{
  m_int = key;
  m_type = VariantKeyType::Integer;
  return *this;
}

VariantKey& VariantKey::operator=(const std::string& key)
{
  m_string = key;
  m_type = VariantKeyType::String;
  return *this;
}

VariantKey& VariantKey::operator=(std::string&& key)
{
  m_string = std::move(key);
  m_type = VariantKeyType::String;
  return *this;
}

VariantKey& VariantKey::operator=(const char key[])
{
  m_string = key;
  m_type = VariantKeyType::String;
  return *this;
}

VariantKey::operator int() const
{
  if(m_type != VariantKeyType::Integer)
  {
    SLIC_ERROR("[Inlet] Attempted to retrieve an integer from a non-integer key");
  }
  return m_int;
}

VariantKey::operator const std::string&() const
{
  if(m_type != VariantKeyType::String)
  {
    SLIC_ERROR("[Inlet] Attempted to retrieve a string from a non-string key");
  }
  return m_string;
}

InletType VariantKey::type() const
{
  if(m_type == VariantKeyType::Integer)
  {
    return InletType::Integer;
  }
  else if(m_type == VariantKeyType::String)
  {
    return InletType::String;
  }
  SLIC_ERROR("[Inlet] VariantKey tagged union is in invalid state");
  return InletType::Nothing;
}

bool VariantKey::operator==(const VariantKey& other) const
{
  if(m_type != other.m_type)
  {
    return false;
  }
  if(m_type == VariantKeyType::Integer)
  {
    return m_int == other.m_int;
  }
  else if(m_type == VariantKeyType::String)
  {
    return m_string == other.m_string;
  }
  SLIC_ERROR("[Inlet] VariantKey tagged union is in invalid state");
  return false;
}

std::ostream& operator<<(std::ostream& out, const VariantKey& key)
{
  if(key.type() == axom::inlet::InletType::Integer)
  {
    out << static_cast<int>(key);
  }
  else
  {
    out << static_cast<std::string>(key);
  }
  return out;
}

}  // end namespace inlet
}  // end namespace axom
