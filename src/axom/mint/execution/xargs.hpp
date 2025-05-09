// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_EXECUTION_ARGS_HPP_
#define MINT_EXECUTION_ARGS_HPP_

/*!
 * \file
 *
 * \brief The xargs::* types are used as template arguments to the execution
 *  model loop/mesh traversal functions. They are used to indicate additional
 *  arguments to the lambda expression that is passed to the loop/mesh traversal
 *  function.
 *
 *  \see interface.hpp
 */

namespace axom
{
namespace mint
{
namespace xargs
{
/*!
 * \brief Indicates that the the lambda expression takes the loop index.
 *
 * \note This is the default argument passed to the lambda expression and
 *  does not need to be specified explicitly.
 */
struct index
{ };

/*!
 * \brief Indicates that the lambda expression takes the corresponding
 *  IJ lattice grid coordinates in addition to the associated linear index.
 *
 * \note This option can be used for node and/or cell loop traversals on
 *  2D Structured meshes.
 */
struct ij
{ };

/*!
 * \brief Indicates that the lambda expression takes the corresponding IJK
 *  lattice grid coordinates in addition to the associated linear index.
 *
 * \note This option can be used for node and/or cell loop traversals on
 *  3D Structured meshes.
 */
struct ijk
{ };

/*!
 * \brief Indicates that the lambda expression takes the corresponding X
 *  coordinate of a mesh node in addition to the associated node index.
 *
 * \note This option can be be use for node mesh traversals with any 1D mesh.
 */
struct x
{ };

/*!
 * \brief Indicates that the lambda expression takes the corresponding X,Y,
 *  coordinates of a mesh node in addition to the associated node index.
 *
 * \note This option can be be used for node mesh traversals with any 2D mesh.
 */
struct xy
{ };

/*!
 * \brief Indicates that the lambda expression takes the corresponding X,Y,Z
 * coordinates of a mesh node in addition to the associated node index.
 *
 * \note This option can be used for node mesh traversals with any 3D mesh.
 */
struct xyz
{ };

/*!
 * \brief Indicates that the lambda expression also takes the node
 *  connectivity information, i.e., the node IDs of the mesh entity being
 *  traversed, e.g., cells or faces and the number of nodes.
 *
 * \note This option can be used for cell/face mesh traversals with any mesh.
 */
struct nodeids
{ };

/*!
 * \brief Indicates that the lambda expression also takes a matrix of the nodal
 *  coordinates and the node connectivity information, i.e., the node IDs of the
 *  mesh entity being traversed, e.g., cells or faces and the number of nodes.
 *
 * \note This option can be used for cell/face mesh traversals with any mesh.
 */
struct coords
{ };

/*!
 * \brief Indicates that the lambda expression also takes the face
 *  connectivity information, i.e., the face IDs of the cell in addition to
 *  the associated cell index. 
 *
 * \note This option can be used for cell mesh traversals with any mesh.
 */
struct faceids
{ };

/*!
 * \brief Indicates that the lambda expression also takes the cell
 *  connectivity information, i.e., the two cell IDs of the face in addition to
 *  the associated face index. 
 *
 * \note This option can be used for face mesh traversals with any mesh.
 */
struct cellids
{ };

} /* namespace xargs */

/*!
 * \brief Traits class used for compile-time checking of xargs types.
 */
template <typename ArgType>
struct xargs_traits
{
  static constexpr bool valid() { return false; };
  static constexpr char* name() { return (char*)""; };
};

template <>
struct xargs_traits<xargs::index>
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::index"); };
};

template <>
struct xargs_traits<xargs::ij>
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::ij"); };
};

template <>
struct xargs_traits<xargs::ijk>
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::ijk"); };
};

template <>
struct xargs_traits<xargs::x>
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::x"); };
};

template <>
struct xargs_traits<xargs::xy>
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::xy"); };
};

template <>
struct xargs_traits<xargs::xyz>
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::xyz"); };
};

template <>
struct xargs_traits<xargs::nodeids>
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::nodeids"); };
};

template <>
struct xargs_traits<xargs::coords>
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::coords"); };
};

template <>
struct xargs_traits<xargs::faceids>
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::faceids"); };
};

template <>
struct xargs_traits<xargs::cellids>
{
  static constexpr bool valid() { return true; };
  static constexpr char* name() { return (char*)("xargs::cellids"); };
};

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_EXEC_ARGS_HPP_ */
