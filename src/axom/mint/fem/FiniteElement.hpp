// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_FINITEELEMENT_HPP_
#define MINT_FINITEELEMENT_HPP_

#include "axom/core/Macros.hpp"  // for disable copy/assignment macros
#include "axom/core/Types.hpp"   // for nullptr definition

#include "axom/core/numerics/Matrix.hpp"  // for Matrix definition

#include "axom/mint/config.hpp"       // for axom::IndexType
#include "axom/mint/fem/FEBasis.hpp"  // for FEBasis traits class

namespace axom
{
namespace mint
{
/*!
 * \enum
 *
 * \brief Enum of return codes from FiniteElement::computeReferenceCoords()
 *
 * \see FiniteElement::computeReferenceCoords()
 */
enum
{
  INVERSE_MAP_FAILED = -1, /*!< INVERSE_MAP_FAILED */
  OUTSIDE_ELEMENT,         /*!< OUTSIDE_ELEMENT    */
  INSIDE_ELEMENT           /*!< INSIDE_ELEMENT     */
};

/*!
 * \brief The FiniteElement object is used to represent a mesh element
 *  \f$ \Omega^e \f$ corresponding to a mesh \f$ \mathcal{M} \f$ .
 *
 *  The FiniteElement is associated with a Finite Element Basis (FEBasis),
 *  which defines a set of shape functions \f$ \mathcal{N}_i( \hat{\xi} ) \f$,
 *  on a corresponding element, \f$ \bar{\Omega}^e \f$ in reference
 *  space. This class provides the following key operations:
 *  <ul>
 *   <li> Evaluate the shape function and its derivatives </li>
 *   <li> Forward/Inverse mapping between physical/reference coordinates</li>
 *   <li> Compute the jacobian \f$ \mathcal{J}(\hat{\xi} ) \f$ </li>
 *  </ul>
 *
 *  The FiniteElement class is the primary object that application codes will
 *  use to perform these operations. Once a FiniteElement object is
 *  instantiated it can then be bound to a Finite Element Basis by invoking
 *  the bind_basis( ) method on the target FiniteElement instance, which is
 *  templated on two enum values: (a) the BasisType, defined in
 *  (\ref FEBasisTypes.hpp) and (b) the CellType defined in(\ref CellType.hpp).
 *  The rationale behind this is to insulate the user from the low-level classes
 *  and provide a simple and unified interface to Finite Element operations.
 *
 * A simple example illustrating how to use the FiniteElement class is given
 * below:
 * \code
 *  using axom;
 *
 *  double xp[2];    // test query point
 *  double xr[2];    // query point reference coordinates
 *  double wgts[4];  // buffer to store interpolation weights
 *
 *  ...
 *
 *  mint::Mesh* m = get_mesh();
 *  const int N   = m->getMeshNumberOfCells();
 *
 *  const bool zero_copy = true;
 *
 *  for ( int idx=0; idx < N; ++idx ) {
 *
 *     ...
 *     // gather cell coordinates in a buffer
 *     double coords[ ] = {
 *         x1[ idx ], y1[ idx ],
 *         x2[ idx ], y2[ idx ],
 *         x3[ idx ], y3[ idx ],
 *         x4[ idx ], y4[ idx ],
 *     };
 *
 *     // put cell coordinates in matrix form
 *     numeric::Matrix< double > m( 2, 4, coords, zero_copy );
 *
 *     // construct finite element and bind to basis
 *     mint::FiniteElement fe( m, mint::QUAD, zero_copy );
 *     mint::bind_basis< MINT_LAGRANGE_BASIS, mint::QUAD >( fe );
 *
 *     // compute reference coordinates
 *     int status = fe.computeReferenceCoords( xp, xr );
 *
 *     if ( status == mint::INSIDE_ELEMENT ) {
 *       fe.evaluateShapeFunction( xr, wgts );
 *       ....
 *     } // END if point is inside fe
 *
 *  } // END for all cells
 *
 *  ...
 * \endcode
 *
 * \see FEBasis
 * \see CellType.hpp
 * \see FEBasisTypes.hpp
 * \see ShapeFunction
 */
class FiniteElement
{
public:
  /*!
   * \brief Default constructor. Disabled.
   */
  FiniteElement() = delete;

  /*!
   * \brief Constructs a FiniteElement instance from a matrix consisting of the
   *  element coordinates and a cell type.
   *
   * \param [in] M (ndims x nnodes) matrix consisting of the element coordinates
   * \param [in] cellType the celltype, e.g., MINT_QUAD, etc.
   * \param [in] useExternal optional argument that indicates whether the
   *  coordinate data will be shallow copied to the internal data-structures.
   *  Defaults to false if argument is not specified.
   *
   * \note The element coordinates are arranged in the supplied matrix, M, such
   *  that the number of rows corresponds to the dimension of the element and
   *  number of columns corresponds to the number of nodes.
   */
  FiniteElement(numerics::Matrix<double>& M, CellType cellType, bool useExternal = false);

  /*!
   * \brief Destructor.
   */
  ~FiniteElement();

  /*!
   * \brief Checks if this instance is pointing to an external, user-supplied,
   *  buffer for the element coordinate information.
   *
   * \return status true if an external buffer is used, otherwise, false.
   */
  bool usesExternalBuffer() const { return m_usingExternal; };

  /*!
   * \brief Overrides the max number of iterations for the Newton-Raphson, used
   *  for the isoparametric inverse mapping from physical coordinates to
   *  reference coordinates.
   *
   * \note Each Basis/Element pair prescribes a default value for the maximum
   *  number of Newton-Raphson iterations. This method provides the flexibility
   *  of overriding this value.
   *
   * \warning If bind_basis() is called after invoking this method, the max
   *  newton iterations would be set to the nominal value. Typically, this
   *  method should be called after bind_basis() is invoked.
   *
   * \param [in] N user-supplied number for
   */
  void setMaxSolverIterations(int numIters) { m_maxNewtonIterations = numIters; };

  /*!
   * \brief Returns the max number of iterations used for the Newton-Raphson.
   * \return N max number of iterations for the Newton-Raphson
   */
  int getMaxSolverIterations() const { return m_maxNewtonIterations; };

  /*!
   * \brief Returns the cell type associated with this FiniteElement instance.
   * \return ctype the corresponding cell type.
   * \post (ctype != MINT_UNDEFINED_CELL) && (ctype < MINT_NUM_CELL_TYPES)
   * \see CellType.hpp for a list of possible return values.
   */
  CellType getCellType() const { return m_ctype; };

  /*!
   * \brief Returns the Basis associated with this FiniteElement instance.
   * \return basis the basis type, may be MINT_UNDEFINED_BASIS if no basis
   *  is bound to this FiniteElement.
   *
   * \see FEBasis
   */
  int getBasisType() const { return m_shape_func_type; };

  /*!
   * \brief Returns the physical dimension for this FiniteElement instance.
   * \return dim the dimension of this FiniteElement.
   */
  int getPhysicalDimension() const { return m_dim; };

  /*!
   * \brief Returns the dimension of the element in reference space.
   * \return ref_dim the reference dimension of the element.
   * \post ref_dim \f$ \in [1,3] \f$
   */
  int getReferenceDimension() const { return m_reference_dim; };

  /*!
   * \brief Returns the number of nodes, i.e., degrees-of-freedom (dofs) of the
   *  FiniteElement.
   * \return N the number of nodes of the element
   */
  int getNumNodes() const { return m_numnodes; };

  /*!
   * \brief Returns the number of degrees of freedom associated with the
   *  Finite Element basis bound to this FiniteElement instance.
   * \return N the number of degrees of freedom, or -1 if no basis is set.
   */
  int getNumDofs() const { return m_numdofs; };

  /*!
   * \brief Returns a pointer to the nodes of the element in physical space,
   *  \f$ \forall\hat{x_i} \in \Omega^e \f$
   *
   * \return ptr pointer to the nodes of the elements in physical space.
   *
   * \note The nodes of the element are arranged in a flat array using a
   *  column-major layout. It is convenient to access the nodes of the element
   *  by wrapping the return pointer in a matrix object with getNumNodes()
   *  columns and getPhysicalDimension() rows. This can be achieved as follows:
   *  \code
   *    ...
   *    double* nodesptr = fe->getPhysicalNodes( );
   *    numerics::Matrix< double > pnodes( ndims, nnodes, nodesptr, true );
   *
   *    for ( int inode=0; inode < nnodes; ++inode ) {
   *       double* node = pnodes.getColumn( inode );
   *
   *       std::cout << "Node [" << inode << "]: ";
   *       for ( int idim=0; idim < ndims << ++idim ) {
   *          std::cout << nodes[ idim ] << " ";
   *       } // END
   *
   *       std::cout << std::endl;
   *
   *    } // END for all nodes
   *    ...
   *  \endcode
   *
   * \see numerics::Matrix
   *
   * \post ptr != nullptr
   */
  /// @{

  double* getPhysicalNodes() { return m_xyz; };
  const double* getPhysicalNodes() const { return m_xyz; };

  /// @}

  /*!
   * \brief Returns a pointer to the nodes of the element in reference space,
   *  \f$ \forall\hat{\xi_i} \in \bar{\Omega}^e \f$
   *
   * \return ptr pointer to the nodes of the element in reference space.
   *
   * \note The nodes of the element are arranged in a flat array using a
   *  column-major layout. It is convenient to access the nodes of the element
   *  by wrapping the return pointer in a matrix object as follows:
   *  \code
   *    ...
   *    double* nodesptr = fe->getReferenceNodes( );
   *    numerics::Matrix< double > pnodes( ndims, nnodes, nodesptr, true );
   *
   *    for ( int inode=0; inode < nnodes; ++inode ) {
   *       double* node = pnodes.getColumn( inode );
   *
   *       std::cout << "Node [" << inode << "]: ";
   *       for ( int idim=0; idim < ndims << ++idim ) {
   *          std::cout << nodes[ idim ] << " ";
   *       } // END
   *
   *       std::cout << std::endl;
   *
   *    } // END for all nodes
   *    ...
   *  \endcode
   *
   * \post ptr == nullptr iff getBasisType()==MINT_UNDEFINED_BASIS
   */
  /// @{

  double* getReferenceNodes() { return m_reference_coords; };
  const double* getReferenceNodes() const { return m_reference_coords; };

  /// @}

  /*!
   * \brief Returns a pointer to the centroid of the reference element
   *  \f$ \xi_c \in \bar{\Omega}^e \f$
   *
   * \return ptr pointer to the reference center.
   *
   * \post ptr == nullptr iff getBasisType()==MINT_UNDEFINED_BASIS
   *
   * \note ptr points to a buffer that is ndims long, where ndims is the
   *  dimension of the reference element.
   */
  /// @{

  double* getReferenceCenter() { return m_reference_center; };
  const double* getReferenceCenter() const { return m_reference_center; };

  /// @}

  /*!
   * \brief Given a point in physical space, \f$ \bar{x_p} \f$, this
   *  method computes the corresponding reference coordinates,
   *  \f$ \bar{xi} \f$ in the reference space of the element.
   *
   * \param [in] xp physical coordinates of the point in query.
   * \param [out] xr computed reference coordinates \f$ \bar{\xi} \f$
   * \param [in] TOL optional tolerance for Newton-Raphson. Default is 1.e-12.
   *
   * \return rc return code
   * <ul>
   *  <li> mint::INVERSE_MAP_FAILED if the Newton-Raphson fails </li>
   *  <li> mint::OUTSIDE_ELEMENT if xp is outside the element </li>
   *  <li> mint::INSIDE_ELEMENT if xp is inside the element  </li>
   * </ul>
   *
   * \pre xp != nullptr
   * \pre xr != nullptr
   * \pre this->getBasisType() != MINT_UNDEFINED_BASIS
   */
  int computeReferenceCoords(const double* xp, double* xr, double TOL = 1.e-12);

  /*!
   * \brief Given reference coordinates, \f$ \xi \in \bar{\Omega}^e \f$,
   *  this method computes the corresponding physical coordinates, given by
   *  \f$ \mathbf{x}(\xi)=\sum\limits_{i=0}^n{N}_i^e({\xi}){x_i} \f$
   *
   * \param [in] xr buffer consisting of the reference coordinates
   * \param [out] pt buffer to store the computed physical coordinates
   *
   * \pre xr != nullptr
   * \pre xp != nullptr
   * \pre this->getBasisType() != MINT_UNDEFINED_BASIS
   */
  void computePhysicalCoords(const double* xr, double* xp);

  /*!
   * \brief Given reference coordinates, \f$ \xi \in \bar{\Omega}^e \f$,
   *  this method computes the jacobian, \f$ \mathcal{J}(\xi) \f$
   *
   * \param [in] xr buffer consisting of the reference coordinates
   * \param [out] J the jacobian matrix evaluated at the given location
   *
   * \pre xr != nullptr
   * \pre this->getBasisType() != MINT_UNDEFINED_BASIS
   */
  void jacobian(const double* xr, numerics::Matrix<double>& J);

  /*!
   * \brief Given reference coordinates, \f$ \xi \in \bar{\Omega}^e \f$,
   *  this method evaluates the shape functions at each degree of freedom,
   *  \f$ \phi_i(\xi)=N_i^e(\xi) \f$
   *
   * \param [in] xr reference coordinates where the shape functions is evaluated
   * \param [out] phi buffer to store the shape functions \f$ \phi_i \f$
   *
   * \note The user-supplied output buffer has to be at least ndofs long
   *
   * \pre xr != nullptr
   * \pre phi != nullptr
   */
  void evaluateShapeFunctions(const double* xr, double* phi);

  /*!
   * \brief Given reference coordinates, \f$ \xi \in \bar{\Omega}^e \f$,
   *  this method evaluates the first derivatives of the shape function at
   *  each degree of freedom, \f$ \partial\phi_i(\xi)=N_i^e(\xi) \f$
   *
   * \param [in] xr reference coordinates where derivatives are computed
   * \param [out] phidot buffer to store the shape function derivatives
   *
   * \note The output buffer
   *
   * \note The derivatives in the output buffer, phidot, are arranged using a
   *  column-major flat array layout. It is therefore convenient to think of
   *  the data as a Matrix, where each row represents a degree of freedom and
   *  each column a dimension.
   *
   * \note The matrix class may be used as a way to access the derivative
   *  information more conveniently as illustrated below:
   * \code
   *    ...
   *    fe->evaluateDerivatives( xr, phidot);
   *    numerics::Matrix< double > shape_derivs( ndofs, ndims, phidot, true );
   *
   *    for ( int idim=0; idim < ndims; ++idim ) {
   *
   *       double* derivs = pnodes.getColumn( idim );
   *
   *       std::cout << "dimension => " << idim << ":\n";
   *
   *       for ( int idof=0; idof < ndofs << ++idof ) {
   *
   *          std::cout << "\t dof[" << idof << "]=" << derivs[ idof ] << "\n";
   *
   *       } // END for each dof
   *
   *    } // END for each dimension
   *    ...
   * \endcode
   *
   * \pre xr != nullptr
   * \pre phidot != nullptr
   */
  void evaluateDerivatives(const double* xr, double* phidot);

  /// \name Friend Functions
  /// @{

  /*!
   * \brief Binds the given FiniteElement instance to a FiniteElement basis.
   *
   * \param [in] fe reference to a finite element object.
   *
   * \tparam BasisType the basis type, e.g., MINT_LAGRANGE_BASIS
   * \tparam CellType  the cell type, e.g., MINT_QUAD, etc.
   *
   * \post fe.getBasisType() != MINT_UNDEFINED_BASIS
   */
  template <int BasisType, CellType CellType>
  friend void bind_basis(FiniteElement& fe);

  /// @}

private:
  /// \name Internal Helper methods
  /// @{

  /*!
   * \brief Helper method used to allocate internal data-structures.
   */
  void setUp();

  /*!
   * \brief Helper method used to deallocate all dynamically allocated memory
   *  and destroy internal data-structures.
   */
  void tearDown();

  /*!
   * \brief Given reference coordinates, \f$ \xi \in \bar{\Omega}^e \f$,
   *  this method checks if the point is inside the reference element.
   *
   * \param [in] xr reference coordinates of the point in query.
   * \param [in] TOL optional user-supplied tolerance. Default is 1.e-12.
   *
   * \note Necessary and sufficient condition to determine if
   *  \f$ \xi \in \bar{\Omega}^e \f$ are:
   *  <ul>
   *    <li>
   *     The shape functions at each degree of freedom should be within the
   *     frame of the reference element \f$[ \xi_0, \xi_1 ]\f$:
   *     \f$ \xi_0 \le N_i^e(\xi) \le \xi_1 \f$, \f$ \forall i \f$ </li>
   *    <li>
   *     The shape functions should sum to unity:
   *     \f$ \sum\limits_{i=0}^n{N}_i^e({\xi}){x_i} = 1 \f$
   *    </li>
   *  </ul>
   *
   * \return status true if inside the reference element, otherwise, false.
   *
   * \pre xr != nullptr
   */
  bool inReferenceElement(const double* xr, double TOL = 1.e-12);

  /// @}

  /// \name Private Data Members
  /// @{

  int m_dim;                 /*!< physical dimension */
  CellType m_ctype;          /*!< cell type, e.g., MINT_QUAD, etc. */
  int m_shape_func_type;     /*!< basis type, e.g., MINT_LAGRANGE */
  int m_maxNewtonIterations; /*!< max newton iterations for inverse map */
  int m_numnodes;            /*!< number of nodes of the element */

  double* m_jac;    /*!< stores jacobian at some point */
  double* m_xyz;    /*!< stores element coordinates */
  double* m_phi;    /*!< stores shape functions at some point */
  double* m_phidot; /*!< stores shape function derivatives */

  bool m_usingExternal; /*!< indicates whether the element coordinates
                                     are pointing to an external buffer. */
  /// @}

  /// \name Reference Element Attributes
  /// @{

  using ShapeFunctionPtr = void (*)(const double* lc, double* phi);
  ShapeFunctionPtr m_shapeFunction;            /*! shape function functor */
  ShapeFunctionPtr m_shapeFunctionDerivatives; /*! derivatives functor */

  double m_reference_min; /*!< min of reference frame \f$ \xi_0 \f$ */
  double m_reference_max; /*!< max of reference frame \f$ \xi_1 \f$ */
  int m_reference_dim;    /*!< dimension of the reference space */
  int m_numdofs;          /*!< number of degrees of freedom */

  double* m_reference_coords; /*!< stores reference coords \f$ \xi_i \f$ */
  double* m_reference_center; /*!< stores reference center \f$ \xi_c \f$ */

  /// @}

  DISABLE_COPY_AND_ASSIGNMENT(FiniteElement);
  DISABLE_MOVE_AND_ASSIGNMENT(FiniteElement);
};

}  // namespace mint
}  // namespace axom

//------------------------------------------------------------------------------
//    IMPLEMENTATION OF TEMPLATED FRIEND METHODS
//------------------------------------------------------------------------------
namespace axom
{
namespace mint
{
template <int BasisType, CellType CELLTYPE>
void bind_basis(FiniteElement& fe)
{
  SLIC_ASSERT(CELLTYPE == fe.getCellType());
  SLIC_ASSERT(fe.m_reference_center != nullptr);
  SLIC_ASSERT(fe.m_reference_coords != nullptr);

  using FEBasisType = mint::FEBasis<BasisType, CELLTYPE>;
  using ShapeType = typename FEBasisType::ShapeFunctionType;

  if(CELLTYPE != fe.getCellType())
  {
    SLIC_WARNING("Inconsistent celltype, cannot bind FEBasis!");
    return;
  }

  fe.m_shape_func_type = FEBasisType::BasisType;
  fe.m_maxNewtonIterations = ShapeType::maxNewtonIters();
  fe.m_reference_dim = ShapeType::dimension();
  fe.m_numdofs = ShapeType::numDofs();
  fe.m_reference_min = ShapeType::min();
  fe.m_reference_max = ShapeType::max();

  ShapeType::coords(fe.m_reference_coords);
  ShapeType::center(fe.m_reference_center);

  fe.m_shapeFunction = &ShapeType::evaluate;
  fe.m_shapeFunctionDerivatives = &ShapeType::derivatives;

  SLIC_ASSERT(fe.m_shape_func_type != MINT_UNDEFINED_BASIS);
  SLIC_ASSERT(fe.m_shapeFunction != nullptr);
  SLIC_ASSERT(fe.m_shapeFunctionDerivatives != nullptr);
}

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_FINITEELEMENT_HPP_ */
