// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file multimat.hpp
 *
 * \brief Contains the MultiMat library header and its template implementation
 *
 */
#ifndef MULTIMAT_H_
#define MULTIMAT_H_

#include "axom/slam.hpp"

#include <vector>
#include <cassert>
#include <stdexcept>
#include <memory>

namespace axom
{
namespace multimat
{
enum class FieldMapping
{
  PER_CELL,
  PER_MAT,
  PER_CELL_MAT
};
enum class DataLayout
{
  CELL_DOM,
  MAT_DOM
};
enum class SparsityLayout
{
  SPARSE,
  DENSE
};
enum class DataTypeSupported
{
  TypeUnknown,
  TypeInt,
  TypeDouble,
  TypeFloat,
  TypeUnsignChar
};

//forward class declarations
template <typename T, typename BiSet>
class MMField2D;
template <typename Field2DType>
class MMSubField2D;
template <typename T, DataLayout D, typename B>
class MMField2DTemplated;

/**
 * \class MultiMat
 *
 * \brief A multimaterial data management class that provides storage in various
 * layouts (dense/sparse, and material-dominant/cell-dominant).
 *
 */
class MultiMat
{
protected:
  // SLAM Set type definitions
  using SetPosType = slam::DefaultPositionType;
  using SetElemType = slam::DefaultPositionType;
  using SetType = slam::Set<SetPosType, SetElemType>;
  using RangeSetType = slam::RangeSet<SetPosType, SetElemType>::ConcreteSet;

public:
  // SLAM Bivariate set type definitions
  using BivariateSetType = slam::BivariateSet<RangeSetType, RangeSetType>;
  using ProductSetType = slam::ProductSet<RangeSetType, RangeSetType>;

private:
  // SLAM Relation typedef
  using IndBufferType = axom::Array<SetPosType>;
  template <typename T>
  using IndViewPolicy = slam::policies::ArrayViewIndirection<SetPosType, T>;

  using VariableCardinality =
    slam::policies::MappedVariableCardinality<SetPosType, IndViewPolicy<SetElemType>>;
  using StaticVariableRelationType =
    slam::StaticRelation<SetPosType, SetElemType, VariableCardinality, IndViewPolicy<SetElemType>, RangeSetType, RangeSetType>;

  using DynamicVariableRelationType = slam::DynamicVariableRelation<RangeSetType, RangeSetType>;

  // SLAM Map type
  using MapStrideType = slam::policies::RuntimeStride<SetPosType>;

  template <typename T>
  using MapType = typename slam::Map<T, RangeSetType, IndViewPolicy<T>, MapStrideType>::ConcreteMap;

  template <typename T, typename BSet = BivariateSetType>
  using BivariateMapType =  //this one has runtime stride
    typename slam::BivariateMap<T, BSet, IndViewPolicy<T>, MapStrideType>::ConcreteMap;

  template <typename T, typename BSet = BivariateSetType>
  using BivariateMapTypeStrideOne =  //this one has compile time stride 1
    typename slam::BivariateMap<T, BSet, IndViewPolicy<T>>::ConcreteMap;

public:
  using SparseRelationType = StaticVariableRelationType;

  // SLAM RelationSet for the set of non-zero cell to mat variables
  using RelationSetType = slam::RelationSet<StaticVariableRelationType>;  //, RangeSetType
  using RelationSetDynType = slam::RelationSet<DynamicVariableRelationType>;

public:
  //1D Field
  template <typename T>
  using Field1D = MapType<T>;

  //2D Field
  //old
  //template <typename T, typename BSet = BivariateSetType>
  //using Field2D = BivariateMapType<T,BSet>;
  template <typename T, typename BSet = BivariateSetType>
  using Field2D = MMField2D<T, BSet>;
  //special
  template <typename T>
  using SparseField2D = MMField2D<T, typename RelationSetType::ConcreteSet>;
  template <typename T>
  using DenseField2D = MMField2D<T, typename ProductSetType::ConcreteSet>;

  template <typename T, DataLayout D, typename B>
  using Field2DTemplated = MMField2DTemplated<T, D, B>;

  template <typename Field2DType>
  using SubField = MMSubField2D<Field2DType>;

  using IndexSet = RangeSetType;                        //For returning set of SparseIndex
  using IdSet = typename BivariateSetType::SubsetType;  //For returning set of DenseIndex

  //Constructors

  /**
   * \brief Constructor to create a new MultiMat object
   *
   * \param data_layout  Select material-dominant or cell-dominant layout to
   *                     store the data in. Default is cell-dominant.
   * \param sparsity_layout Select dense or sparse layout to store the data in.
   *                        Default is sparse.
   */
  MultiMat(DataLayout data_layout = DataLayout::CELL_DOM,
           SparsityLayout sparsity_layout = SparsityLayout::SPARSE);

  /*!
   * \brief Copy constructor for a MultiMat object
   *
   * \param other The MultiMat object to copy from.
   */
  MultiMat(const MultiMat& other);

  /*!
   * \brief Copy assignment operator for a MultiMat object
   *
   * \param other The MultiMat object to copy from.
   */
  MultiMat& operator=(const MultiMat& other)
  {
    if(this != &other)
    {
      MultiMat clone(other);
      *this = std::move(clone);
    }
    return *this;
  }

  /*!
   * \brief Move constructor for MultiMat, defaulted.
   */
  MultiMat(MultiMat&&) = default;

  /*!
   * \brief Move assignment operator for MultiMat, defaulted.
   */
  MultiMat& operator=(MultiMat&&) = default;

  /*!
   * \brief Destructor for MultiMat, defaulted.
   */
  ~MultiMat() = default;

  //Set-up functions
  /**
   * \brief Set the number of materials
   * \pre num_mats > 0
   */
  void setNumberOfMaterials(int num_mats);
  /**
   * \brief Set the number of cells
   * \pre num_cells > 0
   */
  void setNumberOfCells(int num_cells);

  /**
   * \brief Sets the allocator ID to use for allocations. Has the same effect
   *  as calling setSlamAllocatorID() and setFieldAllocatorID() with the same
   *  allocator ID.
   *
   * \param The allocator ID to use
   */
  void setAllocatorID(int alloc_id);

  /**
   * \brief Sets the allocator ID to use for Slam object allocations.
   *  This should point to either a host Umpire pool if only being used on the
   *  host, or a unified memory pool if also being used on the GPU.
   *
   * \param alloc_id the Umpire allocator ID to use
   */
  void setSlamAllocatorID(int alloc_id);

  /**
   * \brief Sets the allocator ID to use for field allocations.
   *
   * \param alloc_id the Umpire allocator ID to use
   */
  void setFieldAllocatorID(int alloc_id);

  /// \brief Returns a pointer to the dense 2d field set
  const ProductSetType* getDense2dFieldSet(DataLayout layout) const
  {
    return &m_denseBivarSet[(int)layout];
  }

  /// \brief Returns a pointer to the sparse 2d field set
  const RelationSetType* getSparse2dFieldSet(DataLayout layout) const
  {
    return &m_sparseBivarSet[(int)layout];
  }

  /**
   * \brief Set the cell-material relation.
   *
   * \detail This function takes the input of boolean vectors and set-up
   * internal data structure. The relation_info vector is indexed by the cell
   * and material index: \n
   *   `idx = mat * num_cell + cell`  \n
   * or `idx = cell * num_mats + mat` \n
   * for material-centric and cell-centric layout, respectively.\n
   * Entry at `idx` should containing 'true' where materials are presented in
   * the cell.\n
   * The number of materials and cell must be set prior to calling this function
   * with setNumberOfMaterials(int) and setNumberOfCells(int)
   *
   * \param A boolean vector of size num_mats * num_cells containing information
   * on if a materials is present in a cell.
   *
   */
  void setCellMatRel(const std::vector<bool>& relation_info, DataLayout layout);

  /*!
   * \brief Set the cell-material relation.
   *
   * \detail This function accepts a compressed sparse row representation of a
   *  relation.
   *
   * \param cardinality The number of secondary elements associated with each
   *  dominant element, i.e.:
   *   * cell-dominant: the number of materials associated with each cell
   *   * material-dominant: the number of cells containing each material
   * \param indices A compressed sparse row array of indices representing the
   *  associated secondary-set elements, i.e.:
   *   * cell-dominant: the index of each material in a cell
   *   * material-dominant: the index of each cell containing the material
   * \param layout The layout of the relation (cell- or material-dominant)
   *
   * \pre The number of materials and cell must be set prior to calling this
   *  function with setNumberOfMaterials(int) and setNumberOfCells(int)
   * \pre If m_slamAllocatorID points to device-accessible memory, cardinality
   *  and indices must be device accessible.
   */
  void setCellMatRel(axom::ArrayView<const SetPosType> cardinality,
                     axom::ArrayView<const SetPosType> indices,
                     DataLayout layout);

  //functions related to fields

  int getNumberOfFields() const { return static_cast<int>(m_fieldNameVec.size()); }

  /**
   * \brief Add a field to the MultiMat object
   *
   *  If a field has the special name 'VolFrac', calling addField will have the
   *  same behavior as calling setVolfracField.
   *  If a field already exists, the call to addExternalField is a no-op, and
   *  the index of the already-existing field will be returned.
   *
   * \tparam T The data type (double, float...) of the field
   * \param field_name The name of the field, used to retrieve the field later
   * \param field_mapping
   * \param data_layout
   * \param sparsity_layout
   * \param data_array The array containing data to the field. The length of the
   *            array should be `num_mats * num_cells * ncomp` if the current
   *            format is Dense, or `num_nonzero * ncomp` if the current format
   *            is Sprase.
   * \param (optional) ncomp The number of component of the field. Default is 1
   * \return int the index of the field, can be used to retrieve the field later
   */
  template <typename T>
  int addField(const std::string& field_name,
               FieldMapping field_mapping,
               DataLayout data_layout,
               SparsityLayout sparsity_layout,
               axom::ArrayView<T> data_array,
               int ncomp = 1);

  /**
   * \brief Add an externally-managed field to the MultiMat object
   *
   *  If a field has the special name 'VolFrac', calling addField will have the
   *  same behavior as calling setVolfracField. A volume fraction field will
   *  always be internally-managed.
   *  If a field already exists, the call to addExternalField is a no-op, and
   *  the index of the already-existing field will be returned.
   *
   * \tparam T The data type (double, float...) of the field
   * \param field_name The name of the field, used to retrieve the field later
   * \param field_mapping
   * \param data_layout
   * \param sparsity_layout
   * \param data_array The array containing data to the field. The length of the
   *            array should be `num_mats * num_cells * ncomp` if the current
   *            format is Dense, or `num_nonzero * ncomp` if the current format
   *            is Sprase.
   * \param (optional) ncomp The number of component of the field. Default is 1
   * \return int the index of the field, can be used to retrieve the field later
   */
  template <typename T>
  int addExternalField(const std::string& field_name,
                       FieldMapping field_mapping,
                       DataLayout data_layout,
                       SparsityLayout sparsity_layout,
                       axom::ArrayView<T> data_array,
                       int ncomp = 1);

  /**
   * \brief Delete a field from the MultiMat object.
   *
   * \param name The name of the field to remove.
   */
  void removeField(const std::string& name);

private:
  template <typename T>
  int addFieldArray_impl(const std::string&,
                         FieldMapping,
                         DataLayout,
                         SparsityLayout,
                         axom::ArrayView<T>,
                         bool,
                         int);

public:
  /**
   * \brief Set the volume fraction field
   * \detail volume fraction field is assumed to be a double. Its field index is
   * always 0, and the name of the field is "Volfrac"
   *
   * \param data_array the array containing the volume fraction information
   * \return int the volume fraction field index, which is always zero.
   */
  template <typename T>
  int setVolfracField(axom::ArrayView<T> data_array, DataLayout layout, SparsityLayout sparsity);

  /// \overload
  int setVolfracField(axom::ArrayView<const double> data_array,
                      DataLayout layout,
                      SparsityLayout sparsity);

  /**
   * \brief Search for and return the field index of the field given its name.
   *
   * \param field_name the name of the field
   * \return if found, the index of the field; otherwise, -1
   */
  int getFieldIdx(const std::string& field_name) const;

  /**
   * \brief Return the name of the field given its index.
   *
   * \param field_idx the index of the field
   * \return if found, the name of the field; otherwise, the empty string
   */
  std::string getFieldName(int field_idx) const;

  /**
   * \brief Search for and return the field given the field name.
   * \detail the field is of type Field1D, containing an entry for each cell
   * or material. To retrieve a field of type Field2D, use get2dField().
   *
   * \tparam T The data type of the field
   * \param field_name the name of the field
   * \return Field1D<T>& the field reference
   */
  template <typename T>
  Field1D<T> get1dField(const std::string& field_name);

  /// \overload
  template <typename T>
  Field1D<const T> get1dField(const std::string& field_name) const;

  /**
   * \brief Search for and return the field given the field name.
   * \detail the field is of type Field2D, containing an entry for each cell and
   * each material. To retrieve a field of type Field1D, use get1dField().
   *
   * \tparam T The data type of the field
   * \param field_name the name of the field
   * \return Field1D<T>& the field reference
   */
  template <typename T>
  Field2D<T> get2dField(const std::string& field_name);

  /// \overload
  template <typename T>
  Field2D<const T> get2dField(const std::string& field_name) const;

  template <typename T, typename BSetType>
  Field2D<T, BSetType> get2dField(const std::string& field_name);

  template <typename T>
  DenseField2D<T> getDense2dField(const std::string& field_name);

  template <typename T>
  SparseField2D<T> getSparse2dField(const std::string& field_name);

  template <typename T, DataLayout D, typename B>
  Field2DTemplated<T, D, B> getTemplated2DField(const std::string& field_name);

  template <typename T, typename BSetType>
  slam::BivariateMap<T, BSetType, IndViewPolicy<T>> get2dFieldAsSlamBivarMap(
    const std::string& field_name);

  /**
   * \brief Get the volume fraction field
   */
  Field2D<double> getVolfracField();
  Field2D<const double> getVolfracField() const;

  /**
   * \brief Get a set of index for a Subfield.
   * \detail for a cell-dominant layout, this is equivalent to getting the set
   * of materials presented in a cell. Vice versa, for a material-dominant
   * layout, this returns a set of cells containing a material.\n
   * getMatInCell() and getCellContainingMat() are layout specific
   * calls for this function.
   *
   * \param idx the index of the subfield
   * \return IdSet A set of index in the subfield
   */
  IdSet getSubfieldIndex(int idx);

  /** Cell-dominant version of getSubfieldIndex() **/
  IdSet getMatInCell(int cell_id);
  /** Material-dominant version of getSubfieldIndex() **/
  IdSet getCellContainingMat(int mat_id);

  /**
   * \brief Get a set to index a Subfield in a Field2D.
   * \detail This returns a set of indices that can be used to directly access
   * a subfield in a Field2D without having to retrieve a Subfield first.\n
   * Accessing the Field2D using this indexing set should be done with the
   * bracket operator.
   *
   * For example:
   * \code
   *   IndexSet indexing_set = getSubfieldIndexingSet( idx );
   *   for( int i : indexing_set );
   *       value = field2d[ i ];
   * \endcode
   *
   * Or to account for multiple components:
   * \code
   *   IndexSet indexing_set = getSubfieldIndexingSet( idx );
   *   for( int i : indexing_set );
   *       for( int c = 0; c < numComp; c++)
   *           value = field2d[ i * numComp + c ];
   * \endcode
   *
   * getIndexingSetOfCell() and getIndexingSetOfMat() are layout specific
   * calls for this function.
   *
   * \param idx The index of the subfield
   * \return IndexSet A Set of index for the subfield
   */
  IndexSet getSubfieldIndexingSet(int idx, DataLayout layout, SparsityLayout sparsity) const;

  /** Cell-dominant version of getSubfieldIndexingSet() **/
  IndexSet getIndexingSetOfCell(int cell_id, SparsityLayout sparsity) const;
  /** Material-dominant version of getSubfieldIndexingSet() **/
  IndexSet getIndexingSetOfMat(int mat_id, SparsityLayout sparsity) const;

  /** Return the number of material this object holds **/
  int getNumberOfMaterials() const { return m_nmats; };
  /** Return the number of cells this object holds **/
  int getNumberOfCells() const { return m_ncells; }

  //Layout modification functions

  void convertFieldLayout(int field_idx, SparsityLayout, DataLayout);

  /*!
   * \brief Converts a field with a given index to a sparse layout.
   *  No-op if the field is already in sparse layout, or if the field is not a
   *  2D cell-material field.
   *
   * \param field_idx the index of the field to convert
   *
   * \post getFieldSparsityLayout(field_idx) == SparsityLayout::SPARSE
   */
  void convertFieldToSparse(int field_idx);

  /*!
   * \brief Converts a field with a given index to a dense layout.
   *  No-op if the field is already in dense layout, or if the field is not a
   *  2D cell-material field.
   *
   * \param field_idx the index of the field to convert
   *
   * \post getFieldSparsityLayout(field_idx) == SparsityLayout::DENSE
   */
  void convertFieldToDense(int field_idx);
  SparsityLayout getFieldSparsityLayout(int field_idx) const;

  void transposeField(int field_idx);
  void convertFieldToMatDom(int field_idx);
  void convertFieldToCellDom(int field_idx);
  DataLayout getFieldDataLayout(int field_idx) const;

  std::string getFieldMappingAsString(int field_idx) const;
  std::string getFieldDataLayoutAsString(int field_idx) const;
  std::string getFieldSparsityLayoutAsString(int field_idx) const;
  std::string getFieldDataTypeAsString(int field_idx) const;

  /** Convert the data to be stored in the specified layout. **/
  void convertLayout(DataLayout, SparsityLayout);
  /** Convert the data to be stored in sparse/compact layout **/
  void convertLayoutToSparse();
  /** Convert the data to be stored in dense layout **/
  void convertLayoutToDense();
  /** Convert the data to be stored in cell-dominant layout. **/
  void convertLayoutToCellDominant();
  /** Convert the data to be stored in material-dominant layout. **/
  void convertLayoutToMaterialDominant();

  /**
   * \brief Get the FieldMapping for a field.
   *
   * A FieldMapping of a field describes if there is an entry for each cell, for
   * each material, or for each cell x material.
   *
   * \param field_idx the index of the field
   */
  FieldMapping getFieldMapping(int field_idx) const { return m_fieldMappingVec[field_idx]; }

  /**
   * \brief Get the DataType for a field.
   *
   * A DataTypeSupported for a field indicates the underlying supported data
   * type (e.g. double).
   *
   * \param field_idx the index of the field
   */
  DataTypeSupported getFieldDataType(int field_idx) const { return m_dataTypeVec[field_idx]; }

  //Dynamic mode functions

  /** Convert the layout to use dynamic layout. In dynamic layout, users can
   *  call addEntry() or removeEntry() to change the materials in a cell.
   */
  void convertToDynamic();
  /** Convert the layout to use static layout. **/
  void convertToStatic();

  /**
   * \brief Add a material to a cell.
   *
   * In a cell-dominant layout, firstIdx is the cell index, secondIdx is the
   * material index.  Conversely, for a material-dominant layout, firstIdx
   * is the material index and secondIdx is the cell index.
   *
   * \return true if the material is added to the cell, false if it already exists.
   */
  bool addEntry(int firstIdx, int secondIdx);
  /**
   * \brief Remove a material from a cell
   *
   * In a cell-dominant layout, firstIdx is the cell index, secondIdx is the
   * material index. Conversely, for a material-dominant layout, firstIdx
   * is the material index and secondIdx is the cell index.
   *
   * \return true if a material is removed. False if it was not in the cell.
   */
  bool removeEntry(int firstIdx, int secondIdx);

  /** Print the detail of this object */
  void print() const;
  /**
   * \brief  Return true if the object is valid, false otherwise.
   *
   * \param verboseOutput (Optional) if true, this function will print out its
   * valid/invalid detail.
   */
  bool isValid(bool verboseOutput = false) const;

protected:
  //Return the Set pointer associated with the given FieldMapping or field idx
  const RangeSetType* getMappedRangeSet(FieldMapping mapping) const
  {
    if(mapping == FieldMapping::PER_CELL)
    {
      return &getCellSet();
    }
    else if(mapping == FieldMapping::PER_MAT)
    {
      return &getMatSet();
    }
    else
    {
      SLIC_ASSERT("Cannot map Cell-Material relation to a RangeSet.");
      return nullptr;
    }
  }

  //Return the BivariateSet used for the specified layouts or the field
  const BivariateSetType* get_mapped_biSet(DataLayout, SparsityLayout) const;
  const BivariateSetType* get_mapped_biSet(int field_idx) const;
  //Return the relation for the specified field idx or the mapping
  StaticVariableRelationType* getRel(int field_idx);

private:  //private functions
  //Given a relation (cell->mat or mat->cell), create the other relation
  void makeOtherRelation(DataLayout layout);

  //helper functions
  template <typename DataType>
  void convertToSparse_helper(int map_i);
  template <typename DataType>
  void convertToDense_helper(int map_i);

  template <typename DataType>
  void transposeField_helper(int field_idx);

  /*!
   * \brief Returns the associated cell set.
   */
  RangeSetType& getCellSet() { return m_sets[0]; }
  const RangeSetType& getCellSet() const { return m_sets[0]; }

  /*!
   * \brief Returns the associated material set.
   */
  RangeSetType& getMatSet() { return m_sets[1]; }
  const RangeSetType& getMatSet() const { return m_sets[1]; }

  /*!
   * \brief Returns a reference to the corresponding array of offsets for a
   *        static relation corresponding to a layout.
   *
   * \param layout The layout type of the relation (cell- or mat-dominant)
   */
  IndBufferType& relBeginVec(DataLayout layout);

  /*!
   * \brief Returns a reference to the corresponding array of indices for
   *        a static relation corresponding to a layout.
   *
   * \param layout The layout type of the relation (cell- or mat-dominant)
   */
  IndBufferType& relIndVec(DataLayout layout);

  /*!
   * \brief Returns a reference to the corresponding array of indices for
   *        a static relation corresponding to a layout.
   *
   * \param layout The layout type of the relation (cell- or mat-dominant)
   */
  IndBufferType& relFirstIndVec(DataLayout layout);

  /*!
   * \brief Returns a reference to the static relation corresponding to a
   *        layout.
   *
   * \param layout The layout type of the relation (cell- or mat-dominant)
   */
  StaticVariableRelationType& relStatic(DataLayout layout);

  /*!
   * \brief Returns a reference to the dynamic relation corresponding to a
   *        layout.
   *
   * \param layout The layout type of the relation (cell- or mat-dominant)
   */
  DynamicVariableRelationType& relDynamic(DataLayout layout);

  /*!
   * \brief Returns a reference to the dominant set of elements in the relation
   *        corresponding to a layout.
   *
   * \param layout The layout type of the relation (cell- or mat-dominant)
   */
  RangeSetType& relDominantSet(DataLayout layout);

  /*!
   * \brief Returns a reference to the secondary set of elements in the relation
   *        corresponding to a layout.
   *
   * \param layout The layout type of the relation (cell- or mat-dominant)
   */
  RangeSetType& relSecondarySet(DataLayout layout);

  /*!
   * \brief Returns a reference to a sparse set corresponding to a relation.
   *
   * \param layout The layout type of the relation (cell- or mat-dominant)
   */
  RelationSetType& relSparseSet(DataLayout layout);
  const RelationSetType& relSparseSet(DataLayout layout) const;

  /*!
   * \brief Returns a reference to a product set corresponding to a relation.
   *
   * \param layout The layout type of the relation (cell- or mat-dominant)
   */
  ProductSetType& relDenseSet(DataLayout layout);
  const ProductSetType& relDenseSet(DataLayout layout) const;

  /*!
   * \brief Returns true if the static relation corresponding to the given data
   *        layout is valid.
   * \param layout The layout type of the relation (cell- or mat-dominant)
   */
  bool hasValidStaticRelation(DataLayout layout) const;

  /*!
   * \brief Returns true if the dynamic relation corresponding to the given data
   *        layout is valid.
   * \param layout The layout type of the relation (cell- or mat-dominant)
   */
  bool hasValidDynamicRelation(DataLayout layout) const;

  template <typename T>
  Field1D<T> get1dFieldImpl(int fieldIdx) const;
  template <typename T>
  Field2D<T> get2dFieldImpl(int fieldIdx) const;

  template <typename BSet>
  BSet getCompatibleBivarSet(int fieldIdx) const;

  template <typename BSet>
  struct ConvertibleTraits;

  template <typename SetType1, typename SetType2, typename InterfaceType>
  struct ConvertibleTraits<slam::ProductSet<SetType1, SetType2, InterfaceType>>
  {
    using Type = ProductSetType;
    constexpr static SparsityLayout Layout = SparsityLayout::DENSE;
  };

  template <typename Relation, typename SetType1, typename SetType2, typename InterfaceType>
  struct ConvertibleTraits<slam::RelationSet<Relation, SetType1, SetType2, InterfaceType>>
  {
    using Type = RelationSetType;
    constexpr static SparsityLayout Layout = SparsityLayout::SPARSE;
  };

private:
  int m_slamAllocatorId;
  int m_fieldAllocatorId;
  unsigned int m_ncells, m_nmats;

  //slam set variables
  axom::Array<RangeSetType> m_sets;
  //slam relation variables (for sparse-layout fields)
  //Depending on the layout of the field, each field can be mapped to different Relations
  //Relation can be nullptr if no field is using said relation
  //cell to mat relation data
  IndBufferType m_cellMatRel_beginsVec;
  IndBufferType m_cellMatRel_indicesVec;
  IndBufferType m_cellMatRel_firstIndicesVec;
  //mat to cell relation data
  IndBufferType m_matCellRel_beginsVec;
  IndBufferType m_matCellRel_indicesVec;
  IndBufferType m_matCellRel_firstIndicesVec;
  //relation objects stored in unified memory
  axom::Array<StaticVariableRelationType> m_staticRelations;
  axom::Array<DynamicVariableRelationType> m_dynamicRelations;

  // sparse layout bivariate sets
  axom::Array<RelationSetType> m_sparseBivarSet;
  // dense layout bivariate sets
  axom::Array<ProductSetType> m_denseBivarSet;

  // Transposition maps for sparse conversions of data layout
  // These map flat indices between cell-dominant and material-dominant layouts
  IndBufferType m_flatCellToMatIndexMap;
  IndBufferType m_flatMatToCellIndexMap;

  struct FieldBacking
  {
  private:
    bool m_isOwned {false};

    axom::Array<std::uint8_t> m_ucharData;
    axom::Array<std::int32_t> m_intData;
    axom::Array<float> m_floatData;
    axom::Array<double> m_dblData;

    axom::ArrayView<std::uint8_t> m_ucharView;
    axom::ArrayView<std::int32_t> m_intView;
    axom::ArrayView<float> m_floatView;
    axom::ArrayView<double> m_dblView;

    template <typename T>
    void setArrayView(axom::ArrayView<T> new_array);

  public:
    template <typename T>
    axom::Array<T>& getArray();

    template <typename T>
    axom::ArrayView<T> getArrayView();

    FieldBacking() = default;

    template <typename T>
    FieldBacking(axom::ArrayView<T> input_array, bool owned, int allocatorID)
    {
      if(owned)
      {
        m_isOwned = true;
        getArray<T>() = axom::Array<T>(input_array, allocatorID);
      }
      else
      {
        m_isOwned = false;
        setArrayView<T>(input_array);
      }
    }

    template <typename T>
    FieldBacking(axom::ArrayView<const T> input_array, bool owned, int allocatorID)
    {
      SLIC_ASSERT(owned == true);
      AXOM_UNUSED_VAR(owned);

      m_isOwned = true;
      getArray<T>() = axom::Array<T>(input_array, allocatorID);
    }

    void moveSpaces(int new_alloc_id)
    {
      if(m_isOwned)
      {
        m_ucharData = axom::Array<unsigned char>(m_ucharData, new_alloc_id);
        m_intData = axom::Array<int>(m_intData, new_alloc_id);
        m_floatData = axom::Array<float>(m_floatData, new_alloc_id);
        m_dblData = axom::Array<double>(m_dblData, new_alloc_id);
      }
      else
      {
        SLIC_ERROR("Cannot move unowned array to a different allocator ID.");
      }
    }

    bool isOwned() const { return m_isOwned; }
  };

  //std::vector of information for each fields
  std::vector<std::string> m_fieldNameVec;
  std::vector<FieldMapping> m_fieldMappingVec;
  std::vector<std::unique_ptr<FieldBacking>> m_fieldBackingVec;
  std::vector<DataTypeSupported> m_dataTypeVec;
  std::vector<DataLayout> m_fieldDataLayoutVec;
  std::vector<SparsityLayout> m_fieldSparsityLayoutVec;
  std::vector<int> m_fieldStrideVec;

  //To store the static layout information when converting to dynamic.
  struct Layout
  {
    DataLayout data_layout;
    SparsityLayout sparsity_layout;
  };

  //Layout used during the static mode
  std::vector<Layout> m_layout_when_static;  //per field
  Layout m_static_layout;

  bool m_dynamic_mode;  //true if in dynamic mode

  template <typename T, typename S>
  friend class MMField2D;  // every type of MMField2D is a friend
  template <typename Field2D>
  friend class MMSubField2D;  // every type of MMSubField2D is a friend
  template <typename T, DataLayout D, typename S>
  friend class MMField2DTemplated;
};  //end MultiMat class

//--------------- MultiMat template function definitions -----------------//

template <class T>
int MultiMat::addField(const std::string& arr_name,
                       FieldMapping arr_mapping,
                       DataLayout data_layout,
                       SparsityLayout sparsity_layout,
                       axom::ArrayView<T> data_arr,
                       int stride)
{
  SLIC_ASSERT(stride > 0);

  //make sure the name does not conflict
  int fieldIdx = getFieldIdx(arr_name);
  if(fieldIdx == 0)
  {  //this is the vol frac array. call setVolfrac instead
    SLIC_ASSERT(arr_mapping == FieldMapping::PER_CELL_MAT);
    SLIC_ASSERT(stride == 1);
    setVolfracField(data_arr, data_layout, sparsity_layout);
    return 0;
  }
  else if(fieldIdx > 0)
  {
    // There is already an array with the current name. And it's not Volfrac.
    // Don't add the new field.
    SLIC_WARNING("Multimat: field with name \"" << arr_name << "\" already exists. Skipping.");

    return fieldIdx;
  }
  else
  {
    //No field with this name. Proceed with adding the field
    return addFieldArray_impl<>(arr_name, arr_mapping, data_layout, sparsity_layout, data_arr, true, stride);
  }
}

template <class T>
int MultiMat::addExternalField(const std::string& arr_name,
                               FieldMapping arr_mapping,
                               DataLayout data_layout,
                               SparsityLayout sparsity_layout,
                               axom::ArrayView<T> data_arr,
                               int stride)
{
  SLIC_ASSERT(stride > 0);

  //make sure the name does not conflict
  int fieldIdx = getFieldIdx(arr_name);
  if(fieldIdx == 0)
  {  //this is the vol frac array. call setVolfrac instead
    SLIC_ASSERT(arr_mapping == FieldMapping::PER_CELL_MAT);
    SLIC_ASSERT(stride == 1);
    SLIC_WARNING(
      "Multimat: attempting to add volume fraction field as an external "
      "field. This will be ignored and the underlying data will be copied.");
    setVolfracField(data_arr, data_layout, sparsity_layout);
    return 0;
  }
  else if(fieldIdx > 0)
  {
    // There is already an array with the current name. And it's not Volfrac.
    // Don't add the new field.
    SLIC_WARNING("Multimat: field with name \"" << arr_name << "\" already exists. Skipping.");

    return fieldIdx;
  }
  else
  {
    //No field with this name. Proceed with adding the field
    return addFieldArray_impl<>(arr_name,
                                arr_mapping,
                                data_layout,
                                sparsity_layout,
                                data_arr,
                                false,
                                stride);
  }
}

template <typename T>
int MultiMat::setVolfracField(axom::ArrayView<T> data_array, DataLayout layout, SparsityLayout sparsity)
{
  axom::Array<double> volfrac_dbls(data_array.size());
  for(int i = 0; i < data_array.size(); i++)
  {
    volfrac_dbls[i] = data_array[i];
  }
  axom::ArrayView<const double> volfrac_view = volfrac_dbls;
  return setVolfracField(volfrac_view, layout, sparsity);
}

template <typename T>
int MultiMat::addFieldArray_impl(const std::string& field_name,
                                 FieldMapping field_mapping,
                                 DataLayout data_layout,
                                 SparsityLayout sparsity_layout,
                                 axom::ArrayView<T> data_arr,
                                 bool owned,
                                 int stride)
{
  unsigned int new_arr_idx = m_fieldNameVec.size();

  m_fieldNameVec.push_back(field_name);
  m_fieldMappingVec.push_back(field_mapping);
  m_fieldBackingVec.emplace_back(new FieldBacking());
  m_fieldDataLayoutVec.push_back(data_layout);
  m_fieldSparsityLayoutVec.push_back(sparsity_layout);
  m_fieldStrideVec.push_back(stride);

  if(std::is_same<T, int>::value)
  {
    m_dataTypeVec.push_back(DataTypeSupported::TypeInt);
  }
  else if(std::is_same<T, double>::value)
  {
    m_dataTypeVec.push_back(DataTypeSupported::TypeDouble);
  }
  else if(std::is_same<T, float>::value)
  {
    m_dataTypeVec.push_back(DataTypeSupported::TypeFloat);
  }
  else if(std::is_same<T, unsigned char>::value)
  {
    m_dataTypeVec.push_back(DataTypeSupported::TypeUnsignChar);
  }
  else
  {
    m_dataTypeVec.push_back(DataTypeSupported::TypeUnknown);
  }

  SLIC_ASSERT(m_fieldNameVec.size() == m_dataTypeVec.size());
  SLIC_ASSERT(m_fieldNameVec.size() == m_fieldMappingVec.size());
  SLIC_ASSERT(m_fieldNameVec.size() == m_fieldDataLayoutVec.size());
  SLIC_ASSERT(m_fieldNameVec.size() == m_fieldSparsityLayoutVec.size());
  SLIC_ASSERT(m_fieldNameVec.size() == m_fieldStrideVec.size());

#ifdef AXOM_DEBUG
  if(field_mapping == FieldMapping::PER_CELL_MAT)
  {
    const BivariateSetType* s = get_mapped_biSet(data_layout, sparsity_layout);
    SLIC_ASSERT(s != nullptr);
    SLIC_ASSERT(s->size() * stride == data_arr.size());
  }
  else
  {
    SLIC_ASSERT(field_mapping == FieldMapping::PER_CELL || field_mapping == FieldMapping::PER_MAT);
    const RangeSetType& s = *getMappedRangeSet(field_mapping);
    SLIC_ASSERT(s.size() * stride == data_arr.size());
  }
#endif

  m_fieldBackingVec.back() = std::make_unique<FieldBacking>(data_arr, owned, m_fieldAllocatorId);

  return new_arr_idx;
}

template <typename T>
MultiMat::Field1D<T> MultiMat::get1dField(const std::string& field_name)
{
  int fieldIdx = getFieldIdx(field_name);

  if(fieldIdx < 0)
  {
    SLIC_ERROR("Multimat: No field with the name \"" + field_name + "\" was found.");
  }

  return get1dFieldImpl<T>(fieldIdx);
}

template <typename T>
MultiMat::Field1D<const T> MultiMat::get1dField(const std::string& field_name) const
{
  int fieldIdx = getFieldIdx(field_name);

  if(fieldIdx < 0)
  {
    SLIC_ERROR("Multimat: No field with the name \"" + field_name + "\" was found.");
  }

  return get1dFieldImpl<const T>(fieldIdx);
}

template <typename T>
MultiMat::Field2D<T> MultiMat::get2dField(const std::string& field_name)
{
  int fieldIdx = getFieldIdx(field_name);

  if(fieldIdx < 0)
  {
    SLIC_ERROR("Multimat: No field with the name \"" + field_name + "\" was found.");
  }

  return get2dFieldImpl<T>(fieldIdx);
}

template <typename T>
MultiMat::Field2D<const T> MultiMat::get2dField(const std::string& field_name) const
{
  int fieldIdx = getFieldIdx(field_name);

  if(fieldIdx < 0)
  {
    throw std::invalid_argument("No field with this name is found");
  }

  return get2dFieldImpl<const T>(fieldIdx);
}

template <typename T>
MultiMat::Field1D<T> MultiMat::get1dFieldImpl(int fieldIdx) const
{
  if(m_fieldMappingVec[fieldIdx] == FieldMapping::PER_CELL ||
     m_fieldMappingVec[fieldIdx] == FieldMapping::PER_MAT)
  {
    return Field1D<T>(*getMappedRangeSet(m_fieldMappingVec[fieldIdx]),
                      m_fieldBackingVec[fieldIdx]->getArrayView<std::remove_const_t<T>>(),
                      m_fieldStrideVec[fieldIdx]);
  }
  else
  {
    SLIC_ASSERT(m_fieldMappingVec[fieldIdx] == FieldMapping::PER_CELL_MAT);

    //Right now we're allowing Field2D (BivariateMap) to be returned as
    // a Field1D (Map) so it can be accessed like a 1d array, but the
    // indexing information would be lost.
    RangeSetType bisetFlat(get_mapped_biSet(fieldIdx)->size());
    return Field1D<T>(bisetFlat,
                      m_fieldBackingVec[fieldIdx]->getArrayView<std::remove_const_t<T>>(),
                      m_fieldStrideVec[fieldIdx]);
  }
}

template <typename T>
MultiMat::Field2D<T> MultiMat::get2dFieldImpl(int fieldIdx) const
{
  SLIC_ASSERT(m_fieldMappingVec[fieldIdx] == FieldMapping::PER_CELL_MAT);

  return Field2D<T>(*this,
                    get_mapped_biSet(fieldIdx),
                    fieldIdx,
                    m_fieldBackingVec[fieldIdx]->getArrayView<std::remove_const_t<T>>(),
                    m_fieldStrideVec[fieldIdx]);
}

template <typename T, typename BSetType>
MultiMat::Field2D<T, BSetType> MultiMat::get2dField(const std::string& field_name)
{
  // Get a reference to the unspecialized BMap
  auto bmap = get2dField<T>(field_name);

  //create instance of that map
  int fi = getFieldIdx(field_name);

  if(fi < 0)
  {
    throw std::invalid_argument("No field with this name is found");
  }

  BSetType bsetValue = getCompatibleBivarSet<BSetType>(fi);

  Field2D<T, BSetType> typedBMap(*this, bsetValue, fi, bmap.getMap()->data(), bmap.stride());

  return typedBMap;
}

template <typename T>
MultiMat::DenseField2D<T> MultiMat::getDense2dField(const std::string& field_name)
{
  // Get a reference to the unspecialized BMap
  auto bmap = get2dField<T>(field_name);

  //create instance of that map
  int fieldIdx = getFieldIdx(field_name);

  if(fieldIdx < 0)
  {
    throw std::invalid_argument("No field with this name is found");
  }

  SLIC_CHECK_MSG(bmap.isDense(),
                 "Attempting to get sparse field \""
                   << field_name << "\" as a "
                   << "dense field. Convert the field to a sparse field first.");

  typename ProductSetType::ConcreteSet prod_set = relDenseSet(m_fieldDataLayoutVec[fieldIdx]);

  DenseField2D<T> typedBMap(*this, prod_set, fieldIdx, bmap.getMap()->data(), bmap.stride());

  return typedBMap;
}

template <typename T>
MultiMat::SparseField2D<T> MultiMat::getSparse2dField(const std::string& field_name)
{
  // Get a reference to the unspecialized BMap
  auto bmap = get2dField<T>(field_name);

  //create instance of that map
  int fieldIdx = getFieldIdx(field_name);

  if(fieldIdx < 0)
  {
    throw std::invalid_argument("No field with this name is found");
  }

  SLIC_CHECK_MSG(bmap.isSparse(),
                 "Attempting to get dense field \""
                   << field_name << "\" as a "
                   << "sparse field. Convert the field to a dense field first.");

  typename RelationSetType::ConcreteSet rel_set = relSparseSet(m_fieldDataLayoutVec[fieldIdx]);

  SparseField2D<T> typedBMap(*this, rel_set, fieldIdx, bmap.getMap()->data(), bmap.stride());

  return typedBMap;
}

template <typename T, DataLayout D, typename B>
MultiMat::Field2DTemplated<T, D, B> MultiMat::getTemplated2DField(const std::string& field_name)
{
  // Get a reference to the unspecialized BMap
  auto bmap = get2dField<T>(field_name);

  int fieldIdx = getFieldIdx(field_name);

  Field2DTemplated<T, D, B> typedBMap(*this, fieldIdx, bmap.getMap()->data(), bmap.stride());

  return typedBMap;
}

// Warning: The return type uses a compile time stride of one!
template <typename T, typename BSetType>
slam::BivariateMap<T, BSetType, MultiMat::IndViewPolicy<T>> MultiMat::get2dFieldAsSlamBivarMap(
  const std::string& field_name)
{
  // Get a reference to the unspecialized BMap
  auto bmap = get2dField<T>(field_name);

  int fieldIdx = getFieldIdx(field_name);

  BSetType bsetValue = getCompatibleBivarSet<BSetType>(fieldIdx);

  // Create instance of templated BivariateMap
  slam::BivariateMap<T, BSetType, IndViewPolicy<T>> typedBMap(bsetValue,
                                                              bmap.getMap()->data(),
                                                              bmap.stride());

  return typedBMap;
}

template <typename BSet>
BSet MultiMat::getCompatibleBivarSet(int fieldIdx) const
{
  using FromBSetType = typename ConvertibleTraits<BSet>::Type;

  SLIC_CHECK_MSG(m_fieldSparsityLayoutVec[fieldIdx] == ConvertibleTraits<BSet>::Layout,
                 "BivariateSet type is incompatible with stored field sparsity layout.");

  const FromBSetType* ptr = static_cast<const FromBSetType*>(get_mapped_biSet(fieldIdx));

  return *ptr;
}

// Output operators.
std::ostream& operator<<(std::ostream& os, FieldMapping mapping);
std::ostream& operator<<(std::ostream& os, DataLayout layout);
std::ostream& operator<<(std::ostream& os, SparsityLayout sparsity);
std::ostream& operator<<(std::ostream& os, DataTypeSupported type);

}  //end namespace multimat
}  //end namespace axom

#include "axom/multimat/mmfield.hpp"

#endif
