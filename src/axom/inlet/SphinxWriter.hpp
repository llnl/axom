// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file SphinxWriter.hpp
 *
 * \brief This file contains the class definition of the SphinxWriter.
 *******************************************************************************
 */

#ifndef INLET_SPHINXWRITER_HPP
#define INLET_SPHINXWRITER_HPP

#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>

#include "axom/sidre.hpp"
#include "axom/inlet/Writer.hpp"

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \class SphinxWriter
 *
 * \brief A Writer that is able write documentation in Sphinx RST format for 
 * a given input file.
 *
 * \see Writer
 *******************************************************************************
 */
class SphinxWriter : public Writer
{
public:
  /*!
  *******************************************************************************
  * \brief A constructor for SphinxWriter.
  * 
  * \param [in] fileName The name of the file the documentation should be written to.
  *******************************************************************************
  */
  SphinxWriter(const std::string& fileName);

  void documentContainer(const Container& container) override;

  void finalize() override;

  virtual ~SphinxWriter() = default;

private:
  /*!
   *****************************************************************************
   * \brief Writes the title in RST syntax.
   *
   * This writes a title to the ostringstream in RST syntax.
   *
   * \param [in] title The title to be written
   *
   *****************************************************************************
   */
  void writeTitle(const std::string& title);

  /*!
   *****************************************************************************
   * \brief Writes the sub-title in RST syntax.
   *
   * This writes a sub-title to the ostringstream in RST syntax.
   *
   * \param [in] sub The sub-title to be written
   *
   *****************************************************************************
   */
  void writeSubtitle(const std::string& sub);

  /*!
   *****************************************************************************
   * \brief Writes a 4 column table in RST syntax.
   *
   * This writes a 4 column table to the ostringstream in RST syntax. The number
   * of rows are determined by the number of Fields found in the vector.
   *
   * \param [in] sub The title of the table written
   * 
   * \param [in] rstTable The 2 dimensional vector containing information to
   * be translated into an RST table
   *
   *****************************************************************************
   */
  void writeTable(const std::string& title, const std::vector<std::vector<std::string>>& rstTable);

  /*!
   *****************************************************************************
   * \brief Writes all tables and their respective titles and descriptions.
   *
   * This parses all of the information from m_rstTables into RST-syntax 
   * documentation and writes it to the ostringstream.
   *
   *****************************************************************************
   */
  void writeAllTables();

  /*!
  *******************************************************************************
  * \struct ContainerData
  *
  * \brief A struct to store data associated with each inlet::Container.
  *
  *******************************************************************************
  */
  struct ContainerData
  {
    /*!
    *******************************************************************************
    * \brief A constructor for the ContainerData struct
    * 
    * This initializes the RST tables's column labels.
    * 
    * \param[in] labels The column labels for the RST table
    *
    *******************************************************************************
    */
    ContainerData(const std::vector<std::string>& fieldLabels,
                  const std::vector<std::string>& functionLabels)
    {
      fieldTable.push_back(fieldLabels);
      functionTable.push_back(functionLabels);
    }

    // Copying shouldn't be needed, these will always be managed in a container
    ContainerData(const ContainerData&) = delete;
    ContainerData(ContainerData&&) = default;

    std::string containerName;
    std::string description;
    bool isSelectedElement = false;
    std::vector<std::vector<std::string>> fieldTable;
    std::vector<std::vector<std::string>> functionTable;
  };

  /*!
  *******************************************************************************
  * \brief Extracts Field information from the given Sidre Group and stores it
  * to be written later.
  * 
  * This extracts information about the Field stored in the given Sidre Group. 
  * This information is stored internally by this class and then written to the
  * document by writeAllTables.
  * 
  * \param [in] sidreGroup The Sidre Group from which Field metadata should be
  * extracted and then stored.
  * \param [inout] currentTable The ContainerData object to write field information to
  *******************************************************************************
  */
  void extractFieldMetadata(const axom::sidre::Group* sidreGroup, ContainerData& currentContainer);

  /*!
  *******************************************************************************
  * \brief Extracts Function information from the given Sidre Group and stores it
  * to be written later.
  * 
  * This extracts information about the Function stored in the given Sidre Group. 
  * This information is stored internally by this class and then written to the
  * document by writeAllTables.
  * 
  * \param [in] sidreGroup The Sidre Group from which Function metadata should be
  * extracted and then stored.
  * \param [inout] currentTable The TableData object to write function information to
  *******************************************************************************
  */
  void extractFunctionMetadata(const axom::sidre::Group* sidreGroup, ContainerData& currentContainer);

  /*!
  *******************************************************************************
  * \brief Gets value information from the given Sidre View and returns
  * it as a string.
  * 
  * \param [in] view The Sidre View containing value information.
  *
  * \return String representation of value information.
  *******************************************************************************
  */
  std::string getValueAsString(const axom::sidre::View* view);

  /*!
  *******************************************************************************
  * \brief Gets range information from the given Sidre View and returns
  * it as a string.
  * 
  * \param [in] view The Sidre View containing range information.
  *
  * \return String representation of range information.
  *******************************************************************************
  */
  std::string getRangeAsString(const axom::sidre::View* view);

  /*!
  *******************************************************************************
  * \brief Gets valid value(s) information from the given Sidre View and returns
  * it as a string.
  * 
  * \param [in] view The Sidre View containing valid value(s) information.
  *
  * \return String representation of valid value(s) information.
  *******************************************************************************
  */
  std::string getValidValuesAsString(const axom::sidre::View* view);

  /*!
  *******************************************************************************
  * \brief Gets valid string value(s) information from the given Sidre Group. 
  * 
  * \param [in] sidreGroup The Sidre Group containing valid string value(s) 
  * information.
  *
  * \return String listing the valid string values.
  *******************************************************************************
  */
  std::string getValidStringValues(const axom::sidre::Group* sidreGroup);

  /*!
  *******************************************************************************
  * \brief Gets function signature information from the given Sidre Group. 
  * 
  * \param [in] sidreGroup The Sidre Group containing function signature information,
  * i.e., one that corresponds to an inlet::Function
  *
  * \return C-style function signature, i.e., Double(Vector, Double)
  *******************************************************************************
  */
  std::string getSignatureAsString(const axom::sidre::Group* sidreGroup);

  std::ofstream m_outFile;
  std::ostringstream m_oss;
  // This is needed to preserve the traversal order of the Inlet::Containers
  std::vector<std::string> m_inletContainerPathNames;
  std::unordered_map<std::string, ContainerData> m_rstTables;
  std::string m_fileName;
  // Used for the RST tables for fields
  std::vector<std::string> m_fieldColLabels;
  // Used for the RST tables for functions
  std::vector<std::string> m_functionColLabels;
};

}  // namespace inlet
}  // namespace axom

#endif
