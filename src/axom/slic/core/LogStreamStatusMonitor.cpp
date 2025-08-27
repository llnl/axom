// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file LogStreamStatusMonitor.cpp
 *
 * \brief Implementation of the LogStreamStatusMonitor class.
 *
 ******************************************************************************
 */

#include "LogStreamStatusMonitor.hpp"

namespace axom
{
namespace slic
{

//------------------------------------------------------------------------------
LogStreamStatusMonitor::LogStreamStatusMonitor()
    : m_streamVec()
{
}

//------------------------------------------------------------------------------
LogStreamStatusMonitor::~LogStreamStatusMonitor()
{
    finalize();
}

//------------------------------------------------------------------------------
void LogStreamStatusMonitor::addStream(LogStream* ls)
{
    m_streamVec.push_back(ls);

}

//------------------------------------------------------------------------------
bool LogStreamStatusMonitor::hasPendingMessages() const
{
    int has_pending_messages = 0;
    for (auto& stream : m_streamVec) {
        has_pending_messages += static_cast<int>(stream->hasPendingMessages());
    }

    return has_pending_messages > 0;
}

//------------------------------------------------------------------------------
void LogStreamStatusMonitor::finalize()
{
}

}
}