// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <cstdlib>  // for free
#include <sstream>  // for std::ostringstream

#ifdef WIN32
  #define NOMINMAX
  #include <Windows.h>
  #include <WinBase.h>
  #include <DbgHelp.h>
#else
  #include <execinfo.h>  // for backtrace()
  #include <ciso646>
  #if !defined(_LIBCPP_VERSION)
    #include <cxxabi.h>  // for abi::__cxa_demangle
  #endif
#endif

constexpr int MAX_FRAMES = 25;

namespace axom
{
namespace slic
{
namespace internal
{
#ifdef WIN32

std::string stacktrace()
{
  void* stack[MAX_FRAMES];
  std::ostringstream oss;

  unsigned short frames;
  SYMBOL_INFO* symbol;
  HANDLE process;

  process = GetCurrentProcess();

  SymInitialize(process, NULL, TRUE);

  frames = CaptureStackBackTrace(0, MAX_FRAMES, stack, NULL);
  symbol = (SYMBOL_INFO*)calloc(sizeof(SYMBOL_INFO) + 256 * sizeof(char), 1);
  symbol->MaxNameLen = 255;
  symbol->SizeOfStruct = sizeof(SYMBOL_INFO);

  oss << "\n** StackTrace of " << frames << " frames **\n";
  for(int i = 0; i < frames; i++)
  {
    char outString[512];
    SymFromAddr(process, (DWORD64)(stack[i]), 0, symbol);

    sprintf_s(outString, "%i: %s - 0x%0X", frames - i - 1, symbol->Name, symbol->Address);
    oss << outString << std::endl;
  }

  free(symbol);
  oss << "=====\n\n";

  return (oss.str());
}

#else /* #ifdef WIN32 */

//------------------------------------------------------------------------------
std::string demangle(char* backtraceString, int frame)
{
  char* mangledName = nullptr;
  char* functionOffset = nullptr;
  char* returnOffset = nullptr;

  #ifdef __APPLE__
  /* On apple machines the mangled function name always starts at the 58th
   * character */
  constexpr int APPLE_OFFSET = 58;
  mangledName = backtraceString + APPLE_OFFSET;
  for(char* p = backtraceString; *p; ++p)
  {
    if(*p == '+')
    {
      functionOffset = p;
    }
    returnOffset = p;
  }
  #else
  for(char* p = backtraceString; *p; ++p)
  {
    if(*p == '(')
    {
      mangledName = p;
    }
    else if(*p == '+')
    {
      functionOffset = p;
    }
    else if(*p == ')')
    {
      returnOffset = p;
      break;
    }
  }
  #endif

  std::ostringstream oss;

  // if the line could be processed, attempt to demangle the symbol
  if(mangledName && functionOffset && returnOffset && mangledName < functionOffset)
  {
    *mangledName = 0;
    mangledName++;
  #ifdef __APPLE__
    *(functionOffset - 1) = 0;
  #endif
    *functionOffset = 0;
    ++functionOffset;
    *returnOffset = 0;
    ++returnOffset;

    int status = false;
  #if !defined(_LIBCPP_VERSION)
    char* realName = abi::__cxa_demangle(mangledName, nullptr, nullptr, &status);
  #else
    char* realName = mangledName;
  #endif

    // if demangling is successful, output the demangled function name
    if(status == 0)
    {
      oss << "Frame " << frame << ": " << realName << std::endl;
    }
    // otherwise, output the mangled function name
    else
    {
      oss << "Frame " << frame << ": " << mangledName << std::endl;
    }

  #if !defined(_LIBCPP_VERSION)
    free(realName);
  #endif
  }
  // otherwise, print the whole line
  else
  {
    oss << "Frame " << frame << ": " << backtraceString << std::endl;
  }

  return (oss.str());
}

std::string stacktrace()
{
  void* array[MAX_FRAMES];

  const int size = backtrace(array, MAX_FRAMES);
  char** strings = backtrace_symbols(array, size);

  // skip first stack frame (points here)
  std::ostringstream oss;
  oss << "\n** StackTrace of " << size - 1 << " frames **\n";
  for(int i = 1; i < size && strings != nullptr; ++i)
  {
    oss << internal::demangle(strings[i], i);
  }

  oss << "=====\n\n";

  free(strings);

  return (oss.str());
}

#endif /* #ifdef WIN32 */

} /* namespace internal */

} /* namespace slic */

} /* namespace axom */
