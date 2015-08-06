/// \file CompileVersion.hh version info about repository and compile
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015

#ifndef COMPILEVERSION_HH
#define COMPILEVERSION_HH

#include <string>
#include <time.h>

namespace aCORN_MPM {
extern const std::string git_sha;       ///< git repository version SHA hash identifier
extern const time_t compile_time;       ///< time of compile
/// print compile info to stdout
void display_version();
}

#endif
