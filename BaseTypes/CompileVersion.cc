/// \file CompileVersion.cc
// This file was produced under the employ of the United States Government,
// and is consequently in the PUBLIC DOMAIN, free from all provisions of
// US Copyright Law (per USC Title 17, Section 105).
//
// -- Michael P. Mendenhall, 2015
#include "CompileVersion.hh"
#include <stdio.h>

namespace aCORN_MPM {
    
const time_t compile_time = time(NULL);

#ifdef GIT_SHA
#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)
const std::string git_sha = STRINGIFY(GIT_SHA);
#else
const std::string git_sha = "unknown";
#endif

void display_version() {
    printf("aCORN_MPM repository %s, compiled %li\n", git_sha.c_str(), compile_time);
}

}
