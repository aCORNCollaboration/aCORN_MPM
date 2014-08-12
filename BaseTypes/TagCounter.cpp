#include "TagCounter.hh"
#include "strutils.hh"

#include <stdlib.h>
#include <cassert>
#include <utility>

template<>
TagCounter<int>::TagCounter(Stringmap m) {
    for(std::multimap< std::string, std::string >::iterator it = m.dat.begin(); it != m.dat.end(); it++)
        add(atoi(it->first.c_str()),atof(it->second.c_str()));
}

template<>
TagCounter<unsigned int>::TagCounter(Stringmap m) {
    for(std::multimap< std::string, std::string >::iterator it = m.dat.begin(); it != m.dat.end(); it++)
        add(atoi(it->first.c_str()),atof(it->second.c_str()));
}

template<>
TagCounter<std::string>::TagCounter(Stringmap m) {
    for(std::multimap< std::string, std::string >::iterator it = m.dat.begin(); it != m.dat.end(); it++)
        add(it->first,atof(it->second.c_str()));
}

template<>
TagCounter< std::pair<unsigned int, unsigned int> >::TagCounter(Stringmap m) {
    for(std::multimap< std::string, std::string >::iterator it = m.dat.begin(); it != m.dat.end(); it++) {
        std::vector<int> v = sToInts(it->first, "/");
        assert(v.size()==2);
        if(v.size()==2)
            add(std::pair<unsigned int, unsigned int>(v[0],v[1]),atof(it->second.c_str()));
    }
}
