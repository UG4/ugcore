
#ifndef DEMANGLE_H_
#define DEMANGLE_H_
#include <string>

namespace ug{
/**
 * demangles C++ function names like _ZZ12ug_backtracev = ug_backtrace().
 * also demangles them when a lot of them appear "in between". make sure they are
 * seperated by ' ', '\n' or '\t' and start with _, like they do in backtrace_symbols.
 * Works only in POSIX, otherwise returns str
 * @sa demangle
 * @param str mangled strings, e.g. _ZZ12ug_backtracev
 * @return the demangled string e.g. ug_backtrace()
 */
std::string demangle_block(const char *str);

/**
 * demangles C++ function and class names
 * Works only in POSIX, otherwise returns str
 * @param str mangled string, containing stuff e.g. 3barI5emptyLi17EE
 * @return the demangled string e.g. bar<empty, 17>
 * */
std::string demangle(const char *str);
}
#endif /* DEMANGLE_H_ */
