
#ifdef untested
#undef untested
#endif

#ifdef itested
#undef itested
#endif

// untested: not reached by any test.
#ifdef TRACE_UNTESTED
#define untested() ( std::cerr <<  "@@#\n@@@:"<< __FILE__ << ":"<< __LINE__ \
          <<":" << __func__ << "\n" )
#else
#define untested()
#endif

// interactively tested: no test required.
#ifdef TRACE_ITESTED
#define itested() ( std::cerr <<  "@@#\n@@@:"<< __FILE__ << ":"<< __LINE__ \
          <<":" << __func__ << "\n" )
#else
#define itested()
#endif

#ifndef incomplete
#define incomplete() ( \
    std::cerr << "@@#\n@@@\nincomplete:" \
              << __FILE__ << ":" << __LINE__ << ":" << __func__ << "\n" )
#endif

#ifndef unreachable
#define unreachable() ( \
    std::cerr << "@@#\n@@@\nunreachable:" \
              << __FILE__ << ":" << __LINE__ << ":" << __func__ << "\n" )
#endif
