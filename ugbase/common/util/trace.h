
#define untested() ( std::cerr <<  "@@#\n@@@:"<< __FILE__ << ":"<< __LINE__ \
          <<":" << __func__ << "\n" )
#define itested()
#define incomplete() ( \
    std::cerr << "@@#\n@@@\nincomplete:" \
              << __FILE__ << ":" << __LINE__ << ":" << __func__ << "\n" )
#define unreachable() ( \
    std::cerr << "@@#\n@@@\nunreachable:" \
              << __FILE__ << ":" << __LINE__ << ":" << __func__ << "\n" )
