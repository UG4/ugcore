
#ifndef COMPATIBILITY_HPP
#define COMPATIBILITY_HPP


#if defined(__cpp_concepts) && __cpp_concepts >= 201907L
    #define USE_CONCEPTS 1
#else
    #define USE_CONCEPTS 0
#endif


#if USE_CONCEPTS
    #define CONCEPT(T) T
#else
    #define CONCEPT(T) typename
#endif



#endif //COMPATIBILITY_HPP
