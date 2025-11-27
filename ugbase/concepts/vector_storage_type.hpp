#ifndef STORAGE_TYPE_HPP
#define STORAGE_TYPE_HPP

#include <concepts>
#include <ostream>

#include "compatibility.hpp"

template<typename V, typename T = typename V::value_type>
concept VectorStorageType =
    requires(V v, const V cv, std::size_t i, std::size_t n) {

    // required member typedefs
    typename V::value_type;
    typename V::size_type;
    typename V::storage_type;

    // constructors
    { V() } -> std::same_as<V>;
    { V(n) } -> std::same_as<V>;
    { V(cv) } -> std::same_as<V>;

    // capacity
    { cv.size() }      -> std::same_as<std::size_t>;
    //{ cv.capacity() }  -> std::same_as<std::size_t>;
    { v.resize(n) }    -> std::same_as<bool>;
    { cv.reserve(n) }  -> std::same_as<bool>;

    // element access
    { cv.at(i) }       -> std::same_as<const T&>;
    { v.at(i) }        -> std::same_as<T&>;

    { cv[i] }          -> std::same_as<const T&>;
    { v[i] }           -> std::same_as<T&>;

    // operator<< must be valid
    { std::declval<std::ostream&>() << cv }
    -> std::same_as<std::ostream&>;
    };


#endif //STORAGE_TYPE_HPP
