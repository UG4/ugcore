#ifndef __H__UG_for_each_in_vec
#define __H__UG_for_each_in_vec


///	Allows iteration over all members of an std::vector compatible type
/**	Use e.g. like this:
 * \code
 * std::vector<T>	vec;
 * //...
 * for_each_in_vec(T& t, vec){
 *    t.memberFunc();
 * }end_for
 * \endcode
 * The '{' and '}' brackets are hereby optional, since a new block is automatically
 * defined between 'for_each_in_vec' and 'end_for'.
 *
 * The specified vector has to feature methods 'T& operator[](size_t i)' and
 * 'size_t size()'.
 * \{ */
#define for_each_in_vec(_vfeDecl, _vfeVec) \
			for(size_t _vfeI = 0; _vfeI < _vfeVec.size(); ++_vfeI){\
				_vfeDecl = _vfeVec[_vfeI];

#define end_for	}
/** \} */

#endif	//__H__UG_for_each_in_vec
