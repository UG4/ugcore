/*
 * schur_complement_inverse_interface.h
 *
 *  Created on: 08.01.2014
 *      Author: mrupp
 */

#ifndef SCHUR_COMPLEMENT_INVERSE_INTERFACE_H_
#define SCHUR_COMPLEMENT_INVERSE_INTERFACE_H_

#include "common/util/smart_pointer.h"

namespace ug{

template <typename TAlgebra>
class SchurComplementOperator;

template<typename TAlgebra>
class ISchurComplementInverse
{
	typedef typename TAlgebra::vector_type vector_type;
public:
	virtual ~ISchurComplementInverse() {}
	virtual bool init(SmartPtr<SchurComplementOperator<TAlgebra> > op) = 0;
	virtual std::string config_string() const = 0;
	virtual bool apply(vector_type& u, const vector_type& f) = 0;
	virtual bool apply_return_defect(vector_type& u, vector_type& f) = 0;
	virtual bool supports_parallel() const = 0;
};

}
#endif /* SCHUR_COMPLEMENT_INVERSE_INTERFACE_H_ */
