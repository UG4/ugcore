/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef IVECTOR_H_
#define IVECTOR_H_

#include "util.h"
#include "te.h"
#include "ivector.h"
#include "vec_functions.h"

namespace ug{
///////////////////////////////////////////////////////////////////
//							Vector
///////////////////////////////////////////////////////////////////


//!
class IVector : public TE_Vector<IVector>
{
public:

// INTERFACE
	//! clones the vector (deep-copy) including values
	virtual SmartPtr<IVector> clone() const = 0;
	//! clones the vector (deep-copy) excluding values
	virtual SmartPtr<IVector> clone_without_values(size_t size) const = 0;

	virtual size_t capacity() const = 0;
	virtual size_t size() const = 0;
	virtual void resize(size_t newSize, bool bCopyValues=true) = 0;
	virtual void clear() = 0;


	// algebra
	// use IVector_AlgebraDownCastTo to cast all these functions to VecAdd etc.
	virtual void vec_add(double a1, double a2, const IVector &iv2) = 0;
	virtual void vec_add(double a1, double a2, const IVector &iv2, double a3, const IVector &iv3) = 0;
	virtual void vec_add(double a1, double a2, const IVector &iv2, double a3, const IVector &iv3, double a4, const IVector &iv4) = 0;
	virtual double vec_prod(const IVector &v) const = 0;
	virtual double norm2() const = 0;

	virtual void set(double d) = 0;

	virtual void print(std::ostream &output) const = 0;
	virtual std::string short_desc() const =0;


// HELPER
	template<typename T>
	T &downcast()
	{
		T *t =dynamic_cast<T*>(this);
		UG_ASSERT(t != nullptr, "could not downcast " << TypeName(this) << " to " << TypeName<T>());
		return *t;
	}
	template<typename T>
	const T &downcast() const
	{
		const T *t =dynamic_cast<const T*>(this);
		UG_ASSERT(t != nullptr, "could not downcast " << TypeName(this) << " to " << TypeName<T>());
		return *t;
	}

// IMPLEMENTATIONS
	double norm() const { return sqrt(norm2()); }

	void operator=(const IVector &v)
	{
		vec_add(0.0, 1.0, v);
	}
	void operator+=(const IVector &v)
	{
		vec_add(1.0, 1.0, v);
	}
	void operator-=(const IVector &v)
	{
		vec_add(1.0, -1.0, v);
	}

	void operator*=(double alpha)
	{
		vec_add(0.0, alpha, *this);
	}

// TEMPLATE EXPRESSIONS IMPLEMENTATIONS
	void operator=(const TE_VecScale<IVector> &t)
	{
		vec_add(0.0, t.scaling(), t.vec());
	}
	void operator=(const TE_VecAdd2<IVector> &t)
	{
		vec_add(0.0, t.a1, t.v1, t.a2, t.v2);
	}
	void operator=(const TE_VecAdd3<IVector> &t)
	{
		vec_add(0.0, t.a1, t.v1, t.a2, t.v2, t.a3, t.v3);
	}
	void operator+=(const TE_VecAdd2<IVector> &t)
	{
		vec_add(1.0, t.a1, t.v1, t.a2, t.v2);
	}
	void operator+=(const TE_VecAdd3<IVector> &t)
	{
		vec_add(1.0, t.a1, t.v1, t.a2, t.v2, t.a3, t.v3);
	}
	void operator-=(const TE_VecAdd2<IVector> &t)
	{
		vec_add(1.0, -t.a1, t.v1, -t.a2, t.v2);
	}
	void operator-=(const TE_VecAdd3<IVector> &t)
	{
		vec_add(1.0, -t.a1, t.v1, -t.a2, t.v2, -t.a3, t.v3);
	}
};

template<typename TVec>
class IVector_AlgebraDownCastTo : public IVector
{
	using this_type = IVector_AlgebraDownCastTo<TVec> ;
	TVec &downcast()
	{
		return *((TVec*) this);
	}
	const TVec &downcast() const
	{
		return *((TVec*) this);
	}
public:
	virtual void vec_add(double a1, double a2, const IVector &iv2)
	{
		VecAdd(a1, downcast(), a2, iv2.downcast<TVec>());
	}
	virtual void vec_add(double a1, double a2, const IVector &iv2, double a3, const IVector &iv3)
	{
		VecAdd(a1, downcast(), a2, iv2.downcast<TVec>(), a3, iv3.downcast<TVec>());
	}
	virtual void vec_add(double a1, double a2, const IVector &iv2, double a3, const IVector &iv3, double a4, const IVector &iv4)
	{
		VecAdd(a1, downcast(), a2, iv2.downcast<TVec>(), a3, iv3.downcast<TVec>(), a4, iv4.downcast<TVec>());
	}
	virtual double vec_prod(const IVector &iv) const
	{
		return VecProd(downcast(), iv.downcast<TVec>());
	}

	virtual double norm2() const
	{
		return VecNorm2(downcast());
	}

	virtual void set(double d)
	{
		TVec &t = downcast();
		for(size_t i=0; i<t.size(); i++)
			t[i] = d;
	}

	virtual void print(std::ostream &output) const
	{
		const TVec &t = downcast();
		output << short_desc() << "\n";
		for(size_t i=0; i<t.size(); i++)
			output << i << " = " << t[i] << "\n";
	}

	virtual std::string short_desc() const
	{
		std::stringstream ss;
		ss << TypeName(downcast()) << "/IVector" <<  "[" << size() << "]";
		return ss.str();
	}

	friend std::ostream &operator<<(std::ostream &output, const this_type &v)
	{
		v.print(output);
		return output;
	}

};

} // namespace ug



#endif /* IVECTOR_H_ */
