/*
 * nedelec.h
 *
 * This file contains implementations of the local shape function set for
 * the so-called Nedelec (or Whitney-1) elements.
 *
 * Created on: 08.02.2013
 * Author: Dmitry Logashenko
 */
#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__NEDELEC__NEDELEC__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__NEDELEC__NEDELEC__

#endif // __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__NEDELEC__NEDELEC__

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Nedelec Set
///////////////////////////////////////////////////////////////////////////////

/**
 * Nedelec (or Whitney-1) base function set for a general element:
 * Not implemented, so this class implements error messages only.
 * For the Nedelec base functions for triangles and tetrahedra cf. the
 * specializations below.
 */
template <typename TRefElement>
class NedelecLSFS
: public
	BaseLocalShapeFunctionSet
		<
			NedelecLSFS<TRefElement>,
			TRefElement::dim, ///< dimensionality of the element
			MathVector<TRefElement::dim>, ///< return type of the shape functions
			MathMatrix<TRefElement::dim, TRefElement::dim> ///< return type of the gradients
		>
{
	public:
	///	Reference Element type
		typedef TRefElement reference_element_type;

	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = reference_element_type::dim;

	private:
	///	Base class
		typedef BaseLocalShapeFunctionSet<NedelecLSFS<reference_element_type>, dim, MathVector<dim>,  MathMatrix<dim,dim> > base_type;

	public:
	///	Domain position type
		typedef typename base_type::position_type position_type;

	///	Shape type
		typedef typename base_type::shape_type shape_type;

	///	Gradient type
		typedef typename base_type::grad_type grad_type;

	protected:
	///	number of shapes (no shapes for the general case)
		static const size_t nsh = 0;

	public:
	///	Constructor
		NedelecLSFS() {};
	
	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline static LFEID type() {return LFEID(LFEID::NEDELEC, 1);}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline static bool continuous() {return false;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			UG_THROW("NedelecLSFS:"
						" Shapes are implemented for triangles and tetrahedra only.");
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline MathVector<dim> shape(const size_t i, const MathVector<dim>& x) const
		{
			UG_THROW("NedelecLSFS:"
						" Shapes are implemented for triangles and tetrahedra only.");
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(MathMatrix<dim,dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			UG_THROW("NedelecLSFS:"
						" Shapes are implemented for triangles and tetrahedra only.");
		}
};

/// Nedelec (or Whitney-1) base function set for triangles
template <>
class NedelecLSFS<ReferenceTriangle>
: public
	BaseLocalShapeFunctionSet
		<
			NedelecLSFS<ReferenceTriangle>,
			ReferenceTriangle::dim, ///< dimensionality of the element
			MathVector<ReferenceTriangle::dim>, ///< return type of the shape functions
			MathMatrix<ReferenceTriangle::dim, ReferenceTriangle::dim> ///< return type of the gradients
		>
{
	public:
	///	Reference Element type
		typedef ReferenceTriangle reference_element_type;

	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = reference_element_type::dim;

	private:
	///	Base class
		typedef BaseLocalShapeFunctionSet<NedelecLSFS<reference_element_type>, dim, MathVector<dim>,  MathMatrix<dim,dim> > base_type;

	public:
	///	Domain position type
		typedef base_type::position_type position_type;

	///	Shape type
		typedef base_type::shape_type shape_type;

	///	Gradient type
		typedef base_type::grad_type grad_type;

	protected:
	///	number of shapes
		static const size_t nsh = reference_element_type::numEdges;

	public:
	///	Constructor
		NedelecLSFS() {};
	
	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline static LFEID type() {return LFEID(LFEID::NEDELEC, 1);}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline static bool continuous() {return false;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{
				case 0:	pos[0] = 0.5;
						pos[1] = 0.0;
						return true;
				case 1:	pos[0] = 0.5;
						pos[1] = 0.5;
						return true;
				case 2:	pos[0] = 0.0;
						pos[1] = 0.5;
						return true;
				default: UG_THROW("NedelecLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline MathVector<dim> shape(const size_t i, const MathVector<dim>& x) const
		{
			MathVector<dim> w1;
			number w0_m, w0_n; // values of the Whitney-0 (linear Lagrange) basis functions
			
			switch (i)
			{
			case 0: // 0->1 (here and below: corners of the edge)
				w0_m = 1.0-x[0]-x[1];
				w0_n = x[0];
				w1[0] = w0_m - w0_n;
				w1[1] =      - w0_n;
				break;
			case 1: // 1->2
				w0_m = x[0];
				w0_n = x[1];
				w1[0] = w0_n;
				w1[1] = w0_m;
				break;
			case 2: // 2->0
				w0_m = x[1];
				w0_n = 1.0-x[0]-x[1];
				w1[0] = - w0_m;
				w1[1] = - w0_m + w0_n;
				break;
			default: UG_THROW("NedelecLSFS: Invalid shape fct index: "<<i);
			}
			
			return w1;
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(MathMatrix<dim,dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			switch(i)
			{
				case 0: // 0->1 (here and below: corners of the edge)
					g(0,0) = -2; g(0,1) = -1;
					g(1,0) = -1; g(1,1) =  0;
					return;
				case 1: // 1->2
					g(0,0) = 0; g(0,1) = 1;
					g(1,0) = 1; g(1,1) = 0;
					return;
				case 2: // 2->0
					g(0,0) =  0; g(0,1) = -1;
					g(1,0) = -1; g(1,1) = -2;
					return;
				default: UG_THROW("NedelecLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}
};

/// Nedelec (or Whitney-1) base function set for tetrahedra
template <>
class NedelecLSFS<ReferenceTetrahedron>
: public
	BaseLocalShapeFunctionSet
		<
			NedelecLSFS<ReferenceTetrahedron>,
			ReferenceTetrahedron::dim, ///< dimensionality of the element
			MathVector<ReferenceTetrahedron::dim>, ///< return type of the shape functions
			MathMatrix<ReferenceTetrahedron::dim, ReferenceTetrahedron::dim> ///< return type of the gradients
		>
{
	public:
	///	Reference Element type
		typedef ReferenceTetrahedron reference_element_type;

	///	Order of Shape functions
		static const size_t order = 1;

	///	Dimension, where shape functions are defined
		static const int dim = reference_element_type::dim;

	private:
	///	Base class
		typedef BaseLocalShapeFunctionSet<NedelecLSFS<reference_element_type>, dim, MathVector<dim>,  MathMatrix<dim,dim> > base_type;

	public:
	///	Domain position type
		typedef base_type::position_type position_type;

	///	Shape type
		typedef base_type::shape_type shape_type;

	///	Gradient type
		typedef base_type::grad_type grad_type;

	protected:
	///	number of shapes
		static const size_t nsh = reference_element_type::numEdges;

	public:
	///	Constructor
		NedelecLSFS() {};
	
	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline static LFEID type() {return LFEID(LFEID::NEDELEC, 1);}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline static bool continuous() {return false;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			switch(i)
			{
				case 0:	pos[0] = 0.5;
						pos[1] = 0.0;
						pos[2] = 0.0;
						return true;
				case 1:	pos[0] = 0.5;
						pos[1] = 0.5;
						pos[2] = 0.0;
						return true;
				case 2:	pos[0] = 0.0;
						pos[1] = 0.5;
						pos[2] = 0.0;
						return true;
				case 3:	pos[0] = 0.0;
						pos[1] = 0.0;
						pos[2] = 0.5;
						return true;
				case 4:	pos[0] = 0.5;
						pos[1] = 0.0;
						pos[2] = 0.5;
						return true;
				case 5:	pos[0] = 0.0;
						pos[1] = 0.5;
						pos[2] = 0.5;
						return true;
				default: UG_THROW("NedelecLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline MathVector<dim> shape(const size_t i, const MathVector<dim>& x) const
		{
			MathVector<dim> w1;
			number w0_m, w0_n; // values of the Whitney-0 (linear Lagrange) basis functions
			
			switch (i)
			{
			case 0: // 0->1 (here and below: corners of the edge)
				w0_m = 1.0-x[0]-x[1]-x[2];
				w0_n = x[0];
				w1[0] = w0_m - w0_n;
				w1[1] =      - w0_n;
				w1[2] =      - w0_n;
				break;
			case 1: // 1->2
				w0_m = x[0];
				w0_n = x[1];
				w1[0] = w0_n;
				w1[1] = w0_m;
				w1[2] = 0.0;
				break;
			case 2: // 2->0
				w0_m = x[1];
				w0_n = 1.0-x[0]-x[1]-x[2];
				w1[0] = - w0_m;
				w1[1] = - w0_m + w0_n;
				w1[2] = - w0_m;
				break;
			case 3: // 0->3
				w0_m = 1.0-x[0]-x[1]-x[2];
				w0_n = x[2];
				w1[0] =      - w0_n;
				w1[1] =      - w0_n;
				w1[2] = w0_m - w0_n;
				break;
			case 4: // 1->3
				w0_m = x[0];
				w0_n = x[2];
				w1[0] = w0_n;
				w1[1] = 0.0;
				w1[2] = w0_m;
				break;
			case 5: // 2->3
				w0_m = x[1];
				w0_n = x[2];
				w1[0] = 0.0;
				w1[1] = w0_n;
				w1[2] = w0_m;
				break;
			default: UG_THROW("NedelecLSFS: Invalid shape fct index: "<<i);
			}
			
			return w1;
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(MathMatrix<dim,dim>& g, const size_t i, const MathVector<dim>& x) const
		{
			switch(i)
			{
				case 0: // 0->1 (here and below: corners of the edge)
					g(0,0) = -2; g(0,1) = -1; g(0,2) = -1;
					g(1,0) = -1; g(1,1) =  0; g(1,2) =  0;
					g(2,0) = -1; g(2,1) =  0; g(2,2) =  0;
					return;
				case 1: // 1->2
					g(0,0) = 0; g(0,1) = 1; g(0,2) = 0;
					g(1,0) = 1; g(1,1) = 0; g(0,2) = 0;
					g(2,0) = 0; g(2,1) = 0; g(2,2) = 0;
					return;
				case 2: // 2->0
					g(0,0) =  0; g(0,1) = -1; g(0,2) =  0;
					g(1,0) = -1; g(1,1) = -2; g(0,2) = -1;
					g(2,0) =  0; g(2,1) = -1; g(2,2) =  0;
					return;
				case 3: // 0->3
					g(0,0) =  0; g(0,1) =  0; g(0,2) = -1;
					g(1,0) =  0; g(1,1) =  0; g(1,2) = -1;
					g(2,0) = -1; g(2,1) = -1; g(2,2) = -2;
					return;
				case 4: // 1->3
					g(0,0) = 0; g(0,1) = 0; g(0,2) = 1;
					g(1,0) = 0; g(1,1) = 0; g(0,2) = 0;
					g(2,0) = 1; g(2,1) = 0; g(2,2) = 0;
					return;
				case 5: // 2->3
					g(0,0) = 0; g(0,1) = 0; g(0,2) = 0;
					g(1,0) = 0; g(1,1) = 0; g(0,2) = 1;
					g(2,0) = 0; g(2,1) = 1; g(2,2) = 0;
					return;
				default: UG_THROW("NedelecLSFS: shape function "<<i<<
									" not found. Only "<<nsh<<" shapes present.");
			}
		}
};

} // namespace ug

/* End of File */
