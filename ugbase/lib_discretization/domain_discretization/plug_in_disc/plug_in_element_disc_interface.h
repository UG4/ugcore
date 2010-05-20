/*
 * plug_in_element_disc_interface.h
 *
 *  Created on: 01.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__PLUG_IN_DISC__PLUG_IN_ELEMENT_DISC_INTERFACE__
#define __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__PLUG_IN_DISC__PLUG_IN_ELEMENT_DISC_INTERFACE__

namespace ug{

enum IPlugInReturn {
	IPlugInReturn_OK = 0,
	IPlugInReturn_ERROR,
	IPlugInReturn_NOT_IMPLEMENTED
};

template <typename TAlgebra>
class IPlugInElementDiscretization{
	public:
		// algebra type
		typedef TAlgebra algebra_type;

		// local matrix type
		typedef typename algebra_type::matrix_type::local_matrix_type local_matrix_type;

		// local vector tyoe
		typedef typename algebra_type::vector_type::local_vector_type local_vector_type;

	public:

		// support assembling on finite elements
		template <typename TElem>
		inline uint num_sh(TElem* elem)
		{ return 0; };

		template <typename TElem>
		inline IPlugInReturn prepare_element_loop(TElem* elem)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

		template <typename TElem>
		inline IPlugInReturn prepare_element(TElem* elem)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

		template <typename TElem>
		inline IPlugInReturn assemble_element_JA(TElem* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

		template <typename TElem>
		inline IPlugInReturn assemble_element_JABnd(TElem* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

		template <typename TElem>
		inline IPlugInReturn assemble_element_JM(TElem* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

		template <typename TElem>
		inline IPlugInReturn assemble_element_A(TElem* elem, local_vector_type& d, const local_vector_type& u, number time=0.0)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

		template <typename TElem>
		inline IPlugInReturn assemble_element_ABnd(TElem* elem, local_vector_type& d, const local_vector_type& u, number time=0.0)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

		template <typename TElem>
		inline IPlugInReturn assemble_element_M(TElem* elem, local_vector_type& d, const local_vector_type& u, number time=0.0)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

		template <typename TElem>
		inline IPlugInReturn assemble_element_f(TElem* elem, local_vector_type& d, number time=0.0)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

		template <typename TElem>
		inline IPlugInReturn finish_element_loop(TElem* elem)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

};




} // namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__PLUG_IN_DISC__PLUG_IN_ELEMENT_DISC_INTERFACE__ */
