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
	IPlugInTerurn_ERROR,
	IPlugInReturn_NOT_IMPLEMENTED
};


class IPlugInElementDiscretization{

	public:
		// support assembling on triangles
		template <typename TElem>
		inline IPlugInReturn prepare_element_discretization(TElem* elem)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

		template <typename TElem>
		inline IPlugInReturn prepare_element(TElem* elem)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

		template <typename TElem>
		inline IPlugInReturn assemble_element_JA(TElem* elem, number mat_values[], number u_values[], const uint num_dofs, number time=0.0)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

		template <typename TElem>
		inline IPlugInReturn assemble_element_JM(TElem* elem, number mat_values[], number u_values[], const uint num_dofs, number time=0.0)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

		template <typename TElem>
		inline IPlugInReturn assemble_element_A(TElem* elem, number def_values[], number u_values[], const uint num_dofs, number time=0.0)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

		template <typename TElem>
		inline IPlugInReturn assemble_element_M(TElem* elem, number def_values[], number u_values[], const uint num_dofs, number time=0.0)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

		template <typename TElem>
		inline IPlugInReturn assemble_element_f(TElem* elem, number def_values[], const uint num_dofs, number time=0.0)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

		template <typename TElem>
		inline IPlugInReturn finish_element_discretization(TElem* elem)
		{ return IPlugInReturn_NOT_IMPLEMENTED; };

	public:
		// support assembling on quadrilaterals
		// TODO: other elements
		// [ ... ]
};




} // namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__PLUG_IN_DISC__PLUG_IN_ELEMENT_DISC_INTERFACE__ */
