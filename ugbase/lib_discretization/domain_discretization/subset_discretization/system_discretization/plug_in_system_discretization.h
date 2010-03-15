/*
 * plug_in_system_discretization.h
 *
 *  Created on: 03.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__SUBSET_DISCRETIZATION__SYSTEM_DISCRETIZATION__PLUG_IN_SYSTEM_DISCRETIZATION__
#define __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__SUBSET_DISCRETIZATION__SYSTEM_DISCRETIZATION__PLUG_IN_SYSTEM_DISCRETIZATION__

namespace ug{

template < typename TDomain, typename TAlgebra, template <typename TDomain, typename TElem> class TElemDisc >
class PlugInSystemDiscretization : public SystemDiscretizationInterface<TDomain, TAlgebra>{
	public:
		virtual bool prepare_element(GeometricObject* elem){
			Trinagle* pTriangle;
			Quadrilateral* pQuadrilateral;

			if((pTriangle = dynamic_cast<Triangle*>(elem)) != NULL)
			{
				return m_elemDiscTriangle.prepare_element(pTriangle);
			}
			else if((pQuadrilateral = dynamic_cast<Quadrilateral*>(elem)) != NULL)
			{
				return m_elemDiscQuadrilateral.prepare_element(pQuadrilateral);
			}
			else
			{
				std::cout << "Type not found" << std::endl;
				return false;
			}
		}


	protected:
		TElemDisc<TDomain, Triangle> m_elemDiscTriangle;
};



} // namespace ug


#endif /* __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__SUBSET_DISCRETIZATION__SYSTEM_DISCRETIZATION__PLUG_IN_SYSTEM_DISCRETIZATION__ */
