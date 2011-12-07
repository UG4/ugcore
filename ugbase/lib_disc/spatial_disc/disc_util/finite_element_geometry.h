/*
 * finite_element_geometry.h
 *
 *  Created on: 08.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__

#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "common/util/provider.h"

#include <cmath>

namespace ug{

template <	typename TElem,	int TWorldDim,
			typename TTrialSpace, typename TQuadratureRule>
class FEGeometry
{
	public:
	///	type of reference element
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	/// reference element dimension
		static const int dim = ref_elem_type::dim;

	/// world dimension
		static const int worldDim = TWorldDim;

	///	type of trial space
		typedef TTrialSpace trial_space_type;

	///	type of quadrature rule
		typedef TQuadratureRule quad_rule_type;

	///	number of shape functions
		static const size_t nsh = trial_space_type::nsh;

	///	number of integration points
		static const size_t nip = quad_rule_type::nip;

	public:
	///	Constructor
		FEGeometry()
			: m_rQuadRule(Provider<quad_rule_type>::get()),
			  m_rTrialSpace(Provider<trial_space_type>::get())
		{
		//	evaluate local shapes and gradients
			for(size_t ip = 0; ip < nip; ++ip)
				for(size_t sh = 0; sh < nsh; ++sh)
				{
					m_vvShape[ip][sh] = m_rTrialSpace.shape(sh, m_rQuadRule.point(ip));
					m_rTrialSpace.grad(m_vvGradLocal[ip][sh], sh, m_rQuadRule.point(ip));
				}
		}

	/// number of integration points
		size_t num_ip() const {return nip;}

	/// number of shape functions
		size_t num_sh() const {return nsh;}

	/// weight for integration point
		number weight(size_t ip) const {return fabs(m_vDetJ[ip]) * m_rQuadRule.weight(ip);}

	/// local integration point
		const MathVector<dim>& local_ip(size_t ip) const {return m_rQuadRule.point(ip);}

	/// global integration point
		const MathVector<worldDim>& global_ip(size_t ip) const
		{
			UG_ASSERT(ip < nip, "Wrong ip.");
			return m_vIPGlobal[ip];
		}

	/// local integration point
		const MathVector<dim>* local_ips() const {return m_rQuadRule.points();}

	/// global integration point
		const MathVector<worldDim>* global_ips() const{return &m_vIPGlobal[0];}

		/// shape function at ip
		number shape(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < nip, "Wrong index"); UG_ASSERT(sh < nsh, "Wrong index");
			return m_vvShape[ip][sh];
		}

	/// local gradient at ip
		const MathVector<dim>& local_grad(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < nip, "Wrong index"); UG_ASSERT(sh < nsh, "Wrong index");
			return m_vvGradLocal[ip][sh];
		}

	/// global gradient at ip
		const MathVector<worldDim>& global_grad(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < nip, "Wrong index"); UG_ASSERT(sh < nsh, "Wrong index");
			return m_vvGradGlobal[ip][sh];
		}

	/// update Geometry for roid
		bool update_local(ReferenceObjectID roid, LFEID lfeID, size_t orderQuad)
		{
			if(roid != geometry_traits<TElem>::REFERENCE_OBJECT_ID)
			{
				UG_LOG("ERROR in 'FEGeometry::update': Geometry only for "
						<<geometry_traits<TElem>::REFERENCE_OBJECT_ID<<", but "
						<<roid<<" requested.\n");
				return false;
			}
			if(lfeID != m_rTrialSpace.type())
			{
				UG_LOG("ERROR in 'FEGeometry::update': Geometry only for "
						<<m_rTrialSpace.type()<<", but "<<lfeID<<" requested.\n");
				return false;
			}
			if(orderQuad > m_rQuadRule.order())
			{
				UG_LOG("ERROR in 'FEGeometry::update': Geometry only for order "
						<< m_rQuadRule.order()<<", but order "<<orderQuad<<" requested.\n");
				return false;
			}

			return true;
		}

	/// update Geometry for corners
		bool update(TElem* pElem, const MathVector<worldDim>* vCorner)
		{
			return update(pElem, vCorner, m_rTrialSpace.type(), m_rQuadRule.order());
		}

	/// update Geometry for corners
		bool update(TElem* pElem, const MathVector<worldDim>* vCorner,
		            LFEID lfeID, size_t orderQuad)
		{
		//	check
			UG_ASSERT(lfeID == m_rTrialSpace.type(), "Wrong type requested.");
			UG_ASSERT(orderQuad <= m_rQuadRule.order(), "Wrong order requested.");

		//	check if element changed
			if(pElem == m_pElem) return true;
			else m_pElem = pElem;

		//	update the mapping for the new corners
			m_mapping.update(vCorner);

		//	compute global integration points
			for(size_t ip = 0; ip < nip; ++ip)
				m_mapping.local_to_global(m_vIPGlobal[ip], local_ip(ip));

		//	evaluate global data
		//	if reference mapping is linear,
			if(ReferenceMapping<ref_elem_type, worldDim>::isLinear)
			{
			// 	compute transformation inverse and determinate at first ip
				m_mapping.jacobian_transposed_inverse(m_vJTInv[0], local_ip(0));
				m_vDetJ[0] = m_mapping.jacobian_det(local_ip(0));

			//	copy values
				for(size_t ip = 1; ip < nip; ++ip)
				{
					m_vJTInv[ip] = m_vJTInv[0];
					m_vDetJ[ip] = m_vDetJ[0];
				}
			}
		//	else compute jacobian for each point
			{
				for(size_t ip = 0; ip < nip; ++ip)
				{
				// 	compute transformation inverse and determinate at ip
					m_mapping.jacobian_transposed_inverse(m_vJTInv[ip], local_ip(ip));

				//	compute determinant
					m_vDetJ[ip] = m_mapping.jacobian_det(local_ip(ip));
				}
			}

		// 	compute global gradients
			for(size_t ip = 0; ip < nip; ++ip)
				for(size_t sh = 0; sh < nsh; ++sh)
					MatVecMult(m_vvGradGlobal[ip][sh], m_vJTInv[ip], m_vvGradLocal[ip][sh]);

		//	we're done
			return true;
		}

	protected:
	///	current element
		TElem* m_pElem;

	///	Quadrature rule
		const quad_rule_type& m_rQuadRule;

	///	Quadrature rule
		const trial_space_type& m_rTrialSpace;

	///	reference mapping
		ReferenceMapping<ref_elem_type, worldDim> m_mapping;

	///	global integration points
		MathVector<worldDim> m_vIPGlobal[nip];

	///	shape functions evaluated at ip
		number m_vvShape[nip][nsh];

	///	global gradient evaluated at ip
		MathVector<dim> m_vvGradLocal[nip][nsh];

	///	local gradient evaluated at ip
		MathVector<worldDim> m_vvGradGlobal[nip][nsh];

	///	jacobian of transformation at ip
		MathMatrix<worldDim,dim> m_vJTInv[nip];

	///	determinate of transformation at ip
		number m_vDetJ[nip];
};



template <int TWorldDim, int TRefDim = TWorldDim>
class DimFEGeometry
{
	public:
	/// reference element dimension
		static const int dim = TRefDim;

	/// world dimension
		static const int worldDim = TWorldDim;

	public:
	///	default Constructor
		DimFEGeometry() :
			m_roid(ROID_UNKNOWN), m_quadOrder(0), m_lfeID(LFEID(LFEID::NONE, LFEID::INVALID)),
			m_vIPLocal(NULL), m_vQuadWeight(NULL)
		{}

	///	Constructor
		DimFEGeometry(size_t order, LFEID lfeid) :
			m_roid(ROID_UNKNOWN), m_quadOrder(order), m_lfeID(lfeid),
			m_vIPLocal(NULL), m_vQuadWeight(NULL)
		{}

	///	Constructor
		DimFEGeometry(ReferenceObjectID roid, size_t order, LFEID lfeid) :
			m_roid(roid), m_quadOrder(order), m_lfeID(lfeid),
			m_vIPLocal(NULL), m_vQuadWeight(NULL)
		{}

	/// number of integration points
		size_t num_ip() const {return m_nip;}

	/// number of shape functions
		size_t num_sh() const {return m_nsh;}

	/// weight for integration point
		number weight(size_t ip) const
		{
			UG_ASSERT(ip < m_nip, "Wrong index");
			return fabs(m_vDetJ[ip]) * m_vQuadWeight[ip];
		}

	/// local integration point
		const MathVector<dim>& local_ip(size_t ip) const
		{
			UG_ASSERT(ip < m_nip, "Wrong index");
			return m_vIPLocal[ip];
		}

	/// global integration point
		const MathVector<worldDim>& global_ip(size_t ip) const
		{
			UG_ASSERT(ip < m_vIPGlobal.size(), "Wrong ip.");
			return m_vIPGlobal[ip];
		}

	/// local integration point
		const MathVector<dim>* local_ips() const {return m_vIPLocal;}

	/// global integration point
		const MathVector<worldDim>* global_ips() const{return &m_vIPGlobal[0];}

	/// shape function at ip
		number shape(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < m_vvShape.size(), "Wrong index");
			UG_ASSERT(sh < m_vvShape[ip].size(), "Wrong index");
			return m_vvShape[ip][sh];
		}

	/// local gradient at ip
		const MathVector<dim>& local_grad(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < m_vvGradLocal.size(), "Wrong index");
			UG_ASSERT(sh < m_vvGradLocal[ip].size(), "Wrong index");
			return m_vvGradLocal[ip][sh];
		}

	/// global gradient at ip
		const MathVector<worldDim>& global_grad(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < m_vvGradGlobal.size(), "Wrong index");
			UG_ASSERT(sh < m_vvGradGlobal[ip].size(), "Wrong index");
			return m_vvGradGlobal[ip][sh];
		}

	/// update Geometry for roid
		bool update_local(ReferenceObjectID roid, LFEID lfeID, size_t orderQuad)
		{
		//	remember current setting
			m_roid = roid;
			m_lfeID = lfeID;
			m_quadOrder = orderQuad;

		//	request for quadrature rule
			try{
			const QuadratureRule<dim>& quadRule
					= QuadratureRuleProvider<dim>::get_rule(roid, orderQuad);

		//	copy quad informations
			m_nip = quadRule.size();
			m_vIPLocal = quadRule.points();
			m_vQuadWeight = quadRule.weights();

			}catch(UG_ERROR_QuadratureRuleNotRegistered& ex){
				UG_LOG("ERROR in FEGeometry::update: " << ex.get_msg() << ".\n");
				return false;
			}

		//	resize for number of integration points
			m_vIPGlobal.resize(m_nip);
			m_vJTInv.resize(m_nip);
			m_vDetJ.resize(m_nip);

			m_vvGradGlobal.resize(m_nip);
			m_vvGradLocal.resize(m_nip);
			m_vvShape.resize(m_nip);

		//	request for trial space
			try{
			const DimLocalShapeFunctionSet<dim>& lsfs
				 = LocalShapeFunctionSetProvider::get<dim>(roid, m_lfeID);

		//	copy shape infos
			m_nsh = lsfs.num_sh();

		//	resize for number of shape functions
			for(size_t ip = 0; ip < m_nip; ++ip)
			{
				m_vvGradGlobal[ip].resize(m_nsh);
				m_vvGradLocal[ip].resize(m_nsh);
				m_vvShape[ip].resize(m_nsh);
			}

		//	get all shapes by on call
			for(size_t ip = 0; ip < m_nip; ++ip)
			{
				lsfs.shapes(&(m_vvShape[ip][0]), m_vIPLocal[ip]);
				lsfs.grads(&(m_vvGradLocal[ip][0]), m_vIPLocal[ip]);
			}

			}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex){
				UG_LOG("ERROR in FEGeometry::update: " << ex.get_msg() << ".\n");
				return false;
			}
			return true;
		}

	/// update Geometry for corners
		bool update(GeometricObject* pElem, const MathVector<worldDim>* vCorner)
		{
			return update(pElem, vCorner, m_lfeID, m_quadOrder);
		}

	/// update Geometry for corners
		bool update(GeometricObject* pElem, const MathVector<worldDim>* vCorner,
					LFEID lfeID, size_t orderQuad)
		{
		//	check if same element
			if(pElem == m_pElem) return true;
			else m_pElem = pElem;

		//	get reference element type
			ReferenceObjectID roid = (ReferenceObjectID)pElem->reference_object_id();

			////////////////////////
			//	local values
			////////////////////////

		//	if already prepared for this roid, skip update of local values
			if(m_roid != roid || lfeID != m_lfeID || (int)orderQuad != m_quadOrder)
				if(!update_local(roid, lfeID, orderQuad))
					return false;

			////////////////////////
			//	global values
			////////////////////////

		//	get reference element mapping
			try{
			DimReferenceMapping<dim, worldDim>& map
				= ReferenceMappingProvider::get<dim, worldDim>(roid);

		//	update the mapping for the new corners
			map.update(vCorner);

		//	compute global integration points
			map.local_to_global(&(m_vIPGlobal[0]), &(m_vIPLocal[0]), m_nip);

		// 	compute transformation inverse and determinate at ip
			map.jacobian_transposed_inverse(&(m_vJTInv[0]), &(m_vIPLocal[0]), m_nip);

		//	compute determinant
			map.jacobian_det(&(m_vDetJ[0]), &(m_vIPLocal[0]), m_nip);

		// 	compute global gradients
			for(size_t ip = 0; ip < m_nip; ++ip)
				for(size_t sh = 0; sh < m_nsh; ++sh)
					MatVecMult(m_vvGradGlobal[ip][sh], m_vJTInv[ip], m_vvGradLocal[ip][sh]);

			}catch(UG_ERROR_ReferenceMappingMissing& ex){
				UG_LOG("ERROR in FEGeometry::update: " << ex.get_msg() << ".\n");
				return false;
			}

		//	we're done
			return true;
		}

	protected:
	///	current reference object id the local values are prepared for
		ReferenceObjectID m_roid;

	///	current element
		GeometricObject* m_pElem;

	///	current integration order
		int m_quadOrder;

	///	current local finite element id
		LFEID m_lfeID;

	///	number of integration point
		size_t m_nip;

	///	local quadrature points
		const MathVector<dim>* m_vIPLocal;

	///	local quadrature weights
		const number* m_vQuadWeight;

	///	global integration points (size = nip)
		std::vector<MathVector<worldDim> > m_vIPGlobal;

	///	jacobian of transformation at ip (size = nip)
		std::vector<MathMatrix<worldDim,dim> > m_vJTInv;

	///	determinant of transformation at ip (size = nip)
		std::vector<number> m_vDetJ;

	///	number of shape functions
		size_t m_nsh;

	///	shape functions evaluated at ip (size = nip x nsh)
		std::vector<std::vector<number> > m_vvShape;

	///	global gradient evaluated at ip (size = nip x nsh)
		std::vector<std::vector<MathVector<dim> > > m_vvGradLocal;

	///	local gradient evaluated at ip (size = nip x nsh)
		std::vector<std::vector<MathVector<worldDim> > > m_vvGradGlobal;
};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__ */
