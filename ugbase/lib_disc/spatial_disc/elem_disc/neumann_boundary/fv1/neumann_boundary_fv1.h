
#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY___NEUMANN_BOUNDARY_FV1__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY___NEUMANN_BOUNDARY_FV1__

// other ug4 modules
#include "common/common.h"

// library intern headers
#include "../neumann_boundary_base.h"

namespace ug{

template<typename TDomain>
class NeumannBoundaryFV1
	: public NeumannBoundaryBase<TDomain>
{
	private:
	///	Base class type
		typedef NeumannBoundaryBase<TDomain> base_type;

	///	Base class type
		typedef NeumannBoundaryFV1<TDomain> this_type;

	/// error estimator type
		typedef SideAndElemErrEstData<TDomain> err_est_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	default constructor
		NeumannBoundaryFV1(const char* function);

	///	add a boundary value
	///	\{
		void add(SmartPtr<CplUserData<number, dim> > data, 			const char* BndSubsets, const char* InnerSubsets);
		void add(SmartPtr<CplUserData<number, dim, bool> > user, 		const char* BndSubsets, const char* InnerSubsets);
		void add(SmartPtr<CplUserData<MathVector<dim>, dim> > user, 	const char* BndSubsets, const char* InnerSubsets);
	/// \}

	protected:
		using typename base_type::Data;

	///	Unconditional scalar user data
		struct NumberData : public base_type::Data
		{
			NumberData(SmartPtr<CplUserData<number, dim> > data,
			           std::string BndSubsets, std::string InnerSubsets)
				: base_type::Data(BndSubsets, InnerSubsets)
			{
				import.set_data(data);
			}

			template<typename TElem, typename TFVGeom>
			void extract_bip(const TFVGeom& geo);

			template <typename TElem, typename TFVGeom>
			void lin_def(const LocalVector& u,
			             std::vector<std::vector<number> > vvvLinDef[],
			             const size_t nip);

			template <int refDim>
			void set_local_ips(const MathVector<refDim>* ips, std::size_t nIPs)
			{import.template set_local_ips<refDim>(ips, nIPs);}

			void set_global_ips(const MathVector<dim>* ips, std::size_t nIPs)
			{import.set_global_ips(ips, nIPs);}

			template <int refDim>
			std::vector<MathVector<refDim> >* local_ips();

			DataImport<number, dim> import;
			std::vector<MathVector<3> > vLocIP_dim3;
			std::vector<MathVector<2> > vLocIP_dim2;	// might have Neumann bnd for lower-dim elements!
			std::vector<MathVector<1> > vLocIP_dim1;
			std::vector<MathVector<dim> > vGloIP;
		};

	///	Conditional scalar user data
		struct BNDNumberData : public base_type::Data
		{
			BNDNumberData(SmartPtr<CplUserData<number, dim, bool> > functor_,
						  std::string BndSubsets, std::string InnerSubsets)
				: base_type::Data(BndSubsets, InnerSubsets), functor(functor_) {}

			SmartPtr<CplUserData<number, dim, bool> > functor;
		};

	///	Unconditional vector user data
		struct VectorData : public base_type::Data
		{
			VectorData(SmartPtr<CplUserData<MathVector<dim>, dim> > functor_,
			           std::string BndSubsets, std::string InnerSubsets)
			: base_type::Data(BndSubsets, InnerSubsets), functor(functor_) {}

			SmartPtr<CplUserData<MathVector<dim>, dim> > functor;
		};

		std::vector<NumberData> m_vNumberData;
		std::vector<BNDNumberData> m_vBNDNumberData;
		std::vector<VectorData> m_vVectorData;

		void update_subset_groups();
		using base_type::update_subset_groups;

	///	current inner subset
		int m_si;

	public:
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	protected:
	///	assembling functions for fv1
	///	\{
		template<typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);
		template<typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);
		template<typename TElem, typename TFVGeom>
		void fsh_elem_loop();
		template<typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	prepares the loop over all elements of one type for the computation of the error estimator
		template <typename TElem, typename TFVGeom>
		void prep_err_est_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for assembling the error estimator
		template <typename TElem, typename TFVGeom>
		void prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	computes the error estimator contribution for one element
		template <typename TElem, typename TFVGeom>
		void compute_err_est_rhs_elem(GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale);

	///	postprocesses the loop over all elements of one type in the computation of the error estimator
		template <typename TElem, typename TFVGeom>
		void fsh_err_est_elem_loop();


	/// \}

		static const int _C_ = 0;

	private:
		void register_all_funcs(bool bHang);
		template <typename TElem, typename TFVGeom>	void register_func();
};

} // end namespac ug

#endif /*__H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY___NEUMANN_BOUNDARY_FV1__*/
