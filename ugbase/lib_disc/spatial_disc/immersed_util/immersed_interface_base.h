/*
 * moving_interface.h
 *
 *  Created on: 15.01.2015
 *      Author: susanne hoellbacher
 */

#ifndef IMMERSED_INTERFACE_BASE_H_
#define IMMERSED_INTERFACE_BASE_H_

#include "lib_disc/assemble_interface.h"
#include "lib_disc/spatial_disc/local_to_global/local_to_global_mapper.h"
#include "interface_handler/interface_handler_base.h"
#include "lib_disc/spatial_disc/disc_util/fv1Cut_geom.h"
#include "lib_disc/spatial_disc/disc_util/fv1FT_geom.h"

namespace ug{

//////////////////////////////////////////////////////////////////////////////////
// class 'IInterfaceMapper': member of class 'IMovingInterface' (see below)
//////////////////////////////////////////////////////////////////////////////////

template <typename TAlgebra>
class IInterfaceMapper : public ILocalToGlobalMapper<TAlgebra>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	///	default Constructor
		IInterfaceMapper(){};

 		IInterfaceMapper(SmartPtr<IInterfaceHandlerLocal> localHandler);

	///	send local entries to global matrix
		virtual void add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat,
				ConstSmartPtr<DoFDistribution> dd);
		virtual void add_local_mat_to_global_interface(matrix_type& mat, const LocalMatrix& lmat,
				ConstSmartPtr<DoFDistribution> dd){};
		virtual void add_local_mat_to_global_interface_for2(matrix_type& mat, const LocalMatrix& lmat,
				ConstSmartPtr<DoFDistribution> dd){};

	///	send local entries to global rhs
		virtual void add_local_vec_to_global(vector_type& vec, const LocalVector& lvec,
				ConstSmartPtr<DoFDistribution> dd);
		virtual void add_local_vec_to_global_interface(vector_type& vec, const LocalVector& lvec,
				ConstSmartPtr<DoFDistribution> dd){}
		virtual void add_local_vec_to_global_interface_for2(vector_type& vec, const LocalVector& lvec,
					ConstSmartPtr<DoFDistribution> dd){}

	///	modifies local solution vector for adapted defect computation
 		virtual void modify_LocalData(LocalMatrix& locJ, LocalVector& locU,
 				ConstSmartPtr<DoFDistribution> dd){};
 		virtual void modify_LocalData(LocalVectorTimeSeries& uT, LocalMatrix& locJ, LocalVector& locU,
 				ConstSmartPtr<DoFDistribution> dd){};

 		virtual void modify_LocalData(LocalVector& locD, LocalVector& tmpLocD, LocalVector& locU,
 				ConstSmartPtr<DoFDistribution> dd){};
 		virtual void modify_LocalData(LocalVectorTimeSeries& uT, LocalVector& locD, LocalVector& tmpLocD, LocalVector& locU,
 				ConstSmartPtr<DoFDistribution> dd, size_t t){};

  	///	modifies global solution vector for adapted defect computation
 		virtual void modify_GlobalSol(vector_type& vecMod, const vector_type& vec,
 				ConstSmartPtr<DoFDistribution> dd){};

 		virtual void modify_GlobalSol(SmartPtr<VectorTimeSeries<vector_type> > vSolMod,
 	      		ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
 	      		ConstSmartPtr<DoFDistribution> dd){};

		virtual void set_identity_mat(matrix_type& mat, const LocalMatrix& lmat,
				ConstSmartPtr<DoFDistribution> dd);

	///	virtual destructor
		virtual ~IInterfaceMapper() {};

	private:
 		SmartPtr<IInterfaceHandlerLocal> m_spInterfaceHandlerLocal;


};

//////////////////////////////////////////////////////////////////////////////////
// class 'IInterfaceBndCond': member of class 'IMovingInterface' (see below)
//////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
class IInterfaceBndCond : public IElemDisc<TDomain>
{
	///	world Dimension
		static const int dim = TDomain::dim;

	public:
	///	defautl constructor
		IInterfaceBndCond(const char* functions, const char* subsets,
							SmartPtr<IInterfaceHandlerLocal> localHandler);
		IInterfaceBndCond(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset,
						  SmartPtr<IInterfaceHandlerLocal> localHandler);
	///	virtual destructor
	 ~IInterfaceBndCond() {};
/*
	public:
		///	type of trial space for each function used
			void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid, bool bStaticRoid)
			{ UG_THROW("'IInterfaceBndCond::prepare_setting' not implemented!\n"); }

 		///	prepares the element loop
			template <typename TElem, typename TFVGeom>
			void prep_elem_loop(const ReferenceObjectID roid, const int si)
			{ UG_THROW("'IInterfaceBndCond::prep_elem_loop' not implemented!\n"); }

		///	prepares the element for evaluation
			template <typename TElem, typename TFVGeom>
			void prep_elem(const LocalVector& u, GridObject* elem,
					const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
			{ UG_THROW("'IInterfaceBndCond::prep_elem' not implemented!\n"); }

		///	finishes the element loop
			template <typename TElem, typename TFVGeom>
			void fsh_elem_loop()
			{ UG_THROW("'IInterfaceBndCond::fsh_elem_loop' not implemented!\n"); }

		///	adds the stiffness part to the local jacobian
			template <typename TElem, typename TFVGeom>
			void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem,
					const MathVector<dim> vCornerCoords[])
			{ UG_THROW("'IInterfaceBndCond::add_jac_A_elem' not implemented!\n"); }

		///	adds the stiffness part to the local defect
			template <typename TElem, typename TFVGeom>
			void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem,
					const MathVector<dim> vCornerCoords[])
			{ UG_THROW("'IInterfaceBndCond::add_def_A_elem' not implemented!\n"); }

			template <typename TElem, typename TFVGeom>
			void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem,
					const MathVector<dim> vCornerCoords[])
			{ UG_THROW("'IInterfaceBndCond::add_jac_M_elem' not implemented!\n"); }

			template <typename TElem, typename TFVGeom>
			void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem,
					const MathVector<dim> vCornerCoords[])
			{ UG_THROW("'IInterfaceBndCond::add_def_M_elem' not implemented!\n"); }

			template <typename TElem, typename TFVGeom>
			void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
			{ UG_THROW("'IInterfaceBndCond::add_rhs_elem' not implemented!\n"); }

*/
	private:
		SmartPtr<IInterfaceHandlerLocal> m_spInterfaceHandlerLocal;

};


//////////////////////////////////////////////////////////////////////////////////
// main class 'IMovingInterface'
//////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
class IMovingInterface
{
	public:
	///	world Dimension
		static const int dim = TDomain::dim;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	// default constructor
		IMovingInterface(){};

		IMovingInterface(SmartPtr<IAssemble<TAlgebra> > ass,
				   	   	 const char* functions, const char* subsets,
				   	 	 SmartPtr<IInterfaceHandlerLocal> localHandler);
		IMovingInterface(SmartPtr<IAssemble<TAlgebra> > ass,
				   	     const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset,
				   	 	 SmartPtr<IInterfaceHandlerLocal> localHandler);

		virtual ~IMovingInterface()	{}


	private:

		SmartPtr<IInterfaceHandlerLocal> m_spInterfaceHandlerLocal;
		SmartPtr<IInterfaceMapper<TAlgebra> > m_spInterfaceMapper;  	// contains member of class 'IInterfaceHandlerLocal'
		SmartPtr<IInterfaceBndCond<TDomain> > m_spInterfaceBndCond;		// contains member of class 'IInterfaceHandlerLocal'


};


}// end namespace ug

#include "immersed_interface_base_impl.h"

#endif /* IMMERSED_INTERFACE_BASE_H_ */
