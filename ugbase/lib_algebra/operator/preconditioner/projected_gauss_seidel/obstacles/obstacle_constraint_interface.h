/*
 * obstacle_constraint_interface.h
 *
 *  Created on: 25.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_CONSTRAINT_INTERFACE__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_CONSTRAINT_INTERFACE__

#include "lib_disc/common/multi_index.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"

#include "lib_disc/spatial_disc/constraints/constraint_interface.h"

using namespace std;

namespace ug{

/// Interface for Obstacle Constraints
/**
 *  The Interface for Obstacle Constraints provides the framework to define obstacle constraints
 *  of the form
 *
 * 			c(u) <= upObs 		(cf. 'set_upper_obstacle')
 *
 * 	and
 *
 * 			c(u) >= lowObs 		(cf. 'set_lower_obstacle')
 *
 * 	where u is the solution vector. Here, 'upObs' and 'lowObs' are user-defined functions,
 * 	which need to be of the same size as the function of unknowns u. The obstacle function
 * 	c is defined by creating a derived class of this interface and by specializing the way
 * 	the correction should be computed (cf. correction_for_lower_obs, correction_for_upper_obs,
 * 	correction_for_lower_and_upper_obs).
 *
 * 	One simple example is the ScalarObstacle-class. Such an obstacle functions can be used in
 * 	combination with projected preconditioners. They should be passed to the preconditioner
 * 	by 'IProjPreconditioner::set_obstacle_constraint'.
 */
template <typename TDomain, typename TAlgebra>
class IObstacleConstraint:
	public IDomainConstraint<TDomain, TAlgebra>
	//	TODO: think about, restructuring the IConstraint-Interface in order to distinguish between
	//	a constraint having an effect on the assembling process and a constraint, which is considered in the
	//	solver (e.g. by means of 'adjust_restriction') only!
{
	public:
	///	Base Type
		typedef IDomainConstraint<TDomain, TAlgebra> base_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Matrix type
		typedef typename algebra_type::matrix_type matrix_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	///	Value type
		typedef typename vector_type::value_type value_type;

	///	Type of domain
		typedef TDomain domain_type;

	///	world Dimension
		static const int dim = domain_type::dim;

	///	Type of position coordinates (e.g. position_type)
		typedef typename domain_type::position_type position_type;

	public:
	/// constructor for an obstacle defined on some subset(s)
		IObstacleConstraint(const GridFunction<TDomain, TAlgebra>& u)
		{
			clear();

			m_spDD = u.dof_distribution();
			m_spDomain = u.domain();
		};

	/// constructor
		IObstacleConstraint(){
			//clear();
			UG_THROW("In 'IObstacleConstraint()': A constructor with a GridFunction as parameter"
					"is needed here!");
		};

	///	adds a lua callback (cond and non-cond)
	#ifdef UG_FOR_LUA
		void add(const char* name, const char* function);
		void add(const char* name, const char* function, const char* subsets);
	#endif

	///	adds a conditional user-defined value as dirichlet condition for a function on subsets and on whole domain
		void add(SmartPtr<UserData<number, dim, bool> > func, const char* function);
		void add(SmartPtr<UserData<number, dim, bool> > func, const char* function, const char* subsets);

	///	adds a user-defined value as dirichlet condition for a function on subsets and on whole domain
		void add(SmartPtr<UserData<number, dim> > func, const char* function);
		void add(SmartPtr<UserData<number, dim> > func, const char* function, const char* subsets);

	///	adds a constant value as dirichlet condition for a function on subsets and on whole domain
		void add(number value, const char* function);
		void add(number value, const char* function, const char* subsets);

	///	adds a user-defined vector as dirichlet condition for a vector-function on subsets and on whole domain
		void add(SmartPtr<UserData<MathVector<dim>, dim> > func, const char* functions);
		void add(SmartPtr<UserData<MathVector<dim>, dim> > func, const char* functions, const char* subsets);


	///	init function calls 'extract_data' and 'obstacle_value'-methods in order to
	///	store the obstacle values set by UserDatas
		void init();

	///	checks if a given dof is in an obstacle subset
		bool is_obs_dof(const DoFIndex& dof);

	///	resets the vector storing the active dofs
		void reset_active_dofs(){m_vActiveDofs.resize(0);}

	///	returns the vector storing the active dofs
		void active_dofs(vector<DoFIndex>& vActiveDoFs)
		{vActiveDoFs = m_vActiveDofs;}


	///	this preprocess-method is called in every init-call of the underlying proj. preconditioner
	///	it is useful to attach e.g. additional data to the obstacle DoFs, like the e.g. the normal vector at
	///	this DoF
		virtual void preprocess(){};

	///	projects the i-th index of the solution onto the admissible set and adjusts the correction
		virtual void adjust_sol_and_cor(value_type& sol_i, value_type& c_i, bool& dofIsAdmissible,
				const DoFIndex& dof) = 0;

	///	the defect needs to be adjusted for the active indices (those indices, which are in contact)
		virtual void adjust_defect_to_constraint(vector_type& d) = 0;


	///	restricts the obstacle values to a coarser grid in a multigrid hierarchy
		virtual void restrict_obs_values() = 0;


	///	Destructor
		virtual ~IObstacleConstraint(){};

	public:
	///////////////////////////////
	// 	Implement Interface
	///////////////////////////////

	/// sets a unity row for all dirichlet indices
		void adjust_jacobian(matrix_type& J, const vector_type& u,
		                     ConstSmartPtr<DoFDistribution> dd, number time = 0.0,
                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = NULL,
							 const number s_a0 = 1.0){
			UG_LOG("IObstacleConstraint::adjust_jacobian() \n");
		};

	/// sets a zero value in the defect for all dirichlet indices
		void adjust_defect(vector_type& d, const vector_type& u,
		                   ConstSmartPtr<DoFDistribution> dd, number time = 0.0,
                           ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = NULL,
						   const std::vector<number>* vScaleMass = NULL,
						   const std::vector<number>* vScaleStiff = NULL){
			UG_LOG("IObstacleConstraint::adjust_defect() \n");
		};

	/// sets the dirichlet value in the solution for all dirichlet indices
		void adjust_solution(vector_type& u,
		                     ConstSmartPtr<DoFDistribution> dd, number time = 0.0){
			UG_LOG("IObstacleConstraint::adjust_solution() \n");
		};

	///	sets unity rows in A and dirichlet values in right-hand side b
		void adjust_linear(matrix_type& A, vector_type& b,
		                   ConstSmartPtr<DoFDistribution> dd, number time = 0.0){
			UG_LOG("IObstacleConstraint::adjust_linear() \n");
		};

	///	sets the dirichlet value in the right-hand side
		void adjust_rhs(vector_type& b, const vector_type& u,
		                ConstSmartPtr<DoFDistribution> dd, number time = 0.0){
			UG_LOG("IObstacleConstraint::adjust_rhs() \n");
		};

	///	sets constraints in prolongation
		virtual void adjust_prolongation(matrix_type& P,
										 ConstSmartPtr<DoFDistribution> ddFine,
										 ConstSmartPtr<DoFDistribution> ddCoarse,
										 number time = 0.0){
			UG_LOG("IObstacleConstraint::adjust_prolongationP() \n");
		};

	///	sets constraints in restriction
		virtual void adjust_restriction(matrix_type& R,
										ConstSmartPtr<DoFDistribution> ddCoarse,
										ConstSmartPtr<DoFDistribution> ddFine,
										number time = 0.0);

	///	returns the type of the constraints
		//TODO: does another type make more sense?
		virtual int type() const {return CT_CONSTRAINTS;}

	private:
	///	extract the UserDatas
		void extract_data();

		template <typename TUserData, typename TScheduledUserData>
		void extract_data(std::map<int, std::vector<TUserData*> >& mvUserDataObsSegment,
		                  std::vector<TScheduledUserData>& vUserData);

	///	clear all UserData-member variables
		void clear();

		void check_functions_and_subsets(FunctionGroup& functionGroup, SubsetGroup& subsetGroup,
				size_t numFct) const;

	///	store the obstacle value set by means of UserDatas
		void init_obstacle_dofs_with_values(number time);

		template <typename TUserData>
		void init_obstacle_dofs_with_values(const std::map<int, std::vector<TUserData*> >& mvUserData, number time);

		template <typename TBaseElem, typename TUserData>
		void init_obstacle_dofs_with_values(const std::vector<TUserData*>& vUserData, int si, number time);

	protected:
	///	map to store obstacle values with its corresponding DoFs
		map<DoFIndex, number> m_mObstacleValues;

	///	stores the subset-indices of the obstacle subsets
		vector<int> m_vObsSubsets;

	///	stores the dofs, which satisfy the constraints with equality
		vector<DoFIndex> m_vActiveDofs;

	private:
	///	pointer to the domain
		ConstSmartPtr<TDomain> m_spDomain;

	///	pointer to the DofDistribution on the whole domain
		ConstSmartPtr<DoFDistribution> m_spDD;

	///	grouping for subset and non-conditional data
		struct NumberData
		{
			const static bool isConditional = false;
			const static size_t numFct = 1;
			typedef MathVector<1> value_type;
			NumberData(SmartPtr<UserData<number, dim> > functor_,
					   std::string fctName_)
				: spFunctor(functor_), fctName(fctName_), bWholeDomain(true)
			{}
			NumberData(SmartPtr<UserData<number, dim> > functor_,
					   std::string fctName_, std::string ssName_)
				: spFunctor(functor_), fctName(fctName_), ssName(ssName_),
				  bWholeDomain(false)
			{}
			bool operator()(MathVector<1>& val, const MathVector<dim> x,
							number time, int si) const
			{
				(*spFunctor)(val[0], x, time, si); return true;
			}

			SmartPtr<UserData<number, dim> > spFunctor;
			std::string fctName;
			std::string ssName;
			size_t fct[numFct];
			SubsetGroup ssGrp;
			bool bWholeDomain;
		};

	///	grouping for subset and conditional data
		struct CondNumberData
		{
			const static bool isConditional = true;
			const static size_t numFct = 1;
			typedef MathVector<1> value_type;
			CondNumberData(SmartPtr<UserData<number, dim, bool> > functor_,
						  std::string fctName_)
				: spFunctor(functor_), fctName(fctName_), bWholeDomain(true)
			{}
			CondNumberData(SmartPtr<UserData<number, dim, bool> > functor_,
						  std::string fctName_, std::string ssName_)
				: spFunctor(functor_), fctName(fctName_), ssName(ssName_),
				  bWholeDomain(true)
			{}
			bool operator()(MathVector<1>& val, const MathVector<dim> x,
							number time, int si) const
			{
				return (*spFunctor)(val[0], x, time, si);
			}

			SmartPtr<UserData<number, dim, bool> > spFunctor;
			std::string fctName;
			std::string ssName;
			size_t fct[numFct];
			SubsetGroup ssGrp;
			bool bWholeDomain;
		};

	///	grouping for subset and conditional data
		struct ConstNumberData
		{
			const static bool isConditional = false;
			const static size_t numFct = 1;
			typedef MathVector<1> value_type;
			ConstNumberData(number value_,
						  std::string fctName_)
				: functor(value_), fctName(fctName_), bWholeDomain(true)
			{}
			ConstNumberData(number value_,
						  std::string fctName_, std::string ssName_)
				: functor(value_), fctName(fctName_), ssName(ssName_),
				  bWholeDomain(false)
			{}
			inline bool operator()(MathVector<1>& val, const MathVector<dim> x,
								   number time, int si) const
			{
				val[0] = functor; return true;
			}

			number functor;
			std::string fctName;
			std::string ssName;
			size_t fct[numFct];
			SubsetGroup ssGrp;
			bool bWholeDomain;
		};

	///	grouping for subset and non-conditional data
		struct VectorData
		{
			const static bool isConditional = false;
			const static size_t numFct = dim;
			typedef MathVector<dim> value_type;
			VectorData(SmartPtr<UserData<MathVector<dim>, dim> > value_,
					   std::string fctName_)
				: spFunctor(value_), fctName(fctName_), bWholeDomain(true)
			{}
			VectorData(SmartPtr<UserData<MathVector<dim>, dim> > value_,
					   std::string fctName_, std::string ssName_)
				: spFunctor(value_), fctName(fctName_), ssName(ssName_),
				  bWholeDomain(false)
			{}
			bool operator()(MathVector<dim>& val, const MathVector<dim> x,
							number time, int si) const
			{
				(*spFunctor)(val, x, time, si); return true;
			}

			SmartPtr<UserData<MathVector<dim>, dim> > spFunctor;
			std::string fctName;
			std::string ssName;
			size_t fct[numFct];
			SubsetGroup ssGrp;
			bool bWholeDomain;
		};

		std::vector<CondNumberData> m_vCondNumberData;
		std::vector<NumberData> m_vNumberData;
		std::vector<ConstNumberData> m_vConstNumberData;

		std::vector<VectorData> m_vVectorData;

	///	non-conditional obstacle values for all subsets
		std::map<int, std::vector<NumberData*> > m_mNumberObsSegment;

	///	constant obstacle values for all subsets
		std::map<int, std::vector<ConstNumberData*> > m_mConstNumberObsSegment;

	///	conditional obstacle values for all subsets
		std::map<int, std::vector<CondNumberData*> > m_mCondNumberObsSegment;

	///	non-conditional obstacle values for all subsets
		std::map<int, std::vector<VectorData*> > m_mVectorObsSegment;
};


} // end namespace ug

// include implementation
#include "obstacle_constraint_interface_impl.h"

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_CONSTRAINT_INTERFACE__ */

