/*
 * domain_disc.h
 *
 *  Created on: 29.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC__

// other ug4 modules
#include "common/common.h"
#include "common/util/string_util.h"

// library intern headers
#include "subset_assemble_util.h"
#include "domain_disc_interface.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_assemble_util.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "disc_item.h"
#include "lib_disc/spatial_disc/domain_disc_interface.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/error_elem_marking_strategy.h"

namespace ug {

/// \ingroup lib_disc_domain_assemble
/// @{

/// generic domain discretization implementing the interface
/**
 * This class template is an implementation of the IDomainDiscretization
 * interface based on the simple groupping of several local (element) 
 * discretizations and constraints.
 *
 * Functions of this class template prepare lists of the local discretizations
 * and elements (where to assemble) for every subset, whereas assembling itself
 * is performed by functions of the so-called 'global assembler class' specified
 * by the TGlobAssembler template parameter. The latter class implements
 * assembling for the generic lists of elements belonging to only one subset.
 * Cf. class StdGlobAssembler for a complete example. Note that all the functions
 * from that example should be implemented (no matter where as regular or static
 * members).
 *
 * \tparam TDomain          domain type
 * \tparam TAlgebra         algebra type
 * \tparam TGlobAssembler   global assembler type
 */
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
class DomainDiscretizationBase
:   public IDomainDiscretization<TAlgebra>,
    protected TGlobAssembler
{
    /// Type of the global assembler
        typedef TGlobAssembler gass_type;
        
	public:
	///	Type of Domain
		typedef TDomain domain_type;

	///	Type of algebra
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	///	Type of approximation space
		typedef ApproximationSpace<TDomain>	approx_space_type;
		
	///	world dimension
		static const int dim = TDomain::dim;
		
	public:
	///	default Constructor
		DomainDiscretizationBase(SmartPtr<approx_space_type> pApproxSpace) :
			m_spApproxSpace(pApproxSpace), m_spAssTuner(new AssemblingTuner<TAlgebra>)
		{};

	/// virtual destructor
		virtual ~DomainDiscretizationBase() {};

	///////////////////////////
	// Time independent part
	///////////////////////////

	/// \copydoc IAssemble::assemble_jacobian()
		virtual void assemble_jacobian(matrix_type& J, const vector_type& u, ConstSmartPtr<DoFDistribution> dd);
		virtual void assemble_jacobian(matrix_type& J, const vector_type& u, const GridLevel& gl)
		{assemble_jacobian(J, u, dd(gl));}

	/// \copydoc IAssemble::assemble_defect()
		virtual void assemble_defect(vector_type& d, const vector_type& u, ConstSmartPtr<DoFDistribution> dd);
		virtual void assemble_defect(vector_type& d, const vector_type& u, const GridLevel& gl)
		{assemble_defect(d, u, dd(gl));}

	/// \copydoc IAssemble::assemble_linear()
		virtual void assemble_linear(matrix_type& A, vector_type& b, ConstSmartPtr<DoFDistribution> dd);
		virtual void assemble_linear(matrix_type& mat, vector_type& rhs, const GridLevel& gl)
		{assemble_linear(mat, rhs, dd(gl));}

	/// \copydoc IAssemble::assemble_rhs()
		virtual void assemble_rhs(vector_type& rhs, const vector_type& u, ConstSmartPtr<DoFDistribution> dd);
		virtual void assemble_rhs(vector_type& rhs, const vector_type& u, const GridLevel& gl)
		{assemble_rhs(rhs, u, dd(gl));}

	/// \copydoc IAssemble::assemble_rhs()
		virtual void assemble_rhs(vector_type& rhs, ConstSmartPtr<DoFDistribution> dd)
		{assemble_rhs(rhs, rhs, dd);}
		virtual void assemble_rhs(vector_type& rhs, const GridLevel& gl)
		{assemble_rhs(rhs, dd(gl));}

	/// \copydoc IAssemble::adjust_solution()
		virtual void adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd);
		virtual void adjust_solution(vector_type& u, const GridLevel& gl)
		{adjust_solution(u, dd(gl));}

	///	wrapper for GridFunction
	/// \{
		void assemble_jacobian(matrix_type& J, GridFunction<TDomain, TAlgebra>& u)
		{assemble_jacobian(J, u, u.dof_distribution());}

		void assemble_defect(vector_type& d, GridFunction<TDomain, TAlgebra>& u)
		{assemble_defect(d, u, u.dof_distribution());}

		void assemble_linear(matrix_type& A, GridFunction<TDomain, TAlgebra>& rhs)
		{assemble_linear(A, rhs, rhs.dof_distribution());}

		void assemble_rhs(vector_type& rhs, GridFunction<TDomain, TAlgebra>& u)
		{assemble_rhs(rhs, u, u.dof_distribution());}

		void assemble_rhs(GridFunction<TDomain, TAlgebra>& b)
		{assemble_rhs(b, b.dof_distribution());}

		void adjust_solution(GridFunction<TDomain, TAlgebra>& u)
		{adjust_solution(u, u.dof_distribution());}
	/// \}

	///////////////////////
	// Time dependent part
	///////////////////////

	/// \copydoc IDomainDiscretization::prepare_timestep()
		virtual void prepare_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                              ConstSmartPtr<DoFDistribution> dd);
		virtual	void prepare_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, const GridLevel& gl)
		{prepare_timestep(vSol, dd(gl));}

	/// \copydoc IDomainDiscretization::assemble_jacobian()
		virtual void assemble_jacobian(matrix_type& J,
									   ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
									   const number s_a0,
									   ConstSmartPtr<DoFDistribution> dd);
		virtual void assemble_jacobian(matrix_type& J,
		                               ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                               const number s_a, const GridLevel& gl)
		{assemble_jacobian(J, vSol, s_a, dd(gl));}

	/// \copydoc IDomainDiscretization::assemble_defect()
		virtual void assemble_defect(vector_type& d,
									 ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
									 const std::vector<number>& vScaleMass,
									 const std::vector<number>& vScaleStiff,
									 ConstSmartPtr<DoFDistribution> dd);
		virtual	void assemble_defect(vector_type& d,
		       	                     ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		       	                     const std::vector<number>& vScaleMass,
		       	                     const std::vector<number>& vScaleStiff,
		       	                     const GridLevel& gl)
		{assemble_defect(d, vSol, vScaleMass, vScaleStiff, dd(gl));}

	/// \copydoc IDomainDiscretization::assemble_linear()
		virtual void assemble_linear(matrix_type& A, vector_type& b,
									 ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
									 const std::vector<number>& vScaleMass,
									 const std::vector<number>& vScaleStiff,
									 ConstSmartPtr<DoFDistribution> dd);
		virtual void assemble_linear(matrix_type& A, vector_type& b,
		                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                             const std::vector<number>& vScaleMass,
		                             const std::vector<number>& vScaleStiff,
		                             const GridLevel& gl)
		{assemble_linear(A, b, vSol, vScaleMass, vScaleStiff, dd(gl));}

	/// \copydoc IDomainDiscretization::assemble_rhs()
		virtual void assemble_rhs(	 vector_type& b,
									 ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
									 const std::vector<number>& vScaleMass,
									 const std::vector<number>& vScaleStiff,
									 ConstSmartPtr<DoFDistribution> dd);
		virtual void assemble_rhs(	 vector_type& b,
		                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                             const std::vector<number>& vScaleMass,
		                             const std::vector<number>& vScaleStiff,
		                             const GridLevel& gl)
		{assemble_rhs(b, vSol, vScaleMass, vScaleStiff, dd(gl));}

	/// \copydoc IDomainDiscretization::adjust_solution()
		virtual void adjust_solution(vector_type& u, number time, ConstSmartPtr<DoFDistribution> dd);
		virtual void adjust_solution(vector_type& u, number time, const GridLevel& gl)
		{adjust_solution(u, time, dd(gl));}

	/// \copydoc IDomainDiscretization::finish_timestep()
		virtual void finish_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, ConstSmartPtr<DoFDistribution> dd);
		virtual void finish_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol, const GridLevel& gl)
		{finish_timestep(vSol, dd(gl));}

	///////////////////////////
	// Mass and Stiffness Matrix
	///////////////////////////

	/// assembles the mass matrix
		virtual void assemble_mass_matrix(matrix_type& M, const vector_type& u,
		                                  ConstSmartPtr<DoFDistribution> dd);
		virtual void assemble_mass_matrix(matrix_type& M, const vector_type& u,
		                                  const GridLevel& gl)
		{assemble_mass_matrix(M, u, dd(gl));}

	/// assembles the stiffness matrix
		virtual void assemble_stiffness_matrix(matrix_type& A, const vector_type& u,
		                                       ConstSmartPtr<DoFDistribution> dd);
		virtual void assemble_stiffness_matrix(matrix_type& A, const vector_type& u,
		                                       const GridLevel& gl)
		{assemble_stiffness_matrix(A, u, dd(gl));}

	///////////////////////////////////////////////////////////
	// Error estimator										///

	/// \copydoc IDomainDiscretization::mark_error()
	// stationary
		virtual void calc_error
		(	const vector_type& u,
			ConstSmartPtr<DoFDistribution> dd,
			vector_type* u_vtk = NULL
		);
		virtual	void calc_error
		(	const vector_type& u,
			const GridLevel& gl,
			vector_type* u_vtk = NULL
		)
		{calc_error(u, dd(gl));}
		virtual	void calc_error(const GridFunction<TDomain,TAlgebra>& u)
		{calc_error(u, u.dd(), NULL);}
		virtual	void calc_error(const GridFunction<TDomain,TAlgebra>& u, vector_type* u_vtk)
		{calc_error(u, u.dd(), u_vtk);}

	// instationary
		virtual void calc_error
		(	ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
			ConstSmartPtr<DoFDistribution> dd,
			std::vector<number> vScaleMass,
			std::vector<number> vScaleStiff,
			vector_type* u_vtk
		);
		virtual void calc_error
		(	ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
			std::vector<number> vScaleMass,
			std::vector<number> vScaleStiff,
			const GridLevel& gl,
			vector_type* u_vtk
		)
		{
			calc_error((ConstSmartPtr<VectorTimeSeries<vector_type> >) vSol, dd(gl),
					   vScaleMass, vScaleStiff, u_vtk);
		}

		virtual void mark_for_refinement
		(	IRefiner& refiner,
			number TOL,
			number refineFrac,
			int maxLevel
		);

		virtual void mark_for_coarsening
		(	IRefiner& refiner,
			number TOL,
			number coarseFrac,
			int maxLevel
		);

		virtual void mark_with_strategy
		(	IRefiner& refiner,
			SmartPtr <IElementMarkingStrategy<TDomain> > strategy
		);

		/// marks error indicators as invalid; in order to revalidate them,
		/// they will have to be newly calculated by a call to calc_error
		virtual void invalidate_error();

		/// returns whether current error values are valid
		virtual bool is_error_valid();

	protected:
		typedef typename domain_traits<dim>::element_type elem_type;
		typedef MultiGrid::AttachmentAccessor<elem_type, Attachment<number> > aa_type;
		Attachment<number> m_aError;
		aa_type m_aaError;
		SmartPtr<MultiGrid> m_pMG;

		bool m_bErrorCalculated;

	// Error estimator										 //
	///////////////////////////////////////////////////////////

	public:
	/// \{
		virtual SmartPtr<AssemblingTuner<TAlgebra> > ass_tuner() {return m_spAssTuner;}
		virtual ConstSmartPtr<AssemblingTuner<TAlgebra> > ass_tuner() const {return m_spAssTuner;}
	/// \}

	public:
	/// adds an element discretization to the assembling process
	/**
	 * This function adds an Element Discretization to the assembling. During
	 * the assembling process the elem disc will assemble a given group of
	 * subsets and add the output to the passed vectors and matrices.
	 * For each Elem Disc one loop over all elements of the subset will be
	 * performed.
	 *
	 * \param[in] 	elem		Element Discretization to be added
	 */
		void add(SmartPtr<IElemDisc<TDomain> > elem)
		{
		//	check that not already registered
			for(size_t i = 0; i < m_vDomainElemDisc.size(); ++i)
				if(m_vDomainElemDisc[i] == elem)
					return;

		//	set approximation space
			elem->set_approximation_space(m_spApproxSpace);

		//	add it
			m_vDomainElemDisc.push_back(elem);
		}

	/// adds a constraint to the assembling process
	/**
	 * This function adds a IConstraint to the assembling. The constraint is
	 * called when all element-wise assembling have been performed.
	 *
	 * \param[in] 	pp		Constraint to be added
	 */
		void add(SmartPtr<IDomainConstraint<TDomain, TAlgebra> > pp)
		{
		//	check that not already registered
			for(size_t i = 0; i < m_vConstraint.size(); ++i)
				if(m_vConstraint[i] == pp)
					return;

		//	add constraint
			m_vConstraint.push_back(pp);
		}

	/// removes a constraint from the assembling process
	/**
	 * This function removes a previously added IConstraint from the assembling.
	 * The constraint's assemblings are no longer called.
	 *
	 * \param[in] 	pp		constraint to be removed
	 */
		void remove(SmartPtr<IDomainConstraint<TDomain, TAlgebra> > pp)
		{
			// check that already registered
			for (size_t i = 0; i < m_vConstraint.size(); i++)
			{
				if (m_vConstraint[i] == pp)
				{
					// remove constraint
					m_vConstraint.erase(m_vConstraint.begin()+i);
					return;
				}
			}

			UG_LOG("Tried to remove DomainConstraint from DomainDisc"
					", but could not find it there.");
		}

	/// adds a disc item to the assembling process
	/**
	 * This function adds a IDiscretizationItem to the assembling. The contained
	 * elem discs and constraints are sorted into the lists
	 *
	 * \param[in] 	di		Disc Item
	 */
		void add(SmartPtr<IDiscretizationItem<TDomain, TAlgebra> > di)
		{
		//	add elem discs
			for(size_t i = 0; i < di->num_elem_disc(); ++i)
				add(di->elem_disc(i));

		//	add constraints
			for(size_t i = 0; i < di->num_constraint(); ++i)
				add(di->constraint(i));
		}

	///	returns number of registered constraints
		virtual size_t num_constraints() const {return m_vConstraint.size();}

	///	returns the i'th constraint
		virtual SmartPtr<IConstraint<TAlgebra> > constraint(size_t i) {return m_vConstraint[i];}

	protected:
	///	set the approximation space in the elem discs and extract IElemDiscs
		void update_elem_discs();
		void update_constraints();
		void update_disc_items();

	protected:
	///	returns the level dof distribution
		ConstSmartPtr<DoFDistribution> dd(const GridLevel& gl) const{return m_spApproxSpace->dof_distribution(gl);}

	protected:
	///	vector holding all registered elem discs
		std::vector<SmartPtr<IElemDisc<TDomain> > > m_vDomainElemDisc;

	///	vector holding all registered elem discs
		std::vector<IElemDisc<TDomain>*> m_vElemDisc;

	//	vector holding all registered constraints
		std::vector<SmartPtr<IDomainConstraint<TDomain, TAlgebra> > > m_vConstraint;

	///	current approximation space
		SmartPtr<approx_space_type> m_spApproxSpace;
		
	///	this object provides tools to adapt the assemble routine
		SmartPtr<AssemblingTuner<TAlgebra> > m_spAssTuner;
	
	private:
	//---- Auxiliary function templates for the assembling ----//
	//	These functions call the corresponding functions from the global assembler for a composed list of elements:
	//-- for stationary problems --//
	template <typename TElem>
	void AssembleMassMatrix(		const std::vector<IElemDisc<domain_type>*>& vElemDisc,
									ConstSmartPtr<DoFDistribution> dd,
									int si, bool bNonRegularGrid,
									matrix_type& M,
									const vector_type& u);
	template <typename TElem>
	void AssembleStiffnessMatrix(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
									ConstSmartPtr<DoFDistribution> dd,
									int si, bool bNonRegularGrid,
									matrix_type& A,
									const vector_type& u);
	template <typename TElem>
	void AssembleJacobian(			const std::vector<IElemDisc<domain_type>*>& vElemDisc,
									ConstSmartPtr<DoFDistribution> dd,
									int si, bool bNonRegularGrid,
									matrix_type& J,
									const vector_type& u);
	template <typename TElem>
	void AssembleDefect( 			const std::vector<IElemDisc<domain_type>*>& vElemDisc,
									ConstSmartPtr<DoFDistribution> dd,
									int si, bool bNonRegularGrid,
									vector_type& d,
									const vector_type& u);
	template <typename TElem>
	void AssembleLinear( 			const std::vector<IElemDisc<domain_type>*>& vElemDisc,
									ConstSmartPtr<DoFDistribution> dd,
									int si, bool bNonRegularGrid,
									matrix_type& A,
									vector_type& rhs);
	template <typename TElem>
	void AssembleRhs(				const std::vector<IElemDisc<domain_type>*>& vElemDisc,
									ConstSmartPtr<DoFDistribution> dd,
									int si, bool bNonRegularGrid,
									vector_type& rhs,
									const vector_type& u);
	template <typename TElem>
	void AssembleErrorEstimator(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
									ConstSmartPtr<DoFDistribution> dd,
									int si, bool bNonRegularGrid,
									const vector_type& u);
	//-- for time-dependent problems --//
	template <typename TElem>
	void PrepareTimestep(			const std::vector<IElemDisc<domain_type>*>& vElemDisc,
									ConstSmartPtr<DoFDistribution> dd,
									int si, bool bNonRegularGrid,
									ConstSmartPtr<VectorTimeSeries<vector_type> > vSol);
	template <typename TElem>
	void AssembleJacobian(			const std::vector<IElemDisc<domain_type>*>& vElemDisc,
									ConstSmartPtr<DoFDistribution> dd,
									int si, bool bNonRegularGrid,
									matrix_type& J,
									ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
									number s_a0);
	template <typename TElem>
	void AssembleDefect( 			const std::vector<IElemDisc<domain_type>*>& vElemDisc,
									ConstSmartPtr<DoFDistribution> dd,
									int si, bool bNonRegularGrid,
									vector_type& d,
									ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
									const std::vector<number>& vScaleMass,
									const std::vector<number>& vScaleStiff);
	template <typename TElem>
	void AssembleLinear( 			const std::vector<IElemDisc<domain_type>*>& vElemDisc,
									ConstSmartPtr<DoFDistribution> dd,
									int si, bool bNonRegularGrid,
									matrix_type& A,
									vector_type& rhs,
									ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
									const std::vector<number>& vScaleMass,
									const std::vector<number>& vScaleStiff);
	template <typename TElem>
	void AssembleRhs(				const std::vector<IElemDisc<domain_type>*>& vElemDisc,
									ConstSmartPtr<DoFDistribution> dd,
									int si, bool bNonRegularGrid,
									vector_type& rhs,
									ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
									const std::vector<number>& vScaleMass,
									const std::vector<number>& vScaleStiff);
	template <typename TElem>
	void AssembleErrorEstimator(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
									ConstSmartPtr<DoFDistribution> dd,
									int si, bool bNonRegularGrid,
									std::vector<number> vScaleMass,
									std::vector<number> vScaleStiff,
									ConstSmartPtr<VectorTimeSeries<vector_type> > vSol);
	template <typename TElem>
	void FinishTimestep(			const std::vector<IElemDisc<domain_type>*>& vElemDisc,
									ConstSmartPtr<DoFDistribution> dd,
									int si, bool bNonRegularGrid,
									ConstSmartPtr<VectorTimeSeries<vector_type> > vSol);
};

/// domain discretization implementing the interface
/**
 * This class template is an implementation of the IDomainDiscretization
 * interface based on the simple groupping of several local (element) 
 * discretizations and constraints.
 *
 * This is the "usual" and "trivial" global discretizations: all the
 * local (element) discretizations are perfomed once for every element, and
 * their contributions are algebraically added to the global data.
 *
 * \tparam TDomain          domain type
 * \tparam TAlgebra         algebra type
 */
template <typename TDomain, typename TAlgebra>
class DomainDiscretization
:	public DomainDiscretizationBase<TDomain, TAlgebra, StdGlobAssembler<TDomain, TAlgebra> >
{
    /// Type of the global assembler
        typedef StdGlobAssembler<TDomain, TAlgebra> gass_type;
        
	public:
	///	Type of Domain
		typedef TDomain domain_type;

	///	Type of algebra
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	///	Type of approximation space
		typedef ApproximationSpace<TDomain>	approx_space_type;
		
	///	world dimension
		static const int dim = TDomain::dim;
		
	public:
	///	default Constructor
		DomainDiscretization(SmartPtr<approx_space_type> pApproxSpace)
		: DomainDiscretizationBase<domain_type, algebra_type, gass_type> (pApproxSpace)
		{};

	/// virtual destructor
		virtual ~DomainDiscretization() {};
};

/// @}

} // end namespace ug

// include documentation
#include "domain_disc_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC__ */
