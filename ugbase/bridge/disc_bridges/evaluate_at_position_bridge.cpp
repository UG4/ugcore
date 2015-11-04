
// extern headers
#include <iostream>
#include <sstream>
#include <string>
#include <limits>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

// lib_disc includes
#include "lib_disc/dof_manager/dof_distribution.h"
#include "lib_disc/function_spaces/grid_function.h"

#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/distributed_grid.h"
#endif

using namespace std;

namespace ug{
namespace bridge{
namespace Evaluate{

template <typename TDomain>
bool CloseVertexExists(const MathVector<TDomain::dim>& globPos,
					   TDomain* dom,
					   const char* subsets,
					   SmartPtr<typename TDomain::subset_handler_type> sh,
					   number maxDist)
{
//	domain type
	typedef TDomain domain_type;
	typedef typename domain_type::grid_type grid_type;
	typedef typename domain_type::subset_handler_type subset_handler_type;
// get position accessor
	grid_type* grid = dom->grid().get();

	const typename domain_type::position_accessor_type& aaPos
										= dom->position_accessor();

	typename subset_handler_type::template traits<Vertex>::const_iterator iterEnd, iter;
	number minDistanceSq = numeric_limits<number>::max();

	#ifdef UG_PARALLEL
		DistributedGridManager* dgm = grid->distributed_grid_manager();
	#endif

	SubsetGroup ssGrp(sh);
	if(subsets != NULL)
		ssGrp.add(TokenizeString(subsets));
	else
		ssGrp.add_all();

	for(size_t i = 0; i < ssGrp.size(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];
	// 	iterate over all elements
		for(size_t lvl = 0; lvl < sh->num_levels(); ++lvl){
			iterEnd = sh->template end<Vertex>(si, lvl);
			iter = sh->template begin<Vertex>(si, lvl);
			for(; iter != iterEnd; ++iter)
			{
			//	get element
			//todo: replace most of the following checks by a spGridFct->contains(...)

				Vertex* vrt = *iter;
				if(grid->has_children(vrt)) continue;

				#ifdef UG_PARALLEL
					if(dgm->is_ghost(vrt))	continue;
					if(dgm->contains_status(vrt, INT_H_SLAVE)) continue;
				#endif
			//	global position
				number buffer = VecDistanceSq(globPos, aaPos[vrt]);
				if(buffer < minDistanceSq)
				{
					minDistanceSq = buffer;
				}
			}
		}
	}
	return 	minDistanceSq < sq(maxDist);
}

/**
 * \defgroup interpolate_bridge Interpolation Bridge
 * \ingroup disc_bridge
 * \{
 */

template <typename TGridFunction>
number EvaluateAtVertex(const MathVector<TGridFunction::dim>& globPos,
                           SmartPtr<TGridFunction> spGridFct,
                           size_t fct,
                           const SubsetGroup& ssGrp,
                           typename TGridFunction::domain_type::subset_handler_type* sh)
{

//	domain type
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::grid_type grid_type;
	typedef typename domain_type::subset_handler_type subset_handler_type;
// get position accessor
	domain_type* dom = spGridFct->domain().get();
	grid_type* grid = dom->grid().get();

	subset_handler_type* domSH = dom->subset_handler().get();

	const typename domain_type::position_accessor_type& aaPos
										= dom->position_accessor();

	std::vector<DoFIndex> ind;
	typename subset_handler_type::template traits<Vertex>::const_iterator iterEnd, iter,chosen;
	double minDistanceSq = 0;

	#ifdef UG_PARALLEL
		DistributedGridManager* dgm = grid->distributed_grid_manager();
	#endif

	bool bInit = false;
	for(size_t i = 0; i < ssGrp.size(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	// 	iterate over all elements
		for(size_t lvl = 0; lvl < sh->num_levels(); ++lvl){
			iterEnd = sh->template end<Vertex>(si, lvl);
			iter = sh->template begin<Vertex>(si, lvl);
			for(; iter != iterEnd; ++iter)
			{
			//	get element
			//todo: replace most of the following checks by a spGridFct->contains(...)
				Vertex* vrt = *iter;
				if(grid->has_children(vrt)) continue;

				#ifdef UG_PARALLEL
					if(dgm->is_ghost(vrt))	continue;
					if(dgm->contains_status(vrt, INT_H_SLAVE)) continue;
				#endif

				int domSI = domSH->get_subset_index(vrt);

			//	skip if function is not defined in subset
				if(!spGridFct->is_def_in_subset(fct, domSI)) continue;

			//	global position
				if(!bInit)
				{
					bInit = true;
					minDistanceSq = VecDistanceSq(globPos, aaPos[vrt]);
					chosen = iter;
				}
				else
				{
					double buffer = VecDistanceSq(globPos, aaPos[vrt]);
					if(buffer < minDistanceSq)
					{
						minDistanceSq = buffer;
						chosen = iter;
					}
				}
			}
		}
	}

	Vertex* vrt = *chosen;
	spGridFct->inner_dof_indices(vrt, fct, ind);
	return 	DoFRef(*spGridFct, ind[0]);

}


template <typename TGridFunction>
number EvaluateAtClosestVertex(const MathVector<TGridFunction::dim>& pos,
                                  SmartPtr<TGridFunction> spGridFct, const char* cmp,
                                  const char* subsets,
                                  SmartPtr<typename TGridFunction::domain_type::subset_handler_type> sh)
{
//	get function id of name
	const size_t fct = spGridFct->fct_id_by_name(cmp);

//	check that function found
	if(fct > spGridFct->num_fct())
		UG_THROW("Evaluate: Name of component '"<<cmp<<"' not found.");


//	create subset group
	SubsetGroup ssGrp(sh);
	if(subsets != NULL)
	{
		ssGrp.add(TokenizeString(subsets));
	}
	else
	{
	//	add all subsets and remove lower dim subsets afterwards
		ssGrp.add_all();
	}

	return EvaluateAtVertex<TGridFunction>(pos, spGridFct, fct, ssGrp, sh.get());
}


/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra dependent parts.
 * All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

//	typedef
	typedef ug::GridFunction<TDomain, TAlgebra> TFct;

//	Interpolate
	{
	//	reg.add_function("Integral", static_cast<number (*)(SmartPtr<UserData<number,dim> >, SmartPtr<TFct>, const char*, number, int)>(&Integral<TFct>), grp, "Integral", "Data#GridFunction#Subsets#Time#QuadOrder");

		//reg.add_function("EvaluateAtClosestVertex", static_cast<number (*)(const std::vector<number>&, SmartPtr<TFct>, const char*, const char*)>(&EvaluateAtClosestVertex<TFct>),grp, "Evaluate_at_closest_vertex", "Position#GridFunction#Component#Subsets");
		reg.add_function("EvaluateAtClosestVertex", &EvaluateAtClosestVertex<TFct>, grp, "Evaluate_at_closest_vertex", "Position#GridFunction#Component#Subsets");
	}
}

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();
	reg.add_function("CloseVertexExists", &CloseVertexExists<TDomain>, grp);
}

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

}

/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

}

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{
}

}; // end Functionality

// end group
/// \}

}// namespace Evaluate

///
void RegisterBridge_Evaluate(Registry& reg, string grp)
{
	grp.append("/Evaluate");
	typedef Evaluate::Functionality Functionality;

	try{
//		RegisterCommon<Functionality>(reg,grp);
//		RegisterDimensionDependent<Functionality>(reg,grp);
		RegisterDomainDependent<Functionality>(reg,grp);
//		RegisterAlgebraDependent<Functionality>(reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}//	end of namespace bridge
}//	end of namespace ug
