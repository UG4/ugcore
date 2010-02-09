/*
 * trialspace.h
 *
 *  Created on: 12.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__TRIALSPACE__
#define __H__LIBDISCRETIZATION__TRIALSPACE__

// extern libraries
#include <cassert>
#include <map>

// other ug4 modules
#include "common/math/ugmath.h"
#include "lib_grid/lib_grid.h"

// library intern headers
#include "referenceelement.h"

namespace ug {


/*
 * This class returns the number of dofs needed in a geometric object.
 */
class ElementDoFPattern
{
	public:
		inline uint total_num_dofs() const
		{
			return m_total_num_dofs;
		}

		inline uint num_dofs(GeometricBaseObject objType) const
		{
			return m_dof_pattern[objType];
		}
		inline uint num_dofs(uint i) const
		{
			return m_dof_pattern[i];
		}

		template <typename TElem>
		inline uint num_dofs() const
		{
			return m_dof_pattern[geometry_traits<TElem>::BASE_OBJECT_TYPE_ID];
		}

		inline void set_num_dofs(GeometricBaseObject objType, uint num)
		{
			assert(objType < NUM_GEOMETRIC_BASE_OBJECTS);
			m_dof_pattern[objType] = num;
			m_total_num_dofs += num;
		}

		template <typename TElem>
		inline void set_num_dofs(uint num)
		{
			m_dof_pattern[geometry_traits<TElem>::BASE_OBJECT_TYPE_ID] = num;
			m_total_num_dofs += num;
		}

		ElementDoFPattern& operator=(const ElementDoFPattern& item)
		{
			for(uint i = 0; i < NUM_GEOMETRIC_BASE_OBJECTS; ++i)
			{
				m_dof_pattern[i] = item.m_dof_pattern[i];
			}
			return *this;
		}
		ElementDoFPattern& operator+=(const ElementDoFPattern& item)
		{
			for(uint i = 0; i < NUM_GEOMETRIC_BASE_OBJECTS; ++i)
			{
				m_dof_pattern[i] += item.m_dof_pattern[i];
			}
			return *this;
		}
		ElementDoFPattern& operator-=(const ElementDoFPattern& item)
		{
			for(uint i = 0; i < NUM_GEOMETRIC_BASE_OBJECTS; ++i)
			{
				m_dof_pattern[i] -= item.m_dof_pattern[i];
			}
			return *this;
		}

		ElementDoFPattern(const ElementDoFPattern& item)
		{
			for(uint i = 0; i < NUM_GEOMETRIC_BASE_OBJECTS; ++i)
			{
				m_dof_pattern[i] = item.m_dof_pattern[i];
			}
		}
		ElementDoFPattern()
		{
			for(uint i = 0; i < NUM_GEOMETRIC_BASE_OBJECTS; ++i)
			{
				m_dof_pattern[i] = 0;
			}
			m_total_num_dofs = 0;
		}

		friend bool operator==(const ElementDoFPattern& lhs, const ElementDoFPattern& rhs);

	private:
		uint m_dof_pattern[NUM_GEOMETRIC_BASE_OBJECTS];
		uint m_total_num_dofs;
};

bool operator==(const ElementDoFPattern& lhs, const ElementDoFPattern& rhs);


enum TrialSpaceType {
	TST_INVALID = 0,
	TST_P1CONFORM
};

/*
 * Trial Function Interface
 *
 * This class can evaluate all trial functions for a given Geometric Object.
 * The use of terms is as follows:
 * We call a 'finite element' a geometric entity, that contains dofs. Specially dofs
 * can be located at the boundary of the geometric object, e.g. in vertices, edges on a triangle.
 * In contrast a geometric object is a grid entity. It can contain dofs, induced by the finite element it
 * belongs to.
 */
template <typename TElem>
class TrialFunctions
{
public:
	typedef TElem geom_obj_type;
	static const std::size_t RefDim = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	typedef MathVector<RefDim> position_type;
	typedef number shape_value_type;
	typedef MathVector<RefDim> grad_value_type;

public:
	// number of dofs on finite element
	virtual uint num_dofs() const = 0;

	// evaluate shape function
	virtual bool evaluate_shape(int nrShapeFct, const position_type& locPos, shape_value_type& value) const = 0;

	// evaluate shape function
	virtual bool evaluate_shape(int nrShapeFct, position_type[], shape_value_type values[], int n) const = 0;

	// evaluate gradient
	virtual bool evaluate_shape_grad(int nrShapeFct, const position_type& locPos, grad_value_type& value) const = 0;

	// local position of dof, returns true if exists, returns false if no meaningful position possible
	virtual bool position_of_dof(int nrShapeFct, position_type& value) const = 0;

	// returns true, if the baseFunctionSpace is p-adaptive
	virtual bool is_adaptive() const = 0;

	// number of dofs per geom object
	virtual const ElementDoFPattern& element_dof_pattern() const = 0;

	// virtual destructor
	virtual ~TrialFunctions()
	{};
};


template <typename TElem>
class P1conform : public TrialFunctions<TElem>{
protected:
		static const std::size_t RefDim = TrialFunctions<TElem>::RefDim;
		static const std::size_t nsh = reference_element_traits<TElem>::NumberCorners;

public:
		bool evaluate_shape(int nrShapeFct, const MathVector< RefDim >& locPos, number& value) const;
		bool evaluate_shape(int nrShapeFct, MathVector< RefDim > locPos[], number values[], int n) const;
		bool evaluate_shape_grad(int nrShapeFct, const MathVector< RefDim >& locPos, MathVector< RefDim >& value) const;
		bool position_of_dof(int nrShapeFct, MathVector< RefDim >& value) const;
		bool is_adaptive() const { return false;}
		uint num_dofs() const { return nsh;	}
		const ElementDoFPattern& element_dof_pattern() const
		{
			return m_ElementDoFPattern;
		}

		static const P1conform<TElem>& inst()
		{
			static P1conform<TElem> instance;
			return instance;
		}

	private:
		// disallow constructor
		P1conform()
		{
			m_ElementDoFPattern.set_num_dofs<Vertex>(1);
		};
		P1conform(const P1conform&);
		P1conform& operator=(const P1conform&);
		~P1conform(){};

		static const uint _order = 1;
		ElementDoFPattern m_ElementDoFPattern;
};

// Singleton, holding all Trial Spaces available
class TrialSpaces {

		// private constructor
		TrialSpaces(){};
		TrialSpaces(const TrialSpaces&);
		TrialSpaces& operator=(const TrialSpaces&);
		~TrialSpaces(){};

		template <typename TElem>
		static inline TrialFunctions<TElem>& get_TrialFunctions(TrialSpaceType type)
		{
			typedef typename ug::TrialFunctions<TElem> trial_functions_type;
			typedef std::map<TrialSpaceType, const trial_functions_type* > CollectionMap;
			CollectionMap& collection_map = get_collection_map<TElem>();

			typename CollectionMap::const_iterator iter = collection_map.find(type);

			if(iter == collection_map.end())
			{
				assert(0 && "Unknown Trial Space Type.\n");
			}

			return *(iter->second);
		}

		template <typename TElem>
		static inline std::map<TrialSpaceType, const TrialFunctions<TElem>* >& get_collection_map()
		{
			typedef typename ug::TrialFunctions<TElem> trial_functions_type;
			typedef std::map<TrialSpaceType, const trial_functions_type* > CollectionMap;
			static CollectionMap m_Collection;

			return m_Collection;
		}

		// remember if a TrialSpace is adaptive
		static std::map<TrialSpaceType, bool> m_bIsAdaptive;
		static std::map<TrialSpaceType, std::map<int, ElementDoFPattern> > m_dimDoFPattern;

	public:
		static TrialSpaces& inst()
		{
			static TrialSpaces myInst;
			return myInst;
		};

		template <typename TElem>
		static bool register_trial_space(TrialSpaceType type, const TrialFunctions<TElem>& collection)
		{
			std::cout << "register_trial_space.\n" << std::flush;

			// remember if trial space remains adaptive or not
			std::map<TrialSpaceType, bool>::iterator iter_adaptive = m_bIsAdaptive.find(type);
			if(iter_adaptive == m_bIsAdaptive.end())
			{
				bool is_adaptive = collection.is_adaptive();
				std::pair<TrialSpaceType, bool> myPair(type, is_adaptive);
				m_bIsAdaptive.insert(myPair);
			}
			else
			{
				iter_adaptive->second = iter_adaptive->second && collection.is_adaptive();
			}

			// remember ElementDoFPattern
			int dim = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
			std::map<TrialSpaceType, std::map<int, ElementDoFPattern> >::iterator iter_pattern = m_dimDoFPattern.find(type);
			if(iter_pattern == m_dimDoFPattern.end())
			{
				std::map<int, ElementDoFPattern> v;
				v.insert(std::pair<int, ElementDoFPattern>(dim, collection.element_dof_pattern()));
				m_dimDoFPattern.insert(std::pair<TrialSpaceType, std::map<int, ElementDoFPattern> >(type, v));
			}
			else
			{
				std::map<int, ElementDoFPattern>::iterator it = (iter_pattern->second).find(dim);
				if(it == (iter_pattern->second).end())
				{
					(iter_pattern->second).insert(std::pair<int, ElementDoFPattern>(dim, collection.element_dof_pattern()));
				}
				else
				{
					for(uint i = 0; i < dim; ++i)
					{
						if((it->second).num_dofs(i) != collection.element_dof_pattern().num_dofs(i))
							assert(0 && "Trial Space for different element does not match.");
					}
				}
			}

			typedef typename ug::TrialFunctions<TElem> trial_functions_type;
			typedef std::map<TrialSpaceType, const trial_functions_type* > CollectionMap;
			CollectionMap& collection_map = get_collection_map<TElem>();

			return collection_map.insert(std::pair<TrialSpaceType, const trial_functions_type*>(type, &collection)).second;
		}

		template <typename TElem>
		static bool unregister_trial_space(TrialSpaceType type)
		{
			typedef const typename ug::TrialFunctions<TElem> trial_functions_type;
			typedef std::map<TrialSpaceType, const trial_functions_type* > CollectionMap;
			CollectionMap& collection_map = get_collection_map<TElem>();

			return collection_map.erase(type) == 1;
		}

		template <typename TElem>
		static const TrialFunctions<TElem>& TrialFunctions(TrialSpaceType type)
		{
			return inst().get_TrialFunctions<TElem>(type);
		}

		static const ElementDoFPattern& get_element_dof_pattern(TrialSpaceType type, uint dim)
		{
			std::map<TrialSpaceType, std::map<int, ElementDoFPattern> >::iterator iter_pattern = m_dimDoFPattern.find(type);
			if(iter_pattern == m_dimDoFPattern.end())
			{
				assert(0 && "Trial Space not found.\n");
			}

			std::map<int, ElementDoFPattern>::iterator it = (iter_pattern->second).find(dim);
			if(it == (iter_pattern->second).end())
			{
				assert(0 && "Dimension not found.\n");
			}

			return it->second;
		}
};

} /* end namespace libDiscretization */


#endif /* __H__LIBDISCRETIZATION__TRIALSPACE__ */
