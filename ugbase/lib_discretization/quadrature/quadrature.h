/*
 * quadrature.h
 *
 *  Created on: 15.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__QUADRATURE__
#define __H__LIBDISCRETIZATION__QUADRATURE__

#include "../reference_element/reference_element.h"

namespace ug{

template <typename TRefElem>
class QuadratureRule{
	public:
		static const int dim = TRefElem::dim;
		typedef MathVector<dim> position_type;
		typedef number weight_type;

	public:
		inline std::size_t num_points() const
		{
			return m_num_points;
		}

		inline position_type point(std::size_t i) const
		{
			assert(i < m_num_points && i >= 0);
			return m_points[i];
		}

		inline position_type* points() const
		{
			return m_points;
		}

		inline weight_type weight(std::size_t i) const
		{
			assert(i < m_num_points && i >= 0);
			return m_weights[i];
		}

		inline weight_type* weights() const
		{
			return m_weights;
		}

		inline int order() const
		{
			return m_order;
		}

	protected:
		position_type* m_points;
		weight_type* m_weights;
		std::size_t m_num_points;
		int m_order;
};

template <typename TRefElem>
class QuadratureRuleFactory;


template <typename TRefElem>
bool RegisterQuadratureRule(QuadratureRuleFactory<TRefElem>& factory);

template <typename TRefElem>
class QuadratureRuleFactory{
	private:
		QuadratureRuleFactory()
		{
			RegisterQuadratureRule<TRefElem>(*this);
		}
		QuadratureRuleFactory(const QuadratureRuleFactory&){};
		QuadratureRuleFactory& operator=(const QuadratureRuleFactory&);

		static const QuadratureRule<TRefElem>& get_rule(int order)
		{
			UG_ASSERT(order < (int)m_rules.size(), "Rule does not exist.");
			if(m_rules[order] == 0)
			{
				for(size_t i = order +1; i < m_rules.size(); ++i)
				{
					if(m_rules[i] != 0) return *m_rules[i];
				}
				UG_ASSERT(0, "Cannot find rule.");
			}
			return *m_rules[order];
		}

		static std::vector<const QuadratureRule<TRefElem>*> m_rules;
		static bool m_initialized;

	public:
		static QuadratureRuleFactory<TRefElem>& instance()
		{
			static QuadratureRuleFactory<TRefElem> inst;
			return inst;
		}

		static bool register_rule(const QuadratureRule<TRefElem>& rule)
		{
			int order = rule.order();
			if((int) m_rules.size() <= order) m_rules.resize(order+1, 0);
			if(m_rules[order] != 0) assert(0 && "Already Quadrature Rule for this order registered.");

			m_rules[order] = &rule;
			return true;
		}

		inline static const QuadratureRule<TRefElem>& get_quadrature_rule(int order)
		{
			return instance().get_rule(order);
		}
};

template <typename TRefElem>
class GaussQuadrature : public QuadratureRule<TRefElem>{
	public:
		GaussQuadrature(int order);
		~GaussQuadrature();
	protected:
		inline bool allocate_memory(std::size_t n);
};


} // namespace ug

#include "quadrature_impl.h"

#endif /* __H__LIBDISCRETIZATION__QUADRATURE__ */
