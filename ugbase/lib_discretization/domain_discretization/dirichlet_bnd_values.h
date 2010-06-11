/*
 * dirichlet_bnd_values.h
 *
 *  Created on: 08.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DIRICHLET_BND_VALUES__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DIRICHLET_BND_VALUES__

namespace ug{

enum BND_TYPE {
	BND_TYPE_NONE = 0,
	BND_TYPE_DIRICHLET,
	BND_TYPE_NEUMANN
};

template <typename TDiscreteFunction>
class DirichletBNDValues {
	public:
		// type of discrete function
		typedef TDiscreteFunction discrete_function_type;

		// type of domain
		typedef typename TDiscreteFunction::domain_type domain_type;

		// type of position coordinates (e.g. position_type)
		typedef typename domain_type::position_type position_type;

		// type of algebra
		typedef typename TDiscreteFunction::algebra_type algebra_type;

		// type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

		// type of local matrix
		typedef typename matrix_type::local_matrix_type local_matrix_type;

		// type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

		// type of algebra vector
		typedef typename vector_type::local_vector_type local_vector_type;

		// type of multi_index used in algebra
		typedef typename matrix_type::index_type index_type;

		// type of local index container
		typedef typename matrix_type::local_index_type local_index_type;

	public:
		DirichletBNDValues(domain_type& domain) :
		  m_domain(domain)
		{
			typename domain_type::subset_handler_type& sh = domain.get_subset_handler();
			int num_sh = sh.num_subsets();
			// TODO: Handle this
			size_t num_fct = 10;
			m_bndfct.resize(num_fct);
			m_bndtype.resize(num_fct);
			for(size_t fct = 0; fct < num_fct; ++fct)
			{
				m_bndfct[fct].resize(num_sh, NULL);
				m_bndtype[fct].resize(num_sh, BND_TYPE_NONE);
			}
		}

		bool clear_dirichlet_jacobian_defect(	 	matrix_type& J, vector_type& d, const discrete_function_type& u, number time = 0.0);
		bool clear_dirichlet_jacobian(				matrix_type& J, const discrete_function_type& u, number time = 0.0);
		bool clear_dirichlet_defect(				vector_type& d, const discrete_function_type& u, number time = 0.0);
		bool set_dirichlet_linear(					matrix_type& mat, vector_type& rhs, const discrete_function_type& u, number time = 0.0);

		bool set_dirichlet_solution(				discrete_function_type& u, number time = 0.0);

	protected:
		template <typename TElem>
		bool clear_dirichlet_jacobian_defect(	typename geometry_traits<TElem>::iterator iterBegin,
												typename geometry_traits<TElem>::iterator iterEnd,
												size_t fct, int si, matrix_type& J, vector_type& d,
												const discrete_function_type& u, number time = 0.0);

		template <typename TElem>
		bool clear_dirichlet_jacobian(			typename geometry_traits<TElem>::iterator iterBegin,
												typename geometry_traits<TElem>::iterator iterEnd,
												size_t fct, int si, matrix_type& J,
												const discrete_function_type& u, number time = 0.0);

		template <typename TElem>
		bool clear_dirichlet_defect(			typename geometry_traits<TElem>::iterator iterBegin,
												typename geometry_traits<TElem>::iterator iterEnd,
												size_t fct, int si, vector_type& d,
												const discrete_function_type& u, number time = 0.0);

		template <typename TElem>
		bool set_dirichlet_solution( 			typename geometry_traits<TElem>::iterator iterBegin,
												typename geometry_traits<TElem>::iterator iterEnd,
												size_t fct, int si, vector_type& x,
												const discrete_function_type& u, number time = 0.0);

		template <typename TElem>
		bool set_dirichlet_linear(				typename geometry_traits<TElem>::iterator iterBegin,
												typename geometry_traits<TElem>::iterator iterEnd,
												size_t fct, int si, matrix_type& mat, vector_type& rhs,
												const discrete_function_type& u, number time = 0.0);

	protected:
		typedef bool (*Boundary_fct)(number&, const position_type&, number);

		struct DirichletSubset
		{
			int si;
			size_t fct;
			Boundary_fct func;
		};

	public:
		inline bool boundary_value(number& value, const position_type& pos, size_t fct, int si, number time)
		{
			return (m_bndfct[fct][si])(value, pos, time);
		}

		// add bndtype and bndfunction to subset s for function fct
		bool add_boundary_value(size_t d, int si, size_t fct, Boundary_fct func, BND_TYPE type)
		{
			std::vector<int>::iterator iter = find(m_bnd_subset[d].begin(), m_bnd_subset[d].end(), si);
			if(iter == m_bnd_subset[d].end())
				m_bnd_subset[d].push_back(si);

			m_bndtype[fct][si] = type;
			m_bndfct[fct][si] = func;
			return true;
		}

		bool is_dirichlet(int si, size_t fct) {return m_bndtype[fct][si] == BND_TYPE_DIRICHLET;}

		size_t num_bnd_subsets(size_t d) {return m_bnd_subset[d].size();}
		int bnd_subset(size_t d, size_t i) {return m_bnd_subset[d][i];}

	protected:
		domain_type& m_domain;

		std::vector<int> m_bnd_subset[domain_type::dim];

		std::vector<std::vector<BND_TYPE> > m_bndtype;
		std::vector<std::vector<Boundary_fct> > m_bndfct;

};


template <typename TDiscreteFunction>
bool
DirichletBNDValues<TDiscreteFunction>::
clear_dirichlet_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u, number time)
{
	for(size_t i = 0; i < num_bnd_subsets(0); ++i)
	{
		int si = bnd_subset(0, i);
		for(size_t fct = 0; fct < u.num_fct(); ++fct)
		{
			if(!is_dirichlet(si, fct)) continue;
			if(clear_dirichlet_jacobian_defect<VertexBase>(u.template begin<VertexBase>(si),u.template end<VertexBase>(si), fct, si, J, d, u, time) == false)
			{
				UG_LOG("Error in set_dirichlet_nodes, aborting.\n");
				return false;
			}
		}
	}
	return true;
}


template <typename TDiscreteFunction>
bool
DirichletBNDValues<TDiscreteFunction>::
clear_dirichlet_jacobian(matrix_type& J, const discrete_function_type& u, number time)
{
	for(size_t i = 0; i < num_bnd_subsets(0); ++i)
	{
		int si = bnd_subset(0, i);
		for(size_t fct = 0; fct < u.num_fct(); ++fct)
		{
			if(!is_dirichlet(si, fct)) continue;
			if(clear_dirichlet_jacobian<VertexBase>(u.template begin<VertexBase>(si),u.template end<VertexBase>(si), fct, si, J, u, time) == false)
			{
				UG_LOG("Error in 'assemble_jacobian' while calling 'clear_dirichlet_jacobian', aborting.\n");
				return false;
			}
		}
	}
	return true;
}

template <typename TDiscreteFunction>
bool
DirichletBNDValues<TDiscreteFunction>::
clear_dirichlet_defect(vector_type& d, const discrete_function_type& u, number time)
{
	for(size_t i = 0; i < num_bnd_subsets(0); ++i)
	{
		int si = bnd_subset(0, i);
		for(size_t fct = 0; fct < u.num_fct(); ++fct)
		{
			if(!is_dirichlet(si, fct)) continue;
			if(clear_dirichlet_defect<VertexBase>(u.template begin<VertexBase>(si),u.template end<VertexBase>(si), fct, si, d, u, time) == false)
			{
				UG_LOG("Error in set_dirichlet_nodes, aborting.\n");
				return false;
			}
		}
	}
	return true;
}

template <typename TDiscreteFunction>
bool
DirichletBNDValues<TDiscreteFunction>::
set_dirichlet_solution(discrete_function_type& u, number time)
{
	for(size_t i = 0; i < num_bnd_subsets(0); ++i)
	{
		int si = bnd_subset(0, i);
		for(size_t fct = 0; fct < u.num_fct(); ++fct)
		{
			if(!is_dirichlet(si, fct)) continue;
			if(set_dirichlet_solution<VertexBase>(u.template begin<VertexBase>(si), u.template end<VertexBase>(si), fct, si, u.get_vector(), u, time) == false)
			{
				UG_LOG("Error in set_dirichlet_nodes, aborting.\n");
				return false;
			}
		}
	}
	return true;
}

template <typename TDiscreteFunction>
bool
DirichletBNDValues<TDiscreteFunction>::
set_dirichlet_linear(		matrix_type& mat, vector_type& rhs, const discrete_function_type& u, number time)
{
	for(size_t i = 0; i < num_bnd_subsets(0); ++i)
	{
		int si = bnd_subset(0, i);
		for(size_t fct = 0; fct < u.num_fct(); ++fct)
		{
			if(!is_dirichlet(si, fct)) continue;
			if(set_dirichlet_linear<VertexBase>(u.template begin<VertexBase>(si), u.template end<VertexBase>(si), fct, si, mat, rhs, u) == false)
			{
				UG_LOG("Error in assemble_linear, aborting.\n");
				return false;
			}
		}
	}
	return true;
}


template <typename TDiscreteFunction>
template <typename TElem>
bool
DirichletBNDValues<TDiscreteFunction>::
clear_dirichlet_jacobian_defect(	typename geometry_traits<TElem>::iterator iterBegin,
									typename geometry_traits<TElem>::iterator iterEnd,
									size_t fct, int si,
									matrix_type& J, vector_type& d,
									const discrete_function_type& u, number time)
{
	local_index_type ind(1);
	local_index_type glob_ind;
	local_vector_type dirichlet_vals;

	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		TElem* elem = *iter;

		if(u.get_multi_indices_of_geom_obj(elem, fct, ind) != 1) assert(0);

		glob_ind.push_back(ind[0]);
		dirichlet_vals.push_back(0.0);
	}

	if(d.set(dirichlet_vals, glob_ind) != true)
		return false;

	if(J.set_dirichlet_rows(glob_ind) != true)
		return false;

	return true;
}

template <typename TDiscreteFunction>
template <typename TElem>
bool
DirichletBNDValues<TDiscreteFunction>::
clear_dirichlet_defect(	typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd,
						size_t fct, int si,
						vector_type& d,
						const discrete_function_type& u, number time)
{
	local_index_type ind(1);
	local_index_type glob_ind;
	local_vector_type dirichlet_vals;

	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		TElem* elem = *iter;

		if(u.get_multi_indices_of_geom_obj(elem, fct, ind) != 1) assert(0);

		glob_ind.push_back(ind[0]);
		dirichlet_vals.push_back(0.0);
	}

	if(d.set(dirichlet_vals, glob_ind) != true)
		return false;

	return true;
}

template <typename TDiscreteFunction>
template <typename TElem>
bool
DirichletBNDValues<TDiscreteFunction>::
clear_dirichlet_jacobian(	typename geometry_traits<TElem>::iterator iterBegin,
							typename geometry_traits<TElem>::iterator iterEnd,
							size_t fct, int si,
							matrix_type& J,
							const discrete_function_type& u, number time)
{
	local_index_type ind(1);
	local_index_type glob_ind;

	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		TElem* elem = *iter;

		if(u.get_multi_indices_of_geom_obj(elem, fct, ind) != 1) assert(0);

		glob_ind.push_back(ind[0]);
	}

	if(J.set_dirichlet_rows(glob_ind) != true)
		return false;

	return true;
}

template <typename TDiscreteFunction>
template <typename TElem>
bool
DirichletBNDValues<TDiscreteFunction>::
set_dirichlet_solution(	typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd,
						size_t fct, int si,
						vector_type& x,
						const discrete_function_type& u, number time)
{
	typename domain_type::position_accessor_type aaPos = u.get_domain().get_position_accessor();

	local_index_type ind(1);
	local_index_type glob_ind;
	local_vector_type dirichlet_vals;

	number val;
	position_type corner;

	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		TElem* elem = *iter;
		// TODO: if TElem != Vertex we have to do something else
		corner = aaPos[elem];

		if(u.get_multi_indices_of_geom_obj(elem, fct, ind) != 1) assert(0);

		if(boundary_value(val, corner, fct, si, time))
		{
			glob_ind.push_back(ind[0]);
			dirichlet_vals.push_back(val);
		}
	}

	if(x.set(dirichlet_vals, glob_ind) != true)
		return false;

	return true;
}


template <typename TDiscreteFunction>
template <typename TElem>
bool
DirichletBNDValues<TDiscreteFunction>::
set_dirichlet_linear(	typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd,
						size_t fct, int si,
						matrix_type& mat, vector_type& rhs,
						const discrete_function_type& u, number time)
{
	typename domain_type::position_accessor_type aaPos = u.get_domain().get_position_accessor();

	local_index_type ind(1);
	local_index_type glob_ind;
	local_vector_type dirichlet_vals;

	number val;
	position_type corner;

	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		TElem* elem = *iter;
		// TODO: if TElem != Vertex we have to do something else
		corner = aaPos[elem];

		if(u.template get_multi_indices_of_geom_obj<TElem>(elem, fct, ind) != 1) assert(0);
		if(boundary_value(val, corner, fct, si, time))
		{
			glob_ind.push_back(ind[0]);
			dirichlet_vals.push_back(val);
		}
	}

	if(rhs.set(dirichlet_vals, glob_ind) != true)
		return false;

	if(mat.set_dirichlet_rows(glob_ind) != true)
		return false;

	return true;
}




} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DIRICHLET_BND_VALUES__ */
