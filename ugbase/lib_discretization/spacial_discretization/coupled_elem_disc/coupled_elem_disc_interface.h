/*
 * coupled_elem_disc_interface.h
 *
 *  Created on: 25.06.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_ELEM_DISC_INTERFACE__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_ELEM_DISC_INTERFACE__

#include <vector>
#include <string>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

#include "../elem_disc/elem_disc_interface.h"
#include "./elem_data/data_items.h"
#include "./elem_data/data_container.h"
#include "lib_discretization/common/local_algebra.h"

namespace ug {

template <typename TAlgebra>
class ICoupledElemDisc : public IElemDisc<TAlgebra> {
	public:
		// algebra type
		typedef TAlgebra algebra_type;

		// local matrix type
		typedef LocalMatrix<typename TAlgebra::matrix_type::entry_type> local_matrix_type;

		// local vector type
		typedef LocalVector<typename TAlgebra::vector_type::entry_type> local_vector_type;

		// local index type
		typedef LocalIndices local_index_type;

	public:
		virtual size_t num_imports() = 0;

		virtual DataImportItem* import(size_t i) = 0;

		virtual bool register_exports(DataContainer& Cont) = 0;

		virtual bool unregister_exports(DataContainer& Cont) = 0;

		virtual bool register_imports(DataContainer& Cont) = 0;

		virtual bool unregister_imports(DataContainer& Cont) = 0;

		virtual bool set_sys_id(size_t sys_id) = 0;

	public:
		template <typename data_type, typename position_type>
		void data_export(	int nr, std::vector<data_type>& val, std::vector<std::vector<data_type> >& deriv,
							const std::vector<position_type>& pos, const local_vector_type& u, bool compDeriv)
		{
			using std::vector;
			typedef void (ICoupledElemDisc<TAlgebra>::*ExportFunc)(vector<data_type>& val, vector<vector<data_type> >& deriv,
															const vector<position_type>& pos, const local_vector_type& u, bool compDeriv);
			vector<vector<ExportFunc> >& vDataExport = get_vDataExport<data_type,position_type>();

			return (this->*(vDataExport[nr][IElemDisc<TAlgebra>::m_id]))(val, deriv, pos, u, compDeriv);
		};

	protected:

		template <typename data_type, typename position_type, typename TFunc>
		bool register_data_export_function(int id, size_t nr, TFunc func)
		{
			using std::vector;
			typedef void (ICoupledElemDisc<TAlgebra>::*ExportFunc)(vector<data_type>&, vector<vector<data_type> >&,
															const vector<position_type>&, const local_vector_type&, bool);
			vector<vector<ExportFunc> >& vDataExport = get_vDataExport<data_type,position_type>();

			if(nr < vDataExport.size() && (size_t)id < vDataExport[nr].size() && vDataExport[nr][id] != 0)
			{
				UG_LOG("Trying to register export function for id "<< id << " and nr " << nr << "twice.\n");
				return false;
			}

			if((size_t)nr >= vDataExport.size())
				vDataExport.resize(nr+1);

			if((size_t)id >= vDataExport[nr].size())
				vDataExport[nr].resize(id+1, 0);

			vDataExport[nr][id] = (ExportFunc)func;

			return 0;
		}

		template <typename data_type, typename position_type>
		bool data_export_function_registered(int id, size_t nr)
		{
			using std::vector;
			typedef void (ICoupledElemDisc<TAlgebra>::*ExportFunc)(vector<data_type>&, vector<vector<data_type> >&,
															const vector<position_type>&, const local_vector_type&, bool);
			vector<vector<ExportFunc> >& vDataExport = get_vDataExport<data_type,position_type>();

			if(nr < vDataExport.size())
			{
				if(id >= 0 && (size_t)id < vDataExport[nr].size())
				{
					if(vDataExport[nr][id] != 0)
						return true;
				}
			}
			return false;
		}

		template <typename data_type, typename position_type>
		std::vector<std::vector<void (ICoupledElemDisc<TAlgebra>::*)(	std::vector<data_type>&, std::vector<std::vector<data_type> >&,
													const std::vector<position_type>&, const local_vector_type&, bool)> >& get_vDataExport()
		{
			using std::vector;
			typedef void (ICoupledElemDisc<TAlgebra>::*ExportFunc)(vector<data_type>&, vector<vector<data_type> >&,
															const vector<position_type>&, const local_vector_type&, bool);
			static vector<vector<ExportFunc> > vDataExport;
			return vDataExport;
		};
};


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_ELEM_DISC_INTERFACE__ */
