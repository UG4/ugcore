/*
 * amg_misc.h
 *
 *  Created on: 19.12.2011
 *      Author: mrupp
 */

#ifndef AMG_MISC_H_
#define AMG_MISC_H_

#include "rsamg/rsamg_impl.h"
namespace ug
{

template<typename TMatrix>
void ReduceToStrongConnections(TMatrix &m1, const TMatrix &const_m2, double theta)
{
	PROFILE_FUNC();
	UG_LOG("Reducing matrix...");
	TMatrix &m2 = const_cast<TMatrix&>(const_m2);
	std::vector<typename TMatrix::connection> con;
	m1.resize(m2.num_rows(), m2.num_cols());

#ifdef UG_PARALLEL
	m1.set_layouts(m2.master_layout(), m2.slave_layout());
	m1.set_communicator(m2.communicator());
	m1.set_process_communicator(m2.process_communicator());
	m1.copy_storage_type(m2);
#endif

	con.reserve(50);
	for(size_t i=0; i<m1.num_rows(); i++)
	{
		con.clear();
		double minConnValue, maxConnValue, diag;
		GetNeighborValues(m2, i, minConnValue, maxConnValue, diag);
		typename TMatrix::connection c;
		double barrier = theta*std::max(maxConnValue, dabs(minConnValue));
		double offLumpedSum=0;
		int iDiag=-1;
		for(typename TMatrix::row_iterator conn = m2.begin_row(i); conn != m2.end_row(i); ++conn)
		{
			if(conn.index() != i)
			{
				if(dabs(amg_offdiag_value(conn.value())) < barrier) // only strong
				{
					offLumpedSum += conn.value();
					continue;
				}
			}
			else
				iDiag = con.size();
			c.dValue = conn.value();
			c.iIndex = conn.index();
			con.push_back(c);
		}
		UG_ASSERT(iDiag != -1, "");
		con[iDiag].dValue += offLumpedSum;
		m1.set_matrix_row(i, &con[0], con.size());
	}
	UG_LOG("done.\n");
}

}

#endif // AMG_MISC_H_
