/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__COMMON__LOCAL_ALGEBRA__
#define __H__UG__LIB_DISC__COMMON__LOCAL_ALGEBRA__

//#define UG_LOCALALGEBRA_ASSERT(cond, exp)
// include define below to assert arrays used in stabilization
#define UG_LOCALALGEBRA_ASSERT(cond, exp) UG_ASSERT((cond), exp)

#include <vector>

#include "./multi_index.h"
#include "./function_group.h"
#include "lib_algebra/small_algebra/small_algebra.h"

namespace ug{


class LocalIndices
{
	public:
	///	Index type used by algebra
		using index_type = size_t;

	///	Component type used by algebra
		using comp_type = index_type;

	public:
	///	Default Constructor
		LocalIndices() {};

	///	sets the number of functions
		void resize_fct(size_t numFct)
		{
			m_vvIndex.resize(numFct);
			m_vLFEID.resize(numFct);
		}

	///	sets the local finite element id for a function
		void set_lfeID(size_t fct, const LFEID& lfeID)
		{
			UG_ASSERT(fct < m_vLFEID.size(), "Invalid index: "<<fct);
			m_vLFEID[fct] = lfeID;
		}

	///	returns the local finite element id of a function
		const LFEID& local_finite_element_id(size_t fct) const
		{
			UG_ASSERT(fct < m_vLFEID.size(), "Invalid index: "<<fct);
			return m_vLFEID[fct];
		}

	///	sets the number of dofs of a function
		void resize_dof(size_t fct, size_t numDoF)
		{
			check_fct(fct);
			m_vvIndex[fct].resize(numDoF);
		}

	///	clears the dofs of a function
		void clear_dof(size_t fct) {resize_dof(fct, 0);}

	/// reserves memory for the number of dofs
		void reserve_dof(size_t fct, size_t numDoF)
		{
			check_fct(fct);
			m_vvIndex[fct].reserve(numDoF);
		}

	///	adds an index (increases size)
		void push_back_index(size_t fct, size_t index) {push_back_multi_index(fct,index,0);}

	///	adds an index (increases size)
		void push_back_multi_index(size_t fct, size_t index, size_t comp)
		{
			check_fct(fct);
			m_vvIndex[fct].push_back(DoFIndex(index,comp));
		}

	///	clears all fct
		void clear() {m_vvIndex.clear();}

	///	number of functions
		size_t num_fct() const {return m_vvIndex.size();}

	/// number of dofs for accessible function
		size_t num_dof(size_t fct) const
		{
			check_fct(fct);
			return m_vvIndex[fct].size();
		}

	/// number of dofs of all accessible (sum)
		size_t num_dof() const
		{
			size_t num = 0;
			for(size_t fct = 0; fct < num_fct(); ++fct) num += num_dof(fct);
			return num;
		}

	/// global algebra multi-index for (fct, dof)
		const DoFIndex& multi_index(size_t fct, size_t dof) const
		{
			check_dof(fct, dof);
			return m_vvIndex[fct][dof];
		}

	/// global algebra index for (fct, dof)
		index_type index(size_t fct, size_t dof) const
		{
			check_dof(fct, dof);
			return m_vvIndex[fct][dof][0];
		}

	/// global algebra index for (fct, dof)
		index_type& index(size_t fct, size_t dof)
		{
			check_dof(fct, dof);
			return m_vvIndex[fct][dof][0];
		}

	/// algebra comp for (fct, dof)
		comp_type comp(size_t fct, size_t dof) const
		{
			check_dof(fct, dof);
			return m_vvIndex[fct][dof][1];
		}

	/// algebra comp for (fct, dof)
		comp_type& comp(size_t fct, size_t dof)
		{
			check_dof(fct, dof);
			return m_vvIndex[fct][dof][1];
		}
	
	///	checks if the local index object references a given index
		bool contains_index(index_type ind)
		{
			for(size_t fct = 0; fct < num_fct(); fct++)
				for(size_t dof = 0; dof < num_dof(fct); dof++)
					if(m_vvIndex[fct][dof][0] == ind)
						return true;
			return false;
		}

	protected:
	///	checks correct fct index in debug mode
		inline void check_fct(size_t fct) const
		{
			UG_LOCALALGEBRA_ASSERT(fct < num_fct(), "Wrong index.");
		}
	///	checks correct dof index in debug mode
		inline void check_dof(size_t fct, size_t dof) const
		{
			check_fct(fct);
			UG_LOCALALGEBRA_ASSERT(dof < num_dof(fct), "Wrong index.");
		}

	protected:
	// 	Mapping (fct, dof) -> local index
		std::vector<std::vector<DoFIndex> > m_vvIndex;

	//	Local finite element ids
		std::vector<LFEID> m_vLFEID;
};

class LocalVector
{
	public:
	///	own type
		using this_type = LocalVector;

	///	Type to store DoF values
		using value_type = number;

	public:
	///	default Constructor
		LocalVector() : m_pIndex(nullptr) {m_vvValue.clear();}

	///	Constructor
		LocalVector(const LocalIndices& ind) {resize(ind);}

	///	resize for current local indices
		void resize(const LocalIndices& ind)
		{
			m_pIndex = &ind;
			m_vvValue.resize(ind.num_fct());
			m_vvValueAcc.resize(m_vvValue.size());
			for(size_t fct = 0; fct < m_vvValue.size(); ++fct)
				m_vvValue[fct].resize(ind.num_dof(fct));
			access_all();
		}

	/// get current local indices
		const LocalIndices& get_indices() const {return *m_pIndex;}

		///////////////////////////
		// vector functions
		///////////////////////////

	this_type& operator=(const this_type& other)
	{
		m_pIndex = other.m_pIndex;
		const size_t numFcts = m_pIndex->num_fct();
		m_vvValue.resize(numFcts);
		for (size_t fct = 0; fct < numFcts; ++fct)
			m_vvValue[fct].resize(m_pIndex->num_dof(fct));
		m_vvValueAcc.resize(numFcts);
		if (other.m_pFuncMap)
			access_by_map(*other.m_pFuncMap);
		else
			access_all();

		return *this;
	}

	/// set all components of the vector
		this_type& operator=(number val)
		{
			for(size_t fct = 0; fct < m_vvValue.size(); ++fct)
				for(size_t dof = 0; dof < m_vvValue[fct].size(); ++dof)
					m_vvValue[fct][dof] = val;
			return *this;
		}

	/// multiply all components of the vector
		this_type& operator*(number val)
		{
			return this->operator*=(val);
		}

	/// multiply all components of the vector
		this_type& operator*=(number val)
		{
			for(size_t fct = 0; fct < m_vvValue.size(); ++fct)
				for(size_t dof = 0; dof < m_vvValue[fct].size(); ++dof)
					m_vvValue[fct][dof] *= val;
			return *this;
		}

	/// add a local vector
		this_type& operator+=(const this_type& rhs)
		{
			UG_LOCALALGEBRA_ASSERT(m_pIndex==rhs.m_pIndex, "Not same indices.");
			for(size_t fct = 0; fct < m_vvValue.size(); ++fct)
				for(size_t dof = 0; dof < m_vvValue[fct].size(); ++dof)
					m_vvValue[fct][dof] += rhs.m_vvValue[fct][dof];
			return *this;
		}

	/// subtract a local vector
		this_type& operator-=(const this_type& rhs)
		{
			UG_LOCALALGEBRA_ASSERT(m_pIndex==rhs.m_pIndex, "Not same indices.");
			for(size_t fct = 0; fct < m_vvValue.size(); ++fct)
				for(size_t dof = 0; dof < m_vvValue[fct].size(); ++dof)
					m_vvValue[fct][dof] -= rhs.m_vvValue[fct][dof];
			return *this;
		}

	///	add a scaled vector
		this_type& scale_append(number s, const this_type& rhs)
		{
			UG_LOCALALGEBRA_ASSERT(m_pIndex==rhs.m_pIndex, "Not same indices.");
			for(size_t fct = 0; fct < m_vvValue.size(); ++fct)
				for(size_t dof = 0; dof < m_vvValue[fct].size(); ++dof)
					m_vvValue[fct][dof] += s * rhs.m_vvValue[fct][dof];
			return *this;
		}

		///////////////////////////
		// restricted DoF access
		///////////////////////////

	/// access only part of the functions using mapping (restrict functions)
		void access_by_map(const FunctionIndexMapping& funcMap)
		{
			m_pFuncMap = &funcMap;
			for(size_t i = 0; i < funcMap.num_fct(); ++i)
			{
				const size_t mapFct = funcMap[i];
				m_vvValueAcc[i] = &(m_vvValue[mapFct][0]);
			}
		}

	///	access all functions
		void access_all()
		{
			m_pFuncMap = nullptr;

			if(m_pIndex==nullptr) {m_vvValueAcc.clear(); return;}

			for(size_t i = 0; i < m_pIndex->num_fct(); ++i)
				m_vvValueAcc[i] = &(m_vvValue[i][0]);
		}

	///	returns the number of currently accessible functions
		size_t num_fct() const
		{
			if(m_pFuncMap == nullptr) return m_vvValue.size();
			return m_pFuncMap->num_fct();
		}

	///	returns the local finite element id of a function
		const LFEID& local_finite_element_id(size_t fct) const
		{
			UG_ASSERT(m_pIndex != nullptr, "No indices present");
			check_fct(fct);
			if(m_pFuncMap == nullptr) return m_pIndex->local_finite_element_id(fct);
			else return m_pIndex->local_finite_element_id((*m_pFuncMap)[fct]);
		}

	///	returns the number of dofs for the currently accessible function
		size_t num_dof(size_t fct) const
		{
			check_fct(fct);
			if(m_pFuncMap == nullptr) return m_vvValue[fct].size();
			else return m_vvValue[ (*m_pFuncMap)[fct] ].size();
		}

	/// access to dof of currently accessible function fct
		number& operator()(size_t fct, size_t dof)
		{
			check_dof(fct,dof);
			return m_vvValueAcc[fct][dof];
		}

	/// const access to dof of currently accessible function fct
		number operator()(size_t fct, size_t dof) const
		{
			check_dof(fct,dof);
			return m_vvValueAcc[fct][dof];
		}

		///////////////////////////
		// all DoF access
		///////////////////////////

	///	returns the number of all functions
		size_t num_all_fct() const {return m_vvValue.size();}

	///	returns the number of dofs for a function (unrestricted functions)
		size_t num_all_dof(size_t fct) const {check_all_fct(fct); return m_vvValue[fct].size();}

	/// access to dof of a fct (unrestricted functions)
		number& value(size_t fct, size_t dof){check_all_dof(fct,dof);return m_vvValue[fct][dof];}

	/// const access to dof of a fct (unrestricted functions)
		const number& value(size_t fct, size_t dof) const{check_all_dof(fct,dof);return m_vvValue[fct][dof];}

	protected:
	///	checks correct fct index in debug mode
		inline void check_fct(size_t fct) const
		{
			UG_LOCALALGEBRA_ASSERT(fct < num_fct(), "Wrong index: fct: "<<fct<<
									" of num_fct: "<<num_fct());
		}
	///	checks correct dof index in debug mode
		inline void check_dof(size_t fct, size_t dof) const
		{
			check_fct(fct);
			UG_LOCALALGEBRA_ASSERT(dof < num_dof(fct), "Wrong index: dof: "
					<<dof<<", num_dof: "<<num_dof(fct)<<" of fct: "<<fct);
		}
	///	checks correct fct index in debug mode
		inline void check_all_fct(size_t fct) const
		{
			UG_LOCALALGEBRA_ASSERT(fct < num_all_fct(), "Wrong index.");
		}
	///	checks correct dof index in debug mode
		inline void check_all_dof(size_t fct, size_t dof) const
		{
			check_all_fct(fct);
			UG_LOCALALGEBRA_ASSERT(dof < num_all_dof(fct), "Wrong index.");
		}

	protected:
	/// Indices
		const LocalIndices* m_pIndex;

	/// Access Mapping
		const FunctionIndexMapping* m_pFuncMap;

	/// Entries (fct, dof)
		std::vector<value_type*> m_vvValueAcc;

	/// Entries (fct, dof)
		std::vector<std::vector<value_type> > m_vvValue;
};

class LocalMatrix
{
	public:
	///	own type
	using this_type = LocalMatrix;

	///	Entry type used by algebra
	using value_type = number;

	public:
	///	Constructor
		LocalMatrix() :
			m_pRowIndex(nullptr), m_pColIndex(nullptr) ,
			m_pRowFuncMap(nullptr), m_pColFuncMap(nullptr)
		{}

	///	Constructor
		LocalMatrix(const LocalIndices& rowInd, const LocalIndices& colInd)
			: m_pRowFuncMap(nullptr), m_pColFuncMap(nullptr)
		{
			resize(rowInd, colInd);
		}

	///	resize for same ind in row and column
		void resize(const LocalIndices& ind) {resize(ind, ind);}

	///	resize for current local indices
		void resize(const LocalIndices& rowInd, const LocalIndices& colInd)
		{
			m_pRowIndex = &rowInd;
			m_pColIndex = &colInd;

			m_CplMat.resize(rowInd.num_fct(), colInd.num_fct());
			m_CplMatAcc.resize(rowInd.num_fct(), colInd.num_fct());

			for(size_t fct1 = 0; fct1 < m_CplMat.num_rows(); ++fct1)
				for(size_t fct2 = 0; fct2 < m_CplMat.num_cols(); ++fct2)
				{
					m_CplMat(fct1, fct2).resize(colInd.num_dof(fct1),
												rowInd.num_dof(fct2));
				}

			access_all();
		}

	/// get current local indices
		const LocalIndices& get_row_indices() const {return *m_pRowIndex;}

	/// get current local indices
		const LocalIndices& get_col_indices() const {return *m_pColIndex;}

		///////////////////////////
		// Matrix functions
		///////////////////////////

	/// set all entries
		this_type& operator=(number val)
		{
			for(size_t fct1=0; fct1 < m_CplMat.num_rows(); ++fct1)
				for(size_t fct2=0; fct2 < m_CplMat.num_cols(); ++fct2)
					m_CplMat(fct1,fct2) = val;
			return *this;
		}

	/// multiply all entries
		this_type& operator*(number val)
		{
			return this->operator*=(val);
		}

	/// multiply matrix
		this_type& operator*=(number val)
		{
			for(size_t fct1=0; fct1 < m_CplMat.num_rows(); ++fct1)
				for(size_t fct2=0; fct2 < m_CplMat.num_cols(); ++fct2)
					m_CplMat(fct1,fct2) *= val;
			return *this;
		}

	/// add matrix
		this_type& operator+=(const this_type& rhs)
		{
			UG_LOCALALGEBRA_ASSERT(m_pRowIndex==rhs.m_pRowIndex &&
			          m_pColIndex==rhs.m_pColIndex, "Not same indices.");
			for(size_t fct1=0; fct1 < m_CplMat.num_rows(); ++fct1)
				for(size_t fct2=0; fct2 < m_CplMat.num_cols(); ++fct2)
					m_CplMat(fct1,fct2) += rhs.m_CplMat(fct1,fct2);
			return *this;
		}

	/// subtract matrix
		this_type& operator-=(const this_type& rhs)
		{
			UG_LOCALALGEBRA_ASSERT(m_pRowIndex==rhs.m_pRowIndex &&
			          m_pColIndex==rhs.m_pColIndex, "Not same indices.");
			for(size_t fct1=0; fct1 < m_CplMat.num_rows(); ++fct1)
				for(size_t fct2=0; fct2 < m_CplMat.num_cols(); ++fct2)
					m_CplMat(fct1,fct2) -= rhs.m_CplMat(fct1,fct2);
			return *this;
		}

	/// add scaled matrix
		this_type& scale_append(number s, const this_type& rhs)
		{
			UG_LOCALALGEBRA_ASSERT(m_pRowIndex==rhs.m_pRowIndex &&
					  m_pColIndex==rhs.m_pColIndex, "Not same indices.");
			for(size_t fct1=0; fct1 < m_CplMat.num_rows(); ++fct1)
				for(size_t fct2=0; fct2 < m_CplMat.num_cols(); ++fct2)
					MatScaleAppend(m_CplMat(fct1,fct2), s, rhs.m_CplMat(fct1,fct2));
			return *this;
		}

		///////////////////////////
		// restricted DoF access
		///////////////////////////

	/// access only part of the functions using mapping (restrict functions)
		void access_by_map(const FunctionIndexMapping& funcMap)
		{
			access_by_map(funcMap, funcMap);
		}

	/// access only part of the functions using mapping (restrict functions)
		void access_by_map(const FunctionIndexMapping& rowFuncMap,
		                   const FunctionIndexMapping& colFuncMap)
		{
			m_pRowFuncMap = &rowFuncMap;
			m_pColFuncMap = &colFuncMap;

			for(size_t i = 0; i < m_pRowFuncMap->num_fct(); ++i)
				for(size_t j = 0; j < m_pColFuncMap->num_fct(); ++j)
				{
					const size_t rowMapFct = rowFuncMap[i];
					const size_t colMapFct = colFuncMap[j];

					m_CplMatAcc(i,j) = &(m_CplMat(rowMapFct, colMapFct));
				}
		}

	///	access all functions
		void access_all()
		{
			m_pRowFuncMap = nullptr;
			m_pColFuncMap = nullptr;

			if(m_pRowIndex==nullptr) {m_CplMatAcc.resize(0,0); return;}

			for(size_t i = 0; i < m_pRowIndex->num_fct(); ++i)
				for(size_t j = 0; j < m_pColIndex->num_fct(); ++j)
					m_CplMatAcc(i,j) = &(m_CplMat(i,j));
		}

	///	returns the number of currently accessible (restricted) functions
		size_t num_row_fct() const
		{
			if(m_pRowFuncMap != nullptr) return m_pRowFuncMap->num_fct();
			return m_CplMatAcc.num_rows();
		}

	///	returns the number of currently accessible (restricted) functions
		size_t num_col_fct() const
		{
			if(m_pColFuncMap != nullptr) return m_pColFuncMap->num_fct();
			return m_CplMatAcc.num_cols();
		}

	///	returns the number of dofs for the currently accessible (restricted) function
		size_t num_row_dof(size_t fct) const
		{
			if(m_CplMat.num_rows()==0) return 0;
			if(m_pRowFuncMap == nullptr) return m_CplMat(fct, 0).num_rows();
			else return m_CplMat( (*m_pRowFuncMap)[fct], 0).num_rows();
		}

	///	returns the number of dofs for the currently accessible (restricted) function
		size_t num_col_dof(size_t fct) const
		{
			if(m_CplMat.num_cols()==0) return 0;
			if(m_pRowFuncMap == nullptr) return m_CplMat(0, fct).num_cols();
			else return m_CplMat( 0, (*m_pColFuncMap)[fct]).num_cols();
		}

	/// access to (restricted) coupling (rowFct, rowDoF) x (colFct, colDoF)
		number& operator()(size_t rowFct, size_t rowDoF,
		                   size_t colFct, size_t colDoF)
		{
			check_dof(rowFct, rowDoF, colFct, colDoF);
			return (*m_CplMatAcc(rowFct,colFct))(rowDoF, colDoF);
		}

	/// const access to (restricted) coupling (rowFct, rowDoF) x (colFct, colDoF)
		number operator()(size_t rowFct, size_t rowDoF,
		                        size_t colFct, size_t colDoF) const
		{
			check_dof(rowFct, rowDoF, colFct, colDoF);
			return (*m_CplMatAcc(rowFct,colFct))(rowDoF, colDoF);
		}

		///////////////////////////
		// all DoF access
		///////////////////////////

	///	returns the number of all functions
		size_t num_all_row_fct() const{	return m_CplMat.num_rows();}

	///	returns the number of all functions
		size_t num_all_col_fct() const{return m_CplMat.num_cols();}

	///	returns the number of dofs for a function
		size_t num_all_row_dof(size_t fct) const {return m_CplMat(fct, 0).num_rows();}

	///	returns the number of dofs for a function
		size_t num_all_col_dof(size_t fct) const {return m_CplMat(0, fct).num_cols();}

	/// access to coupling (rowFct, rowDoF) x (colFct, colDoF)
		number& value(size_t rowFct, size_t rowDoF,
		              size_t colFct, size_t colDoF)
		{
			check_all_dof(rowFct, rowDoF, colFct, colDoF);
			return (m_CplMat(rowFct,colFct))(rowDoF, colDoF);
		}

	/// const access to coupling (rowFct, rowDoF) x (colFct, colDoF)
		number value(size_t rowFct, size_t rowDoF,
		                   size_t colFct, size_t colDoF) const
		{
			check_all_dof(rowFct, rowDoF, colFct, colDoF);
			return (m_CplMat(rowFct,colFct))(rowDoF, colDoF);
		}

	protected:
	///	checks correct (fct1,fct2) index in debug mode
		inline void check_fct(size_t rowFct, size_t colFct) const
		{
			UG_LOCALALGEBRA_ASSERT(rowFct < num_row_fct(), "Wrong index.");
			UG_LOCALALGEBRA_ASSERT(colFct < num_col_fct(), "Wrong index.");
		}
	///	checks correct (dof1,dof2) index in debug mode
		inline void check_dof(size_t rowFct, size_t rowDoF,
		                      size_t colFct, size_t colDoF) const
		{
			check_fct(rowFct, colFct);
			UG_LOCALALGEBRA_ASSERT(rowDoF < num_row_dof(rowFct), "Wrong index.");
			UG_LOCALALGEBRA_ASSERT(colDoF < num_col_dof(colFct), "Wrong index.");
		}
	///	checks correct (fct1,fct2) index in debug mode
		inline void check_all_fct(size_t rowFct, size_t colFct) const
		{
			UG_LOCALALGEBRA_ASSERT(rowFct < num_all_row_fct(), "Wrong index.");
			UG_LOCALALGEBRA_ASSERT(colFct < num_all_col_fct(), "Wrong index.");
		}
	///	checks correct (dof1,dof2) index in debug mode
		inline void check_all_dof(size_t rowFct, size_t rowDoF,
							  size_t colFct, size_t colDoF) const
		{
			check_all_fct(rowFct, colFct);
			UG_LOCALALGEBRA_ASSERT(rowDoF < num_all_row_dof(rowFct), "Wrong index.");
			UG_LOCALALGEBRA_ASSERT(colDoF < num_all_col_dof(colFct), "Wrong index.");
		}

	protected:
	// 	Row indices
		const LocalIndices* m_pRowIndex;

	//	Column indices
		const LocalIndices* m_pColIndex;

	/// Row Access Mapping
		const FunctionIndexMapping* m_pRowFuncMap;

	/// Column Access Mapping
		const FunctionIndexMapping* m_pColFuncMap;

	//	\todo: Think of a better (faster) storage (one plain array?)

	//	type of cpl matrices between two functions
		using LocalCplMatrix = DenseMatrix<VariableArray2<value_type> >;

	//	type of Func-Coupling matrices
		using FctCplMatrix = DenseMatrix<VariableArray2<LocalCplMatrix> >;

	//	type of Func-Coupling pointer matrices
		using FctCplAccMatrix = DenseMatrix<VariableArray2<LocalCplMatrix*> >;

	// 	Entries (fct1, fct2, dof1, dof2)
		FctCplMatrix m_CplMat;

	// 	Entries (fct1, fct2, dof1, dof2)
		FctCplAccMatrix m_CplMatAcc;
};

inline
std::ostream& operator<< (std::ostream& outStream, const ug::LocalMatrix& mat)
{
	for(size_t fct1 = 0; fct1 < mat.num_row_fct(); ++fct1)
		for(size_t fct2 = 0; fct2 < mat.num_col_fct(); ++fct2)
			for(size_t dof1 = 0; dof1 < mat.num_row_dof(fct1); ++dof1)
				for(size_t dof2 = 0; dof2 < mat.num_col_dof(fct2); ++dof2)
				{
					outStream << "[("<< fct1 << "," << dof1 << ") x ("
							  	     << fct2 << "," << dof2 << ")]: "
							  	     << mat(fct1, dof1, fct2, dof2) << "\n";
				}
	return outStream;
}

inline
std::ostream& operator<< (std::ostream& outStream, const ug::LocalVector& vec)
{
	for(size_t fct = 0; fct < vec.num_fct(); ++fct)
		for(size_t dof = 0; dof < vec.num_dof(fct); ++dof)
		{
			outStream << "["<<fct<<","<<dof<<"]:" << vec(fct, dof) << "\n";
		}
	return outStream;
}

template <typename TVector>
void GetLocalVector(LocalVector& lvec, const TVector& vec)
{
	const LocalIndices& ind = lvec.get_indices();

	for(size_t fct=0; fct < lvec.num_all_fct(); ++fct)
		for(size_t dof=0; dof < lvec.num_all_dof(fct); ++dof)
		{
			const size_t index = ind.index(fct,dof);
			const size_t comp = ind.comp(fct,dof);
			lvec.value(fct,dof) = BlockRef(vec[index], comp);
		}
}

template <typename TVector>
void AddLocalVector(TVector& vec, const LocalVector& lvec)
{
	const LocalIndices& ind = lvec.get_indices();

	for(size_t fct=0; fct < lvec.num_all_fct(); ++fct)
		for(size_t dof=0; dof < lvec.num_all_dof(fct); ++dof)
		{
			const size_t index = ind.index(fct,dof);
			const size_t comp = ind.comp(fct,dof);
			BlockRef(vec[index], comp) += lvec.value(fct,dof);
		}
}

template <typename TMatrix>
void AddLocalMatrixToGlobal(TMatrix& mat, const LocalMatrix& lmat)
{
	//PROFILE_FUNC_GROUP("algebra")
	const LocalIndices& rowInd = lmat.get_row_indices();
	const LocalIndices& colInd = lmat.get_col_indices();

	for(size_t fct1=0; fct1 < lmat.num_all_row_fct(); ++fct1)
		for(size_t dof1=0; dof1 < lmat.num_all_row_dof(fct1); ++dof1)
		{
			const size_t rowIndex = rowInd.index(fct1,dof1);
			const size_t rowComp = rowInd.comp(fct1,dof1);

			for(size_t fct2=0; fct2 < lmat.num_all_col_fct(); ++fct2)
				for(size_t dof2=0; dof2 < lmat.num_all_col_dof(fct2); ++dof2)
				{
					const size_t colIndex = colInd.index(fct2,dof2);
					const size_t colComp = colInd.comp(fct2,dof2);

					BlockRef(mat(rowIndex, colIndex), rowComp, colComp)
								+= lmat.value(fct1,dof1,fct2,dof2);
				}
		}
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__COMMON__LOCAL_ALGEBRA__ */
