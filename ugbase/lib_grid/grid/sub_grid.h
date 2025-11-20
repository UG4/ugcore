/*
 * Copyright (c) 2016:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG_sub_grid
#define __H__UG_sub_grid

namespace ug{

///	Instances represent a part of a grid.
/** One can iterate over the objects in a sub-grid and one may check whether
 * elements are contained in a sub-grid.
 *
 * \sa SubGrid
 */
class ISubGrid {
public:
	virtual ~ISubGrid() = default;

	virtual const GridObjectCollection& goc () const = 0;
	virtual bool is_contained(Vertex*) const = 0;
	virtual bool is_contained(Edge*) const = 0;
	virtual bool is_contained(Face*) const = 0;
	virtual bool is_contained(Volume*) const = 0;
};


///	specializes ISubGrid for general callback classes.
/**
 *\tparam TCallbackCls	A class that features the following methods
 *						- bool operator () (Vertex* v) const;
 *						- bool operator () (Edge* e) const;
 *						- bool operator () (Face* f) const;
 *						- bool operator () (Volume* v) const;
 */
template <typename TCallbackCls>
class SubGrid : public ISubGrid{
public:
	SubGrid(GridObjectCollection goc, const TCallbackCls& cb) :
		m_goc (goc),
		m_callbacks (cb)
	{
	}

	const GridObjectCollection& goc () const override {return m_goc;}
	bool is_contained (Vertex* e) const override {return m_callbacks(e);}
	bool is_contained (Edge* e) const override {return m_callbacks(e);}
	bool is_contained (Face* e) const override {return m_callbacks(e);}
	bool is_contained (Volume* e) const override {return m_callbacks(e);}

private:
	GridObjectCollection m_goc;
	TCallbackCls m_callbacks;
};


///	Callbacks that return true if an element is contained in a sub-grid.
/**	Please make sure that the associated sub-grid exists until the callback-class
 * was destroyed.*/
class IsInSubGrid
{
	public:
		IsInSubGrid(const ISubGrid& subGrid) :
			m_subGrid(subGrid)	{}

		bool operator () (Vertex* v)	{return callback(v);}
		bool operator () (Edge* e)	{return callback(e);}
		bool operator () (Face* f)	{return callback(f);}
		bool operator () (Volume* v)	{return callback(v);}

	private:
		template <typename TElem>
		bool callback(TElem* e)		{return m_subGrid.is_contained(e);}

	private:
		const ISubGrid&	m_subGrid;
};


///	Callbacks that return true if an element is not contained in a sub-grid.
/**	Please make sure that the associated sub-grid exists until the callback-class
 * was destroyed.*/
class IsNotInSubGrid
{
	public:
		IsNotInSubGrid(const ISubGrid& subGrid) :
			m_subGrid(subGrid)	{}

		bool operator () (Vertex* v)	{return callback(v);}
		bool operator () (Edge* e)	{return callback(e);}
		bool operator () (Face* f)	{return callback(f);}
		bool operator () (Volume* v)	{return callback(v);}

	private:
		template <typename TElem>
		bool callback(TElem* e)		{return !m_subGrid.is_contained(e);}

	private:
		const ISubGrid&	m_subGrid;
};

}//	end of namespace

#endif