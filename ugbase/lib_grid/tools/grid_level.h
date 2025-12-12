/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_GRID__TOOLS__GRID_LEVEL__
#define __H__UG__LIB_GRID__TOOLS__GRID_LEVEL__

#include <string>

#include "common/error.h"

namespace ug {


class GridLevel
{
	public:
	///	indicates that always top level should be used
		enum{TOP = -1};

	///	type of view
		enum class ViewType : bool {LEVEL = false, SURFACE = true};


	public:
	///	constructor creation surface grid
		GridLevel() : m_level(TOP), m_type(ViewType::SURFACE), m_bWithGhosts(false) {}

	///	constructor
		GridLevel(int level, ViewType type, bool bWithGhosts = false)
			: m_level(level), m_type(type), m_bWithGhosts(bWithGhosts)
		{}

	///	constructor
	explicit GridLevel(int level)
			: m_level(level), m_type(ViewType::SURFACE), m_bWithGhosts(false)
		{}

	///	constructor
		GridLevel(int level, const std::string& type)
			: m_level(level), m_bWithGhosts(false)
		{
			if(type == "top") {m_type = ViewType::LEVEL;}
			else if(type == "surf") {m_type = ViewType::SURFACE;}
			else UG_THROW("Grid Level Type not in ['top' | 'surf'].");
		}

	///	returns the level
		[[nodiscard]] int level() const {return m_level;}

	///	returns the type
		[[nodiscard]] ViewType type() const {return m_type;}

	///	returns if ghosts are considered as part of the level
		[[nodiscard]] bool ghosts() const {return m_bWithGhosts;}

	///	returns if top level
		[[nodiscard]] bool top() const {return level() == TOP;}

	///	returns if type is level
		[[nodiscard]] bool is_level() const {return type() == ViewType::LEVEL;}

	///	returns if type is surface
		[[nodiscard]] bool is_surface() const {return type() == ViewType::SURFACE;}

	///	operator ==
		bool operator == (const GridLevel& rhs) const {
			return (this->level() == rhs.level() && this->type() == rhs.type()
					&& this->ghosts() == rhs.ghosts());
		}

	///	operator !=
		bool operator != (const GridLevel& rhs) const {
			return !(this->operator == (rhs));
		}

	///	operator <
		bool operator < (const GridLevel& rhs) const
		{
			if(this->type() != rhs.type()) return this->type() < rhs.type();
			if(this->ghosts() != rhs.ghosts()) return !this->ghosts();
			if(this->level() == rhs.level()) return false;
			if(this->level() == TOP) return false;
			if(rhs.level() == TOP) return true;
			return this->level() < rhs.level();
		}

	///	operator >
		bool operator > (const GridLevel& rhs) const
		{
			if(this->type() != rhs.type()) return this->type() > rhs.type();
			if(this->ghosts() != rhs.ghosts()) return this->ghosts();
			if(this->level() == rhs.level()) return false;
			if(this->level() == TOP) return true;
			if(rhs.level() == TOP) return false;
			return this->level() > rhs.level();
		}

	///	operator <=
		bool operator <= (const GridLevel& rhs) const
		{
			return (*this < rhs || *this == rhs);
		}

	///	operator >=
		bool operator >= (const GridLevel& rhs) const
		{
			return (*this > rhs || *this == rhs);
		}

	protected:
		int m_level; ///< the grid level
		ViewType m_type;	 ///< type (i.e. surface or level view)
		bool m_bWithGhosts;	///< with ghosts (only senseful in parallel)
};

/// writes to the output stream
std::ostream& operator << (std::ostream& out, const GridLevel& v);

std::ostream& operator << (std::ostream& stream, const GridLevel::ViewType & type);

/// returns appendix for a grid level
std::string GridLevelAppendix(const GridLevel& gl, int minfill = 2);

} // end namespace ug

#endif