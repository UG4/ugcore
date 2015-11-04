
#ifndef __H__UG__LIB_GRID__TOOLS__GRID_LEVEL__
#define __H__UG__LIB_GRID__TOOLS__GRID_LEVEL__

#include <string>
#include <iostream>

namespace ug{

class GridLevel
{
	public:
	///	indicates that always top level should be used
		enum{TOP = -1};

	///	type of view
		enum ViewType{LEVEL = 0, SURFACE = 1};

	public:
	///	constructor creation surface grid
		GridLevel() : m_level(TOP), m_type(SURFACE), m_bWithGhosts(false) {}

	///	constructor
		GridLevel(int level, ViewType type, bool bWithGhosts = false)
			: m_level(level), m_type(type), m_bWithGhosts(bWithGhosts)
		{}

	///	constructor
		GridLevel(int level)
			: m_level(level), m_type(SURFACE), m_bWithGhosts(false)
		{}

	///	constructor
		GridLevel(int level, const std::string& type)
			: m_level(level), m_bWithGhosts(false)
		{
			if(type == "top") {m_type = LEVEL;}
			else if(type == "surf") {m_type = SURFACE;}
			else UG_THROW("Grid Level Type not in ['top' | 'surf'].");
		}

	///	returns the level
		int level() const {return m_level;}

	///	returns the type
		ViewType type() const {return m_type;}

	///	returns if ghosts are considered as part of the level
		bool ghosts() const {return m_bWithGhosts;}

	///	returns if top level
		bool top() const {return level() == TOP;}

	///	returns if type is level
		bool is_level() const {return type() == LEVEL;}

	///	returns if type is surface
		bool is_surface() const {return type() == SURFACE;}

	///	operator ==
		bool operator==(const GridLevel& rhs) const {
			return (this->level() == rhs.level() && this->type() == rhs.type()
					&& this->ghosts() == rhs.ghosts());
		}

	///	operator !=
		bool operator!=(const GridLevel& rhs) const {
			return !(this->operator==(rhs));
		}

	///	operator <
		bool operator<(const GridLevel& rhs) const
		{
			if(this->type() != rhs.type()) return this->type() < rhs.type();
			if(this->ghosts() != rhs.ghosts()) return !this->ghosts();
			if(this->level() == rhs.level()) return false;
			if(this->level() == TOP) return false;
			if(rhs.level() == TOP) return true;
			return this->level() < rhs.level();
		}

	///	operator >
		bool operator>(const GridLevel& rhs) const
		{
			if(this->type() != rhs.type()) return this->type() > rhs.type();
			if(this->ghosts() != rhs.ghosts()) return this->ghosts();
			if(this->level() == rhs.level()) return false;
			if(this->level() == TOP) return true;
			if(rhs.level() == TOP) return false;
			return this->level() > rhs.level();
		}

	///	operator <=
		bool operator<=(const GridLevel& rhs) const
		{
			return (*this < rhs || *this == rhs);
		}

	///	operator >=
		bool operator>=(const GridLevel& rhs) const
		{
			return (*this > rhs || *this == rhs);
		}

	protected:
		int m_level; ///< the grid level
		ViewType m_type;	 ///< type (i.e. surface or level view)
		bool m_bWithGhosts;	///< with ghosts (only senseful in parallel)
};

/// writes to the output stream
std::ostream& operator<<(std::ostream& out,	const GridLevel& v);

/// returns appendix for a grid level
std::string GridLevelAppendix(const GridLevel& gl, int minfill = 2);

} // end namespace ug

#endif /* __H__UG__LIB_GRID__TOOLS__GRID_LEVEL__ */
