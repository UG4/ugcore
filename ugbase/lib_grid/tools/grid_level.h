/*
 * grid_level.h
 *
 *  Created on: 15.12.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_GRID__TOOLS__GRID_LEVEL__
#define __H__UG__LIB_GRID__TOOLS__GRID_LEVEL__


namespace ug{

class GridLevel
{
	public:
	///	indicates that always top level should be used
		enum{TOPLEVEL = -1};

	///	type of view
		enum ViewType{SURFACE = 0, LEVEL = 1};

	public:
	///	constructor creation surface grid
		GridLevel() : m_level(TOPLEVEL), m_type(SURFACE), m_bWithGhosts(false) {}

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
		bool with_ghosts() const {return m_bWithGhosts;}

	///	operator ==
		bool operator==(const GridLevel& rhs) const {
			return (this->level() == rhs.level() && this->type() == rhs.type()
					&& this->with_ghosts() == rhs.with_ghosts());
		}

	///	operator !=
		bool operator!=(const GridLevel& rhs) const {
			return !(this->operator==(rhs));
		}

	///	operator <
		bool operator<(const GridLevel& rhs) const
		{
			if(this->with_ghosts() != rhs.with_ghosts())
				return this->with_ghosts();
			else if(this->type() != rhs.type())
				return this->type() < rhs.type();
			else return this->level() < rhs.level();
		}

	///	operator >
		bool operator>(const GridLevel& rhs) const
		{
			if(this->with_ghosts() != rhs.with_ghosts())
				return !this->with_ghosts();
			else if(this->type() != rhs.type())
				return this->type() > rhs.type();
			else return this->level() > rhs.level();
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
inline std::ostream& operator<<(std::ostream& out,	const GridLevel& v)
{
	if(v.type() == GridLevel::SURFACE) out << "(surface, ";
	else if(v.type() == GridLevel::LEVEL) out << "(level, ";
	else UG_THROW("type of GridLevel not found.");

	if(v.level() == GridLevel::TOPLEVEL) out << "toplevel)";
	else out << v.level() << ')';
	return out;
}


} // end namespace ug

#endif /* __H__UG__LIB_GRID__TOOLS__GRID_LEVEL__ */
