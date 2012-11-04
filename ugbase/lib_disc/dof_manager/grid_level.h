/*
 * grid_level.h
 *
 *  Created on: 15.12.2011
 *      Author: andreasvogel
 */

#ifndef GRID_LEVEL_H_
#define GRID_LEVEL_H_


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
		GridLevel() : m_level(TOPLEVEL), m_type(SURFACE) {}

	///	constructor
		GridLevel(int level, ViewType type) : m_level(level), m_type(type) {}

	///	constructor
		GridLevel(int level) : m_level(level), m_type(SURFACE) {}

	///	constructor
		GridLevel(int level, std::string type) : m_level(level)
		{
			if(type == "top") {m_type = LEVEL;}
			else if(type == "surf") {m_type = SURFACE;}
			else UG_THROW("Grid Level Type not in ['top' | 'surf'].");
		}

	///	returns the level
		int level() const {return m_level;}

	///	returns the type
		ViewType type() const {return m_type;}

	///	operator ==
		bool operator==(const GridLevel& rhs) const {
			return (this->level() == rhs.level() && this->type() == rhs.type());
		}

	///	operator !=
		bool operator!=(const GridLevel& rhs) const {
			return !(this->operator==(rhs));
		}

	protected:
		int m_level; ///< the grid level
		ViewType m_type;	 ///< type (i.e. surface or level view)
};

/// writes to the output stream
inline std::ostream& operator<<(std::ostream& out,	const GridLevel& v)
{
	if(v.type() == GridLevel::SURFACE) out << "(surface, ";
	else if(v.type() == GridLevel::LEVEL) out << "(level, ";
	else UG_THROW("type of GridLevel not found.");

	if(v.level() == GridLevel::TOPLEVEL) out << "toplevel)";
	else out << v.level() << ")";
	return out;
}


} // end namespace ug

#endif /* GRID_LEVEL_H_ */
