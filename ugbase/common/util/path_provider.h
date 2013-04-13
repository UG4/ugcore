// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 22.03.2011 (m,d,y)

#ifndef __H__UG__path_provider__
#define __H__UG__path_provider__

#include <string>
#include <map>
#include <stack>

#include "common/util/file_util.h"

namespace ug
{

/// \addtogroup ugbase_common_util
/// \{

////////////////////////////////////////////////////////////////////////
///	Constants used by PathProvider
enum PathTypes
{
	APP_PATH = 0,	///< path in which the binary lies
	DATA_PATH,
	SCRIPT_PATH,
	ROOT_PATH,
	PLUGIN_PATH,
	GRID_PATH,
	APPS_PATH,		///< path in which the application-scripts lie

//	always last
	MAX_PATH_CONSTANT
};

////////////////////////////////////////////////////////////////////////
///	Singleton which stores common paths and a stack of current paths.
/**	All paths are initially set to "".
 *
 * Note that all public methods of PathProvider are static. That means
 * you have to call them through the :: operator. E.g.
 *
 * \code
 * std::string appPath = PathProvider::get_path(APP_PATH);
 * \endcode
 */
class PathProvider
{
	public:
	///	sets the path for the given constant.
	/**
	 * \param pathType	should be one of the constants enumerated in PathTypes
	 * 					or a used defined constant starting from
	 * 					MAX_PATH_CONSTANT + 1.
	 */
		static inline void set_path(PathTypes pathType, const std::string& path)
		{inst().m_map[pathType] = path;}

	///	returns the path associated with the given constant.
	/**
	 * \param pathType	should be one of the constants enumerated in PathTypes
	 * 					or a used defined constant starting from
	 * 					MAX_PATH_CONSTANT + 1.
	 */
		static inline const std::string& get_path(PathTypes pathType)
		{return inst().m_map[pathType];}

	///	returns true, if the path associated with the given constant exists.
	/**
	 * \param pathType	should be one of the constants enumerated in PathTypes
	 * 					or a used defined constant starting from
	 * 					MAX_PATH_CONSTANT + 1.
	 */
		static inline bool has_path(PathTypes pathType)
		{return inst().m_map.find(pathType) != inst().m_map.end();}

	///	returns the current path
	/**	current paths are stored in a stack. The top of the stack is considered
	 * to be the most current path and is returned by this method.
	 *
	 * \param defPath	(optional) If the stack is empty, the path associated with
	 * 					defPath is returned. By default defPath is set to APP_PATH.
	 */
		static inline const std::string& get_current_path(PathTypes defPath = APP_PATH)
		{
			if(inst().m_curPaths.empty())
				return get_path(defPath);
			return inst().m_curPaths.top();
		}

	///	returns true if a current path exists, false if not.
		static inline bool has_current_path()
		{return !inst().m_curPaths.empty();}

	///	pushes a path to the stack of current paths
		static inline void push_current_path(const std::string& path)
		{inst().m_curPaths.push(path);}

	///	pops a path from the stack of current paths
		static inline void pop_current_path()
		{inst().m_curPaths.pop();}

	///	clears the stack of current paths. This makes sense if an error was catched.
		static inline void clear_current_path_stack()
		{while(has_current_path()) {pop_current_path();}}

	/**
	 * @param relativeFilename (in) relative filename
	 * @param absoluteFilename (out) absolute filename
	 * @return true if file exists relative to current path
	 */
		static inline bool get_filename_relative_to_current_path(const std::string &relativeFilename, std::string &absoluteFilename)
		{
			if(inst().has_current_path() == false) return false;
			absoluteFilename = inst().get_current_path() + "/" + relativeFilename;
			return FileExists(absoluteFilename.c_str());
		}

	/**
	 * @param relativeFilename (in) relative filename
	 * @param absoluteFilename (out) absolute filename
	 * @return true if file exists relative to current path
	 */
		static inline bool get_filename_relative_to_path(PathTypes pathType, const std::string &relativeFilename, std::string &absoluteFilename)
		{
			if(inst().has_path(pathType) == false) return false;
			absoluteFilename = inst().get_path(pathType) + "/" + relativeFilename;
			return FileExists(absoluteFilename.c_str());
		}

	private:
		PathProvider()	{}
		PathProvider(const PathProvider&)	{}

		static PathProvider& inst()
		{
			static PathProvider pp;
			return pp;
		}

	private:
		std::map<PathTypes, std::string>	m_map;
		std::stack<std::string>		m_curPaths;
};

// end group ugbase_common_util
/// \}

}//	end of namespace

#endif
