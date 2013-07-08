/*
 * profilenode_management.h
 *
 */

#ifndef PROFILENODE_MANAGEMENT_H_
#define PROFILENODE_MANAGEMENT_H_

#include <stack>

class AutoProfileNode;

class ProfileNodeManager
{
	public:
		static void add(AutoProfileNode* node);
		static void release_latest();

	private:
		ProfileNodeManager();
		~ProfileNodeManager();
		static ProfileNodeManager& inst();

	private:
		std::stack<AutoProfileNode*>	m_nodes;
};

class AutoProfileNode
{
	friend class ProfileNodeManager;

	public:
#ifdef UG_PROFILER_SHINY
		AutoProfileNode();
#endif
#if defined(UG_PROFILER_SCALASCA) || defined(UG_PROFILER_VAMPIR)
		AutoProfileNode(const char* name);
#endif
#ifdef UG_PROFILER_SCOREP
		AutoProfileNode(SCOREP_User_RegionHandle name);
#endif
		~AutoProfileNode();

	private:
		void release();
		inline bool is_active()		{return m_bActive;}

	private:
		bool m_bActive;
#if defined(UG_PROFILER_SCALASCA) || defined(UG_PROFILER_VAMPIR)
		const char* m_pName;
#endif
#ifdef UG_PROFILER_SCOREP
		SCOREP_User_RegionHandle m_pHandle;
#endif
};


#endif /* PROFILENODE_MANAGEMENT_H_ */
