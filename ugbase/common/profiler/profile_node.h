/*
 * profile_node.h
 *
 *  Created on: May 22, 2012
 *      Author: Martin Rupp
 */

#ifndef __H_UG__PROFILE_NODE__
#define __H_UG__PROFILE_NODE__

#include <string>
#include <vector>
#include <sstream>

namespace ug
{

/**
 * UGProfilerNode class for more information about Shiny's ProfileNode.
 *
 * \note do NOT introduce variables or virtual functions to this class.
 * Shiny::ProfileNode are cast directly to UGProfilerNode and therefore are
 * assumed to have exactly the same size and format.
 * If you really need to change that you'd have to change the whole GetProfileNode-process.
 */
class UGProfilerNode
#if SHINY_PROFILER
: public Shiny::ProfileNode
#endif
{
public:
	/// \return number of entries in this profiler node
	double get_avg_entry_count() const;

	/// \return time in milliseconds spend in this node excluding subnodes
	double get_avg_self_time_ms() const;

	/// \return time in milliseconds spend in this node including subnodes
	double get_avg_total_time_ms() const;

	/// \return time in seconds spend in this node excluding subnodes
	double get_avg_self_time() const;

	/// \return time in seconds spend in this node including subnodes
	double get_avg_total_time() const;

	std::string call_tree(double dSkipMarginal) const;

	std::string call_tree() const;

	std::string child_self_time_sorted(double dSkipMarginal) const;

	std::string child_self_time_sorted() const;

	std::string total_time_sorted(double dSkipMarginal) const;

	std::string total_time_sorted() const;

	std::string entry_count_sorted(double dSkipMarginal) const;

	std::string entry_count_sorted() const;


	/// \return true if node has been found
	bool valid() const;

#if SHINY_PROFILER
private:
	std::string print_node(double full, size_t offset=0) const;
	const UGProfilerNode *get_first_child() const;
	const UGProfilerNode *get_last_child() const;
	const UGProfilerNode *get_next_sibling() const;

	void rec_print(double full, std::stringstream &s, size_t offset, double dSkipMarginal) const;
	void add_nodes(std::vector<const UGProfilerNode*> &nodes) const;
	std::string child_sorted(const char *name, bool sortFunction(const UGProfilerNode *a, const UGProfilerNode *b),
			double dSkipMarginal) const;
	static void log_header(std::stringstream &s, const char *name);
	static bool self_time_sort(const UGProfilerNode *a, const UGProfilerNode *b);
	static bool total_time_sort(const UGProfilerNode *a, const UGProfilerNode *b);
	static bool entry_count_sort(const UGProfilerNode *a, const UGProfilerNode *b);
#endif

	// do NOT add variables or virtual functions here (see above).
};


const UGProfilerNode *GetProfileNode(const char *name);
bool GetProfilerAvailable();

}


#endif /* __H_UG__PROFILE_NODE__ */
