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
 * UGProfileNode class for more information about Shiny's ProfileNode.
 *
 * \note do NOT introduce variables or virtual functions to this class.
 * Shiny::ProfileNode are cast directly to UGProfileNode and therefore are
 * assumed to have exactly the same size and format.
 * If you really need to change that you'd have to change the whole GetProfileNode-process.
 */
class UGProfileNode
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

	double get_self_mem() const;
	double get_total_mem() const;

	std::string call_tree(double dSkipMarginal) const;

	std::string call_tree() const;

	std::string child_self_time_sorted(double dSkipMarginal) const;

	std::string child_self_time_sorted() const;

	std::string total_time_sorted(double dSkipMarginal) const;

	std::string total_time_sorted() const;

	std::string entry_count_sorted(double dSkipMarginal) const;

	std::string entry_count_sorted() const;

	std::string groups() const;

	/// \return true if node has been found
	bool valid() const;

	static const UGProfileNode *get_root();

#if SHINY_PROFILER
public:
	void add_nodes(std::vector<const UGProfileNode*> &nodes) const;
	void write_node(std::ostream &s) const;
	const UGProfileNode *get_first_child() const;
	const UGProfileNode *get_last_child() const;
	const UGProfileNode *get_next_sibling() const;
	const UGProfileNode *find_next_in_tree() const;
	std::string print_node(double full, double fullMem, size_t offset=0) const;
	static void log_header(std::stringstream &s, const char *name);
private:
	std::string get_mem_info(double fullMem) const;


	void rec_print(double full, double fullMem, std::stringstream &s, size_t offset, double dSkipMarginal) const;
	std::string child_sorted(const char *name, bool sortFunction(const UGProfileNode *a, const UGProfileNode *b),
			double dSkipMarginal) const;

	static bool self_time_sort(const UGProfileNode *a, const UGProfileNode *b);
	static bool total_time_sort(const UGProfileNode *a, const UGProfileNode *b);
	static bool entry_count_sort(const UGProfileNode *a, const UGProfileNode *b);	
#endif

	// do NOT add variables or virtual functions here (see above).
};


const UGProfileNode *GetProfileNode(const char *name);
bool GetProfilerAvailable();
void WriteProfileData(const char *filename);
}


#endif /* __H_UG__PROFILE_NODE__ */
