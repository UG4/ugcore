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

	/**
	 * @param dSkipMarginal 	nodes with full*dSkipMarginal > node->full[ms or mem] are skipped
	 * @return call tree profile information
	 */
	std::string call_tree(double dSkipMarginal) const;
	std::string call_tree() const { return call_tree(0.0);	}

	/**
	 * @param dSkipMarginal 	nodes with full*dSkipMarginal > node->full[ms or mem] are skipped
	 * @return a table with all profile nodes information, sorted by self time
	 */
	std::string child_self_time_sorted(double dSkipMarginal) const;
	std::string child_self_time_sorted() const { return child_self_time_sorted(0.0); }

	/**
	 * @param dSkipMarginal 	nodes with full*dSkipMarginal > node->full[ms or mem] are skipped
	 * @return a table with all profile nodes information, sorted by total time
	 */
	std::string total_time_sorted(double dSkipMarginal) const;
	std::string total_time_sorted() const { return total_time_sorted(0.0); }


	/**
	 * @param dSkipMarginal 	nodes with full*dSkipMarginal > node->full[ms or mem] are skipped
	 * @return a table with all profile nodes information, sorted by self memory
	 */
	std::string child_self_memory_sorted(double dSkipMarginal) const;
	std::string child_self_memory_sorted() const { return child_self_memory_sorted(0.0); }

	/**
	 * @param dSkipMarginal 	nodes with full*dSkipMarginal > node->full[ms or mem] are skipped
	 * @return a table with all profile nodes information, sorted by total memory
	 */
	std::string total_memory_sorted(double dSkipMarginal) const;
	std::string total_memory_sorted() const { return total_memory_sorted(0.0); }


	/**
	 * @param dSkipMarginal 	nodes with full*dSkipMarginal > node->full[ms or mem] are skipped
	 * @return a table with all profile nodes information, sorted by entry count
	 */
	std::string entry_count_sorted(double dSkipMarginal) const;
	std::string entry_count_sorted() const { return entry_count_sorted(0.0); }

	/**
	 * @return Profiling group information
	 */
	std::string groups() const;

	/// \return true if node has been found
	bool valid() const;

	static const UGProfileNode *get_root();

	static void CheckForTooSmallNodes();
#if SHINY_PROFILER

public:
	const UGProfileNode *get_first_child() const;
	const UGProfileNode *get_last_child() const;
	const UGProfileNode *get_next_sibling() const;
	const UGProfileNode *find_next_in_tree() const;


	/**
	 * @brief writes this node and its subnodes in PDXML format to an ostream buffer
	 * @param s ostream buffer
	 */
	void PDXML_rec_write(std::ostream &s) const;

	/**
	 * @brief prints the node information into a string
	 * @param fullMs	full time to calculate percentage of
	 * @param fullMem	full allocated mem to calculate percentage of
	 * @param offset	number of spaces to add from the left to get a tree-like structure (only applies to the name)
	 * @return tabular node information
	 */
	std::string print_node(double fullMs, double fullMem, size_t offset=0) const;

	/**
	 * recursively adds this node and its subnodes to the array 'nodes'
	 * @param nodes	vector to push nodes into
	 */
	void rec_add_nodes(std::vector<const UGProfileNode*> &nodes) const;

	static void log_header(std::stringstream &s, const char *name);
private:
	void check_for_too_small_nodes(double fullMs, std::map<std::string, const UGProfileNode *> &list) const;

	/**
	 * @brief prints the memory information to a node
	 * @param fullMem	full allocated mem to calculate percentage of
	 */
	std::string get_mem_info(double fullMem) const;


	/**
	 * @brief recursive print this node and its subnodes into stringstream s
	 * @param fullMs		full time to calculate percentage of
	 * @param fullMem		full allocated mem to calculate percentage of
	 * @param s				stringstream to print information into
	 * @param offset		number of spaces to add from the left to get a tree-like structure (only applies to the name)
	 * @param dSkipMarginal	if != 0.0, only nodes which have full*dSkipMarginal < this->total[ms or mem] are printed
	 * @return
	 */
	void rec_print(double fullMs, double fullMem, std::stringstream &s, size_t offset, double dSkipMarginal) const;


	/**
	 * @param name 			the name of the table to print
	 * @param sortFunction	how to sort the table
	 * @param dSkipMarginal	if != 0.0, only nodes which have full*dSkipMarginal < this->total[ms or mem] are printed
	 * @return table with sorted profiling information
	 * @sa self_time_sort, total_time_sort, entry_count_sort, self_memory_sort, total_memory_sort
	 */
	std::string print_child_sorted(const char *name, bool sortFunction(const UGProfileNode *a, const UGProfileNode *b),
			double dSkipMarginal) const;

	// sort functions
	static bool self_time_sort(const UGProfileNode *a, const UGProfileNode *b);
	static bool total_time_sort(const UGProfileNode *a, const UGProfileNode *b);
	static bool self_memory_sort(const UGProfileNode *a, const UGProfileNode *b);
	static bool total_memory_sort(const UGProfileNode *a, const UGProfileNode *b);
	static bool entry_count_sort(const UGProfileNode *a, const UGProfileNode *b);	
#endif

	// do NOT add variables or virtual functions here (see above).
};


/// This singleton represents a UGProfileNode that has not been found.
class UGProfileNodeNull
: public UGProfileNode
{
	public:
       static UGProfileNodeNull* getInstance()
       {
           static UGProfileNodeNull instance;
           return &instance;
       }

   private:
       UGProfileNodeNull() {};
       UGProfileNodeNull(UGProfileNodeNull const&); // do not implement
       void operator=(UGProfileNodeNull const&); // do not implement
};

#define PROFILER_NULL_NODE UGProfileNodeNull::getInstance()


const UGProfileNode *GetProfileNode(const char *name);
const UGProfileNode *GetProfileNode(const char *name, const UGProfileNode *node);
bool GetProfilerAvailable();

///	Writes profile data of process 0 to the specified file
void WriteProfileDataXML(const char *filename);

///	Writes profile data to the specified file
/**	Writes profile data of all procs (procId == -1) or of specified proc
 * (procId >= 0) to the specified file
 */
void WriteProfileDataXML(const char *filename, int procId);

void WriteCallLog(const char *filename);
void WriteCallLog(const char *filename, int procId);

}


#endif /* __H_UG__PROFILE_NODE__ */
