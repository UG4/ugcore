/*
 * pcl_distributed_domain_info.h
 *
 *  Created on: 27.01.2011
 *      Author: sreiter, avogel, iheppner
 */

#ifndef __H__PCL__DISTRIBUTED_DOMAIN_DECOMPOSITION__
#define __H__PCL__DISTRIBUTED_DOMAIN_DECOMPOSITION__

namespace pcl
{

class IDomainDecompositionInfo
{

    public:
	/// 	mapping method "proc-id" ==> "subdomain-id"
	/**
	 * This functions determines the subdomain a processor lives
	 * \param[out] 	procID		id of processor
	 * \return 		int			id of subdomain the processor operates on
	 */
		virtual int map_proc_id_to_subdomain_id(int procID) const = 0;

		virtual int get_num_subdomains() const = 0;

		virtual int get_num_procs_per_subdomain() const = 0;

	public:
		// destructor
		virtual ~IDomainDecompositionInfo() {};

}; /* end class 'IDomainDecompositionInfo' */

class StandardDomainDecompositionInfo : public IDomainDecompositionInfo
{
	///	constructors
    public:
		StandardDomainDecompositionInfo() :
			m_num_subdomains(1),
			m_num_procs_per_subdomain(1)
			{}

		StandardDomainDecompositionInfo(int numSubdomains) :
			m_num_subdomains(numSubdomains),
			m_num_procs_per_subdomain(1)
			/*m_pDebugWriter(NULL)*/
		{
			if(numSubdomains > 0 && pcl::GetNumProcesses() > 1)
				m_num_procs_per_subdomain = pcl::GetNumProcesses() / numSubdomains;
		}

    public:
	/// 	mapping method "proc-id" ==> "subdomain-id"
	/**
	 * This functions determines the subdomain a processor lives
	 * \param[out] 	procID		id of processor
	 * \return 		int			id of subdomain the processor operates on
	 */
		int map_proc_id_to_subdomain_id(int procID) const
        {
			return procID / m_num_procs_per_subdomain; 
        }

		void set_num_subdomains(int numSubdomains) {m_num_subdomains = numSubdomains;}

		int get_num_subdomains() const {return m_num_subdomains;}

		int get_num_procs_per_subdomain() const {return m_num_procs_per_subdomain;}

	///	set debug output - needed?
		/* 
		void set_debug(IDebugWriter<algebra_type>* debugWriter)
		{
			m_pDebugWriter = debugWriter;
		}
	protected:
		bool write_debug(const vector_type& vec, const char* filename)
		{
		//	if no debug writer set, we're done
			if(m_pDebugWriter == NULL) return true;

		//	write
			return m_pDebugWriter->write_vector(vec, filename);
		}
		*/

	public:
		// destructor
		~StandardDomainDecompositionInfo() {};

	protected:
	// 	number of subdomains
		int m_num_subdomains;

	// 	number of procs per subdomains
		//std::vector<int> m_vnum_procs_per_subdomain; // as vector, for variable distribution of processors over subdomains
		int m_num_procs_per_subdomain;

	//	Debug Writer - needed?
			/* 
		IDebugWriter<algebra_type>* m_pDebugWriter;
			*/
}; /* end class 'StandardDomainDecompositionInfo' */


}//	end of namespace
#endif
