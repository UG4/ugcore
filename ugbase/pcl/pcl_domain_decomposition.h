
#ifndef __H__PCL__DOMAIN_DECOMPOSITION__
#define __H__PCL__DOMAIN_DECOMPOSITION__

namespace pcl
{

/// \addtogroup pcl
/// \{

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

		//virtual int get_num_procs_per_subdomain() const = 0;
		virtual int get_num_spatial_dimensions() const = 0;

		virtual void get_subdomain_procs(std::vector<int>& procsOut,
										int subdomIndex) = 0;
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
			m_num_spatial_dimensions(2),
			m_num_procs_per_subdomain(1)
			{}

		StandardDomainDecompositionInfo(int numSubdomains) :
			m_num_subdomains(numSubdomains),
			m_num_procs_per_subdomain(-1)
			/*m_pDebugWriter(NULL)*/
		{
			if(numSubdomains > 0 && pcl::NumProcs() > 1)
				m_num_procs_per_subdomain = pcl::NumProcs() / numSubdomains;
			else
				m_num_procs_per_subdomain = -1;
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

		void set_num_subdomains(int numSubdomains) {
			m_num_subdomains = numSubdomains;

			if(numSubdomains > 0 && pcl::NumProcs() > 1)
				m_num_procs_per_subdomain = pcl::NumProcs() / numSubdomains;
			else
				m_num_procs_per_subdomain = -1;
		}

		int get_num_subdomains() const {return m_num_subdomains;}

		void set_num_spatial_dimensions(int dim) {
			m_num_spatial_dimensions = dim;
		}

		int get_num_spatial_dimensions() const {return m_num_spatial_dimensions;}

		//int get_num_procs_per_subdomain() const {return m_num_procs_per_subdomain;}

		virtual void get_subdomain_procs(std::vector<int>& procsOut,
										int subdomIndex)
		{
			procsOut.resize(m_num_procs_per_subdomain);
			int first = subdomIndex * m_num_procs_per_subdomain;
			for(int i = 0; i < m_num_procs_per_subdomain; ++i)
				procsOut[i] = first + i;
		}

	public:
		// destructor
		~StandardDomainDecompositionInfo() {};

	protected:
	// 	number of subdomains
		int m_num_subdomains;

	// 	number of spatial dimensions
		int m_num_spatial_dimensions;

	// 	number of procs per subdomains
		//std::vector<int> m_vnum_procs_per_subdomain; // as vector, for variable distribution of processors over subdomains
		int m_num_procs_per_subdomain;
}; /* end class 'StandardDomainDecompositionInfo' */

// end group pcl
/// \}

}//	end of namespace
#endif
