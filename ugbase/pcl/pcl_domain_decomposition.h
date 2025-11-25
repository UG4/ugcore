/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Sebastian Reiter, Andreas Vogel, Ingo Heppner
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

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
		[[nodiscard]] virtual int map_proc_id_to_subdomain_id(int procID) const = 0;

		[[nodiscard]] virtual int get_num_subdomains() const = 0;

		//virtual int get_num_procs_per_subdomain() const = 0;
		[[nodiscard]] virtual int get_num_spatial_dimensions() const = 0;

		virtual void get_subdomain_procs(std::vector<int>& procsOut,
										 int subdomIndex) = 0;
	public:
		// destructor
		virtual ~IDomainDecompositionInfo() = default;

};

class StandardDomainDecompositionInfo : public IDomainDecompositionInfo
{
	///	constructors
    public:
		StandardDomainDecompositionInfo() :
			m_num_subdomains(1),
			m_num_spatial_dimensions(2),
			m_num_procs_per_subdomain(1)
			{}

		explicit StandardDomainDecompositionInfo(int numSubdomains) :
			m_num_subdomains(numSubdomains),
			m_num_procs_per_subdomain(-1)
			/*m_pDebugWriter(nullptr)*/
		{
			if(numSubdomains > 0 && NumProcs() > 1)
				m_num_procs_per_subdomain = NumProcs() / numSubdomains;
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
		[[nodiscard]] int map_proc_id_to_subdomain_id(int procID) const override {
			return procID / m_num_procs_per_subdomain;
        }

		void set_num_subdomains(int numSubdomains) {
			m_num_subdomains = numSubdomains;

			if(numSubdomains > 0 && NumProcs() > 1)
				m_num_procs_per_subdomain = NumProcs() / numSubdomains;
			else
				m_num_procs_per_subdomain = -1;
		}

		[[nodiscard]] int get_num_subdomains() const override {
			return m_num_subdomains;
		}

		void set_num_spatial_dimensions(int dim) {
			m_num_spatial_dimensions = dim;
		}

		[[nodiscard]] int get_num_spatial_dimensions() const override {
			return m_num_spatial_dimensions;
		}

		//int get_num_procs_per_subdomain() const {return m_num_procs_per_subdomain;}

		void get_subdomain_procs(std::vector<int>& procsOut,
	                         int subdomIndex) override {
			procsOut.resize(m_num_procs_per_subdomain);
			int first = subdomIndex * m_num_procs_per_subdomain;
			for(int i = 0; i < m_num_procs_per_subdomain; ++i)
				procsOut[i] = first + i;
		}

	public:
		// destructor
		~StandardDomainDecompositionInfo() override = default;

	protected:
	// 	number of subdomains
		int m_num_subdomains;

	// 	number of spatial dimensions
		int m_num_spatial_dimensions;

	// 	number of procs per subdomains
		//std::vector<int> m_vnum_procs_per_subdomain; // as vector, for variable distribution of processors over subdomains
		int m_num_procs_per_subdomain;
};

// end group pcl
/// \}

}//	end of namespace
#endif
