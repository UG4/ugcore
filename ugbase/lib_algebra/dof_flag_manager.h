/*
 * Copyright (c) 2026:  CEMSE, KAUST
 * Author: Dmitry Logashenko
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

#ifndef __H__UG__LIB_ALGEBRA__DOF_FLAG_MANAGER__
#define __H__UG__LIB_ALGEBRA__DOF_FLAG_MANAGER__

#include <vector>
#include <string>

#include "algebra_type.h"
#include "cpu_algebra/vector.h"

namespace ug{

/// DoF flag name and allocation manager
/**
 * This class assignes names to particular flag blocks in flag sets
 * and determines the number of the flag sets. It provides tools for
 * allocating the flags.
 *
 * Note that new flag blocks can be added but they cannot be removed
 * separately. (Otherwise, one had to restructure the flags in all the
 * vectors.)
 */
template <typename TAlgebra>
class DoFFlagManager
{
public:
	typedef TAlgebra algebra_type;
	typedef typename algebra_type::vector_type vector_type;
	typedef typename vector_type::flag_unit_type flag_unit_type;
	
	static const int blockSize = algebra_type::blockSize;
	static const int flagBlocksPerUnit = algebra_type::flagBlocksPerUnit;
	static const flag_unit_type flagBlockFrame
		= (blockSize <= 0)? ~0 : (1 << (blockSize + 1)) - 1; // 1's, blockSize times

/// Returns the static (single) flag manager for particular algebra
	static DoFFlagManager* get();

private:

	/// Private constructor: This class is used only for singletons
	DoFFlagManager () {};

	struct flag_block_alloc_type
	{
		std::string block_name; //< name of the flag block
		size_t set_idx; //< index of the set
		size_t block_idx; //< index of the block in the set (in block size units)
		
		flag_block_alloc_type(const char* name, size_t s, size_t i)
		:	block_name(name), set_idx(s), block_idx(i) {};
	};
	
///	finds the allocation structure by the name of the block
	flag_block_alloc_type * find_alloc_by_name (const char* name)
	{
		for(size_t i = 0; i < m_flag_blocks.size(); i++)
			if(m_flag_blocks[i].block_name == name)
				return &(m_flag_blocks[i]);
		return NULL;
	}
	
public:
	
///	allocates a new flag block (and possibly a new flag set)
	void add
	(
		const char* name ///< name for the new flag block
	)
	{
		if(m_flag_blocks.size() == 0) // no flags at all up to now
		{
			m_flag_blocks.push_back(flag_block_alloc_type(name, 0, 0));
			return;
		}
		
		if(find_alloc_by_name(name) != NULL)
		{
			UG_THROW("DoFFlagManager::add: Attempt to add a flag block '" << name
				<< "', but such name is already used.");
		}
		
		size_t last_alloc = m_flag_blocks.size() - 1;
		size_t set = m_flag_blocks[last_alloc].set_idx;
		size_t block = m_flag_blocks[last_alloc].block_idx;
		
		if (block + 1 >= flagBlocksPerUnit)
		{
		//	allocate a new set
			set++; block = 0;
		}
		else block++;
		
		m_flag_blocks.push_back(flag_block_alloc_type(name, set, block));
	}
	
///	gets the set number and the minor bit number of the flag block
	bool get_flag_offsets ///< the function returns false if the name is not found
	(
		const char* name, ///< name of the block
		size_t& set_idx, ///< index of the set
		size_t& bit_idx ///< index of the minor bit
	) const
	{
		flag_block_alloc_type * a = find_alloc_by_name(name);
		
		if(a == NULL) return false;
		
		set_idx = a->set_idx;
		bit_idx = a->block_idx * blockSize;
	}
	
	//! returns the flag block in the minor bits of the return value
	/**
	 * Remark: Note that the major bits may be arbitrary (contain the other flags).
	 * We do not care about them for efficiency reasons. Use flagBlockFrame
	 * to clean it.
	 */
	inline static flag_unit_type flag
	(
		const vector_type& vec, ///< the DoF vector to read the flag from
		size_t i, ///< algebra index of the dof block
		size_t s, ///< index of the flag set
		size_t f_offset ///< index of the minor bit
	)
	{
		return vec.flag(i,s) >> f_offset;
	}	

	//! sets a flag block at a given algebra index
	/**
	 * WARNING! Note that all the flags for a given block are set simultaneously.
	 *
	 * WARNING! We assume that only the last blockSize bits in 'flag_block'
	 * may be non-zero! Use flagBlockFrame to clean the other bits!
	 */
	inline static void set_flag
	(
		vector_type& vec, ///< the DoF vector to set the flag in
		size_t i, ///< algebra index of the DoF block
		size_t s, ///< index of the flag set
		size_t f_offset, ///< index of the minor bit of the flag block
		flag_unit_type flag_block ///< value of the block (in the minor bits) to set
	)
	{
		if(blockSize <= 0) // note that his is a constant expression
		{
		// Variable block size: We assume, there is a single flag block per flag unit
			vec.flag(i,s) = flag_block;
		}
		else
		{
			const flag_unit_type mask = flagBlockFrame << f_offset;
			
			vec.flag(i,s) = (vec.flag(i,s) & (~mask)) | (flag_block << f_offset);
		}
	}
///	creates the flags for a particular vector
	void create_flags
	(
		vector_type& v ///< the DoF vector to create the flags for
	) const
	{
		if(m_flag_blocks.size() == 0) // no flags at all up to now
			return;
		
		size_t num_sets = m_flag_blocks[m_flag_blocks.size()-1].set_idx + 1;
		
		v.create_flags(num_sets);
	}

///	prints the flag block allocations
	void print()
	{
		UG_LOG("Allocated DoF flags:\n");
		UG_LOG("--------------------\n");
		for(size_t i = 0; i < m_flag_blocks.size(); i++)
		{
			flag_block_alloc_type& a = m_flag_blocks[i];
			UG_LOG("Flag block '" << a.block_name << "': "
				"index " << a.block_idx << " in set " << a.set_idx << "\n");
		}
		UG_LOG("--------------------\n");
		UG_LOG("(DoF block size: " << blockSize << ")\n");
	}

private:

	std::vector<flag_block_alloc_type> m_flag_blocks;
};

/// Returns the static (single) flag manager for particular algebra
template <typename TAlgebra>
DoFFlagManager<TAlgebra>* DoFFlagManager<TAlgebra>::get()
{
	static DoFFlagManager<TAlgebra> dof_flag_manager; // ONLY ONE INSTANCE for this template param!
	
	return &dof_flag_manager;
}

/* The following functions and classes are implemented for the brige */

///	Adds flags to the dof flag manager
template <typename TAlgebra>
void AddDoFFlag
(
	const char* name ///< name for the new flag block
)
{
	DoFFlagManager<TAlgebra>::get()->add(name);
}

///	Prints flags added to the dof flag manager
template <typename TAlgebra>
void PrintDoFFlags()
{
	DoFFlagManager<TAlgebra>::get()->print();
}

/// Creates flags for a given DoF vector
template <typename TAlgebra>
void CreateDoFFlags
(
	SmartPtr<typename TAlgebra::vector_type> spVec ///< the DoF vector to create the flags for
)
{
	DoFFlagManager<TAlgebra>::get()->create_flags(*spVec);
}

/**
 * The following class is an auxiliary tool to be able to automatically
 * match the calls for the current algebra type in the scripts (like LUA).
 * It does not contain any fields but serves like a function itself.
 */
template <typename TAlgebra>
class DeclareDoFFlags
{
public:

///	Dummy constructor
	DeclareDoFFlags() {};

///	This constructor simply adds a new dof flag to the dof flag manager
	DeclareDoFFlags
	(
		const char* name ///< name for the new flag block
	)
	{
		DoFFlagManager<TAlgebra>::get()->add(name);
	}

///	This constructor simply adds new dof flags to the dof flag manager
	DeclareDoFFlags
	(
		std::vector<std::string> names ///< name for the new flag block
	)
	{
		for(size_t i = 0; i < names.size(); i++)
			DoFFlagManager<TAlgebra>::get()->add(names[i].c_str());
	}

///	Prints the dof flag names
	void print()
	{
		DoFFlagManager<TAlgebra>::get()->print();
	}
};

}
#endif /* __H__UG__LIB_ALGEBRA__DOF_FLAG_MANAGER__ */
