#include "./p1conform.h"


namespace ug{

///////////////////////////////////////
// P1StorageManager
///////////////////////////////////////

void
P1StorageManager::
set_subset_handler(ISubsetHandler& sh)
{
	if(m_pSH != NULL) clear();

	m_pSH = &sh;
	m_pSH->enable_subset_attachments(true);
}

void
P1StorageManager::
clear_subset_handler()
{
	if(m_pSH != NULL) clear();
	m_pSH->enable_subset_attachments(false);
	m_pSH = NULL;
}

void
P1StorageManager::
clear()
{
	if(m_pSH == NULL) return;

	for(size_t si = 0; si < m_vSubsetInfo.size(); ++si)
	{
		m_pSH->detach_from<VertexBase>(m_vSubsetInfo[si].aDoF, si);
		m_vSubsetInfo[si].aaDoFVRT.invalidate();
	}
	m_vSubsetInfo.clear();
}

void
P1StorageManager::
update_attachments()
{
	if(m_pSH == NULL) return;

	size_t num_subsets =  m_pSH->num_subsets();

	// Create level dof distributors
	for(size_t si = m_vSubsetInfo.size(); si < num_subsets; ++si)
	{
		m_vSubsetInfo.push_back(SubsetInfo());
		m_pSH->attach_to<VertexBase>(m_vSubsetInfo[si].aDoF, si);
		m_vSubsetInfo[si].aaDoFVRT.access(*m_pSH, m_vSubsetInfo[si].aDoF, si);
	}
}


///////////////////////////////////////
// P1ConformFunctionPattern
///////////////////////////////////////

bool
P1ConformFunctionPattern::
add_discrete_function(const char* name, LocalShapeFunctionSetID id, int dim)
{
	// for a P1 dof manager only Lagrange P1 function space is permitted
	if(id != LSFS_LAGRANGEP1)
		{UG_LOG("P1ConformDoFDistributor: Only LSFS_LAGRANGEP1 functions are supported.\n"); return false;}

	return FunctionPattern::add_discrete_function(name, id, dim);
}

bool
P1ConformFunctionPattern::
add_discrete_function(const char* name, LocalShapeFunctionSetID id, const SubsetGroup& SubsetIndices, int dim)
{
	// for a P1 dof manager only Lagrange P1 function space is permitted
	if(id != LSFS_LAGRANGEP1)
		{UG_LOG("P1ConformDoFDistributor: Only LSFS_LAGRANGEP1 functions are supported.\n"); return false;}

	return FunctionPattern::add_discrete_function(name, id, SubsetIndices, dim);
}


///////////////////////////////////////
// P1ConformDoFDistribution
///////////////////////////////////////

size_t
P1ConformDoFDistribution::
num_indices(ReferenceObjectID refID, int si, const FunctionGroup& funcGroup) const
{
	const ReferenceElement& refElem = ReferenceElementFactory::get_reference_element(refID);

	size_t numFct = 0;
	for(size_t fct = 0; fct < funcGroup.num_fct(); ++fct)
	{
		if(is_def_in_subset(funcGroup[fct], si))
			numFct++;
	}

	return numFct * refElem.num_obj(0);
}

bool
P1ConformDoFDistribution::
prepare_indices(ReferenceObjectID refID, int si, LocalIndices& ind, bool withHanging) const
{
	const ReferenceElement& refElem = ReferenceElementFactory::get_reference_element(refID);

	if(!withHanging)
	{
		ind.clear();
		size_t numInd = 0;
		size_t numFct = 0;
		for(size_t fct = 0; fct < ind.num_fct(); ++fct)
		{
			if(!is_def_in_subset(ind.fct_id(fct), si)) continue;
			for(size_t dof = 0; dof < refElem.num_obj(0); ++dof)
			{
				LocalIndices::multi_index_type dof_ind;
				dof_ind[0] = dof + numFct * refElem.num_obj(0);
				dof_ind[1] = 0;
				ind.add_dof(fct, dof_ind);
			}
			numInd += refElem.num_obj(0);
			numFct++;
		}
		ind.set_num_indices(numInd);
	}
	return true;
}


size_t
P1ConformDoFDistribution::
num_inner_indices(ReferenceObjectID refID, int si, const FunctionGroup& funcGroup) const
{
	if(refID != ROID_VERTEX) return 0;
	else return num_indices(refID, si, funcGroup);
}

bool
P1ConformDoFDistribution::
prepare_inner_indices(ReferenceObjectID refID, int si, LocalIndices& ind) const
{
	ind.clear();
	if(refID != ROID_VERTEX) return true;
	else return prepare_indices(refID, si, ind);
}

///////////// creation /////////////////

bool
P1ConformDoFDistribution::
distribute_dofs()
{
	if(m_pStorageManager == NULL)
	{
		UG_LOG("In 'P1ConformDoFDistribution::distribute_dofs:"
				"Storage Manager not set. Aborting.\n");
		return false;
	}

	if(m_pFunctionPattern == NULL)
	{
		UG_LOG("In 'P1ConformDoFDistribution::distribute_dofs:"
				"Function Pattern not set. Aborting.\n");
		return false;
	}

	// Attach indices
	m_pStorageManager->update_attachments();

	// iterators
	geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

	// loop subsets
	size_t i = 0;
	size_t i_per_subset = 0;

	// reset number of dofs
	m_vNumDoFs.clear(); m_vNumDoFs.resize(num_subsets(), 0);
	m_vvOffsets.clear(); m_vvOffsets.resize(num_subsets());

	// loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
		// create offsets
		size_t count = 0;
		m_vvOffsets[si].resize(m_pFunctionPattern->num_fct());
		for(size_t fct = 0; fct < m_pFunctionPattern->num_fct(); ++fct)
		{
			if(!is_def_in_subset(fct, si)) m_vvOffsets[si][fct] = (size_t) -1;
			else m_vvOffsets[si][fct] = count++;
		}

		// skip if no dofs to be distributed
		if(!(m_pFunctionPattern->num_fct(si)>0)) continue;

		iterBegin = m_goc.begin<VertexBase>(si);
		iterEnd =  m_goc.end<VertexBase>(si);

		// remember number of functions
		size_t num_fct =  m_pFunctionPattern->num_fct(si);

		// loop Verices
		i_per_subset = 0;
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
			// get vertex
			VertexBase* vrt = *iter;

			// write index
			m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] = i;
			i += num_fct;
			i_per_subset += num_fct;
		}
		m_vNumDoFs[si] = i_per_subset;
	}
	m_numDoFs = i;

	UG_LOG(std::setw(8) << m_numDoFs <<" |     " << 1 << "     | " );

	for(int si = 0; si < num_subsets(); ++si)
	{
		UG_LOG( " (" << si << ","<<1<<"," << std::setw(8) << m_vNumDoFs[si] << ") ");
	}
	UG_LOG(std::endl);

	return true;
}

///////////// Help functions /////////////////
VertexBase* P1ConformDoFDistribution::get_vertex(VertexBase* vrt, size_t i) const
{
	UG_ASSERT(i < 1, "A Vertex has only one vertex");
	return vrt;
}

VertexBase* P1ConformDoFDistribution::get_vertex(EdgeBase* edge, size_t i) const
{
	UG_ASSERT(i < edge->num_vertices(), "Wrong number of vertex");
	return edge->vertex(i);
}

VertexBase* P1ConformDoFDistribution::get_vertex(Face* face, size_t i) const
{
	UG_ASSERT(i < face->num_vertices(), "Wrong number of vertex");
	return face->vertex(i);
}

VertexBase* P1ConformDoFDistribution::get_vertex(Volume* vol, size_t i) const
{
	UG_ASSERT(i < vol->num_vertices(), "Wrong number of vertex");
	return vol->vertex(i);
}


///////////////////////////////////////
// GroupedP1ConformDoFDistribution
///////////////////////////////////////


size_t
GroupedP1ConformDoFDistribution::
num_indices(ReferenceObjectID refID, int si, const FunctionGroup& funcGroup) const
{
	const ReferenceElement& refElem = ReferenceElementFactory::get_reference_element(refID);

	for(size_t fct = 0; fct < funcGroup.num_fct(); ++fct)
	{
		if(is_def_in_subset(funcGroup[fct], si))
			return refElem.num_obj(0);
	}

	return 0;
}


size_t
GroupedP1ConformDoFDistribution::
num_inner_indices(ReferenceObjectID refID, int si, const FunctionGroup& funcGroup) const
{
	if(refID != ROID_VERTEX) return 0;
	else return num_indices(refID, si, funcGroup);
}

bool
GroupedP1ConformDoFDistribution::
prepare_indices(ReferenceObjectID refID, int si, LocalIndices& ind, bool withHanging) const
{
	const ReferenceElement& refElem = ReferenceElementFactory::get_reference_element(refID);

	ind.clear();
	size_t numInd = 0;
	for(size_t fct = 0; fct < ind.num_fct(); ++fct)
	{
		if(!is_def_in_subset(ind.fct_id(fct), si)) continue;

		for(size_t dof = 0; dof < refElem.num_obj(0); ++dof)
		{
			LocalIndices::multi_index_type dof_ind;
			dof_ind[0] = dof;
			dof_ind[1] = ind.fct_id(fct);
			ind.add_dof(fct, dof_ind);
		}
		numInd = refElem.num_obj(0);
	}
	ind.set_num_indices(numInd);
	return true;
}

bool
GroupedP1ConformDoFDistribution::
prepare_inner_indices(ReferenceObjectID refID, int si, LocalIndices& ind) const
{
	ind.clear();
	if(refID != ROID_VERTEX) return 0;
	else return prepare_indices(refID, si, ind);
}

///////// Creation /////////////

bool
GroupedP1ConformDoFDistribution::
distribute_dofs()
{
	// Attach indices
	m_pStorageManager->update_attachments();

	// iterators
	geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

	// loop subsets
	size_t i = 0;

	// reset number of dofs
	m_vNumDoFs.clear(); m_vNumDoFs.resize(num_subsets());

	// loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
		// skip if no dofs to be distributed
		if(!(m_pFunctionPattern->num_fct(si)>0)) continue;

		iterBegin = m_goc.begin<VertexBase>(si);
		iterEnd =  m_goc.end<VertexBase>(si);

		// loop Verices
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
			// get vertex
			VertexBase* vrt = *iter;

			// write index
			m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] = i;
			i ++;
		}
		if(si==0) m_vNumDoFs[si] = i;
		else m_vNumDoFs[si] = i - m_vNumDoFs[si-1];
	}
	m_numDoFs = i;

	UG_LOG(std::setw(8) << m_numDoFs <<" | " << "variable" << "  | " );

	for(int si = 0; si < num_subsets(); ++si)
	{
		size_t num_fct =  m_pFunctionPattern->num_fct(si);
		UG_LOG( " (" << si << "," <<num_fct<<","<< std::setw(8) << m_vNumDoFs[si] << ") ");
	}
	UG_LOG(std::endl);

	return true;
}

///////////// Help functions /////////////////
VertexBase* GroupedP1ConformDoFDistribution::get_vertex(VertexBase* vrt, size_t i) const
{
	UG_ASSERT(i < 1, "A Vertex has only one vertex");
	return vrt;
}

VertexBase* GroupedP1ConformDoFDistribution::get_vertex(EdgeBase* edge, size_t i) const
{
	UG_ASSERT(i < edge->num_vertices(), "Wrong number of vertex");
	return edge->vertex(i);
}

VertexBase* GroupedP1ConformDoFDistribution::get_vertex(Face* face, size_t i) const
{
	UG_ASSERT(i < face->num_vertices(), "Wrong number of vertex");
	return face->vertex(i);
}

VertexBase* GroupedP1ConformDoFDistribution::get_vertex(Volume* vol, size_t i) const
{
	UG_ASSERT(i < vol->num_vertices(), "Wrong number of vertex");
	return vol->vertex(i);
}





} // end namespace ug
