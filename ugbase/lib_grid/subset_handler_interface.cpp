//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d24   (reworked y09 m12 d15)


#include <cassert>
#include "subset_handler_interface.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	SubsetInfo implementation
SubsetInfo::SubsetInfo()
{
	name = "";
	materialIndex = -1;
	subsetState = SS_NONE;
}

////////////////////////////////////////////////////////////////////////
//	ISubsetHandler implementation
ISubsetHandler::
ISubsetHandler(uint supportedElements) : m_aSubsetIndex(false), m_aIterator(false)
{
	m_pGrid = NULL;
	m_supportedElements = supportedElements;
	m_defaultSubsetIndex = -1;
	m_bSubsetInheritanceEnabled = true;
	m_bSubsetAttachmentsEnabled = false;
}
/*
ISubsetHandler::
ISubsetHandler(GridType& grid, uint supportedElements) : m_aSubsetIndex(false), m_aIterator(false)
{
	m_pGrid = NULL;
	m_supportedElements = supportedElements;
	m_defaultSubsetIndex = -1;
	m_bSubsetInheritanceEnabled = true;
	assign_grid(grid);
}
*/
ISubsetHandler::
ISubsetHandler(const ISubsetHandler& sh) : m_aSubsetIndex(false), m_aIterator(false)
{
/*
	m_pGrid = NULL;
	assert(!"WARNING in SubsetHandler::SubsetHandler(const SubsetHandler& sel): Copy-Constructor not yet implemented!");
	LOG("WARNING in SubsetHandler::SubsetHandler(const SubsetHandler& sel): Copy-Constructor not yet implemented! Expect unexpected behaviour!");
*/
	m_pGrid = NULL;
	m_supportedElements = SHE_NONE;
	m_defaultSubsetIndex = -1;
	m_bSubsetInheritanceEnabled = true;
	m_bSubsetAttachmentsEnabled = false;
	
	Grid* pGrid = sh.get_assigned_grid();

	if(pGrid){
		assign_grid(*pGrid);		
		assign_subset_handler(sh);
	}
}

ISubsetHandler::
~ISubsetHandler()
{
	if(m_pGrid != NULL)
		m_pGrid->unregister_observer(this);
}

ISubsetHandler& ISubsetHandler::operator = (const ISubsetHandler& sh)
{
	clear();
	assign_subset_handler(sh);
	return *this;
}

template <class TIterator> static void
CopySubsetIndices(ISubsetHandler& dest, const ISubsetHandler& src,
					TIterator destIterBegin, TIterator destIterEnd,
					TIterator srcIterBegin, TIterator srcIterEnd)
{
	while((srcIterBegin != srcIterEnd) && (destIterBegin != destIterEnd))
	{
		dest.assign_subset(*destIterBegin, src.get_subset_index(*srcIterBegin));
		++destIterBegin;
		++srcIterBegin;
	}
}

void ISubsetHandler::assign_subset_handler(const ISubsetHandler& sh)
{
//	get the source grid
	Grid* srcGrid = sh.get_assigned_grid();
	Grid* destGrid = m_pGrid;
	
//	copy status and options
	set_supported_elements(sh.m_supportedElements);
	enable_subset_inheritance(sh.m_bSubsetInheritanceEnabled);
	set_default_subset_index(sh.m_defaultSubsetIndex);

//TODO: enable attachment support based on the source-hanlders attachment support!?!

//	make sure that both accessors have a valid grid
	if(srcGrid && destGrid){
	//	assign the subsets for each element-type
		CopySubsetIndices(*this, sh, destGrid->begin<VertexBase>(), destGrid->end<VertexBase>(),
							srcGrid->begin<VertexBase>(), srcGrid->end<VertexBase>());

		CopySubsetIndices(*this, sh, destGrid->begin<EdgeBase>(), destGrid->end<EdgeBase>(),
							srcGrid->begin<EdgeBase>(), srcGrid->end<EdgeBase>());

		CopySubsetIndices(*this, sh, destGrid->begin<Face>(), destGrid->end<Face>(),
							srcGrid->begin<Face>(), srcGrid->end<Face>());

		CopySubsetIndices(*this, sh, destGrid->begin<Volume>(), destGrid->end<Volume>(),
							srcGrid->begin<Volume>(), srcGrid->end<Volume>());
							
//TODO:	copy attachments!?!
	}
}

void ISubsetHandler::
create_required_subset_infos(int index)
{
	m_subsetInfos.resize(index+1);

	if(subset_attachments_are_enabled())
		resize_attachment_pipes(index+1);

	add_required_subset_lists(index);
}

void
ISubsetHandler::
assign_grid(Grid& grid)
{
	grid.register_observer(this, OT_GRID_OBSERVER | OT_VERTEX_OBSERVER | OT_EDGE_OBSERVER |
									OT_FACE_OBSERVER | OT_VOLUME_OBSERVER);
}

void ISubsetHandler::
set_default_subset_index(int subsetIndex)
{
	m_defaultSubsetIndex = subsetIndex;
	if(subsetIndex < 0)
		m_defaultSubsetIndex = -1;
}


Grid* ISubsetHandler::
get_assigned_grid() const
{
	return m_pGrid;
}

bool ISubsetHandler::
elements_are_supported(uint shElements)
{
	return (m_supportedElements & shElements) == shElements;
}

void ISubsetHandler::
set_supported_elements(uint shElements)
{
//	do this in two steps:
//	1: disable the element-support that is no longer required.
//	2: enable the element-support that was not already enabled.
	
//	disable the elements that shall be disabled.
//	(the ones which shall not be set, but are currently active.)
	disable_element_support((!shElements) & m_supportedElements);

//	enable the elements that are not already enabled
	enable_element_support(shElements & (!m_supportedElements));
}

void ISubsetHandler::
enable_element_support(uint shElements)
{
//	if no grid is assigned, we can't do anything.
	if(m_pGrid)
	{
	//	check for each option whether it should be enabled.
	//	to reduce unnecessary operations, we have to make sure that
	//	that option hasn't already been enabled.

		if((shElements & SHE_VERTEX) &&
			(!elements_are_supported(SHE_VERTEX)))
		{
//LOG("enabling vertex support\n");
		//	enable vertex-support.
			m_pGrid->attach_to_vertices(m_aSubsetIndex);
			m_pGrid->attach_to_vertices(m_aIterator);
			m_aaSubsetIndexVRT.access(*m_pGrid, m_aSubsetIndex);
			m_aaIteratorVRT.access(*m_pGrid, m_aIterator);
			m_supportedElements |= SHE_VERTEX;
			reset_subset_indices(SHE_VERTEX);
		}

		if((shElements & SHE_EDGE) &&
			(!elements_are_supported(SHE_EDGE)))
		{
//LOG("enabling edge support\n");
		//	enable edge support
			m_pGrid->attach_to_edges(m_aSubsetIndex);
			m_pGrid->attach_to_edges(m_aIterator);
			m_aaSubsetIndexEDGE.access(*m_pGrid, m_aSubsetIndex);
			m_aaIteratorEDGE.access(*m_pGrid, m_aIterator);
			m_supportedElements |= SHE_EDGE;
			reset_subset_indices(SHE_EDGE);
		}

		if((shElements & SHE_FACE) &&
			(!elements_are_supported(SHE_FACE)))
		{
//LOG("enabling face support\n");
		//	enable face support
			m_pGrid->attach_to_faces(m_aSubsetIndex);
			m_pGrid->attach_to_faces(m_aIterator);
			m_aaSubsetIndexFACE.access(*m_pGrid, m_aSubsetIndex);
			m_aaIteratorFACE.access(*m_pGrid, m_aIterator);
			m_supportedElements |= SHE_FACE;
			reset_subset_indices(SHE_FACE);
		}

		if((shElements & SHE_VOLUME) &&
			(!elements_are_supported(SHE_VOLUME)))
		{
//LOG("enabling volume support\n");
		//	enable volume support
			m_pGrid->attach_to_volumes(m_aSubsetIndex);
			m_pGrid->attach_to_volumes(m_aIterator);
			m_aaSubsetIndexVOL.access(*m_pGrid, m_aSubsetIndex);
			m_aaIteratorVOL.access(*m_pGrid, m_aIterator);
			m_supportedElements |= SHE_VOLUME;
			reset_subset_indices(SHE_VOLUME);
		}
	}
}

void ISubsetHandler::
disable_element_support(uint shElements)
{
//	if no grid is assigned, we can't do anything.
	if(m_pGrid)
	{
	//	check for each option whether it should be disabled.
	//	to reduce unnecessary operations, we have to make sure that
	//	that option hasn't already been disabled.

		if((shElements & SHE_VERTEX) && elements_are_supported(SHE_VERTEX))
		{
//LOG("disabling vertex support\n");
			m_pGrid->detach_from_vertices(m_aSubsetIndex);
			m_pGrid->detach_from_vertices(m_aIterator);
		}

		if((shElements & SHE_EDGE) && elements_are_supported(SHE_EDGE))
		{
//LOG("disabling edge support\n");
			m_pGrid->detach_from_edges(m_aSubsetIndex);
			m_pGrid->detach_from_edges(m_aIterator);
		}

		if((shElements & SHE_FACE) && elements_are_supported(SHE_FACE))
		{
//LOG("disabling face support\n");
			m_pGrid->detach_from_faces(m_aSubsetIndex);
			m_pGrid->detach_from_faces(m_aIterator);
		}

		if((shElements & SHE_VOLUME) && elements_are_supported(SHE_VOLUME))
		{
//LOG("disabling volume support\n");
			m_pGrid->detach_from_volumes(m_aSubsetIndex);
			m_pGrid->detach_from_volumes(m_aIterator);
		}
	}

//	remove the disabled elements from the set of currently supported elements.
	m_supportedElements &= (!shElements);
}

void ISubsetHandler::
enable_subset_inheritance(bool bEnable)
{
	m_bSubsetInheritanceEnabled = bEnable;
}

void ISubsetHandler::
set_subset_info(int subsetIndex, const SubsetInfo& subsetInfo)
{
	subset_info_required(subsetIndex);
	m_subsetInfos[subsetIndex] = subsetInfo;
}

SubsetInfo& ISubsetHandler::
subset_info(int subsetIndex)
{
	subset_info_required(subsetIndex);
	return m_subsetInfos[subsetIndex];
}

 
const SubsetInfo& ISubsetHandler::
subset_info(int subsetIndex) const
{
	assert(((subsetIndex >= 0) && (subsetIndex < (int)num_subset_infos())) && "ERROR in SubsetHandler::subset_info(..) const: bad subset index. Use non-const version to avoid this Problem.");

	return m_subsetInfos[subsetIndex];
}

void ISubsetHandler::
clear()
{
//	erase subsets
	erase_subset_lists();

//	erase subset indfos
	m_subsetInfos.clear();

//	reset all element subset-indices.
	reset_subset_indices();
}

void ISubsetHandler::
clear_subset(int subsetIndex)
{
	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subset_infos()) &&
			"ERROR in SubsetHandler::clear_subset(): bad subset index.");

	m_subsetInfos[subsetIndex] = SubsetInfo();

	change_subset_indices(subsetIndex, -1);
	clear_subset_lists(subsetIndex);
}

void ISubsetHandler::
clear_subsets()
{
	for(int i = 0; i < (int)num_subset_infos(); ++i)
		clear_subset(i);
}


template <class TElem, class TAAInd>
static void
ResetSubsetIndices(Grid* pGrid, TAAInd& aaInd)
{
	typedef typename geometry_traits<TElem>::iterator iterator;
//	in the given subset to newInd.
	for(iterator iter = pGrid->begin<TElem>();
		iter != pGrid->end<TElem>(); iter++)
		aaInd[*iter] = -1;
}

void ISubsetHandler::
reset_subset_indices(uint shElements)
{
	if(m_pGrid)
	{
		if((shElements & SHE_VERTEX) && elements_are_supported(SHE_VERTEX))
			ResetSubsetIndices<VertexBase>(m_pGrid, m_aaSubsetIndexVRT);
		if((shElements & SHE_EDGE) && elements_are_supported(SHE_EDGE))
			ResetSubsetIndices<EdgeBase>(m_pGrid, m_aaSubsetIndexEDGE);
		if((shElements & SHE_FACE) && elements_are_supported(SHE_FACE))
			ResetSubsetIndices<Face>(m_pGrid, m_aaSubsetIndexFACE);
		if((shElements & SHE_VOLUME) && elements_are_supported(SHE_VOLUME))
			ResetSubsetIndices<Volume>(m_pGrid, m_aaSubsetIndexVOL);
	}
}

int ISubsetHandler::
get_subset_index(GeometricObject* elem) const
{
	uint type = elem->base_object_type_id();
	switch(type)
	{
		case VERTEX:
		{
			VertexBase* pVrt = dynamic_cast<VertexBase*>(elem);
			assert((pVrt != NULL) && "ERROR in SubsetHandler::get_subset_index(GeometricObject* elem): elem type and base_type_id do not match!");
			return get_subset_index(pVrt);
		}
		case EDGE:
		{
			EdgeBase* pEdge = dynamic_cast<EdgeBase*>(elem);
			assert((pEdge != NULL) && "ERROR in SubsetHandler::get_subset_index(GeometricObject* elem): elem type and base_type_id do not match!");
			return get_subset_index(pEdge);
		}
		case FACE:
		{
			Face* pFace = dynamic_cast<Face*>(elem);
			assert((pFace != NULL) && "ERROR in SubsetHandler::get_subset_index(GeometricObject* elem): elem type and base_type_id do not match!");
			return get_subset_index(pFace);
		}
		case VOLUME:
		{
			Volume* pVol = dynamic_cast<Volume*>(elem);
			assert((pVol != NULL) && "ERROR in SubsetHandler::get_subset_index(GeometricObject* elem): elem type and base_type_id do not match!");
			return get_subset_index(pVol);
		}
	}

//	we should never arrive at this point
	assert(!"ERROR in SubsetHandler::get_subset_index(GeometricObject* elem): Program should never reach this point!");
	return -1;
}

void ISubsetHandler::
insert_subset(int subsetIndex)
{
//	append required subsets
	if(subsetIndex >= 0)
	{
	//	make sure that the subset-infos
		subset_info_required(num_subset_infos());
		if(subsetIndex < (int)num_subset_infos()-1)
			move_subset(num_subset_infos() - 1, subsetIndex);
	}
}

void ISubsetHandler::
erase_subset(int subsetIndex)
{
//	assign all elements of this subset to -1
//	delete the subset
//	move all subsets with higher index one entry up
//	and correct indices of assigned elements.
	if((subsetIndex >= 0) && (subsetIndex < (int)num_subset_infos()))
	{
		change_subset_indices(subsetIndex, -1);

	//	clear and delete pipes of erased subset
		if(subset_attachments_are_enabled())
		{
			clear_attachment_pipes(subsetIndex);
			delete m_vertexAttachmentPipes[subsetIndex];
			delete m_edgeAttachmentPipes[subsetIndex];
			delete m_faceAttachmentPipes[subsetIndex];
			delete m_volumeAttachmentPipes[subsetIndex];
		}

		for(uint i = subsetIndex + 1; i < num_subset_infos(); ++i)
		{
		//	alter indices
			change_subset_indices(i, i-1);

		//	move the subset
			m_subsetInfos[i-1] = m_subsetInfos[i];
			
		//	move the pipes
			if(subset_attachments_are_enabled())
			{
				m_vertexAttachmentPipes[i-1] = m_vertexAttachmentPipes[i];
				m_edgeAttachmentPipes[i-1] = m_edgeAttachmentPipes[i];
				m_faceAttachmentPipes[i-1] = m_faceAttachmentPipes[i];
				m_volumeAttachmentPipes[i-1] = m_volumeAttachmentPipes[i];
			}
		}

	//	resize the subset vector
		uint numNewSubsets = num_subset_infos() - 1;
		erase_subset_lists(subsetIndex);
		m_subsetInfos.resize(numNewSubsets);
		resize_attachment_pipes(numNewSubsets);
	}
	else
		LOG("WARNING in SubsetHandler::erase_subset(...): bad subset index: " << subsetIndex << endl);
}

void ISubsetHandler::
swap_subsets(int subsetIndex1, int subsetIndex2)
{
	if((subsetIndex1 >= 0) && (subsetIndex1 < (int)num_subset_infos())
		&& (subsetIndex2 >= 0) && (subsetIndex2 < (int)num_subset_infos()))
	{
	//	set all indices of subset1 to subsetIndex2 and vice versa.
	//	swap pointers of subsets afterwards
		change_subset_indices(subsetIndex1, subsetIndex2);
		change_subset_indices(subsetIndex2, subsetIndex1);

	//	swap pointers
		SubsetInfo tmpSI = m_subsetInfos[subsetIndex1];
		m_subsetInfos[subsetIndex1] = m_subsetInfos[subsetIndex2];
		m_subsetInfos[subsetIndex2] = tmpSI;

	//	store from-attachment-pipes
		if(subset_attachments_are_enabled())
		{
			VertexAttachmentPipe* apFromVrt = m_vertexAttachmentPipes[subsetIndex1];
			m_vertexAttachmentPipes[subsetIndex1] = m_vertexAttachmentPipes[subsetIndex2];
			m_vertexAttachmentPipes[subsetIndex2] = apFromVrt;
			EdgeAttachmentPipe* apFromEdge = m_edgeAttachmentPipes[subsetIndex1];
			m_edgeAttachmentPipes[subsetIndex1] = m_edgeAttachmentPipes[subsetIndex2];
			m_edgeAttachmentPipes[subsetIndex2] = apFromEdge;
			FaceAttachmentPipe* apFromFace = m_faceAttachmentPipes[subsetIndex1];
			m_faceAttachmentPipes[subsetIndex1] = m_faceAttachmentPipes[subsetIndex2];
			m_faceAttachmentPipes[subsetIndex2] = apFromFace;
			VolumeAttachmentPipe* apFromVol = m_volumeAttachmentPipes[subsetIndex1];
			m_volumeAttachmentPipes[subsetIndex1] = m_volumeAttachmentPipes[subsetIndex2];
			m_volumeAttachmentPipes[subsetIndex2] = apFromVol;
		}

	//	swap the lists
		swap_subset_lists(subsetIndex1, subsetIndex2);
	}
	else
		LOG("WARNING in SubsetHandler::swap_subsets(...): bad indices: " << subsetIndex1 << ", " << subsetIndex2 << endl);
}

void ISubsetHandler::
move_subset(int indexFrom, int indexTo)
{
	if((indexFrom >= 0) && (indexFrom < (int)num_subset_infos()) && (indexTo >= 0))
	{
		int moveDir = 0;
	//	we have to distinguish two cases
		if(indexFrom < indexTo)
			moveDir = 1;
		else if(indexTo < indexFrom)
			moveDir = -1;

		if(moveDir != 0)
		{
		//	check if we have to append subsets
			subset_info_required(indexTo);

		//	store the from-subset-info
			SubsetInfo siFrom = m_subsetInfos[indexFrom];

		//	store from-attachment-pipes
			VertexAttachmentPipe* apFromVrt;
			EdgeAttachmentPipe* apFromEdge;
			FaceAttachmentPipe* apFromFace;
			VolumeAttachmentPipe* apFromVol;
			if(subset_attachments_are_enabled())
			{
				apFromVrt = m_vertexAttachmentPipes[indexFrom];
				apFromEdge = m_edgeAttachmentPipes[indexFrom];
				apFromFace = m_faceAttachmentPipes[indexFrom];
				apFromVol = m_volumeAttachmentPipes[indexFrom];
			}

		//	assign new indices to elements in from subset (no iterators are changed)
			change_subset_indices(indexFrom, indexTo);

		//	move the subsets between indexFrom and indexTo
			for(int i = indexFrom; i != indexTo; i+= moveDir)
			{
				int iNext = i+moveDir;
			//	assign new indices
				change_subset_indices(iNext, i);

			//	move info
				m_subsetInfos[i] = m_subsetInfos[iNext];
			}

		//	assign stored info
			m_subsetInfos[indexTo] = siFrom;
			
		//	assign stored attachment pipes
			if(subset_attachments_are_enabled())
			{
				m_vertexAttachmentPipes[indexTo] = apFromVrt;
				m_edgeAttachmentPipes[indexTo] = apFromEdge;
				m_faceAttachmentPipes[indexTo] = apFromFace;
				m_volumeAttachmentPipes[indexTo] = apFromVol;
			}

		//	move the subsets lists
			move_subset_lists(indexFrom, indexTo);
		}
	}
	else
		LOG("WARNING in SubsetHandler::move_subset(...): bad indices: " << indexFrom << ", " << indexTo << endl);
}

////////////////////////////////////////////////////////////////////////
//	attachments

void ISubsetHandler::
resize_attachment_pipes(size_t newSize)
{
	while(m_vertexAttachmentPipes.size() < newSize)
		m_vertexAttachmentPipes.push_back(new VertexAttachmentPipe(this));
	while(m_edgeAttachmentPipes.size() < newSize)
		m_edgeAttachmentPipes.push_back(new EdgeAttachmentPipe(this));
	while(m_faceAttachmentPipes.size() < newSize)
		m_faceAttachmentPipes.push_back(new FaceAttachmentPipe(this));
	while(m_volumeAttachmentPipes.size() < newSize)
		m_volumeAttachmentPipes.push_back(new VolumeAttachmentPipe(this));
}

void ISubsetHandler::
clear_attachment_pipes()
{
	for(size_t i = 0; i < m_vertexAttachmentPipes.size(); ++i){
	//	all pipes have the same size
		delete m_vertexAttachmentPipes[i];
		delete m_edgeAttachmentPipes[i];
		delete m_faceAttachmentPipes[i];
		delete m_volumeAttachmentPipes[i];
	}
	m_vertexAttachmentPipes.clear();
	m_edgeAttachmentPipes.clear();
	m_faceAttachmentPipes.clear();
	m_volumeAttachmentPipes.clear();
}

void ISubsetHandler::
clear_attachment_pipes(int subsetIndex)
{
	m_vertexAttachmentPipes[subsetIndex]->clear();
	m_edgeAttachmentPipes[subsetIndex]->clear();
	m_faceAttachmentPipes[subsetIndex]->clear();
	m_volumeAttachmentPipes[subsetIndex]->clear();
}

void ISubsetHandler::
enable_subset_attachments(bool bEnable)
{
	if(bEnable &! subset_attachments_are_enabled())
	{
		LOG("  enabling...\n");
	//	enable subset-attachments.
		m_bSubsetAttachmentsEnabled = true;
	//	attach data-indices to the elements
//TODO:	allow subset-attachments for vertices, edges, faces and volumes separately.
		if(m_pGrid)
		{
			LOG("  attaching data...\n");
			m_pGrid->attach_to_vertices(m_aDataIndex);
			m_pGrid->attach_to_edges(m_aDataIndex);
			m_pGrid->attach_to_faces(m_aDataIndex);
			m_pGrid->attach_to_volumes(m_aDataIndex);
			m_aaDataIndVRT.access(*m_pGrid, m_aDataIndex);
			m_aaDataIndEDGE.access(*m_pGrid, m_aDataIndex);
			m_aaDataIndFACE.access(*m_pGrid, m_aDataIndex);
			m_aaDataIndVOL.access(*m_pGrid, m_aDataIndex);
			
			resize_attachment_pipes(num_subset_infos());
		//	tell the derived class that it should register all
		//	of its elements at the pipe.
			register_subset_elements_at_pipe();
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	grid callbacks
void ISubsetHandler::
registered_at_grid(Grid* grid)
{
	if(m_pGrid != NULL)
		m_pGrid->unregister_observer(this);

	m_pGrid = grid;
	if(m_pGrid != NULL)
	{
	//	initialise attachments and accessors.
	//	do this whith a little trick:
	//	set the supported-element-options to SHE_NONE,
	//	then call enable for all element-types that should be supported.
		uint tmpOpts = m_supportedElements;
		m_supportedElements = SHE_NONE;
		enable_element_support(tmpOpts);

	//	enable attachments - if required
		if(subset_attachments_are_enabled())
		{
			m_bSubsetAttachmentsEnabled = false;
			enable_subset_attachments(true);
		}

	//DEBUG: log m_supportedElements
		//LOG("supported elements after registration: " << m_supportedElements << "\n");
	}
}

void ISubsetHandler::
unregistered_from_grid(Grid* grid)
{
	assert(m_pGrid && "this method should only be called while the handler is registered at a grid.");
	
	if(m_pGrid)
	{
		assert((m_pGrid == grid) && "ERROR in SubsetHandler::unregistered_from_grid(...): Grids do not match.");

	//	clear the subsets
		m_subsetInfos.clear();

	//	disable all currently supported elements (this will remove any attachments)
		disable_element_support(m_supportedElements);

	//	clear attachment data
		if(subset_attachments_are_enabled())
		{
			m_pGrid->detach_from_vertices(m_aDataIndex);
			m_pGrid->detach_from_edges(m_aDataIndex);
			m_pGrid->detach_from_faces(m_aDataIndex);
			m_pGrid->detach_from_volumes(m_aDataIndex);
			
		//	clear pipes
			clear_attachment_pipes();
		}
	
	//DEBUG: log m_supportedElements
		//LOG("supported elements after deregistration: " << m_supportedElements << "\n");
		m_pGrid = NULL;
	}
}

void ISubsetHandler::
elements_to_be_cleared(Grid* grid)
{
	for(uint si = 0; si < num_subset_infos(); si++)
		clear_subset_lists(si);
}

//	vertex callbacks
void ISubsetHandler::
vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent)
{
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::vertex_created(...): Grids do not match.");
//LOG("new vertex...\n");
	m_aaSubsetIndexVRT[vrt] = -1;
//LOG("si_before assignement: " << get_subset_index(vrt) << endl);
	if((pParent != NULL) && m_bSubsetInheritanceEnabled)
		assign_subset(vrt, get_subset_index(pParent));
	else if(m_defaultSubsetIndex != -1)
		assign_subset(vrt, m_defaultSubsetIndex);
//LOG("vertex creation done\n");
}

void ISubsetHandler::
vertex_to_be_erased(Grid* grid, VertexBase* vrt)
{
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::vertex_to_be_erased(...): Grids do not match.");
	if(m_aaSubsetIndexVRT[vrt] != -1)
		assign_subset(vrt, -1);

}

//	edge callbacks
void ISubsetHandler::
edge_created(Grid* grid, EdgeBase* edge, GeometricObject* pParent)
{
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::edge_created(...): Grids do not match.");
	m_aaSubsetIndexEDGE[edge] = -1;

	if((pParent != NULL) && m_bSubsetInheritanceEnabled)
		assign_subset(edge, get_subset_index(pParent));
	else if(m_defaultSubsetIndex != -1)
		assign_subset(edge, m_defaultSubsetIndex);
}

void ISubsetHandler::
edge_to_be_erased(Grid* grid, EdgeBase* edge)
{
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::edge_to_be_erased(...): Grids do not match.");
	if(m_aaSubsetIndexEDGE[edge] != -1)
		assign_subset(edge, -1);

}

//	face callbacks
void ISubsetHandler::
face_created(Grid* grid, Face* face, GeometricObject* pParent)
{
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::face_created(...): Grids do not match.");
	m_aaSubsetIndexFACE[face] = -1;

	if((pParent != NULL) && m_bSubsetInheritanceEnabled)
		assign_subset(face, get_subset_index(pParent));
	else if(m_defaultSubsetIndex != -1)
		assign_subset(face, m_defaultSubsetIndex);
}

void ISubsetHandler::
face_to_be_erased(Grid* grid, Face* face)
{
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::face_to_be_erased(...): Grids do not match.");
	if(m_aaSubsetIndexFACE[face] != -1)
		assign_subset(face, -1);
}

//	volume callbacks
void ISubsetHandler::
volume_created(Grid* grid, Volume* vol, GeometricObject* pParent)
{
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::volume_created(...): Grids do not match.");
	m_aaSubsetIndexVOL[vol] = -1;

	if((pParent != NULL) && m_bSubsetInheritanceEnabled)
		assign_subset(vol, get_subset_index(pParent));
	else if(m_defaultSubsetIndex != -1)
		assign_subset(vol, m_defaultSubsetIndex);
}

void ISubsetHandler::
volume_to_be_erased(Grid* grid, Volume* vol)
{
	assert((m_pGrid == grid) && "ERROR in SubsetHandler::volume_to_be_erased(...): Grids do not match.");
	if(m_aaSubsetIndexVOL[vol] != -1)
		assign_subset(vol, -1);
}

}//	end of namespace

