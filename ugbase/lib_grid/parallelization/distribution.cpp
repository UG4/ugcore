// created by Sebastian Reiter
// s.b.reiter@gmail.com
// 29.11.2012 (d,m,y)

#include <sstream>
#include "common/static_assert.h"
#include "common/util/table.h"
#include "distribution.h"
#include "distributed_grid.h"
#include "lib_grid/tools/selector_multi_grid.h"
#include "lib_grid/algorithms/selection_util.h"
#include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "parallelization_util.h"
#include "lib_grid/file_io/file_io.h"

//#define LG_DISTRIBUTION_DEBUG
//#define LG_DISTRIBUTION_Z_OUTPUT_TRANSFORM 40


using namespace std;

namespace ug{

enum InterfaceStates{
	IS_UNASSIGNED = 0,
	IS_NORMAL = 1,
	IS_VMASTER = 1<<1,
	IS_VSLAVE = 1<<2,
	IS_DUMMY = 1<<3,
	HAS_PARENT = 1<<4 // only used for dummies currently
};


struct TargetProcInfo
{
	TargetProcInfo() {}
	TargetProcInfo(int pID, byte intfcState) :
		procID(pID), interfaceState(intfcState) {}

	int procID;
	byte interfaceState; // or-combinations of constants from InterfaceStates
};

typedef Attachment<vector<TargetProcInfo> >	ADistInfo;

///	Automatically attaches ADistInfo to all elements of a grid.
/**	On destruction, the attached dist info attachments are removed again.
 * You may access the dist-info in an element through the get method.
 * The get method returns a reference to the attached DistInfo object.
 *
 * Make sure that the given grid is valid while the DistInfoSupplier exists.
 */
class DistInfoSupplier{
	public:
		DistInfoSupplier(Grid& grid) : m_grid(grid), m_aDistInfo("distribution-info")
		{
			m_grid.attach_to_all(m_aDistInfo);
			m_aaDistInfoVRT.access(grid, m_aDistInfo);
			m_aaDistInfoEDGE.access(grid, m_aDistInfo);
			m_aaDistInfoFACE.access(grid, m_aDistInfo);
			m_aaDistInfoVOL.access(grid, m_aDistInfo);
		}

		~DistInfoSupplier()
		{
			m_grid.detach_from_all(m_aDistInfo);
		}

		vector<TargetProcInfo>& get(VertexBase* vrt)	{return m_aaDistInfoVRT[vrt];}
		vector<TargetProcInfo>& get(EdgeBase* edge)		{return m_aaDistInfoEDGE[edge];}
		vector<TargetProcInfo>& get(Face* face)			{return m_aaDistInfoFACE[face];}
		vector<TargetProcInfo>& get(Volume* vol)		{return m_aaDistInfoVOL[vol];}
		vector<TargetProcInfo>& get(GeometricObject* obj)
		{
			int objType = obj->base_object_id();
			switch(objType){
				case VERTEX:	return get(static_cast<VertexBase*>(obj));
				case EDGE:		return get(static_cast<EdgeBase*>(obj));
				case FACE:		return get(static_cast<Face*>(obj));
				case VOLUME:	return get(static_cast<Volume*>(obj));
				default:	UG_THROW("Unknown geometric object base type."); break;
			}
		}

		ADistInfo dist_info_attachment()	{return m_aDistInfo;}

		template <class TElem>
		StringStreamTable get_debug_info(TElem* e)
		{
			const vector<TargetProcInfo>& di = get(e);
			StringStreamTable t;
			for(size_t i = 0; i < di.size(); ++i){
				t(0, i+1) << "p" << di[i].procID;
			}

			size_t ri = 1;
			t(ri, 0) << "normal";
			for(size_t i = 0; i < di.size(); ++i)
				t(ri, i+1) << ((di[i].interfaceState & IS_NORMAL) != 0);

			ri = 2;
			t(ri, 0) << "vmaster";
			for(size_t i = 0; i < di.size(); ++i)
				t(ri, i+1) << ((di[i].interfaceState & IS_VMASTER) != 0);

			ri = 3;
			t(ri, 0) << "vslave";
			for(size_t i = 0; i < di.size(); ++i)
				t(ri, i+1) << ((di[i].interfaceState & IS_VSLAVE) != 0);

			ri = 4;
			t(ri, 0) << "dummy";
			for(size_t i = 0; i < di.size(); ++i)
				t(ri, i+1) << ((di[i].interfaceState & IS_DUMMY) != 0);

			ri = 5;
			t(ri, 0) << "has parent";
			for(size_t i = 0; i < di.size(); ++i)
				t(ri, i+1) << ((di[i].interfaceState & HAS_PARENT) != 0);

			return t;
		}

	private:
	//	copy construction unsupported.
		DistInfoSupplier(const DistInfoSupplier& di) : m_grid(di.m_grid) {}

		Grid& 		m_grid;
		ADistInfo	m_aDistInfo;
		Grid::AttachmentAccessor<VertexBase, ADistInfo> m_aaDistInfoVRT;
		Grid::AttachmentAccessor<EdgeBase, ADistInfo> m_aaDistInfoEDGE;
		Grid::AttachmentAccessor<Face, ADistInfo> m_aaDistInfoFACE;
		Grid::AttachmentAccessor<Volume, ADistInfo> m_aaDistInfoVOL;
};


////////////////////////////////////////////////////////////////////////////////
///	Communicates the distribution infos through existing interfaces
/**	Distribution infos are packed into the send buffer for each node and are
 * either merged with existing entries or existing entries are simply overwritten.
 * The merge/overwrite behavior can be chosen through the member method enable_merge.
 */
template <class TLayout>
class ComPol_SynchronizeDistInfos : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout								Layout;
		typedef typename Layout::Type				GeomObj;
		typedef typename Layout::Element			Element;
		typedef typename Layout::Interface			Interface;
		typedef typename Interface::const_iterator	InterfaceIter;

		ComPol_SynchronizeDistInfos(DistInfoSupplier& distInfos, bool merge) :
			m_distInfos(distInfos), m_mergeEnabled(merge)	{}

		virtual ~ComPol_SynchronizeDistInfos()	{}

		void enable_merge(bool enable)	{m_mergeEnabled = enable;}
		bool merge_enabled()			{return m_mergeEnabled;}

		virtual int
		get_required_buffer_size(const Interface& interface)		{return -1;}

	///	write target processes and move-flag
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& intfc)
		{
			for(InterfaceIter iter = intfc.begin(); iter != intfc.end(); ++iter){
				Element elem = intfc.get_element(iter);
				Serialize(buff, m_distInfos.get(elem));
			}
			return true;
		}

	///	read target processes and move-flag
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& intfc)
		{
			if(m_mergeEnabled){
				vector<TargetProcInfo> tpInfo;
				for(InterfaceIter iter = intfc.begin(); iter != intfc.end(); ++iter){
					tpInfo.clear();
					Deserialize(buff, tpInfo);

					Element elem = intfc.get_element(iter);
					vector<TargetProcInfo>& tpInfoDest = m_distInfos.get(elem);
					size_t initialInfoSize = tpInfoDest.size();

					for(size_t i_src = 0; i_src < tpInfo.size(); ++i_src){
						int procID = tpInfo[i_src].procID;
						bool gotOne = false;

					//	we only have to check entries up to initialInfoSize, since
					//	all following entries have been added during this operation.
					//	Since there are no double entries in tpInfo, there's no
					//	need to check against those new entries in tpInfoDest.
						for(size_t i = 0; i < initialInfoSize; ++i){
							if(procID == tpInfoDest[i].procID){
								tpInfoDest[i].interfaceState
												|= tpInfo[i_src].interfaceState;

								gotOne = true;
								break;
							}
						}

						if(!gotOne)
							tpInfoDest.push_back(tpInfo[i_src]);
					}
				}
			}
			else{
				for(InterfaceIter iter = intfc.begin(); iter != intfc.end(); ++iter){
					Element elem = intfc.get_element(iter);
					vector<TargetProcInfo>& tpInfo = m_distInfos.get(elem);
					tpInfo.clear();
					Deserialize(buff, tpInfo);
				}
			}

			return true;
		}

	protected:
		DistInfoSupplier& 	m_distInfos;
		bool				m_mergeEnabled;
};


////////////////////////////////////////////////////////////////////////////////
template <class TElem>
static void SynchronizeDistInfos(MultiGrid& mg, DistInfoSupplier& distInfos)
{
	typedef typename GridLayoutMap::Types<TElem>::Layout	ElemLayout;
	GridLayoutMap& glm = mg.distributed_grid_manager()->grid_layout_map();
	pcl::InterfaceCommunicator<ElemLayout>	com;
	ComPol_SynchronizeDistInfos<ElemLayout>	compolSync(distInfos, true);

	compolSync.enable_merge(true);
	com.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, compolSync);
	com.communicate();

	compolSync.enable_merge(false);
	com.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, compolSync);
	com.communicate();

//	if an element is a mulit-v-master, all v-master copies are connected to
//	all v-slave copies. However, multi-v-masters are not connected with one another
//	through h-interfaces. That's a pitty and requires this triple communication
	compolSync.enable_merge(true);
	com.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, compolSync);
	com.communicate();

	compolSync.enable_merge(true);
	com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compolSync);
	com.communicate();

	compolSync.enable_merge(false);
	com.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, compolSync);
	com.communicate();
}

#ifdef LG_DISTRIBUTION_DEBUG
////////////////////////////////////////////////////////////////////////////////
static void SaveDistSelectorToFile(MGSelector& msel, const char* filename)
{
//	create a subset handler which holds different subsets for the different selection states
	MultiGrid& mg = *msel.multi_grid();
	SubsetHandler sh(mg);

	for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
		for(MGSelector::traits<Volume>::level_iterator iter = msel.begin<Volume>(lvl);
			iter != msel.end<Volume>(lvl); ++iter)
		{
			sh.assign_subset(*iter, msel.get_selection_status(*iter));
		}

		for(MGSelector::traits<Face>::level_iterator iter = msel.begin<Face>(lvl);
			iter != msel.end<Face>(lvl); ++iter)
		{
			sh.assign_subset(*iter, msel.get_selection_status(*iter));
		}

		for(MGSelector::traits<EdgeBase>::level_iterator iter = msel.begin<EdgeBase>(lvl);
			iter != msel.end<EdgeBase>(lvl); ++iter)
		{
			sh.assign_subset(*iter, msel.get_selection_status(*iter));
		}

		for(MGSelector::traits<VertexBase>::level_iterator iter = msel.begin<VertexBase>(lvl);
			iter != msel.end<VertexBase>(lvl); ++iter)
		{
			sh.assign_subset(*iter, msel.get_selection_status(*iter));
		}
	}

	const char* subsetNames[] = {"unassigned", "normal", "vmaster", "normal+vmaster",
								 "vslave", "normal+vslave", "vmaster+vslave",
								 "normal+vmaster+vslave", "dummy", "normal+dummy",
								 "vmaster+dummy", "normal+vmaster+dummy", "vslave+dummy",
								 "normal+vslave+dummy", "vmaster+vslave+dummy",
								 "normal+vmaster+vslave+dummy"};

	for(int i = 0; i < 16; ++i)
		sh.subset_info(i).name = subsetNames[i];

	AssignSubsetColors(sh);
	EraseEmptySubsets(sh);
	SaveGridHierarchyTransformed(mg, sh, filename, LG_DISTRIBUTION_Z_OUTPUT_TRANSFORM);
}

////////////////////////////////////////////////////////////////////////////////
static void SaveDistInfosToFile(MultiGrid& mg, DistInfoSupplier& infoSupplier,
								const char* filename)
{
//	create a subset handler which holds different subsets for the different selection states
	SubsetHandler sh(mg);

//	write a file for each adressed process
	for(int pi = 0; pi < pcl::NumProcs(); ++pi){
		sh.clear();

		for(MultiGrid::traits<Volume>::iterator iter = mg.begin<Volume>();
			iter != mg.end<Volume>(); ++iter)
		{
			vector<TargetProcInfo>& infos = infoSupplier.get(*iter);
			for(size_t i = 0; i < infos.size(); ++i){
				if(infos[i].procID == pi)
					sh.assign_subset(*iter, infos[i].interfaceState);
			}
		}

		for(MultiGrid::traits<Face>::iterator iter = mg.begin<Face>();
			iter != mg.end<Face>(); ++iter)
		{
			vector<TargetProcInfo>& infos = infoSupplier.get(*iter);
			for(size_t i = 0; i < infos.size(); ++i){
				if(infos[i].procID == pi)
					sh.assign_subset(*iter, infos[i].interfaceState);
			}
		}

		for(MultiGrid::traits<EdgeBase>::iterator iter = mg.begin<EdgeBase>();
			iter != mg.end<EdgeBase>(); ++iter)
		{
			vector<TargetProcInfo>& infos = infoSupplier.get(*iter);
			for(size_t i = 0; i < infos.size(); ++i){
				if(infos[i].procID == pi)
					sh.assign_subset(*iter, infos[i].interfaceState);
			}
		}

		for(MultiGrid::traits<VertexBase>::iterator iter = mg.begin<VertexBase>();
			iter != mg.end<VertexBase>(); ++iter)
		{
			vector<TargetProcInfo>& infos = infoSupplier.get(*iter);
			for(size_t i = 0; i < infos.size(); ++i){
				if(infos[i].procID == pi)
					sh.assign_subset(*iter, infos[i].interfaceState);
			}
		}

		const char* subsetNames[] = {"unassigned", "normal", "vmaster", "normal+vmaster",
									 "vslave", "normal+vslave", "vmaster+vslave",
									 "normal+vmaster+vslave", "dummy", "normal+dummy",
									 "vmaster+dummy", "normal+vmaster+dummy", "vslave+dummy",
									 "normal+vslave+dummy", "vmaster+vslave+dummy",
									 "normal+vmaster+vslave+dummy"};

		for(int i = 0; i < 16; ++i)
			sh.subset_info(i).name = subsetNames[i];

		AssignSubsetColors(sh);
		EraseEmptySubsets(sh);
		if(sh.num_subsets() > 0){
			stringstream ss;
			ss << filename << "_p" << pcl::ProcRank() << "_for_p" << pi << ".ugx";
			SaveGridHierarchyTransformed(mg, sh, ss.str().c_str(), LG_DISTRIBUTION_Z_OUTPUT_TRANSFORM);
		}
	}
}

template <class TElem>
static void WriteDistInfosToTextFile(MultiGrid& mg, DistInfoSupplier& infoSupplier,
									 const char* filename)
{
	typedef typename MultiGrid::traits<TElem>::iterator TElemIter;

	Table<std::stringstream> table(mg.num<TElem>() + 1, 3);
	table(0, 0) << "lvl";	table(0, 1) << "center";	table(0, 2) << "interface states";

	int row = 1;
	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		for(TElemIter iter = mg.begin<TElem>(lvl); iter != mg.end<TElem>(lvl);
			++iter, ++row)
		{
			TElem* e = *iter;

			table(row, 0) << lvl;
			table(row, 1) << GetGeometricObjectCenter(mg, e);

			vector<TargetProcInfo>& infos = infoSupplier.get(e);

			for(size_t i = 0; i < infos.size(); ++i){
				table(row, 2) << "p" << infos[i].procID << ": ";
				byte is = infos[i].interfaceState;
				if(is & IS_NORMAL)	table(row, 2) << "normal ";
				if(is & IS_VMASTER)	table(row, 2) << "vmaster ";
				if(is & IS_VSLAVE)	table(row, 2) << "vslave ";
				if(is & IS_DUMMY)	table(row, 2) << "dummy ";
			}

			table(row, 2) << "| ";
		}
	}

	ofstream out(filename);
	if(!out){
		UG_THROW("Couldn't open file " << filename << " for output.");
	}
	out << table;
	out.close();
}

template <class TElem>
static string LocateElement(MultiGrid& mg, TElem* e)
{
	stringstream ssLocator;
	ssLocator << "at " << GetGeometricObjectCenter(mg, e)
			  << " on level " << mg.get_level(e);
	return ssLocator.str();
}

template <class TElem>
static bool PerformValidityCheck(DistributedGridManager& dgm)
{
	typedef typename Grid::traits<TElem>::iterator	TElemIter;

	bool isValid = true;

	UG_LOG("DEBUG: Performing validity check on distributed grid ");
	switch(TElem::BASE_OBJECT_ID){
	case VERTEX:	UG_LOG("for vertices:\n"); break;
	case EDGE:		UG_LOG("for edges:\n"); break;
	case FACE:		UG_LOG("for faces:\n"); break;
	case VOLUME:	UG_LOG("for volumes:\n"); break;
	}

	MultiGrid& mg = *dgm.get_assigned_grid();
	for(TElemIter iter = mg.begin<TElem>(); iter != mg.end<TElem>(); ++iter)
	{
		TElem* e = *iter;
	//	Make sure that pure vertical masters (ghosts) do not have children
		if(dgm.is_ghost(e)){
			if(mg.has_children(e)){
				UG_LOG("  Ghost has child " << LocateElement(mg, e) << endl);
				isValid = false;
			}
		}
	//	Make sure that pure vertical slaves do not have parents
		else if(dgm.contains_status(e, ES_V_SLAVE)
				&& (!dgm.is_in_horizontal_interface(e)))
		{
			if(mg.get_parent(e)){
				UG_LOG("  Pure vertical slave has parent " << LocateElement(mg, e) << endl);
				isValid = false;
			}
		}
	}

	UG_LOG("DEBUG: Validity check done with result ");
	if(isValid){
		UG_LOG("SUCCESS\n");
	}
	else{
		UG_LOG("FAIL\n");
	}

	return isValid;
}

static bool PerformValidityCheck(DistributedGridManager& dgm)
{
	bool isValid = true;
	isValid &= PerformValidityCheck<VertexBase>(dgm);
	isValid &= PerformValidityCheck<EdgeBase>(dgm);
	isValid &= PerformValidityCheck<Face>(dgm);
	isValid &= PerformValidityCheck<Volume>(dgm);
	return isValid;
}

#endif //LG_DISTRIBUTION_DEBUG

////////////////////////////////////////////////////////////////////////////////
// ATTENTION - THIS DOESN'T REALLY WORK!
//	mpirun -n 4 ugshell -ex adaptive_mg/moving_front.lua -redistributionSteps 2 -redistributionProcs 2
//template <class TElem>
//void AdjustGhostSelection(MGSelector& msel, ISelector::status_t status)
//{
//	DistributedGridManager& dgm = *msel.grid()->distributed_grid_manager();
//	GridLayoutMap& glm = dgm.grid_layout_map();
//
//	if(!glm.has_layout<TElem>(INT_V_MASTER))
//		return;
//
//	typedef typename GridLayoutMap::Types<TElem>::Layout	Layout;
//	typedef typename Layout::iterator						LIter;
//	typedef typename Layout::Interface						Interface;
//	typedef typename Interface::iterator					IIter;
//
//	Layout& layout = glm.get_layout<TElem>(INT_V_MASTER);
//	for(size_t lvl = 0; lvl < layout.num_levels(); ++lvl){
//		for(LIter liter = layout.begin(lvl); liter != layout.end(lvl); ++liter){
//			Interface& intfc = layout.interface(liter);
//			for(IIter iiter = intfc.begin(); iiter != intfc.end(); ++iiter){
//				TElem* e = intfc.get_element(iiter);
//				if(dgm.is_ghost(e)){
//					msel.select(e, status);
//				}
//			}
//		}
//	}
//}
//
//void AdjustGhostSelection(MGSelector& msel, ISelector::status_t status)
//{
//	UG_LOG("DEBUG: AdjustGhostSelection, #sel-vrts before: " << msel.num<VertexBase>() << endl);
//	Grid& g = *msel.grid();
//	if(g.num_vertices())
//		AdjustGhostSelection<VertexBase>(msel, status);
//	if(g.num_edges())
//		AdjustGhostSelection<EdgeBase>(msel, status);
//	if(g.num_faces())
//		AdjustGhostSelection<Face>(msel, status);
//	if(g.num_volumes())
//		AdjustGhostSelection<Volume>(msel, status);
//	UG_LOG("DEBUG: AdjustGhostSelection, #sel-vrts after: " << msel.num<VertexBase>() << endl);
//}


////////////////////////////////////////////////////////////////////////////////
///	Recursively selects unselected sides.
template <class TElem>
static void SelectAssociatedSides(MGSelector& msel, TElem* e,
								  ISelector::status_t status = ISelector::SELECTED)
{
	//UG_DLOG(LIB_GRID, 1, "dist-start: SelectAssociatedSides\n");
	GDIST_PROFILE_FUNC();

	UG_ASSERT(msel.multi_grid(), "");
	MultiGrid& mg = *msel.multi_grid();

	typedef typename TElem::side TSide;
	typename MultiGrid::traits<TSide>::secure_container sides;

	mg.associated_elements(sides, e);
	for(size_t i = 0; i < sides.size(); ++i){
		TSide* s = sides[i];
		//if(!msel.is_selected(sides[i])){
			ISelector::status_t nstate = status | msel.get_selection_status(s);
			msel.select(s, nstate);
			if(TElem::HAS_SIDES)
				SelectAssociatedSides(msel, s, nstate);
		//}
	}

	//UG_DLOG(LIB_GRID, 1, "dist-stop: SelectAssociatedSides\n");
}


////////////////////////////////////////////////////////////////////////////////
/**	selects unselected constrained elements of all selected constraining elements
 * and associated unselected low-dim elems. An exception is made for constraining
 * elements which are pure vertical masters. Associated constrained elements won't
 * be selected in this case, since pure vertical masters mustn't have children.*/
static void SelectAssociatedConstrainedElements(MGSelector& msel,
								ISelector::status_t status = ISelector::SELECTED)
{
	UG_DLOG(LIB_GRID, 1, "dist-start: SelectAssociatedConstrainedElements\n");
	GDIST_PROFILE_FUNC();

	const bool selectAll = true;

//	constraining triangles
	{
		typedef ConstrainingTriangle TElem;
		typedef MGSelector::traits<TElem>::level_iterator TIter;
		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
			for(TIter iter = msel.begin<TElem>(lvl);
				iter != msel.end<TElem>(lvl); ++iter)
			{
				ConstrainingFace* e = *iter;
			//	we won't select constrained elements of pure v-masters!
				if(msel.get_selection_status(e) == IS_VMASTER)
					continue;

				for(size_t i = 0; i < e->num_constrained_vertices(); ++i){
					VertexBase* cd = e->constrained_vertex(i);
					ISelector::status_t nstate = status | msel.get_selection_status(cd);
					if(selectAll || !msel.is_selected(cd)){
						msel.select(cd, nstate);
					}
				}
				for(size_t i = 0; i < e->num_constrained_edges(); ++i){
					EdgeBase* cd = e->constrained_edge(i);
					ISelector::status_t nstate = status | msel.get_selection_status(cd);
					if(selectAll || !msel.is_selected(cd)){
						msel.select(cd, nstate);
						SelectAssociatedSides(msel, cd, nstate);
					}
				}
				for(size_t i = 0; i < e->num_constrained_faces(); ++i){
					Face* cd = e->constrained_face(i);
					ISelector::status_t nstate = status | msel.get_selection_status(cd);
					if(selectAll || !msel.is_selected(cd)){
						msel.select(cd, nstate);
						SelectAssociatedSides(msel, cd, nstate);
					}
				}
			}
		}
	}

//	constraining quadrilaterals
	{
		typedef ConstrainingQuadrilateral TElem;
		typedef MGSelector::traits<TElem>::level_iterator TIter;
		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
			for(TIter iter = msel.begin<TElem>(lvl);
				iter != msel.end<TElem>(lvl); ++iter)
			{
				ConstrainingFace* e = *iter;
			//	we won't select constrained elements of pure v-masters!
				if(msel.get_selection_status(e) == IS_VMASTER)
					continue;

				for(size_t i = 0; i < e->num_constrained_vertices(); ++i){
					VertexBase* cd = e->constrained_vertex(i);
					ISelector::status_t nstate = status | msel.get_selection_status(cd);
					if(selectAll || !msel.is_selected(cd)){
						msel.select(cd, nstate);
					}
				}
				for(size_t i = 0; i < e->num_constrained_edges(); ++i){
					EdgeBase* cd = e->constrained_edge(i);
					ISelector::status_t nstate = status | msel.get_selection_status(cd);
					if(selectAll || !msel.is_selected(cd)){
						msel.select(cd, nstate);
						SelectAssociatedSides(msel, cd, nstate);
					}
				}
				for(size_t i = 0; i < e->num_constrained_faces(); ++i){
					Face* cd = e->constrained_face(i);
					ISelector::status_t nstate = status | msel.get_selection_status(cd);
					if(selectAll || !msel.is_selected(cd)){
						msel.select(cd, nstate);
						SelectAssociatedSides(msel, cd, nstate);
					}
				}
			}
		}
	}

//	constraining edges
	{
		typedef ConstrainingEdge TElem;
		typedef MGSelector::traits<TElem>::level_iterator TIter;
		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
			for(TIter iter = msel.begin<TElem>(lvl);
				iter != msel.end<TElem>(lvl); ++iter)
			{
				ConstrainingEdge* e = *iter;
			//	we won't select constrained elements of pure v-masters!
				if(msel.get_selection_status(e) == IS_VMASTER)
					continue;

				for(size_t i = 0; i < e->num_constrained_vertices(); ++i){
					VertexBase* cd = e->constrained_vertex(i);
					ISelector::status_t nstate = status | msel.get_selection_status(cd);
					if(selectAll || !msel.is_selected(cd)){
						msel.select(cd, nstate);
					}
				}
				for(size_t i = 0; i < e->num_constrained_edges(); ++i){
					EdgeBase* cd = e->constrained_edge(i);
					ISelector::status_t nstate = status | msel.get_selection_status(cd);
					if(selectAll || !msel.is_selected(cd)){
						msel.select(cd, nstate);
						SelectAssociatedSides(msel, cd, nstate);
					}
				}
			}
		}
	}
	UG_DLOG(LIB_GRID, 1, "dist-stop: SelectAssociatedConstrainedElements\n");
}

////////////////////////////////////////////////////////////////////////////////
/**	selects unselected constraining elements of all selected constrained elements
 * and associated unselected low-dim elems.*/
//static void SelectAssociatedConstrainingElements(MGSelector& msel,
//								ISelector::status_t status = ISelector::SELECTED)
//{
//	UG_DLOG(LIB_GRID, 1, "dist-start: SelectAssociatedConstrainingElements\n");
//	GDIST_PROFILE_FUNC();
//
//	const bool selectAll = true;
//
////	constrained triangles
//	{
//		typedef ConstrainedTriangle TElem;
//		typedef MGSelector::traits<TElem>::level_iterator TIter;
//		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
//			for(TIter iter = msel.begin<TElem>(lvl);
//				iter != msel.end<TElem>(lvl); ++iter)
//			{
//				ConstrainedFace* e = *iter;
//				if(GeometricObject* cg = e->get_constraining_object()){
//				//	we won't select pure v-masters!
////					if(msel.get_selection_status(cg) == IS_VMASTER)
////						continue;
//
//					ISelector::status_t nstate = status | msel.get_selection_status(cg);
//					if(selectAll || !msel.is_selected(cg)){
//						msel.select(cg, nstate);
//						UG_ASSERT(dynamic_cast<ConstrainingFace*>(cg),
//								  "constraining object of a face has to be a "
//								  "ConstrainingFace!");
//						SelectAssociatedSides(msel, static_cast<Face*>(cg), nstate);
//					}
//				}
//			}
//		}
//	}
//
////	constrained quadrilaterals
//	{
//		typedef ConstrainedQuadrilateral TElem;
//		typedef MGSelector::traits<TElem>::level_iterator TIter;
//		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
//			for(TIter iter = msel.begin<TElem>(lvl);
//				iter != msel.end<TElem>(lvl); ++iter)
//			{
//				ConstrainedFace* e = *iter;
//				if(GeometricObject* cg = e->get_constraining_object()){
//				//	we won't select pure v-masters!
////					if(msel.get_selection_status(cg) == IS_VMASTER)
////						continue;
//
//					ISelector::status_t nstate = status | msel.get_selection_status(cg);
//					if(selectAll || !msel.is_selected(cg)){
//						msel.select(cg, nstate);
//						UG_ASSERT(dynamic_cast<Face*>(cg),
//								  "constraining object of a face has to be a "
//								  "Face!");
//						SelectAssociatedSides(msel, static_cast<Face*>(cg), nstate);
//					}
//				}
//			}
//		}
//	}
//
////	constrained edges
//	{
//		typedef ConstrainedEdge TElem;
//		typedef MGSelector::traits<TElem>::level_iterator TIter;
//		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
//			for(TIter iter = msel.begin<TElem>(lvl);
//				iter != msel.end<TElem>(lvl); ++iter)
//			{
//				ConstrainedEdge* e = *iter;
//				if(GeometricObject* cg = e->get_constraining_object()){
//				//	we won't select pure v-masters!
////					if(msel.get_selection_status(cg) == IS_VMASTER)
////						continue;
//
//					ISelector::status_t nstate = status | msel.get_selection_status(cg);
//					if(selectAll || !msel.is_selected(cg)){
//						msel.select(cg, nstate);
//						switch(cg->base_object_id()){
//						case EDGE:
//							SelectAssociatedSides(msel, static_cast<EdgeBase*>(cg), nstate);
//							break;
//						case FACE:
//							SelectAssociatedSides(msel, static_cast<Face*>(cg), nstate);
//							break;
//						}
//					}
//				}
//			}
//		}
//	}
//
////	constrained vertices
//	{
//		typedef ConstrainedVertex TElem;
//		typedef MGSelector::traits<TElem>::level_iterator TIter;
//		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
//			for(TIter iter = msel.begin<TElem>(lvl);
//				iter != msel.end<TElem>(lvl); ++iter)
//			{
//				ConstrainedVertex* e = *iter;
//				if(GeometricObject* cg = e->get_constraining_object()){
//				//	we won't select pure v-masters!
////					if(msel.get_selection_status(cg) == IS_VMASTER)
////						continue;
//
//					ISelector::status_t nstate = status | msel.get_selection_status(cg);
//					if(selectAll || !msel.is_selected(cg)){
//						msel.select(cg, nstate);
//						switch(cg->base_object_id()){
//						case EDGE:
//							SelectAssociatedSides(msel, static_cast<EdgeBase*>(cg), nstate);
//							break;
//						case FACE:
//							SelectAssociatedSides(msel, static_cast<Face*>(cg), nstate);
//							break;
//						}
//					}
//				}
//			}
//		}
//	}
//	UG_DLOG(LIB_GRID, 1, "dist-stop: SelectAssociatedConstrainingElements\n");
//}


////////////////////////////////////////////////////////////////////////////////
///	Recursively selects all children of selected vertices
/**	This method is required, since if a distributed vertex has a child which
 * is not connected to other distributed elements, then the child wouldn't be selected.
 * Note that for edges and faces the methods SelectAssociatedConstrainedElements
 * takes care of this.
 * Children of pure vertical masters won't be selected, since those mustn't have
 * children.
 */
static void SelectChildrenOfSelectedShadowVertices(MGSelector& msel,
								ISelector::status_t status = ISelector::SELECTED)
{
	UG_DLOG(LIB_GRID, 1, "dist-start: SelectChildrenOfSelectedShadowVertices\n");
	GDIST_PROFILE_FUNC();

	UG_ASSERT(msel.multi_grid(), "The selector has to operate on a grid!");
	MultiGrid& mg = *msel.multi_grid();
	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();

	Grid::edge_traits::secure_container edges;

	for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
		for(MGSelector::traits<VertexBase>::level_iterator iter = msel.begin<VertexBase>(lvl);
			iter != msel.end<VertexBase>(lvl); ++iter)
		{
			VertexBase* vrt = *iter;
			if(msel.get_selection_status(vrt) == IS_VMASTER)
				continue;

			VertexBase* child = mg.get_child_vertex(vrt);
			if(!child)
				continue;

		//	check whether vrt has an associated edge which does not have children
			mg.associated_elements(edges, vrt);
			for(size_t i = 0; i < edges.size(); ++i){
				if(!(distGridMgr.is_ghost(edges[i]) || mg.has_children(edges[i]))){
					msel.select(child, msel.get_selection_status(child) | status);
					break;
				}
			}
		}
	}
	UG_DLOG(LIB_GRID, 1, "dist-stop: SelectChildrenOfSelectedShadowVertices\n");
}


////////////////////////////////////////////////////////////////////////////////
/**	The method operates on selected entries only. Make sure that all elements
 * of type TElem which are being sent to a process are selected.
 *
 * If a selected element has no children and if it is a vertical master, it will
 * be marked as vertical master again.
 *
 * If a selected element has unselected children, then those children will be
 * selected as vertical master.
 *
 * This method only works correctly if called for the elements of highest dimension.
 */
template <class TElem>
static void AssignVerticalMasterAndSlaveStates(MGSelector& msel, bool partitionForLocalProc)
{
	UG_DLOG(LIB_GRID, 1, "dist-start: AssignVerticalMasterAndSlaveStates\n");
	GDIST_PROFILE_FUNC();

	UG_ASSERT(msel.multi_grid(), "Selector has to operate on a MultiGrid");
	MultiGrid& mg = *msel.multi_grid();
	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();

//	we start on the highest level and go downwards to avoid side
//	effects from repeated selection adjustment.
	typedef typename MGSelector::traits<TElem>::level_iterator TIter;
	for(int lvl = (int)msel.num_levels() - 1; lvl >= 0 ; --lvl){
		for(TIter iter = msel.begin<TElem>(lvl);
			iter != msel.end<TElem>(lvl);)
		{
			TElem* e = *iter;
			++iter;

			if((msel.get_selection_status(e) & IS_VMASTER)
				|| (msel.get_selection_status(e) & IS_VSLAVE))
			{
			//	nothing to do here...
				continue;
			}

		//	assign vertical master states first
			size_t numChildren = mg.num_children<TElem>(e);
			GeometricObject* parent = mg.get_parent(e);
			bool parentIsSelected = false;
			if(parent)
				parentIsSelected = msel.is_selected(parent);

			if(numChildren){
				for(size_t i = 0; i < numChildren; ++i){
					TElem* c = mg.get_child<TElem>(e, i);
					if(!msel.is_selected(c))
						msel.select(c, IS_VMASTER);
				}
			}
			else if(distGridMgr.contains_status(e, ES_V_MASTER)){
				if(parentIsSelected || (partitionForLocalProc && (lvl == 0))){
					msel.select(e, IS_VMASTER);
					continue;
				}
				else{
					msel.deselect(e);
					continue;
				}
			}

		//	and now slave states
			if(parent){
				if(!msel.is_selected(parent))
					msel.select(e, IS_VSLAVE);
			}
			else{
				if(distGridMgr.contains_status(e, ES_V_SLAVE))
					msel.select(e, IS_VSLAVE);
			}
		}
	}

	UG_DLOG(LIB_GRID, 1, "dist-stop: AssignVerticalMasterAndSlaveStates\n");
}

////////////////////////////////////////////////////////////////////////////////
/**	VSlaves will be ignored.*/
template <class TElem>
static void SelectUnselectedRootElementsAsVMasters(MGSelector& msel)
{
	UG_DLOG(LIB_GRID, 1, "dist-start: SelectUnselectedRootElementsAsVMasters\n");
	GDIST_PROFILE_FUNC();

	typedef typename Grid::traits<TElem>::iterator TIter;

	UG_ASSERT(msel.multi_grid(), "Selector has to operate on a MultiGrid");
	MultiGrid& mg = *msel.multi_grid();
	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();

	for(TIter iter = mg.begin<TElem>(0); iter != mg.end<TElem>(0); ++iter){
		if(!msel.is_selected(*iter)){
			if(!distGridMgr.contains_status(*iter, ES_V_SLAVE))
				msel.select(*iter, IS_VMASTER);
		}
	}
	UG_DLOG(LIB_GRID, 1, "dist-stop: SelectUnselectedRootElementsAsVMasters\n");
}

////////////////////////////////////////////////////////////////////////////////
/**	VMasters will be ignored.*/
template <class TElem>
static void SelectSelectedRootElementsAsVSlaves(MGSelector& msel)
{
	UG_DLOG(LIB_GRID, 1, "dist-start: SelectSelectedRootElementsAsVSlaves\n");
	GDIST_PROFILE_FUNC();

	typedef typename Grid::traits<TElem>::iterator TIter;

	UG_ASSERT(msel.multi_grid(), "Selector has to operate on a MultiGrid");
	MultiGrid& mg = *msel.multi_grid();
	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();

	for(TIter iter = msel.begin<TElem>(0); iter != msel.end<TElem>(0); ++iter){
		if(!distGridMgr.contains_status(*iter, ES_V_MASTER))
			msel.select(*iter, IS_VSLAVE);
	}
	UG_DLOG(LIB_GRID, 1, "dist-stop: SelectSelectedRootElementsAsVSlaves\n");
}

////////////////////////////////////////////////////////////////////////////////
static void SelectElementsForTargetPartition(MGSelector& msel,
								SubsetHandler& shPartition, int partitionIndex,
								bool partitionForLocalProc,
								bool createVerticalInterfaces)
{
	UG_DLOG(LIB_GRID, 1, "dist-start: SelectElementsForTargetPartition\n");
	GDIST_PROFILE_FUNC();

//	elements which do not have parents (so called root-elements), and which are
//	not v-slaves have to have a copy on the local proc.
//	If they are not contained in the partition for the local proc, we'll add a
//	copy on the local proc and make it the v-master copy.
//	Note that this should only affect elements in the base-level.
//	The assignment is performed in SelectUnselectedRootElementsAsVMasters.
//todo	assert that only elements in the base-level do not have parents
//		(regarding the global grid)

	UG_ASSERT(msel.multi_grid(), "Selector has to operate on a MultiGrid");
	MultiGrid& mg = *msel.multi_grid();

	if(mg.num<Volume>() > 0){
		if(partitionIndex >= 0)
			SelectSubsetElements<Volume>(msel, shPartition, partitionIndex, IS_NORMAL);
		if(createVerticalInterfaces){
			if(partitionForLocalProc)
				SelectUnselectedRootElementsAsVMasters<Volume>(msel);
			else
				SelectSelectedRootElementsAsVSlaves<Volume>(msel);
		}
	}
	else if(mg.num<Face>() > 0){
		if(partitionIndex >= 0)
			SelectSubsetElements<Face>(msel, shPartition, partitionIndex, IS_NORMAL);
		if(createVerticalInterfaces){
			if(partitionForLocalProc)
				SelectUnselectedRootElementsAsVMasters<Face>(msel);
			else
				SelectSelectedRootElementsAsVSlaves<Face>(msel);
		}
	}
	else if(mg.num<EdgeBase>() > 0){
		if(partitionIndex >= 0)
			SelectSubsetElements<EdgeBase>(msel, shPartition, partitionIndex, IS_NORMAL);
		if(createVerticalInterfaces){
			if(partitionForLocalProc)
				SelectUnselectedRootElementsAsVMasters<EdgeBase>(msel);
			else
				SelectSelectedRootElementsAsVSlaves<EdgeBase>(msel);
		}
	}
	else if(mg.num<VertexBase>() > 0){
		if(partitionIndex >= 0)
			SelectSubsetElements<VertexBase>(msel, shPartition, partitionIndex, IS_NORMAL);
		if(createVerticalInterfaces){
			if(partitionForLocalProc)
				SelectUnselectedRootElementsAsVMasters<VertexBase>(msel);
			else
				SelectSelectedRootElementsAsVSlaves<VertexBase>(msel);
		}
	}

	if(mg.num<Volume>() > 0){
		if(createVerticalInterfaces)
			AssignVerticalMasterAndSlaveStates<Volume>(msel, partitionForLocalProc);
		AssignSelectionStateToSides<Volume>(msel, true);
	}
	else if(mg.num<Face>() > 0){
		if(createVerticalInterfaces)
			AssignVerticalMasterAndSlaveStates<Face>(msel, partitionForLocalProc);
		AssignSelectionStateToSides<Face>(msel, true);
	}
	else if(mg.num<EdgeBase>() > 0){
		if(createVerticalInterfaces)
			AssignVerticalMasterAndSlaveStates<EdgeBase>(msel, partitionForLocalProc);
		AssignSelectionStateToSides<EdgeBase>(msel, true);
	}
	else if(mg.num<VertexBase>() > 0){
		if(createVerticalInterfaces)
			AssignVerticalMasterAndSlaveStates<VertexBase>(msel, partitionForLocalProc);
	//	no sides to assign...
	}

//	select associated constraining elements first, since they may reference
//	additional unselected constrained elements.
//	UG_LOG("DEBUG: SELECTING CONSTRAINING ELEMENTS...\n");
//	SelectAssociatedConstrainingElements(msel, IS_DUMMY);
	SelectAssociatedConstrainedElements(msel, IS_DUMMY | HAS_PARENT);
	SelectChildrenOfSelectedShadowVertices(msel, IS_DUMMY | HAS_PARENT);
	UG_DLOG(LIB_GRID, 1, "dist-stop: SelectElementsForTargetPartition\n");
}

////////////////////////////////////////////////////////////////////////////////
template <class TElem>
static void AddTargetProcToDistInfos(MGSelector& msel,
									DistInfoSupplier& distInfos, int targetProc)
{
	UG_DLOG(LIB_GRID, 1, "dist-start: AddTargetProcToDistInfos\n");
	GDIST_PROFILE_FUNC();

	typedef typename Grid::traits<TElem>::iterator	TElemIter;


	for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
		for(TElemIter iter = msel.begin<TElem>(lvl);
			iter != msel.end<TElem>(lvl); ++iter)
		{
			TElem* e = *iter;
			byte selState = msel.get_selection_status(e);

			distInfos.get(e).push_back(
					TargetProcInfo(targetProc, selState));
		}
	}

	UG_DLOG(LIB_GRID, 1, "dist-stop: AddTargetProcToDistInfos\n");
}


////////////////////////////////////////////////////////////////////////////////
///	DistInfos are post-processed and some values are adjusted (primarily missing vslaves-marks are added)
/**	In some situations a copy of an element may be marked as vmaster but some
 * associated copies are neither marked as vmaster or vslave. This would be invalid
 * and we have to mark those copies as vslaves in those situations.
 *
 * This occurs in situations where a low-dim element with copies
 * on p1 and p2 (no v-interface) is distributed from p1 to a third process.
 *
 * This method should be called after the distribution infos have been synchronized.
 * Since it performs the exactly same actions on all processes for synchronized
 * dist-infos, no further communication is required afterwards.
 */
template <class TElem>
static void PostProcessDistInfos(MultiGrid& mg, DistInfoSupplier& distInfos)
{
//	iterate over all elements and check for each whether a copy is marked as vmaster.
//	If this is the case, all other elements have to be in v-interfaces, too.
//	If a copy isn't in a v-interface, it will be marked as vslave.
	for(typename MultiGrid::traits<TElem>::iterator iter = mg.begin<TElem>();
		iter != mg.end<TElem>(); ++iter)
	{
		TElem* e = *iter;
		vector<TargetProcInfo>& di = distInfos.get(e);
		if(di.size() < 2)
			continue;

		bool gotVMaster = false;
		bool gotNeither = false;
		for(size_t i = 0; i < di.size(); ++i){
			TargetProcInfo& tpi = di[i];
			if(tpi.interfaceState & IS_VMASTER){
				gotVMaster = true;
			}
			else if(!(tpi.interfaceState & IS_VSLAVE)){
				gotNeither = true;
			}
		}

		if(gotVMaster && gotNeither){
		//	those which have neither a vmaster or vslave mark have to be marked as vslaves.
			for(size_t i = 0; i < di.size(); ++i){
				TargetProcInfo& tpi = di[i];
				if(!(tpi.interfaceState & (IS_VMASTER | IS_VSLAVE))){
					tpi.interfaceState |= IS_VSLAVE;
				}
			}
		}
	}

	#ifdef LG_DISTRIBUTION_DEBUG
		UG_LOG("DEBUG: DUMMY CHECK\n");
		for(typename MultiGrid::traits<TElem>::iterator iter = mg.begin<TElem>();
			iter != mg.end<TElem>(); ++iter)
		{
			TElem* e = *iter;
			vector<TargetProcInfo>& di = distInfos.get(e);

			if(di.size() < 2)
				continue;

			bool allDummies = true;
			for(size_t i = 0; i < di.size(); ++i){
				TargetProcInfo& tpi = di[i];
				bool isNormal = ((tpi.interfaceState & IS_NORMAL) != 0);
				bool isVMaster = ((tpi.interfaceState & IS_VMASTER) != 0);
				bool isVSlave = ((tpi.interfaceState & IS_VSLAVE) != 0);
				bool isDummy = ((tpi.interfaceState & IS_DUMMY) != 0);
				if(isNormal || isVMaster || isVSlave){
					allDummies = false;
					break;
				}
				else{
					if(!isDummy){
						UG_THROW("Element doesn't have a valid interface state: "
								 << ElementDebugInfo(mg, e));
					}
				}
			}
			if(allDummies){
				UG_THROW("The element (" << ElementDebugInfo(mg, e) << ") has only dummy marks:\n"
						<< distInfos.get_debug_info(e));
			}
		}
	#endif
}


////////////////////////////////////////////////////////////////////////////////
/**	\param partitionIsEmpty	If no elements are selected for a target-partition, the
 * 							corresponding entry is set to false. This can happen even
 * 							if shPartition contains elements for that partition:
 * 							vmaster elements whose parents are not sent to the same
 * 							partition don't have to be sent either and are thus ignored...
 */
static void FillDistInfos(MultiGrid& mg, SubsetHandler& shPartition, MGSelector& msel,
						DistInfoSupplier& distInfos, const std::vector<int>* processMap,
						const pcl::ProcessCommunicator& procComm,
						bool createVerticalInterfaces,
						vector<bool>& partitionIsEmpty)
{
	UG_DLOG(LIB_GRID, 1, "dist-start: FillDistInfos\n");
	GDIST_PROFILE_FUNC();

	partitionIsEmpty.resize(shPartition.num_subsets());

	for(int i_part = 0; i_part < shPartition.num_subsets(); ++i_part){

		int targetProc = i_part;
		if(processMap)
			targetProc = (*processMap)[i_part];

		bool localPartition = (targetProc == pcl::ProcRank());

		msel.clear();
		SelectElementsForTargetPartition(msel, shPartition, i_part,
									 localPartition, createVerticalInterfaces);

		partitionIsEmpty[i_part] = msel.empty();

		if(!partitionIsEmpty[i_part]){
		//DEBUG:	temporarily save selection to a file
			#ifdef LG_DISTRIBUTION_DEBUG
			{
				stringstream ss;
				ss << "dist-selection-p" << pcl::ProcRank() << "for-p"<< i_part << ".ugx";
				SaveDistSelectorToFile(msel, ss.str().c_str());
			}
			#endif

			AddTargetProcToDistInfos<Volume>(msel, distInfos, targetProc);
			AddTargetProcToDistInfos<Face>(msel, distInfos, targetProc);
			AddTargetProcToDistInfos<EdgeBase>(msel, distInfos, targetProc);
			AddTargetProcToDistInfos<VertexBase>(msel, distInfos, targetProc);
		}
	}

#ifdef LG_DISTRIBUTION_DEBUG
	{
		//stringstream ss;
		//ss << "dist_infos_vrt_before_sync_p" << pcl::ProcRank() << ".ugx";
		//WriteDistInfosToTextFile<VertexBase>(mg, distInfos, ss.str().c_str());
		SaveDistInfosToFile(mg, distInfos, "dist_infos_before_sync");
	}
#endif

	SynchronizeDistInfos<VertexBase>(mg, distInfos);
	SynchronizeDistInfos<EdgeBase>(mg, distInfos);
	SynchronizeDistInfos<Face>(mg, distInfos);
	SynchronizeDistInfos<Volume>(mg, distInfos);

	PostProcessDistInfos<VertexBase>(mg, distInfos);
	PostProcessDistInfos<EdgeBase>(mg, distInfos);
	PostProcessDistInfos<Face>(mg, distInfos);
	PostProcessDistInfos<Volume>(mg, distInfos);

#ifdef LG_DISTRIBUTION_DEBUG
	{
		//stringstream ss;
		//ss << "dist_infos_vrt_after_sync_p" << pcl::ProcRank() << ".ugx";
		//WriteDistInfosToTextFile<VertexBase>(mg, distInfos, ss.str().c_str());
		SaveDistInfosToFile(mg, distInfos, "dist_infos_after_sync");
	}
#endif

	UG_DLOG(LIB_GRID, 1, "dist-stop: FillDistInfos\n");
}

////////////////////////////////////////////////////////////////////////////////
/**	Based on the list of target processes given by DistInfoSupplier, layouts and
 * interfaces are generated in the given GridLayoutMap.
 *
 * \todo	Think about caching interfaces to speed up this method.
 */
template <class TElem>
static void CreateLayoutsFromDistInfos(MultiGrid& mg, GridLayoutMap& glm,
										DistInfoSupplier& distInfos,
										AGeomObjID& aGID)
{
	UG_DLOG(LIB_GRID, 1, "dist-start: CreateLayoutsFromDistInfos\n");
	GDIST_PROFILE_FUNC();

	typedef typename MultiGrid::traits<TElem>::iterator	TIter;


	int localProcID = pcl::ProcRank();

	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		for(TIter iter = mg.begin<TElem>(lvl); iter != mg.end<TElem>(lvl); ++iter)
		{
			TElem* e = *iter;
			vector<TargetProcInfo>& di = distInfos.get(e);

			if(di.size() < 2)
				continue;

		//	get the process with the lowest rank, on which a normal copy of this
		//	element lies (ignore pure vertical masters)
		//	this lowest rank is required to decide, which process a horizontal
		//	master should reside on
			byte localInterfaceState = 0;
			int minProc = pcl::NumProcs();
			int minVMasterProc = pcl::NumProcs();
			int minVMasterNoVSlave = pcl::NumProcs();
			int minNormalProc = pcl::NumProcs();

		//	the lowest proc which holds a v-slave or a normal entry.
		//	dummies are ignored here, since we don't want them to be h-masters.
			int minRegularHMasterProc = pcl::NumProcs();
			bool isVMaster = false;
			bool isVSlave = false;
			bool isNormal = false;
			bool isDummy = false;

			bool vMasterExists = false;
			bool dummyExists = false;
			//bool isDummy = false;
			//bool isNormal = false;
			bool createNormalHInterface = false;
//			bool forceHInterface = false;

			int numVSlaveProcs = 0;
			for(size_t i = 0; i < di.size(); ++i){
				TargetProcInfo& tpi = di[i];
				if(tpi.procID == localProcID)
					localInterfaceState = tpi.interfaceState;

				if(tpi.interfaceState & IS_VMASTER){
				//	if there is more than one vmaster, then we have to build
				//	h interfaces
//					if(vMasterExists){
//						createNormalHInterface = true;
//						forceHInterface = true;
//					}

					vMasterExists = true;
					if(tpi.procID == localProcID){
						isVMaster = true;
					}
					if(tpi.procID < minVMasterProc)
						minVMasterProc = tpi.procID;
					if(!(tpi.interfaceState & IS_VSLAVE)){
						if(tpi.procID < minVMasterNoVSlave)
							minVMasterNoVSlave = tpi.procID;
					}
				}

				if(tpi.interfaceState & IS_VSLAVE){
					if(tpi.procID == localProcID)
						isVSlave = true;
					if(tpi.procID < minRegularHMasterProc)
						minRegularHMasterProc = tpi.procID;
					++numVSlaveProcs;
				}

				if(tpi.interfaceState & (IS_NORMAL)){
					createNormalHInterface = true;
					if(tpi.procID < minRegularHMasterProc){
						minRegularHMasterProc = tpi.procID;
						minNormalProc = tpi.procID;
					}
					if(tpi.procID == localProcID)
						isNormal = true;
				}

				if(tpi.interfaceState & (IS_DUMMY)){
					createNormalHInterface = true;
					dummyExists = true;
					if(tpi.procID == localProcID)
						isDummy = true;

				//	if you don't want to have dummies, which are h-masters, then
				//	remove the following lines
//					if(tpi.procID < minRegularHMasterProc)
//						minRegularHMasterProc = tpi.procID;

				}
				if(tpi.procID < minProc)
					minProc = tpi.procID;
			}

			UG_ASSERT((!createNormalHInterface)
					  || (minRegularHMasterProc < pcl::NumProcs()),
					  "invalid h-master process. The local node (" << ElementDebugInfo(mg, e)
					  << ") has the following flags:\n"
					  << distInfos.get_debug_info(e) << "\n");

		//	if one process is marked as vmaster but not as a vslave, then we have
		//	to be careful if we adjust states on processes which are marked as
		//	both vmaster and vslave.
			if(minVMasterNoVSlave < pcl::NumProcs())
				minVMasterProc = minVMasterNoVSlave;

		//	in some situations, lower dimensional elements can be marked as a
		//	vmaster and as vslave at the same time. We currently allow for elements
		//	to be part of both interfaces at the same time.
			if(isVMaster && isVSlave){
//				if(localProcID == minVMasterProc)
//					isVSlave = false;
//				else
//					isVMaster = false;

//				isVMaster = isVSlave = false;

			//	adjacent normal full-dimensional elements should thus exist and a
			//	horizontal interface has to be built.
				createNormalHInterface = true;
				UG_ASSERT(minRegularHMasterProc < pcl::NumProcs(), "invalid h-master process");
			}
//			else if((!(isVMaster || isVSlave)) && vMasterExists){
//			//	check whether a vmaster copy exists. If this is the case,
//			//	the element itself has to be a vslave.
//			//	This occurs in situations where a low-dim element with copies
//			//	on p1 and p2 (no v-interface) is distributed from p1 to a third process.
//			//	a v-interface on p1 will then be created and
//				isVSlave = true;
//			}

		//	dummies are only required where no normal or slave state is set
			if(isDummy && (isNormal || isVSlave))
				isDummy = false;

		//	there only may be one v-master copy
//			if(isVMaster && (localProcID != minVMasterProc)){
//				isVMaster = false;
//				isVSlave = true;
//			}

		//	if this condition is fulfilled, some kind of h-interface will be built
			//bool createHInterface = createNormalHInterface || (numVSlaveProcs > 1);

			for(size_t i = 0; i < di.size(); ++i){
				TargetProcInfo& tpi = di[i];
				if(tpi.procID == localProcID)
					continue;

				bool tpIsVMaster = (tpi.interfaceState & IS_VMASTER);
				bool tpIsVSlave = (tpi.interfaceState & IS_VSLAVE);
				bool tpIsNormal = (tpi.interfaceState & IS_NORMAL);
				bool tpIsDummy = (tpi.interfaceState & IS_DUMMY);


				if(tpIsVMaster && tpIsVSlave){
//					if(tpi.procID == minVMasterProc)
//						tpIsVSlave = false;
//					else
//						tpIsVMaster = false;

//					tpIsVMaster = tpIsVSlave = false;

					createNormalHInterface = true;
					UG_ASSERT(minRegularHMasterProc < pcl::NumProcs(), "invalid h-master process");
				}
//				else if((!(tpIsVMaster || tpIsVSlave)) && vMasterExists){
//					tpIsVSlave = true;
//				}

				if(tpIsDummy && (tpIsNormal || tpIsVSlave))
					tpIsDummy = false;
//				if(tpIsVMaster && (tpi.procID != minVMasterProc)){
//					tpIsVMaster = false;
//					tpIsVSlave = true;
//				}

				bool interfaceCreated = false;

			//	add entry to vertical interface if necessary
				//if(isVSlave && (tpi.procID == minVMasterProc)){
				if(isVSlave && tpIsVMaster){
					glm.get_layout<TElem>(INT_V_SLAVE).
						interface(tpi.procID, lvl).push_back(e);
				}
				if(isVMaster && tpIsVSlave){
					glm.get_layout<TElem>(INT_V_MASTER).
						interface(tpi.procID, lvl).push_back(e);

				//	we have to destroy parent-child relations to possibly existing
				//	children. Those children are now indirectly connected through
				//	copies on other processes.
					//if(!createHInterface && mg.has_children(e))
					if(!createNormalHInterface && mg.has_children(e))
						mg.clear_child_connections(e);
				}
				if(isVSlave && tpIsVSlave){
					UG_ASSERT(minRegularHMasterProc < pcl::NumProcs(), "invalid h-master process");
				//	we still have to build a horizontal interface, this time
				//	however only between vertical slaves
//					if(tpIsVSlave && (!tpWasVMaster)){
//						if(!(isVMaster || wasVMaster)){
//					if(isVSlave && tpIsVSlave){
//					if(tpIsVSlave){
//						if(!(isVMaster)){

					if(localProcID == minRegularHMasterProc){
					//	horizontal master
						interfaceCreated = true;
						glm.get_layout<TElem>(INT_H_MASTER).
							interface(tpi.procID, lvl).push_back(e);
					}
					else if(tpi.procID == minRegularHMasterProc){
					//	horizontal slave
						interfaceCreated = true;
						glm.get_layout<TElem>(INT_H_SLAVE).
							interface(tpi.procID, lvl).push_back(e);
					}

//						}
//					}
				}

			//	add entry to horizontal interface if necessary
				if(!interfaceCreated && createNormalHInterface){
					UG_ASSERT(minRegularHMasterProc < pcl::NumProcs(), "invalid h-master process");

				//	check whether the target process would also create a normal h interface
					if(localProcID == minRegularHMasterProc){
					//	horizontal master
					//	only build the interface if the process is not a pure
					//	v-master
						if(tpi.interfaceState != IS_VMASTER){
						//if((tpi.procID != minVMasterProc) || (tpi.interfaceState != IS_VMASTER)){
						//if(forceHInterface || (tpi.interfaceState != IS_VMASTER)){
							glm.get_layout<TElem>(INT_H_MASTER).
								interface(tpi.procID, lvl).push_back(e);
						}
					}
					else if(tpi.procID == minRegularHMasterProc){
					//	horizontal slave
					//	only build the interface if the process is not a pure
					//	v-master
						if(localInterfaceState != IS_VMASTER){
						//if((localProcID != minVMasterProc) || (localInterfaceState != IS_VMASTER)){
						//if(forceHInterface || (localInterfaceState != IS_VMASTER)){
							glm.get_layout<TElem>(INT_H_SLAVE).
								interface(tpi.procID, lvl).push_back(e);
						}
					}
				}

			//	finally we have to make sure, that dummies which do not have a parent
			//	are v-slaves.
				if(dummyExists && (!vMasterExists)){
				//	the lowest normal process will be transformed to a v-master
					UG_ASSERT(minNormalProc < pcl::NumProcs(), "invalid minNormalProc!");

					if(isDummy && (!(localInterfaceState & HAS_PARENT))
						&& (tpi.procID == minNormalProc))
					{
						glm.get_layout<TElem>(INT_V_SLAVE).
								interface(tpi.procID, lvl).push_back(e);
					}
					else if(tpIsDummy && (!(tpi.interfaceState & HAS_PARENT))
							&& (localProcID == minNormalProc))
					{
						glm.get_layout<TElem>(INT_V_MASTER).
								interface(tpi.procID, lvl).push_back(e);
					}
				}
			}
		}
	}

//	Now sort the interface entries in the different layouts
	CompareByAttachment<TElem, AGeomObjID> gidCmp(mg, aGID);
	if(glm.has_layout<TElem>(INT_H_MASTER))
		glm.get_layout<TElem>(INT_H_MASTER).sort_interface_entries(gidCmp);
	if(glm.has_layout<TElem>(INT_H_SLAVE))
		glm.get_layout<TElem>(INT_H_SLAVE).sort_interface_entries(gidCmp);
	if(glm.has_layout<TElem>(INT_V_MASTER))
		glm.get_layout<TElem>(INT_V_MASTER).sort_interface_entries(gidCmp);
	if(glm.has_layout<TElem>(INT_V_SLAVE))
		glm.get_layout<TElem>(INT_V_SLAVE).sort_interface_entries(gidCmp);

	UG_DLOG(LIB_GRID, 1, "dist-stop: CreateLayoutsFromDistInfos\n");
}

////////////////////////////////////////////////////////////////////////////////
bool DistributeGrid(MultiGrid& mg,
					SubsetHandler& shPartition,
					GridDataSerializationHandler& serializer,
					bool createVerticalInterfaces,
					const std::vector<int>* processMap,
					const pcl::ProcessCommunicator& procComm)
{
	GDIST_PROFILE_FUNC();
	PCL_DEBUG_BARRIER(procComm);

	UG_STATIC_ASSERT(IS_DUMMY < 256, RedistributeGrid_IS_DUMMY_too_big);

	UG_DLOG(LIB_GRID, 1, "dist-start: DistributeGrid\n");
	const char* errprefix = "ERROR in DistributeGrid: ";

	if(!mg.is_parallel()){
		UG_THROW(errprefix << "Can't distribute a serial grid! Compile ug with -DPARALLEL=ON");
	}


	UG_DLOG(LIB_GRID, 2, "dist: Informing msg-hub that distribution starts\n");
	GDIST_PROFILE(gdist_distStartsCallback);
	GridDataSerializationHandler	userDataSerializer;
	SPMessageHub msgHub = mg.message_hub();
	msgHub->post_message(GridMessage_Distribution(GMDT_DISTRIBUTION_STARTS, userDataSerializer));
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();


	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();
	GridLayoutMap& glm = distGridMgr.grid_layout_map();

//	The selector will be of frequent use to speed up some algorithms
	MGSelector msel(mg);

	#ifdef LG_DISTRIBUTION_DEBUG
		PerformValidityCheck(distGridMgr);
	#endif

//	Since we will change huge parts of the underlying grid and the grid-layout-map,
//	we'll disable auto-insertion of elements in the distributed-grid-manager.
//	This means we carefully have to take care of all interface changes.
	distGridMgr.enable_interface_management(false);

////////////////////////////////
//	GLOBAL IDS
//todo:	only create global ids if they aren't already present
	GDIST_PROFILE(gdist_CreateGlobalIDs);
	UG_DLOG(LIB_GRID, 2, "dist-DistributeGrid: Create global vertex ids\n");
	CreateAndDistributeGlobalIDs<VertexBase>(mg, glm);
	CreateAndDistributeGlobalIDs<EdgeBase>(mg, glm);
	CreateAndDistributeGlobalIDs<Face>(mg, glm);
	CreateAndDistributeGlobalIDs<Volume>(mg, glm);
	MultiElementAttachmentAccessor<AGeomObjID> aaID(mg, aGeomObjID);
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();


	#ifdef LG_DISTRIBUTION_DEBUG
	{
		UG_LOG("DEBUG: WRITING GLOBAL VERTEX IDS TO FILE\n");
		stringstream ss;
		ss << "global_ids_vrt_p" << pcl::ProcRank() << ".txt";
		WriteDebugValuesToFile<VertexBase>(ss.str().c_str(), mg, aGeomObjID, false);
	}
	{
		UG_LOG("DEBUG: WRITING GLOBAL EDGE IDS TO FILE\n");
		stringstream ss;
		ss << "global_ids_edge_p" << pcl::ProcRank() << ".txt";
		WriteDebugValuesToFile<EdgeBase>(ss.str().c_str(), mg, aGeomObjID, false);
	}
	{
		UG_LOG("DEBUG: WRITING GLOBAL FACE IDS TO FILE\n");
		stringstream ss;
		ss << "global_ids_face_p" << pcl::ProcRank() << ".txt";
		WriteDebugValuesToFile<Face>(ss.str().c_str(), mg, aGeomObjID, false);
	}
	#endif

////////////////////////////////
//	FILL THE DISTRIBUTION INFOS (INVOLVES COMMUNICATION...)
	GDIST_PROFILE(gdist_FillDistInfos);
	UG_DLOG(LIB_GRID, 2, "dist-DistributeGrid: Fill distribution infos\n");
	vector<bool> partitionIsEmpty;
	DistInfoSupplier distInfos(mg);
	FillDistInfos(mg, shPartition, msel, distInfos, processMap, procComm,
				  createVerticalInterfaces, partitionIsEmpty);
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();

//	DEBUG: output distInfos...
	#ifdef LG_DISTRIBUTION_DEBUG
	{
		SaveDistInfosToFile(mg, distInfos, "dist_infos_before_distribution");
	}
	#endif

////////////////////////////////
//	COMMUNICATE INVOLVED PROCESSES
	GDIST_PROFILE(gdist_CommunicateInvolvedProcs);
	UG_DLOG(LIB_GRID, 2, "dist-DistributeGrid: CommunicateInvolvedProcesses\n");

//	each process has to know with which other processes it
//	has to communicate.
	vector<int> sendToRanks, recvFromRanks, sendPartitionInds;

	if(processMap && (shPartition.num_subsets() > (int)processMap->size())){
		UG_THROW("process-map is too small for the given number of partitions!");
	}

//	for each subset which is not emtpy we'll have to send data to
//	the associated process.
	for(int si = 0; si < shPartition.num_subsets(); ++si){
	//	instead of simply querying shPartition.empty(si), we'll check partitionIsEmpty[si],
	//	since this array tells whether data is actually sent to a partition.
	//	E.g. vmasters which are contained in shPartition are not necessarily sent
	//	to the associated target process...
		if(!partitionIsEmpty[si]){
			int toProc = si;
		//	if a process map exists, we'll use the associated process
			if(processMap)
				toProc = processMap->at(si);

			sendToRanks.push_back(toProc);
			sendPartitionInds.push_back(si);
		}
	}

	pcl::CommunicateInvolvedProcesses(recvFromRanks, sendToRanks, procComm);

	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();


////////////////////////////////
//	SERIALIZE THE GRID, THE GLOBAL IDS AND THE DISTRIBUTION INFOS.
	GDIST_PROFILE(gdist_Serialization);
	UG_DLOG(LIB_GRID, 2, "dist-DistributeGrid: Serialization\n");
	AInt aLocalInd("distribution-tmp-local-index");
	mg.attach_to_all(aLocalInd);
	MultiElementAttachmentAccessor<AInt> aaInt(mg, aLocalInd);

//	out and sendSegSizes will be used to distribute the grid.
	BinaryBuffer out;
	vector<int> outSegSizes;

//	the magic number is used for debugging to make sure that the stream is read correctly
	int magicNumber1 = 75234587;
	int magicNumber2 = 560245;

	ADistInfo aDistInfo = distInfos.dist_info_attachment();

	GridDataSerializationHandler distInfoSerializer;
	distInfoSerializer.add(GeomObjAttachmentSerializer<VertexBase, ADistInfo>::create(mg, aDistInfo));
	distInfoSerializer.add(GeomObjAttachmentSerializer<EdgeBase, ADistInfo>::create(mg, aDistInfo));
	distInfoSerializer.add(GeomObjAttachmentSerializer<Face, ADistInfo>::create(mg, aDistInfo));
	distInfoSerializer.add(GeomObjAttachmentSerializer<Volume, ADistInfo>::create(mg, aDistInfo));

//	now perform the serialization
	int localPartitionInd = -1;
	for(size_t i_to = 0; i_to < sendPartitionInds.size(); ++i_to){
		int partInd = sendPartitionInds[i_to];
		bool localPartition = (sendToRanks[i_to] == pcl::ProcRank());
		if(localPartition)
			localPartitionInd = partInd;

	//	the last size is required to calculate the size of the new segment
		size_t oldSize = out.write_pos();

	//	don't serialize the local partition since we'll keep it here on the local
	//	process anyways.
		if(!localPartition){
		//	write a magic number for debugging purposes
			out.write((char*)&magicNumber1, sizeof(int));

		//	select the elements of the current partition
			msel.clear();
			SelectElementsForTargetPartition(msel, shPartition, partInd,
										 localPartition, createVerticalInterfaces);
			//AdjustGhostSelection(msel, ISelector::DESELECTED);

			SerializeMultiGridElements(mg, msel.get_geometric_objects(), aaInt, out, &aaID);


		//	serialize associated data
			distInfoSerializer.write_infos(out);
			distInfoSerializer.serialize(out, msel.get_geometric_objects());
			serializer.write_infos(out);
			serializer.serialize(out, msel.get_geometric_objects());
			userDataSerializer.write_infos(out);
			userDataSerializer.serialize(out, msel.get_geometric_objects());

		//	write a magic number for debugging purposes
			out.write((char*)&magicNumber2, sizeof(int));
		}

	//	size of the segment we just wrote to out
		outSegSizes.push_back((int)(out.write_pos() - oldSize));
	}
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();



////////////////////////////////
//	COMMUNICATE SERIALIZED DATA
	GDIST_PROFILE(gdist_CommunicateSerializedData);
	UG_DLOG(LIB_GRID, 2, "dist-DistributeGrid: Distribute data\n");
//	now distribute the packs between involved processes
	BinaryBuffer in;
	vector<int> inSegSizes(recvFromRanks.size());

	procComm.distribute_data(in, GetDataPtr(inSegSizes),
							GetDataPtr(recvFromRanks), (int)recvFromRanks.size(),
							out.buffer(), GetDataPtr(outSegSizes),
							GetDataPtr(sendToRanks), (int)sendToRanks.size());
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();


////////////////////////////////
//	INTERMEDIATE CLEANUP
	GDIST_PROFILE(gdist_IntermediateCleanup);
	UG_DLOG(LIB_GRID, 2, "dist-DistributeGrid: Intermediate cleanup\n");

	msgHub->post_message(GridMessage_Creation(GMCT_CREATION_STARTS));

//	we have to remove all elements which won't stay on the local process.
//	To do so, we'll first select all elements that stay, invert that selection
//	and erase all elements which are selected thereafter.
	if(createVerticalInterfaces || (localPartitionInd != -1)){
		msel.clear();
		SelectElementsForTargetPartition(msel, shPartition, localPartitionInd,
									 	 true, createVerticalInterfaces);
		InvertSelection(msel);

	//	make sure that constrained/constraining connections won't be harmed
	//	this is a little cumbersome in the moment. Ideally constrained/constraining
	//	elements should unregister from each other automatically on destruction.
		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
			for(ConstrainedVertexIterator iter = msel.begin<ConstrainedVertex>(lvl);
				iter != msel.end<ConstrainedVertex>(lvl); ++iter)
			{
				GeometricObject* co = (*iter)->get_constraining_object();
				if(co && !msel.is_selected(co)){
					switch(co->base_object_id()){
						case EDGE:{
							if(ConstrainingEdge* ce = dynamic_cast<ConstrainingEdge*>(co))
								ce->unconstrain_object(*iter);
						}break;
						case FACE:{
							if(co->reference_object_id() == ROID_TRIANGLE){
								if(ConstrainingTriangle* ce = dynamic_cast<ConstrainingTriangle*>(co))
									ce->unconstrain_object(*iter);
							}
							else{
								if(ConstrainingQuadrilateral* ce = dynamic_cast<ConstrainingQuadrilateral*>(co))
									ce->unconstrain_object(*iter);
							}
						}break;
						default: break;
					}
				}
			}

			for(ConstrainedEdgeIterator iter = msel.begin<ConstrainedEdge>(lvl);
				iter != msel.end<ConstrainedEdge>(lvl); ++iter)
			{
				GeometricObject* co = (*iter)->get_constraining_object();
				if(co && !msel.is_selected(co)){
					switch(co->base_object_id()){
						case EDGE:{
							if(ConstrainingEdge* ce = dynamic_cast<ConstrainingEdge*>(co))
								ce->unconstrain_object(*iter);
						}break;
						case FACE:{
							if(co->reference_object_id() == ROID_TRIANGLE){
								if(ConstrainingTriangle* ce = dynamic_cast<ConstrainingTriangle*>(co))
									ce->unconstrain_object(*iter);
							}
							else{
								if(ConstrainingQuadrilateral* ce = dynamic_cast<ConstrainingQuadrilateral*>(co))
									ce->unconstrain_object(*iter);
							}
						}break;
						default: break;
					}
				}
			}

			for(ConstrainingEdgeIterator iter = msel.begin<ConstrainingEdge>(lvl);
				iter != msel.end<ConstrainingEdge>(lvl); ++iter)
			{
				ConstrainingEdge* e = *iter;
				for(size_t i = 0; i < e->num_constrained_vertices(); ++i){
					ConstrainedVertex* cv = dynamic_cast<ConstrainedVertex*>(e->constrained_vertex(i));
					UG_ASSERT(cv, "Constrained vertices have to be of the type ConstrainedVertex");
					cv->set_constraining_object(NULL);
				}

				for(size_t i = 0; i < e->num_constrained_edges(); ++i){
					ConstrainedEdge* cde = dynamic_cast<ConstrainedEdge*>(e->constrained_edge(i));
					UG_ASSERT(cde, "Constrained edges have to be of the type ConstrainedEdge");
					cde->set_constraining_object(NULL);
				}
			}


			for(ConstrainedTriangleIterator iter = msel.begin<ConstrainedTriangle>(lvl);
				iter != msel.end<ConstrainedTriangle>(lvl); ++iter)
			{
				GeometricObject* co = (*iter)->get_constraining_object();
				if(co && !msel.is_selected(co)){
					if(ConstrainingTriangle* ce = dynamic_cast<ConstrainingTriangle*>(co))
						ce->unconstrain_object(*iter);
				}
			}

			for(ConstrainedQuadrilateralIterator iter = msel.begin<ConstrainedQuadrilateral>(lvl);
				iter != msel.end<ConstrainedQuadrilateral>(lvl); ++iter)
			{
				GeometricObject* co = (*iter)->get_constraining_object();
				if(co && !msel.is_selected(co)){
					if(ConstrainingQuadrilateral* ce = dynamic_cast<ConstrainingQuadrilateral*>(co))
						ce->unconstrain_object(*iter);
				}
			}

			for(ConstrainingTriangleIterator iter = msel.begin<ConstrainingTriangle>(lvl);
				iter != msel.end<ConstrainingTriangle>(lvl); ++iter)
			{
				ConstrainingFace* e = *iter;
				for(size_t i = 0; i < e->num_constrained_vertices(); ++i){
					ConstrainedVertex* cv = dynamic_cast<ConstrainedVertex*>(e->constrained_vertex(i));
					UG_ASSERT(cv, "Constrained vertices have to be of the type ConstrainedVertex");
					cv->set_constraining_object(NULL);
				}

				for(size_t i = 0; i < e->num_constrained_edges(); ++i){
					ConstrainedEdge* cde = dynamic_cast<ConstrainedEdge*>(e->constrained_edge(i));
					UG_ASSERT(cde, "Constrained edges have to be of the type ConstrainedEdge");
					cde->set_constraining_object(NULL);
				}

				for(size_t i = 0; i < e->num_constrained_faces(); ++i){
					ConstrainedFace* cdf = dynamic_cast<ConstrainedFace*>(e->constrained_face(i));
					UG_ASSERT(cdf, "Constrained faces have to be of the type ConstrainedFace");
					cdf->set_constraining_object(NULL);
				}
			}

			for(ConstrainingQuadrilateralIterator iter = msel.begin<ConstrainingQuadrilateral>(lvl);
				iter != msel.end<ConstrainingQuadrilateral>(lvl); ++iter)
			{
				ConstrainingFace* e = *iter;
				for(size_t i = 0; i < e->num_constrained_vertices(); ++i){
					ConstrainedVertex* cv = dynamic_cast<ConstrainedVertex*>(e->constrained_vertex(i));
					UG_ASSERT(cv, "Constrained vertices have to be of the type ConstrainedVertex");
					cv->set_constraining_object(NULL);
				}

				for(size_t i = 0; i < e->num_constrained_edges(); ++i){
					ConstrainedEdge* cde = dynamic_cast<ConstrainedEdge*>(e->constrained_edge(i));
					UG_ASSERT(cde, "Constrained edges have to be of the type ConstrainedEdge");
					cde->set_constraining_object(NULL);
				}

				for(size_t i = 0; i < e->num_constrained_faces(); ++i){
					ConstrainedFace* cdf = dynamic_cast<ConstrainedFace*>(e->constrained_face(i));
					UG_ASSERT(cdf, "Constrained faces have to be of the type ConstrainedFace");
					cdf->set_constraining_object(NULL);
				}
			}
		}

		GDIST_PROFILE(gdist_ErasingObjects);
		EraseSelectedObjects(msel);
		GDIST_PROFILE_END();
	}
	else{
	//	nothing remains on the local process...
		GDIST_PROFILE(gdist_ClearGeometry);
		mg.clear_geometry();
		GDIST_PROFILE_END();
	}

	{
		GDIST_PROFILE(gdist_ClearLayoutMap);
	//	the grid layout map will be rebuilt from scratch
		glm.clear();
		GDIST_PROFILE_END();
	}
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();

////////////////////////////////
//	DESERIALIZE INCOMING GRIDS
	GDIST_PROFILE(gdist_Deserialize);
	distInfoSerializer.deserialization_starts();
	serializer.deserialization_starts();
	userDataSerializer.deserialization_starts();

	vector<VertexBase*>	vrts;
	vector<EdgeBase*> edges;
	vector<Face*> faces;
	vector<Volume*> vols;

	for(size_t i = 0; i < recvFromRanks.size(); ++i){
	//	there is nothing to serialize from the local rank
		if(recvFromRanks[i] == pcl::ProcRank())
			continue;

		UG_DLOG(LIB_GRID, 2, "Deserializing from rank " << recvFromRanks[i] << "\n");

	//	read the magic number and make sure that it matches our magicNumber
		int tmp = 0;
		in.read((char*)&tmp, sizeof(int));
		if(tmp != magicNumber1){
			UG_THROW("ERROR in RedistributeGrid: "
					 "Magic number mismatch before deserialization.\n");
		}

		DeserializeMultiGridElements(mg, in, &vrts, &edges, &faces, &vols, &aaID);

	//	deserialize the associated data (global ids have already been deserialized)
		distInfoSerializer.read_infos(in);
		distInfoSerializer.deserialize(in, vrts.begin(), vrts.end());
		distInfoSerializer.deserialize(in, edges.begin(), edges.end());
		distInfoSerializer.deserialize(in, faces.begin(), faces.end());
		distInfoSerializer.deserialize(in, vols.begin(), vols.end());

		serializer.read_infos(in);
		serializer.deserialize(in, vrts.begin(), vrts.end());
		serializer.deserialize(in, edges.begin(), edges.end());
		serializer.deserialize(in, faces.begin(), faces.end());
		serializer.deserialize(in, vols.begin(), vols.end());

		userDataSerializer.read_infos(in);
		userDataSerializer.deserialize(in, vrts.begin(), vrts.end());
		userDataSerializer.deserialize(in, edges.begin(), edges.end());
		userDataSerializer.deserialize(in, faces.begin(), faces.end());
		userDataSerializer.deserialize(in, vols.begin(), vols.end());

	//	read the magic number and make sure that it matches our magicNumber
		tmp = 0;
		in.read((char*)&tmp, sizeof(int));
		if(tmp != magicNumber2){
			UG_THROW("ERROR in RedistributeGrid: "
					 "Magic number mismatch after deserialization.\n");
		}

		UG_DLOG(LIB_GRID, 2, "Deserialization from rank " << recvFromRanks[i] << " done\n");
	}
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();

//	DEBUG: output distInfos...
	#ifdef LG_DISTRIBUTION_DEBUG
	{
		SaveDistInfosToFile(mg, distInfos, "dist_infos_after_distribution");
	}
	#endif

////////////////////////////////
//	CREATE LAYOUTS
	GDIST_PROFILE(gdist_CreateLayouts);
	CreateLayoutsFromDistInfos<VertexBase>(mg, glm, distInfos, aGeomObjID);
	CreateLayoutsFromDistInfos<EdgeBase>(mg, glm, distInfos, aGeomObjID);
	CreateLayoutsFromDistInfos<Face>(mg, glm, distInfos, aGeomObjID);
	CreateLayoutsFromDistInfos<Volume>(mg, glm, distInfos, aGeomObjID);
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();

////////////////////////////////
//	UPDATE THE DISTRIBUTED GRID MANAGER
	GDIST_PROFILE(gdist_UpdateDistGridManager);
	UG_DLOG(LIB_GRID, 2, "dist-DistributeGrid: Update DistributedGridManager\n");
	glm.remove_empty_interfaces();
	distGridMgr.enable_interface_management(true);
	distGridMgr.grid_layouts_changed(false);
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();

	mg.detach_from_all(aLocalInd);

	#ifdef LG_DISTRIBUTION_DEBUG
		PerformValidityCheck(distGridMgr);
	#endif

////	DEBUGGING...
//	{
//		static int counter = 0;
//		stringstream ss;
//		ss << "parallel-grid-layout-after-redist-" << counter << "-p" << pcl::ProcRank() << ".ugx";
//		UG_LOG("DEBUG SAVE OF PARALLEL GRID LAYOUT IN DistributeGrid\n");
//		SaveParallelGridLayout(mg, ss.str().c_str(), 0.1);
//		++counter;
//
//		if(!TestGridLayoutMap(mg, glm)){
//			UG_THROW("TestGridLayoutMap failed after redistribution!");
//		}
//	}



//	execute callbacks for external postprocessing
	GDIST_PROFILE(gdist_ExternalPostProcessing);
	UG_DLOG(LIB_GRID, 2, "dist: Informing msg-hub that distribution stops\n");

	msgHub->post_message(GridMessage_Creation(GMCT_CREATION_STOPS));
	//msgHub->post_message(GridMessage_Distribution(GMDT_GRID_SERIALIZATION_DONE));

//	we'll inform deserializers now, that deserialization is complete.
	distInfoSerializer.deserialization_done();
	serializer.deserialization_done();
	userDataSerializer.deserialization_done();
	msgHub->post_message(GridMessage_Distribution(GMDT_DISTRIBUTION_STOPS, userDataSerializer));
	//msgHub->post_message(GridMessage_Distribution(GMDT_DATA_SERIALIZATION_DONE));
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();

	UG_DLOG(LIB_GRID, 1, "dist-stop: DistributeGrid\n");
	return true;
}

}// end of namespace
