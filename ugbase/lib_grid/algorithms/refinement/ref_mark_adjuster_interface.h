#ifndef __H__UG__refmark_adjuster_interface__
#define __H__UG__refmark_adjuster_interface__

#include "refiner_interface.h"
#include "common/util/smart_pointer.h"

namespace ug{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

/**	Used to steer the mark-adjustment process in classes derived from IRefiner.
 * Instances can be registered at IRefiner and are then used to mark new elements
 * based on already marked elements. Different specializations can thereby implement
 * different rules (e.g. for parallel refinement or anisotropic refinement).
 * The different specialization should be constructed orthogonal to each other,
 * so that they can be chained.
 *
 * The callback-methods 'ref_marks_changed' and 'coarsen_marks_changed' are called
 * with lists of elements which have been marked since the last adjustment-round.
 * The callbacks can then use IRefiner::mark to adjust marks as required.
 * If a callback changes a ref-mark on an element,
 * this will automatically trigger a new round of adjustments
 * (then only considering elements which have been marked in the current round).
 *
 * Adjustment normally stops as soon as no adjuster marked any new elements.
 */
class IRefMarkAdjuster
{
	public:
		IRefMarkAdjuster() :
			m_nodeDependencyOrder1(true)
		{}

		virtual ~IRefMarkAdjuster()	{}

		virtual void ref_marks_changed(IRefiner& ref,
										const std::vector<Vertex*>& vrts,
										const std::vector<Edge*>& edges,
										const std::vector<Face*>& faces,
										const std::vector<Volume*>& vols)
		{}

		virtual void coarsen_marks_changed(IRefiner& ref,
										const std::vector<Vertex*>& vrts,
										const std::vector<Edge*>& edges,
										const std::vector<Face*>& faces,
										const std::vector<Volume*>& vols)
		{return;}

	///	enables or disables node-dependency-order-1.
	/**	\{
	 * If enabled, hanging nodes may only depend on non-hanging nodes.
	 * An edge containing a hanging node thus will not have a hanging-node
	 * as a corner vertex.
	 *
	 * Enabled by default.*/
		void enable_node_dependency_order_1(bool bEnable)	{m_nodeDependencyOrder1 = bEnable;}
		bool node_dependency_order_1_enabled()				{return m_nodeDependencyOrder1;}
	/**	\} */

	private:
		bool	m_nodeDependencyOrder1;
};

typedef SmartPtr<IRefMarkAdjuster>	SPIRefMarkAdjuster;

/// @}

}// end of namespace

#endif
