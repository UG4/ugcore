//	created by Sebastian Reiter
//	s.b.reiter@gmail.com
//	april 2013

#ifndef __H__UG__REFINEMENT_PROJECTION_HANDLER__
#define __H__UG__REFINEMENT_PROJECTION_HANDLER__

#include "../refinement_callbacks.h"
#include "lib_grid/tools/subset_handler_interface.h"

namespace ug{

///	Register different refinement callbacks for different subsets
/**	Using set_callback you may specify a callback which will be executed,
 * if the parent element is contained in the specified subset.
 * You may also specify a default callback which is executed if no parent is
 * present, if the parent isn't assigned to a subset or if no callback has
 * been assigned to the parent's subset. By default the default callback is
 * RefinementCallbackLinear.
 *
 * Use IRefiner::set_refinement_callback to set an instance of this class as
 * the refinement callback of a refiner.
 */
template <class TAPosition>
class RefinementProjectionHandler : public IRefinementCallback {
	public:
		RefinementProjectionHandler(SmartPtr<ISubsetHandler> sh, TAPosition aPos);

		virtual ~RefinementProjectionHandler();

		void set_default_callback(SmartPtr<IRefinementCallback> callback);
		void set_callback(int subsetIndex, SmartPtr<IRefinementCallback> callback);
		void set_callback(std::string subsetName, SmartPtr<IRefinementCallback> callback);

	////////////////////////////////////////
	//	IMPLEMENTATION OF IRefinementCallback
	///	called when a new vertex was created from an old vertex.
		virtual void new_vertex(VertexBase* vrt, VertexBase* parent);
	///	called when a new vertex was created from an old edge.
		virtual void new_vertex(VertexBase* vrt, EdgeBase* parent);
	///	called when a new vertex was created from an old face.
		virtual void new_vertex(VertexBase* vrt, Face* parent);
	///	called when a new vertex was created from an old volume.
		virtual void new_vertex(VertexBase* vrt, Volume* parent);

	///	callback for vertices in flat grids.
		virtual void flat_grid_vertex_encountered(VertexBase* vrt);

	///	returns the position of the given vertex.
		virtual int current_pos(number* coordsOut, VertexBase* vrt, int maxCoords);

	private:
		template <class TParent>
		void handle_new_vertex(VertexBase* vrt, TParent* parent);

		SmartPtr<ISubsetHandler>					m_sh;
		Grid::VertexAttachmentAccessor<TAPosition>	m_aaPos;
		std::vector<SmartPtr<IRefinementCallback> >	m_callbacks;
		SmartPtr<IRefinementCallback>				m_defaultCallback;
};

}	//	end of namespace

////////////////////
#include "refinement_projection_handler_impl.hpp"

#endif /* REFINEMENT_PROJECTION_HANDLER_H_ */
