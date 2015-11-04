#ifndef __H__UG__fractured_media_refiner__
#define __H__UG__fractured_media_refiner__

#include <queue>
#include "hanging_node_refiner_t.h"

namespace ug
{

/**	This class takes special care for degenerated fractures during refinement.
 * Degenerated faces are refined using anisotropic refinement, so that their
 * degenerated sides are not refined.
 *
 * Currently this only works in 2d.
 *
 * \todo	Add support for degenerated volumes
 * \todo	Use a IsDegenerated callback instead of thresholds
 */
template <class TGrid, class TAPosition>
class FracturedMediaRefiner : public THangingNodeRefiner<TGrid>
{
	typedef THangingNodeRefiner<TGrid>	BaseClass;

	public:
		using BaseClass::mark;

	public:
		FracturedMediaRefiner(IRefinementCallback* refCallback = NULL);
		FracturedMediaRefiner(TGrid& g,
							  IRefinementCallback* refCallback = NULL);

		virtual ~FracturedMediaRefiner();

	///	if the aspect ratio is smaller then the given threshold, the element is considered a fracture element.
		void set_aspect_ratio_threshold(number threshold);

	///	\todo: replace this with a callback
		void set_position_attachment(TAPosition& aPos);

	///	Marks a face for refinement.
	/**	Uses the degenerated-edge-threshold to determine, whether the face
	 * is a fracture-face or not.
	 * Make sure to specify a position attachment before marking any elements.*/
		virtual bool mark(Face* f, RefinementMark refMark = RM_REFINE);

	protected:
		number aspect_ratio(Face* f);
		virtual void collect_objects_for_refine();

	private:
		Grid::VertexAttachmentAccessor<TAPosition>	m_aaPos;
		std::queue<Face*>	m_queDegeneratedFaces;
		number				m_aspectRatioThreshold;
};

}//	end of namespace

#endif
