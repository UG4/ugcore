// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_horizontal_anisotropy_adjuster
#define __H__UG_horizontal_anisotropy_adjuster

#include "../ref_mark_adjuster_interface.h"

namespace ug{

///	Selects additional edges to preserve layers in a grid with horizontal anisotropies
/** Selects additional edges to preserve layers in a grid with horizontal anisotropies
 * and builds a closure so that no hanging nodes are generated.
 * 
 * \note	This adjuster regards the grid as a serial grid. If the grid represents a part
 * 			of a distributed grid, then the additional use of a parallel adjuster is required.
 */
template <class TAPos>
class HorizontalAnisotropyAdjuster : public IRefMarkAdjuster
{
	public:
		typedef TAPos									position_attachment_t;
		typedef Grid::VertexAttachmentAccessor<TAPos>	position_accessor_t;

		HorizontalAnisotropyAdjuster(TAPos aPos) : m_aPos(aPos) {}

		virtual ~HorizontalAnisotropyAdjuster()	{}

		virtual void ref_marks_changed(IRefiner& ref,
										const std::vector<Vertex*>& vrts,
										const std::vector<Edge*>& edges,
										const std::vector<Face*>& faces,
										const std::vector<Volume*>& vols);
	private:
		TAPos	m_aPos;
};

}//	end of namespace

#endif	//__H__UG_horizontal_anisotropy_adjuster
