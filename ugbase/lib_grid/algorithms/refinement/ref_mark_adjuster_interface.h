// created by Sebastian Reiter
// s.b.reiter@gmail.com
// 15.11.2012 (d,m,y)

#ifndef __H__UG__refmark_adjuster_interface__
#define __H__UG__refmark_adjuster_interface__

namespace ug{

class IRefMarkAdjuster
{
	public:
		virtual ~IRefMarkAdjuster()	{}

		virtual void ref_marks_changed(IRefiner& ref,
									   const std::vector<VertexBase*>& vrts,
									   const std::vector<EdgeBase*>& edge,
									   const std::vector<Face*>& face,
									   const std::vector<Volume*>& vol);
};

}// end of namespace

#endif
