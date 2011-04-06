// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 24.01.2011 (m,d,y)

#ifndef __H__UG__REFINER_INTERFACE__
#define __H__UG__REFINER_INTERFACE__

#include "refinement_callbacks.h"

namespace ug
{

///	The refiner interface allows to mark elements for refinement and to call refine.
class IRefiner
{
	public:
		IRefiner(IRefinementCallback* refCallback = NULL) :
			m_refCallback(refCallback)	{}

		virtual ~IRefiner()	{}

		void set_refinement_callback(IRefinementCallback* refCallback)
			{m_refCallback = refCallback;}

		IRefinementCallback* get_refinement_callback()
			{return m_refCallback;}

	///	has to return the associated grid. Pure virtual
		virtual Grid* get_associated_grid() = 0;

	///	clears all marks. Default implementation is empty
		virtual void clear_marks()	{}

	///	Marks a vertex for refinement. Default implementation is empty
		virtual void mark_for_refinement(VertexBase* v)	{}

	///	Marks an edge for refinement. Default implementation is empty
		virtual void mark_for_refinement(EdgeBase* e)	{}

	///	Marks a face for refinement. Default implementation is empty
		virtual void mark_for_refinement(Face* f)		{}

	///	Marks a volume for refinement. Default implementation is empty
		virtual void mark_for_refinement(Volume* v)		{}

	///	marks all elements between iterBegin and iterEnd.
	/**	the value-type of TIterator has to be a pointer to a type derived
	 * 	from either EdgeBase, Face or Volume.*/
		template <class TIterator>
		void mark_for_refinement(const TIterator& iterBegin, const TIterator& iterEnd)
			{
				TIterator iter = iterBegin;
				while(iter != iterEnd){
					mark_for_refinement(*iter);
					++iter;
				}
			}

	/// Performs refinement on the marked elements. Pure virtual.
		virtual void refine() = 0;

	protected:
		IRefinementCallback*	m_refCallback;
};

}//	end of namespace

#endif
