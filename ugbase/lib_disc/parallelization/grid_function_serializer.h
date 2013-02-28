// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 27, 2013

#ifndef __H__UG__grid_function_serializer__
#define __H__UG__grid_function_serializer__

#include "lib_grid/algorithms/serialization.h"

namespace ug{

template <class TGridFct>
class GridFunctionSerializer : public GridDataSerializer
{
	public:
		static SPGridDataSerializer create(SmartPtr<TGridFct> fct)
		{return SPGridDataSerializer(new GridFunctionSerializer(fct));}

		GridFunctionSerializer()									{}

		GridFunctionSerializer(SmartPtr<TGridFct> fct) : m_fct(fct)	{}

		void set_function(SmartPtr<TGridFct> fct)					{m_fct = fct;}

		//virtual void write_info(BinaryBuffer& out) const;

		//virtual void read_info(BinaryBuffer& in);

		virtual void write_data(BinaryBuffer& out, VertexBase* o) const	{write(out, o);}
		virtual void write_data(BinaryBuffer& out, EdgeBase* o) const	{write(out, o);}
		virtual void write_data(BinaryBuffer& out, Face* o) const		{write(out, o);}
		virtual void write_data(BinaryBuffer& out, Volume* o) const		{write(out, o);}

		virtual void read_data(BinaryBuffer& in, VertexBase* o)			{read(in, o);}
		virtual void read_data(BinaryBuffer& in, EdgeBase* o)			{read(in, o);}
		virtual void read_data(BinaryBuffer& in, Face* o)				{read(in, o);}
		virtual void read_data(BinaryBuffer& in, Volume* o)				{read(in, o);}

	private:
		template <class TElem>
		void write(BinaryBuffer& out, TElem* e) const
		{
			std::vector<size_t>	indices;
			m_fct->inner_algebra_indices(e, indices);

			for(size_t i = 0; i < indices.size(); ++i){
				Serialize(out, (*m_fct)[indices[i]]);
			}
		}

		template <class TElem>
		void read(BinaryBuffer& in, TElem* e)
		{
			std::vector<size_t>	indices;
			m_fct->inner_algebra_indices(e, indices);

			for(size_t i = 0; i < indices.size(); ++i){
				Deserialize(in, (*m_fct)[indices[i]]);
			}
		}

	private:
		SmartPtr<TGridFct> m_fct;
};

}// end of namespace

#endif
