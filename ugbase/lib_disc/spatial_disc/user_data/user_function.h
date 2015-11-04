
#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__USER_FUNCTION__
#define __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__USER_FUNCTION__

namespace ug{

template <typename TData, typename TDataIn = TData>
class IFunction
{
	public:
	///	evaluates the data
		virtual void operator() (TData& out, int numArgs, ...) = 0;

	///	virtual destructor
		virtual ~IFunction() {}
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__USER_FUNCTION__ */
