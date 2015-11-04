
#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON_UPDATE_INTERFACE__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON_UPDATE_INTERFACE__

namespace ug {

/// general interface for data updates during Newton process
class INewtonUpdate
{
	public:
		virtual void update() = 0;
		virtual ~INewtonUpdate() {};
};

}

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON_UPDATE_INTERFACE__ */
