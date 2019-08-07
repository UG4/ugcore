/*
 * interface_handler_local.h
 *
 *  Created on: 15.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_PROVIDER_H_
#define INTERFACE_PROVIDER_H_

namespace ug{

    
class IInterfaceProvider
{
	public:

/// default constructor:
	IInterfaceProvider(){};

/// destructor

	~IInterfaceProvider() {}

};


class IInterfaceProvider
{
	public:

/// default constructor:
	IInterfaceProvider(){};

/// destructor

	~IInterfaceProvider() {}

};


template <int TWorldDim>
class InterfaceProviderBase : public IInterfaceProvider
{

	public:
///	world Dimension
	static const int dim = TWorldDim;

/// default constructor:
	InterfaceProviderBase()
	{
		clear();
		UG_LOG("InterfaceProviderBase constructor\n");
	};

/// destructor
	~InterfaceProviderBase() {}

};


}// end namespace ug


#endif /* INTERFACE_PROVIDER_H_ */
