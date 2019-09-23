/*
 * interface_handler_local.h
 *
 *  Created on: 15.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_PROVIDER_H_
#define INTERFACE_PROVIDER_H_

namespace ug{


template <int TWorldDim>
class IInterfaceProvider
{
    public:
///	world Dimension
    static const int dim = TWorldDim;
    
/// default constructor:
	IInterfaceProvider(){};

/// destructor

	virtual ~IInterfaceProvider() {}

    virtual number get_LSvalue_byPosition(MathVector<dim> vrtPos, const int prtIndex) = 0;
    virtual const int get_orientation() const= 0;
    virtual void set_orientation(const int orientation) = 0;

};



}// end namespace ug


#endif /* INTERFACE_PROVIDER_H_ */
