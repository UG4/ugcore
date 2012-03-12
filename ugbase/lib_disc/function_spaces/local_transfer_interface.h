/*
 * local_transfer_interface.h
 *
 *  Created on: 07.03.2012
 *      Author: andreasvogel
 */

#ifndef LOCAL_TRANSFER_INTERFACE_H_
#define LOCAL_TRANSFER_INTERFACE_H_

namespace ug{

// predeclaration
class MGDoFDistribution;

class ILocalTransfer
{
	public:
	///	constructor
		ILocalTransfer() {}

	///	virtual destructor
		virtual ~ILocalTransfer() {}

	///	returns if prolongation is performed on type
		virtual bool prolongation_needed(GeometricBaseObject gbo) const = 0;

	///	returns if restriction is performed on type
		virtual bool restriction_needed(GeometricBaseObject gbo) const = 0;

	///	prolongate
		virtual void prolongate_values(VertexBase* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const = 0;
		virtual void prolongate_values(EdgeBase* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const = 0;
		virtual void prolongate_values(Face* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const = 0;
		virtual void prolongate_values(Volume* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const = 0;

	///	restrict
		virtual void restrict_values(VertexBase* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const = 0;
		virtual void restrict_values(EdgeBase* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const = 0;
		virtual void restrict_values(Face* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const = 0;
		virtual void restrict_values(Volume* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const = 0;

};

} // end namespace ug

#endif /* LOCAL_TRANSFER_INTERFACE_H_ */
