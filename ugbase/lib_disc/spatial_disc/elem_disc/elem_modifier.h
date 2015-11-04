
#ifndef ELEM_MODIFIER_H_
#define ELEM_MODIFIER_H_

#include "elem_disc_interface.h"

namespace ug{

template <typename TDomain>
class IElemDisc;


template <typename TDomain>
class IElemDiscModifier
{

	protected:
	///	own type
		typedef IElemDiscModifier<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = TDomain::dim;

	public:
	///	Constructor (setting default values)
	/// \{
		IElemDiscModifier(): m_pElemDisc(NULL){};
		IElemDiscModifier(IElemDisc<TDomain>* myElemDisc) : m_pElemDisc(myElemDisc) {};
	 /// \}

	/// Virtual destructor
		virtual ~IElemDiscModifier(){}


	/// virtual initiates pre-computations before the standard element assembling
		virtual void preprocess(LocalVector& u, LocalVector& d, LocalVector& tmpD, GridObject* elem,
								MathVector<dim> vCornerCoords[], LocalIndices& ind);

	/// virtual initiates pre-computations before the standard element assembling
		virtual void preprocess(LocalVector& u, LocalMatrix& J, GridObject* elem,
								MathVector<dim> vCornerCoords[], LocalIndices& ind);

	/// virtual initiates post-computations after the standard element assembling
		virtual void postprocess(const LocalVector& u, LocalVector& d, LocalIndices& ind);

	/// virtual initiates post-computations after the standard element assembling
		virtual void postprocess(const LocalVector& u, LocalMatrix& J, LocalIndices& ind);

        void set_elem_disc(IElemDisc<TDomain>* myElemDisc){ m_pElemDisc = myElemDisc; }


      protected:
        IElemDisc<TDomain>* m_pElemDisc;
};

/*
 IElemDiscModifier_Local : IElemDiscModifier

 u.acces_by_map(...);
*/
} // end name space ug

#include "elem_modifier_impl.h"


#endif /* ELEM_MODIFIER_H_ */
