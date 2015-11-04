#ifndef __H__UG__owned_pointer__
#define __H__UG__owned_pointer__

namespace ug
{

/// \addtogroup ugbase_common_util
/// \{

///	Holds and automatically deletes a pointer, similar to std::auto_ptr. USE WITH CARE!
/**	WARNING: USE WITH CARE!
 * Class shows uncommon behavior when copied.
 *
 * The behavior of this class is very similar to std::auto_ptr and somewhat
 * similar to boost::scoped_ptr. It takes a pointer in its constructor and
 * automatically deletes the associated object in its destructor or in a call to reset.
 * The most notable property of OwnedPtr is that an OwnedPtr 'owns' the associated
 * object. That means that two OwnedPtrs should never point to the same instance.
 * The copy-constructor thus shows the somewhat irritating behavior of invalidating
 * the original copy and transferring the ownership of the associated object to
 * the new copy of OwnedPtr. This behavior again is the same in std::auto_ptr,
 * however, in contrary to std::auto_ptr, instances of OwnedPtr can be used
 * in standard containers like std::vector. This is explicitly not advised, since
 * many algorithms working with iterators or standard containers will not work
 * with OwnedPtrs.
 * Storing OwnedPtrs in lib_grids attachments may sometimes be useful. Therefore
 * this slightly dangerous class was introduced.
 */
template <class T>
class OwnedPtr
{
	public:
		typedef T* TPtr;
		typedef T& TRef;

		OwnedPtr(TPtr p = 0) : m_p(p)	{}

	///	Transfers ownership of the associated object from ap to this.
	/**	Note that ap looses ownership of the associated object.
	 * ap thus is not const at all, despite beeing declared const.*/
		OwnedPtr(const OwnedPtr& op)
		{
			m_p = op.get();
			op.invalidate();
		}

	///	Transfers ownership of the associated object from ap to this.
	/**	Note that ap looses ownership of the associated object.
	 * ap thus is not const at all, despite beeing declared const.*/
		OwnedPtr& operator=(const OwnedPtr& op)
		{
			reset(op.get());
			op.invalidate();
			return *this;
		}

		void reset(TPtr p = 0)				{if(m_p) delete m_p; m_p = p;}

		TPtr operator ->() const			{return m_p;}
		TRef operator *() const				{return *m_p;}

		TPtr& get() const					{return m_p;}

		operator bool()	const				{return m_p != 0;}

	private:
	///	This method is a hack, to allow that only one instance owns the pointer
		void invalidate() const				{m_p = 0;}

	private:
		mutable TPtr m_p;
};

// end group ugbase_common_util
/// \}

}//	end of namespace

#endif
