#ifndef __H__UNIQUE_IDENTITY__
#define __H__UNIQUE_IDENTITY__

namespace ug
{

/// \addtogroup ugbase_common_util
/// \{

////////////////////////////////////////////////////////////////////////////////////////////////
//	UID
///	supplies a unique ID.
/**
 * Derivatives of this class can be identified by their unique ID.
 * Useful for hashing and other applications, where the need for identification is given.
 */
class UID
{
	public:
		UID()
			{
				static unsigned int ID = 1;
				m_uID = ID;
				ID++;
			}

		UID(const UID& uid)					{m_uID = uid.id();}
		inline unsigned int id() const		{return m_uID;}

		inline bool operator == (const UID& uid)
			{
				if(uid.id() == this->id())
					return true;
				return false;
			}

		virtual ~UID() {};

	private:
		unsigned int m_uID;
};

// end group ugbase_common_util
/// \}

}
#endif

