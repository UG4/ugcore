//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	year 2007
#ifndef	__H__HASH__
#define	__H__HASH__

#include <vector>
#include <list>
#include <cassert>

namespace ug
{

///	this template function creates a hash key for a number value. Specify to support other types
template <typename TKey> unsigned long hash_key(const TKey& key)
{
	return (unsigned long)key;
}

///	a generic Hash class
template <class TVal, class TKey> class Hash
{
	private:
		class HashEntry;
		typedef typename std::list<HashEntry>::iterator HashEntryIterator;
		typedef typename std::list<HashEntry>::const_iterator ConstHashEntryIterator;

	public:
		class Iterator;
		class ConstIterator;

		friend class ConstIterator;

	public:
		Hash()	{init(97);}
		Hash(unsigned long listSize)	{init(listSize);}

	///	warning: set_hash_size deletes all content
	//TODO: allow hash-resize without content-loss
		void set_hash_size(int hashSize)	{m_v.clear(); init(hashSize);}
		int get_hash_size()	const			{return m_listSize;}
		
		void clear()
		{
			m_v.clear();
			init(m_listSize);
		}

		void add(const TVal& val, const TKey& key)
		{
			unsigned long hKey = hash_key(key);
			m_v[hKey % m_listSize].push_back(HashEntry(hKey, val, key));
		}

		void erase(const TKey& key)
		{
			unsigned long hKey = hash_key(key);
			unsigned long hIndex = hKey % m_listSize;
			typename std::list<HashEntry>::iterator iter = m_v[hIndex].begin();
			typename std::list<HashEntry>::iterator iEnd = m_v[hIndex].end();

			while(iter != iEnd)
			{
				typename std::list<HashEntry>::iterator tIter = iter++;
				if((*tIter).hashKey == hKey)
				{
					if((*tIter).key == key)
						m_v[hIndex].erase(tIter);
				}
			}
		}

		bool has_entries(const TKey& key) const
		{
			unsigned long hKey = hash_key(key);
			unsigned long hIndex = hKey % m_listSize;
			return find_next_valid(m_v[hIndex].begin(), key, hIndex) != m_v[hIndex].end();
		}

		TVal& first(const TKey& key)
		{
			unsigned long hKey = hash_key(key);
			unsigned long hIndex = hKey % m_listSize;
			HashEntryIterator iter = find_next_valid(m_v[hIndex].begin(), key, hIndex);
			assert((iter != m_v[hIndex].end()) && "ERROR in Hash::first(...): Hash does not contain an entry with the given key.");
			return (*iter).value;
		}

		const TVal& first(const TKey& key) const
		{
			unsigned long hKey = hash_key(key);
			unsigned long hIndex = hKey % m_listSize;
			ConstHashEntryIterator iter = find_next_valid(m_v[hIndex].begin(), key, hIndex);
			assert((iter != m_v[hIndex].end()) && "ERROR in Hash::first(...): Hash does not contain an entry with the given key.");
			return (*iter).value;
		}

		Iterator begin(const TKey& key)
		{
			Iterator iter;
			unsigned long hKey = hash_key(key);
			unsigned long hIndex = hKey % m_listSize;
			iter.m_entryIter = find_next_valid(m_v[hIndex].begin(), key, hIndex);
			iter.m_key = key;
			iter.m_hIndex = hIndex;
			iter.m_pHash = this;
			return iter;
		}

		Iterator end(const TKey& key)
		{
			Iterator iter;
			unsigned long hKey = hash_key(key);
			unsigned long hIndex = hKey % m_listSize;
			iter.m_entryIter = m_v[hIndex].end();
			iter.m_key = key;
			iter.m_hIndex = hIndex;
			iter.m_pHash = this;
			return iter;
		}

		void advance_iterator(Iterator& iter)
		{
			iter.m_entryIter++;
			iter.m_entryIter = find_next_valid(iter.m_entryIter, iter.m_key, iter.m_hIndex);
		}

	//	const-version
		ConstIterator begin(const TKey& key) const
		{
			unsigned long hKey = hash_key(key);
			unsigned long hIndex = hKey % m_listSize;
			return ConstIterator(this, find_next_valid(m_v[hIndex].begin(), key, hIndex), key, hIndex);
		}

	//	const-version
		ConstIterator end(const TKey& key) const
		{
			unsigned long hKey = hash_key(key);
			unsigned long hIndex = hKey % m_listSize;
			return ConstIterator(this, m_v[hIndex].end(), key, hIndex);
		}

	private:
		void init(unsigned long listSize)
		{
			m_listSize = listSize;
			m_v.resize(listSize);
		}

		///	returns iter if iter is valid
		HashEntryIterator find_next_valid(HashEntryIterator iter, const TKey& key, unsigned long hIndex)
		{
			HashEntryIterator iEnd = m_v[hIndex].end();
			while(iter != iEnd)
			{
				if((*iter).key == key)
					return iter;
				iter++;
			}
			return iEnd;
		}

	//	const-version
		ConstHashEntryIterator find_next_valid(ConstHashEntryIterator iter, const TKey& key, unsigned long hIndex) const
		{
			ConstHashEntryIterator iEnd = m_v[hIndex].end();
			while(iter != iEnd)
			{
				if((*iter).key == key)
					return iter;
				iter++;
			}
			return iEnd;
		}

	private:
		std::vector<std::list<HashEntry> >			m_v;
		int	m_listSize;

	////////////////////////////////
	//	class definitions
	private:
		class HashEntry
		{
			public:
				HashEntry()	{};
				HashEntry(unsigned long hk, const TVal& v, const TKey& k) : hashKey(hk), key(k), value(v)	{}

				unsigned long		hashKey;
				TKey	key;
				TVal	value;
		};

	public:
		class Iterator
		{
			//friend class Hash;
			public:
				Iterator operator ++()	{m_pHash->advance_iterator(*this); return *this;}
				Iterator operator ++(int unused)	{Iterator i = *this; m_pHash->advance_iterator(*this); return i;}

				bool operator ==(const Iterator& iter) const
				{
					if(m_entryIter == iter.m_entryIter
						&& m_pHash == iter.m_pHash
						&& m_key == iter.m_key)
						return true;
					return false;
				}

				bool operator !=(const Iterator& iter) const
				{
					if(m_entryIter != iter.m_entryIter
						|| m_pHash != iter.m_pHash
						|| m_key != iter.m_key)
						return true;
					return false;
				}

				TVal& operator *()	{return (*m_entryIter).value;}
				const TVal& operator *()const	{return (*m_entryIter).value;}

			protected:
				HashEntryIterator		m_entryIter;
				Hash*					m_pHash;
				unsigned long			m_key;
				unsigned long			m_hIndex;
		};

		class ConstIterator
		{
			friend class Hash;
			public:
				ConstIterator()	{}
				ConstIterator(const Hash* hash, ConstHashEntryIterator entryIter, unsigned long key, unsigned long hIndex) :
							m_pHash(hash), m_entryIter(entryIter), m_key(key), m_hIndex(hIndex)
							{}

				ConstIterator operator ++() const			{advance(); return *this;}
				ConstIterator operator ++(int unused) const	{ConstIterator i = *this; advance(); return i;}

				void advance() const
				{
					m_entryIter++;
					m_entryIter = m_pHash->find_next_valid(m_entryIter, m_key, m_hIndex);
				}

				bool operator ==(const Iterator& iter) const
				{
					if(m_entryIter == iter.m_entryIter
						&& m_pHash == iter.m_pHash
						&& m_key == iter.m_key)
						return true;
					return false;
				}

				const TVal& operator *()const	{return (*m_entryIter).value;}

				bool operator !=(const ConstIterator& iter) const
				{
					if(m_entryIter != iter.m_entryIter
						|| m_pHash != iter.m_pHash
						|| m_key != iter.m_key)
						return true;
					return false;
				}

			protected:
				const Hash*						m_pHash;
				mutable ConstHashEntryIterator	m_entryIter;
				unsigned long					m_key;
				unsigned long					m_hIndex;
		};
};

}//	end of namespace
#endif
