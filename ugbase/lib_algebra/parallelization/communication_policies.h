/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Andreas Vogel, Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__COMMUNICATION_POLICIES__
#define __H__LIB_ALGEBRA__PARALLELIZATION__COMMUNICATION_POLICIES__

#include "pcl/pcl.h"
#include "common/serialization.h"
#include "parallel_index_layout.h"
#include "algebra_id.h"


namespace ug{

// predeclaration of block_traits
/**	The block-traits are required by communication policies to determine whether buffers
 *	have a fixed size or whether they have to be treated in a more flexible way.
 *	such traits could look like this (default implementation):
 *	\code
 *	template <> struct block_traits<double>
 *	{
 *		enum{
 *			is_static = 1
 *		};
 *	};
 *	\endcode
 */
template <typename t> struct block_traits
{
	enum{
		is_static = 1
	};
};

/**
 * \brief Communication Policies for parallel Algebra
 *
 * Algebra Communication Policies are the basic building blocks for the
 * parallel communication between Vectors and Matrices.
 *
 * \defgroup lib_algebra_parallelization_policies Parallel Algebra Communication Policies
 * \ingroup lib_algebra_parallelization
 */

/// \addtogroup lib_algebra_parallelization_policies
/// @{

/// Communication Policy to copy values of a vector
/**
 * This class is used as a policy to copy values on the interfaces of a
 * parallel vector. The collecting interfaces pack the values on the interface
 * into a stream. The extracting interfaces receive the stream and overwrites
 * the values of the vector with the sent data.
 *
 * \tparam	TVector		Vector type
 */
template <class TVector>
class ComPol_VecCopy : public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
	///	Default constructor
		ComPol_VecCopy() : m_pVecDest(NULL), m_pVecSrc(NULL) {}

	///	Constructor setting the vector
		ComPol_VecCopy(TVector* pVec): m_pVecDest(pVec), m_pVecSrc(pVec)	{}

	///	Constructor setting the vector
		ComPol_VecCopy(TVector* pVecDest, const TVector* pVecSrc)
			: m_pVecDest(pVecDest), m_pVecSrc(pVecSrc)	{}

	///	set the vector to work on
		void set_vector(TVector* pVec)	{m_pVecDest = pVec; m_pVecSrc = pVec;}

	///	set the vector to work on
		void set_vector(TVector* pVecDest, const TVector* pVecSrc)
			{m_pVecDest = pVecDest; m_pVecSrc = pVecSrc;}

	/// returns the buffer size
	/**
	 * This function returns the size of the buffer needed for the communication
	 * of passed interface. If the vector has fixed size entries this is just
	 * the number of interface entries times the size of the entry. In case
	 * of a variable size entry type a negative value is returned to indicate
	 * that no buffer size can be determined in advanced.
	 *
	 * \param[in]	interface	Interface that will communicate
	 */
		virtual int
		get_required_buffer_size(const Interface& interface)
		{
			if(block_traits<typename TVector::value_type>::is_static)
				return interface.size() * sizeof(typename TVector::value_type);
			else
				return -1;
		}

	///	writes the interface values into a buffer that will be sent
	/**
	 * This function collects all entries of the vector into a buffer that
	 * are part of the interface.
	 *
	 * \param[out]		buff		Buffer
	 * \param[in]		interface	Interface that will communicate
	 */
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_VecCopy_collect, "algebra parallelization");
		//	check that vector has been set
			if(m_pVecSrc == NULL) return false;

		//	rename for convenience
			const TVector& v = *m_pVecSrc;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	write entry into buffer
				Serialize(buff, v[index]);
			}
			return true;
		}

	///	writes values from a buffer into the interface
	/**
	 * This function writes the buffer values into the vector.
	 *
	 * \param[out]		buff		Buffer
	 * \param[in]		interface	Interface that communicates
	 */
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_VecCopy_extract, "algebra parallelization");
		//	check that vector has been set
			if(m_pVecDest == NULL) return false;

		//	rename for convenience
			TVector& v = *m_pVecDest;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	write entry into vector
				Deserialize(buff, v[index]);
			}
			return true;
		}

	private:
	//	pointer to current vector
		TVector* m_pVecDest;
		const TVector* m_pVecSrc;
};

/// Communication Policy to copy scaled values of a vector
/**
 * This class is used as a policy to copy scaled values on the interfaces of a
 * parallel vector. The collecting interfaces pack the values on the interface
 * into a stream. The extracting interfaces receive the stream and overwrites
 * the values of the vector with the sent data times the scaling.
 *
 * \tparam	TVector		Vector type
 */
template <class TVector>
class ComPol_VecScaleCopy : public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
	///	Default Constructor
		ComPol_VecScaleCopy() : m_pVec(NULL), m_scale(1.0)	{}

	///	Constructor setting vector and scaling factor
		ComPol_VecScaleCopy(TVector* pVec, number scale)
			: m_pVec(pVec), m_scale(scale)	{}

	///	sets the vector that we be used for communication
		void set_vector(TVector* pVec)	{m_pVec = pVec;}

	///	sets the scaling factor
		void set_scale(number scale)	{m_scale = scale;}

	/// returns the buffer size
	/**
	 * This function returns the size of the buffer needed for the communication
	 * of passed interface. If the vector has fixed size entries this is just
	 * the number of interface entries times the size of the entry. In case
	 * of a variable size entry type a negative value is returned to indicate
	 * that no buffer size can be determined in advanced.
	 *
	 * \param[in]	interface	Interface that will communicate
	 */
		virtual int
		get_required_buffer_size(const Interface& interface)
		{
			if(block_traits<typename TVector::value_type>::is_static)
				return interface.size() * sizeof(typename TVector::value_type);
			else
				return -1;
		}

	///	writes the interface values into a buffer that will be sent
	/**
	 * This function collects all entries of the vector into a buffer that
	 * are part of the interface.
	 *
	 * \param[out]		buff		Buffer
	 * \param[in]		interface	Interface that will communicate
	 */
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_VecScaleCopy_collect, "algebra parallelization");
		//	check that vector has been set
			if(m_pVec == NULL) return false;

		//	rename for convenience
			TVector& v = *m_pVec;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	copy value into buffer
				Serialize(buff, v[index]);
			}
			return true;
		}

	///	scales values of a buffer and writes the into the interface
	/**
	 * This function writes the scaled buffer values into the vector.
	 *
	 * \param[out]		buff		Buffer
	 * \param[in]		interface	Interface that communicates
	 */
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_VecScaleCopy_extract, "algebra parallelization");
		//	check that vector has been set
			if(m_pVec == NULL) return false;

		//	rename for convenience
			TVector& v = *m_pVec;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	copy value from buffer
				Deserialize(buff, v[index]);

			//	scale value
				v[index] *= m_scale;
			}
			return true;
		}

	private:
		TVector* m_pVec;

		number m_scale;
};

/// Communication Policy to add values of a vector
/**
 * This class is used as a policy to add values on the interfaces of a
 * parallel vector. The collecting interfaces pack the values on the interface
 * into a stream. The extracting interfaces receive the stream and adds
 * the values of the vector.
 *
 * \tparam	TVector		Vector type
 */
template <class TVector>
class ComPol_VecAdd : public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
	///	Default constructor
		ComPol_VecAdd() : m_pVecDest(NULL), m_pVecSrc(NULL)	{}

	///	Constructor setting the values
		ComPol_VecAdd(TVector* pVec) : m_pVecDest(pVec),  m_pVecSrc(pVec)	{}

	///	Constructor setting the values
		ComPol_VecAdd(TVector* pVecDest, const TVector* pVecSrc)
			: m_pVecDest(pVecDest),  m_pVecSrc(pVecSrc)	{}

	///	sets the vector used in communication
		void set_vector(TVector* pVec)	{m_pVecDest = pVec; m_pVecSrc = pVec;}

	///	sets the vector used in communication
		void set_vector(TVector* pVecDest, const TVector* pVecSrc)
			{m_pVecDest = pVecDest; m_pVecSrc = pVecSrc;}

	/// returns the buffer size
	/**
	 * This function returns the size of the buffer needed for the communication
	 * of passed interface. If the vector has fixed size entries this is just
	 * the number of interface entries times the size of the entry. In case
	 * of a variable size entry type a negative value is returned to indicate
	 * that no buffer size can be determined in advanced.
	 *
	 * \param[in]	interface	Interface that will communicate
	 */
		virtual int
		get_required_buffer_size(const Interface& interface)
		{
			if(block_traits<typename TVector::value_type>::is_static)
				return interface.size() * sizeof(typename TVector::value_type);
			else
				return -1;
		}

	///	writes the interface values into a buffer that will be sent
	/**
	 * This function collects all entries of the vector into a buffer that
	 * are part of the interface.
	 *
	 * \param[out]		buff		Buffer
	 * \param[in]		interface	Interface that will communicate
	 */
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_VecAdd_collect, "algebra parallelization");
		//	check that vector has been set
			if(m_pVecSrc == NULL) return false;

		//	rename for convenience
			const TVector& v = *m_pVecSrc;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			// copy entry
				Serialize(buff, v[index]);
			}
			return true;
		}

	///	adds values of a buffer to the interface values
	/**
	 * This function adds the buffer values to the vector values.
	 *
	 * \param[out]		buff		Buffer
	 * \param[in]		interface	Interface that communicates
	 */
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_VecAdd_extract, "algebra parallelization");
		//	check that vector has been set
			if(m_pVecDest == NULL) return false;

		//	rename for convenience
			TVector& v = *m_pVecDest;

		//	entry
			typename TVector::value_type entry;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	copy entry
				Deserialize(buff, entry);

			//	add entry
				v[index] += entry;
			}
			return true;
		}

	private:
		TVector* m_pVecDest;
		const TVector* m_pVecSrc;
};

/// Communication Policy to add values of a vector
/**
 * This class is used as a policy to add values on the interfaces of a
 * parallel vector. The collecting interfaces pack the values on the interface
 * into a stream. The extracting interfaces receive the stream and adds
 * the values of the vector.
 *
 * \tparam	TVector		Vector type
 */
template <class TVector>
class ComPol_VecScaleAdd : public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
	///	Default constructor
		ComPol_VecScaleAdd() : m_pVec(NULL), m_scale(1.0) {}

	///	Constructor setting the values
		ComPol_VecScaleAdd(TVector* pVec) : m_pVec(pVec)	{}

	///	sets the vector used in communication
		void set_vector(TVector* pVec)	{m_pVec = pVec;}

	///	sets the scaling factor
		void set_scale(number scale)	{m_scale = scale;}

	/// returns the buffer size
	/**
	 * This function returns the size of the buffer needed for the communication
	 * of passed interface. If the vector has fixed size entries this is just
	 * the number of interface entries times the size of the entry. In case
	 * of a variable size entry type a negative value is returned to indicate
	 * that no buffer size can be determined in advanced.
	 *
	 * \param[in]	interface	Interface that will communicate
	 */
		virtual int
		get_required_buffer_size(const Interface& interface)
		{
			if(block_traits<typename TVector::value_type>::is_static)
				return interface.size() * sizeof(typename TVector::value_type);
			else
				return -1;
		}

	///	writes the interface values into a buffer that will be sent
	/**
	 * This function collects all entries of the vector into a buffer that
	 * are part of the interface.
	 *
	 * \param[out]		buff		Buffer
	 * \param[in]		interface	Interface that will communicate
	 */
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_VecScaleAdd_collect, "algebra parallelization");
		//	check that vector has been set
			if(m_pVec == NULL) return false;

		//	rename for convenience
			TVector& v = *m_pVec;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			// copy entry
				Serialize(buff, v[index]);
			}
			return true;
		}

	///	scales and adds values of a buffer to the interface values
	/**
	 * This function sclaes and adds the buffer values to the vector values.
	 *
	 * \param[out]		buff		Buffer
	 * \param[in]		interface	Interface that communicates
	 */
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_VecScaleAdd_extract, "algebra parallelization");
		//	check that vector has been set
			if(m_pVec == NULL) return false;

		//	rename for convenience
			TVector& v = *m_pVec;

		//	entry
			typename TVector::value_type entry;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	copy entry
				Deserialize(buff, entry);

			//	add entry
				v[index] += entry * m_scale;
			}
			return true;
		}

	private:
		TVector* m_pVec;

		number m_scale;
};


/// Communication Policy to add values of a vector and reset value to zero on sending interface
/**
 * This class is used as a policy to add values on the interfaces of a
 * parallel vector and set the values to zero on the sending interface.
 * The collecting routine packs the interface values of the given vector
 * into a stream and and sets those entries to zero immediately after packing.
 * The extracting interfaces receive the stream and add the values to those of
 * already existing in the vector.
 *
 * \note	If an entry is sent over multiple interfaces, only the first connected
 * 			proc will receive the full value. All others will receive 0. This helps
 * 			to keep the vector additive.
 *
 * \tparam	TVector		Vector type
 */
template <class TVector>
class ComPol_VecAddSetZero : public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
	///	Default Constructor
		ComPol_VecAddSetZero() : m_pVec(0)	{}

	///	Constructor setting vector
		ComPol_VecAddSetZero(TVector* pVec) : m_pVec(pVec)	{}

	///	sets the vector that we be used for communication
		void set_vector(TVector* pVec)	{m_pVec = pVec;}

	/// returns the buffer size
	/**
	 * This function returns the size of the buffer needed for the communication
	 * of passed interface. If the vector has fixed size entries this is just
	 * the number of interface entries times the size of the entry. In case
	 * of a variable size entry type a negative value is returned to indicate
	 * that no buffer size can be determined in advanced.
	 *
	 * \param[in]	interface	Interface that will communicate
	 */
		virtual int
		get_required_buffer_size(const Interface& interface)
		{
			if(block_traits<typename TVector::value_type>::is_static)
				return interface.size() * sizeof(typename TVector::value_type);
			else
				return -1;
		}

	///	writes the interface values into a buffer that will be sent and then sets
	/// the value to zero on the interface
	/**
	 * This function collects all entries of the vector into a buffer that
	 * are part of the interface. The values on the interface are then set to
	 * zero.
	 *
	 * \param[out]		buff		Buffer
	 * \param[in]		interface	Interface that will communicate
	 */
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_VecAddSetZero_collect, "algebra parallelization");
		//	check that vector has been set
			if(m_pVec == NULL) return false;

		//	rename for convenience
			TVector& v = *m_pVec;

			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			// copy values
				Serialize(buff, v[index]);
				v[index] *= 0;
			}

			return true;
		}

	///	adds the values of a buffer to the values on the interface
	/**
	 * This function adds the values of the buffer to the interface values.
	 *
	 * \param[out]		buff		Buffer
	 * \param[in]		interface	Interface that communicates
	 */
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_VecAddSetZero_extract, "algebra parallelization");
		//	check that vector has been set
			if(m_pVec == NULL) return false;

		//	rename for convenience
			TVector& v = *m_pVec;

		//	entry
			typename TVector::value_type entry;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	copy values
				Deserialize(buff, entry);

			//	add entry
				v[index] += entry;
			}
			return true;
		}

	private:
		TVector* m_pVec;
};

/// Communication Policy to subtract values of a vector
/**
 * This class is used as a policy to subtract values on the interfaces of a
 * parallel vector. The collecting interfaces pack the values on the interface
 * into a stream. The extracting interfaces receive the stream and subtracts
 * the values of the vector.
 *
 * \tparam	TVector		Vector type
 */
template <class TVector>
class ComPol_VecSubtract : public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
	///	Default constructor
		ComPol_VecSubtract() : m_pVec(NULL)	{}

	///	Constructor setting the values
		ComPol_VecSubtract(TVector* pVec) : m_pVec(pVec)	{}

	///	sets the vector used in communication
		void set_vector(TVector* pVec)	{m_pVec = pVec;}

	/// returns the buffer size
	/**
	 * This function returns the size of the buffer needed for the communication
	 * of passed interface. If the vector has fixed size entries this is just
	 * the number of interface entries times the size of the entry. In case
	 * of a variable size entry type a negative value is returned to indicate
	 * that no buffer size can be determined in advanced.
	 *
	 * \param[in]	interface	Interface that will communicate
	 */
		virtual int
		get_required_buffer_size(const Interface& interface)
		{
			if(block_traits<typename TVector::value_type>::is_static)
				return interface.size() * sizeof(typename TVector::value_type);
			else
				return -1;
		}

	///	writes the interface values into a buffer that will be sent
	/**
	 * This function collects all entries of the vector into a buffer that
	 * are part of the interface.
	 *
	 * \param[out]		buff		Buffer
	 * \param[in]		interface	Interface that will communicate
	 */
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_VecSubtract_collect, "algebra parallelization");
		//	check that vector has been set
			if(m_pVec == NULL) return false;

		//	rename for convenience
			TVector& v = *m_pVec;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			// get index
				const size_t index = interface.get_element(iter);

			// copy value
				Serialize(buff, v[index]);
			}
			return true;
		}

	///	subtracts values of a buffer to the interface values
	/**
	 * This function subtracts the buffer values to the vector values.
	 *
	 * \param[out]		buff		Buffer
	 * \param[in]		interface	Interface that communicates
	 */
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_VecSubtract_extract, "algebra parallelization");
		//	check that vector has been set
			if(m_pVec == NULL) return false;

		//	rename for convenience
			TVector& v = *m_pVec;

		// entry
			typename TVector::value_type entry;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	copy vector
				Deserialize(buff, entry);

			//	subtract entry
				v[index] -= entry;
			}
			return true;
		}

	private:
		TVector* m_pVec;
};

/// Communication Policy to check consistency of a vector
template <class TVector>
class ComPol_CheckConsistency : public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
	///	Default constructor
		ComPol_CheckConsistency() : m_pVec(NULL)	{}

	///	Constructor setting the values
		ComPol_CheckConsistency(const TVector* pVec) : m_pVec(pVec)	{}

	///	sets the vector used in communication
		void set_vector(const TVector* pVec)	{m_pVec = pVec;}

	/// returns the buffer size
		virtual int
		get_required_buffer_size(const Interface& interface)
		{
			if(block_traits<typename TVector::value_type>::is_static)
				return interface.size() * sizeof(typename TVector::value_type);
			else
				return -1;
		}

	///	writes the interface values into a buffer that will be sent
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_CheckConsistency_collect, "algebra parallelization");
		//	check that vector has been set
			if(m_pVec == NULL) return false;

		//	rename for convenience
			const TVector& v = *m_pVec;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			// get index
				const size_t index = interface.get_element(iter);

			// copy value
				Serialize(buff, v[index]);
			}
			return true;
		}

	///	checks if recieved values are equal to process local values
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_CheckConsistency_extract, "algebra parallelization");
		//	check that vector has been set
			if(m_pVec == NULL) return false;

		//	rename for convenience
			const TVector& v = *m_pVec;

		// entry
			typename TVector::value_type entry, diff;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	copy vector
				Deserialize(buff, entry);

				diff = v[index];
				diff -= entry;

			//	check
				if(VecNormSquared(diff) > 1e-20 * VecNormSquared(entry)){
					UG_LOG_ALL_PROCS("Slave and master value at index "<<index<<
					                 " differ by: "<<VecNormSquared(diff)<<"\n");
				}
			}

			return true;
		}

	private:
		const TVector* m_pVec;
};

/// Communication Policy to subtract only one slave value per master of a vector
/**
 * This class is used as a policy to subtract values on the interfaces of a
 * parallel vector. The collecting interfaces - supposed to consist of slave
 * dofs - pack the values on the interface into a stream. 
 * The extracting interfaces - supposed to consist of master dofs - receive
 * the stream and subtracts only one slave value (the first) per master.
 *
 * \tparam	TVector		Vector type
 */
template <class TVector>
class ComPol_VecSubtractOnlyOneSlave : public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
	///	Default constructor
		ComPol_VecSubtractOnlyOneSlave() : m_pVec(NULL)	{}

	///	Constructor setting the values
		ComPol_VecSubtractOnlyOneSlave(TVector* pVec) 
	  	{
			set_vector(pVec);
		}

	///	sets the vector used in communication
		void set_vector(TVector* pVec) 
		{
			m_pVec = pVec;
			m_vProcessed.resize(m_pVec->size(), false);
		}

	/// clear processed flag
		void clear()
		{
			m_vProcessed.clear();
			m_vProcessed.resize(m_pVec->size(), false);
		}

	/// returns the buffer size
	/**
	 * This function returns the size of the buffer needed for the communication
	 * of passed interface. If the vector has fixed size entries this is just
	 * the number of interface entries times the size of the entry. In case
	 * of a variable size entry type a negative value is returned to indicate
	 * that no buffer size can be determined in advanced.
	 *
	 * \param[in]	interface	Interface that will communicate
	 */
		virtual int
		get_required_buffer_size(const Interface& interface)
		{
			if(block_traits<typename TVector::value_type>::is_static)
				return interface.size() * sizeof(typename TVector::value_type);
			else
				return -1;
		}

	///	writes the interface values into a buffer that will be sent
	/**
	 * This function collects all entries of the vector into a buffer that
	 * are part of the interface.
	 *
	 * \param[out]		buff		Buffer
	 * \param[in]		interface	Interface that will communicate
	 */
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{

			PROFILE_BEGIN_GROUP(ComPol_VecSubtractOnlyOneSlave_collect, "algebra parallelization");
		//	check that vector has been set
			if(m_pVec == NULL) return false;

		//	rename for convenience
			TVector& v = *m_pVec;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			// get index
				const size_t index = interface.get_element(iter);

			// copy value
				Serialize(buff, v[index]);
			}
			return true;
		}

	///	subtracts values of a buffer to the interface values
	/**
	 * This function subtracts the buffer values to the vector values.
	 *
	 * \param[out]		buff		Buffer
	 * \param[in]		interface	Interface that communicates
	 */
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_VecSubtractOnlyOneSlave_extract, "algebra parallelization");
		//	check that vector has been set
			if(m_pVec == NULL) return false;

		//	rename for convenience
			TVector& v = *m_pVec;

		// entry
			typename TVector::value_type entry;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	copy vector
				Deserialize(buff, entry);

			//	subtract entry
				if(m_vProcessed[index] == false)
				{
					v[index] -= entry;
					m_vProcessed[index] = true;
				}
			}
			return true;
		}

	private:
		TVector* m_pVec;
		std::vector<bool> m_vProcessed;
};



/// Communication Policy sending fractions to
/**
 * This Communication policy assumes that the matrix is stored additively with
 * no overlap. It then copies the row entries of the slaves indices to the
 * corresponding master rows. Only those values are copied that have if the index
 * the slave couples to has also a representation on the other process. This means
 * that couplings to inner indices or to other slave/master, that have only
 * connections to other processes, are not taken into account.
 *
 * \tparam TMatrix	matrix type
 */
template <class TAlgebra>
class ComPol_MatDistributeDiag
	: public pcl::ICommunicationPolicy<IndexLayout>
{
	public:

	typedef typename TAlgebra::matrix_type TMatrix;
	typedef typename TAlgebra::vector_type TVector;

	///	Constructor setting the vector
	/**
	 * vGlobalID must have size >= mat.num_rows()
	 */
	ComPol_MatDistributeDiag(TMatrix& rMat, TVector &rWeight, double theta=1.0)
			: m_rMat(rMat), m_rWeight(rWeight), m_dTheta(theta)
		{
			UG_ASSERT( m_rMat.num_rows()==m_rWeight.size(), "sizes should match");
		}

	/// returns the buffer size
		/**
		 * This function returns the size of the buffer needed for the communication
		 * of passed interface. If the vector has fixed size entries this is just
		 * the number of interface entries times the size of the entry. In case
		 * of a variable size entry type a negative value is returned to indicate
		 * that no buffer size can be determined in advanced.
		 *
		 * \param[in]	interface	Interface that will communicate
		 */
		virtual int
		get_required_buffer_size(const Interface& interface)
			{
				if(block_traits<typename TMatrix::value_type>::is_static)
					return interface.size() * sizeof(typename TMatrix::value_type);
				else
					return -1;
			}


	///	writes the interface values into a buffer that will be sent
		virtual bool collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(CMatDistributeDiag_collect, "algebra parallelization");
			typedef typename TMatrix::value_type block_type;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);
				block_type& a_ii = m_rMat(index, index);

			// send rescaled diagonal entry
				block_type out_ii = a_ii;
				out_ii = out_ii * (1.0/BlockRef(m_rWeight[index], 0));
				//	write diagonal entry to stream
				Serialize(buff, out_ii);

			//	std::cerr << "Sending(" << index << ")=" << out_ii <<  " "<< BlockRef(m_rWeight[index], 0) << std::endl;
			// store average
				a_ii = out_ii;

			// after the first call, we need to reset the weight (on master)
			// in order to communicate correct entries in subsequent if calls
				m_rWeight[index] = 1.0;

			}

		///	done
			return true;
		}

		virtual bool
		begin_layout_extraction(const Layout* pLayout)
		{ return true; }

	///	writes values from a buffer into the interface
		virtual bool extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(MatDistributeDiag_extract, "algebra parallelization");
		//	block type of associated matrix
			typedef typename TMatrix::value_type block_type;

			block_type in_ii;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				//	get index
				const size_t index = interface.get_element(iter);

				//	read incoming diagonal
				Deserialize(buff, in_ii);

				// set matrix entry
			//	std::cerr << "Recv'd(" << index << ")=" << in_ii << std::endl;
				m_rMat(index, index) += in_ii;
			}

		///	done
			return true;
		}

	private:
	//	pointer to current vector
		TMatrix& m_rMat;
		TVector& m_rWeight;
		double m_dTheta;


};

/// Communication Policy to copy slave couplings to master row
/**
 * This Communication policy assumes that the matrix is stored additively with
 * no overlap. It then copies the row entries of the slaves indices to the
 * corresponding master rows. Only those values are copied that have if the index
 * the slave couples to has also a representation on the other process. This means
 * that couplings to inner indices or to other slave/master, that have only
 * connections to other processes, are not taken into account.
 *
 * \tparam TMatrix	matrix type
 */
template <class TMatrix>
class ComPol_MatAddRowsOverlap0
	: public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
	///	Constructor setting the vector
	/**
	 * vGlobalID must have size >= mat.num_rows()
	 */
		ComPol_MatAddRowsOverlap0(TMatrix& rMat, AlgebraIDVec& vGlobalID)
			: m_rMat(rMat), m_vGlobalID(vGlobalID)
		{
			UG_ASSERT(vGlobalID.size() >= m_rMat.num_rows(), "too few Global ids");
		}

	///	writes the interface values into a buffer that will be sent
		virtual bool collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_MatAddRowsOverlap0_collect, "algebra parallelization");
			typedef typename TMatrix::row_iterator row_iterator;
			typedef typename TMatrix::value_type block_type;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	count number of row entries
				const row_iterator rowEnd = m_rMat.end_row(index);
				size_t numRowEntry = 0;
				for(row_iterator it_k = m_rMat.begin_row(index); it_k != rowEnd; ++it_k)
					numRowEntry++;

			//	write number of row entries to stream
				Serialize(buff, numRowEntry);

			//	write entries and global id to stream
				for(row_iterator it_k = m_rMat.begin_row(index); it_k != rowEnd; ++it_k)
				{
					const size_t k = it_k.index();
					block_type& a_ik = it_k.value();

				//	write global entry to buffer
					Serialize(buff, m_vGlobalID[k]);

				//	write entry into buffer
					Serialize(buff, a_ik);
				}
			}

		///	done
			return true;
		}

		virtual bool
		begin_layout_extraction(const Layout* pLayout)
		{
		//	fill the map global->local
			GenerateAlgebraIDHashList(m_algIDHash, m_vGlobalID);
			return true;
		}

	///	writes values from a buffer into the interface
		virtual bool extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_MatAddRowsOverlap0_extract, "algebra parallelization");
		//	block type of associated matrix
			typedef typename TMatrix::value_type block_type;

		//	we'll read global ids into this variable
			AlgebraID gID;

		//	we'll read blocks into this var
			block_type block;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	read the number of connections
				size_t numConnections = 0;
				Deserialize(buff, numConnections);

			//	read each connection
				for(size_t i_conn = 0; i_conn < numConnections; ++i_conn){
					Deserialize(buff, gID);
					Deserialize(buff, block);

				//	if gID exists on this process, then add the connection to
				//	the matrix.
					size_t conInd;
					if(m_algIDHash.get_entry(conInd, gID)){
					//	add connection between index and conInd to matrix
						m_rMat(index, conInd) += block;
					}
				}
			}

		///	done
			return true;
		}

	private:
	//	pointer to current vector
		TMatrix& m_rMat;

	//	map localID->globalID
		AlgebraIDVec& m_vGlobalID;

	//	map globalID->localID
		AlgebraIDHashList	m_algIDHash;

};


/// Communication Policy to copy couplings between interfaces
/**
 * This Communication policy copies the row entries of the source interface
 * corresponding master interface. Connections are only copied if they connect indices
 * which lie in the same interface. This means
 * that couplings to inner indices or to other slave/master interfaces,
 * are not taken into account.
 *
 * \tparam TMatrix	matrix type
 */
template <class TMatrix>
class ComPol_MatCopyRowsOverlap0
	: public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
	///	Constructor setting the vector
	/**
	 * vGlobalID must have size >= mat.num_rows()
	 */
		ComPol_MatCopyRowsOverlap0(TMatrix& rMat, AlgebraIDVec& vGlobalID)
			: m_rMat(rMat), m_vGlobalID(vGlobalID)
		{
			UG_ASSERT(vGlobalID.size() >= m_rMat.num_rows(), "too few Global ids");
		}

	///	writes the interface values into a buffer that will be sent
		virtual bool collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_MatAddRowsOverlap0_collect, "algebra parallelization");
			typedef typename TMatrix::row_iterator row_iterator;
			typedef typename TMatrix::value_type block_type;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	count number of row entries
				const row_iterator rowEnd = m_rMat.end_row(index);
				size_t numRowEntry = 0;
				for(row_iterator it_k = m_rMat.begin_row(index); it_k != rowEnd; ++it_k)
					numRowEntry++;

			//	write number of row entries to stream
				Serialize(buff, numRowEntry);

			//	write entries and global id to stream
				for(row_iterator it_k = m_rMat.begin_row(index); it_k != rowEnd; ++it_k)
				{
					const size_t k = it_k.index();
					block_type& a_ik = it_k.value();

				//	write global entry to buffer
					Serialize(buff, m_vGlobalID[k]);

				//	write entry into buffer
					Serialize(buff, a_ik);
				}
			}

		///	done
			return true;
		}

		virtual bool
		begin_layout_extraction(const Layout* pLayout)
		{
		//	fill the map global->local
			GenerateAlgebraIDHashList(m_algIDHash, m_vGlobalID);
			return true;
		}

	///	writes values from a buffer into the interface
		virtual bool extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_MatAddRowsOverlap0_extract, "algebra parallelization");
		//	block type of associated matrix
			typedef typename TMatrix::value_type block_type;

		//	we'll read global ids into this variable
			AlgebraID gID;

		//	we'll read blocks into this var
			block_type block;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	read the number of connections
				size_t numConnections = 0;
				Deserialize(buff, numConnections);

			//	read each connection
				for(size_t i_conn = 0; i_conn < numConnections; ++i_conn){
					Deserialize(buff, gID);
					Deserialize(buff, block);

				//	if gID exists on this process, then set the connection to
				//	the new value.
					size_t conInd;
					if(m_algIDHash.get_entry(conInd, gID)){
						m_rMat(index, conInd) = block;
					}
				}
			}

		///	done
			return true;
		}

	private:
	//	pointer to current vector
		TMatrix& m_rMat;

	//	map localID->globalID
		AlgebraIDVec& m_vGlobalID;

	//	map globalID->localID
		AlgebraIDHashList	m_algIDHash;

};


/// Comm-Pol to add matrix rows of inner-interface couplings
/**
 * \tparam TMatrix	matrix type
 */
template <class TMatrix>
class ComPol_MatAddSetZeroInnerInterfaceCouplings
	: public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
	///	Constructor setting the matrix
		ComPol_MatAddSetZeroInnerInterfaceCouplings(TMatrix& rMat)
			: m_rMat(rMat)
		{}

	///	writes the interface values into a buffer that will be sent
		virtual bool collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_MatAddInnerInterfaceCouplings_collect, "algebra parallelization");
			typedef typename TMatrix::row_iterator row_iterator;
			typedef typename TMatrix::value_type block_type;

		//	build map
			Hash<size_t, size_t> hash((size_t)(interface.size() * 1.1));
			size_t interfaceIndex = 0;

			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter, ++interfaceIndex){
				hash.insert(interface.get_element(iter), interfaceIndex);
			}

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	count number of row entries
				const row_iterator rowEnd = m_rMat.end_row(index);
				size_t numRowEntry = 0;
				for(row_iterator it_k = m_rMat.begin_row(index); it_k != rowEnd; ++it_k){
					if(!hash.has_entry(it_k.index())) continue;
					numRowEntry++;
				}

			//	write number of row entries to stream
				Serialize(buff, numRowEntry);

			//	write entries and global id to stream
				for(row_iterator it_k = m_rMat.begin_row(index); it_k != rowEnd; ++it_k)
				{
					const size_t target = it_k.index();
					if(!hash.get_entry(interfaceIndex, target)) continue;
					block_type& a_ik = it_k.value();

				//	write global entry to buffer
					Serialize(buff, interfaceIndex);

				//	write entry into buffer
					Serialize(buff, a_ik);

				//	set entry to zero
					a_ik *= 0;
				}
			}

			return true;
		}

	///	writes values from a buffer into the interface
		virtual bool extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_MatAddRowsOverlap0_extract, "algebra parallelization");
		//	block type of associated matrix
			typedef typename TMatrix::value_type block_type;

		//	build map
			std::vector<size_t> vMapInterfaceToGlob(interface.size());
			int k = 0;
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter, ++k){
				vMapInterfaceToGlob[k] = interface.get_element(iter);
			}

		//	we'll read blocks into this var
			block_type block;
			size_t locID;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	read the number of connections
				size_t numConnections = 0;
				Deserialize(buff, numConnections);

			//	read each connection
				for(size_t i_conn = 0; i_conn < numConnections; ++i_conn){
					Deserialize(buff, locID);
					Deserialize(buff, block);

					UG_ASSERT(locID < vMapInterfaceToGlob.size(), "Invalid"
					          " id: "<<locID<<", size: "<<vMapInterfaceToGlob.size());
					m_rMat(index, vMapInterfaceToGlob[locID]) += block;
				}
			}

		///	done
			return true;
		}

	private:
	//	pointer to current vector
		TMatrix& m_rMat;
};



/// @}

}//	end of namespace


#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__COMMUNICATION_POLICIES__ */
