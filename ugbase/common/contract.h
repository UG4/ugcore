/**
 * \file contract.h
 * \addtogroup ugbase_common
 * \{
 *
 * \brief A Design by Contract(TM) system for C++
 *
 * This file provides a basic Design by Contract(TM) system as described in
 * <i>Meyer, B., Applying �Design by Contract�, IEEE Computer 25(10), 40-51, October 1992.</i>
 * It is inspired by the article <i>Welch,D./Strong, S., An Exception-Based Assertion
 * Mechanism for C++, J. of Object-Oriented Programming, July/August 1998, pp. 50-60</i>
 * and by <i>Isernhagen, R., Softwaretechnik in C und C++, Carl Hanser Verlag M�nchen Wien, pp. 694-706</i>
 *
 * A contract for a class constist of a set of pre- and post-condition as well as a
 * set of invariants.
 * <i>Pre-conditions</i> are the conditions which must be met for the correct
 * function of a method.
 * The method in turn guarantees certain <i>post-conditions</i> after its correct use.
 * <i>Invariants</i> are conditions which have to be <code>true</code> during the lifetime
 * of an object.
 * To check post-conditions it is sometimes necessary to have a copy of the object in its
 * state before the method was called. The contract system supports this with the statement
 * <code>uses_old(Type)</code>. It creates a copy of the object named <code>old</code>.
 * <code>old</code> variables should only be used in the contract statements.
 * The checking of the contract can be controlled by the <code>CONTRACT_LEVEL</code> macro.
 * A value of <code>0</code> means no contract checks, a value of <code>1</code> the checking
 * of all pre-conditions, a value of <code>2</code> the checking of post-conditions and the
 * instanciation of the <code>old</code> variables. A value of <code>3</code> enables the
 * checking of the invariants. Below you can find a simple example of a contract-checked
 * stack which consumes integer values.
 *
 * \code
 *
 * #define CONTRACT_LEVEL 2
 * #include "contract.h"
 *
 * class Stack
 * {
 * private:
 * 	int* _data;
 * 	std::size_t _n;
 * 	std::size_t _size;
 *
 * 	Stack& operator= (const Stack& rhs);
 * 	begin_invariants
 * 		invariant(count() >= 0 && count() <= _size);
 * 		invariant(_size >= 0);
 * 		invariant(implies(_size > 0, _data != NULL));
 * 	end_invariants
 *
 * public:
 *
 * 	Stack(std::size_t size = 10) : _n(0), _size(size)
 * 	{
 * 		pre(size > 0);
 * 		_data = new int[size];
 * 		post(count() == _n);
 * 		check_invariants();
 * 	}
 *
 * 	Stack(const Stack& orig) : _n(orig._n), _size(orig._size)
 * 	{
 * 		_data = new int[orig._size];
 * 		for(std::size_t i = 0; i < orig._n; i++)
 * 		{
 * 			_data[i] = orig._data[i];
 * 		}
 * 		post(_n == orig._n && _size == orig._size);
 * 	}
 *
 * 	~Stack()
 * 	{
 * 		delete [] _data;
 * 	}
 *
 * 	const std::size_t& count() const
 * 	{
 * 		return _n;
 * 	}
 *
 * 	bool is_empty() const
 * 	{
 * 		return _n == 0;
 * 	}
 *
 * 	bool is_full() const
 * 	{
 * 		return _n == _size;
 * 	}
 *
 * 	void push(const int& x)
 * 	{
 * 		std::cout << x << std::endl;
 * 		pre(!is_full());
 * 		uses_old(Stack);
 * 		_data[_n++] = x;
 * 		post(!is_empty());
 * 		post(top() == x);
 * 		post(count() == old.count() + 1);
 * 		check_invariants();
 * 	}
 *
 * 	int pop()
 * 	{
 * 		pre(!is_empty());
 * 		uses_old(Stack);
 * 		int result = _data[--_n];
 * 		post(!is_full());
 * 		post(count() == old.count() - 1);
 * 		check_invariants();
 * 		return result;
 * 	}
 *
 * 	int top() const
 * 	{
 * 		pre(!is_empty());
 * 		uses_old(Stack);
 * 		int result = _data[_n-1];
 * 		post(!is_empty());
 * 		post(count() == old.count());
 * 		check_invariants();
 * 		return result;
 * 	}
 * };
 *
 * \endcode
 *
 * \note Doxygen supplies special tags for contracts, i.e. <code>\\pre</code>, <code>\\post</code>
 * and <code>\\invariant</code>
 *
 * \author Alexander Heusel
 *
 */

#ifndef CONTRACT_H_
#define CONTRACT_H_

#include <cstddef>
#include <stdexcept>
#include <ostream>

#define CHECK_NOTHING 0
#define CHECK_PRE 1
#define CHECK_ENSURE 2
#define CHECK_INVARIANT 3

/**
 * \def CONTRACT_LEVEL
 *
 * \brief This macro defines the level of contract-checking.
 */
#ifndef CONTRACT_LEVEL
	#define CONTRACT_LEVEL CHECK_NOTHING
#endif


#define CONTRACT_ERROR(type, exp, file, line) \
	throw ug::contract_error(#type"("exp")", ug::contract_error::type, file, line)

/**
 * \def pre(exp)
 *
 * \brief The pre-condition macro.
 */
#if(CONTRACT_LEVEL >= CHECK_PRE)
	#define pre(exp)												\
		if(!(exp))													\
		{															\
			CONTRACT_ERROR(pre, #exp, __FILE__, __LINE__);			\
		}
#else
	#define pre(exp)
#endif

/**
 * \def post(exp)
 *
 * \brief The post-condition macro.
 */
#if(CONTRACT_LEVEL >= CHECK_POST)
	#define post(exp)												\
		if(!(exp))													\
		{															\
			CONTRACT_ERROR(post, #exp, __FILE__, __LINE__);			\
		}
#else
	#define post(exp)
#endif

/**
 * \def implies(exp1, exp2)
 *
 * \brief A macro to check implications.
 */
#define implies(exp1, exp2) (!(exp1) || (exp2))

/**
 * \def uses_old(type)
 *
 * \brief A macro which generates a copy of the current object.
 */
#if(CONTRACT_LEVEL >= CHECK_POST)
	#define uses_old(type) type old(*this)
#else
	#define uses_old(type)
#endif

/**
 * \def begin_invariants
 *
 * \brief Start marker for the list of invariants
 */
/**
 * \def end_invariants
 *
 * \brief End marker for the list of invariants
 */
/**
 * \def invariant(exp)
 *
 * \brief A macro to check an invariant.
 */

#if(CONTRACT_LEVEL >= CHECK_INVARIANT)
	#define begin_invariants protected: virtual void _check_invariants(const char* file, std::size_t line) const {
	#define end_invariants }
	#define invariant(exp)											\
		if(!(exp))													\
		{															\
			CONTRACT_ERROR(invariant, #exp, file, line);			\
		}
	#define check_invariants() _check_invariants(__FILE__, __LINE__)
#else
	#define begin_invariants
	#define end_invariants
	#define invariant(exp)
	#define check_invariants()
#endif

// end group ugbase_common
/// \}

namespace ug
{

/// \addtogroup ugbase_common
/// \{

	/**
	 * \class contract_error
	 *
	 * \brief The error class which is thrown if any condition
	 * of a contract is violated.
	 *
	 */
	class contract_error : public std::runtime_error
	{
	public:

		/**
		 * \enum Type
		 *
		 * The type of the breach of contract.
		 */
		enum Type
		{
			pre, /**< Constant for the failure of a pre-condition*/
			post, /**< Constant for the failure of a post-condition*/
			invariant /**< Constant for the failure of an invariant*/
		};

	private:
		std::string _file;
		std::size_t _line;
		Type _type;

	public:
		explicit contract_error(const std::string& __what,
								Type __type,
								const std::string& __file,
								std::size_t __line)
		: runtime_error(__what), _file(__file), _line(__line), _type(__type)
		{}

		virtual ~contract_error() throw()
		{}

		/**
		 * Returns a string containing the path to the source-file
		 * in which the exception occurred.
		 *
		 * \return The source-file in which the exception occurred.
		 */
		const char* file() const
		{
			return _file.c_str();
		}

		/**
		 * Returns the line-number in which the exception occurred.
		 *
		 * \return The line-number in which the exception occurred.
		 */
		std::size_t line() const
		{
			return _line;
		}

		/**
		 * Returns the type of the breach of contract.
		 *
		 * \return The type of the breach of contract.
		 */
		Type type() const
		{
			return _type;
		}

		/**
		 * Writes the exception content to the given <code>std::ostream</code>.
		 *
		 * \param s The <code>std::ostream</code> to which the exception is dumped.
		 * \param e The exception to be dumped.
		 * \return A reference to the passed <code>std::ostream</code>.
		 */
		friend std::ostream& operator<< (std::ostream& s, const contract_error& e)
		{
			return s << e.what() << " in file " << e._file << ", in line " << e._line;
		}

	};

// end group ugbase_common
/// \}

}

#endif /* CONTRACT_H_ */
