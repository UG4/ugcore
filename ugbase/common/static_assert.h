//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m9 d11

#ifndef __H__STATIC_ASSERT__
#define __H__STATIC_ASSERT__

//	STATIC_ASSERT is only active during debug-mode.
#ifndef NDEBUG

template <bool> struct CompileTimeAssertion;
template <> struct CompileTimeAssertion<true>
{
	void do_assert()	{};
};

////////////////////////////////////////////////////////////////////////
//	STATIC_ASSERT
///	Checks an expression at compile-time and raises a compile-error if the expression equals 0.
/**
 * STATIC_ASSERT is only active if NDEBUG is not defined.
 *
 * \param expr: an arbitrary expression that can be converted to bool.
 * \param msg: this message will be part of the error-message that the
 * 				compiler prints to the screen.
 * 				Only intended to help the user to identify the error-source.
 * 				The format of the message has to adhere to the rules
 * 				that apply for class-names or variable-names.
 *
 * (bad) example:	STATIC_ASSERT(sizeof(int) == 4, size_of_int_has_to_be_4);
 */
#define STATIC_ASSERT(expr, msg) \
	{CompileTimeAssertion<expr ? true : false> STATIC_ASSERT_ERROR_##msg;\
	STATIC_ASSERT_ERROR_##msg.do_assert();}

#else
//	define an empty STATIC_ASSERT.
#define STATIC_ASSERT(expr, msg)

#endif

#endif
