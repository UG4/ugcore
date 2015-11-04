
#ifndef __H__UG__COMMON__TRANSPOSE_H_
#define __H__UG__COMMON__TRANSPOSE_H_


template<typename T>
class TE_TRANSPOSED
{
public:
	typedef typename T::value_type value_type;
	TE_TRANSPOSED(const T&_t) : t(_t) {}
	inline size_t num_rows() const
	{
		return t.num_cols();
	}

	inline size_t num_cols() const
	{
		return t.num_rows();
	}

	const value_type &operator() (size_t r, size_t c) const
	{
		return t(c, r);
	}

	value_type &operator() (size_t r, size_t c)
	{
		return t(c, r);
	}
private:
	const T &t;
};

template<typename T>
inline TE_TRANSPOSED<T> transpose(const T &t)
{
	return TE_TRANSPOSED<T>(t);
}

inline const double &transpose(const double &t)
{
	return t;
}

inline double &transpose(double &t)
{
	return t;
}

#endif /* TRANSPOSE_H_ */
