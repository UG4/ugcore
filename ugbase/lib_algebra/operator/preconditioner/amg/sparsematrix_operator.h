
namespace ug
{

template <typename matrix_type, typename vector_type>
class SparseMatrixOperator : public virtual IMatrixOperator<vector_type, vector_type, matrix_type >
{
	public:
		SparseMatrixOperator() {}
		SparseMatrixOperator(matrix_type *pMatrix) : m_pMatrix(pMatrix) { }
		void setmatrix(matrix_type *pMatrix) { m_pMatrix = pMatrix; }
		bool init(const vector_type &) { return true; }
		bool init() { return true; }
		bool apply(vector_type &c, const vector_type &d)
		{
			MatMult(c, 1.0, *m_pMatrix, d);
			return true;
		}
		bool apply_sub(vector_type &c, const vector_type &d)
		{
			MatMultAdd(c, 1.0, c, -1.0, *m_pMatrix, d);
			return true;
		}

		virtual ~SparseMatrixOperator(){}

	public:
	// 	Access to matrix
		virtual matrix_type& get_matrix() { return *m_pMatrix; };

	private:
		matrix_type *m_pMatrix;
};


}
