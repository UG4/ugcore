// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// Somewhen back in 2005-2007

#ifndef __UG__LIB_GRID__SPARSE_MATRIX__
#define __UG__LIB_GRID__SPARSE_MATRIX__

namespace ug{
namespace libgrid_simplealg{

///	This is a very basic sparse matrix structure.
/**	This matrix class should only be used for very simple problems, with a small
 * number of non-zero entries per row.
 * It is only used in some of lib_grids algorithms, to keep lib_grid independend of
 * lib_algebra.
 * If you're working with ug, but not inside lib_grid, please consider using the
 * more advanced algebra structures from lib_algebra.
 */
template <class Number>
class SparseMatrix
{
public:
	SparseMatrix()
		{
			m_iNumRows = 0;
			m_iNumCols = 0;
			m_iDataCount = 0;
		}
	~SparseMatrix()
		{
			clear();
		}

	void clear()
		{
			if(m_iDataCount)
			{
				delete m_pData;
				delete m_pIndex;
			}

			m_iNumRows = 0;
			m_iNumCols = 0;
			m_iDataCount = 0;
		}

	void clear_entries()
		{
			for(int i = 0; i < m_iDataCount; ++i)
			{
				m_pData[i] = 0;
				m_pIndex[i] = -1;
			}
		}
		
	void create(int NumRows, int NumCols, int NumNonZeroPerRow)
		{
			if(m_iDataCount)
			{
				delete[] m_pData;
				delete[] m_pIndex;
			}
			m_iNumRows = NumRows;
			m_iNumCols = NumCols;
			m_iNumNonZero = NumNonZeroPerRow;
			m_iDataCount = m_iNumRows * m_iNumNonZero;

			m_pData = new Number[m_iDataCount];
			m_pIndex = new int[m_iDataCount];

			for(int i = 0; i < m_iDataCount; ++i)
			{
				m_pData[i] = 0;
				m_pIndex[i] = -1;
			}
		}
		
	void set(int RowIndex, int ColIndex, Number value)
		{
			int BaseIndex = RowIndex * m_iNumNonZero;
			for(int i = 0; i < m_iNumNonZero; ++i)
			{
				if((m_pIndex[BaseIndex + i] == -1) || (m_pIndex[BaseIndex + i] == ColIndex))
				{
					m_pIndex[BaseIndex + i] = ColIndex;
					m_pData[BaseIndex + i] = value;
					break;
				}
			}
		}
		
	Number get(int RowIndex, int ColIndex)
		{
			int BaseIndex = RowIndex * m_iNumNonZero;
			for(int i = 0; i < m_iNumNonZero; ++i)
			{
				if(m_pIndex[BaseIndex + i] == ColIndex)
					return m_pData[BaseIndex + i];
			}
			return 0;
		}
		
	void vmult(Number* pVecOut, Number* pVecIn)
		{
			for(int i = 0; i < m_iNumRows; ++i)
			{
				pVecOut[i] = 0;
				int Index = i*m_iNumNonZero;
				for(int j = 0; j < m_iNumNonZero; ++j)
				{
					if(m_pIndex[Index] != -1)
						pVecOut[i] += (m_pData[Index] * pVecIn[m_pIndex[Index]]);
					++Index;
				}
			}
		}
		
	void vmultT(Number* pVecOut, Number* pVecIn)
		{
			int i;
			for(i = 0; i < m_iNumCols; ++i)
				pVecOut[i] = 0;
			for(i = 0; i < m_iNumRows; ++i)
			{
				int Index = i*m_iNumNonZero;
				for(int k = 0; k < m_iNumNonZero; ++k)
				{
					if(m_pIndex[Index] != -1)
						pVecOut[m_pIndex[Index]] += m_pData[Index] * pVecIn[i];
					++Index;
				}
			}
		}

		//void print();
	private:
		Number*	m_pData;
		int*	m_pIndex;

		int		m_iNumRows;
		int		m_iNumCols;
		int		m_iNumNonZero;
		int		m_iDataCount;
};

}//	end of namespace
}//	end of namespace
#endif
