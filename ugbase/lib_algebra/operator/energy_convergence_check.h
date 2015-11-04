
#ifndef __H__LIB_ALGEBRA__OPERATOR__ENERGY_CONVERGENCE_CHECK__
#define __H__LIB_ALGEBRA__OPERATOR__ENERGY_CONVERGENCE_CHECK__

#include "convergence_check.h"

namespace ug{


template <typename TVector>
class EnergyConvCheck : public StdConvCheck<TVector>
{
	typedef StdConvCheck<TVector> base_type;
	public:
	EnergyConvCheck() : base_type() {}
	EnergyConvCheck(int maxSteps, number minDefect, number relReduction) : base_type(maxSteps, minDefect, relReduction) {}
	EnergyConvCheck(int maxSteps, number minDefect, number relReduction, bool verbose)
	 	 : base_type(maxSteps, minDefect, relReduction, verbose) {}
	EnergyConvCheck(int maxSteps, number minDefect, number relReduction, bool verbose,bool suppressUnsuccessful)
	 	 : base_type(maxSteps, minDefect, relReduction, verbose, suppressUnsuccessful) {}

	virtual SmartPtr<IConvergenceCheck<TVector> > clone()
	{
		SmartPtr<EnergyConvCheck<TVector> > newInst(new EnergyConvCheck<TVector>);
		// use std assignment (implicit member-wise is fine here)
		*newInst = *this;
		return newInst;
	}

	void start(const TVector& d)
	{
		base_type::start_defect(energy_norm(d));
	}
	void update(const TVector& d)
	{
		base_type::update_defect(energy_norm(d));
	}

	double energy_norm(const TVector &d)
	{
		if(tmp.valid() == false || tmp->size() != d.size())
		{
			tmp = d.clone_without_values();
			tmp2 = d.clone_without_values();
		}
		TVector &t = *tmp;
		TVector &t2 = *tmp2;
		t = d;
#ifdef UG_PARALLEL
		t.change_storage_type(PST_CONSISTENT);
#endif
		m_op->apply(t2, t);
		return sqrt(VecProd(t, t2));
	}

	virtual std::string config_string() const
	{
		std::stringstream ss;
		ss << "EnergyConvCheck( max steps = " << base_type::m_maxSteps << ", min defect = " << base_type::m_minDefect <<
				", relative reduction = " << base_type::m_relReduction << ")";
		return ss.str();
	}

	void set_linear_operator(SmartPtr<ILinearOperator<TVector> > op)
	{
		m_op = op;
	}

private:
	SmartPtr<TVector> tmp, tmp2;
	SmartPtr<ILinearOperator<TVector> > m_op;

};


}

#endif /* __H__LIB_ALGEBRA__OPERATOR__ENERGY_CONVERGENCE_CHECK__ */
