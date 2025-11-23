/*
 * Copyright (c) 2010-2020:  G-CSC, Goethe University Frankfurt
 * Author: Tim Schön
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


#ifndef __H__UG__LIB_DISC__TIME_DISC__FINISHED_CONDITIONS
#define __H__UG__LIB_DISC__TIME_DISC__FINISHED_CONDITIONS


#include <vector>
#include "common/util/smart_pointer.h"
#include "common/types.h"

namespace ug {
	

class IFinishedCondition
{
    public:
        virtual bool check_finished(number time, int step) { return false; }
		virtual ~IFinishedCondition() = default;
};

class FinishedTester
{
	public:
		using cond_type = SmartPtr<IFinishedCondition>;
		FinishedTester()= default;

		bool is_finished(number time, int step)
		{
			for(size_t i = 0; i < m_conditions.size(); i++)
			{
				if(m_conditions[i]->check_finished(time, step))
				{
					return true;
				}
			}
			return false;
		}

		void add_condition(cond_type condition)
		{
			m_conditions.push_back(condition);
		}

	
	private:
		std::vector<cond_type> m_conditions;
};

class MaxStepsFinishedCondition : public IFinishedCondition
{	
	public:
		MaxStepsFinishedCondition(int max_timesteps) : m_max_timesteps(max_timesteps)
		{}

		bool check_finished(number time, int step) override {
			return step >= m_max_timesteps;
		}
	private:
		int m_max_timesteps;
};

class TemporalFinishedCondition : public IFinishedCondition
{
	public:
		TemporalFinishedCondition(number end_time, number max_step_size, number relative_precision_bound) : 
			m_end_time(end_time), m_max_step_size(max_step_size), m_relative_precision_bound(relative_precision_bound)
		{}

		bool check_finished(number time, int step) override {
			return (time >= m_end_time) || ((m_end_time-time)/m_max_step_size <= m_relative_precision_bound);
		}

		void set_max_step_size(number max_step_size)
		{
			m_max_step_size = max_step_size;
		}

	private:
		number m_end_time;
		number m_max_step_size;
		number m_relative_precision_bound;
};

}
#endif
