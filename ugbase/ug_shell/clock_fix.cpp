/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

/*
 this fixes a linker problem on hermit (and possibly on other clusters):

 /opt/gcc/4.8.2/snos/lib/gcc/x86_64-suse-linux/4.8.2/../../../../lib64/libgomp.a(time.o): In function `gomp_ialias_omp_get_wtime':
time.c:(.text+0xd): undefined reference to `clock_gettime'
time.c:(.text+0x3e): undefined reference to `clock_gettime'
/opt/gcc/4.8.2/snos/lib/gcc/x86_64-suse-linux/4.8.2/../../../../lib64/libgomp.a(time.o): In function `gomp_ialias_omp_get_wtick':
time.c:(.text+0x5d): undefined reference to `clock_getres'
time.c:(.text+0x8e): undefined reference to `clock_getres'

it has something to do with lib rt  http://ubuntuforums.org/showthread.php?t=1870586
so -lrt should solve the problem, but it does NOT, as is written here
http://stackoverflow.com/questions/17150075/undefined-reference-to-clock-gettime-although-lrt-is-given
i couldn't get -Wl,--as-needed to work, so i just added a reference to clock_gettime and clock_getres here.
 */


#include <time.h>
#include <unistd.h>
#include <stdio.h>

bool clock_fix_dummy()
{
	struct timespec tps;
	struct timespec res;
	return clock_gettime(CLOCK_REALTIME, &tps) && clock_getres(CLOCK_MONOTONIC, &res);

}
