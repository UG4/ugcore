
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
