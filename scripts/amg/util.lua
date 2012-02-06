function PrintParallelProfileNode(name)
	pn = GetProfileNode(name)
	t = pn:get_avg_total_time_ms()/to100 * 100
	tmin = ParallelMin(t)
	tmax = ParallelMax(t)
	printf("%s:\n%.2f %%, min: %.2f %%, max: %.2f %%", name, t, tmin, tmax)
end