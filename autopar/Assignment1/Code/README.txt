Instructions:
	To serial and parallel run of the benchmarks, inside "NPB_OMP_REMOVE" folder open terminal and run the following command: "python3 script.py"
	
* "script.py" uses "run_test.py" and "autoPar.py" scripts to compile the benchmarks and run them in serial with sbatch.
	Then, it creates folder name "parallel_benchamrks" and all the parallel benchmarks will be there.
	It makes them parallel, compiles them and executes them with sbatch.
	
* all the results will be at NPB_OMP_REMOVE/results.