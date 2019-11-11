-------------------- Compar stage 1 - par4all compiler --------------------

Running the benchmarks before and after par4all compilation and compare the runtime to the article results.

Steps:

1. Move to Nas Benchmarks folder - NPB_OMP_REMOVE:
    $ cd /path/to/NPB_OMP_REMOVE/

2. Run the following commands in the same order:
    $ module load intel
    $ module load par4all
    $ source $set_p4a_env

3. Compile the files by running the script compile_and_run_all_p4a.py:
    The arguments of the script are:
        -dir    Path to the Nas Benchmarks directory - NPB_OMP_REMOVE (the current directory by default)
        -t      The test name, or 'all' to compile all the tests at the same time ('all' by default)
        -c      The test class, in capital letter ('C' by default)
    
    $ python3 compile_and_run_all_p4a.py -dir "/path/to/NPB_OMP_REMOVE/" -t all -c C

4. Run the tests:
    If the compilation from step 2 ended successfully, you will see in the NPB_OMP_REMOVE directory the following folders: 
        - serial_outputs
        - parallel_outputs
    The executable files will be in each folder in a folder named 'bin'.

    If it does, run the tests by running the script twice, one for serials and second for parallels:
        $ python3 run_tests.py -dir "/path/to/NPB_OMP_REMOVE/serial_results/bin" -sbatch -p grid
        $ python3 run_tests.py -dir "/path/to/NPB_OMP_REMOVE/parallel_results/bin" -sbatch -p grid

5. Check the results:
    When the running will end, the results will save in a folder named 'results'.
    The 'results' folder will be both at the 'serial_outputs' and 'parallel_outputs' folders.
