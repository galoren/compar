## This directory contains directories:
* Retults - Contains graphs with comparisons between the benchmark tests run time results and the output results of their run
* NPB_OMP_REMOVE - Contains the Cetus compiler script (cetus_script.py), benchmark test (BT,LU,MG,CG,EP,SP,UA) .

# Cetus compiler script
[cetus_script.py](https://github.com/yoelv92/cetus_project/blob/master/NPB_OMP_REMOVE/cetus_script.py)  is an script that the main purpose is to compile sourse to sourse using Cetus to make a parallel progrm and compile with icc Intel compiler and run benchmark test of the NAS problems (BT,LU,MG,CG,EP,SP,UA), as a result before using the script must to load the Cetus compiler and Intel icc compiler.

#### the directives to use Cetus and icc compilers are:
```sh
$ module load intel

$ module load cetus
```
When choosing the parallel option (-p) in the sctipt, the script is opening new directory(if not exist) and copy the chosen benchmark test and all the following directories and files they needed. 
All the parallel benchmark test are in the Cetus_output directory .

### Collaborative directories:
* Cetus_output - Contains copy of benchmark tests with openMP directives(parallel code), the directory is created after using the cetus_script.py on -p option.
* bin - Contains all the files that after the compiling.
* results - Contains all the results thet return from the grid.

### Collaborative files:
#### The files must be in the same directory with cetus_script.py.
* [Makefile](https://github.com/yoelv92/cetus_project/blob/master/NPB_OMP_REMOVE/Makefile) - Use to compile the tests
* [run_test.py](https://github.com/yoelv92/cetus_project/blob/master/NPB_OMP_REMOVE/run_tests.py) - The script ment for runig the tests after compiling .

### Script options: 
| Option | Input | Description | Default |
| ------ | ------ | ------ | ------ |
| -dir | path to the tests directory | Path to the directory containing the benchmark tests directories | Current directory |
| -c | C or W | The size for the input class C or W  | C |
| -t | BT, LU, MG, CG, EP, SP, UA  |Choosing which benchmark test to run (BT,LU,MG,CG,EP,SP,UA)  | Run all |
| -p | Nathing | Compile sourse to sourse using Cetus to make benchmark test run in parallel and compile with Intel icc compiler| Only serialy Intel icc compiler |
| -g | Nathing | the benchmark test to run in the grid  | Not sending to grid |

#### To execute:
```sh
$ python3 cetus_test.py  -dir <path> -c <class> -t <test> 
```

