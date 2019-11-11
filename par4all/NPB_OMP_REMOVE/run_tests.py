import argparse
import os
import subprocess
import distutils
import time
from distutils import dir_util


def make_batch_file(test_name, test_dir):
    print('\n@@ Making batch script')
    batch_file_path = os.path.join(test_dir, 'batch_job.sh')
    batch_file = open(batch_file_path,'w')
    batch_file.write(
        '#!/bin/bash\n'
        + '#SBATCH --exclusive\n'
        + '$@\n'
    ) 
    batch_file.close()
    return batch_file_path


def run_test(test, dest_dir, start=100, end=100, stride=100, times=1, sbatch=False, slurm_partition='nxtscl'):
    test_name = os.path.splitext(os.path.basename(test))[0]
    if sbatch:
        batch_file = make_batch_file(test_name, dest_dir)     
        for n in range(start, end+1, stride):
            for i in range(0,times):
                job_name = test_name+'_'+str(n)+'_'+str(i)
                print('\n@ running job: '+job_name)
                sub_proc = subprocess.Popen(['sbatch', '-p', slurm_partition, '-o', job_name, '-J', job_name,
                                             batch_file, test, '-n', str(n), '-j', job_name], cwd=dest_dir)
                sub_proc.wait()
    else:
        for n in range(start, end+1, stride):
            for i in range(0,times):
                job_name = test_name+'_'+str(n)+'_'+str(i)
                print('\n@ running job: '+job_name)
                sub_proc = subprocess.Popen([test, '-n', str(n), '-j', job_name, '|', 'tee',job_name], cwd=dest_dir)
                sub_proc.wait()


def main(dir, start=100, end=200, stride=100, times=1, sbatch=False, slurm_partition='nxtscl'):
    print('\n@@ Got dir: '+dir)
    print('\n@@ Got times: '+str(times))

    dir = os.path.abspath(dir)
    dir_name = os.path.basename(dir)
    timestr = time.strftime("%y%m%d%H%M%S")
    dest_dir = os.path.abspath(os.path.join('results',dir_name+'_results_'+timestr))
    os.makedirs(dest_dir, exist_ok=True)
    print('\n@@ destination directory: '+dest_dir)

    for root, dirs, files in os.walk(dir):
        for name in files:
            split_ext = os.path.splitext(name)
            if split_ext[1] == '.x':
                test = os.path.join(root,name)
                print("\n@@ runing test: "+test)
                test_dest_dir = os.path.join(dest_dir, split_ext[0])
                os.makedirs(test_dest_dir, exist_ok=True)
                run_test(test, test_dest_dir, start, end, stride, times, sbatch, slurm_partition)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Runs tests with varying input sizes.')
    parser.add_argument('-dir',
        dest='dir',
        help='Path to the directory containing the tests.')
    parser.add_argument('-start',
        nargs='?',
        const=1,
        default=1,
        type=int,
        dest='start',
        help='Maximum input size.')
    parser.add_argument('-end',
        nargs='?',
        const=1,
        default=1,
        type=int,
        dest='end',
        help='Maximum input size.')
    parser.add_argument('-stride',
        nargs='?',
        const=1,
        default=1,
        type=int,
        dest='stride',
        help='Maximum input size.')
    parser.add_argument('-t',
        nargs='?',
        const=1,
        default=1,
        type=int,
        dest='times',
        help='Number of times each test will be executed.')
    parser.add_argument('-sbatch',
        action='store_true',
        dest='sbatch',
        help='Run with sbatch.')
    parser.add_argument('-p',
        nargs='?',
        const=1,
        default='nxtscl',
        type=str,
        dest='slurm_partition',
        help='Slurm partition to use with sbatch.')
    args = parser.parse_args()

    main(args.dir,
         start=args.start,
         end=args.end,
         stride=args.stride,
         times=args.times,
         sbatch=args.sbatch,
         slurm_partition=args.slurm_partition)
    print("\n@@ Done")
