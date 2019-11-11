import argparse
import os
import subprocess
import shutil
import re

SERIAL_DIR_NAME = 'serial_outputs'
PARALLEL_DIR_NAME = 'parallel_outputs'

tests = ['bt', 'cg', 'ep', 'lu', 'mg', 'sp', 'ua']
tests.remove('ua') # problem with compilation

TESTS_FILES_DICT = {}

def make_dict_of_lines_to_replace(main_dir):
    tests_files = {}

    # {
    #   "test name capital" :
    #   [   -->   list of all files of the test that need to be changed
    #       (
    #           path from main dir of the file,
    #           [   -->   list of all lines that need to be replaced
    #               (original line, replacement line),
    #           ]
    #       ),
    #   ],
    # }
    t = (
        os.path.join(main_dir, 'BT', 'bt.c'),
        [
            (
                'printf(" Number of available threads: %5d\\n", omp_get_max_threads());',
                '// printf(" Number of available threads: %5d\\n", omp_get_max_threads());'
            ),
        ]
    )
    l = []
    l.append(t)
    tests_files['BT'] = l
    t = (
        os.path.join(main_dir, 'CG', 'cg.c'),
        [
            (
                'printf(" Number of available threads: %5d\\n", omp_get_max_threads());',
                '// printf(" Number of available threads: %5d\\n", omp_get_max_threads());'
            ),
            (
                'num_threads = omp_get_num_threads();',
                'num_threads = 32;'
            ),
            (
                'myid = omp_get_thread_num();',
                'myid = 0;'
            ),
        ]
    )
    l = []
    l.append(t)
    tests_files['CG'] = l
    t = (
        os.path.join(main_dir, 'EP', 'ep.c'),
        [
            (
                'printf("\\n Number of available threads:          %13d\\n", omp_get_max_threads());',
                '// printf("\\n Number of available threads:          %13d\\n", omp_get_max_threads());'
            ),
        ]
    )
    l = []
    l.append(t)
    tests_files['EP'] = l
    t = (
        os.path.join(main_dir, 'LU', 'read_input.c'),
        [
            (
                'printf(" Number of available threads: %5d\\n", omp_get_max_threads());',
                '// printf(" Number of available threads: %5d\\n", omp_get_max_threads());'
            ),
        ]
    )
    l = []
    l.append(t)
    t = (
        os.path.join(main_dir, 'LU', 'ssor.c'),
        [
            (
                'mthreadnum = omp_get_num_threads() - 1;',
                '// mthreadnum = omp_get_num_threads() - 1;'
            ),
            (
                'iam = omp_get_thread_num();',
                '// iam = omp_get_thread_num();'
            ),
        ]
    )
    l.append(t)
    tests_files['LU'] = l
    t = (
        os.path.join(main_dir, 'MG', 'mg.c'),
        [
            (
                'printf(" Number of available threads: %5d\\n", omp_get_max_threads());',
                '// printf(" Number of available threads: %5d\\n", omp_get_max_threads());'
            ),
            (
                'myid = omp_get_thread_num();',
                '// myid = omp_get_thread_num();'
            ),
            (
                'num_threads = omp_get_num_threads();',
                'num_threads = 32;'
            ),
        ]
    )
    l = []
    l.append(t)
    tests_files['MG'] = l
    t = (
        os.path.join(main_dir, 'SP', 'sp.c'),
        [
            (
                'printf(" Number of available threads: %5d\\n", omp_get_max_threads());',
                '// printf(" Number of available threads: %5d\\n", omp_get_max_threads());'
            ),
        ]
    )
    l = []
    l.append(t)
    tests_files['SP'] = l
    return tests_files


def replace_lines_for_p4a(test, replace_back=False):
    test = test.upper()
    files = TESTS_FILES_DICT[test]
    if replace_back:
        files = [(
            os.path.splitext(file_element[0])[0] + '.p4a' + os.path.splitext(file_element[0])[1],
            file_element[1]
        ) 
        for file_element in files]

    for file_element in files:
        print('\n\nopening', file_element[0], '...\n\n')
        with open(file_element[0], 'r+') as f:
            filedata = f.read()
            for line_element in file_element[1]:
                if replace_back:
                    print('replacing', line_element[1], 'with', line_element[0], '...\n\n')
                    filedata = filedata.replace(line_element[1], line_element[0])
                else:
                    print('replacing', line_element[0], 'with', line_element[1], '...\n\n')
                    filedata = filedata.replace(line_element[0], line_element[1])
            f.seek(0)
            f.write(filedata)
            f.truncate()


def remove_bswap_function(main_dir):
    common_files = [
        'c_print_results.c',
        'c_timers.c',
        'print_results.c',
        'randdp.c',
    ]
    files_path = []
    common_path = os.path.join(main_dir, 'common')
    for file_name in common_files:
        file_path = os.path.join(common_path, file_name)
        files_path.append(file_path)
    file_path = os.path.join(main_dir, 'BT', 'error.c')
    files_path.append(file_path)
    file_path = os.path.join(main_dir, 'BT', 'verify.c')
    files_path.append(file_path)
    file_path = os.path.join(main_dir, 'LU', 'read_input.c')
    files_path.append(file_path)
    file_path = os.path.join(main_dir, 'LU', 'domain.c')
    files_path.append(file_path)
    file_path = os.path.join(main_dir, 'LU', 'l2norm.c')
    files_path.append(file_path)
    file_path = os.path.join(main_dir, 'LU', 'error.c')
    files_path.append(file_path)
    file_path = os.path.join(main_dir, 'LU', 'verify.c')
    files_path.append(file_path)
    file_path = os.path.join(main_dir, 'SP', 'set_constants.c')
    files_path.append(file_path)
    file_path = os.path.join(main_dir, 'SP', 'rhs.c')
    files_path.append(file_path)
    file_path = os.path.join(main_dir, 'SP', 'error.c')
    files_path.append(file_path)
    file_path = os.path.join(main_dir, 'SP', 'verify.c')
    files_path.append(file_path)
    for test in tests:
        file_path = os.path.join(main_dir, test.upper(), '{}.c'.format(test.lower()))
        files_path.append(file_path)

    bswap_regex = re.compile(r'static __uint64_t __bswap_64[^\}]*\}', re.DOTALL)
    for path in files_path:
        try:
            with open(path, 'r+') as f:
                file_data = f.read()
                file_data = re.sub(bswap_regex, '', file_data)
                f.seek(0)
                f.write(file_data)
                f.truncate()
        except Exception as e:
            print(path, 'doesn\'t exist')
            print(e)


def remove_p4a_suffix(directory):
    print('\n@@ renaming all files in {} with suffix .p4a.c to .c\n'.format(directory))
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.p4a.c'):
                subprocess.run(['mv -v {} {}'.format(file, file[0:-6]+'.c')], shell=True, cwd=directory)


def run_make_veryclean(main_dir):
    print('\n\n@@@ running make clean\n')
    subprocess.run(['make', 'veryclean'], cwd=main_dir)


def compile_serial_test(main_serial_dir, test, test_class):
    test_dir = os.path.join(main_serial_dir, test.upper())
    common_dir = os.path.join(main_serial_dir, 'common')

    run_make_veryclean(main_serial_dir)
    print('\n\n@@@ making original {} serial test with class {}\n'.format(test, test_class))
    subprocess.run(['make', test, 'CLASS={}'.format(test_class)], cwd=main_serial_dir)


def copy_and_move_back(main_dir, src_dir, dest_dir):
    subprocess.run(['cp -r {}/* {}'.format(src_dir, dest_dir)], shell=True, cwd=main_dir)
    subprocess.run(['mv {} {}'.format(dest_dir, src_dir)], shell=True, cwd=main_dir)


def remove_wtime_sgi64_file(main_dir):
    wtime_sgi64_file = os.path.join(main_dir, 'common', 'wtime_sgi64.c')
    os.remove(wtime_sgi64_file) if os.path.exists(wtime_sgi64_file) else None


def rename_wtime_sgi64_file(main_dir, rename_back=False):
    added_extension = '.bak'
    wtime_sgi64_file = os.path.join(main_dir, 'common', 'wtime_sgi64.c')
    if rename_back:
        current_name = wtime_sgi64_file + added_extension
        os.rename(current_name, wtime_sgi64_file) if os.path.exists(current_name) else None
    else:
        new_name = wtime_sgi64_file + added_extension
        os.rename(wtime_sgi64_file, new_name) if os.path.exists(wtime_sgi64_file) else None


def create_new_dir(main_dir, dir_name):
    new_dir = os.path.join(main_dir, dir_name)
    os.makedirs(new_dir, exist_ok=True)
    return new_dir


def move_results_to_bin(src_dir, dest_dir):
    shutil.rmtree(os.path.join(src_dir, 'bin'), ignore_errors=True)
    try:
        os.rename(dest_dir, os.path.join(src_dir, 'bin'))
    except Exception as e:
        print('Error in renaming bin directory')
        print(e)


def compile_parallel_test(main_parallel_dir, test, test_class):
    test_dir = os.path.join(main_parallel_dir, test.upper())
    common_dir = os.path.join(main_parallel_dir, 'common')

    run_make_veryclean(main_parallel_dir)
    subprocess.run(['make', test, 'CLASS={}'.format(test_class)], cwd=main_parallel_dir)

    replace_lines_for_p4a(test)
    rename_wtime_sgi64_file(main_parallel_dir)

    print('\n\n@@@ running p4a')
    subprocess.run(['p4a -v --no-pointer-aliasing --log -O *.c ../common/*.c -I ../common -I../sys'], shell=True, cwd=test_dir)

    replace_lines_for_p4a(test, replace_back=True)
    rename_wtime_sgi64_file(main_parallel_dir, rename_back=True)

    remove_p4a_suffix(test_dir)
    remove_p4a_suffix(common_dir)

    remove_bswap_function(main_parallel_dir)
    run_make_veryclean(main_parallel_dir)

    print('\n\n@@@ making p4a parallel {} test with class {}\n'.format(test, test_class))
    subprocess.run(['make', test, 'CLASS={}'.format(test_class)], cwd=main_parallel_dir)

def compile_test(main_serial_dir, temp_serial_bin_result, main_parallel_dir, temp_parallel_bin_result, test, test_class):
    compile_serial_test(main_serial_dir, test, test_class)
    exec_serial_file_path = os.path.join(main_serial_dir, 'bin', '{}.C.x'.format(test.lower()))
    try:
        shutil.move(exec_serial_file_path, temp_serial_bin_result)
    except Exception as e:
        print('There is no serial executable file for', test)
        print(e)
    compile_parallel_test(main_parallel_dir, test, test_class)
    exec_parallel_file_path = os.path.join(main_parallel_dir, 'bin', '{}.C.x'.format(test.lower()))
    try:
        shutil.move(exec_parallel_file_path, temp_parallel_bin_result)
    except Exception as e:
        print('There is no parallel executable file for', test)
        print(e)

def main(main_dir, test, test_class):
    global TESTS_FILES_DICT

    # creates all needed directories
    main_dir = os.path.abspath(main_dir)
    parent_dir = os.path.abspath(os.pardir)
    serial_dir = os.path.join(main_dir, SERIAL_DIR_NAME)
    parallel_dir = os.path.join(main_dir, PARALLEL_DIR_NAME)

    print('\n\n@@@ deleting old directories...')
    shutil.rmtree(serial_dir, ignore_errors=True)
    shutil.rmtree(parallel_dir, ignore_errors=True)

    print('\n\n@@@ creating new directories...')
    temp_serial_dir = create_new_dir(parent_dir, SERIAL_DIR_NAME)
    temp_parallel_dir = create_new_dir(parent_dir, PARALLEL_DIR_NAME)

    print('\n\n@@@ copy all files...')
    copy_and_move_back(main_dir, main_dir, temp_serial_dir)
    copy_and_move_back(main_dir, main_dir, temp_parallel_dir)
    shutil.rmtree(os.path.join(parallel_dir, SERIAL_DIR_NAME), ignore_errors=True)  # fix it later with more elegant way

    # delete wtime_sgi64.c file (p4a crashed because of this file)
    # remove_wtime_sgi64_file(serial_dir)
    # remove_wtime_sgi64_file(parallel_dir)

    # create temp results directory
    temp_serial_bin_result = create_new_dir(serial_dir, 'temp_bin')
    temp_parallel_bin_result = create_new_dir(parallel_dir, 'temp_bin')

    TESTS_FILES_DICT = make_dict_of_lines_to_replace(parallel_dir)

    if test == 'all':
        for test in tests:
            compile_test(serial_dir, temp_serial_bin_result, parallel_dir, temp_parallel_bin_result, test, test_class)
    else:
        compile_test(serial_dir, temp_serial_bin_result, parallel_dir, temp_parallel_bin_result, test, test_class)
    
    move_results_to_bin(serial_dir, temp_serial_bin_result)
    move_results_to_bin(parallel_dir, temp_parallel_bin_result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Runs tests with varying input sizes.')
    parser.add_argument('-dir',
        nargs='?',
        const=1,
        default='.',
        type=str,
        dest='main_dir',
        help='Path to the directory containing the tests.')
    parser.add_argument('-t',
        nargs='?',
        const=1,
        default='all',
        type=str,
        dest='test',
        help='test benchmark to rum')
    parser.add_argument('-c',
        nargs='?',
        const=1,
        default='C',
        type=str,
        dest='test_class',
        help='Test class = problem size')

    args = parser.parse_args()
    try:
        main(args.main_dir, args.test, args.test_class)
    except Exception as e:
        print('Exception:', e)
