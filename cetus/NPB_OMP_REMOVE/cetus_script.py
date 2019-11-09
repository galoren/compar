import argparse
import os
import subprocess
#done test 'lu', 'mg'
#missing test ua

TESTS = ['bt', 'lu', 'mg' ,'cg', 'ep', 'sp']

#---------------------------------------------------------------- 
def creating_cetus_directory(main_dir):
    cetus_dir = main_dir+"/Cetus_output"
    os.makedirs(cetus_dir, exist_ok=True)#creating cetus directory
    subprocess.run(['cp -r common sys config {}' .format(cetus_dir)], shell=True, cwd=main_dir)
    subprocess.run(['cp Makefile {}' .format(cetus_dir)], shell=True, cwd=main_dir)
    os.makedirs(cetus_dir+"/bin", exist_ok=True)#creating bin directory
    
#----------------------------------------------------------------          
    
#----------------------------------------------------------------     
def cetus_compiler(main_dir,cetus_dir,test):
    print("\n==========================================  Cetus Compiling :"+test.upper()+" ==========================================\n")
    test_dir='{}/{}'.format(main_dir, test.upper())
    subprocess.run(['cp -r {} {}'.format(test_dir, cetus_dir)], shell=True, cwd=main_dir)#coping all test to new one
    cetus_test_dir='{}/{}'.format(cetus_dir, test.upper())
    
    subprocess.run(['cp -r ../common/*.h .' ], shell=True, cwd=cetus_test_dir)

    for root, dirs, files in os.walk(cetus_test_dir):
        for file in files:
          if os.path.splitext(file)[1] == '.c': 
              
              subprocess.run([' cetus -alias=3 {}  '.format(file)], shell=True, cwd=cetus_test_dir)
             
            
    subprocess.run(['mv cetus_output/* {}'.format(cetus_test_dir)], shell=True, cwd=cetus_test_dir)        
    subprocess.run(['rm -r cetus_output'], shell=True, cwd=cetus_test_dir)
#---------------------------------------------------------------- 
  
#---------------------------------------------------------------- 
def compile_test(des_dir,test,test_class):  
    print("\n========================= Compiling : "+test.upper()+" ======================================\n")
    
    subprocess.run(['make', test.upper(), 'CLASS={}'.format(test_class)], cwd=des_dir)   
#---------------------------------------------------------------- 

#---------------------------------------------------------------- 
def main(main_dir, test, test_class,parallel,grid):
    
    main_dir = os.path.abspath(main_dir)
    des_dir=main_dir
    
    if parallel:
        creating_cetus_directory(main_dir)
        cetus_dir = main_dir+"/Cetus_output"
        des_dir=cetus_dir
        
        
    if test == 'suite':
      subprocess.run(['rm bin/*'], shell=True, cwd=des_dir)#clean bin directory
      for test in TESTS:
        if parallel:
          cetus_compiler(main_dir,cetus_dir,test)
          
        compile_test(des_dir,test,test_class)  
                
    else:
      if parallel:
            cetus_compiler(main_dir,cetus_dir,test)
      subprocess.run(['rm bin/*'], shell=True, cwd=des_dir)#clean bin directory 
      compile_test(des_dir,test,test_class)
            
 
    print("========================================== Starting tests ==========================================")
    if grid :
      subprocess.run(['python3 run_tests.py -dir {} -p grid -sbatch '.format(des_dir+"/bin"  )], shell=True, cwd=main_dir)
    else:
      subprocess.run(['python3 run_tests.py -dir {} '.format(des_dir+"/bin"  )], shell=True, cwd=main_dir)     
    
    
    print("DONE!")
#---------------------------------------------------------------- 

#---------------------------------------------------------------- 
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
        default='suite',
        type=str,
        dest='test',
        help='test benchmark to run')
    parser.add_argument('-p',
        nargs='?',
        const=1,
        default=False,
        type=bool,
        dest='parallel',
        help='compile the problem in parallel')
    parser.add_argument('-g',
        nargs='?',
        const=1,
        default=False,
        type=bool,
        dest='grid',
        help='send to grid ')
    parser.add_argument('-c',
        nargs='?',
        const=1,
        default='C',
        type=str,
        dest='test_class',
        help='Test class = problem size')

    args = parser.parse_args()

    main(args.main_dir, args.test, args.test_class ,args.parallel,args.grid)
#----------------------------------------------------------------     
    
    
    
    
    

