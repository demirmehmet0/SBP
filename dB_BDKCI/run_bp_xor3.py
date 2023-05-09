import argparse
import os
import glob
from subprocess import Popen
import _thread as thread

parser = argparse.ArgumentParser(description='Run Boyar-Peralta with XOR2 / XOR3 / XOR4')
parser.add_argument('-xor4c', type=float, required=False, help='Cost of XOR4 gates')
parser.add_argument('-xor3c', type=float, required=False, help='Cost of XOR3 gates')
parser.add_argument('-iwsec', action='store_true', help='Run IWSEC algorithm for XOR3 gates')
parser.add_argument('-xor2c', type=float, required=True, help='Cost of XOR2 gates')
parser.add_argument('-iterations', type=int, required=True, help='Number of iterations')
parser.add_argument('-path', type=str, required=False, help='Path for matrix file', default='.')
parser.add_argument('-matrix', type=str, required=True, help='Name of matrix file')

args = parser.parse_args()

files = glob.glob(args.path + "/*")

for file in files:
    if args.matrix.lower() == file.split("/")[-1].lower():
        args.matrix = file.split("/")[-1]

if args.iwsec:
    if args.xor3c is not None:
        print("IWSEC algorithm does not support dynamic XOR3 costs")
        exit (1)
    elif args.xor4c is not None:
        print("IWSEC algorithm does not support XOR4 gates")
        exit (1)

    rnbp =  'g++ -O2 -o bp_xor3_iwsec_{}_x2-{} -std=c++11 -DTIME_LIMIT={} -DXOR2C={} -DIWSEC main_globalopt_rowcol.cpp'.format(args.iterations, args.xor2c, args.iterations, args.xor2c)
    rnbp_exe = 'bp_xor3_iwsec_{}_x2-{}'.format(args.iterations, args.xor2c)
    
    rnbp_less_than =  'g++ -O2 -o bp_xor3_iwsec_{}_x2-{}_less_than -std=c++11 -DTIME_LIMIT={} -DXOR2C={} -DIWSEC -DLESS_THAN main_globalopt_rowcol.cpp'.format(args.iterations, args.xor2c, args.iterations, args.xor2c)
    rnbp_less_than_exe = 'bp_xor3_iwsec_{}_x2-{}_less_than'.format(args.iterations, args.xor2c)

elif args.xor4c is not None:
    if args.xor3c is not None:
        if args.xor2c==0 or args.xor3c==0 or args.xor4c==0:
            print ("XOR-2 / XOR-3 / XOR-4 gate costs cannot be 0")
            exit (1)
        rnbp =  'g++ -O2 -o bp_xor4_{}_x4-{}_x3-{}_x2-{} -std=c++11 -DTIME_LIMIT={} -DXOR4C={} -DXOR3C={} -DXOR2C={} main_globalopt_rowcol.cpp'.format(args.iterations, args.xor4c, args.xor3c, args.xor2c, args.iterations, args.xor4c, args.xor3c, args.xor2c)
        rnbp_exe = 'bp_xor4_{}_x4-{}_x3-{}_x2-{}'.format(args.iterations, args.xor4c, args.xor3c, args.xor2c)
        
        rnbp_less_than =  'g++ -O2 -o bp_xor4_{}_x4-{}_x3-{}_x2-{}_less_than -std=c++11 -DTIME_LIMIT={} -DXOR4C={} -DXOR3C={} -DXOR2C={} -DLESS_THAN main_globalopt_rowcol.cpp'.format(args.iterations, args.xor4c, args.xor3c, args.xor2c, args.iterations, args.xor4c, args.xor3c, args.xor2c)
        rnbp_less_than_exe = 'bp_xor4_{}_x4-{}_x3-{}_x2-{}_less_than'.format(args.iterations, args.xor4c, args.xor3c, args.xor2c)
    else:
        print ("Please pass XOR-3 gate cost")
        exit (1)

elif args.xor3c is not None:
    if args.xor2c==0 or args.xor3c==0:
        print ("XOR-2 / XOR-3 gate costs cannot be 0")
        exit (1)
    if args.xor3c > 2 * args.xor2c:
        print("WARN: XOR3 cost is greater than 2 times XOR2 cost. Solution may not be optimal")
    rnbp =  'g++ -O2 -o bp_xor3_{}_x3-{}_x2-{} -std=c++11 -DTIME_LIMIT={} -DXOR3C={} -DXOR2C={} main_globalopt_rowcol.cpp'.format(args.iterations, args.xor3c, args.xor2c, args.iterations, args.xor3c, args.xor2c)
    rnbp_exe = 'bp_xor3_{}_x3-{}_x2-{}'.format(args.iterations, args.xor3c, args.xor2c)
    
    rnbp_less_than =  'g++ -O2 -o bp_xor3_{}_x3-{}_x2-{}_less_than -std=c++11 -DTIME_LIMIT={} -DXOR3C={} -DXOR2C={} -DLESS_THAN main_globalopt_rowcol.cpp'.format(args.iterations, args.xor3c, args.xor2c, args.iterations, args.xor3c, args.xor2c)
    rnbp_less_than_exe = 'bp_xor3_{}_x3-{}_x2-{}_less_than'.format(args.iterations, args.xor3c, args.xor2c)

else:
    if args.xor2c==0:
        print ("XOR-2 gate cost cannot be 0")
        exit (1)
    rnbp =  'g++ -O2 -o bp_xor2_{}_x2-{} -std=c++11 -DTIME_LIMIT={} -DXOR2C={} main_globalopt_rowcol.cpp'.format(args.iterations, args.xor2c, args.iterations, args.xor2c)
    rnbp_exe = 'bp_xor2_{}_x2-{}'.format(args.iterations, args.xor2c)
    
    rnbp_less_than =  'g++ -O2 -o bp_xor2_{}_x2-{}_less_than -std=c++11 -DTIME_LIMIT={} -DXOR2C={} -DLESS_THAN main_globalopt_rowcol.cpp'.format(args.iterations, args.xor2c, args.iterations, args.xor2c)
    rnbp_less_than_exe = 'bp_xor2_{}_x2-{}_less_than'.format(args.iterations, args.xor2c)


out = os.system(rnbp)
if out != 0:
    exit(1)
out = os.system(rnbp_less_than)
if out != 0:
    exit(1)

commands = [
    './' + rnbp_exe + ' ' + args.matrix + ' < ' + args.path + "/" + args.matrix,
    './' + rnbp_less_than_exe + ' ' + args.matrix + ' < ' + args.path + "/" + args.matrix,
]

processes = [Popen(cmd, shell=True) for cmd in commands]

for p in processes: p.wait()
