import numpy as np
import os

# 遍历当前文件夹下所有的out文件
def get_out_files():
    files = os.listdir()
    out_files = []
    for file in files:
        if file.endswith('.out\uf00d'):
            out_files.append(file)
    return out_files

indexs = []
coordinates = []

# 读取out文件
def read_out_file(file, num):
    with open(file, 'r') as f:
        lines = f.readlines()
        normal_termination = False
        line = lines[-1]

        if 'Normal termination of Gaussian 16' in line:
            normal_termination = True
            print(f'{file} has been read successfully!')
            return
            
        if not normal_termination:
            print(f'{file} did not terminate normally.Please check the file!')

            lines = [line.strip() for line in lines]
            # 遍历每一行，如果包含了'Standard orientation' 则记录下标
            for i, line in enumerate(lines):
                if 'Standard orientation' in line:
                    indexs.append(i)
            for i in range(indexs[-1] + 5, len(lines)):
                if '---------------------' in lines[i]:
                    break
                coordinates.append(lines[i])
            # C:6 O:8 H:1
            with open(f'out_pdb{num}.pdb', 'w') as f:
                for line in coordinates:
                    line = line.split()
                    atom_type = line[1]
                    x, y, z = map(float, line[3:6])
                    if atom_type == '1':
                        f.write(f'HETATM    {line[0]}  H  UNL     1      {x:.3f}  {y:.3f}  {z:.3f}  1.00  0.00           H\n')
                    elif atom_type == '6':
                        f.write(f'HETATM    {line[0]}  C  UNL     1      {x:.3f}  {y:.3f}  {z:.3f}  1.00  0.00           C\n')
                    elif atom_type == '8':
                        f.write(f'HETATM    {line[0]}  O  UNL     1      {x:.3f}  {y:.3f}  {z:.3f}  1.00  0.00           O\n')
        # 清空coordinates和indexs
        coordinates.clear()
        indexs.clear()
            
files  = get_out_files()
num = 0
for file in files:
    read_out_file(file, num)
    num += 1

# 遍历当前文件夹中所有的PDB文件
def get_pdb_files():
    pdb_files = []
    for file in os.listdir():
        if file.endswith('.pdb'):
            pdb_files.append(file)
    return pdb_files

# 逐一读取PDB文件中的原子坐标
def read_pdb_file(file):
    pos = []
    with open(file, 'r') as f:
        pdb_data = f.readlines()
        pdb_data = [line.strip() for line in pdb_data if line.startswith('HETATM')]
    for line in pdb_data:

        line = line.split()
        atom_type = line[2]
        # 根据实际的PDB来修改
        if len(line) == 11:
            x, y, z = map(float, line[5:8])
        elif len(line) == 9:
            x, y, z = map(float, line[5:8])
        elif len(line) == 8:
            x, y, z = map(float, line[4:7])
        else:
            continue

        if len(atom_type) == 1: 
            pos.append([atom_type, x, y, z])
        elif len(atom_type) == 2:
            # atom_type[1] 大写转小写
            pos.append([atom_type[0], x, y, z])
        elif len(atom_type) > 2:
            if atom_type[0].isalpha() and atom_type[1].isdigit():
                pos.append([atom_type[0], x, y, z])
            elif atom_type[0].isalpha() and atom_type[1].isalpha():
                pos.append([atom_type[0] + atom_type[1].lower(), x, y, z])
        else: continue
    return np.array(pos, dtype=object)

# 开一个子文件夹, 存放即将要写入的gjf文件
def make_dir():
    if not os.path.exists('gjf_files'):
        os.mkdir('gjf_files')

# 写入gjf文件       opt freq  主要用于算HOMO LUMO或几何优化  后期还要再做一个求RESP2电荷的  在气相和液相中继续做几何优化和算单点
def write_gjf():
    pdb_files = get_pdb_files()
    
    make_dir()
    for file in pdb_files:
        pos = read_pdb_file(file)
        
        with open('gjf_files/' + file[:-4] + '.gjf', 'w') as f:
            f.write('%nprocshared=14\n')  
            f.write('%mem=4GB\n')
            f.write('%chk=' + file[:-4] + '.chk\n')
            f.write('# opt freq b3lyp/6-31g(d) em=gd3bj  \n\n')
            f.write('Title Card Required\n\n')
            f.write('0 1\n')       # 电荷与自选多重度 自己改写
            
            # 按照注释中的数据模式写入
            for atom in pos:
                f.write(' {:<2} {:>12.6f} {:>12.6f} {:>12.6f}\n'.format(atom[0], float(atom[1]), float(atom[2]), float(atom[3])))
            f.write('\n\n')
            f.write('--link1--\n')
            f.write(f'%oldchk=' + file[:-4] + '.chk\n')
            f.write(f'%chk=' + file[:-4] + 'energy' + '.chk\n')
            f.write('# b3lyp/def2tzvp em=gd3bj geom=allcheck\n\n')

# 生成一个bash文件, 用于批量运行g16
def write_bash():
    with open('gjf_files/'+ 'run_g16.sh', 'w') as f:
        f.write('#!/bin/bash\n')
        for file in get_pdb_files():
            f.write('g16 < ' + file[:-4] + '.gjf |tee ' + file[:-4] + '.out\n')

write_gjf()
write_bash()

