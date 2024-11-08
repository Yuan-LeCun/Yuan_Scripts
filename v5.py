from time import time
import numpy as np
import multiprocessing as mp

"""
DMM:  Represent Solvent Molecules
FSI:  Negative ions
对于单溶剂来说, 需要修改的地方 已打上Tag
Tag(1) atom_type in {} : 溶剂分子中带负电的原子类型
Tag(2) atom_type in {} : FSI的类型
Tag(3) 中心C原子是sp3杂化, 对于此类sp3杂化分子而言仅考虑O配位情况; 若是sp2杂化, 则考虑带负电且有一个明显的近邻峰的原子配位情况;
       typeDMM[idsDMM[j]] == atom_type(O) and LDMM <= RDF_(K-DMM(O))
Tag(4) 问题同Tag(3) 
       typeFSI[idsFSI[z]] == atom_type(O) and LFSI <= RDF_(K-FSI(O)) or\
         (typeFSI[idsFSI[z]] == atom_type(F) and LFSI <= K-FSI(F)) or \
            (typeFSI[idsFSI[z]] == atom_type(N) and LFSI <= K-FSI(N)
"""

# anion: K+, Li+, Na+
# solvent: E3, N7, H6, etc.
# cation: FSI-



def read_file(file_dir, atoms, box_size):
    with open(file_dir, 'r') as f:

        info = f.readlines() 
        pos_anion, ids_anion, mol_anion = [], [], {1: {}, 2: {}, 3: {}}
        pos_sol_e3, ids_sol_e3, mol_sol_e3, type_sol_e3 = [], [], [], {}
        posFSI, idsFSI, molFSI, typeFSI, AGGFSI = [], [], [], {}, {}
        
        # 读取文件中的原子信息
        for i in range(atoms):
            line = info[i].split()
            atom_id, mol, atom_type = int(line[0]), int(line[1]), int(line[2])
            x, y, z = map(float, line[3:6])

            if atom_type == 1:
                pos_anion.append([x, y, z])
                ids_anion.append(atom_id)
                for key in mol_anion:
                    mol_anion[key].setdefault(str(mol), 0)
            # Tag(1) 
            elif atom_type in {6, 8, 9}:
                pos_sol_e3.append([x, y, z])
                ids_sol_e3.append(atom_id)
                mol_sol_e3.append(mol)
                type_sol_e3[atom_id] = atom_type
            # Tag(2)
            elif atom_type in {2, 3, 5}:
                posFSI.append([x, y, z])
                idsFSI.append(atom_id)
                molFSI.append(mol)
                typeFSI[atom_id] = atom_type
                AGGFSI.setdefault(str(mol), 0)

    return pos_anion, ids_anion, mol_anion, pos_sol_e3, ids_sol_e3, mol_sol_e3, posFSI, idsFSI, molFSI, typeFSI, type_sol_e3, AGGFSI

def process_file(file, atoms, box_size):
    pos_anion, ids_anion, mol_anion, pos_sol_e3, ids_sol_e3, mol_sol_e3, posFSI, idsFSI, molFSI, typeFSI, type_sol_e3, AGGFSI = read_file(file, atoms, box_size)

    # 转换为数组
    pos_anion = np.array(pos_anion)
    pos_sol_e3 = np.array(pos_sol_e3)
    posFSI = np.array(posFSI)

    total = len(ids_anion)
    ratioCN, ratioAGGFSI = {}, {}
    ratioAGG = {'AGG': 0, 'CIP': 0, 'SSIP': 0}

    # 遍历K与DMM和FSI的配位数
    for i in range(len(ids_anion)):
        num_sol_e3, numFSI, numCN_sol_e3, numCNFSI = 0, 0, 0, 0
        mol_sol_e3_2, molFSI2 = {}, {}

        for j in range(len(ids_sol_e3)):
            x_sol_e3, y_sol_e3, z_sol_e3 = pos_anion[i] - pos_sol_e3[j]
            x_sol_e3 = x_sol_e3 - box_size if x_sol_e3 > box_size * 0.7 else x_sol_e3 + box_size if x_sol_e3 <= -box_size * 0.7 else x_sol_e3
            y_sol_e3 = y_sol_e3 - box_size if y_sol_e3 > box_size * 0.7 else y_sol_e3 + box_size if y_sol_e3 <= -box_size * 0.7 else y_sol_e3
            z_sol_e3 = z_sol_e3 - box_size if z_sol_e3 > box_size * 0.7 else z_sol_e3 + box_size if z_sol_e3 <= -box_size * 0.7 else z_sol_e3
            L_sol_e3_dist = np.sqrt(x_sol_e3**2 + y_sol_e3**2 + z_sol_e3**2)

            # Tag(3)
            if (type_sol_e3[ids_sol_e3[j]] == 6 and L_sol_e3_dist <= 4.116):
                if mol_sol_e3[j] not in mol_sol_e3_2:
                    num_sol_e3 += 1
                    numCN_sol_e3 += 1
                    mol_sol_e3_2[mol_sol_e3[j]] = 'exist'
                else:
                    numCN_sol_e3 += 1

        for z in range(len(idsFSI)):
            xFSI, yFSI, zFSI = pos_anion[i] - posFSI[z]
            xFSI = xFSI - box_size if xFSI > box_size * 0.7 else xFSI + box_size if xFSI <= -box_size * 0.7 else xFSI
            yFSI = yFSI - box_size if yFSI > box_size * 0.7 else yFSI + box_size if yFSI <= -box_size * 0.7 else yFSI
            zFSI = zFSI - box_size if zFSI > box_size * 0.7 else zFSI + box_size if zFSI <= -box_size * 0.7 else zFSI

            LFSI = np.sqrt(xFSI**2 + yFSI**2 + zFSI**2)

            # RDF第一个峰后的峰谷  
            # Tag(4)
            if (typeFSI[idsFSI[z]] == 2 and LFSI <= 4.14) or (typeFSI[idsFSI[z]] == 3 and LFSI <= 3.804) or (typeFSI[idsFSI[z]] == 5 and LFSI <= 3.66):
                if molFSI[z] not in molFSI2:
                    numFSI += 1
                    numCNFSI += 1
                    molFSI2[molFSI[z]] = 'exist'
                else:
                    numCNFSI += 1

        key = f'Solv_e3={num_sol_e3},FSI={numFSI},CNDMM={numCN_sol_e3},CNFSI={numCNFSI}'
        ratioCN[key] = ratioCN.get(key, 0) + 1

    # 遍历FSI与K的配位情况
    for d in range(len(idsFSI)):
        for e in range(len(ids_anion)):
            xFSI, yFSI, zFSI = pos_anion[e] - posFSI[d]
            xFSI = xFSI - box_size if xFSI > box_size * 0.7 else xFSI + box_size if xFSI <= -box_size * 0.7 else xFSI
            yFSI = yFSI - box_size if yFSI > box_size * 0.7 else yFSI + box_size if yFSI <= -box_size * 0.7 else yFSI
            zFSI = zFSI - box_size if zFSI > box_size * 0.7 else zFSI + box_size if zFSI <= -box_size * 0.7 else zFSI
            LFSI = np.sqrt(xFSI**2 + yFSI**2 + zFSI**2)

#    pos_anion, ids_anion, mol_anion, pos_sol_e3, ids_sol_e3, mol_sol_e3, posFSI, idsFSI, molFSI, typeFSI, type_sol_e3, AGGFSI = read_file(file, atoms, box_size)

            # 同上Tag(4)
            if (typeFSI[idsFSI[d]] == 2 and LFSI <= 4.14) or (typeFSI[idsFSI[d]] == 3 and LFSI <= 3.804) or (typeFSI[idsFSI[d]] == 5 and LFSI <= 3.66): 
                if mol_anion[1].get(str(ids_anion[e]), 0) == 0:
                    mol_anion[1][str(ids_anion[e])] = str(molFSI[d])
                    AGGFSI[str(molFSI[d])] += 1
                elif mol_anion[1].get(str(ids_anion[e])) != str(molFSI[d]):
                    if mol_anion[2].get(str(ids_anion[e]), 0) == 0:
                        mol_anion[2][str(ids_anion[e])] = str(molFSI[d])
                        AGGFSI[str(molFSI[d])] += 1
                    elif mol_anion[2].get(str(ids_anion[e])) != str(molFSI[d]):
                        mol_anion[3][str(ids_anion[e])] = str(molFSI[d])
                        AGGFSI[str(molFSI[d])] += 1

    # 统计SSIP, CIP, AGG
    for key in mol_anion[1]:
        if mol_anion[1][key] == 0:
            ratioAGG['SSIP'] += 1
        elif mol_anion[2][key] == 0:
            ratioAGG['CIP'] += 1
        else:
            ratioAGG['AGG'] += 1

    for key in AGGFSI:
        agg_key = f'{AGGFSI[key]}K-FSI'
        ratioAGGFSI[agg_key] = ratioAGGFSI.get(agg_key, 0) + 1

    return total, ratioCN, ratioAGG, ratioAGGFSI

# 多进程合并
def collect_results(result, total, ratioCN, ratioAGG, ratioAGGFSI):
    # result: 单个进程结果 
    total_res, ratioCN_res, ratioAGG_res, ratioAGGFSI_res = result
    total.value += total_res
    for key, value in ratioCN_res.items():
        ratioCN[key] = ratioCN.get(key, 0) + value
    for key, value in ratioAGG_res.items():
        ratioAGG[key] = ratioAGG.get(key, 0) + value
    for key, value in ratioAGGFSI_res.items():
        ratioAGGFSI[key] = ratioAGGFSI.get(key, 0) + value

if __name__ == '__main__':
    start = time()
    # 提前定义好参数 file.ipynb会给出
    atoms = 12068
    box_size = 49.99699348576955
    manager = mp.Manager() # 多进程管理器，用于在多进程间共享数据
    # 共享变量 total，赋值0 ，共享ratioCN、ratioAGGFSI、ratioAGG 字典集
    total = manager.Value('i', 0) 
    ratioCN = manager.dict()
    ratioAGGFSI = manager.dict()
    ratioAGG = manager.dict({'AGG': 0, 'CIP': 0, 'SSIP': 0})

    with open('file_list.txt', 'r') as f:
        csv_files = [file.strip() for file in f.readlines()]

    # process pool
    pool = mp.Pool(mp.cpu_count())
    for file in csv_files:
        pool.apply_async(process_file, args=(file, atoms, box_size), callback=lambda res: collect_results(res, total, ratioCN, ratioAGG, ratioAGGFSI))

    pool.close() # close the pool
    pool.join()  # wait all process done

    print('-------end reading files-------')
    print(dict(ratioAGG)) 
    print(dict(ratioAGGFSI)) 
    print(dict(ratioCN)) 

    print('time:', time() - start)
    
    # 将ratioCN 写入文件 如'DMM=3,FSI=0,CNDMM=7,CNFSI=0' 仅写入3 0 7 0
    with open('ratioCN.txt', 'w') as f:
        for key, value in ratioCN.items():
            f.write(f'{key.split(",")[0].split("=")[1]} {key.split(",")[1].split("=")[1]} {key.split(",")[2].split("=")[1]} {key.split(",")[3].split("=")[1]} {value}\n')


