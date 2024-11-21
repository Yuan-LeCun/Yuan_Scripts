import numpy as np
import time
from multiprocessing import Pool, cpu_count

"""
溶剂化壳层中离子键平均寿命!!!
该代码用于跟踪钾离子与周围分子(H2O、DMSO、FSI)之间的连接关系
通过比较当前和之前的连接分子集合，检测连接的变化，并计算连接的持续时间，最终求得平均连接寿命
使用多进程加速计算 temp_v2.py 耗时425s temp_v3.py 耗时67s
简化版的并行计算
"""

def calculate_distance(pos1, pos2, boxsize):
    delta = pos1 - pos2
    delta = np.where(delta > boxsize * 0.7, delta - boxsize, delta)
    delta = np.where(delta < -boxsize * 0.7, delta + boxsize, delta)
    return np.linalg.norm(delta, axis=1)

def read_file(filename, atoms):
    with open(filename, 'r') as f:
        info = f.readlines()
        posK, IdsK = [], []
        posH2O, IdsH2O = [], []
        posDMSO, IdsDMSO = [], []
        posFSI, IdsFSI, molFSI, typeFSI = [], [], [], []

        for i in range(atoms):
            line = info[i].split()
            Id, Mol, Type = int(line[0]), int(line[1]), int(line[2])
            x, y, z = map(float, line[3:6])
            if Type == 7:
                posK.append([x, y, z])
                IdsK.append(Id)
            elif Type == 1:
                posH2O.append([x, y, z])
                IdsH2O.append(Id)
            elif Type == 5:
                posDMSO.append([x, y, z])
                IdsDMSO.append(Id)
            elif Type in {8, 9, 11}:
                posFSI.append([x, y, z])
                IdsFSI.append(Id)
                typeFSI.append(Type)
                molFSI.append(Mol)

    return np.array(posK), IdsK, np.array(posH2O), IdsH2O, np.array(posDMSO), IdsDMSO, np.array(posFSI), IdsFSI, molFSI, np.array(typeFSI)

def process_file(args):
    file_dir, atoms, boxsize = args
    posK, IdsK, posH2O, IdsH2O, posDMSO, IdsDMSO, posFSI, IdsFSI, molFSI, typeFSI = read_file(file_dir, atoms)
    step = int(file_dir.split('.')[-2])
    time_step = step/500     # timestep=2  250steps totaltime:500fs
    
    results = []
    for idx, IdK in enumerate(IdsK):
        pos_K = posK[idx]
        dist_H2O = calculate_distance(posH2O, pos_K, boxsize)
        dist_DMSO = calculate_distance(posDMSO, pos_K, boxsize)
        dist_FSI = calculate_distance(posFSI, pos_K, boxsize)

        connected_H2O = set(np.array(IdsH2O)[dist_H2O <= 3.564])
        connected_DMSO = set(np.array(IdsDMSO)[dist_DMSO <= 3.684])

        # np.where(condition, x, y)  condition为True返回x,否则返回y
        cutoff_FSI = np.where(
            typeFSI == 8, 3.684,
            np.where(
                typeFSI == 9, 3.516,
                np.where(
                    typeFSI == 11, 3.66,
                    0
                )
            )
        )
        connected_FSI = set(np.array(IdsFSI)[dist_FSI <= cutoff_FSI])
        connected_all = connected_H2O | connected_DMSO | connected_FSI
        
        results.append((IdK, time_step, connected_all, connected_H2O, connected_DMSO, connected_FSI))
    return results

start = time.time()

def main():
    boxsize = 52.68883463461055
    atoms = 11965
    with open('file_list2.txt', 'r') as f:
        csv_files = [line.strip() for line in f]

    # 初始化字典和计数器
    K_all, K_H2O, K_DMSO, K_FSI = {}, {}, {}, {}
    time_all, time_H2O, time_DMSO, time_FSI = 0, 0, 0, 0
    N_all, N_H2O, N_DMSO, N_FSI = 0, 0, 0, 0
    
    # 创建进程池
    num_cores = cpu_count()
    pool = Pool(processes=num_cores)
    
    # 准备参数
    args = [(file, atoms, boxsize) for file in csv_files]
    
    # 并行处理文件 逻辑基本不变, 都是把process_file先封装好, 然后用pool.imap并行处理
    # 按理来说, file_results中已经有了所有的连接信息和寿命信息, 但是我懒得再开数组了,直接复制上面内容
    for i, file_results in enumerate(pool.imap(process_file, args)):
        step = int(csv_files[i].split('.')[-2])
        
        if step == 0:
            # 初始化所有K离子的连接信息
            for IdK, time_step, connected_all, connected_H2O, connected_DMSO, connected_FSI in file_results:
                K_all[IdK] = [time_step] + list(connected_all)
                K_H2O[IdK] = [time_step] + list(connected_H2O)
                K_DMSO[IdK] = [time_step] + list(connected_DMSO)
                K_FSI[IdK] = [time_step] + list(connected_FSI)
        else:
            for IdK, time_step, connected_all, connected_H2O, connected_DMSO, connected_FSI in file_results:
                prev_all = set(K_all[IdK][1:]) if len(K_all[IdK]) > 1 else set()
                prev_H2O = set(K_H2O[IdK][1:]) if len(K_H2O[IdK]) > 1 else set()
                prev_DMSO = set(K_DMSO[IdK][1:]) if len(K_DMSO[IdK]) > 1 else set()
                prev_FSI = set(K_FSI[IdK][1:]) if len(K_FSI[IdK]) > 1 else set()

                # 溶剂化壳层是否完全转变?
                if prev_all.isdisjoint(connected_all):   
                    N_all += 1
                    time_all += (time_step - K_all[IdK][0])
                    K_all[IdK] = [time_step] + list(connected_all)
                if prev_H2O != connected_H2O:
                    N_H2O += 1
                    time_H2O += (time_step - K_H2O[IdK][0])
                    K_H2O[IdK] = [time_step] + list(connected_H2O)
                if prev_DMSO != connected_DMSO:
                    N_DMSO += 1
                    time_DMSO += (time_step - K_DMSO[IdK][0])
                    K_DMSO[IdK] = [time_step] + list(connected_DMSO)
                if prev_FSI != connected_FSI:
                    N_FSI += 1
                    time_FSI += (time_step - K_FSI[IdK][0])
                    K_FSI[IdK] = [time_step] + list(connected_FSI)
        
        print(f'{i+1} completed')
    
    pool.close()
    pool.join()

    stop = time.time()
    # 计算平均寿命
    ave_all = time_all / N_all if N_all else 0
    ave_H2O = time_H2O / N_H2O if N_H2O else 0
    ave_DMSO = time_DMSO / N_DMSO if N_DMSO else 0
    ave_FSI = time_FSI / N_FSI if N_FSI else 0

    print(f'life time of all is {ave_all} ps')
    print(f'life time of H2O is {ave_H2O} ps')
    print(f'life time of DMSO is {ave_DMSO} ps')
    print(f'life time of FSI is {ave_FSI} ps')
    print(f'total time: {stop - start} s')

if __name__ == '__main__':
    main()
