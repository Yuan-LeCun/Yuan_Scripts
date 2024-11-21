import numpy as np
import time
import multiprocessing as mp
'''
ratioCN: 分子中配位数情况, 大多数时间都不做考虑
ratioFSI: FSI周围的K, TMP, FSI的数量 (以阴离子为中心)
ratioAGG: 金属阳离子周围的FSI, TMP, FSI的数量 (以金属阳离子为中心)

单溶剂模板 K+ FSI- TMP

ITEM: ATOMS id mol type x y z 
2827 214 7 12.3826 10.3825 13.5417 
id 2827: 原子的ID, 是一个唯一的标识符,用于区分模拟中的不同原子
mol 214: 分子的ID, 表示该原子所属的分子, 用于区分不同的分子或分子片段
type  7: 子类型, 是一个数字，用于区分不同类型的原子, 例如不同的元素或者同一元素的不同状态

mol ID是否是按照packmol的顺序来的? K[80] , FSI[80], 那么K的mol ID: 1-80, FSI的mol ID: 81-160
'''

ratioFSI = {}
total, cycle = 0, 1
start = time.time()

# 计算两个坐标之间的距离 
def calculate_distance(pos1, pos2, boxsize):
    delta = pos1 - pos2
    delta = np.where(delta > boxsize * 0.7, delta - boxsize, delta)
    delta = np.where(delta < -boxsize * 0.7, delta + boxsize, delta)
    return np.linalg.norm(delta, axis=1)

def read_file(file_dir, atoms):
    with open(file_dir, 'r') as f:
        info = f.readlines()
        pos_anion, mol_anion = [], []
        pos_FSI, mol_FSI, type_FSI = [], [], []
        pos_Solv, mol_Solv, type_Solv = [], [], []

        for i in range(atoms):
            line = info[i].split()
            mol_id, mol_type = int(line[1]), int(line[2])
            x, y, z = map(float, line[3:6])
            if mol_type == 1:
                pos_anion.append([x, y, z])
                mol_anion.append(mol_id)
            elif mol_type in {2, 3, 5}:
                pos_FSI.append([x, y, z])
                mol_FSI.append(mol_id)
                type_FSI.append(mol_type)
                
            elif mol_type in {7, 8, 9}:
                pos_Solv.append([x, y, z])
                mol_Solv.append(mol_id)
                type_Solv.append(mol_type)
    return (np.array(pos_anion), np.array(mol_anion), np.array(pos_FSI), np.array(mol_FSI), np.array(type_FSI), \
             np.array(pos_Solv), np.array(mol_Solv), np.array(type_Solv))

def process_file_FSI(file_dir, atoms, boxsize):
    posK, molK, posFSI, molFSI, typeFSI, posTMP, molTMP, typeTMP = read_file(file_dir, atoms)
    
    """
    num = 0
    for i, mol in enumerate(molFSI):
        print(f'{i}: {mol}')
        num += 1
    # num = 80 i=542
    """
    # 每次处理数据时 需要指定FSI mol的范围
    for a in range(81, 161):     # 这里是顺序遍历FSI mol ,从81开始一直到160
        # 在每一个dump or atom.csv中 molFSI都是数组的形式, 里面存储了每一个FSI的mol ID
        mask_a = molFSI == a
        print(molFSI)
        print(len(molFSI))   # 80*7=560

        # print(mask_a[0])

        # np.any(mask_a) 判断mask_a中是否有True, 如果没有, 则跳过这个循环
        if not np.any(mask_a):
            continue

        molFSIK, molFSITMP, molFSIFSI = set(), set(), set()
        countFSI = 0
        # np.where(condition) 返回的是一个tuple, tuple[0]是一个数组, 里面存储了condition为True的索引
        for b in np.where(mask_a)[0]:
            # print(f'{b} {mask_a[b]}')

            countFSI += 1
            if countFSI > 7:
                break
            # print(posFSI[b]) # True-FSI的坐标
            # 利用mask_a 和 np.where(condition)[0]就不需要每次都取遍历了, 直接算距离就行了, 因为大部分都是False但np.where把False给去掉了
            # 数组的广播机制, posK - posFSI[b] 返回的是一个数组, 里面存储了每一个K与FSI[b]的坐标差
            FSI_K_distances = calculate_distance(posK, posFSI[b], boxsize)
            FSI_K_valid = ((typeFSI[b] == 2) & (FSI_K_distances <= 3.972)) | \
                    ((typeFSI[b] == 3) & (FSI_K_distances <= 3.348)) | \
                    ((typeFSI[b] == 5) & (FSI_K_distances <= 3.468))
            molFSIK.update(molK[FSI_K_valid])

            FSI_FSI_dist = calculate_distance(posFSI, posFSI[b], boxsize)
            FSI_FSI_valid = molFSI[(FSI_FSI_dist <= 4) & (molFSI != a)]
            molFSIFSI.update(FSI_FSI_valid)

            FSI_TMP_dist = calculate_distance(posTMP, posFSI[b], boxsize)
            FSI_TMP_valid = ((typeTMP == 7) & (FSI_TMP_dist <= 3.6)) | \
                             ((typeTMP == 8) & (FSI_TMP_dist <= 4)) | \
                             ((typeTMP == 9) & (FSI_TMP_dist <= 4))
            molFSITMP.update(molTMP[FSI_TMP_valid])

        key = f'K={len(molFSIK)},TMP={len(molFSITMP)},FSI={len(molFSIFSI)}'
        ratioFSI[key] = ratioFSI.get(key, 0) + 1

    return ratioFSI

# 元素和范围需要确定! 还需要修改!
def process_file_Solv(file_dir, atoms, boxsize):
    posK, molK, posFSI, molFSI, typeFSI, posTMP, molTMP, typeTMP = read_file(file_dir, atoms)

    for b in range(161, 845):
        mask_b  = molTMP == b
        if not np.any(mask_b):
            continue
        molFSIK, molFSITMP, molFSIFSI = set(), set(), set()
        countTMP = 0

        for c in np.where(mask_b)[0]:
            countTMP += 1
            # 限制TMP的数量 <= 7
            if countTMP > 7:
                break

            TMP_K_distances = calculate_distance(posK, posTMP[c], boxsize)
            TMP_K_valid = ((typeTMP[c] == 7) & (TMP_K_distances <= 3.6)) | \
                          ((typeTMP[c] == 8) & (TMP_K_distances <= 4)) | \
                          ((typeTMP[c] == 9) & (TMP_K_distances <= 4))
            molFSIK.update(molK[TMP_K_valid])

            TMP_FSI_dist = calculate_distance(posFSI, posTMP[c], boxsize)
            TMP_FSI_valid = ((typeFSI == 2) & (TMP_FSI_dist <= 3.972)) | \
                            ((typeFSI == 3) & (TMP_FSI_dist <= 3.348)) | \
                            ((typeFSI == 5) & (TMP_FSI_dist <= 3.468))
            molFSIFSI.update(molFSI[TMP_FSI_valid])

            TMP_TMP_dist = calculate_distance(posTMP, posTMP[c], boxsize)
            TMP_TMP_valid = ((typeTMP == 7) & (TMP_TMP_dist <= 3.6)) | \
                            ((typeTMP == 8) & (TMP_TMP_dist <= 4)) | \
                            ((typeTMP == 9) & (TMP_TMP_dist <= 4))
            molFSITMP.update(molTMP[TMP_TMP_valid])


with open('file_list.txt', 'r') as f:
    start = time.time()
    csv_files = [file.strip() for file in f.readlines()]
    for file in csv_files:
        ratioFSI = process_file_FSI(file, 12428, 51.83015687378338)
        print(f'{cycle} completed')
        cycle += 1
    stop = time.time()
    print(ratioFSI)
    print(f'Execution time: {stop - start}')
     # 写入文件
    # ratioFSI从大到小排序
    ratioFSI = dict(ratioFSI)
    ratioFSI = dict(sorted(ratioFSI.items(), key=lambda x: x[1], reverse=True))
    with open('ratioFSI2.txt', 'w') as f:
        for key, value in ratioFSI.items():
            f.write(f'{key}: {value}\n')



  

