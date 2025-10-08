import re
from ase.atoms import Atoms
def read_xyz_path(file_path, index):
    # 读取多个原子结构的信息存储在images列表，通过指定index返回相应原子结构
    # This function reads first all atoms and then yields based on the index.
    # Perfomance could be improved, but this serves as a simple reference.
    # It'd require more code to estimate the total number of images
    # without reading through the whole file (note: the number of atoms
    # can differ for every image).
    with open(file_path, 'r') as fp:
        lines = fp.readlines()
        images = []
        while len(lines) > 0:
            symbols = []
            positions = []
            natoms = int(lines.pop(0))
            lines.pop(0)  # Comment line; ignored
            for _ in range(natoms):
                line = lines.pop(0)
                symbol, x, y, z = line.split()[:4]
                symbol = symbol.lower().capitalize()
                symbols.append(symbol)
                positions.append([float(x), float(y), float(z)])
            images.append(Atoms(symbols=symbols, positions=positions))
        # Atoms类型是可迭代对象，将对象的值逐个返回给调用方
        yield from images[index]
from descriptor.cnPkg import calc_CN_g
def struc_reader(file_path):
    ls_coord = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        try:
            natoms = int(lines.pop(0))
        except:
            print("Invalid XYZ file !!!")
    for index in range(natoms):
        struc = read_xyz_path(file_path, index)
        res = calc_CN_g(struc, natoms, index)
        print(next(res))
        if index==30:
            break
    #     each_line = xyz_raw[i]
    #     x,y,z = map(eval, pat.findall(each_line))
    #     ls.append([x,y,z])
    # return [ls, n]
# struc_reader("../bulk_xyz/AgAl.xyz")
# 读取团簇的xyz文件 file_path=xyz文件路径
# return: 团簇坐标ls，团簇原子数n
def clusterReader(file_path):
    with open(file_path, 'r') as f:
        lines = f.read()
    raw = lines.split('\n')
    n = eval(raw[0])
    xyz_raw = lines.split('\n')[2:]
    ls = []
    pat = re.compile(r"[-+]?\d*\.\d+")
    for i in range(n):
        each_line = xyz_raw[i]
        x,y,z = map(eval, pat.findall(each_line))
        ls.append([x,y,z])
    return [ls, n]
def xyz_reader(file_path):
    # 读取团簇的xyz文件
    return clusterReader(file_path)


def string2numbers(str):
    ls = list(map(eval, str.split()))
    return ls
def poscar_reader(file_path):
    ls = dataReader(file_path)
    basis_vector = []
    for i in range(2, 5):
         basis_vector.append(string2numbers(ls[i]))
    num = eval(ls[5])
    direct = ls[7:]
    pos = []
    for r in direct:
        pos.append(string2numbers(r))
    return [pos,num,basis_vector]
    
# 死循环 读文件内容到文件尾，如实返回文件内容
def dataReader(file_path):
    ls = []
    with open(file_path, 'r', encoding="utf-8") as f:
        while True:
            line = f.readline()
            if line == '':
                break
            ls.append(line)
    return ls
# Pre:
# Can:读取文件的倒数n行内容
# Tip:
def tail(file_path, n, block=-1024):
    with open(file_path, 'rb') as fp:
        # seek(offset, whence)
        # offset=移动位置，whence=参考位置
        # 2代表末尾
        fp.seek(0, 2)
        # tell()返回指针在文件的位置
        file_size = fp.tell()
        while True:
            if file_size >= abs(block):
                fp.seek(block, 2)
                # readlines读取起始位置应该是文件指针
                s = fp.readlines()
                if len(s) > n:
                    return s[-n:]
                else:
                    block *= 2
            else:
                block = -file_size
        
    
# Pre:有广义配位数文件在CN_path
# Can:读取CN_path,返回内容字典
def readCN(CN_path):
    db_CN = {}
    lines=dataReader(CN_path)
    for line in lines:
        ls=line.split()
        key= int(ls[0])
        db_CN[key]=list(map(eval, ls[1:]))
    return db_CN

# 在readCN的基础上，读取多条团簇CN记录
def readCN_pro(CN_path):
    pat_cluster_num=re.compile(r"(?<=>)[0-9]*")
    db_CN={}
    lines=dataReader(CN_path)
    for i in range(len(lines)):
        line=lines[i].strip()
        if line=='':
            continue
        match=pat_cluster_num.search(line)
        if match:
            cur=eval(match.group(0))
            if cur not in db_CN.keys():
                db_CN[cur] = []
            else:
                print("重复团簇")
        else:
            vestaID, val_CN=list(map(eval, line.split()))
            db_CN[cur].append([vestaID, val_CN])
    return db_CN
# 在readCN_pro的基础上，读取标准结构化的数据集合信息(led的标准
def read_std_struct_info(file_path):
    pat_cluster_num=re.compile(r"(?<=>)[0-9]*")
    dt_res={}
    lines=dataReader(file_path)
    for i in range(len(lines)):
        line=lines[i].strip()
        if line=='':
            continue
        match=pat_cluster_num.search(line)
        if match:
            cur=eval(match.group(0))
            if cur not in dt_res.keys():
                dt_res[cur] = {}
            else:
                print("Duplicate KEY")
        else:
            # 此处异常捕获是应对非数字的字符
            # 问题可能出在：输入数据不合适vestaID有多条
            try:
                vestaID, val_CN=list(map(eval, line.split()))
                dt_res[cur][vestaID] = val_CN
            except:
                pass
    return dt_res
# 在readCN_pro的基础上，读取多条团簇allInOne数据的记录
def read_cluster_info(cluster_info_path):
    pat_cluster_num=re.compile(r"(?<=>)[0-9]*")
    db_cluster_info={}
    lines=dataReader(cluster_info_path)
    for i in range(len(lines)):
        line=lines[i].strip()
        if line=='':
            continue
        match=pat_cluster_num.search(line)
        if match:
            cur=eval(match.group(0))
            if cur not in db_cluster_info.keys():
                db_cluster_info[cur] = []
            else:
                print("重复团簇")
        else:
            each_atom = list(map(eval, line.split()))
            db_cluster_info[cur].append(each_atom)
    return db_cluster_info

# Pre:有相对曲率角文件在CA_r_path
# Can:读取相对曲率角文件内容，返回内容字典
def readCA_r(CA_r_path):
    lines=dataReader(CA_r_path)
    db_CA_r={}
    for line in lines:
        ls=line.split()
        # key保留int类型吧
        key = int(ls[0])
        db_CA_r[key]=list(map(eval, ls[1:]))
    return db_CA_r
# 从CN、CA_r的字典中，提取指定位点(依据vestaID)的数据
# dt_sites_vestaID存储所有团簇的期望位点, keys_wanted代表想处理的团簇
# dt_desc是描述符数据的字典，key代表团簇，value代表该团簇的描述符数据列表
def seek_dt_data(dt_desc, dt_sites_vestaID, keys_wanted):
    res = {}
    for num in keys_wanted:
        ls=[]
        for vestaID in dt_sites_vestaID[num]:
            codeID = vestaID-1
            ls.append(dt_desc[num][codeID])
        res[num]=ls
    return res
# Pre:有近邻表文件在neighbor_path (vasp的OUTCAR中明确出现nearest neighbor table)
# Can:读取近邻表文件内容，返回内容字典
def readNeighbor(neighbor_path):
    db_neighbor = {}
    lines = dataReader(neighbor_path)
    for line in lines:
        ls = line.split()
        key = ls[0]
        db_neighbor[key]=list(map(eval, ls[1:]))
    return db_neighbor