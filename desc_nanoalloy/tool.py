coordinate_out_fmt = "    {:18.12f}    {:18.12f}    {:18.12f}"
def sFloatPoint(num,point=4,width=''):
    if point != -1:
        fmt_str = "{{:<{}.{}f}}".format(width, point)
        string = fmt_str.format(num)
    else:
        string=str(num)
    return string
def averageBondLen(dis, n):
    BOND = 3.3
    lenRes = 0
    cnt = 0
    for i in range(n):
        for j in range(i):
            if dis[i][j] <= BOND:
                lenRes += dis[i][j]
                cnt += 1
    return lenRes/cnt
def get_cluster_diameter(dis):
    # dis是二维数组
    dia = max(x for x in max(dis))
    return dia
def cbm_bond_prop(dis, n):
    BOND=3.3
    bond_min = 999
    cnt = 0
    for i in range(n):
        for j in range(i):
            if dis[i][j]<=BOND:
                cnt+=1
                if bond_min > dis[i][j]:
                    bond_min = dis[i][j]
    return [bond_min,cnt]
def writeList(file_handle, List, point, sep='\n'):
    for i in range(len(List)):
        file_handle.write(sFloatPoint(List[i], point)+sep)
def writeMatrix(file_handle, Matrix, point, pre=''):
    for i in range(len(Matrix)):
        List = Matrix[i]
        file_handle.write(pre)
        for j in range(len(List)):
            file_handle.write(sFloatPoint(List[j],point)+' ')
        file_handle.write('\n')
def writeTwoList(file_handle, ls1, ls2, sep1='\n', \
                 sep2='\n', p1=0, p2=3):
    if len(ls1) != len(ls2):
        print("Two lists must have same len")
        return -1
    for i in range(len(ls1)):
        file_handle.write(sFloatPoint(ls1[i],p1)+sep1)
        file_handle.write(sFloatPoint(ls2[i],p2)+sep2)
    return 0

def genDataFile(file_handle, dt, order_of_access=None):
    with open(file_handle, 'w') as f:
        if order_of_access == None:
            ls = dt.keys()
        else:
            ls = order_of_access
        for key in ls:
            print(key,file=f)
            print(dt[key],file=f)
def find_sites(ls_cluster_xyz, n, alpha_shape):
    points=alpha_shape.vertices
    flag=[0]*n
    epsilon=0.001
    ls_alpha2vesta = [0]*len(points)
    for i in range(len(points)):
        x, y, z=[points[i][k] for k in range(3)]
        for j in range(n):
            if flag[j]==1:
                continue
            x1, y1, z1=[ls_cluster_xyz[j][k] for k in range(3)]
            if (abs(x-x1)+abs(y-y1)+abs(z-z1))<=epsilon:
                flag[j]=1
                ls_alpha2vesta[i] = j+1
    sites_codeID=[i for i in range(n) if flag[i]==1]
    return [sites_codeID,ls_alpha2vesta]
def id2xyz(id, ls_xyz):
    return ls_xyz[id]

def calcEad_old(tag, CN, CA_r):
    k, b=0, 0
    if CN >= 4:
        if tag == "O2":
            k=0.19
            b=-1.16
        elif tag == "CO":
            k=0.17
            b=-1.55
        return k*CN+b
    else:
        if tag=="O2":
            k=0.74
            b=-0.56
        elif tag=="CO":
            k=1.25
            b=-1.65
        return k*CA_r+b

def calcEad_old(tag, CN_g, CA_r):
    k, b=0, 0
    if CN_g >= 4:
        if tag == "O2":
            k=0.14
            b=-0.83
        elif tag == "CO":
            k=0.15
            b=-1.47
        return k*CN_g+b
    else:
        if tag=="O2":
            k= 1.58
            b=-1.02
        elif tag=="CO":
            k=1.59
            b=-1.82
        return k*CA_r+b

def calcEad(adsorbate, subtrate, GCN, RCA):
    k, b=0, 0
    seg_flag = False
    dt_k = {
    "Au": {"O2":[1.58,0.14],"CO":[1.59,0.15]},
    "Pt": {"O":[0.118], "O2":[0.227], "OH":[0.166], "OOH":[0.152], "H2O":[0.092], "H2O2":[0.087]},
    "Cu": {"C3H8O":[0.16], "CH3COCH3":[0.27], "CO":[0.0932], "CHO":[0.108]},
    "Co": {"C3H8O":[0.16], "CH3COCH3":[0.27]}
    }
    dt_b = {
    "Au": {"O2":[-1.02,-0.83],"CO":[-1.82,-1.47]},
    "Pt": {"O":[-1.936], "O2":[-2.418], "OH":[-3.811], "OOH":[-2.303], "H2O":[-0.886], "H2O2":[-0.906]},
    "Cu": {"C3H8O":[-1.76], "CH3COCH3":[-1.43], "CO":[-0.5543], "CHO":[0.4624]},
    "Co": {"C3H8O":[-1.93], "CH3COCH3":[-1.98]}
    }
    slope = dt_k[subtrate][adsorbate]
    intercept = dt_b[subtrate][adsorbate]
    if seg_flag == True:
        if GCN >= 4:
            k, b = slope[1], intercept[1]
            Ead = k*GCN + b
        else:
            k, b = slope[0], intercept[0]
            Ead = k*RCA + b
    else:
        k,b = slope[0], intercept[0]
        Ead = k*GCN + b
    return Ead


import numpy as np
def mov_to_pos(arr_cart):
    for col in range(3):
        mi = min(arr_cart[:,col])
        if mi < 0:
            delta = abs(mi)+1
            # np允许对array进行广播操作
            arr_cart[:,col] += delta
    return arr_cart
def get_lat_from_xyz(arr_cart):
    cart_1d=arr_cart.flatten()
    ma = max(cart_1d)
    cell_len=ma*2
    lat_vec=[[cell_len, 0, 0], [0, cell_len, 0], [0, 0, cell_len]]
    return lat_vec
def cart2frac(arr_cart):
    arr_cart = mov_to_pos(arr_cart)
    lat_vec=get_lat_from_xyz(arr_cart)
    conv_mat = np.linalg.inv(np.array(lat_vec))
    mat_frac=[]
    for i in range(len(arr_cart)):
        mat_frac.append((conv_mat*arr_cart[i].T).T)
    ls_frac = []
    for mat in mat_frac:
        ls = []
        for i in range(3):
            ls.append(mat[i][i])
        ls_frac.append(ls)
    return [ls_frac,lat_vec]
spaces = lambda num: " "*num
def xyz2poscar(poscar_file, lat_vec, ls_coordinate,num, sys_name="sys",scale_factor=1.0,symbol=None):
    with open(poscar_file, "w", encoding="utf-8") as fp:
        fp.write(sys_name+'\n')
        fp.write(str(scale_factor)+'\n')
        for i in range(3):
            fp.write(coordinate_out_fmt.format(*lat_vec[i])+'\n')
        fp.write(spaces(4)+symbol+'\n')
        fp.write(spaces(4)+str(num)+'\n')
        fp.write("Direct"+'\n')
        for i in range(num):
            fp.write(coordinate_out_fmt.format(*ls_coordinate[i])+'\n')

import re
def sort_path(raw_paths):
    pat=re.compile(r"[0-9]+")
    ls = []
    for path in raw_paths:
        tag = int(pat.search(path).group())
        ls.append([tag,path])
    ls.sort(key=lambda x:x[0])
    sorted_paths = [el[1] for el in ls]
    return sorted_paths
def get_atom_num_from_vesta(vesta_file):
    pat_start=re.compile(r"THERI")
    pat_atom = re.compile(r"\d+[\s\t]+[a-zA-Z0-9]{2,}[\s\t]+[0-9\.]+")
    flag = False
    with open(vesta_file, 'r', encoding='utf-8') as fp:
        cnt=0
        for line in fp:
            if re.search(pat_start, line)!=None:
                flag=True
                continue
            if flag:
                if re.search(pat_atom, line):
                    cnt += 1
                else:
                    break
    return cnt

def get_path_dict(dir):
    xyz_and_vesta = ""
    dt = {}
    pat=re.compile(r"[0-9]+")
    for path in dir:
        if re.search(r"vesta$", path) != None:
            idx=1
        elif re.search(r"xyz$", path)!= None:
            idx=0
        size=int(pat.search(path).group())
        if size not in dt.keys():
            dt[size] = ['','']
        dt[size][idx] = path
    dt={size: value for size, value in dt.items()
        if value[0]!='' and value[1]!=''}
        
    return dt

def get_size_from_path(path, extension):
    pat=re.compile(r"[0-9]+")
    if re.search(r"{}$".format(extension), path)!=None:
        size=int(pat.search(path).group())
        return size
    else:
        print("Find nothing")

# 用在调试中
def buyer(op, select, cur_tag):
    if op==1 and cur_tag!=select:
        return True
    else:
        return False

import glob
import os
def get_xyz_path(path):
    pat = re.compile(r"(?<=[\/\\]).*(?=\.xyz)")
    ls = glob.glob(os.path.join(path,"*xyz"))
    ls.sort()
    ls_xyz_tag = []
    for xyz_file_path in ls:
        xyz_tag = pat.search(xyz_file_path).group(0)
        ls_xyz_tag.append(xyz_tag)
    return [ls,ls_xyz_tag]

def get_fypte_path(path, ftype):
    pat_name = re.compile(r"(?<=[\/\\])[^\/\\]+(?=\."+re.escape(ftype)+"$)")
    ls_file = glob.glob(os.path.join(path, "*{}".format(ftype)))
    ls_file.sort()
    ls_tag = []
    for file_path in ls_file:
        match = pat_name.search(file_path)
        file_tag = match.group(0)
        if match:
            ls_tag.append(file_tag)
    return [ls_file, ls_tag]


def classify_number(num, ls_interval):
    for i in range(len(ls_interval)-1):
        if ls_interval[i]<=num<ls_interval[i+1]:
            return ls_interval[i]
    print("上限不合适")
def get_type_of_site(val_CN):
    ls_interval=[0, 4, 6, 8]
    lower_bound = classify_number(val_CN, ls_interval)
    type_geo = ls_interval.index(lower_bound)
    return type_geo



    
