import numpy as np
BOND = 3.3
cn_max = 12
def seek_dis(struc, i, j):
    s = 0
    for k in range(3):
        s += (struc[i][k]-struc[j][k])**2
    yield s**0.5
def seek_neigh(struc, i, natoms):
    for j in range(natoms):
        for dis in seek_dis(struc, i, j):
            if dis < BOND:
                yield j
def calc_CN(struc, natoms, i):
    gen_neigh=seek_neigh(struc, natoms, i)
    s = 0
    for j in gen_neigh:
        gen_dis=seek_dis(struc, i, j)
        if next(gen_dis)<=BOND:
            s += 1
    yield s

def calc_CN_g(struc, natoms, i):
    gen_neigh = seek_neigh(struc, natoms, i)
    s = 0
    for j in gen_neigh:
        gen_CN = calc_CN(struc, natoms, j)
        s += next(gen_CN)
    yield s
    


def seekDis(cluster, n):
    dis = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i):
            s = 0
            for k in range(3):
                s += (cluster[i][k]-cluster[j][k])**2
            dis[i][j] = s**0.5
            dis[j][i] = dis[i][j]
    return dis
def seekNeighbor(dis, n):
    neighbor=[[] for _ in range(n)]
    global BOND
    for i in range(n):
        for j in range(i):
            if dis[i][j] <= BOND:
                neighbor[i].append(j)
                neighbor[j].append(i)
    return neighbor
def calc_cn(n, dis):
    cn_max=12
    global BOND
    neighbor = seekNeighbor(dis, n)
    cn = [0]*n
    CN = [0]*n
    for i in range(n):
        for j in range(i):
            if dis[i][j] <= BOND:
                cn[i] += 1
                cn[j] += 1
    for i in range(n):
        s = 0
        for j in neighbor[i]:
            s += cn[j]
        CN[i] = eval("%.2f" % (s/cn_max))
    # cn为int列表， CN默认保留两位小数浮点数
    return [cn, CN]