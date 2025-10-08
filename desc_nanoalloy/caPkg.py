import math

import numpy as np

MAX = 0xfffff
BOND = 3.3
def seekDeg(v1, v2):
    x1,y1,z1 = [v1[i] for i in range(3)]
    x2,y2,z2 = [v2[i] for i in range(3)]
    val = (x1*x2+y1*y2+z1*z2)/((x1*x1+y1*y1+z1*z1)**0.5 * (x2*x2+y2*y2+z2*z2)**0.5)
    # Return the arc cosine (measured in radians) of val
    rad = math.acos(val)
    deg = math.degrees(rad)
    return deg

def seekMid(p, q):
    x1,y1,z1 = [p[i] for i in range(3)]
    x2,y2,z2 = [q[i] for i in range(3)]
    res = [(x1+x2)/2, (y1+y2)/2, (z1+z2)/2]
    return res
def dot2Vec(p, q):
    x1,y1,z1 = [p[i] for i in range(3)]
    x2,y2,z2 = [q[i] for i in range(3)]
    v = [x1-x2, y1-y2, z1-z2]
    return v
def calc_vec_mag(p, q):
    x1,y1,z1 = [p[i] for i in range(3)]
    x2,y2,z2 = [q[i] for i in range(3)]
    mag = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5
    return mag
def calc_ca(cluster, n, dis,neighbor, debug=False):
    ca=[0.0]*n
    caR=[0.0]*n
    for i in range(n):
        n_s = neighbor[i]
        vertex = cluster[i]
        if len(n_s) <= 1:
            if debug:
                print("case 1 CA")
            continue
        elif len(n_s) == 2:
            if debug:
                print("case 2 CA")
            first = n_s[0]
            second = n_s[1]
            v1 = dot2Vec(vertex,cluster[first])
            v2 = dot2Vec(vertex,cluster[second])
            ca[i] =seekDeg(v1,v2)
            ca[i] = eval("%.2f"%ca[i])
        else:
            mi = MAX
            s = 0
            cnt = 0
            first = 0
            for idx in n_s:
                if dis[i][idx] < mi:
                    mi = dis[i][idx]
                    first = idx
            v1 = dot2Vec(vertex, cluster[first])
            for j in range(len(n_s)):
                second = n_s[j]
                if second == first:
                    continue
                for k in range(j):
                    third = n_s[k]
                    if third == first:
                        continue
                    if dis[second][third] > BOND:
                        continue
                    cnt += 1
                    mid = seekMid(cluster[second],cluster[third])
                    v2 = dot2Vec(vertex, mid)
                    s += seekDeg(v1,v2)
            #     -led有坑
            if s == 0:
                if debug:
                    print("OH, 000!!!")
                continue
            ca[i] = s/cnt
        ca[i] = eval("%.2f" % ca[i])
        caR[i] = eval("%.2f" % (ca[i]/180))
    # 默认ca caR保留两位小数浮点数
    return [ca, caR]
