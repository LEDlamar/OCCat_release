from math import sqrt,log10,exp,log
k = 2.316e-4
h = 4.136e-15
T = 298
# J/(mol K) >> eV/K
delta_S1=-2.05e-3
delta_S2=-2.14e-3
# delta_S1=-197.7
# delta_S2=-205.1


def calcNumDen(E, delta_S):
    global k, T
    num = -(E - T*delta_S)
    den = k*T
    # print("-(E - T*delta_S)={}, k*T={}".format(num,den))
    return (num,den)
def calc_k(E, delta_S):
    global k,h,T
    num, den = calcNumDen(E, delta_S)
    # print(k*T/h)
    return exp(num/den)*k*T/h
def calc_K(delta_E,delta_S):
    global k,T
    num, den=calcNumDen(delta_E, delta_S)
    return exp(num/den)
def calc_E(E_ad_CO, E_ad_O2):
    keys = ["co_ad","OOCO","TS_1","TS_2","a3+","a3_","a4"]
    E = {key:0 for key in keys}
    # 图5-16
    E["co_ad"] = 1.51*(E_ad_CO+E_ad_O2)+0.08
    E["OOCO"] = 1.49*E["co_ad"]+0.85
    E["TS_1"] = 0.75*E["co_ad"]-0.12
    E["TS_2"] = 0.71*E["OOCO"]-0.16
    # 5.2.2 
    E["a3_pos"] = E["TS_1"]-E["co_ad"]
    E["a3_neg"] = E["TS_1"]-E["OOCO"]
    E["a4_pos"] = E["TS_2"]-E["OOCO"]

    return E
def solve(a, b, c):
    if b*b - 4*a*c > 0:
        x1,x2 = -b+sqrt(b*b-4*a*c),-b-sqrt(b*b-4*a*c)
        x1 /= (2*a)
        x2 /= (2*a)
        return (x1, x2)
    return (None, None)
def calc_cov(Ead_CO,Ead_O2, debug=False, p=None):
    if p == None:
        p = {"CO": 0.01, "O2": 0.21}
    E=calc_E(Ead_CO, Ead_O2)
    K1=calc_K(Ead_CO, delta_S1)
    K2=calc_K(Ead_O2, delta_S2)
    k3_pos=calc_k(E["a3_pos"], 0)
    k3_neg=calc_k(E["a3_neg"], 0)
    k4=calc_k(E["a4_pos"], 0)
    theta={"CO": 0, "O2": 0, "OOCO": 0, "*": 0}
    a=k3_pos*K1*K2*p["CO"]*p["O2"]/(k3_neg+k4)
    b=1+K1*p["CO"]+K2*p["O2"]
    c=-1
    x1, x2=solve(a, b, c)
    if x1>0 and x1<1:
        theta["*"]=x1
    elif x2>0 and x2<1:
        theta["*"]=x2
    elif x1==0 or x2==0:
        theta["*"] = 0
    else:
        if debug:
            print("calc Rate ERROR")
        return theta
    theta["CO"]=p["CO"]*theta["*"]*K1
    theta["O2"]=p["O2"]*theta["*"]*K2
    return theta
# unit as bar
def finalRate(Ead_CO,Ead_O2, debug=False, p=None):
    if p == None:
        p = {"CO": 0.01, "O2": 0.21}
    E=calc_E(Ead_CO, Ead_O2)
    K1=calc_K(Ead_CO, delta_S1)
    K2=calc_K(Ead_O2, delta_S2)
    k3_pos=calc_k(E["a3_pos"], 0)
    k3_neg=calc_k(E["a3_neg"], 0)
    k4=calc_k(E["a4_pos"], 0)
    theta={"CO": 0, "O2": 0, "OOCO": 0, "*": 0}
    a = k3_pos*K1*K2*p["CO"]*p["O2"]/(k3_neg+k4)
    b = 1+K1*p["CO"]+K2*p["O2"]
    c = -1
    x1, x2 = solve(a, b, c)
    # print(Ead_CO,Ead_O2)
    # print('K1', K1)
    # print('K2', K2)
    # print('k3+',k3_pos)
    # print('k3-',k3_neg)
    # print('k4',k4)
    # print(x1, x2)
    if x1 > 0 and x1 < 1:
        theta["*"] = x1
    elif x2 > 0 and x2 < 1:
        theta["*"] = x2
    else:
        if debug:
            print("calc Rate ERROR")
        return 1
    theta["CO"]=p["CO"]*theta["*"]*K1
    theta["O2"]=p["O2"]*theta["*"]*K2
    theta["OOCO"]=k3_pos/(k3_neg+k4)*theta["CO"]*theta["O2"]
    r4=k4*theta["OOCO"]
    Rate = r4
    return Rate


def calc_MA(N, ls_rate):
    NA=6.02e23
    AW_au=196.967
    rate_tot = sum(ls_rate)
    MA = (NA*rate_tot)/(N*AW_au)
    return MA

def calc_barrier(Ead_CO, Ead_O2):
    # 计算活化能能垒
    Ea1 = -0.47 *(Ead_CO+Ead_O2) - 0.26
    Ea2 = -0.67*(Ead_CO+Ead_O2)-0.45
    return [Ea1,Ea2]


# ORR on pure Pt cluster
# Ref: Finding optimal surface sites on heterogeneous catalysts by counting nearest neighbors
# Fig 2, manually fit line
def calc_U_limit(CN_g, sys="Pt_ORR"):
    if CN_g < 0 or CN_g > 12:
        return "invalid GCN"
    if sys=="Pt_ORR":
        if CN_g<8.283:
            return 0.19*CN_g-0.72
        elif CN_g >= 8.283:
            return -0.17*CN_g+2.26
# Use CN_g to predict adsorption energy E_ads
# Model is: adsorbate sys

def calc_E_ads(GCN, tar="", cluster_metal="Pt"):
    if tar == "OOH":
        k,b=0.18,2.81
    elif tar == "OH":
        k,b = 0.19,-0.26
    return k*GCN+b

# HER by Pt clusters
def HER_Pt(GCN):
    if GCN < 7.7:
        i_log=0.91*GCN-8.73
    else:
        i_log=-0.94*GCN+5.53
    return i_log

# Isopropanol Dehydrogenation on Co cluster
def calc_E_delta(GCN,type="iPrO"):
    if type == "iPrO":
        a,b = 0.1,-1.2
    elif type == "H":
        a,b = 0.02,-0.36
    elif type == "TS_OHCH":
        a,b = 0.25,-1.28
    else:
        return None
    return a*GCN+b
def calc_G_delta(GCN,type):
    TS_H2 = 0.61
    TS_CH3COCH3 = 1.18
    print(type)
    if type == "iPrO":
        G=calc_E_delta(GCN, type)-(TS_H2/2-TS_CH3COCH3)
    elif type == "H":
        G = calc_E_delta(GCN, type)+TS_H2/2
    elif type == "TS_OHCH":
        G = calc_E_delta(GCN, type)-(TS_H2/2-TS_CH3COCH3)
    else:
        return None
    return G
def iso_dehydro_Co(GCN):
    # eV/K
    k_B=8.617e-5
    # eV·s
    h=4.136e-15
    T=423
    G_3=calc_G_delta(GCN,"iPrO")+calc_G_delta(GCN,"H")
    G_4=calc_G_delta(GCN,"TS_OHCH")+calc_G_delta(GCN,"H")
    delta_G = G_4-G_3
    # print(G_4, G_3, delta_G)
    R=8.617e-5
    tof = k_B*T/h*exp(-delta_G/R/T)
    return tof
    # R=8.617e-5 eV·mol⁻¹·K⁻¹
    # print(tof)
    

# CO2 electro reduction, limiting potential U_L
def calc_U_L(GCN):
    if 3.1 <= GCN <= 8.4:
        return -0.067*GCN-0.416
    elif GCN < 3.1:
        return 0.16*GCN-1.12
    elif GCN > 8.4:
        return -1.2
def overpotential_CO2(GCN):
    return -calc_U_L(GCN)
