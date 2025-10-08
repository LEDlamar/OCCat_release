import sympy as sp
from sympy import exp, symbols, Eq, solve
from descriptor.tool import calcEad


# list of Free energy
def get_pot_rev(ls_G):
    res=[-G_free for G_free in ls_G]
    return res


def tof_ORR(U3=0.8171, U4=1.7521, U5=1.2229, U6=1.128, U9=0.5509):
    SA, SB, S1, S2, S3, S4, S5, S6=symbols('SA SB S1 S2 S3 S4 S5 S6')
    # x1=O2
    x1=0.0000234
    x3 = 0
    A=1e9
    k_h=2.083e10
    kT_h=6.21e12
    e=1.6e-19
    T=298.15
    U=1
    kB=1.38e-23
    U30=U3
    U40=U4
    U50=U5
    U60=U6
    U90=U9
    
    # barrier=0.26
    # k3=kT_h*exp(-barrier*e/(kB*T))*exp(-e*0.5*(U-U30)/(kB*T))
    # k4=kT_h*exp(-barrier*e/(kB*T))*exp(-e*0.5*(U-U40)/(kB*T))
    # k5=kT_h*exp(-barrier*e/(kB*T))*exp(-e*0.5*(U-U50)/(kB*T))
    # k6=kT_h*exp(-barrier*e/(kB*T))*exp(-e*0.5*(U-U60)/(kB*T))
    # k9=kT_h*exp(-barrier*e/(kB*T))*exp(-e*0.5*(U-U90)/(kB*T))
    #
    # K3=exp(-(-U30+U)*e/(kB*T))
    # K4=exp(-(-U40+U)*e/(kB*T))
    # K5=exp(-(-U50+U)*e/(kB*T))
    # K6=exp(-(-U60+U)*e/(kB*T))
    # K9=exp(-(-U90+U)*e/(kB*T))
    #
    
    k3=(6.21*1e12)*exp(-0.26*e/T/(1.38*1e-23))*exp(-e*0.5*(U-U30)/(1.38*1e-23)/T)
    k4=(6.21*1e12)*exp(-0.26*e/T/(1.38*1e-23))*exp(-e*0.5*(U-U40)/(1.38*1e-23)/T)
    k5=(6.21*1e12)*exp(-0.26*e/T/(1.38*1e-23))*exp(-e*0.5*(U-U50)/(1.38*1e-23)/T)
    k6=(6.21*1e12)*exp(-0.26*e/T/(1.38*1e-23))*exp(-e*0.5*(U-U60)/(1.38*1e-23)/T)
    k9=(6.21*1e12)*exp(-0.26*e/T/(1.38*1e-23))*exp(-e*0.5*(U-U90)/(1.38*1e-23)/T)
    
    K3=exp(-(-U30+U)*e/kB/T)
    K4=exp(-(-U40+U)*e/kB/T)
    K5=exp(-(-U50+U)*e/kB/T)
    K6=exp(-(-U60+U)*e/kB/T)
    K9=exp(-(-U90+U)*e/kB/T)
    
    k32=k3/K3
    k42=k4/K4
    k52=k5/K5
    k62=k6/K6
    k92=k9/K9
    

    L4=Eq((k3*SA*x1-k32*S2)-(k4*S2-k42*S3)-(k9*S2-k92*SA*x3),0)
    # L4=Eq((k4*S2-k42*S3+k9*S2-k92*SA)-(k3*SA*x1-k32*S2), 0)
    L5=Eq((k4*S2-k42*S3)-(k5*S3-k52*S4), 0)
    L6=Eq((k5*S3-k52*S4)-(k6*S4-k62*SA), 0)
    L8=Eq(SA+S2+S3+S4, 1)
    sol=solve([L4, L5, L6, L8], [SA, S2, S3, S4], dict=True)
    solutions=sol[0]
    SA_sol=solutions[SA]
    S2_sol=solutions[S2]
    S3_sol=solutions[S3]
    S4_sol=solutions[S4]
    
    r3=k3*x1*SA_sol-k32*S2_sol
    r4=k4*S2_sol-k42*S3_sol
    r5=k5*S3_sol-k52*S4_sol
    r6=k6*S4_sol-k62*SA_sol
    r9=k9*S2_sol-k92*SA_sol*x3
    T1=r3+r4+r5+r6
    # print(f"k3 = {k3.evalf(6)}")
    # print(f"k4 = {k4.evalf(6)}")
    # print(f"k5 = {k5.evalf(6)}")
    # print(f"k6 = {k6.evalf(6)}")
    # print(f"k9 = {k9.evalf(6)}")
    #
    # print()
    #
    # print(f"k32 = {k32.evalf(6)}")
    # print(f"k42 = {k42.evalf(6)}")
    # print(f"k52 = {k52.evalf(6)}")
    # print(f"k62 = {k62.evalf(6)}")
    # print(f"k92 = {k92.evalf(6)}")
    
    # print(f"K3 = {K3.evalf(6)}")
    # print(f"K4 = {K4.evalf(6)}")
    # print(f"K5 = {K5.evalf(6)}")
    # print(f"K6 = {K6.evalf(6)}")
    # print(f"K9 = {K9.evalf(6)}")
    
    # print()
    # print(f"SA = {SA_sol.evalf(4)}")
    # print(f"S2 = {S2_sol.evalf(4)}")
    # print(f"S3 = {S3_sol.evalf(4)}")
    # print(f"S4 = {S4_sol.evalf(4)}")
    # print(f"r3 = {r3.evalf(6)}")
    # print(f"r4 = {r4.evalf(6)}")
    # print(f"r5 = {r5.evalf(6)}")
    # print(f"r6 = {r6.evalf(6)}")
    # print(f"r9 = {r9.evalf(6)}")
    return T1.evalf(5)



def get_orr_ads(GCN=4):
    G_OH=-289.2+279.084
    G_O=-283.54+279.084
    G_OOH=-293.672+279.084
    
    G_H2O=-14.22
    G_H2=-6.8
    
    RCA=None
    E_O2=calcEad("O2", "Pt", GCN, RCA)
    E_OOH=calcEad("OOH", "Pt", GCN, RCA)
    E_O=calcEad("O", "Pt", GCN, RCA)
    E_H2O=calcEad("H2O", "Pt", GCN, RCA)
    E_OH=calcEad("OH", "Pt", GCN, RCA)
    E_H2O2=calcEad("H2O2", "Pt", GCN, RCA)
    
    H3=E_OOH
    H4=E_O+E_H2O-E_OOH
    H4_new=E_O+G_O+G_H2O-E_OOH-G_OOH-G_H2/2
    H5=E_OH-E_O
    H5_new=E_OH+G_OH-E_O-G_O
    
    H6=G_H2O-E_OH-G_OH
    H9=E_H2O2-E_OOH
    return [H3, H4, H5, H6, H9]
# print(tof_ORR())
