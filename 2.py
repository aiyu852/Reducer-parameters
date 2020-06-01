# -*-coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

r_p = 25.5
z_p = 16
z_c = z_p-1
i_h = z_p/z_c
r_rp = 2.5
a = 0.9
k1 = a*z_p/r_p
b = 4.5


# 修形参数
# D_Rp = 0.0572
# D_Rrp = -0.0472

# 最佳修形量公式

# D_j = 0.01
# D_Rp = D_j/(1-(1-k1**2)**(1/2))
# D_Rrp = -D_j*(1-k1**2)**(1/2)/(1-(1-k1**2)**(1/2))
d_Rp = [0]*1000
d_Rrp = [0]*1000
d_Rp[0] = 0.01
d_Rrp[0] = 0

# print(D_Rp)
# print(D_Rrp)
print("\n")


def calc_delta_max(r_p, k1, z_p, r_rp, b, F_max):
    rou = r_p*(1+k1**2-2*k1**2)**(3/2) / \
        np.abs(k1*(z_p+1)*k1-1-z_p*k1**2)+r_rp  # +D_Rrp
    miu = 0.3
    E = 206000
    c = 4.99/1000*np.sqrt(2*(1-miu**2)/E*F_max/b*2*rou*r_rp/(r_rp+rou))
    w_max = 2*(1-miu**2)/E*F_max/np.pi/b*(2/3+np.log(16*r_rp*rou/c**2))

    # f_max
    J = np.pi * 3**4/64
    f_max = F_max*31*16**3/48/64/E/J

    delta_max = f_max + w_max
    return delta_max


def clac_delta_i(delta_max, k1):
    delta = [0]*180
    for i in range(180):
        j = np.pi/180*i
        delta[i] = delta_max*np.sin(j)/np.sqrt(1+k1**2-2*k1*np.cos(j))
    return delta


def find_phi_mn(delta, Q):
    phi_mn = []
    for i in range(180):
        if delta[i] > Q[i]:
            phi_mn.append(i)
    return phi_mn


def find_mANDn(phi_mn):
    for i in range(48):
        phi_m = i * 360/z_p
        if phi_m >= phi_mn[0]:
            for j in range(len(phi_mn)):
                if (phi_m+j*360/z_p) > phi_mn[-1]:
                    m_n = j
                    break
            break
    return phi_m, m_n


def calc_F_max(m_n, a, z_c, phi_m, z_p, k1, D_Rp, D_Rrp, delta_max):
    U = 0
    for i in range(m_n):
        r_c = a*z_c
        j = phi_m+i*360/z_p
        k = j*np.pi/180
        S_i = 1/np.sqrt(1+k1**2-2*k1*np.cos(k))
        # 计算l_i
        l_i = r_c*np.sin(k)*S_i
        # delta_phi_i
        delta_phi_i = D_Rrp*(1-np.sin(k)*S_i)+D_Rp*S_i * \
            (1-k1*np.cos(k)-np.sqrt(1-k1**2)*np.sin(k))
        U += l_i*(l_i/r_c-delta_phi_i/delta_max)
    F_max = 0.55 * T / U
    return F_max


def calc_Fi(j):
    k = j*np.pi/180
    # delta_i
    delta_i = delta_max*np.sin(k)/np.sqrt(1+k1**2-2*k1*np.cos(k))
    # Q_i
    S_i = 1/np.sqrt(1+k1**2-2*k1*np.cos(k))
    Q_i = D_Rrp*(1-np.sin(k)*S_i)+D_Rp*S_i * \
        (1-k1*np.cos(k)-np.sqrt(1-k1**2)*np.sin(k))
    # Fi
    if (delta_i-Q_i) < 0:
        Fi = 0
    else:
        Fi = (delta_i-Q_i)/delta_max*F_max[-1]
    phi_i = j
    Fii = Fi / 14000
    return Fi, phi_i, Fii, delta_i, Q_i


def clac_sigma(r_p, k1, z_p, r_rp, b_c, Fi):
    rou = r_p*(1+k1**2-2*k1**2)**(3/2) / \
        np.abs(k1*(z_p+1)*k1-1-z_p*k1**2)+r_rp
    rou_ei = np.abs(rou*r_rp/(rou-r_rp))
    sigma = 0.418*(206000*Fi/b_c/rou_ei)**(1/2)
    return sigma


# 计算初始间隙
'''
该初始间隙是u型曲线
'''
sigma_min = 1000
d_rp = 0

for n in range(1000):
    if n >= 1:
        d_Rp[n] = d_Rp[n-1]+0.0001
        d_Rrp[n] = d_Rrp[n-1]-0.0001
    D_Rp = d_Rp[n]
    D_Rrp = d_Rrp[n]

    Q = [0]*180
    for i in range(180):
        j = i*np.pi/180
        S_i = 1/np.sqrt(1+k1**2-2*k1*np.cos(j))
        Q[i] = D_Rrp*(1-np.sin(j)*S_i)+D_Rp*S_i * \
            (1-k1*np.cos(j)-np.sqrt(1-k1**2)*np.sin(j))
    # plt.plot(Q, label=r'$\Delta\phi_i$')
    # plt.grid(True)
    # plt.show()

    # 待啮合点的法线方向位移
    '''
    该曲线应该是n型的
    phi0 = 0*np.pi/180
    '''

    # 计算F_max[0]
    T = 8460  # N·mm
    F_a = 2.2*T/k1/z_c/r_p
    F_b = 0.55*T/a/z_c
    F_max = [0]*100
    F_max[0] = (F_a + F_b) / 2

    # phi_m = 0
    # phi_n = 0
    # m_n = 0

    for p in range(1, 100):

        delta_max = calc_delta_max(r_p, k1, z_p, r_rp, b, F_max[p-1])

        delta = clac_delta_i(delta_max, k1)

        phi_mn = find_phi_mn(delta, Q)

        phi_m, m_n = find_mANDn(phi_mn)

        F_max[p] = calc_F_max(m_n, a, z_c, phi_m, z_p,
                              k1, D_Rp, D_Rrp, delta_max)
        delta_F_max = np.abs(F_max[p]-F_max[p-1])
        # print(p)
        # print(F_max[p])
        # print(delta_F_max)
        # print(F_max[p]/1000)
        if delta_F_max < F_max[p]/1000:
            F_max[-1] = (F_max[p] + F_max[p-1]) / 2
            break

    # print(F_max[-1])
    delta_max = calc_delta_max(r_p, k1, z_p, r_rp, b, F_max[-1])
    delta = clac_delta_i(delta_max, k1)
    phi_mn = find_phi_mn(delta, Q)
    phi_m, m_n = find_mANDn(phi_mn)

    Fi = []
    Fii = []
    phii = []

    F_i = []
    F_ii = []
    phi_i = []

    sigmai = []
    sigmaii = []
    sigma_i = []
    sigma_ii = []

    # print(phi_mn)
    for i in range(180):
        o, p, q, r, s = calc_Fi(i)
        F_i.append(o)
        phi_i.append(p)
        F_ii.append(q)
        sigma_i.append(clac_sigma(r_p, k1, z_p, r_rp, b, F_i[i]))
        sigma_ii.append(clac_sigma(r_p, k1, z_p, r_rp, b, F_i[i])/60000)

    for i in range(m_n):
        j = phi_m+i*360/z_p
        o, p, q, delta_i, Q_i = calc_Fi(j)
        # plot
        Fi.append(o)
        # print(i+1)
        # print(o)
        phii.append(p)
        Fii.append(q)
        sigmai.append(clac_sigma(r_p, k1, z_p, r_rp, b, Fi[i]))
        # print(sigmai[i])
        sigmaii.append(clac_sigma(r_p, k1, z_p, r_rp, b, Fi[i])/60000)

        # plt.scatter(j, delta_i)
        # plt.scatter(j, Q_i)

    # plt.plot(delta, label=r'$\Delta\delta_i$')
    # plt.xlim([0, 180])
    # plt.ylim([0, 0.01])
    # plt.xlabel("角参量"+r'$\phi$/rad')
    # plt.ylabel("初始间隙"+r'$\Delta\phi_i$/mm'+'\n'+"变形量"+r'$\Delta\delta_i$/mm')
    # plt.legend()

    # plt.figure()
    # plt.plot(phi_i, F_ii)
    # plt.plot(phi_i, sigma_i)
    # plt.scatter(phii, Fii)
    # plt.scatter(phii, sigmai)
    # plt.show()
    sum = 0
    for m in range(m_n):
        sum += sigmai[m]
    if sigma_min > sum/m_n:
        sigma_min = sum/m_n
        d_rp = D_Rp
        print(sigmai)
print(sigma_min)
print(d_rp)