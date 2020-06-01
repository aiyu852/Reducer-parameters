# -*-coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号


# Alpha = []
# for i in range(100, 850):
#     alpha_i = 90.0-i/10+55.561
#     # 转换为弧度制
#     Alpha.append(alpha_i*np.pi/180)

# # 余弦定理
# L_a = 260
# L_b = 84.876
# L_c = []
# for i in range(750):
#     Num = L_a**2+L_b**2-2*L_a*L_b*np.cos(Alpha[i])
#     L_c.append(np.sqrt(Num))
# # print(L_c)

# # 求出L_c上的高即为执行器力臂
# h = []
# Ave = 0
# Max = 0
# Max_x = 0
# Min = 50
# Min_x = 0
# for i in range(750):
#     S_2 = L_a*L_b*np.sin(Alpha[i])
#     h.append(S_2/L_c[i])
#     Ave = Ave + h[i]
#     if h[i] > Max:
#         Max = h[i]
#         Max_x = i/10+10
#     if h[i] < Min:
#         Min = h[i]
#         Min_x = i/10+10
# Ave = Ave/750
# # print(Ave)

# # 创建x轴
# x = []
# for i in range(100, 850):
#     x.append(i/10)


# plt.plot(x, h)
# plt.grid(True)

# plt.xlabel("大腿角度(°)")
# plt.ylabel("力臂长度(mm)")

# plt.text(7,85,"平均力臂长度为 %f" % Ave)

# plt.annotate('最大力臂长度为 \n%f' % Max, xy=(Max_x, Max), xytext=(65, 70),
#              arrowprops=dict(facecolor='black'))
# plt.annotate('最小力臂长度为 \n%f' % Min, xy=(Min_x, Min), xytext=(30, 52),
#              arrowprops=dict(facecolor='black'))

# plt.show()

# Alpha = []
# for i in range(1050):
#     alpha_i = 26.667+i/10
#     # 转换为弧度制
#     Alpha.append(alpha_i*np.pi/180)

# # 余弦定理
# L_a = 340
# L_b = 80
# L_c = []
# for i in range(1050):
#     Num = L_a**2+L_b**2-2*L_a*L_b*np.cos(Alpha[i])
#     L_c.append(np.sqrt(Num))
# # print(L_c)

# # 求出L_c上的高即为执行器力臂
# h = []
# Ave = 0
# Max = 0
# Max_x = 0
# Min = 50
# Min_x = 0
# for i in range(1050):
#     S_2 = L_a*L_b*np.sin(Alpha[i])
#     h.append(S_2/L_c[i])
#     Ave = Ave + h[i]
#     if h[i] > Max:
#         Max = h[i]
#         Max_x = i/10
#     if h[i] < Min:
#         Min = h[i]
#         Min_x = i/10
# Ave = Ave/1050
# # print(Ave)

# # 创建x轴
# x = []
# for i in range(1050):
#     x.append(i/10)


# plt.plot(x, h)
# plt.grid(True)

# plt.xlabel("小腿角度(°)")
# plt.ylabel("力臂长度(mm)")

# plt.text(-4,78,"平均力臂长度 \n 为 %f" % Ave)

# plt.annotate('最大力臂长度为 \n%f' % Max, xy=(Max_x, Max), xytext=(50, 70),
#              arrowprops=dict(facecolor='black'))
# plt.annotate('最小力臂长度为 \n%f' % Min, xy=(Min_x, Min), xytext=(20, 50),
#              arrowprops=dict(facecolor='black'))

# plt.show()

###
# 三
###

# Alpha = []
# for i in range(100, 850):
#     alpha_i = 90.0-i/10+55.561
#     # 转换为弧度制
#     Alpha.append(alpha_i*np.pi/180)
# Alpha1 = []
# for i in range(1050):
#     alpha_i1 = 26.667+i/10
#     # 转换为弧度制
#     Alpha1.append(alpha_i1*np.pi/180)
# # 余弦定理
# L_a = 260
# L_b = 84.876
# L_thigh = []
# L_a1 = 340
# L_b1 = 80
# L_shank = []
# for i in range(750):
#     Num = L_a**2+L_b**2-2*L_a*L_b*np.cos(Alpha[i])
#     L_thigh.append(np.sqrt(Num)-140)
# for i in range(1050):
#     Num = L_a1**2+L_b1**2-2*L_a1*L_b1*np.cos(Alpha1[i])
#     L_shank.append(np.sqrt(Num)-208)
# # 求出L_c上的高即为执行器力臂
# F0 = 1854
# F_max = []
# turningx = 0
# for i in range(750):
#     P1 = 1.2*7.8**4/L_thigh[i]**2*10**4
#     if F0 < P1:
#         F_max.append(F0)
#         if F_max[i-1] != F_max[i]:
#             turningx = i/10
#             print(turningx)
#     else:
#         F_max.append(P1)
# F1_max = []
# turningx1 = 0
# for i in range(1050):
#     P1 = 1.2*7.8**4/L_shank[i]**2*10**4
#     if F0 < P1:
#         F1_max.append(F0)
#     else:
#         F1_max.append(P1)
#         if F1_max[i-1] == F0:
#             turningx1 = 105-i/10
# # F_max.reverse()
# F1_max.reverse()
# # 创建x轴
# x = []
# for i in range(750):
#     x.append(i/10)
# x1 = []
# for i in range(1050):
#     x1.append(i/10)

# plt.figure()
# plt.plot(x, F_max)
# plt.xticks(range(0, 75, 5))
# plt.xlim([0, 75])
# plt.ylim([1200, 1900])
# plt.plot(turningx, F0, 'om')
# plt.annotate('转折点\n%.1f°' % turningx, xy=(
#     turningx, F0), xytext=(turningx, F0-150), arrowprops=dict(arrowstyle='->'))
# plt.text(0, 1250, "最小值为%.1f N" % F_max[0])
# plt.text(50, 1860, "最大值为%.1f N" % F_max[749])
# plt.grid(True)
# plt.xlabel("大腿角度(°)")
# plt.ylabel("最大容许推力(N)")

# plt.figure()
# plt.plot(x1, F1_max)
# plt.xticks(range(0, 105, 5))
# plt.xlim([0, 105])
# plt.ylim([1200, 1900])
# plt.plot(turningx1, F0, 'om')
# plt.annotate('转折点\n%.1f°' % turningx1, xy=(
#     turningx1, F0), xytext=(turningx1, F0-150), arrowprops=dict(arrowstyle='->'))
# plt.text(5, 1230, "最小值为%.1f N" % F1_max[0])
# plt.text(65, 1860, "最大值为%.1f N" % F1_max[1049])
# plt.grid(True)
# plt.xlabel("小腿角度(°)")
# plt.ylabel("最大容许推力(N)")

# plt.show()

Alpha = []
for i in range(100, 850):
    alpha_i = 90.0-i/10+55.561
    # 转换为弧度制
    Alpha.append(alpha_i*np.pi/180)
Alpha1 = []
for i in range(1050):
    alpha_i1 = 26.667+i/10
    # 转换为弧度制
    Alpha1.append(alpha_i1*np.pi/180)
# 余弦定理
L_a = 260
L_b = 84.876
L_thigh = []
L_a1 = 340
L_b1 = 80
L_shank = []
for i in range(750):
    Num = L_a**2+L_b**2-2*L_a*L_b*np.cos(Alpha[i])
    L_thigh.append(np.sqrt(Num)-140)
for i in range(1050):
    Num = L_a1**2+L_b1**2-2*L_a1*L_b1*np.cos(Alpha1[i])
    L_shank.append(np.sqrt(Num)-208)
# 求出L_c上的高即为执行器力臂
v0 = 11007
v_max = []
turningx = 0
for i in range(750):
    n1 = 3.4*7.8/L_thigh[i]**2*10**7
    if v0 < n1:
        v_max.append(v0)
        if v_max[i-1] != v_max[i]:
            turningx = i/10
            print(turningx)
    else:
        v_max.append(n1)
v1_max = []
turningx1 = 0
for i in range(1050):
    n1 = 3.4*7.8/L_shank[i]**2*10**7
    if v0 < n1:
        v1_max.append(v0)
    else:
        v1_max.append(n1)
        if v1_max[i-1] == v0:
            turningx1 = 105-i/10
# F_max.reverse()
v1_max.reverse()
# 创建x轴
x = []
for i in range(750):
    x.append(i/10)
x1 = []
for i in range(1050):
    x1.append(i/10)

plt.figure()
plt.plot(x, v_max)
plt.xticks(range(0, 75, 5))
plt.xlim([0, 75])
plt.plot(turningx, v0, 'om')
plt.annotate('转折点\n%.1f°' % turningx, xy=(
    turningx, v0), xytext=(turningx, v0-500), arrowprops=dict(arrowstyle='->'))
plt.text(0, v_max[0], "最小值为%.1f r/min" % v_max[0])
plt.text(50, v0+40, "最大值为%.1f r/min" % v_max[749])
plt.grid(True)
plt.xlabel("大腿角度(°)")
plt.ylabel("电机最大转速(r/min)")

plt.figure()
plt.plot(x1, v1_max)
plt.xticks(range(0, 105, 5))
plt.xlim([0, 105])
plt.plot(turningx1, v0, 'om')
plt.annotate('转折点\n%.1f°' % turningx1, xy=(
    turningx1, v0), xytext=(turningx1, v0-500), arrowprops=dict(arrowstyle='->'))
plt.text(5, v1_max[0], "最小值为%.1f r/min" % v1_max[0])
plt.text(65, v0+40, "最大值为%.1f r/min" % v1_max[1049])
plt.grid(True)
plt.xlabel("小腿角度(°)")
plt.ylabel("电机最大转速(r/min)")

plt.show()