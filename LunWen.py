import math
import sympy
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import json


def drawlimit():
    # 限制坐标
    plt.xlim((0, 1))
    plt.ylim((0, 1))
    # x=y
    x = np.linspace(0, 1)
    plt.plot(x, x)
    # 转化条件


def init1():
    # 输入数据
    feed = input("筛板塔：进料量  吨/年；")  # 55000
    if feed == "":
        print("\033[0;31;40m使用调试数据30000，40，97\033[0m")
        feed = 30000
        feed_composition = 40
        distillate = 97
    else:
        feed_composition = input("进料组成  %；")  # 40
        distillate = input("塔顶馏出液组成  %")  # 96
    return feed, feed_composition, distillate


def init3():
    # 输入数据
    xq = input("xq")
    if xq == "":
        print("\033[0;31;40m使用调试数据0.44,0.657\033[0m")
        xq = 0.44
        yq = 0.657
    else:
        yq = input("yq")
    return xq, yq


def init2():
    # 输入数据
    print("由图解法可知yd，yw")
    YD = input("yd")  # 95.34
    if YD == "":
        print("\033[0;31;40m使用调试数据0.988,0.024,6,10\033[0m")
        YD = 0.988
        YW = 0.024
        Nj = 6
        Nt = 10
    else:
        YW = input("yW")
        Nj = input("N精")
        Nt = input("N提")
    return YD, YW, Nj, Nt


def init4():
    # 输入数据
    td = input("最小的温度")  # 55000
    if td == "":
        print("\033[0;31;40m使用调试数据81.023，110.215\033[0m")
        td = 81.023
        tw = 110.215
    else:
        tw = input("最大的温度")  #
    return td, tw


def uuu(td, tw, xd, xw):
    td = float(td)
    tw = float(tw)
    list = [[80, 0.308, 0.311], [90, 0.279, 0.286], [100, 0.255, 0.264], [110, 0.233, 0.254], [120, 0.215, 0.228]]
    for i in range(len(list) - 1):
        if list[i + 1][0] + 10 > td > list[i][0]:
            befor_info = list[i]
            after_info = list[i + 1]
            duben = befor_info[1] + ((after_info[1] - befor_info[1]) / 10) * (td - befor_info[0])
            dujben = befor_info[2] + ((after_info[2] - befor_info[2]) / 10) * (td - befor_info[0])
        if list[i + 1][0] + 10 > tw > list[i][0]:
            befor_info = list[i]
            after_info = list[i + 1]
            wuben = befor_info[1] + ((after_info[1] - befor_info[1]) / 10) * (tw - befor_info[0])
            wujben = befor_info[2] + ((after_info[2] - befor_info[2]) / 10) * (tw - befor_info[0])
    fduben = float(f'{duben:.3f}')
    fdujben = float(f'{dujben:.3f}')
    fwuben = float(f'{wuben:.3f}')
    fwujben = float(f'{wujben:.3f}')
    uld = fdujben * (1 - xd) + fduben * xd
    ulw = fwujben * (1 - xw) + fwuben * xw
    fuld = float(f'{uld:.3f}')
    fulw = float(f'{ulw:.3f}')
    u = (fuld + fulw) / 2
    print(f'原料液与塔顶、塔底产品的摩尔分率\ntd:u苯{fduben}u甲苯{fdujben},tw:u苯{fwuben}u甲苯{fwujben}，uld:{fuld},ulw:{fulw},')
    return u


def draw1(Xw, Xd, Xf):
    # 图1
    drawlimit()
    tixingtux = [1.0, 0.78], [0.78, 0.581], [0.581, 0.412], [0.412, 0.258], [0.258, 0.130], [0.130, 0.0], [Xf, Xf]
    tixingtuy = [1.0, 0.9], [0.9, 0.777], [0.777, 0.633], [0.633, 0.456], [0.456, 0.262], [0.262, 0.0], [Xf, 1]
    for i in range(len(tixingtux)):
        plt.plot(tixingtux[i], tixingtuy[i], 'k-')
    plt.title('Schematic method')
    plt.show()


def draw2(Xw, Xd, Xf, R):
    # 图2
    drawlimit()
    # 算点求位置
    yb = Xd / (R + 1)
    dy = (Xd - yb) * Xf / Xd + yb
    tujiex = [0, Xd], [Xf, Xf], [Xf, Xw]
    tujiey = [yb, Xd], [Xf, dy], [dy, Xw]
    for i in range(len(tujiex)):
        plt.plot(tujiex[i], tujiey[i])
    print(f"由化工原理下册书p24页图1-20可得相对坐标位置\nb点{(0, yb)}a点{(Xd, Xd)}c点{(Xw, Xw)}e点{(Xf, Xf)}d点{(Xf, dy)}")
    tixingtux = [1.0, 0.78], [0.78, 0.581], [0.581, 0.412], [0.412, 0.258], [0.258, 0.130], [0.130,
                                                                                             0.0],  # ,[Xd,Xd],[Xw,Xw]
    tixingtuy = [1.0, 0.9], [0.9, 0.777], [0.777, 0.633], [0.633, 0.456], [0.456, 0.262], [0.262, 0.0],  # , [0,1],[0,1]
    for i in range(len(tixingtux)):
        plt.plot(tixingtux[i], tixingtuy[i], 'k-')
    plt.title('Schematic method')
    plt.show()


def draw3(Xw, Xd, Xf, R):
    # 图2
    drawlimit()
    tixingtux = [1.0, 0.78], [0.78, 0.581], [0.581, 0.412], [0.412, 0.258], [0.258, 0.130], [0.130, 0.0], [Xd, Xd], [Xw,
                                                                                                                     Xw]
    tixingtuy = [1.0, 0.9], [0.9, 0.777], [0.777, 0.633], [0.633, 0.456], [0.456, 0.262], [0.262, 0.0], [0, 1], [0, 1]
    for i in range(len(tixingtux)):
        plt.plot(tixingtux[i], tixingtuy[i], 'k-')
    plt.title('Schematic method')
    plt.show()


def draw4(yd, yw, xd, xw):
    plt.xlim((0, 1))
    plt.ylim((0, 120))
    benx = [0.00, 1.00, 3.00, 5.00, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0,
            80.0, 85.0, 90.0, 95.0, 97.0, 99.0, 100.0]
    jiabeny = [0.00, 2.50, 7.11, 11.2, 20.8, 29.4, 37.2, 44.2, 50.7, 56.6, 61.9, 66.7, 71.3, 75.5, 79.1, 82.5, 85.7,
               88.5, 91.2, 93.6, 95.9, 98.0, 98.8, 99.61, 100.0]
    t = [110.56, 109.91, 108.79, 107.61, 105.05, 102.79, 100.75, 98.84, 97.13, 95.58, 94.09, 92.69, 91.40, 90.11,
         80.80, 87.63, 86.52, 85.44, 84.40, 83.33, 82.25, 81.11, 80.66, 80.21, 80.01]
    for i in range(len(benx) - 1):
        xtx = [benx[i] * 0.01, benx[i + 1] * 0.01]
        xty = [t[i], t[i + 1]]
        plt.plot(xtx, xty, 'k-')
        xtx = [jiabeny[i] * 0.01, jiabeny[i + 1] * 0.01]
        xty = [t[i], t[i + 1]]
        plt.plot(xtx, xty, 'k-')
    x = [xd, xw]
    y = [yd, yw]
    for i in range(len(x)):
        xtx = [x[i], x[i]]
        xty = [y[i], 120]
        plt.plot(xtx, xty, 'k-')
    plt.title('t-x-y')
    plt.show()


number_list = init1()
# 转化数据
FloatFeed = float(str(number_list[0])[0:2]) * 0.1
FloatFeedComposition = float(int(number_list[1]) * 0.01)
FloatDistillate = float(int(number_list[2]) * 0.01)
# 原料液与塔顶、塔底产品的摩尔分率
Xf = (FloatFeedComposition / 78.11) / ((FloatFeedComposition / 78.11) + ((1 - FloatFeedComposition) / 92.13))
Xd = (FloatDistillate / 78.11) / ((FloatDistillate / 78.11) + ((1 - FloatDistillate) / 92.13))
Xw = (0.01 / 78.11) / ((0.01 / 78.11) + ((1 - 0.01) / 92.13))
f_Xf = float(f'{Xf:.3f}')
f_Xd = float(f'{Xd:.3f}')
f_Xw = float(f'{Xw:.3f}')
print(f'原料液与塔顶、塔底产品的摩尔分率\nXf:{f_Xf},Xd:{f_Xd},Xw:{f_Xw}')
# 原料液与塔顶、塔底产品的平均摩尔分率
Mf = f_Xf * 78.11 + (1 - f_Xf) * 92.13
Md = f_Xd * 78.11 + (1 - f_Xd) * 92.13
Mw = f_Xw * 78.11 + (1 - f_Xw) * 92.13
f_Mf = float(f'{Mf:.2f}')
f_Md = float(f'{Md:.2f}')
f_Mw = float(f'{Mw:.2f}')
print(f'原料液与塔顶、塔底产品的平均摩尔分率\nMf:{f_Mf},Md:{f_Md},Mw:{f_Mw}')
# 物料衡算
F_ = (FloatFeed * math.pow(10, 7)) / 7200
f_F_ = float(f'{F_:.0f}')
F = f_F_ / f_Mf
f_F = float(f'{F:.2f}')
D, W = sympy.symbols('D W')
number_dict = sympy.solve([D + W - f_F, f_Xd * D + Xw * W - f_F * f_Xf], [D, W])
f_D = float(f'{number_dict[D]:.2f}')
f_W = float(f'{number_dict[W]:.2f}')
print(f'物料衡算\nF‘:{f_F_},F:{f_F},D:{f_D},W:{f_W}')
draw1(f_Xw, f_Xd, f_Xf)
number_list = init3()
# 求最小回流比及操作线回流比
f_yq = float(number_list[0])
R_min = (f_Xd - f_yq) / (float(number_list[1]) - f_Xf)
f_R_min = float(f'{R_min:.2f}')
R = 1.7 * f_R_min
f_R = float(f'{R:.2f}')
print(f'求最小回流比及操作线回流比\nR_min:{f_R_min},R:{f_R}')
# 求精馏塔的气、液相负荷
L = f_R * f_D
f_L = float(f'{L:.2f}')
V = (f_R + 1) * f_D
f_V = float(f'{V:.2f}')
L_ = f_L + f_F
f_L_ = float(f'{L_:.2f}')
f_V_ = f_V
print(f'求精馏塔的气、液相负荷\nL:{f_L},V:{f_V},L’:{f_L_},V‘:{f_V_}')
# 求操作线方程
y = "{:.3}x+{:.3}".format(f_L / f_V, f_D / f_V * f_Xd)
y_ = "{:.3}x'-{:.3}".format(f_L_ / f_V_, f_W / f_V_ * 0.012)
print(f'精馏段操作线方程为\ny={y}')
print(f'提馏段操作线方程为\ny’={y_}')
# 图解法理论板数
draw2(f_Xw, f_Xd, f_Xf, f_R)
# 实际板层数的求取
draw3(f_Xw, f_Xd, f_Xf, f_R)
number_list = init2()
f_yd = float(number_list[0])
f_yw = float(number_list[1])
ad = (f_yd * (1 - f_Xd)) / ((1 - f_yd) * f_Xd)
aw = (f_yw * (1 - f_Xw)) / ((1 - f_yw) * f_Xw)
a = math.sqrt(ad * aw)
f_a = float(f'{a:.3f}')
draw4(f_yd, f_yw, f_Xd, f_Xw)
other_number_list = init4()
f_td = float(other_number_list[0])
f_tw = float(other_number_list[1])
u = uuu(f_td, f_tw, f_Xd, f_Xw)
f_u = float(f'{u:.3f}')
Et = 0.49 * math.pow(f_u * f_a, -0.245)
f_Et = float(f'{Et:.3f}')
Nj = float(number_list[2]) // Et
f_Nj = float(f'{Nj + 1:.3f}')
Nt = float(number_list[3]) // Et
f_Nt = float(f'{Nt + 1:.3f}')
print(f'实际板层数的求取\nad={ad},aw={aw},a={f_a},U={f_u},Et={f_Et},Nj={f_Nj},Nt={f_Nt},总板数={f_Nt + f_Nj}')
# 进料板与塔顶、塔底平均摩尔质量计算
Mvdm = f_Xd * 78.11 + (1 - f_Xd) * 92.13
f_Mvdm = float(f'{Mvdm:.2f}')
Mldm = f_yd * 78.11 + (1 - f_yd) * 92.13
f_Mldm = float(f'{Mldm:.2f}')
Mvfm = f_yq * 78.11 + (1 - f_yq) * 92.13
f_Mvfm = float(f'{Mvfm:.2f}')
Mlfm = f_Xf * 78.11 + (1 - f_Xf) * 92.13
f_Mlfm = float(f'{Mlfm:.2f}')
Mvwm = f_yw * 78.11 + (1 - f_yw) * 92.13
f_Mvwm = float(f'{Mvwm:.2f}')
Mlwm = f_Xw * 78.11 + (1 - f_Xw) * 92.13
f_Mlwm = float(f'{Mlwm:.2f}')
Mvm = (f_Mvdm + f_Mvfm) / 2
f_Mvm = float(f'{Mvm:.2f}')
Mlm = (f_Mldm + f_Mlfm) / 2
f_Mlm = float(f'{Mlm:.2f}')
M_vm = (f_Mvwm + f_Mvfm) / 2
f_M_vm = float(f'{M_vm:.2f}')
M_lm = (f_Mlwm + f_Mlfm) / 2
f_M_lm = float(f'{M_lm:.2f}')
print(f'进料板与塔顶、塔底平均摩尔质量计算'
      f'\nMVDm={f_Mvdm},MLDm={f_Mldm},MVFm={f_Mvfm},MLFm={f_Mlfm},'
      f'MVWm={f_Mvwm},MLWm={f_Mlwm},MVm={f_Mvm},MLm={f_Mlm},'
      f'M’Vm={f_M_vm},M’Lm={f_M_lm}')
# 塔内各段操作条件和物性数据表
f_Pd = 105.3
f_SP = 0.7
Pw = f_Pd + (f_Nt + f_Nj) * f_SP
f_Pw = float(f'{Pw:.2f}')
Pf = f_Pd + f_SP * f_Nj
f_Pf = float(f'{Pf:.2f}')
Pm = (f_Pd + f_Pf) / 2
f_Pm = float(f'{Pm:.2f}')
pwm = (f_Pw + f_Pf) / 2
f_pwm = float(f'{pwm:.2f}')
print(f'塔内各段操作条件和物性数据表\nPd={f_Pd},三角P={f_SP},Pw={f_Pw},Pf={f_Pf},Pm={f_Pm},Pwm={f_pwm}')


# 操作温度计算
def wendu(zufen, yali):  # 组分，压力
    ba = 6.023
    bb = 1206.35
    bc = 220.24
    ja = 6.078
    jb = 1343.94
    jc = 219.58
    t = 80
    while 1:
        pa = 10 ** (ba - (bb / (t + bc)))
        pb = 10 ** (ja - (jb / (t + jc)))
        p = pa * zufen + pb * (1 - zufen)
        if p > yali:
            t = float(f'{t:.1f}')
            break
        else:
            t += 0.1
    return t


p = 283.2871 * f_Xd + 123.645 * (1 - f_Xd)
p = float(f'{p:.1f}')
print(f'操作温度计算\nPd={p}')
td = wendu(f_Xd, f_Pd)
tf = wendu(f_Xf, f_Pf)
tw = wendu(f_Xw, f_Pw)
importtd = td
print(f'td塔顶温度={td},tf进料板温度={tf},tw塔釜温度={tw}')
f_tdm = (td + tf) / 2
f_twm = (tw + tf) / 2
print(f'tdm精馏段平均温度={f_tdm},twm提馏段平均温度={f_twm}')
# 气相平均密度计算
pvm = f_Pm * f_Mvm / (8.314 * (f_tdm + 273.15))
p_vm = f_pwm * f_M_vm / (8.314 * (f_twm + 273.15))
f_pvm = float(f'{pvm:.2f}')
f_p_vm = float(f'{p_vm:.2f}')
print(f'气相平均密度计算\npvm精馏段={f_pvm},p“vm提馏段={f_p_vm}')


# 液相相平均密度计算
def yyy(any):
    any = float(any)
    list = [[80, 815, 810], [90, 803.9, 800.2], [100, 792.5, 790.3], [110, 780.3, 780.3], [120, 768.9, 770]]
    for i in range(len(list) - 1):
        if list[i + 1][0] + 10 > any > list[i][0]:
            befor_info = list[i]
            after_info = list[i + 1]
            pa = befor_info[1] + ((after_info[1] - befor_info[1]) / 10) * (td - befor_info[0])
            pb = befor_info[2] + ((after_info[2] - befor_info[2]) / 10) * (td - befor_info[0])
    pa = float(f'{pa:.1f}')
    pb = float(f'{pb:.1f}')
    return pa, pb


plist = yyy(td)
pldm = 1 / ((f_Xd / plist[0]) + ((1 - f_Xd) / plist[1]))
f_pldm = float(f'{pldm:.1f}')
print(f'液相平均密度计算\n塔顶液相\npldm={f_pldm},pa={plist[0]},pb={plist[1]}')
plist = yyy(tf)
plfm = 1 / ((f_Xf / plist[0]) + ((1 - f_Xf) / plist[1]))
f_plfm = float(f'{plfm:.1f}')
print(f'进料板液\nplfm={f_plfm},pa={plist[0]},pb={plist[1]}')
plist = yyy(tw)
plwm = 1 / ((f_Xw / plist[0]) + ((1 - f_Xw) / plist[1]))
f_plwm = float(f'{plwm:.1f}')
print(f'提馏段液相\nplwm={f_plwm},pa={plist[0]},pb={plist[1]}')
plm = (f_pldm + f_plfm) / 2
p_lm = (f_plwm + f_plfm) / 2
f_plm = float(f'{plm:.1f}')
f_p_lm = float(f'{p_lm:.1f}')
print(f'plm塔釜液相={f_plm},p’lm提馏段液相={f_p_lm}')


# 液体表面张力计算
def zzz(any):
    any = float(any)
    list = [[80, 21.7, 21.69], [90, 20.06, 20.59], [100, 18.85, 19.94], [110, 17.66, 18.41], [120, 16.49, 17.31]]
    for i in range(len(list) - 1):
        if list[i + 1][0] + 10 > any > list[i][0]:
            befor_info = list[i]
            after_info = list[i + 1]
            pa = befor_info[1] + ((after_info[1] - befor_info[1]) / 10) * (td - befor_info[0])
            pb = befor_info[2] + ((after_info[2] - befor_info[2]) / 10) * (td - befor_info[0])
    pa = float(f'{pa:.2f}')
    pb = float(f'{pb:.2f}')
    return pa, pb


plist = zzz(td)
aldm = ((f_Xd * plist[0]) + ((1 - f_Xd) * plist[1]))
f_aldm = float(f'{aldm:.1f}')
print(f'液体表面张力计算\n塔顶液相\naldm={f_aldm},pa={plist[0]},pb={plist[1]}')
plist = zzz(tf)
alfm = ((f_Xf * plist[0]) + ((1 - f_Xf) * plist[1]))
f_alfm = float(f'{alfm:.1f}')
print(f'进料板液\nalfm={f_alfm},pa={plist[0]},pb={plist[1]}')
plist = zzz(tw)
alwm = ((f_Xw * plist[0]) + ((1 - f_Xw) * plist[1]))
f_alwm = float(f'{alwm:.1f}')
print(f'塔釜液相\nalwm={f_alwm},pa={plist[0]},pb={plist[1]}')
alm = (f_aldm + f_alfm) / 2
a_lm = (f_alwm + f_alfm) / 2
f_alm = float(f'{alm:.1f}')
f_a_lm = float(f'{a_lm:.1f}')
print(f'alm精馏段液相={f_alm},a‘lm提馏段液相={f_a_lm}')


def f3(a):
    b = float(f'{a:3f}')
    return b


# 塔径的计算
vs = f_V * f_Mvm / (3600 * f_pvm)
ls = f_L * f_Mlm / (3600 * f_plm)
f_vs = f3(vs)
f_ls = f3(ls)
print(f'精馏段的气、液相体积流量为\nvs={f_vs},ls={f_ls}')
v_s = f_V_ * f_M_vm / (3600 * f_p_vm)
l_s = f_L_ * f_M_lm / (3600 * f_p_lm)
f_v_s = f3(v_s)
f_l_s = f3(l_s)
print(f'提馏段的气、液相体积流量为\nv‘s={f_v_s},l’s={f_l_s}')
jsmith = ((f_ls * 3600) / (f_vs * 3600)) * ((f_plm / f_pvm) ** 0.5)
tsmith = ((f_l_s * 3600) / (f_v_s * 3600)) * ((f_p_lm / f_p_vm) ** 0.5)
f_jsmith = f3(jsmith)
f_tsmith = f3(tsmith)
print(f'smith数据\njsmith={f_jsmith},tsmith={f_tsmith}')


def init5():
    print("查阅下册p161的smith关系图")
    jc20 = input("精馏c20")  # 55000
    if jc20 == "":
        print("\033[0;31;40m使用调试数据0.082，0.081\033[0m")
        jc20 = 0.082
        tc20 = 0.081
    else:
        tc20 = input("提馏c20")
    return jc20, tc20


number_list = init5()
jc = float(number_list[0]) * ((f_alm / 20) ** 0.2)
tc = float(number_list[1]) * ((f_a_lm / 20) ** 0.2)
f_jc = f3(jc)
f_tc = f3(tc)
print(f'c数据\njc={f_jc},tc={f_tc}')
jumax = f_jc * (((f_plm - f_pvm) / f_pvm) ** 0.5)
tumax = f_tc * (((f_p_lm - f_p_vm) / f_p_vm) ** 0.5)
f_jumax = f3(jumax)
f_tumax = f3(tumax)
print(f'umax数据\njumax={f_jumax},tumax={f_tumax}')
ju = 0.6 * f_jumax
tu = 0.6 * f_tumax
f_ju = f3(ju)
f_tu = f3(tu)
print(f'u数据\nju={f_ju},tu={f_tu}')
jd = ((4 * f_vs) / (3.14 * f_ju)) ** 0.5
td = ((4 * f_v_s) / (3.14 * f_tu)) ** 0.5
f_jd = f3(jd)
f_td = f3(td)
print(f'D数据\njd={f_jd},td={f_td}')
if f_jd >= f_td:
    tand = float(f'{f_jd:.1f}') + 0.1
else:
    tand = float(f'{f_td:.1f}') + 0.1
print(f'按标准塔径圆整后\nD={tand}')
A = 3.14 / 4 * (tand ** 2)
print(f'塔截面积\nA={A}')
sju = f_vs / A
stu = f_v_s / A
f_sju = f3(sju)
f_stu = f3(stu)
print(f'u数据\n实际u={f_sju},实际u_={f_stu}')
# 精馏塔有效高度计算
ans = (f_Nt + f_Nj) // 6 + 1
zj = (f_Nj - 1) * 0.45
zt = (f_Nt - 1) * 0.45
z = ((f_Nt + f_Nj) - ans - 1) * 0.45 + (ans + 1) * 0.8
print(f'Z精={zj},Z提={zt}\n开{ans}个人孔，一个进料板，高度均为0.8m，\n有效高度为z{z}')
# 溢流堰长Lw
f_lw = f3(0.65 * tand)
print(f'溢流堰长Lw={f_lw}')
# ·出口堰高hw
f_how = f3(0.00284 * ((f_ls * 3600 / f_lw) ** (2 / 3)))
f_h_ow = f3(0.00284 * ((f_l_s * 3600 / f_lw) ** (2 / 3)))
print(f'how数据\nhow={f_how},h‘ow={f_h_ow}')
f_hw = f3(0.05 - f_how)
f_h_w = f3(0.08 - f_h_ow)
print(f'出口堰高hw\nhw={f_hw},h’w={f_h_w}')


# 降液管的宽度Wd和降液管的面积Af


def init6():
    print("查阅下册p166的弓形降液管关系图")
    WdD = input("Wd/D=")  # 55000
    if WdD == "":
        print("\033[0;31;40m使用调试数据0.082，0.081\033[0m")
        WdD = 0.124
        AfAt = 0.09
    else:
        AfAt = input("Af/At=")
    return WdD, AfAt


number_list = init6()
f_wd = f3(float(number_list[0]) * tand)
f_af = f3(float(number_list[1]) * A)
print(f'wd={f_wd},af={f_af}')
# 液体在降液管内停留的时间为
f_seita = f3(f_af * 0.45 * 3600 / 3600 / f_ls)
f_seita_ = f3(f_af * 0.45 * 3600 / 3600 / f_l_s)
print(f'液体在降液管内停留的时间为Θ={f_seita},Θ‘={f_seita_}')
# 降液管的底隙高度ho
f_ho = f3(f_ls / f_lw / 0.12)
f_h_o = f3(f_l_s / f_lw / 0.12)
print(f'降液管的底隙高度ho={f_ho},ho‘={f_h_o}')
m = f_hw - f_ho
m_ = f_h_w - f_h_o
print(f'hw-ho={m},h’w-h‘o={m_}')
# 塔板布置
f_uo = f3(10 / (f_pvm ** 0.5))
f_u_o = f3(10 / (f_p_vm ** 0.5))
print(f'本塔采用分块式塔板uo={f_uo},u’o={f_u_o}')
f_n = f3(4 * f_vs / 3.14 / (0.039 ** 2) / f_uo)
f_n = float(f'{f_n:.0f}') + 1
f_n_ = f3(4 * f_v_s / 3.14 / (0.039 ** 2) / f_u_o)
f_n_ = float(f'{f_n_:.0f}') + 1
print(f'n={f_n},n’={f_n_}')
# 开孔区面积Aa
X = f3(tand / 2 - (f_wd + 0.075))
R = f3(tand / 2 - 0.06)
Aa = f3(2 * (X * ((R ** 2 - X ** 2) ** 0.5) + math.pi / 180 * (R ** 2) * math.degrees(math.asin(X / R))))
print(f'X={X},R={R},Aa={Aa}')
# 排间距
f_jt_ = Aa / (f_n * 0.075)
f_tt_ = Aa / (f_n_ * 0.075)
print(f'精馏段={f_jt_},提馏段’={f_tt_}')
t_list = [65, 80, 100]
if f_jt_ <= f_tt_:
    any = f_tt_
    number = 0  # 按照精馏段计算
else:
    any = f_jt_
    number = 1  # 按照提馏段计算
t_ = 0
for i in range(len(t_list) - 1):
    if t_list[i] < any * 1000 < t_list[i + 1]:
        t_ = t_list[i]
    else:
        if any * 1000 > 100:
            t_ = 100
        else:
            t_ = 65
print(f'故取t’={t_}')


def init7():
    print("画图计算圈数")
    n = input("阀数N=")  #
    if n == "":
        print("\033[0;31;40m使用调试数据76\033[0m")
        n = 76
    return n


n = float(init7())
importn = n
if number == 0:
    any = "精馏段"
else:
    any = "提馏段"
print(f"\033[0;31;40m请注意,下面的数据都将以{any}的数据计算\033[0m")
if number == 0:
    # 精馏
    dict = {
        'vs': f_vs, 'pvm': f_pvm, 'su': f_sju, 'plm': f_plm, 'how': f_how,
        'hw': f_hw, 'ls': f_ls, 'ho': f_ho, 'hl': 0.05
    }
else:
    # 提馏
    dict = {
        'vs': f_v_s, 'pvm': f_p_vm, 'su': f_stu, 'plm': f_p_lm, 'how': f_h_ow,
        'hw': f_h_w, 'ls': f_l_s, 'ho': f_h_o, 'hl': 0.08
    }

importdict = {}


def qvfen(dict, n, lw, tand, wd, at, af):
    uo = f3(4 * dict["vs"] / 3.14 / (0.039 ** 2) / n)
    f = f3(uo * (dict["pvm"] ** 0.5))
    print(f'Uo={uo},f={f}')
    kkl = f3(dict["su"] / uo)
    print(f'开孔率={kkl}')
    # 临界孔速
    f_uoc = f3((73 / dict["pvm"]) ** (1 / 1.825))
    if f_uoc < uo:
        f_hc = f3(5.34 * (dict["pvm"] * (uo ** 2) / 2 / 9.81 / dict["plm"]))
        print("uo>uoc")
    else:
        f_hc = f3(19.9 * (uo ** 0.175 / dict["plm"]))
        print("uc<uoc")
    print(f'临界孔速uoc={f_uoc}，hc={f_hc}')
    f_h1 = f3((dict["hw"] + dict["how"]) * 0.5)
    print(f'液层的压降h1={f_h1}')
    f_hp = f3(f_hc + f_h1)
    f_pp = f3(dict["plm"] * 9.81 * f_hp)
    print(f'气体通过筛板的压降（单板压降）hf和△Pf\nhp={f_hp},△Pp={f_pp}')
    # 降液管液泛校核
    f_hd = f3(0.153 * ((dict["ls"] / lw / dict["ho"]) ** 2))
    print(f"\033[0;33;40m请注意,下面的数据是计算提馏段与精馏段的Hd\033[0m")
    if dict["hl"] == 0.05:
        f_Hd = f3(f_hc + f_h1 + dict["hl"] + f_hd)
        f_H_d = f3(f_hc + f_h1 + 0.08 + f_hd)
    else:
        f_Hd = f3(f_hc + f_h1 + 0.05 + f_hd)
        f_H_d = f3(f_hc + f_h1 + dict["hl"] + f_hd)
    print(f'降液管液泛校核\nhd={f_hd},Hd={f_Hd},H‘d={f_H_d}')
    fai = f3((0.45 + dict["hw"]) * 0.5)
    faidict = {"fai": fai}
    importdict.update(faidict)
    print(f'Ф={fai}')
    # 雾沫夹带量校核
    f_Z = f3(tand - 2 * wd)
    f_ab = f3(at - 2 * af)
    print(f"板上液流长度\nZ={f_Z},ab={f_ab}")
    print(f"\033[0;33;40m请注意,下面的数据是计算提馏段与精馏段的泛点率\033[0m")
    dict2 = {"z": f_Z, "ab": f_ab}
    dict.update(dict2)


qvfen(dict, n, f_lw, tand, f_wd, A, f_af)


def init10():
    jcf = input("Cf=")  #
    if jcf == "":
        print("\033[0;31;40m使用调试数据0.13,0.13\033[0m")
        jcf = 0.13
        tcf = 0.13
    else:
        tcf = input("C'f=")
    return jcf, tcf


numberlist = init10()
F1 = f3((100 * f_vs * ((f_pvm / (f_plm - f_pvm)) ** 0.5) + 1.36 * f_ls * dict["z"]) / (
        dict["ab"] * float(numberlist[0]) * 1))
F1_ = f3(
    (100 * f_v_s * ((f_p_vm / (f_p_lm - f_p_vm)) ** 0.5) + 1.36 * f_l_s * dict["z"]) / (
            dict["ab"] * float(numberlist[1]) * 1))
print(f"第一组公式\nF={F1},F'={F1_}")
F2 = f3((100 * f_vs * ((f_pvm / (f_plm - f_pvm)) ** 0.5)) / (0.78 * 1 * A * float(numberlist[0])))
F2_ = f3((100 * f_vs * ((f_pvm / (f_plm - f_pvm)) ** 0.5)) / (0.78 * 1 * A * float(numberlist[1])))
print(f"第二组公式\nF={F2},F'={F2_}")
# 性能图
jk1 = f3(((0.8 * 1 * float(numberlist[0]) * dict["ab"]) ** 2 / (f_pvm / (f_plm - f_pvm))) ** 0.5)
jb1 = f3(((1.36 * dict["z"]) ** 2 / (f_pvm / (f_plm - f_pvm))) ** 0.5)
print(f"Vs={jk1}-{jb1}Ls")


def check_json(any):
    with open("../date.json", 'r') as read:
        file = read.read()
        if len(file) > 0:
            result = json.loads(file)
            if any not in result:
                return "yes"
            else:
                return "no"
        else:
            data = {}
            with open("../date.json", 'w') as write:
                json.dump(data, write)


def read_json(any):
    with open("../date.json", 'r') as read:
        file = read.read()
        result = json.loads(file)
        return result[any]


def write(anyfloat, anyname):
    anyfloat.sort()
    dict = {anyname: anyfloat}
    print(f"{anyname}{anyfloat}")
    with open("../date.json", 'r') as read:
        result = json.load(read)
        result.update(dict)
    with open("../date.json", 'w') as write:
        json.dump(result, write)


tk1 = f3(((0.8 * 1 * float(numberlist[1]) * dict["ab"]) ** 2 / (f_p_vm / (f_p_lm - f_p_vm))) ** 0.5)
tb1 = f3(((1.36 * dict["z"]) ** 2 / (f_p_vm / (f_p_lm - f_p_vm))) ** 0.5)
print(f"V's={tk1}-{tb1}L's")
# （2）降液管液泛线（气相负荷上限线）
ja = f3(0.4416 * f_pvm / f_plm / (0.039 ** 4) / (importn ** 2))
jb = f3(0.153 / (f_lw ** 2) / (f_ho ** 2))
jc = f3(1.5 * f_hw - 0.5 * (0.45 + f_hw))
jd = f3(1 / (f_lw ** (2 / 3)))
print(f"如果出现--情况即为+号")
print(f"vs2=-{f3(jc / ja)}-{f3(jb / ja)}ls2-{f3(jd / ja)}ls2/3")
ta = f3(0.4416 * f_p_vm / f_p_lm / (0.039 ** 4) / (importn ** 2))
tb = f3(0.153 / (f_lw ** 2) / (f_h_o ** 2))
tc = f3(1.5 * f_h_w - 0.5 * (0.45 + f_h_w))
td = f3(1 / (f_lw ** (2 / 3)))
print(f"如果出现--情况即为+号")
print(f"vs2=-{f3(tc / ta)}-{f3(tb / ta)}ls2-{f3(td / ta)}ls2/3")
# 漏液线（气相负荷下限线)
jv = f3(3.14 / 4 * 0.039 * 0.039 * importn * 5 / (f_pvm ** 0.5))
tv = f3(3.14 / 4 * 0.039 * 0.039 * importn * 5 / (f_p_vm ** 0.5))
print(f"Vs={jv}V's={tv}")
lsmin = f3(((0.006 * 1000 / 2.84) ** (2 / 3)) * f_lw / 3600)
print(f"lsmin={lsmin}")
lsmax = f3(f_af * 0.45 / 5)
print(f"lsmax={lsmax}")
jyqb = f3(f_vs / f_ls)
tyqb = f3(f_v_s / f_l_s)
print(f"精馏段气液比={jyqb}提馏段气液比={tyqb}")
print("表格数据一览，不满意的可以重新取")
ls_str = []
ls_float = []
vs_float = []
ans = check_json("ls")
if ans == "yes":
    while 1:
        number = f3(random.uniform(lsmin, lsmax * 1.2))
        if number not in ls_str:
            ls_str.append(str(number))
        else:
            continue
        if len(ls_str) == 30:
            write(ls_str, "ls")
            break
ans = check_json("v_s")
if ans == "yes":
    for i in range(30):
        number = float(read_json("ls")[i])
        ls_float.append(number)
        ans = tk1 - tb1 * number
        vs_float.append(ans)
    write(vs_float, "v_s")
vs_float = []
ans = check_json("vs")
if ans == "yes":
    for i in range(30):
        number = float(read_json("ls")[i])
        ls_float.append(number)
        ans = jk1 - jb1 * number
        vs_float.append(ans)
    write(vs_float, "vs")
vs_float = []
ans = check_json("2vs")
if ans == "yes":
    for i in range(30):
        number = float(read_json("ls")[i])
        ans = -f3(jc / ja) - f3(jb / ja) * (number ** 2) - f3(jd / ja) * (number ** (2 / 3))
        vs_float.append(ans)
    write(vs_float, "2vs")
vs_float = []
ans = check_json("2v_s")
if ans == "yes":
    for i in range(30):
        number = float(read_json("ls")[i])
        ans = -f3(tc / ja) - f3(tb / ta) * (number ** 2) - f3(td / ta) * (number ** (2 / 3))
        vs_float.append(ans)
    write(vs_float, "2v_s")
print("表格数据抄下面的，自动取得，30个点取8个")
list = [1, 5, 10, 15, 18, 20, 25, 29]
for i in range(8):
    ls = float(read_json("ls")[list[i]])
    vs = float(read_json("vs")[-list[i]])
    v_s = float(read_json("v_s")[-list[i]])
    vs2 = float(read_json("2vs")[-list[i]])
    v_s2 = float(read_json("2v_s")[-list[i]])
    print(f"第{i + 1}列ls={ls},vs={vs},v's={v_s},vs={vs2},v's={v_s2}")
print("操作弹性计算自己算")


def draw5(vsmin, lsmin, lsmax, jyqb):
    vs1 = float(read_json("vs")[-1])
    vs2 = float(read_json("2vs")[-1])
    if vs1 > vs2:
        ans = vs1
    else:
        ans = vs2
    plt.xlim((0, lsmax * 1.7))
    plt.ylim((0, ans * 1.1))
    ls = read_json("ls")
    vs = read_json("vs")
    vs2 = read_json("2vs")
    for i in range(len(ls) - 1):
        xtx = [float(ls[i]), float(ls[i + 1])]
        xty = [vs[-1 - i], vs[-2 - i]]
        plt.plot(xtx, xty, 'r-')
        xtx1 = [float(ls[i]), float(ls[i + 1])]
        xty1 = [vs2[-1 - i], vs2[-2 - i]]
        plt.plot(xtx1, xty1, 'y-')
    plt.plot(xtx, xty, c='r', lw=3., label='过量雾沫夹带线')
    plt.plot(xtx1, xty1, c='y', lw=3., label='降液管液泛线')
    # vsmin
    xtx = [0, lsmax * 1.1]
    xty = [vsmin, vsmin]
    plt.plot(xtx, xty, 'b-')
    plt.plot(xtx, xty, c='b', lw=3., label='漏液线')
    # lsmin
    xtx = [lsmin, lsmin]
    xty = [0, ans]
    plt.plot(xtx, xty, 'c-')
    plt.plot(xtx, xty, c='c', lw=3., label='液相负荷下限线')
    # lsmax
    xtx = [lsmax, lsmax]
    xty = [0, ans]
    plt.plot(xtx, xty, 'm-')
    plt.plot(xtx, xty, c='m', lw=3., label='液相负荷上限线')
    # p点
    ans = jyqb * lsmin
    xtx = [0, lsmin]
    xty = [0, ans]
    plt.plot(xtx, xty, 'k-')
    plt.plot(xtx1, xty1, c='k', lw=3., label='操作线')
    ans2 = jyqb * lsmax * 0.9
    xtx = [lsmin, lsmin + lsmax * 0.9]
    xty = [ans, ans2]
    plt.plot(xtx, xty, 'k-')
    print(f"交点坐标({lsmin},{ans})")
    plt.rcParams['font.sans-serif'] = ['SimSun']
    plt.title('Load performance graph')
    plt.legend()
    plt.show()


draw5(jv, lsmin, lsmax, jyqb)
def draw6(vsmin, lsmin, lsmax, jyqb):
    vs1 = float(read_json("v_s")[-1])
    vs2 = float(read_json("2v_s")[-1])
    if vs1 > vs2:
        ans = vs1
    else:
        ans = vs2
    plt.xlim((0, lsmax * 1.7))
    plt.ylim((0, ans * 1.1))
    ls = read_json("ls")
    vs = read_json("v_s")
    vs2 = read_json("2v_s")
    for i in range(len(ls) - 1):
        xtx = [float(ls[i]), float(ls[i + 1])]
        xty = [vs[-1 - i], vs[-2 - i]]
        plt.plot(xtx, xty, 'r-')
        xtx1 = [float(ls[i]), float(ls[i + 1])]
        xty1 = [vs2[-1 - i], vs2[-2 - i]]
        plt.plot(xtx1, xty1, 'y-')
    plt.plot(xtx, xty, c='r', lw=3., label='过量雾沫夹带线')
    plt.plot(xtx1, xty1, c='y', lw=3., label='降液管液泛线')
    # vsmin
    xtx = [0, lsmax * 1.1]
    xty = [vsmin, vsmin]
    plt.plot(xtx, xty, 'b-')
    plt.plot(xtx, xty, c='b', lw=3., label='漏液线')
    # lsmin
    xtx = [lsmin, lsmin]
    xty = [0, ans]
    plt.plot(xtx, xty, 'c-')
    plt.plot(xtx, xty, c='c', lw=3., label='液相负荷下限线')
    # lsmax
    xtx = [lsmax, lsmax]
    xty = [0, ans]
    plt.plot(xtx, xty, 'm-')
    plt.plot(xtx, xty, c='m', lw=3., label='液相负荷上限线')
    # p点
    ans = jyqb * lsmin
    xtx = [0, lsmin]
    xty = [0, ans]
    plt.plot(xtx, xty, 'k-')
    plt.plot(xtx1, xty1, c='k', lw=3., label='操作线')
    ans2 = jyqb * lsmax * 0.9
    xtx = [lsmin, lsmin + lsmax * 0.9]
    xty = [ans, ans2]
    plt.plot(xtx, xty, 'k-')
    print(f"交点坐标({lsmin},{ans})")
    plt.rcParams['font.sans-serif'] = ['SimSun']
    plt.title('Load performance graph')
    plt.legend()
    plt.show()


draw6(tv, lsmin, lsmax, tyqb)
# 1、塔顶冷凝器的试算与初选
sanjiaot1 = f3(importtd - 30)
sanjiaot2 = f3(importtd - 40)
sanjiaotm = f3((sanjiaot1 - sanjiaot2) / (math.log(sanjiaot1 / sanjiaot2, math.e)))
print(f"△t1={sanjiaot1}△t2={sanjiaot2}△tm={sanjiaotm}")


def rrr(any):
    any = float(any)
    list = [[80, 394.1, 379.9], [90, 386.9, 373.8], [100, 379.3, 367.6], [110, 371.5, 361.2], [120, 363.2, 354.6]]
    for i in range(len(list) - 1):
        if list[i + 1][0] + 10 > any > list[i][0]:
            befor_info = list[i]
            after_info = list[i + 1]
            pa = befor_info[1] + ((after_info[1] - befor_info[1]) / 10) * (td - befor_info[0])
            pb = befor_info[2] + ((after_info[2] - befor_info[2]) / 10) * (td - befor_info[0])
    pa = float(f'{pa:.1f}')
    pb = float(f'{pb:.1f}')
    return pa, pb


rrrlist = rrr(importtd)
print(f"γ苯={rrrlist[0]}γ甲苯={rrrlist[1]}")
r = f3(FloatDistillate * rrrlist[0] + (1 - FloatDistillate) * rrrlist[1])
q = f3(r * f_V * f_Mvm / 3600)
print(f"γ={r}Q={q}")
s = f3(q * 1000 / 550 / sanjiaotm)
print(f"S={s}")
print(f"查化工原理上册书p361页内导流浮头式")


def init12():
    n = input("管根数=")  #
    if n == "":
        print("\033[0;31;40m使用调试数据164,2,3000\033[0m")
        n = 164
        fai = 2
        l = 3000
    else:
        fai = input("外径=mm")
        l = input("管长mm")
    return n, fai, l


huanreqilist = init12()
s = f3(float(huanreqilist[0]) * 3.14 * float(huanreqilist[1]) * 0.001 * (float(huanreqilist[2]) / 1000 - 0.1))
k = f3(q * 1000 / s / sanjiaotm)
print(f"S0={s}K0={k}")


def d(vs, u):
    ans = f3((4 * vs / 3.14 / u) ** 0.5)
    return ans


def u(vs, d):
    ans = f3(4 * vs / 3.14 / (d ** 2))
    return ans


def init11():
    n = input("直径=mm")  #
    if n == "":
        print("\033[0;31;40m使用调试数据100\033[0m")
        n = 100
    return float(n) / 1000

print("查书化工原理上次350页")
v = f3(f_F * f_Mf)
vs = f3(v / f_plm / 3600)
jld = d(vs, 0.5)
print(f"V={v}Vs={vs}D={jld}")
su = u(vs, init11())
print(f"u={su}")
hld = d(f_ls, 0.5)
print(f"D={hld}")
su = u(f_ls, init11())
print(f"u={su}")
v = f3(f_W * f_Mf)
vs = f3(v / f_p_lm / 3600)
jld = d(vs, 0.5)
print(f"V={v}Vs={vs}D={jld}")
su = u(vs, init11())
print(f"u={su}")
jld = d(f_v_s, 20)
print(f"D={jld}")
su = u(f_v_s, init11())
print(f"u={su}")
jld = d(f_vs, 20)
print(f"D={jld}")
su = u(f_vs, init11())
print(f"u={su}")
os.system("pause")