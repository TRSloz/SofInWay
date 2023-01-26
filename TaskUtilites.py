import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP


def get_psat(EOS, T, mf):
    """
    Расчёт давления жидкой и паровой фаз при заданном составе
    :param EOS: CoolProp.AbstractState;
    :param T: Температура раствора, (K);
    :param mf: [x1, x2] массовые концентрации компонентов раствора (kg/kg);
    :return: [T (K), x1 (kg/kg), PsatLiquid (bar), PsatVapor (bar)]
    """

    # Set the mass fractions of the mixture
    EOS.set_mass_fractions(mf)
    # Liquid saturated pressure
    EOS.update(CP.QT_INPUTS, 0, T)
    PsatL = EOS.saturated_vapor_keyed_output(CP.iP)
    # Vapor saturated pressure
    EOS.update(CP.QT_INPUTS, 1, T)
    PsatV = EOS.saturated_vapor_keyed_output(CP.iP)

    return np.array([[mf[0], PsatL*1e-5, PsatV*1e-5]])


def get_tsat(EOS, P, mf):
    """
    Расчёт температуры жидкой и паровой фаз при заданном составе
    :param EOS: CoolProp.AbstractState;
    :param P: Давление раствора, (Pa);
    :param mf: [x1, x2] массовые концентрации компонентов раствора (kg/kg);
    :return: [x1 (kg/kg), TsatLiquid (C), TsatVapor (C)]
    """

    # Set the mass fractions of the mixture
    EOS.set_mass_fractions(mf)
    # Liquid saturated pressure
    EOS.update(CP.PQ_INPUTS, P, 0)
    TsatL = EOS.T()
    # Vapor saturated pressure
    EOS.update(CP.PQ_INPUTS, P, 1)
    TsatV = EOS.T()

    return np.array([[mf[0], TsatL-273.15, TsatV-273.15]])


def get_py(EOS, T, Lmf):
    """
    Расчёт давления насыщения и состава паровой фазы
    :param EOS: CoolProp.AbstractState;
    :param T: Температура раствора (K);
    :param Lmf: [x1, x2] массовые концентрации компонентов раствора в жидкой фазе, (kg/kg);
    :return: [y1 (kg/kg), Psat (bar)]
    """

    # Set the mass fractions of the mixture
    EOS.set_mass_fractions(Lmf)
    # Saturated pressure
    EOS.update(CP.QT_INPUTS, 0, T)
    #Psat = EOS.saturated_vapor_keyed_output(CP.iP)
    Psat = EOS.p()
    # Vapor phase mass fractions
    EOS.set_mole_fractions(EOS.mole_fractions_vapor())
    Vmf = EOS.get_mass_fractions()

    return Vmf[0], Psat*1e-5


def get_conductivity(Sub1, Sub2, T, Lmf):
    """
    Расчёт теплопроводности жидкой фазы раствора
    По данным о чистых компонентах. Модель прогнозирования Ли [Рид 1982]
    :param Sub1: Нименование компонента 1;
    :param Sub2: Нименование компонента 2;
    :param T: Температура раствора (C);
    :param Lmf: [x1, x2] массовые концентрации компонентов раствора в жидкой фазе, (kg/kg);
    :return:  thermal conductivity (W/m/K)
    """
    Tsat=T+273.15
    Mix = Sub1 + '&' + Sub2
    REFPROP = CP.AbstractState("REFPROP", Mix)
    # Задаём массовую концентрацию компонентов в жидкой фазе
    REFPROP.set_mass_fractions(Lmf)
    REFPROP.update(CP.QT_INPUTS, 0, Tsat)
    # Соответствующая мольная концентрация компонентов в жидкой фазе
    xmole = REFPROP.mole_fractions_liquid()
    # Свойства чистых компонентов
    # - вода
    REFPROP_W = CP.AbstractState("REFPROP", "Water")
    REFPROP_W.update(CP.QT_INPUTS, 0, Tsat)
    V1    = 1/REFPROP_W.rhomolar()
    Lmbd1 = REFPROP_W.conductivity()
    # - этиленгликоль
    REFPROP_EG = CP.AbstractState("REFPROP", "EthyleneGlycol")
    REFPROP_EG.update(CP.QT_INPUTS, 0, Tsat)
    V2    = 1/REFPROP_EG.rhomolar()
    Lmbd2 = REFPROP_EG.conductivity()
    # Объемные доли
    F1 = xmole[0]*V1 / (xmole[0]*V1 + xmole[1]*V2)
    F2 = xmole[1]*V2 / (xmole[0]*V1 + xmole[1]*V2)
    # Теплопроводность
    Lmbd12 = 2 / (1/Lmbd1+1/Lmbd2)
    Lmbd_mix = F1**2*Lmbd1 + 2*F1*F2*Lmbd12 + F2**2*Lmbd2
    return np.array([[Lmf[0], Lmbd_mix]])


def get_viscosity(Sub1, Sub2, T, Lmf):
    """
    Расчёт кинематической вязкости жидкой фазы раствора
    По данным о чистых компонентах. Модель прогнозирования Лобе [Рид 1982]
    :param Sub1: Нименование компонента 1;
    :param Sub2: Нименование компонента 2;
    :param T: Температура раствора (C);
    :param Lmf: [x1, x2] массовые концентрации компонентов раствора в жидкой фазе, (kg/kg);
    :return: кинематическая вязкость (cSt)
    """
    Tsat=T+273.15
    Mix = Sub1 + '&' + Sub2
    REFPROP = CP.AbstractState("REFPROP", Mix)

    # Задаём массовую концентрацию компонентов в жидкой фазе
    REFPROP.set_mass_fractions(Lmf)
    REFPROP.update(CP.QT_INPUTS, 0, Tsat)
    # Соответствующая мольная концентрация компонентов в жидкой фазе
    xmole = REFPROP.mole_fractions_liquid()

    # Свойства чистых компонентов
    # - вода
    REFPROP_W = CP.AbstractState("REFPROP", "Water")
    REFPROP_W.update(CP.QT_INPUTS, 0, Tsat)
    V1 = 1/REFPROP_W.rhomolar()
    Nu1 = REFPROP_W.viscosity()/REFPROP_W.rhomass()
    # - этиленгликоль
    REFPROP_EG = CP.AbstractState("REFPROP", "EthyleneGlycol")
    REFPROP_EG.update(CP.QT_INPUTS, 0, Tsat)
    V2 = 1/REFPROP_EG.rhomolar()
    Nu2 = REFPROP_EG.viscosity()/REFPROP_EG.rhomass()

    # Объемные доли
    F1 = xmole[0]*V1 / (xmole[0]*V1 + xmole[1]*V2)
    F2 = xmole[1]*V2 / (xmole[0]*V1 + xmole[1]*V2)

    # Характеристический параметр
    LnAB = np.log(Nu2/Nu1)
    Alpha1 = -1.7*LnAB
    Alpha2 = 0.27*LnAB + np.sqrt(1.3*LnAB)

    # Кинематическая вязкость
    Nu_mix = F1*Nu1*np.exp(F2*Alpha2) + F2*Nu2*np.exp(F1*Alpha1)
    return np.array([[Lmf[0], Nu_mix*1e6]])


def plotdata(Tmix, x1, y1, Psat, dt):
    """
    График результатов (Изотерма)
    """

    fig, ax = plt.subplots(figsize=(7, 5), dpi=100)

    ax.plot(dt[:,0], dt[:, 1], color='r', linewidth = 2) # Liquid saturation line
    ax.plot(dt[:,0], dt[:, 2], color='b', linewidth = 2) # Vapor saturation line
    ax.fill_between(dt[:, 0], dt[:, 1], dt[:, 2],
                    facecolor='deepskyblue',
                    alpha=0.2)

    # Result points
    ax.plot(x1, Psat, 'o', color='limegreen', markeredgecolor='0')
    ax.plot(y1, Psat, 'o', color='limegreen', markeredgecolor='0')

    ax.set_xlim([0, 1])
    ax.set_xlabel('Mass fraction of Propane, kg/kg', fontstyle='italic')
    ax.set_ylabel('Saturated pressure, bar', fontstyle='italic')
    ax.set_title('Propane/n-Butane', fontweight='black')
    ax.text(0.05,9.5, 'T = '+str(Tmix)+'K', fontweight = 'bold' )

    # X liquid arrow
    ax.annotate(
        'X propane',
        xy=(x1, Psat), xycoords='data',
        xytext=(x1, 1), textcoords='data',
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="->", shrinkB=3, edgecolor='forestgreen'),
        horizontalalignment='center',
        verticalalignment='top',
    )

    # Y liquid arrow
    ax.annotate(
        'Y propane',
        xy=(y1, Psat), xycoords='data',
        xytext=(y1, 1.5), textcoords='data',
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="<-", shrinkB=3, edgecolor='forestgreen'),
        horizontalalignment='center',
        verticalalignment='top',
    )

    # Psat arrows
    ax.annotate(
        '',
        xy=(x1, Psat), xycoords='data',
        xytext=(y1, Psat), textcoords='data',
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="<-", shrinkA=0, shrinkB=3, edgecolor='black'),
        horizontalalignment='center',
        verticalalignment='top',
    )
    ax.annotate(
        'P sat',
        xy=(x1, Psat), xycoords='data',
        xytext=(-0.07, Psat), textcoords='data',
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="<-", shrinkA=0, shrinkB=3, edgecolor='black'),
        horizontalalignment='center',
        verticalalignment='top',
    )

    plt.show()


def plotdata2(Pmix, x1, y1, Tsat, dt):
    """
    График результатов (Изобара)
    """

    fig, ax = plt.subplots(figsize=(7, 5), dpi=100)

    ax.plot(dt[:,0], dt[:, 1], color='r', linewidth=2) # Liquid saturation line
    ax.plot(dt[:,0], dt[:, 2], color='b', linewidth=2) # Vapor saturation line
    ax.fill_between(dt[:, 0], dt[:, 1], dt[:, 2],
                    facecolor='lightgreen',
                    alpha=0.2)

    # Result points
    ax.plot(x1, Tsat, 'o', color='deepskyblue', markeredgecolor='0')
    ax.plot(y1, Tsat, 'o', color='deepskyblue', markeredgecolor='0')

    ax.set_xlim([0, 1])
    ax.set_xlabel('Mass fraction of Water, kg/kg', fontstyle='italic')
    ax.set_ylabel('Saturated temperature, C', fontstyle='italic')
    ax.set_title('Water/EtylenGlycole', fontweight='black')
    ax.text(0.7,220, 'P = '+str(np.round(Pmix*1e-5, 2))+' bar', fontweight = 'bold' )

    # X liquid arrow
    ax.annotate(
        'X water',
        xy=(x1, Tsat), xycoords='data',
        xytext=(x1, 113), textcoords='data',
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="->", shrinkA=4, shrinkB=3, edgecolor='forestgreen'),
        horizontalalignment='center',
        verticalalignment='top',
    )

    # Y liquid arrow
    ax.annotate(
        'Y water',
        xy=(y1, Tsat), xycoords='data',
        xytext=(y1, 107), textcoords='data',
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="<-", shrinkA=16, shrinkB=3, edgecolor='forestgreen'),
        horizontalalignment='center',
        verticalalignment='top',
    )

    # Tsat arrows
    ax.annotate(
        '',
        xy=(x1, Tsat), xycoords='data',
        xytext=(y1, Tsat), textcoords='data',
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="<-", shrinkA=0, shrinkB=3, edgecolor='black'),
        horizontalalignment='center',
        verticalalignment='top',
    )
    ax.annotate(
        'T sat',
        xy=(x1, Tsat), xycoords='data',
        xytext=(-0.11, Tsat), textcoords='data',
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="<-", shrinkA=26, shrinkB=3, edgecolor='black'),
        horizontalalignment='center',
        verticalalignment='center',
    )

    plt.show()


def plotdata3(Pmix, x1, y1, Tsat, dt):
    """
    График результатов (Изобара)
    """

    fig, ax = plt.subplots(figsize=(7, 5), dpi=100)

    ax.plot(dt[:,0], dt[:, 1], color='r', linewidth=2) # Liquid saturation line
    ax.plot(dt[:,0], dt[:, 2], color='b', linewidth=2) # Vapor saturation line
    ax.fill_between(dt[:, 0], dt[:, 1], dt[:, 2],
                    facecolor='lightgreen',
                    alpha=0.2)

    # Result points
    ax.plot(x1, Tsat, 'o', color='deepskyblue', markeredgecolor='0')
    ax.plot(y1, Tsat, 'o', color='deepskyblue', markeredgecolor='0')

    ax.set_xlim([0, 1])
    ax.set_xlabel('Mass fraction of Water, kg/kg', fontstyle='italic')
    ax.set_ylabel('Saturated temperature, C', fontstyle='italic')
    ax.set_title('Water/EtylenGlycole', fontweight='black')
    ax.text(0.7,220, 'P = '+str(np.round(Pmix*1e-5,2))+' bar', fontweight='bold' )

    # X liquid arrow
    ax.annotate(
        'X water',
        xy=(x1, Tsat), xycoords='data',
        xytext=(x1, 107), textcoords='data',
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="<-", shrinkA=16, shrinkB=3, edgecolor='forestgreen'),
        horizontalalignment='center',
        verticalalignment='top',
    )

    # Y liquid arrow
    ax.annotate(
        'Y water',
        xy=(y1, Tsat), xycoords='data',
        xytext=(y1, 113), textcoords='data',
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="->", shrinkA=4, shrinkB=3, edgecolor='forestgreen'),
        horizontalalignment='center',
        verticalalignment='top',
    )

    # Tsat arrows
    ax.annotate(
        '',
        xy=(x1, Tsat), xycoords='data',
        xytext=(y1, Tsat), textcoords='data',
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="->", shrinkA=0, shrinkB=3, edgecolor='black'),
        horizontalalignment='center',
        verticalalignment='top',
    )
    ax.annotate(
        'T sat',
        xy=(x1, Tsat), xycoords='data',
        xytext=(-0.15, Tsat), textcoords='data',
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="<-", shrinkA=43, shrinkB=3, edgecolor='black'),
        horizontalalignment='center',
        verticalalignment='center',
    )

    plt.show()


def plotdata4(x1, Tsat, Lmbd_mix1, Lmbd_mix2, dt1, dt2):
    """
    График результатов (теплопроводность)
    """

    fig, ax = plt.subplots(figsize=(7, 5), dpi=100)

    ax.plot(dt1[:,0], dt1[:, 1], color='r', linewidth=1)
    ax.plot(dt2[:,0], dt2[:, 1], color='b', linewidth=1)

    # Result points
    ax.plot(x1[0], Lmbd_mix1, 'o', color='pink', markeredgecolor='0')
    ax.plot(x1[1], Lmbd_mix2, 'o', color='deepskyblue', markeredgecolor='0')

    ax.set_xlim([0, 1])
    ax.set_xlabel('Mass fraction of Water, kg/kg', fontstyle='italic')
    ax.set_ylabel('Thermal conductivity, W/m/K', fontstyle='italic')
    ax.set_title('Water/EtylenGlycole', fontweight='black')
    ax.text(0.1, 0.63, 'T = '+str(np.round(Tsat[0], 2))+' C', fontweight='bold', color='r')
    ax.text(0.7, 0.30, 'T = '+str(np.round(Tsat[1], 2))+' C', fontweight='bold', color='b')

    # arrow 1
    ax.annotate(
        'Task 1.1',
        xy=(x1[0], Lmbd_mix1), xycoords='data',
        xytext=(x1[0]-0.1, Lmbd_mix1+0.1), textcoords='data',
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="->", shrinkA=0, shrinkB=3, edgecolor='0.2'),
        horizontalalignment='center',
        verticalalignment='top',
    )

    # arrow 2
    ax.annotate(
        'Task 1.2',
        xy=(x1[1], Lmbd_mix2), xycoords='data',
        xytext=(x1[1]+0.1, Lmbd_mix2+0.15), textcoords='data',
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="->", shrinkA=0, shrinkB=3, edgecolor='0.2'),
        horizontalalignment='center',
        verticalalignment='top',
    )

    plt.show()


def plotdata5(x1, Tsat, Nu_mix1, Nu_mix2, dt1, dt2):
    """
    График результатов (кинематическая вязкость)
    """

    fig, ax = plt.subplots(figsize=(7, 5), dpi=100)

    ax.plot(dt1[:,0], dt1[:, 1], color='b', linewidth=1)
    ax.plot(dt2[:,0], dt2[:, 1], color='r', linewidth=1)

    # Result points
    ax.plot(x1[0], Nu_mix1, 'o', color='deepskyblue', markeredgecolor='0')
    ax.plot(x1[1], Nu_mix2, 'o', color='pink', markeredgecolor='0')

    ax.set_xlim([0, 1])
    ax.set_xlabel('Mass fraction of Water, kg/kg', fontstyle='italic')
    ax.set_ylabel('Kinematic viscosity, cSt', fontstyle='italic')
    ax.set_title('Water/EtylenGlycole', fontweight='black')
    ax.text(0.2, 1.0, 'T = '+str(np.round(Tsat[0], 2))+' C', fontweight='bold', color='b')
    ax.text(0.1, 0.2, 'T = '+str(np.round(Tsat[1], 2))+' C', fontweight='bold', color='r')

    # arrow 1
    ax.annotate(
        'Task 1.1',
        xy=(x1[0], Nu_mix1), xycoords='data',
        xytext=(x1[0], Nu_mix1+0.2), textcoords='data',
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="->", shrinkA=0, shrinkB=3, edgecolor='0.2'),
        horizontalalignment='center',
        verticalalignment='top',
    )

    # arrow 2
    ax.annotate(
        'Task 1.2',
        xy=(x1[1], Nu_mix2), xycoords='data',
        xytext=(x1[1]+0.05, Nu_mix2-0.15), textcoords='data',
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="->", shrinkA=0, shrinkB=3, edgecolor='0.2'),
        horizontalalignment='center',
        verticalalignment='top',
    )

    plt.show()