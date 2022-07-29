import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
import makecharts as mc
from copy import copy
import resistances as rs

#Исходные данные (все единицы в СИ, если не указано иное)
Q_in=1. #объемный расход на входе в систему, м3/с
P_in=150000 #давление на входе в систему, Па
d_temp=10. #подогрев воды в теплообменнике, С
t_out_target=20 #целевая темпеартура на выходе из системы в точке измерения, С
d=0.05 #диаметр всех труб в системе, м
S=np.pi*(d/2)**2 #площадь сечения труб
Rho=1000 #плотность воды, кг/м3
Nu=1.006e-6 #вязкость жидкости, м2/с
g=9.81
L1=10 #длина основной трубы (систему) через котрую течет основной поток с нагревом
L2=20 #длина "обходной" трубы с холодной водой, подмешиваемой в основной поток через регулируемый клапан
dH_friction1=rs.dH_pipe #функция для вычисления потерь трения в основной системе
dH_friction2=rs.dH_pipe #функция для вычисления потерь трения в обходной трубе
dH_valve=rs.dH_valve3 #функция для расчета потерь в регулируемом клапане

KP=1 #пропорциональный коэффициент регулятора
KI=1 #интегральный коэффициент регулятора

#временнЫе настройки расчета
T=100 #общее время расчета
dt=0.01 #шаг расчета по времени
t_in_by_time= {0:10,#изменение температуры воды (С) на входе в систему по времени
               40:30,
               60:15,
               80:20}

#функция для расчета площади проходного сечения регулируемого клапана в зависимости от величины управляющего сигнала: от 0 до 1 (0 - клапан закрыт, 1 - открыт)
def S_valve(t):
    return np.pi*(d/2)**2*t

#функция расчета управляющего сигнала ПИ-регулятора
def PI_controller(error, prev_error):
    return prev_error + KP*(error - prev_error) + KI*error

#первоначальные приближения
t_in=10 #температура на входе в систему
valve_position=0.5

#Массивы для записи результатов:
time=[t_in]

#вспомогательная функция для расчета потоков по всем трубам системы, здесь x=Q1/Q_in - относительный расход через трубу 1
def iteration(x, valve_position):
    #уравнение Бернулли для трубы 1 от точки рассоединения с трубой 2 до точки соединения с ней же: P_in/Rho/g +  V1**2/2/g = P_out/Rho/g +  V1**2/2/g  + dHfriction_pipe1
    #для трубы 2 от точки рассоединения с трубой 1 до точки соединения с ней же: P_in/Rho/g +  V2**2/2/g = P_out/Rho/g +  V2**2/2/g  + dHfriction_pipe2 + dH_valve
    V1=Q_in*x/S #находим в первом приближении скорость исходя из значения на предыдущем шаге
    dHfriction_pipe1 = dH_friction1(V1, d, Nu, L1, g) #потери давления на трение в основной системе
    V2 = Q_in * (1-x) / S #скорость потока в обходной трубе
    dHfriction_pipe2 = dH_friction2(V2, d, Nu, L2, g)  # потери давления на трение в обходной системе
    dHvalve = dH_valve(valve_position,V2,g)  # потери давления на трение в обходной системе
    P_out_1 = (P_in / Rho / g - dHfriction_pipe1) * Rho * g
    P_out_2 = (P_in / Rho / g - dHfriction_pipe2 - dHvalve) * Rho * g
    return P_out_2-P_out_1

#расчет
for i in range(1,int(T/dt)):
    t=i*dt
    time.append(t)
    if i%1000==0:
        print(f"step calculated: {i}")

    #как температура меняется на входе в систему
    for time_, temperature_ in t_in_by_time.items():
        if (abs(time_ - t) < dt * 0.5):
            t_in = temperature_

    #находим давление на выходе из системы
    res=root_scalar(iteration,x0=valve_position,x1=valve_position*1.0001,method='secant',args=(valve_position))
    P_out=res.root
    print(P_out)




