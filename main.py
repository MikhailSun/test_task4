import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
import makecharts as mc
from copy import copy
import resistances as rs

#Исходные данные (все единицы в СИ, если не указано иное)
Q_in=0.2 #объемный расход на входе в систему, м3/с
P_in=150000 #давление на входе в систему, Па
d_temp=10. #подогрев воды в теплообменнике, С
t_out_target=17.5 #целевая темпеартура на выходе из системы в точке измерения, С
d=0.1 #диаметр всех труб в системе, м
S=np.pi*(d/2)**2 #площадь сечения труб
Rho=1000 #плотность воды, кг/м3
Nu=1.006e-6 #вязкость жидкости, м2/с
g=9.81
L1=10 #длина основной трубы (систему) через котрую течет основной поток с нагревом
L2=20 #длина "обходной" трубы с холодной водой, подмешиваемой в основной поток через регулируемый клапан
dH_friction1=rs.dH_pipe #функция для вычисления потерь трения в основной системе
dH_friction2=rs.dH_pipe #функция для вычисления потерь трения в обходной трубе
dH_valve=rs.dH_valve3 #функция для расчета потерь в регулируемом клапане

KP=0.05 #пропорциональный коэффициент регулятора
KI=0.01 #интегральный коэффициент регулятора

#временнЫе настройки расчета
T=10 #общее время расчета
dt=0.1 #шаг расчета по времени
t_in_by_time= {0:10,#изменение температуры воды (С) на входе в систему по времени
               40:10,
               60:10,
               80:10}

#функция для расчета площади проходного сечения регулируемого клапана в зависимости от величины управляющего сигнала: от 0 до 1 (0 - клапан закрыт, 1 - открыт)
def S_valve(t):
    return np.pi*(d/2)**2*t

#функция расчета управляющего сигнала ПИ-регулятора
def PI_controller(error, sum_of_prev_error, KP, KI, dt):
    # return sum_of_prev_error + KP*(error - sum_of_prev_error) + KI*error
    proportional=KP*error
    integral=KI*(sum_of_prev_error+error*dt)
    return proportional+integral

#первоначальные приближения
t_in=10 #температура на входе в систему
t_out=10. #на выходе
valve_position=0.0001
sum_of_prev_error=0.
prev_x=0.99

#Массивы для записи результатов:
time=[t_in]
res_t_out=[np.nan]
res_G1=[np.nan]
res_G2=[np.nan]
res_signal=[np.nan]
res_valve_position=[np.nan]

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
    # res=root_scalar(iteration,x0=prev_x,x1=prev_x*1.0001,method='secant',args=(valve_position))
    res = root_scalar(iteration, bracket=[0.001,0.999999], method='bisect', args=(valve_position))
    x=res.root
    prev_x=x
    #расходы потоков через основной контур и обходную трубу
    V1 = Q_in * x / S
    G1 = Rho*V1*S
    V2 = Q_in * (1 - x) / S
    G2 = Rho * V2 * S
    #температуры потоков до смешения:
    t1=t_in+d_temp
    t2=t_in
    #после смешения:
    t_out=(G1*t1+G2*t2)/Q_in/Rho
    #ошибка регулирования:
    error=t_out-t_out_target
    #считаем через ПИ-регулятор управляющий сигнал
    signal=PI_controller(error, sum_of_prev_error, KP, KI, dt)
    sum_of_prev_error+=error
    # if len(sum_of_prev_error)>10:
    #     del sum_of_prev_error[0]
    #полученный сигнал пересчитываем в вид от 0 до 1 для управления клапаном
    if signal>1:
        valve_position=1
    elif signal<0:
        valve_position=0
    else:
        valve_position=signal

    res_t_out.append(t_out)
    res_G1.append(G1)
    res_G2.append(G2)
    res_signal.append(signal)
    res_valve_position.append(valve_position)
t_out_target
fig1=mc.Chart(points_for_plot=[{'x':time,'y':res_t_out,'label':'t_out','c':'red'},{'x':time,'y':[t_out_target]*len(time),'label':'t_target','c':'black'}],xlabel='time, s',ylabel='t_out, C', title='Temperature', dpi=150,figure_size=(5,5))
fig2=mc.Chart(points_for_plot=[{'x':time,'y':res_G1,'label':'res_G1','c':'red'},{'x':time,'y':res_G2,'label':'res_G2','c':'blue'}],xlabel='time, s',ylabel='G, kg/s', title='Massflow', dpi=150,figure_size=(5,5))
fig3=mc.Chart(points_for_plot=[{'x':time,'y':res_signal,'label':'signal','c':'red'},{'x':time,'y':res_valve_position,'label':'valve_position','c':'blue'}],xlabel='time, s',ylabel='signal, valve_pos', title='PI controller', dpi=150,figure_size=(5,5))

#проверим характеристику регулируемого клапана
v_array=[0.1,1,5]
x_array=np.linspace(0,1,50)
dH1_array=[dH_valve(x,v_array[0],g) for x in x_array]
dH2_array=[dH_valve(x,v_array[1],g) for x in x_array]
dH3_array=[dH_valve(x,v_array[2],g) for x in x_array]
fig4=mc.Chart(points_for_plot=[{'x':x_array,'y':dH1_array,'label':'v=0.1','c':'red'},
                               {'x':x_array,'y':dH2_array,'label':'v=1','c':'blue'},
                               {'x':x_array,'y':dH3_array,'label':'v=5','c':'green'},],xlabel='valve position',ylabel='Hydraulic resistance, m', title='Hydraulic resistance of control valve', dpi=150,figure_size=(5,5))
plt.show()




