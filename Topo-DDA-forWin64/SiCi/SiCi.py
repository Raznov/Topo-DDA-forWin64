import numpy as np
import scipy as sp
import math
import matplotlib
import matplotlib.pyplot as plt


def Si(x):
    dis = 0.1
    sum = 0.0
    for i in range(1, int(round(x/dis))):
        sum += dis*math.sin(i*dis)/(i*dis)
    end
    return sum


def Ci(x):
    gamma = 0.577216
    dis = 0.1
    sum = 0.0
    for i in range(1, int(round(x/dis))):
        sum += dis*(1-math.cos(i*dis))/(i*dis)
    end
    return gamma+math.log(x)-sum


def Sid(x):
    return math.sin(x)/x


def Cid(x):
    return (1-math.cos(x))/x


number = 10000
dis = 0.1
gamma = 0.577216

pos = np.linspace(dis, dis*number, number)
pos_plot = np.linspace(dis, dis*(number-1), number-1)
Si_plot = np.zeros(number-1)
Ci_plot = np.zeros(number-1)

sum_Si = 0
sum_Ci = 0
for i in range(pos.shape[0]-1):
    sum_Si += dis*(Sid(pos[i])+Sid(pos[i+1]))/2
    sum_Ci += dis*(Cid(pos[i])+Cid(pos[i+1]))/2
    Si_plot[i] = sum_Si
    Ci_plot[i] = gamma+math.log(pos[i])-sum_Ci

fSi = open("Si.txt", "w")
fCi = open("Ci.txt", "w")
fSi.write(str(number-1)+"\n"+str(dis)+"\n")
fCi.write(str(number-1)+"\n"+str(dis)+"\n")
for i in range(number-1):
    fSi.write(str(Si_plot[i])+"\n")
    fCi.write(str(Ci_plot[i])+"\n")
fSi.close()
fCi.close()


##plt.plot(pos_plot,Si_plot)
##plt.plot(pos_plot,Ci_plot)
##plt.show()