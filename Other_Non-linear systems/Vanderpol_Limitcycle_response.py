
"""
##Vanderpol Oscillator
from scipy import linspace
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def vdp(t, z):
    x, y = z
    return [y, mu*(1 - x**2)*y - x]
#############inputs############################
a, b = 0, 200######time
mus = [0.3]
x1=0.1; y1=0
Z=[]
###############################
styles = ["-", "--", ":"]
t = linspace(a, b, 5000)
for mu, style in zip(mus, styles):
    sol = solve_ivp(vdp, [a, b], [x1, y1], t_eval=t)
    plt.plot(sol.y[0], sol.y[1], style)
# make a little extra horizontal room for legend
plt.xlim([-4,4])    
plt.legend([f"$\mu={m}$" for m in mus])
#plt.axes().set_aspect(1)
plt.plot(sol.y[0], sol.y[1], style)
plt.show()
plt.plot(sol.t, sol.y[0], style)
plt.show()
