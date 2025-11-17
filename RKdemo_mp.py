import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def RK1Solve(f, y0, nsteps, x0, xmax) -> np.array:
    h = (xmax - x0) / nsteps  # step size
    x = x0                     # independent variable
    y = y0                     # dependent variable to plot vs x
    points = [(x0, y0)]        # store points for plotting

    for i in range(nsteps - 1):
        k1 = h * f(x, y)
        y = y + k1
        x += h
        points.append((x, y))
    return np.array(points)

def RK2Solve(f, y0, nsteps, x0, xmax) -> np.array:
    h = (xmax - x0) / nsteps  # step size
    x = x0                     # independent variable
    y = y0                     # dependent variable to plot vs x
    points = [(x0, y0)]        # store points for plotting

    for i in range(nsteps - 1):
        k1 = h * f(x, y)
        k2 = h * f(x + h / 2, y + k1 / 2)
        y = y + k2
        x += h
        points.append((x, y))
    return np.array(points)

def RK4Solve(f, y0, nsteps, x0, xmax) -> np.array:
    h = (xmax - x0) / nsteps
    x = x0
    y = y0
    points = [(x0, y0)]
    for i in range(nsteps - 1):
        k1 = h * f(x, y)
        k2 = h * f(x + h / 2, y + k1 / 2)
        k3 = h * f(x + h / 2, y + k2 / 2)
        k4 = h * f(x + h, y + k3)
        y = y + (k1 + 2*k2 + 2*k3 + k4) / 6
        x += h
        points.append((x, y))
    return np.array(points)

# The differential equation to be solved
def fun1(x, y):
    return 2*x -1-3* y  
def exact_solution(x):
    return (2*x)/3 -(5/9)+(32/9) * np.exp(-3*x)
y0 = 3
nsteps = 30
x0 = 0
xmax = 3

# Solve our DEQ using RK1 or RK2 methods!
tg1 = RK1Solve(fun1, 3, 30, 0, 3)  # initial condition y(0)=3
tg2 = RK2Solve(fun1, 3, 30, 0, 3)
tg4 = RK4Solve(fun1, 3, 30, 0, 3)
x_exact = np.linspace(0, 3, 300)
y_exact = (2*x_exact)/3 -(5/9)+(32/9) * np.exp(-3*x_exact)  # exact solution

err1 = np.abs(tg1[:, 1] - exact_solution(tg1[:, 0]))
err2 = np.abs(tg2[:, 1] - exact_solution(tg2[:, 0]))
err4 = np.abs(tg4[:, 1] - exact_solution(tg4[:, 0]))

with PdfPages('RK4.pdf') as pdf:
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.5, 0.85, 'Runge-Kutta Methods Comparison', 
             ha='center', fontsize=20, weight='bold')
    info_text = f"""
ODE Problem:
    y' + 3y = 2x - 1
    Also written as: dy/dx = 2x - 1 - 3y
    Initial condition: y(0) = 3
    Domain: x=∈ [0, 3]
    Number of steps: {nsteps}
    Step size h = {(xmax-x0)/nsteps:.4f}

Exact Solution:
    y(x) = (2x)/3 - 5/9 + (32/9)·exp(-3x)

Methods:
    • RK1: First-order method
    • RK2: Second-order method
    • RK4: Fourth-order method

Error:
    RK4 vs RK1: {err1[-1]/err4[-1]:.1f}x more accurate
    RK4 vs RK2: {err2[-1]/err4[-1]:.1f}x more accurate

Conclusion:
    RK4 demonstrates far better accuracy compared to both
    RK1 and RK2 methods
    """
    
    fig.text(0.5, 0.45, info_text, fontsize=11, family='monospace',
             verticalalignment='center', horizontalalignment='center')
    plt.axis('off')
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()
    
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(tg1[:, 0], tg1[:, 1], 'r^', markersize=8, label='RK1 Solution', 
            markeredgewidth=1.5)
    ax.plot(tg2[:, 0], tg2[:, 1], 'gv', markersize=8, label='RK2 Solution', 
            markeredgewidth=1.5)
    ax.plot(tg4[:, 0], tg4[:, 1], 'bo', markersize=7, label='RK4 Solution', 
            fillstyle='none', markeredgewidth=1.5)
    ax.plot(x_exact, y_exact, 'k--', label='Exact Solution')
    
    ax.set_title('ODE Solutions: y\' + 3y = 2x - 1, y(0) = 3', fontsize=16, weight='bold')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend(fontsize=12, loc='best')
    ax.grid(True, alpha=0.3, linestyle='--')
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()
    
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.semilogy(tg1[:, 0], err1, 'r^-', markersize=7, linewidth=1.5, 
                 label='RK1 Error', markeredgewidth=1.5)
    ax.semilogy(tg2[:, 0], err2, 'gv-', markersize=7, linewidth=1.5, 
                 label='RK2 Error', markeredgewidth=1.5)
    ax.semilogy(tg4[:, 0], err4, 'bo-', markersize=6, linewidth=1.5, 
                 label='RK4 Error', fillstyle='none', markeredgewidth=1.5)
    ax.set_title('Absolute Error', fontsize=14, weight='bold')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('Absolute Error', fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()
