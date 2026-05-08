#!/usr/bin/python3 env

import math
import numpy as np
#TASK 1

def f(x):
    return x**3 - 3*x**2 + x - 1 #general function

def df(x):
    return 3*x**2 - 6*x + 1    #for Newton's method


#BISECTION METHOD

def bisection(f, a, b, tol = 1e-6, no_iter = 50):
    results = []

    if f(a)*f(b) < 0:
        i = 0
        while i < no_iter:
            xi = (a + b) / 2
            results.append(f(xi))
            if abs(f(xi)) < tol:
                break

            if (f(xi) * f(a) < 0):
                b = xi
            else:
                a = xi

            i += 1
        return xi, results

    else:
        print("Bisection method failed. f(a) and f(b) must have opposite signs")



#SECANT METHOD
def secant(f, x0, x1, tol = 1e-6, no_iter = 50):
    results = []
    f0 = f(x0)
    f1 = f(x1)
    i = 0
    while (i <= no_iter):
        sub_f = f(x1) - f(x0)

        if sub_f == 0:
            print ("Division by 0 error.")
        sub_x = x1 - x0
        del_f = (sub_f) / (sub_x) 
        xi = x1 - f(x1)/del_f
        results.append(f(xi))
        
        if abs(f(xi)) < tol:
            break
        else:
            x0 = x1
            x1 = xi

        i += 1

    return xi, results 



#NEWTON RAPHTON'S METHOD
def newton_r(f, df, x0, tol = 1e-6, no_iter = 50):
    results = []
    i = 0
    while i<= no_iter:
        f0 = f(x0)
        df0 = df(x0)
        xi = x0 - f0 / df0
        results.append(f(xi))

        if abs(f(xi)) < tol:
            break
        else:
            x0 = xi

        if df0 == 0:
            print("Division by 0 error")

        i += 1

    return xi, results


#FIXED-POINT ITERATION METHOD

#Making each x the subject and defining them as g1, g2 and g3 to test which one converges
#better with my initial point

def g1(x):
    return 1 + 3*x**2 - x**3 #I made x the subject of f(x) = 0

def dg1(x):
    return 6*x - 3*x**2  #Its derivative

def g2(x):
    return math.sqrt((1 / 3) * (x**3 - x - 1)) #I made x^2 the subject of f(x) = 0

def dg2(x):
    return (3*x**2 + 1) / (6 * math.sqrt((x**3 + x - 1) / 3))  #Its derivative
    print

def g3(x): 
    return math.cbrt(3*x**2 - x + 1) #I made x^3 the subject 

def dg3(x):
    return (6*x - 1) / (3 * ((3*x**2 - x + 1)**(2 / 3)))  #Its derivative



  #Begin method  

def fixed_p(f, g1, g2, g3, dg1, dg2, dg3, x0, tol = 1e-6, no_iter = 50):
    results = []
    if abs(dg1(x0)) < 1:
        g = g1
        
    elif ((x0**3 + x0 - 1) / 3) >= 0 and abs((dg2(x0))) < 1:
        g = g2
        
    elif abs(dg3(x0)) < 1:
        g = g3

    else:
        return None

    i = 0
    while i < no_iter:
        xi = g(x0)
        results.append(f(xi))
        if abs(f(xi)) < tol:
            break
        else:
            x0 = xi
            i += 1
    return xi, results


#CALLING AND MAKING THE FILE data.dat
bisection_root, bisection_data = bisection(f, -15, 5)
secant_root, secant_data = secant(f, 0.3, 1, 1e-6, 50)
newton_root, newton_data = newton_r(f, df, 0.5)
fixed_root, fixed_data = fixed_p(f, g1, g2, g3, dg1, dg2, dg3, 0.5)

max_iter = max(len(bisection_data), len(secant_data), len(newton_data), len(fixed_data))

with open("data.dat", "w") as file:
    file.write("Iteration  fx_Bisection  fx_Secant      fx_Newton    fx_Fixed\n")
    for i in range(max_iter):
        b = bisection_data[i] if i < len(bisection_data) else bisection_data[-1]
        s = secant_data[i] if i < len(secant_data) else secant_data[-1]
        n = newton_data[i] if i < len(newton_data) else newton_data[-1]
        fp = fixed_data[i] if i < len(fixed_data) else fixed_data[-1]
        file.write(f"{i:<5} {b:<15.7f} {s:<15.6f} {n:<15.7f} {fp:<15.7f}\n")

print("data.dat made successfully.")



#TASK 2
#Euler's Method
def f(t, y):
    return (-2.2067 * 1e-12) * (y**4 - 81*1e8)

def euler(f, t0, y0, h, no_iter = 50):
    while t0 < 480:
        yi = y0 + f(t0, y0) * h
        t0 = t0 + h
        y0 = yi

    return yi

e480 = euler(f, 0, 1200, 480)
e240 = euler(f, 0, 1200, 240)
e120 = euler(f, 0, 1200, 120)
e60 = euler(f, 0, 1200, 60)
e30 = euler(f, 0, 1200, 30)

step_size = [480, 240, 120, 60, 30]
euler_r = [e480, e240, e120, e60, e30]
exact = [1635.4, 537.6, 100.80, 32.607, 14.806]

with open("euler.dat", "w") as file:
    file.write("StepSize    Theta(480)        Exact\n")
    for i in range(5):
        file.write(f"{step_size[i]:<10} {euler_r[i]:<15.6f} {exact[i]:<15.6f}\n")

print("euler.dat made successfully.")



#NAIVE GAUSSIAN ELIMINATION

def gauss_elim(R1, R2, R3):
    NGE = np.array([R1, R2, R3])

    # Forward elimination

    #Step 1
    wing = ( R2[0] / R1[0]) * R1
    R2 = R2 - wing
    NGE = np.array([R1, R2, R3])

    #Step 2
    wing = (R3[0] / R1[0]) * R1
    R3 = R3 - wing
    NGE = np.array([R1, R2, R3])

    #Step 3
    wing = (R3[1] / R2[1]) * R2
    R3 = R3 - wing
    NGE = np.array([R1, R2, R3])
    print("Augmented upper triangular matrix after elimination = \n", NGE)

    #BACK SUBSTITUTING
    print("Back substituting: ")

    #Step 4
    x3 = R3[3] / R3[2]
    x2 = (R2[3] - R2[2]*x3) / R2[1]
    x1 = (R1[3] - R1[2]*x3 - R1[1]*x2) / R1[0]


    return np.array([x1, x2, x3])

R1 = np.array([17, 14, 23, 24.5], dtype=float)
R2 = np.array([-7.54, -3.54, 2.7, 2.352], dtype=float)
R3 = np.array([6, 1, 3, 14], dtype=float)

solution = gauss_elim(R1, R2, R3)
print("Solution:", solution)

#Verifying Solution
R1 = [17, 14, 23]
R2 = [-7.54, -3.54, 2.7]
R3 = [6, 1, 3]

A = np.array([R1, R2, R3], dtype=float)

b = np.array([24.5, 2.352, 14], dtype=float)

x = np.linalg.solve(A, b)

print("Verifying answer using np.linalg.solve(A, b): ")
print("x1 = ", x[0])
print("x2 = ", x[1])
print("x3 = ", x[2])

