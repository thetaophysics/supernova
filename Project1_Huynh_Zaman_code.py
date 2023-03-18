#Declare libraries 
import numpy as np
from scipy import integrate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#Declare constants and take input values
e_Ni = 3.3696*pow(10,15)   #[erg*g/day]
e_Co = 5.8752*pow(10,14)    #[erg*g/day]
ra = pow(10,13)
# ra0 = 0.
# ra1 = 1.e11
# ra2 = 1.e12
# ra3 = 1.e14
ma_Ni = 0.1*1.99*pow(10,33) #[g]
# ma_Ni_min = 0.05*1.99*pow(10,33) #[g]
# ma_Ni_max = 1.99*pow(10,33) #[g]
# b = float(input("Enter M_ej:[solar mass unit) "))
ma_ej = 7*(1.99*pow(10,33)) #[g]
# ma_ej_min = 3*(1.99*pow(10,33)) #[g]
# ma_ej_max = 15*(1.99*pow(10,33)) #[g]
# c = float(input("Enter SN ejecta velocity: "))
velocity = 7000*86400*10**5  #[km/day]
l = []
# lo = []

#func0: diffusion time with input M_ej and velocity
def t_d(m_ej, velo):    
    beta = 13.7 #integration constant
    c = 3*10**10*86400  #cm/day
    k = 0.2 #cm/g optical opacity 
    # diff_time = np.sqrt((10./3.)*(k/(beta*c))*m_ej*(1/velo))
    diff_time = np.sqrt((10.*k*m_ej)/(3.*beta*c*velo))
    return diff_time

#func1: first integral for Ni with input R_o, velocity, and time
def func1(r0, v, M_ej, t):
    t_Ni = 113  #characteristic half-life time scales f or radiactive Ni [days]
    dif1 = t_d(M_ej,v)
    x2 = (t**2/dif1**2)+((2.*r0*t)/(v*dif1**2))  #x2 is inside integral
    f1 = (r0/(v*dif1)+(t/dif1))*np.exp(x2)*np.exp(-t/t_Ni)
    return f1
def func2(r0, v, M_ej, t):
    t_Co = 8.8  #characteristic half-life time scales f or radiactive Ni [days]
    dif2 = t_d(M_ej,v)
    x2 = (t**2/dif2**2)+((2.*r0*t)/(v*dif2**2))  #x2 is inside integral
    f2 = (r0/(v*dif2)+(t/dif2))*np.exp(x2)*np.exp(-t/t_Co)
    return f2

#Changing masses for #4
# ma_Ni_1 = 0.01*1.99*pow(10,33) #[g]
# ma_Ni_2 = 0.5*pow(10,33) #[g]

#Changing only M_ej
# ma_ej_1 = 3.*(1.99*pow(10,33)) #[g]
# ma_ej_2 = 5.*(1.99*pow(10,33)) #[g]
# ma_ej_3 = 7.*(1.99*pow(10,33)) #[g]
# ma_ej_4 = 9.*(1.99*pow(10,33)) #[g]
# ma_ej_5 = 11.*(1.99*pow(10,33)) #[g]
# ma_ej_6 = 15.*(1.99*pow(10,33)) #[g]
# ma_ej_7 = 20.*(1.99*pow(10,33)) #[g]
# l1 = []
# l2 = []
# l3 = []
# l4 = []
# l5 = []
# l6 = []
# l7 = []
# dif_t1 = t_d(ma_ej_1, velocity)
# dif_t2 = t_d(ma_ej_2, velocity)
# dif_t3 = t_d(ma_ej_3, velocity)
# dif_t4 = t_d(ma_ej_4, velocity)
# dif_t5 = t_d(ma_ej_5, velocity)
# dif_t6 = t_d(ma_ej_6, velocity)
# dif_t7 = t_d(ma_ej_7, velocity)


# Call diffusion time funct
dif_t = t_d(ma_ej, velocity)


#Perform integration using simpson rule
gridsize = 1000   #using Simpson requires odd number of intervals
t = np.linspace(0.,400., 401)

for i in t:
    #Original Graph
    grid = np.linspace(0., i, gridsize)
    outs_x = (i**2/dif_t**2) + ((2*ra*i)/(velocity*dif_t**2))
    f1 = func1(ra, velocity, ma_ej, grid)
    f2 = func2(ra, velocity, ma_ej, grid)
    f1_integrate = integrate.simps(f1)    #integration for first integral
    f2_integrate = integrate.simps(f2)    #integration for 2nd integral
    lum = ((2*ma_Ni)/dif_t)*np.exp(-outs_x)*((e_Ni-e_Co)*f1_integrate + e_Co*f2_integrate)
    l.append(lum/86400.)

    # outs_x0 = (i**2/dif_t**2) + ((2*ra0*i)/(velocity*dif_t**2))
    # f1o = func1(ra0, velocity, ma_ej, grid)
    # f2o = func2(ra0, velocity, ma_ej, grid)
    # f1o_integrate = integrate.simps(f1)    #integration for first integral
    # f2o_integrate = integrate.simps(f2)    #integration for 2nd integral
    # lumo = ((2*ma_Ni)/dif_t)*np.exp(-outs_x0)*((e_Ni-e_Co)*f1o_integrate + e_Co*f2o_integrate)
    # lo.append(lumo/86400.)

    # outs_x1 = (i**2/dif_t**2) + ((2*ra1*i)/(velocity*dif_t**2))
    # f11 = func1(ra1, velocity, ma_ej, grid)
    # f21 = func2(ra1, velocity, ma_ej, grid)
    # f11_integrate = integrate.simps(f11)    #integration for first integral
    # f21_integrate = integrate.simps(f21)    #integration for 2nd integral
    # lum1 = ((2*ma_Ni)/dif_t)*np.exp(-outs_x1)*((e_Ni-e_Co)*f11_integrate + e_Co*f21_integrate)
    # l1.append(lum1/86400.)

    # outs_x2 = (i**2/dif_t**2) + ((2*ra2*i)/(velocity*dif_t**2))
    # f12 = func1(ra2, velocity, ma_ej, grid)
    # f22 = func2(ra2, velocity, ma_ej, grid)
    # f12_integrate = integrate.simps(f12)    #integration for first integral
    # f22_integrate = integrate.simps(f22)    #integration for 2nd integral
    # lum2 = ((2*ma_Ni)/dif_t)*np.exp(-outs_x2)*((e_Ni-e_Co)*f12_integrate + e_Co*f22_integrate)
    # l2.append(lum2/86400.)

    # outs_x3 = (i**2/dif_t**2) + ((2*ra3*i)/(velocity*dif_t**2))
    # f13 = func1(ra3, velocity, ma_ej, grid)
    # f23 = func2(ra3, velocity, ma_ej, grid)
    # f13_integrate = integrate.simps(f13)    #integration for first integral
    # f23_integrate = integrate.simps(f23)    #integration for 2nd integral
    # lum3 = ((2*ma_Ni)/dif_t)*np.exp(-outs_x3)*((e_Ni-e_Co)*f13_integrate + e_Co*f23_integrate)
    # l3.append(lum3/86400.)


    # #Calculations for chhanging M_ej only
    # outs_x1 = (i**2/dif_t1**2) + ((2*ra*i)/(velocity*dif_t1**2))
    # outs_x2 = (i**2/dif_t2**2) + ((2*ra*i)/(velocity*dif_t2**2))
    # outs_x3 = (i**2/dif_t3**2) + ((2*ra*i)/(velocity*dif_t3**2))
    # outs_x4 = (i**2/dif_t4**2) + ((2*ra*i)/(velocity*dif_t4**2))
    # outs_x5 = (i**2/dif_t5**2) + ((2*ra*i)/(velocity*dif_t5**2))
    # outs_x6 = (i**2/dif_t6**2) + ((2*ra*i)/(velocity*dif_t6**2))
    # outs_x7 = (i**2/dif_t7**2) + ((2*ra*i)/(velocity*dif_t7**2))
    # f11 = func1(ra, velocity, ma_ej_1, grid)
    # f21 = func2(ra, velocity, ma_ej_1, grid)
    # f1_integrate1 = integrate.simps(f11)    #integration for first integral
    # f2_integrate1 = integrate.simps(f21)    #integration for 2nd integral
    # lum1 = ((2*ma_Ni)/dif_t1)*np.exp(-outs_x1)*((e_Ni-e_Co)*f1_integrate1 + e_Co*f2_integrate1)
    # l1.append(lum1/86400.)

    # f12 = func1(ra, velocity, ma_ej_2, grid)
    # f22 = func2(ra, velocity, ma_ej_2, grid)
    # f1_integrate2 = integrate.simps(f12)    #integration for first integral
    # f2_integrate2 = integrate.simps(f22)    #integration for 2nd integral
    # lum2 = ((2*ma_Ni)/dif_t2)*np.exp(-outs_x2)*((e_Ni-e_Co)*f1_integrate2 + e_Co*f2_integrate2)
    # l2.append(lum2/86400.)

    # f13 = func1(ra, velocity, ma_ej_3, grid)
    # f23 = func2(ra, velocity, ma_ej_3, grid)
    # f1_integrate3 = integrate.simps(f13)    #integration for first integral
    # f2_integrate3 = integrate.simps(f23)    #integration for 2nd integral
    # lum3 = ((2*ma_Ni)/dif_t3)*np.exp(-outs_x3)*((e_Ni-e_Co)*f1_integrate3 + e_Co*f2_integrate3)
    # l3.append(lum3/86400.)

    # f14 = func1(ra, velocity, ma_ej_4, grid)
    # f24 = func2(ra, velocity, ma_ej_4, grid)
    # f1_integrate4 = integrate.simps(f14)    #integration for first integral
    # f2_integrate4 = integrate.simps(f24)    #integration for 2nd integral
    # lum4 = ((2*ma_Ni)/dif_t4)*np.exp(-outs_x4)*((e_Ni-e_Co)*f1_integrate4 + e_Co*f2_integrate4)
    # l4.append(lum4/86400.)

    # f15 = func1(ra, velocity, ma_ej_5, grid)
    # f25 = func2(ra, velocity, ma_ej_5, grid)
    # f1_integrate5 = integrate.simps(f15)    #integration for first integral
    # f2_integrate5 = integrate.simps(f25)    #integration for 2nd integral
    # lum5 = ((2*ma_Ni)/dif_t5)*np.exp(-outs_x5)*((e_Ni-e_Co)*f1_integrate5 + e_Co*f2_integrate5)
    # l5.append(lum5/86400.)

    # f16 = func1(ra, velocity, ma_ej_6, grid)
    # f26 = func2(ra, velocity, ma_ej_6, grid)
    # f1_integrate6 = integrate.simps(f16)    #integration for first integral
    # f2_integrate6 = integrate.simps(f26)    #integration for 2nd integral
    # lum6 = ((2*ma_Ni)/dif_t6)*np.exp(-outs_x6)*((e_Ni-e_Co)*f1_integrate6 + e_Co*f2_integrate6)
    # l6.append(lum6/86400.)

    # f17 = func1(ra, velocity, ma_ej_7, grid)
    # f27 = func2(ra, velocity, ma_ej_7, grid)
    # f1_integrate7 = integrate.simps(f17)    #integration for first integral
    # f2_integrate7 = integrate.simps(f27)    #integration for 2nd integral
    # lum7 = ((2*ma_Ni)/dif_t7)*np.exp(-outs_x7)*((e_Ni-e_Co)*f1_integrate7 + e_Co*f2_integrate7)
    # l7.append(lum7/86400.)


    #Calculations for various masses
    # outs_x1 = (i**2/dif_t1**2) + ((2*ra*i)/(velocity*dif_t1**2))
    # outs_x2 = (i**2/dif_t2**2) + ((2*ra*i)/(velocity*dif_t2**2))
    # #First graph for M_Ni1 and M_ej1
    # f1A = func1(ra, velocity, ma_ej_1, grid)
    # f2A = func2(ra, velocity, ma_ej_1, grid)
    # f1_integrateA = integrate.simps(f1A)    #integration for first integral
    # f2_integrateA = integrate.simps(f2A)    #integration for 2nd integral
    # lum1 = ((2*ma_Ni_1)/dif_t1)*np.exp(-outs_x1)*((e_Ni-e_Co)*f1_integrateA + e_Co*f2_integrateA)
    # l1.append(lum1/86400.)

    # #Graph 2 for M_ni1 and M_ej2
    # f1B = func1(ra, velocity, ma_ej_2, grid)
    # f2B = func2(ra, velocity, ma_ej_2, grid)
    # f1_integrateB = integrate.simps(f1B)    #integration for first integral
    # f2_integrateB = integrate.simps(f2B)    #integration for 2nd integral
    # lum2 = ((2*ma_Ni_1)/dif_t2)*np.exp(-outs_x2)*((e_Ni-e_Co)*f1_integrateB + e_Co*f2_integrateB)
    # l2.append(lum2/86400.)

    # #Graph 3 for M_ni2 and  M_ej1
    # f1C = func1(ra, velocity, ma_ej_1, grid)
    # f2C = func2(ra, velocity, ma_ej_1, grid)
    # f1_integrateC = integrate.simps(f1C)    #integration for first integral
    # f2_integrateC = integrate.simps(f2C)    #integration for 2nd integral
    # lum3 = ((2*ma_Ni_2)/dif_t1)*np.exp(-outs_x1)*((e_Ni-e_Co)*f1_integrateC + e_Co*f2_integrateC) 
    # l3.append(lum3/86400.)

    # #Graph 4 for M_ni2 and max M_ej2
    # f1D = func1(ra, velocity, ma_ej_2, grid)
    # f2D = func2(ra, velocity, ma_ej_2, grid)
    # f1_integrateD = integrate.simps(f1D)    #integration for first integral
    # f2_integrateD = integrate.simps(f2D)    #integration for 2nd integral
    # lum4 = ((2*ma_Ni_2)/dif_t2)*np.exp(-outs_x2)*((e_Ni-e_Co)*f1_integrateD + e_Co*f2_integrateD)
    # l4.append(lum4/86400.)




# Part 3
# Getting the max index to cut the lum list from the max_value
max_lum = np.max(l)
max_index = l.index(max_lum)
print(max_index, max_lum)
new_lum = l[max_index:]
new_t = t[max_index:]

#Define function for fitting f(t) = Ae^(−Bt)
def fit_funct(t,a,b):
    return (a*np.exp(-b*t))

popt,pcov = curve_fit(fit_funct ,new_t,new_lum,bounds=((1e44,0.01),(2e44,0.02)))

print (popt)

#Make plots using matplotlib
# Set a title for the plot
plt.title('Supernova Lightcurve Models')
# plt.title('Peak Luminosities with changing Initial Radii')
#Set label for x-axis 
plt.xlabel('Time (days)')
# plt.ylabel('Luminosity(erg/day)')
plt.ylabel('log [Luminosity(erg/s)]')
plt.plot(t,l,'r', label = 'M_Ni = 0.1M_sun, M_ej = 7.0M_sun')
# plt.plot(t,l,'r', label = 'Radii = 10^13cm')
# plt.plot(t,lo,'g', label = 'Radii = 0cm')
# plt.plot(t,l1,'b', label = 'Radii = 10^11cm')
# plt.plot(t,l2,'y', label = 'Radii = 10^12cm')
# plt.plot(t,l3,'c', label = 'Radii = 10^14cm')

# plt.plot(t,l,'r', label = 'M_Ni = 0.1M_sun, M_ej = 7.0M_sun')
#Plot the cut at max luminosity
# only one line may be specified; full height
plt.axvline(x=max_index)
# Plot fit with A and B guesses
plt.plot(new_t, fit_funct(new_t, *popt), 'b', label = 'fit')
# plt.plot(t, np.log10(l1),'b' , label = 'M_Ni=0.1M_sun and M_ej=3M_sun' )
# plt.plot(t, np.log10(l2),'g',label = 'M_Ni=0.1M_sun and M_ej= 5M_sun' )
# plt.plot(t, np.log10(l3),'r',label = 'M_Ni=0.1M_sun and M_ej=7M_sun' )
# plt.plot(t, np.log10(l4),'c',label = 'M_Ni=0.1M_sun and M_ej=9M_sun' )
# plt.plot(t, np.log10(l5),'m',label = 'M_Ni=0.1M_sun and M_ej=11M_sun' )
# plt.plot(t, np.log10(l6),'y',label = 'M_Ni=0.1M_sun and M_ej=15M_sun' )
# plt.plot(t, np.log10(l7),'k',label = 'M_Ni=0.1M_sun and M_ej=20M_sun' )

plt.legend()
plt.xlim(0., 400.)
# plt.xlim(20., 90.)
# plt.ylim(3.7e43, 4.6e43)
plt.show()


