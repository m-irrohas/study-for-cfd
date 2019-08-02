import numpy as np
import matplotlib.pyplot as plt

# 定数
gamma = 1.4
# 初期条件
x_min = -2.
x_max = 2.
dx = 0.01
dt = 0.001
nx = int((x_max-x_min)/dx)
nt = 200
pl = 1.0
pr = 0.1
rhol = 1.0
rhor = 0.125
Ul = 0.0
Ur = 0.0

X = np.linspace(x_min, x_max, nx)

def get_e(p, rho, u):
    return p/(gamma-1.)+1./2.*rho*u**2

def get_p(e, rho, u):
    return (gamma-1.)*(e-1./2.*rho*u**2)

def get_E(Q):
    E = np.empty_like(Q)
    rho = Q[:,0]
    m = Q[:,1]
    e = Q[:,2]
    E[:,0] = m
    E[:,1] = (gamma-1)*e+1/2*(3-gamma)*m**2/rho
    E[:,2] = gamma*e*m/rho-1/2*(gamma-1)*m**3/rho**2
    return E

def get_E_vec(Q_vec):
    rho = Q_vec[0]
    m = Q_vec[1]
    e = Q_vec[2]
    E = np.array([m,
                (gamma-1)*e+1/2*(3-gamma)*m**2/rho,
                gamma*e*m/rho-1/2*(gamma-1)*m**3/rho**2])
    return E

# 初期化
Q = np.empty((nx, 3))
for i in range(nx):
    if X[i]<=0:
        q0 = rhol
        q1 = rhol*Ul
        q2 = get_e(pl, rhol, Ul)
    elif X[i] > 0:
        q0 = rhor
        q1 = rhor*Ur
        q2 = get_e(pr, rhor, Ur)
    else:
        raise Exception("okasii")
    Q[i,:] = np.array([q0,q1,q2])

print(Q.shape)
E = get_E(Q)
print(E.shape)
print("ここから数値解析")

print("通常の方法")
eps = 3.3
for t in range(nt):
    Q_old = Q
    E = get_E(Q)
    for i in range(1,nx-1):
        Q_mid = Q_old[i,:] - dt/dx*(E[i,:]-E[i-1,:])
        Q_mid_next = Q_old[i+1,:] - dt/dx*(E[i+1,:]-E[i,:])

        E_mid = get_E_vec(Q_mid)
        E_mid_next = get_E_vec(Q_mid_next)

        p_pre = get_p(Q_old[i-1,2], Q_old[i-1,0], Q_old[i-1,1]/Q_old[i-1,0])
        p = get_p(Q_old[i,2], Q_old[i,0], Q_old[i,1]/Q_old[i,0])
        p_next = get_p(Q_old[i+1,2], Q_old[i+1,0], Q_old[i+1,1]/Q_old[i+1,0])

        visc_co = eps*abs((p_next-2*p+p_pre))/abs((p_next+2*p+p_pre))
        visc = visc_co*(Q_old[i+1,:]-2*Q_old[i,:]+Q_old[i-1,:])
        Q[i,:] = 1/2*(Q[i,:]+Q_mid)-1/2*dt/dx*(E_mid_next-E_mid) + visc

e_result = Q[:,2]
rho_result = Q[:,0]
u_result = Q[:,1]/Q[:,0]
P_result = get_p(e_result, rho_result, u_result)
plt.clf()
plt.plot(X,P_result, label="Pressure", color="pink")
plt.legend()
plt.savefig("./output/maccormac/pressure.png")
plt.show()
plt.clf()
plt.plot(X,rho_result, label="density", color="violet")
plt.legend()
plt.savefig("./output/maccormac/density.png")
plt.show()
plt.clf()
plt.plot(X,u_result, label="velocity", color="lime")
plt.legend()
plt.savefig("./output/maccormac/velocity.png")
plt.show()
plt.clf()
plt.plot(X,e_result, label="hi-energy", color="coral")
plt.legend()
plt.savefig("./output/maccormac/energy.png")
plt.show()

# 初期化
Q = np.empty((nx, 3))
for i in range(nx):
    if X[i]<=0:
        q0 = rhol
        q1 = rhol*Ul
        q2 = get_e(pl, rhol, Ul)
    elif X[i] > 0:
        q0 = rhor
        q1 = rhor*Ur
        q2 = get_e(pr, rhor, Ur)
    else:
        raise Exception("okasii")
    Q[i,:] = np.array([q0,q1,q2])

print(Q.shape)
E = get_E(Q)
print(E.shape)



print("Roeスキームで解く。")
def get_flux_from_roe(Q, Q_next):
    """
    Q_jとQ_{j+1}から中間の数値流速をもとめる。
    Q_j ndarray(float) [rho, rho*u, e]で構成される
    Q_{j+1} ndarray(float) 同じ。
    """
    rho_pre = Q[0]
    u_pre = Q[1]/rho_pre
    e_pre = Q[2]
    p_pre = get_p(e_pre, rho_pre, u_pre)
    H_pre = (e_pre+p_pre)/rho_pre
    E_pre = np.array([Q[1],
                    (gamma-1)*e_pre+1/2*(3-gamma)*Q[1]**2/rho_pre,
                    gamma*e_pre*Q[1]/rho_pre-1/2*(gamma-1)*Q[1]**3/rho_pre**2])

    rho_next = Q_next[0]
    u_next = Q_next[1]/rho_next
    e_next = Q_next[2]
    p_next = get_p(e_next, rho_next, u_next)
    H_next = (e_next+p_next)/rho_next
    E_next = np.array([Q_next[1],
                    (gamma-1)*e_next+1/2*(3-gamma)*Q_next[1]**2/rho_next,
                    gamma*e_next*Q_next[1]/rho_next-1/2*(gamma-1)*Q_next[1]**3/rho_next**2])

    r_rho_pre = np.sqrt(rho_pre)
    r_rho_next = np.sqrt(rho_next)
    rho_ave = r_rho_pre*r_rho_next
    u_ave = (r_rho_pre*u_pre + r_rho_next*u_next)/(r_rho_pre + r_rho_next)
    H_ave = (r_rho_pre*H_pre + r_rho_next*H_next)/(r_rho_pre + r_rho_next)
    c_ave = np.sqrt((gamma-1.)*(H_ave-1./2.*u_ave**2))

    # 対角化
    Lambda = np.diag(np.array([abs(u_ave-c_ave), abs(u_ave), abs(u_ave+c_ave)]))
    R = np.array([[1, 1, 1],
                [u_ave-c_ave, u_ave, u_ave+c_ave],
                [H_ave-u_ave*c_ave, 1./2.*u_ave**2, H_ave+u_ave*c_ave]])
    b1 = u_ave**2*(gamma-1)/2/c_ave**2
    b2 = (gamma-1)/c_ave**2
    R_inv = np.array([[1./2.*(b1+u_ave/c_ave), -1./2.*(1/c_ave+b2*u_ave), 1./2.*b2],
                    [1-b1, b2*u_ave, -b2],
                    [1./2.*(b1-u_ave/c_ave), 1./2.*(1./c_ave-b2*u_ave), 1./2.*b2]])
    A = R.dot(Lambda).dot(R_inv)
    Q_next = Q_next.reshape((3,1))
    Q = Q.reshape((3,1))
    E_next = E_next.reshape((3,1))
    E_pre = E_pre.reshape((3,1))
    E_flux = 1./2.*(E_next+E_pre-A.dot(Q_next-Q))
    return E_flux


for t in range(nt):
    Q_old = Q
    for i in range(1, nx-1):
        #目標 Q_i
        E_flux_pre = get_flux_from_roe(Q_old[i-1, :], Q_old[i, :])
        E_flux_next = get_flux_from_roe(Q_old[i, :], Q_old[i+1, :])
        Q[i, :] = Q_old[i, :] - dt/dx*(E_flux_next-E_flux_pre).reshape(1,3)

print(Q)
e_result = Q[:,2]
rho_result = Q[:,0]
u_result = Q[:,1]/Q[:,0]
P_result = get_p(e_result, rho_result, u_result)
plt.clf()
plt.plot(X,P_result, label="Pressure", color="pink")
plt.legend()
plt.savefig("./output/roe/pressure.png")
plt.show()
plt.clf()
plt.plot(X,rho_result, label="density", color="violet")
plt.legend()
plt.savefig("./output/roe/density.png")
plt.show()
plt.clf()
plt.plot(X,u_result, label="velocity", color="lime")
plt.legend()
plt.savefig("./output/roe/velocity.png")
plt.show()
plt.clf()
plt.plot(X,e_result, label="hi-energy", color="coral")
plt.legend()
plt.savefig("./output/roe/energy.png")
plt.show()
