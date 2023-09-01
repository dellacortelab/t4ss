import numpy as np

n_lambdas = 31
n_iterations = 100

u_kln = np.array([])
for i in range(0,n_lambdas):
  print(f'Grabbing u_{i}ln')
  a = np.loadtxt(f'lambda_{i}/u_kln.txt')
  a = a.reshape(n_lambdas,n_lambdas,n_iterations)
  a = a[i].reshape(1,n_lambdas,n_iterations)
  if len(u_kln) == 0:
    u_kln = a
  else:
    u_kln = np.concatenate((u_kln,a),axis=0)
print('u_kln shape: ',u_kln.shape)

newarr = u_kln.reshape(u_kln.shape[0], (u_kln.shape[1]*u_kln.shape[2]))
np.savetxt("u_kln.txt",newarr)