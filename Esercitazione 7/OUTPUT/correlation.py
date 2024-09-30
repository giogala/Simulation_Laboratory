import numpy as np


dir = ['Gas','Liquid','Solid']

chi=np.zeros([500000,3])
for k in range(0,3):
    df_temp= np.genfromtxt(dir[k]+'/Statistic/potential_energy.dat')
    #df_temp = pd.read_csv(dir[k]+'/Statistic/potential_energy.dat',sep='\t')
    x=df_temp[:,0]
    y=df_temp[:,1]
    err=df_temp[:,3]
    file = open(dir[k]+'/Statistic/correlation.txt',"w")
    
    ym_q = (sum(y)/len(x))**2
    yq_m = sum(y**2)/len(x)
    for i in range(0,200):
        mm=0
        m2=0
        m=0
        for j in range(0,len(x)-i):
            mm+=y[j]*y[j+i]/(len(x)-i)
            m2+=y[j+i]/(len(x)-i)
            m+=y[j]/(len(x)-i)
        chi[i,k]= (mm - m2*m)
        print(i)
    chi[:,k] /= (yq_m - ym_q)
    for i in range(0,200):
        file.write(str(i)+'\t'+str(chi[i,k])+'\n')


