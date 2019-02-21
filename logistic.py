from scipy import stats
from scipy.stats import logistic

f=open("list.txt","r")
ls=f.readlines()
f.close()
for i in range(len(ls)):
    ls[i]=float(ls[i])

param=logistic.fit(ls)
mode=param[0]
scale=param[1]
mean=logistic.mean(loc=param[0], scale=param[1])
print(mean)
print(scale)
