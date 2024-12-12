

x0=0
xi = .1
x1=3
function0= x1**x1
a=0
b=0

def func(x0):
    return x0
while x0<=x1:
    #iterate x
    x0 +=xi
    #function statement
    a= func(x0)*xi
    b+=a
print(b)