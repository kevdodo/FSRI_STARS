
def func(x):
    return x**2
x0=0
xi = 1.5
x1=3

n =int((x1-x0)/xi)
function0= x1**x1
trap=0
trapr = 0
b=0

tot = 0
st = func(x0)
end = func(x1)


con = 3*xi/8
for i in range (1,n):
    x=x0+xi*i
    if i%3 ==0:
        sim =2*func(x)
    else:
        sim = 3*func(x)
    tot += sim 
    
print(con*(tot+st+end))