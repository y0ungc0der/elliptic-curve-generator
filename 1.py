import random
def Diskr():
    D = random.randint(-10**7, -10**6 - 1)
    K = QuadraticField(D) #(D,t)?
    N = K.class_number()
    while ((N<=3600) or ((D%3)==0 or (D%8)!=1)):
        D = random.randint(-10**7, -10**6 - 1)
        K = QuadraticField(D) #(D,t)?
        N = K.class_number()
    if ((N>=3600) and ((D%3)!=0 or (D%8)==1)):
        return (D)
        
# рассчет р, вычисление разложения N/c=R, где N=p+1+a или N=p+1-a
def P(D):
    a=2**120;
    while (a<2**125-11111):
        b=2**115
        while(b<2**120):
            p=((a**2)*(abs(D) + 1) // 4) + abs(D)*a*b + abs(D)*(b**2)
            #(abs(D)+1)/4*a**2+abs(D)*a*b+abs(D)*b**2
            
            if p in Primes():
                break
            b+=1
        if p in Primes():
            break
        a+=1     
    
    if p in Primes():
        print(p)
    else:
        print("error, not prime")
    print(len(bin(p))-2)
    
    N1=p+1+a
    N2=p+1-a
    c=1
    while True:
        r1=N1//c
        r2=N2//c
        if((N1%c) == 0 and r1 in Primes()):
            print("N1",N1%c,c)
            return N1,c,r1,p
        elif((N2%c) == 0 and r2 in Primes()):
            print("N2",N2%c,c)
            return N2,c,r2,p
        else:
            c=c+1
            if(c>50):
                print ("a lot of")
                return 0, 0, 0,0
Dk=Diskr()
N,C,R,p=P(Dk)
while (p==0):
    Dk=Diskr()
    N,C,R,p=P(Dk)
print(N,C,R,Dk)

#Вычисление символа Лежандра для нахождения координаты x
def fun(K,A,B,p):
    x=1
    expr=K(x**(3)+A*x+B)
    l=kronecker_symbol(expr,p)
    while(l!=1):
        x=x+1
        expr=K(x**(3)+A*x+B)
        l=kronecker_symbol(expr,p)
    return x
#проверка на принадлежность точки кривой
def ch(x,y,A,B,p):
    e=(x**3+A*x+B)%p
    if (e==y**(2)%p):
        return True
    return False
# рассчет коэффициентов A,B, построение кривой,проверка корректности параметров
def coffi(p,j,C,R):
    K=GF(p)
    A = K(3*j / (1728 - j))
    B = K(2*j / (1728 - j))

    E = EllipticCurve(K,[A,B])
    x=fun(K,A,B,p)
    e=K(x**(3)+A*x+B)
    y=sqrt(e)%p
    
    print("res",ch(x,y,A,B,p))
    P=E([x,y])
    rup=R*P
    
    rupc=C*P*R
    if rup[0]==0 and rup[1]==1 and rup[2]==0:
        print(R,"is order")
        print(E,P)
    elif(rupc[0]==0 and rupc[1]==1 and rupc[2]==0):
        print(R,"is orderc")
        print(E,P)
        
    else:
        print("is not orderrrr")
        t=2
        while(t<20):
            
            l=kronecker_symbol(t,p)
            if(l==-1):
                break
            else:
                t+=1
        A=K(A*(t**2))
        B=K(B*(t**3))
        E = EllipticCurve(K,[A,B])
        x=fun(K,A,B,p)
        e=K(x**(3)+A*x+B)
        y=sqrt(e)%p
        print(x,y)
        print("res",ch(x,y,A,B,p))
        P=E([x,y])
        rup=R*P
        rupc=C*P*R
        if rup[0]==0 and rup[1]==1 and rup[2]==0:
            print(R,"is order")
            print(E,P)
        elif(rupc[0]==0 and rupc[1]==1 and rupc[2]==0):
            print(R,"is orderc")
            print(E,P)
        else:
            print("=(((")        
Disk = -6094079
R = 94963127712575456382870935540246696255536192341071950639817830665100231338089
C = 32 
p = 3038820086802414604251869937287894280177159484142298205390043485090267683163423
j = 531432412335783608259183730437010937439245660086266581867668012322656101674228
coffi(p,j,C,R)

import random

# поиск простого p размера 356 бит
def PrimeP():
    p = randint(2**255, 2**256 -1)
    while p not in Primes():
        p = randint(2**255, 2**256 -1)
    print ("p =", p)
    return p

# поиск простого D
def PrimeD(p):
    for D in Primes():
        if D % 4 == 3 and kronecker_symbol(-D, p) == -1:
            print("D =", D)
            return D

# генерация суперсингулярной кривой
def SSCurve(p, D):
    if p == 2:
        print("y^2 + y = x^3")
        return
    
    elif p == 3:
        print("y^2 = x^3 - x")
        return
    
    # вычисление полинома Гильберта
    Gilbert = hilbert_class_polynomial(-D) 
    print("полином Гильберта:", Gilbert)
    # поиск корней
    j = int(Gilbert.roots()[0][0])
    print("корень полинома:", j)
    
    if D == 3:
        print("y^2 = x^3 - 1")
        return
    
    a = (27 * j) / 4 * (1728 - j)
    print("x^3 + {} * x + {}".format(a, a))

def main():    
    p = PrimeP()
    D = PrimeD(p)
    SSCurve(p, D)
    
if __name__ == "__main__":
    main()