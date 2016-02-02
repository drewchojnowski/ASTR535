import math

def combinemags(m1=0,m2=0):
    f1=math.pow(10,-0.4*m1)
    f2=math.pow(10,-0.4*m2)
    fsum=f1+f2
    msum=-2.5*math.log10(fsum)
    print 'f1 = ',f1
    print 'f2 = ',f2
    print 'f_tot = ',fsum
    print 'm_tot =  ',msum

    return msum
