function f(z, x=1.0, zeta=0.05)
   
    w = z[1]
    k0 = z[2]
    m0 = 10.0/9.0
    p = 3.0
    F = 10.0

    F/sqrt(4*k0*m0*w^2*x^(p+1)*zeta^2-2*k0*m0*w^2*x^(p+1)+k0^2*x^(2*p)+m0^2*w^4*x^2)
    
 end
 
 
 #
 # Derivada de f em relação a variável de projeto (x)
 # 
 function dfx(z ,x=1.0, zeta=0.05)

    w = z[1]
    k0 = z[2]
    m0 = 10.0/9.0
    p = 3.0
    F = 10.0
    
    -(F*(4*k0*m0*(p+1)*w^2*x^p*zeta^2+2*k0^2*p*x^(2*p-1)-2*k0*m0*(p+1)*w^2*x^p+2*m0^2*w^4*x))/(2*(4*k0*m0*w^2*x^(p+1)*zeta^2-2*k0*m0*w^2*x^(p+1)+k0^2*x^(2*p)+m0^2*w^4*x^2)^(3/2))    
    
 end
 
#
# References
#
function Referencias(x=1.0,zeta=0.05)

   @assert zeta==0.05
   @assert x==1.0

   E = 5.866107742312822

   Var = 6.869719177052227

   dE = -13.35903021749422

   dVar = -43.44923714683614

   return E,Var,dE,dVar

end