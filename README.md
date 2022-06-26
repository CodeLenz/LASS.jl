# LASS
Local Averaged Stratified Sampling

Valentini, F., Silva, O.M., Torii, A.J. et al. Local averaged stratified sampling method. J Braz. Soc. Mech. Sci. Eng. 44, 294 (2022). https://doi.org/10.1007/s40430-022-03589-6

```julia
using LASS
using  Distributions, Statistics

# Lets load QuadGK to evaluate some reference values
load QuadGK

#
# Lets define a function of z. The argument must z
# be a vector
#
f(z::Vector,x=1.0) = (z[1]-5*x)^2 + z[1]*x;

# and assume that z is N(5,1.0)
pz = Normal(5.0,1.0);

# Lets obtain some reference values. 
# The expected value E[f] is 
pE(z)=f([z])*pdf(pz,z);
E = quadgk(pE,-Inf,Inf)[1];

# And the variance is
pV(z)=(f([z])-E)^2*pdf(pz,z);
Var = quadgk(pV,-Inf,Inf)[1];

#
#
# Now lets use the Local Averaged Stratified Sampling Method
#
# generate a large number of realizations. Each realization must 
# be stored in a column of a nz x n matriz.
n = 1_000_000;
distrib = rand(pz,1,n);

# Generate the bins
Nb = [50];
bins = Generate_bins(distrib,Nb);

# Evaluate E and Var using only  50 evaluations
EL,VarL = Lass(bins,f);

println("E[f] $E - reference")
println("E[f] $EL obtained with $(Nb[1]) evaluations")
println()
println("Var[f] $Var - reference")
println("Var[f] $VarL obtained with $(Nb[1]) evaluations")


#
# The same thing can be done to evaluate E[df] and Var[df]
# where df is the derivative of f w.r.t a deterministic 
# design variable x
#
dfx(z::Vector,x=1.0) = -10*(z[1]-5*x) + z[1];

# Lets obtain some reference values. 
# The derivative of the expected value E[f] is 
pdE(z)=dfx([z])*pdf(pz,z);
dE = quadgk(pdE,-Inf,Inf)[1];

# And the derivative of the variance is
# d(f-E)^2 ==> 2(f-E)*(df-dE)
#
pdV(z)=2*(f([z])-E)*(dfx([z])-dE)*pdf(pz,z);
dVar = quadgk(pdV,-Inf,Inf)[1];

# Evaluate dE and dVar using only  50 evaluations
dEL,dVarL = dLass(bins,f,dfx,1);


println("dE[f] $dE - reference")
println("dE[f] $dEL obtained with $(Nb[1]) evaluations")
println()
println("dVar[f] $dVar - reference")
println("dVar[f] $dVarL obtained with $(Nb[1]) evaluations")
```
The results for this problem are

```julia
E[f] 6.000000000000621 - reference
E[f] 5.997049691574353 obtained with 50 evaluations
```
``` julia
Var[f] 3.0000000000000417 - reference
Var[f] 2.98478094281721 obtained with 50 evaluations
```
and 

``` julia
dE[f] 4.999999999999871 - reference
dE[f] [4.9985604747676415] obtained with 50 evaluations
``` 
``` julia
dVar[f] 5.999999999960915 - reference
dVar[f] [5.9715595256430305] obtained with 50 evaluations
``` 



