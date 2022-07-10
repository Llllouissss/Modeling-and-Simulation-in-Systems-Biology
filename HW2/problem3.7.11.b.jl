using Plots
using DifferentialEquations
using ModelingToolkit

Plots.gr(lw = 2,fmt = :png)

# Convenience functions


# Y ~ (K * x / Kt * ( 1 + x/Kt) + x/Kr*(1 +x/Kr)) / (K*(((1 + x/Kt)^2) +(1 + x/Kr)^2))
Y1(x) = (0*x/1000*(1 + x/1000) + x*(1+x)) / (0*(1 + x/1000)^2 +(1 + x/1)^2)
Y2(x) = (x/10)/(1+x/10)

#params1 = [Kt => 1000,Kr => 1,K => 500]
#params2 = [Kt => 1000,Kr => 1,K => 1000]
#params3 = [Kt => 1000,Kr => 1,K => 1500]
plot(x->Y1(x), 0.0, 200.0, labels = "K = 0")
plot!(x->Y2(x),0.0,200.0,labels="KR = KT = 10" )
# plot!(x->Y3(x),0.0,200.0,labels="K = 1500" )


