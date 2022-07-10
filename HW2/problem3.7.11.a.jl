using Plots
using DifferentialEquations
using ModelingToolkit

Plots.gr(lw = 2,fmt = :png)

# Convenience functions


# Y ~ (K * x / Kt * ( 1 + x/Kt) + x/Kr*(1 +x/Kr)) / (K*(((1 + x/Kt)^2) +(1 + x/Kr)^2))
Y1(x) = (500 * x / 1000 * ( 1 + x/1000) + x/1*(1 +x/1)) / (500*(((1 + x/1)^2) +(1 + x/1)^2))
Y2(x) = (1000 * x / 1000 * ( 1 + x/1000) + x/1*(1 +x/1)) / (1000*(((1 + x/1)^2) +(1 + x/1)^2))
Y3(x) = (1500 * x / 1000 * ( 1 + x/1000) + x/1*(1 +x/1)) / (1500*(((1 + x/1)^2) +(1 + x/1)^2))
#params1 = [Kt => 1000,Kr => 1,K => 500]
#params2 = [Kt => 1000,Kr => 1,K => 1000]
#params3 = [Kt => 1000,Kr => 1,K => 1500]
plot(x->Y1(x), 0.0, 100.0, labels = "K = 500")
plot!(x->Y2(x),0.0,100.0,labels="K = 1000" )
plot!(x->Y3(x),0.0,100.0 labels="K = 1500" )


