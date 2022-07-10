using Plots
using DifferentialEquations
using ModelingToolkit

Plots.gr(lw = 2,fmt = :png)

# Convenience functions
hill(x, k) = x / (x + k)
hill(x, k, n) = hill(x^n, k^n)
reducedmm(x, k) = x / k

@parameters v_0 v_m1 v_m2 v_m3 Km_1 Km_2 Km_3
@variables t S1(t) S2(t) S3(t) v1(t) v2(t) v3(t)
D = Differential(t)

eqsFull = [ v1 ~ v_m1 * hill(S1, Km_1),
            v2 ~ v_m2 * hill(S2, Km_2),
            v3 ~ v_m3 * hill(S3, Km_3),
            D(S1) ~ v_0 - v1,
            D(S2) ~ v1 - v2,
            D(S3) ~ v2 - v3]

@named fullsys = ODESystem(eqsFull)
fullSys = structural_simplify(fullsys)

eqsRe =   [ v1 ~ v_m1 * S1 / Km_1,
            v2 ~ v_m2 * S2 / Km_2,
            v3 ~ v_m3 * S3 / Km_3,
            D(S1) ~ v_0 - v1,
            D(S2) ~ v1 - v2,
            D(S3) ~ v2 - v3]
@named reSys = ODESystem(eqsRe)
reSys = structural_simplify(reSys)

u0 = [S1 => 6.0, S2 => 4.0, S3 => 4.0]
tend = 10.0
params = [v_0 => 2, v_m1 => 9, v_m2 => 12, v_m3 => 15, Km_1 => 1.0, Km_2 => 0.4, Km_3 => 3.0]

sol1 = solve(ODEProblem(fullSys, u0, tend, params))

plot(sol1, ylims=(0.0, 1.0),
     title="Problem 3.7.5 (1) for (s1, s2, s3) = (6.0, 4.0, 4.0)",
     xlabel="Time (arbitrary units)",
     ylabel="Concentration (arbitrary units)")

     