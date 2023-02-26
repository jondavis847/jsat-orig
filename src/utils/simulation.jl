
function simulate(ic::InitialConditions,p,tspan; nruns = 1, dt = 0.1)
    x = defineModel(ic)
    simulate(x,p,tspan;nruns = nruns, dt = dt)
end


function simulate(x, p, tspan; nruns=1, dt=0.1)
    #Run nominal 
    xn,pn = initModel(x,p)
    nominal_prob = ODEProblem(model!, xn, tspan, pn, callback=model_cb)

    sol = solve(nominal_prob, Tsit5(), saveat=dt)
    sol = convertSol(sol.t, sol.u)

    #Run monte carlo
    if nruns > 1
        p.config.montecarlo = true
        mc_prob = ODEProblem(model!, x, tspan, p, callback=model_cb)
        eprob = EnsembleProblem(mc_prob, prob_func=prob_func, output_func=output_func)
        sols = solve(eprob, Tsit5(), trajectories=nruns, saveat=dt)
        sol_out = (nominal=sol,dispersed=sols)
    else
        sol_out = sol
    end
    return sol_out
end

function prob_func(prob, i, repeat)
    print("Simulating run $i\n")    
    xm,pm = initModel(prob.u0,prob.p)
    remake(prob, p=pm, u0 = xm)
end
function output_func(sol, i)
    sol = convertSol(sol.t, sol.u)
    return (sol, false)
end

function convertSol(t, u)
    u_out = digdeeper(u)
    return (t=t, u=u_out)
end
function digdeeper(in)
    d = Dict()
    for k in keys(in[1])
        if in[1][k] isa ComponentVector
            merge!(d, Dict(k => digdeeper(map(x -> x[k], in))))
        else
            merge!(d, Dict(k => map(x -> x[k], in)))
        end
    end
    return NamedTuple(d)
end


#Seems OBE but leaving for future reference in case it's not
function recursive_sim(prob, tspan, interval)
    n = (tspan[2] - tspan[1]) / interval
    t = Float64[tspan[1]]
    u = [prob.u0]
    t0 = tspan[1]
    tf = t0 + interval
    uf = prob.u0
    for i in 1:n
        prob = remake(prob, u0=uf, tspan=(t0, tf))
        sol = solve(prob, Tsit5(), callback=model_cb)
        t0 = sol.t[end]
        tf = t0 + interval
        push!(t, t0)
        uf = sol.u[end]
        push!(u, uf)
    end
    return convertSol(t, u)
end


