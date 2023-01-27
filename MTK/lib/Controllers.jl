struct Step <: Controller
    time::Float64
    duration::Float64
    height::Float64
end

function step(x,t,p)
        dx = ((t > p.time) & (t < p.time + p.duration)) ? p.height : 0
end