function sensors!(dx, x, p, t)
    tam!(dx,x,p,t)
    iru!(dx,x,p,t)
    startrackers!(dx,x,p,t)
end

function sensors_cb!(integrator)
end