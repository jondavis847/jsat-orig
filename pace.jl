

includet("lib.jl")

using Revise

#Mass properties
mass = 1506.27
inertia = [1529.097 -58.0883 -26.71023
            -58.0883  1400.682 83.51491
            -26.71023  83.51491 2320.778]
cg = [-1.311708,     -0.122079,    -0.0302493]

#Rigid Body
rb = make(RigidBody(inertia,zeros(3),zeros(3)))


#RWA
r1 = ReactionWheel(0.231, 1, [0.76604,0.64279,0], 50)
r2 = ReactionWheel(0.231, 1, [0,0.64279,0.76604], 50)
r3 = ReactionWheel(0.231, 1, [-0.76604,0.64279,0], 50)
r4 = ReactionWheel(0.231,1,[0,0.64279,-0.76604], 50)
rw = make([r1,r2,r3,r4])

@named rw_command = Step3(times = [1.,4.,7.])
@named rw_sys = ODESystem([
        connect(rw_command.output,rw.input_u)
    ], t, systems = [rw,rw_command] )
