include("jsat_mtk.jl")

#Make the Rigid Body
rb = make(RigidBody(1000*I(3),zeros(3),zeros(3)))

#Make the Reaction Wheels
r1 = ReactionWheel(0.25, 1, [ 1, 0, 0], 20)
r2 = ReactionWheel(0.35, 1.1, [ 0, 1, 0], 30)
r3 = ReactionWheel(0.45, 1.2, [ 0, 0, 1], 40)
rw = make([r1,r2,r3])

#Make the reaction wheel command
@named s = step(times = [1,3,5], durations = ones(3), values = ones(3))

#Fake the Thrusters
@named th = step(times = [Inf,Inf,Inf], durations = ones(3),values = ones(3))

eqs = [    
    s.u ~ rw.u
    th.u ~ rb.Te
    rw.T ~ rb.Ti
    rw.H ~ rb.Hi
]
sys = compose(ODESystem(eqs,t;name=:sys),rb,rw,s,th)
sys_ = structural_simplify(sys)
sol = test(sys_)
