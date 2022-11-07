include("jsat_mtk.jl")

#Make the Rigid Body
rb = make(RigidBody(1000*I(3),zeros(3),zeros(3)))
rb_ = structural_simplify(rb)

#Make the Reaction Wheels
r1 = ReactionWheel(0.25, 1, [ 1, 0, 0], 20)
r2 = ReactionWheel(0.35, 1.1, [ 0, 1, 0], 30)
r3 = ReactionWheel(0.45, 1.2, [ 0, 0, 1], 40)
rw = make([r1,r2,r3])
rw_ = structural_simplify(rw)

#Make the reaction wheel command
@named s = step(times = [1,3,5], durations = ones(3), values = ones(3))
s_ = structural_simplify(s)

#Fake the Thrusters
@named th = step(times = [Inf,Inf,Inf], durations = ones(3),values = ones(3))
th_ = structural_simplify(th)

eqs = [    
    s_.u ~ rw_.u
    th_.u ~ rb_.Te
    rw_.T ~ rb_.Ti
    rw_.H ~ rb_.Hi
]
sys = compose(ODESystem(eqs,t;name=:sys),rb_,rw_,s_,th_)
sys_ = structural_simplify(sys)
sol = test(sys_)
