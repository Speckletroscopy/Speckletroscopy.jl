using DrWatson
@quickactivate :Speckletroscopy

light_source_dict = Dict(
    :n  => 100,
    :Em => [1.0,1.0,1.0],
    :νm => [456810,456813,456815], #GHz
    :σ  => 15.0,
    :fγ => 2.0e6
)
light_source = LightSource(light_source_dict)

e_field = eField(light_source,42)
time_vals = collect(0:0.001:1) # time in ns
efield_vals = map(t->eFieldT(t,e_field),time_vals)

plot(time_vals, abs.(efield_vals), label = false)

xlabel!("time (ns)")
ylabel!("E(t)")
copy_ticks()

plot!()
