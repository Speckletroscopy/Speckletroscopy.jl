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


fig = Figure()
ax = Axis(fig[1, 1], 
    xlabel = "time (ns)",
    ylabel = L"I(t)",
    xgridvisible = false,
    xticksmirrored=false,
    xminorticksvisible = true,
    xtickalign = true,
    xminortickalign = true,
    xlabelsize = 20,
    ygridvisible = false,
    yticksmirrored=true,
    yminorticksvisible = true,
    ytickalign = true,
    yminortickalign = true,
    ylabelsize = 20
)

e_field = eField(light_source,18)
time_vals = collect(0:0.001:1) # time in ns
efield_vals = map(t->eFieldT(t,e_field),time_vals)
intensity_vals = map(et->intensity(et),efield_vals)

lines!(time_vals, intensity_vals, label = "Instantaneous")
hlines!(ax,mean(intensity_vals), label = L"\mu", color = :red)
hlines!(ax,std(intensity_vals), label = L"\sigma", color= :green)
axislegend()

fig
