using DrWatson
@quickactivate :Speckletroscopy

light_source_dict = Dict(
    :n  => 100,
    :Em => [1.0,1.0,1.0],
    :νm => [456810,456813,456815], #GHz
    :σ  => 55.0,
    :fγ => 2.0e6
)
light_source = LightSource(light_source_dict)

e_field = eField(light_source,42)

time_vals = collect(0:0.001:20) # time in ns
τ_vals = collect(0:0.001:1) # offset in ns

efield_vals = map(t->eFieldT(t,e_field),time_vals)
intensity_vals = map(et->intensity(et),efield_vals)

ibar = mean(intensity_vals)

correlation_avg = ones(length(τ_vals))
correlation_std = ones(length(τ_vals))

sample_size = 10000
for i in 0:(length(τ_vals)-1)
    correlation = map(j->intensity_vals[j]*intensity_vals[j+i],1:sample_size)
    correlation_avg[i+1] = mean(correlation)
    correlation_std[i+1] = std(correlation)
end

g2_vals = correlation_avg/ibar^2
std_g2_vals = correlation_std/ibar^2
g2_calc = map(τ->g2Calc(τ,light_source),τ_vals)

cut_low = 100
mean_g2 = mean(g2_vals[cut_low:end])
mean_std_g2 = mean(std_g2_vals[cut_low:end])

plot(τ_vals, g2_vals, label = "\$g^{(2)}(\\tau)\$")
plot!(τ_vals, std_g2_vals, label = "\$\\delta g^{(2)}(\\tau)\$")
plot!(τ_vals, g2_calc, label = "\$g^{(2)}_{calc}(\\tau)\$", color = :black, linestyle = :dash)
hline!([sqrt(3)], label = "\$\\sqrt{3}\$", linestyle = :dot)
hline!([mean_g2], label = L"\mu", linestyle = :dot)
hline!([mean_std_g2], label = L"\sigma", linestyle = :dot)
xlabel!(L"\tau \mathrm{~(ns)}")
ylabel!("Normalized values (unitless)")
copy_ticks()

savefig(plotsdir("g2.svg"))
plot!()