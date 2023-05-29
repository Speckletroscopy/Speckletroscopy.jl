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


fig = Figure()
ax = Axis(fig[1, 1], 
    xlabel = L"\log(N)",
    ylabel = L"\frac{Var[I(t)I(t+τ)]}{\bar{I}}",
    xticksmirrored=false,
    xgridvisible = false,
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

e_field = eField(light_source,43)

time_vals = collect(0:0.001:20) # time in ns
τ_vals = collect(0:0.001:1) # offset in ns

efield_vals = map(t->eFieldT(t,e_field),time_vals)
intensity_vals = map(et->intensity(et),efield_vals)

# mean intensity
ibar = mean(intensity_vals)


log_sample_size = collect(1.0:0.1:4.0)
sample_size = 10 .^log_sample_size
sample_size = round.(Int,sample_size)

cut_low = 100 # low cut on τ to get beyond standard correlated regime

normalized_mean_correlation_var = ones(length(sample_size))
# iterate over sample size
for (k,nsamples) in enumerate(sample_size)
    # correlation_avg = ones(length(τ_vals))
    correlation_var = ones(length(τ_vals))
    for i in 0:(length(τ_vals)-1) # iterate over τ
        correlation = map(j->intensity_vals[j]*intensity_vals[j+i],1:nsamples)
        # correlation_avg[i+1] = mean(correlation)
        correlation_var[i+1] = var(correlation)
    end
    # take the mean value of the variance after low cut on τ
    normalized_mean_correlation_var[k] = mean(correlation_var[cut_low:end])/ibar^4
end

lines!(ax, log_sample_size, normalized_mean_correlation_var)

fig

# g2_vals = correlation_avg/ibar^2
# std_g2_vals = correlation_std/ibar^2
# g2_calc = map(τ->g2Calc(τ,light_source),τ_vals)

# cut_low = 100
# mean_g2 = mean(g2_vals[cut_low:end])
# mean_std_g2 = mean(std_g2_vals[cut_low:end])
