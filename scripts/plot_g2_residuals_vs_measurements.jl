using Distributed
@everywhere using DrWatson

@everywhere begin
    @quickactivate :Speckletroscopy
    using SharedArrays
end

light_source_dict = Dict(
    :n  => 100,
    :Em => [1.0,1.0,1.0],
    :νm => [456810,456813,456815], #GHz
    :σ  => 55.0,
    :fγ => 2.0e6
)
light_source = LightSource(light_source_dict)

# mean intensity
ibar = light_source.n*sum(x->real(x*conj(x)),light_source.Em)


time_vals = collect(0:0.001:20) # time in ns
τ_vals = collect(0:0.001:1) # offset in ns

g2_calc = map(τ->g2Calc(τ,light_source),τ_vals)


# define sample sizes
log_sample_size = collect(1.0:0.1:3.0)
sample_size = 10 .^log_sample_size
sample_size = round.(Int,sample_size)

cut_low = 100 # low cut on τ to get beyond standard correlated regime

ninstances = 500 
normalized_mean_correlation_residual_std_array = SharedArray{Float64}(length(sample_size),ninstances)

@info "Calculating $(ninstances) instances"
@sync @distributed for inst in 1:ninstances
    @info "Instance $(inst)"
    e_field = eField(light_source)
    efield_vals = map(t->eFieldT(t,e_field),time_vals)
    intensity_vals = map(et->intensity(et),efield_vals)

    normalized_mean_correlation_residuals_std = ones(length(sample_size))
    # iterate over sample size
    for (k,nsamples) in enumerate(sample_size)
        correlation_avg = ones(length(τ_vals))
        for i in 0:(length(τ_vals)-1) # iterate over τ
            correlation = map(j->intensity_vals[j]*intensity_vals[j+i],1:nsamples)
            correlation_avg[i+1] = mean(correlation)
        end
        # calculate residuals compared to g2_calc
        residuals = correlation_avg/ibar^2 .- g2_calc
        normalized_mean_correlation_residuals_std[k]= std(residuals)
    end
    normalized_mean_correlation_residual_std_array[:,inst] = normalized_mean_correlation_residuals_std
end

g2_residual_mean_by_samples = map(k->mean(normalized_mean_correlation_residual_std_array[k,:]),1:length(sample_size))
g2_residual_std_by_samples= map(k->std(normalized_mean_correlation_residual_std_array[k,:]),1:length(sample_size))

plot(log_sample_size, g2_residual_mean_by_samples, ribbon = g2_residual_std_by_samples, label=false)
hline!([1/light_source.n], label = L"\frac{1}{n}", color=:red)
xlabel!(L"\log_{10}(\mathrm{sample~size})")
ylabel!(L"\mathrm{residuals~of~}g^{(2)}(\tau)")
copy_ticks()

savefig(plotsdir("g2_residual_scaling.svg"))

plot!()
