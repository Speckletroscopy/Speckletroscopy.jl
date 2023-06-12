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

time_vals = collect(0:0.001:20) # time in ns
τ_vals = collect(0:0.001:1) # offset in ns



# define sample sizes
log_sample_size = collect(1.0:0.1:3.0)
sample_size = 10 .^log_sample_size
sample_size = round.(Int,sample_size)

# mean intensity
ibar = light_source.n*sum(x->real(x*conj(x)),light_source.Em)

cut_low = 100 # low cut on τ to get beyond standard correlated regime

ninstances = 100 

correlation_var_array = SharedArray{Float64}(length(sample_size),ninstances)

@info "Calculating $(ninstances) instances"
@sync @distributed for inst in 1:ninstances
    @info "Instance $(inst)"
    e_field = eField(light_source)
    efield_vals = map(t->eFieldT(t,e_field),time_vals)
    intensity_vals = map(et->intensity(et),efield_vals)


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
    correlation_var_array[:,inst] = normalized_mean_correlation_var
end

mean_var_correlation = map(k->mean(correlation_var_array[k,:]),1:length(sample_size))
std_var_correlation = map(k->std(correlation_var_array[k,:]),1:length(sample_size))

plot(log_sample_size, mean_var_correlation, ribbon = std_var_correlation, label=false)
xlabel!(L"\log_{10}(\mathrm{sample~size})")
ylabel!(L"\mathrm{Normalized~variance~of~}I(t)I(t+\tau)")
copy_ticks()

savefig(plotsdir("classical_variance_scaling.svg"))

plot!()
