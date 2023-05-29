using DrWatson
@quickactivate :Speckletroscopy

νHα2 = [456810,456813,456815] #GHz

# dictionary of all parameters we wish to test
paramDict = Dict(
    :n    => [10,20,30,40,80,100,160], # number of atoms
    :νm   => [νHα2], # line frequencies in GHz
    :Em   => ["ones"], # relative line magnitudes
    :σ    => [10.0], # Doppler broadening in GHz
    :fγ   => [2.0e6],#,"shot10%","shot50%",10.0,1.0,0.16], # mean photon count rate in GHz
    :deadtime   => [0.0], # detector deadtime in nanoseconds
    :resolution => [0.01],#,0.10], # detector resolution in nanoseconds
    :jitter     => [0.015], # detector timing jitter in nanoseconds 
    :efficiency => [0.9], # detector efficiency
    :darkcounts => [1.0e-8], # detector dark count rate in GHz
    :duration   => [40.0], # duration of each correlation measurement in nanoseconds
    :window     => ["halfwindow"], # time over which to average correlations in nanoseconds
    :repeat     => [100], # number of times to repeat correlation measurement
    :reinstance => [true] # control whether or not frequencies and phases should be reinstanced between measurements
)

# Take all combinations of parameters above and place into individual simulation dictionaries
paramVec = SpeckleParamsVector(paramDict)

nn = 1

instance_natoms100 = SpeckleInstance(paramVec[nn])

avg_counts = instance_natoms100.:γint.γvec

nbar_per_bin = instance_natoms100.:γint.nbar/length(avg_counts)

window = 3000
max_offset = 1000
nbar_per_window = window*nbar_per_bin

bigG2_pairs = map(0:max_offset) do offset
    i_t = view(avg_counts, 1:window)
    i_t_tau = view(avg_counts, (1+offset):(window+offset))
    prod_ii_ttau = i_t .* i_t_tau
    mean_prod = mean(prod_ii_ttau)
    std_prod = std(prod_ii_ttau)
    [mean_prod, std_prod] ./ nbar_per_bin^2
end

bigG2_mean = map(x->x[1], bigG2_pairs)
bigG2_std = map(x->x[2], bigG2_pairs)

fig = Figure()
ax = Axis(fig[1, 1], 
    xlabel = L"\tau (ps)",
    ylabel = L"g^{(2)}(\tau)",
    xgridvisible = false,
    xticksmirrored=true,
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

τvals = collect(0:max_offset)

low_cut = 10
high_cut = 0

lines!(τvals[low_cut:(end-high_cut)],bigG2_mean[low_cut:(end-high_cut)], label=L"g^{(2)}(\tau)")
lines!(τvals[low_cut:(end-high_cut)],bigG2_std[low_cut:(end-high_cut)], label=L"$Δg^{(2)}(\tau)$") 
lines!(τvals[low_cut:(end-high_cut)], map(τ->g2Calc(τ,paramVec[nn]),τvals[low_cut:(end-high_cut)]), label=L"$g^{(2)}_{calc}(\tau)$") 




hlines!(ax,sqrt(3.0),color=:red, label = L"\sqrt{3}")

axislegend()

Makie.save(plotsdir("speckle_constant_noise.svg"),fig)

# Display the plot
fig