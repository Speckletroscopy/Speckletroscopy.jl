### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ e40b2d2a-f889-11ec-3c35-3fe8ddb46a6c
using Pkg

# ╔═╡ d70cfb4b-c9a6-4161-bddf-79c73b5ee70a
Pkg.activate(".")

# ╔═╡ 7babf9f0-87cc-4448-84a9-48daae1a83e4
using Speckletroscopy;

# ╔═╡ 694e95ed-6a86-49b2-ae3c-683ae7fd4e14
νHα2 = [456810,456813,456815] #GHz

# ╔═╡ 7b1d7ad9-4e2b-401f-825a-204f62aa82f0
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
				:duration   => [20.0], # duration of each correlation measurement in nanoseconds
				:window     => ["halfwindow"], # time over which to average correlations in nanoseconds
				:repeat     => [100], # number of times to repeat correlation measurement
				:reinstance => [true] # control whether or not frequencies and phases should be reinstanced between measurements
				)

# ╔═╡ 44160a16-9d42-41ea-8f2e-e0c3b8e4fa06
begin
	shifts = convert.(Float64,lineShifts(νHα2))
	nyquistSampling = 2*shifts
	nyquistRes = 1.0 ./ nyquistSampling
end

# ╔═╡ 3c029b45-1c9a-42ab-bc8d-6efb741587fd
shifts

# ╔═╡ 3c371f49-6c13-4931-99d0-b985b4374a6a
paramVec = SpeckleParamsVector(paramDict)

# ╔═╡ b13eec73-ed66-44e0-804d-574db45e0ed2
run_vec = Speckletroscopy.run(paramVec[7])

# ╔═╡ 2c2db19e-3ba5-43c2-90ae-8ada6dae0631
this_sim = run_vec

# ╔═╡ 4f146129-8db5-49c3-b7fb-8f501e23654c
sumcorr = sum(this_sim.corr)

# ╔═╡ 8f50ee7b-f121-488f-bfff-3fbd520ec9f9
τ = this_sim.params.resolution*collect(0:length(sumcorr.data)-1)

# ╔═╡ 88cc30eb-cfc9-4c17-8a01-04806aa59338
let
	# percentStartCut = 40
	# percentEndCut = 50
	
	# percentStartCut /=100
	# percentEndCut /= 100
	# lengthTot = length(τ)
	# cutStart = convert(Int64,floor(lengthTot*percentStartCut))
	# cutEnd = convert(Int64,floor(lengthTot*percentEndCut))

	# plot()
	# plot!(τ[1+cutStart:end-cutEnd],sumcorr.data[1+cutStart:end-cutEnd]/length(this_sim.corr),label="sim")
	# plot!(τ[1+cutStart:end-cutEnd],map(t->g2Calc(t,this_sim.params),τ[1+cutStart:end-cutEnd]), label="model")
end

# ╔═╡ 7866fc91-0cf4-4ffe-9445-dc83c1c31a2f
1/this_sim.params.σ

# ╔═╡ 1864a410-fee7-497c-9fb4-8c794c6461bd
# generate an instance of an electric field for testing
begin
	instance = SpeckleInstance(this_sim.params)
	ef = eField(instance.source,3)
end

# ╔═╡ f285faa8-3485-4079-84b8-994b4c0cd1a4
let
histogram(collect(ivec(transpose(ef.νn) .+ Δm(ef.source))),nbins=20)
histogram!(ef.νn,nbins=20)
end

# ╔═╡ 7aaa5eb2-124f-4f6e-af62-9d5d6b5d9b31
histogram(ef.νn)

# ╔═╡ 6713a4dd-9bb1-403c-8723-6c826a52eb18


# ╔═╡ f56f6a4f-0162-421f-affc-21dd462c6647
begin
	# Δt = [0.01]#,0.05,0.1,0.166667,0.25]
	# tmax = 2.0
	# tseries = map(dt->collect(0:dt:tmax),Δt)
	# tdata = collect(zip(Δt,tseries))
	# meanIntensitySeries = map(tdat->
	# map(tt->meanIntensity(tt,tdat[1],ef),tdat[2]),tdata)
	# ttseries = collect(0:0.001:tmax)
	# intensitySeries = map(tt->intensity(tt,ef),ttseries)
end

# ╔═╡ 6eb5b6c9-f2c9-4b0f-8ea1-ac413d3bb173
# plot some intensities at different timescales
# begin
# 	plot()
# 	for i=1:length(tseries)
# 		plot!(tseries[i]*paramDict[:σ][1] ,meanIntensitySeries[i]/90)
# 	end
# 	plot!(ttseries*paramDict[:σ][1],intensitySeries/90)
# 	# vline!(tseries[1]*paramDict[:σ][1])
# 	plot!(title="Mean Intensity")
# end

# ╔═╡ 73e6073b-5c3c-4db2-b63c-f242e45dfec2
begin

	# tfine = ttseries
	# intensityFine = intensitySeries
	# Δtfine = 0.001
	# intensityFineFT = meanFFT(intensityFine,(1,length(intensityFine)))
	# intensityFineFT = abs.(intensityFineFT)
	# intensityFineFT = convert.(Float64,intensityFineFT)
	# maxIntensityFineFT = max(intensityFineFT...)
	# scaledIntensityFineFT = intensityFineFT/maxIntensityFineFT

	# fineFreq = fftFreq(Δtfine,(1,length(intensityFine)))

	# fineFreq,fineCoeff = fftPositiveFreq(fineFreq,scaledIntensityFineFT)
end


# ╔═╡ 4c95ac68-3470-4a91-ab0f-ccd1f21a2df0
let
	# percentStartCut = 0
	# percentEndCut = 90
	
	# percentStartCut /=100
	# percentEndCut /= 100
	# lengthTot = length(fineFreq)
	# cutStart = convert(Int64,floor(lengthTot*percentStartCut))
	# cutEnd = convert(Int64,floor(lengthTot*percentEndCut))

	
	# plot()
	# plot!(fineFreq[1+cutStart:end-cutEnd],fineCoeff[1+cutStart:end-cutEnd])
	# vline!([shifts...])
end

# ╔═╡ Cell order:
# ╠═e40b2d2a-f889-11ec-3c35-3fe8ddb46a6c
# ╠═d70cfb4b-c9a6-4161-bddf-79c73b5ee70a
# ╠═7babf9f0-87cc-4448-84a9-48daae1a83e4
# ╠═694e95ed-6a86-49b2-ae3c-683ae7fd4e14
# ╠═7b1d7ad9-4e2b-401f-825a-204f62aa82f0
# ╠═44160a16-9d42-41ea-8f2e-e0c3b8e4fa06
# ╠═3c029b45-1c9a-42ab-bc8d-6efb741587fd
# ╠═3c371f49-6c13-4931-99d0-b985b4374a6a
# ╠═b13eec73-ed66-44e0-804d-574db45e0ed2
# ╠═2c2db19e-3ba5-43c2-90ae-8ada6dae0631
# ╠═4f146129-8db5-49c3-b7fb-8f501e23654c
# ╠═8f50ee7b-f121-488f-bfff-3fbd520ec9f9
# ╠═88cc30eb-cfc9-4c17-8a01-04806aa59338
# ╠═7866fc91-0cf4-4ffe-9445-dc83c1c31a2f
# ╠═1864a410-fee7-497c-9fb4-8c794c6461bd
# ╠═f285faa8-3485-4079-84b8-994b4c0cd1a4
# ╠═7aaa5eb2-124f-4f6e-af62-9d5d6b5d9b31
# ╠═6713a4dd-9bb1-403c-8723-6c826a52eb18
# ╠═f56f6a4f-0162-421f-affc-21dd462c6647
# ╠═6eb5b6c9-f2c9-4b0f-8ea1-ac413d3bb173
# ╠═73e6073b-5c3c-4db2-b63c-f242e45dfec2
# ╠═4c95ac68-3470-4a91-ab0f-ccd1f21a2df0
