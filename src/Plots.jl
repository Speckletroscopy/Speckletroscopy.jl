function copy_ticks(sp::Plots.Subplot)
	ptx = twinx(sp)
	plot!(ptx,xlims=xlims(sp),ylims=ylims(sp),xformatter=_->"",yformatter=_->"")
	pty = twiny(sp)
	plot!(pty,xlims=xlims(sp),ylims=ylims(sp),xformatter=_->"",yformatter=_->"")
end

"""
	copy_ticks(plt::Plots.Plot = current())
	
Call `copy_ticks()` in same cell as plot function in order to get ticks on all sides
"""
copy_ticks(plt::Plots.Plot = current()) = copy_ticks(plt[1])

export copy_ticks


# #-------------------------------------------------------------------------------
# """
#     γIntensityPlot(times::Array,γint::Array,prefix::String)

# Plots average photon counts and saves to file. Returns relative path to resulting plot.
# """
# function γIntensityPlot(times::Array,γint::Array,prefix::String)
#     γplotName = string(prefix,"avg-photons-vs-time.svg")
#     γplot = plot(times,γint,label=false)
#     xlabel!("time (ns)")
#     ylabel!("Avg photon counts")
#     title!("Average Photon Counts vs Time")
#     savefig(γplot,γplotName)
#     @info "Saved $γplotName"
#     return γplotName
# end

# export γIntensityPlot
# #-------------------------------------------------------------------------------

# """
#     γCountPlot(times::Array,γint::Array,prefix::String)

# Plots 'measured' photon counts and saves to file. Returns relative path to resulting plot.
# """
# function γCountPlot(times::Array,γcount::Array, γint::Array,prefix::String)
#     γplotName = string(prefix,"photon-counts-vs-time.svg")
#     γplot = bar(times,γcount,label="Simulated Counts")
#     sig1 = sqrt.(γint)
#     sig1low = map(x -> x < 1 ? x : sqrt(x),γint)
#     plot!(γplot,times,γint,label = "Average Counts"; ribbon = (sig1low,sig1))
#     xlabel!("time (ns)")
#     ylabel!("photon counts")
#     title!("Photon Counts vs Time")
#     savefig(γplot,γplotName)
#     @info "Saved $γplotName"
#     return γplotName
# end

# export γCountPlot
# """
#     γCorrTimePlot(timeVec::Vector, g2Vec::Vector, params::Dict, prefix::String)

# Saves plot of photon correlation time series
# """
# function γCorrTimePlot(timeDF::IndexedTable, params::Dict, prefix::String)
#     out = []
#     γplotName = string(prefix,"time-domain-photon-correlation.svg")
#     inzcorr  = timeDF[!,:corr1] .> 0
#     nzcorr = timeDF[!,:corr1][inzcorr]
#     nztime = timeDF[!,:time][inzcorr]
# 	γplot = plot(nztime,nzcorr,label = false)
#     xlabel!(L"\tau \textrm{ (ns)}")
# 	ylabel!("\$g^{(2)}(\\tau)\$")
#     title!("Photon Correlations vs Time Offset")
#     savefig(γplot,γplotName)
#     push!(out,γplotName)
#     @info "Saved $(γplotName)"
#     if params[:repeat] > 1
#         γplotName = string(prefix,"time-domain-photon-correlation-sum.svg")
#         nzcorr = timeDF[!,:sum][inzcorr]
#         γplot = plot(nztime,nzcorr,label = false)
#         xlabel!(L"\tau \textrm{ (ns)}")
#         ylabel!("\$\\sum_{i=1}^{$(params[:repeat])}g_i^{(2)}(\\tau)\$")
#         title!("Sum of Photon Correlations vs Time Offset")
#         savefig(γplot,γplotName)
#         push!(out,γplotName)
#         @info "Saved $(γplotName)"
#     end
#     return out
# end

# function γCorrTimePlot(corr::Vector{T},params::SpeckleParams) where {T<:Number}
#     times = collect(0:(length(corr)-1))*params.resolution
#     out = plot(times,corr,label="simulation")
#     plot!(out,times,map(time->g2Calc(time,params),times),label = "calculated")
#     xlabel!(L"\tau \textrm{ (ns)}")
# 	ylabel!("\$g^{(2)}(\\tau)\$")
#     title!("Photon Correlations vs Time Offset")
#     return out
# end

# function γCorrTimePlot(sim::SpeckleSim; ind::Union{String,Int} = "avg")
#     if ind isa Int
#         @assert 0 < ind < length(sim) "Index out of bounds"
#     elseif ind == "avg"
#         cavg = +(sim.corr...)
#         return γCorrTimePlot(cavg.data/length(sim.corr),sim.params)

#     end
#     @assert false "Invalid ind in γCorrTimePlot"
# end

# export γCorrTimePlot
# #-------------------------------------------------------------------------------

# """
#     γCorrFreqPlot(freqVec::Vector, fftVec::Vector, params::Dict, prefix::String)

# Saves plot of photon correlation Fourier transform
# """
# function γCorrFreqPlot(freqDF::IndexedTable, params::Dict, prefix::String)
#     out = []
#     γplotName = string(prefix,"frequency-domain-photon-correlation.svg")
# 	γplot = plot(freqDF[!,:freq],freqDF[!,:corr1],label = false)
# 	xlabel!("frequency (GHz)")
# 	ylabel!("\$\\hat{g}^{(2)}(\\nu)\$")
# 	title!("Photon Correlations vs Frequency")
#     if length(params[:νm]) > 1
#         shifts = map(x->abs(x[2]-x[1]),subsets(params[:νm],2))
#         shiftname = length(params[:νm]) > 2 ? "frequency shifts" : "frequency shift"
#         vline!(shifts,label=shiftname,ls=:dash)
#     end
#     savefig(γplot,γplotName)
#     push!(out,γplotName)
#     @info "Saved $(γplotName)"
#     if params[:repeat] > 1
#         γplotName = string(prefix,"frequency-domain-photon-correlation-sum.svg")
#         γplot = plot(freqDF[!,:freq],freqDF[!,:sum],label = false)
#         xlabel!("frequency (GHz)")
#         ylabel!("\$\\sum_{i=1}^{$(params[:repeat])}\\hat{g}^{(2)}(\\nu)\$")
#         title!("Sum of Photon Correlations vs Frequency")
#         if length(params[:νm]) > 1
#             shifts = map(x->abs(x[2]-x[1]),subsets(params[:νm],2))
#             shiftname = length(params[:νm]) > 2 ? "frequency shifts" : "frequency shift"
#             vline!(shifts,label=shiftname,ls=:dash)
#         end
#         savefig(γplot,γplotName)
#         push!(out,γplotName)
#         @info "Saved $(γplotName)"
#     end
#     return out
# end

# export γCorrFreqPlot
# #-------------------------------------------------------------------------------
