################################################################################
# File tree manipulation
################################################################################
#-------------------------------------------------------------------------------
"""
    resultsDir()

Returns the name of the directory where results are being stored
"""
function resultsDir()
    return results_directory
end

"""
    resultsDir(newdir::String)

Sets the directory where results are stored to newdir
"""
function resultsDir(newdir::String)
    global results_directory = newdir
    return results_directory
end
export resultsDir

#-------------------------------------------------------------------------------

"""
    plotsDir(sim::SpeckleSim)

Returns the path to the plots directory for the input simulation
"""
function plotsDir(sim::SpeckleSim)
    return plotsDir(sim.id)
end

function plotsDir(id::UUID)
    return joinpath(resultsDir(),string(id),"plots")
end

export plotsDir


#-------------------------------------------------------------------------------

"""
    dataDir(sim::SpeckleSim)

Returns the path to the data directory for the input simulation
"""
function dataDir(sim::SpeckleSim)
    return dataDir(sim.id)
end

function dataDir(id::UUID)
    return joinpath(resultsDir(),string(id),"data")
end

export dataDir

#-------------------------------------------------------------------------------
################################################################################
# Database management
################################################################################
"""
    mergeall(a::IndexedTable,b::IndexedTable)

Merges tables a and b. Unmatched columns are filled in with missing values.
"""
function mergeall(a::IndexedTable,b::IndexedTable)
    # make column names into sets for set operations
    colseta = Set(colnames(a))
    colsetb = Set(colnames(b))
    allcols = union(colseta, colsetb)
    commoncols = intersect(colseta,colsetb)
    # return just the regular merge if all columns are shared
    if commoncols == allcols
        return merge(a,b)
    end

    # get column names not common between the tables
    cols_anotb = setdiff(colseta,colsetb)
    cols_bnota = setdiff(colsetb,colseta)

    # merge parts of table with common columns
    commona = select(a,Tuple(commoncols))
    commonb = select(b,Tuple(commoncols))
    commontbl = merge(commona,commonb)
    # join uncommon columns into the table
    if length(cols_anotb) > 0
        onlya = select(a,(:id,cols_anotb...))
        commontbl = join(commontbl,onlya; how=:outer)
    end
    if length(cols_bnota) > 0
        onlyb = select(b,(:id,cols_bnota...))
        commontbl = join(commontbl,onlyb; how=:outer)
    end
    return commontbl
end

#-------------------------------------------------------------------------------
function tabulate(sim::SpeckleSim)
    simDict = Dict{Symbol,AbstractVector}()
    # iterate through fields in sim.params and expand arrays to individual columns
    for key in fieldnames(typeof(sim.params))
        if typeof(getfield(sim.params,key)) <: AbstractArray
            keyStr = string(key)
            for (i,val) in enumerate(getfield(sim.params,key))
                iKeyStr = string(keyStr,i)
                simDict[Symbol(iKeyStr)] = [val]
            end
        else
            simDict[key] = [getfield(sim.params,key)]
        end
    end
    simDict[:dt]  = [sim.dt]
    simDict[:id]  = [sim.id]
    simDict[:bst] = [sim.bs.t]
    simDict[:bsr] = [sim.bs.r]
    simDict[:elapsed] = [sim.elapsed]
    return DataFrame(simDict; pkey=:id)
end

function tabulate(simvec::Vector{T}) where {T<:SpeckleSim}
    out = tabulate(simvec[1])
    if length(simvec) == 1
        return out
    else
        for i=2:length(simvec)
            simtbl = tabulate(simvec[i])
            out = mergeall(out,simtbl)
        end
    end
    return out
end
export tabulate
#-------------------------------------------------------------------------------

################################################################################
# Save functions
################################################################################
function save(sim::SpeckleSim)
    # save the simulation to the results directory
    for i = 1:length(sim.readout)
        beamDict = Dict(:b1=>sim.readout[i].beam1, :b2=>sim.readout[i].beam2)
        beamTbl  = DataFrame(beamDict)
        beamName = string("counts",i,".csv")
        beamPath = joinpath(dataDir(sim),beamName)
        CSV.write(beamPath,beamTbl)

        corrDict = Dict{Symbol,Vector}()
        if typeof(sim.corr[i]) <: CorrelationVector
            corrDict[Symbol(sim.corr[i].n)] = sim.corr[i].data
        else
            for j=1:size(sim.corr[i].data)[2]
                corrDict[Symbol(j)] = sim.corr[i].data[:,j]
            end
        end
        corrTbl = DataFrame(corrDict)
        corrName = string("correlation",i,".csv")
        corrPath = joinpath(dataDir(sim),corrName)
        CSV.write(corrPath,corrTbl)
    end
    @info "Saved data for simulation id: $(sim.id)"
    return nothing
end

#-------------------------------------------------------------------------------

function save(sfft::SpeckleFFT)
    for (i,single) in enumerate(sfft.singles)
        singleDict = Dict(:frequency=>sfft.freqs, :coefficient=>single)
        singleTable = DataFrame(singleDict)
        singleName = string("singleft",i,".csv")
        singlePath = joinpath(dataDir(sfft.id),singleName)
        DataFrame.write(singleTable,singlePath)
    end

    sumFFTdict = Dict(:frequency=>sfft.freqs,:coefficient=>sfft.sumFFT)
    sumFFTtable = DataFrame(sumFFTdict)
    sumFFTname = "sumFFT.csv"
    sumFFTpath = joinpath(dataDir(sfft.id),sumFFTname)
    DataFrame.write(sumFFTtable,sumFFTpath)

    FFTsumdict = Dict(:frequency=>sfft.freqs,:coefficient=>sfft.FFTsum)
    FFTsumtable = DataFrame(FFTsumdict)
    FFTsumname = "FFTsum.csv"
    FFTsumpath = joinpath(dataDir(sfft.id),FFTsumname)
    DataFrame.write(FFTsumtable,FFTsumpath)
end

#-------------------------------------------------------------------------------
function save(plt::Plots.Plot,name::String,sim::SpeckleSim)
    figPath = joinpath(plotsDir(sim),name)
    savefig(plt,figPath)
    @info "Saved $name in $(plotsDir(sim))"
end
