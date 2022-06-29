module SimV1
using Speckletroscopy

function runsim()

    
    # νHα1 = [456813] #GHz
    νHα2 = [456810,456813,456815] #GHz
    paramDict = Dict(
                     :n    => [10,20,30],#,40,80],#,100,160], # number of atoms
                    :νm   => [νHα2], # line frequencies in GHz
                    :Em   => ["ones"], # relative line magnitudes
                    :σ    => [20.0], # Doppler broadening in GHz
                    :fγ   => [2.0e6],#,"shot10%","shot50%",10.0,1.0,0.16], # mean photon count rate in GHz
                    :deadtime   => [0.0], # detector deadtime in nanoseconds
                    :resolution => [0.010],#,0.10], # detector resolution in nanoseconds
                    :jitter     => [0.015], # detector timing jitter in nanoseconds 
                    :efficiency => [0.9], # detector efficiency
                    :darkcounts => [1.0e-8], # detector dark count rate in GHz
                    :duration   => [20.0], # duration of each correlation measurement in nanoseconds
                    :window     => ["halfwindow"], # time over which to average correlations in nanoseconds
                    :repeat     => [100], # number of times to repeat correlation measurement
                    :reinstance => [true] # control whether or not frequencies and phases should be reinstanced between measurements
                    )

    test = SpeckleParamsVector(paramDict)
    println(test[1])
    Speckletroscopy.run(test[1])
end
export runsim

end


using .SimV1
SimV1.runsim()
