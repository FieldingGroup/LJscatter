#!/usr/bin/env julia

import Pkg
Pkg.activate(".")

using HDF5
using Plots
using DelimitedFiles

h5path = raw"bases/allData.h5"
fileName = "bases/basisConcProfiles.h5"
amplitudes = [0.,0.1,0.2,.5,1.,2.,5.,10.,20.,50.,100.]

"""
    load_input(file)

Loads the input file. It should have the following format *variable name*, *value*, *(description)*:

The variable names are 

*cs_file*,
*r_grid start*,
*r_grid step*,
*r_grid stop*,
*e0*,
*num_electrons*,
*jet_radius*,
*crossSectionSetID*,
*extraInfo*

and can be in any order.
"""
function load_input(file)
    #Load input csv
    inputs = readdlm(file,',')
    #Place them in a dictionary
    inputs = Dict(zip(inputs[:,1],inputs[:,2]))
    # Assign each value to a global variable
    global cs_file = inputs["cs_file"]
    global r_gridStart = inputs["r_grid start"]
    global r_gridStep = inputs["r_grid step"] 
    global r_gridStop = inputs["r_grid stop"] 
    global r_grid = [x for x in inputs["r_grid start"]:inputs["r_grid step"]:inputs["r_grid stop"]] #depth grid
    global e0 = inputs["e0"] #Escape threshold
    global num_electrons = inputs["num_electrons"] # Number of random walks to perform
    global jet_radius = inputs["jet_radius"] # 10 micron radius liquid microjet
    global crossSectionSetID = inputs["crossSectionSetID"] #ID of cross-section set ID
    global extraInfo = inputs["extraInfo"]
    return nothing
end


"""
    gaussian(x,amp,centre,FWHM)

A gauss amp function with FWHM"""
function gaussian(x,amp,centre,FWHM)
    wid = FWHM * 0.60056
    y = amp * exp( -((x - centre)/wid)^2)
    return y
end

"""
    concProfile(x,amp,centre,FWHM,offset)

Defines a concentration profile in terms of a Gaussian with a vertical offset.
"""
function concProfile(x,amp,centre,FWHM,offset)
    return gaussian(x,amp,centre,FWHM) + offset
end

"""
    standProfile(x, amp, offset)

Defines a standard concentratin profile with its gaussian centred
at 0.1 nm and with a FWHM of 0.4 nm.
"""
standProfile(x, amp, offset) = concProfile(x, amp, 0.1, 0.4, offset)

"""
    normaliseProfile(r_grid,concProfile)

Normalises the concentration profile using rectangular integration.

```r_grid``` must be equally spaced by ```r_gridStep```"""
function normaliseProfile(r_gridStep,concProfile)
    normFactor = r_gridStep * sum(concProfile)
    return 1/normFactor .* concProfile
end
"""
    applyConcProfile(basis_mat,concProfile)

Multiplies each depth slice by the values in concProfile.
"""
function applyConcProfile(basis_mat,concProfile)
    for i in axes(basis_mat,2)
        basis_mat[:,i,:] = basis_mat[:,i,:] .* concProfile[i]
    end
    return basis_mat
end

"""
    integrateAxis(scaledBasis; dims = 2)

Integrates a 3D array over one of its axes and return an array with one fewer dimension.

`dims = 2`: integrate over depth

`dims = 1`: integrate over initial energy

`dims = 3`: integrate over 
"""
function integrateAxis(scaledBasis; dims = 2)
    basis_integrated = sum(scaledBasis, dims=dims)
    basis_integrated = dropdims(basis_integrated, dims=dims)
end

#open HDF5 file
function main()
    load_input("PostUpgrade\\hdf5 test\\input_scattering_test.csv")

    basis_mat = h5open(h5path,"r") do fid
        read(fid["data"]["bases"])
    end

    #Create the file if it doesn't exist
    h5open(fileName, "cw") do fileID
    end

    #Overwrite the file and open in write mode
    h5open(fileName, "w") do fileID
        #Copy metadata, cross-sections and ELPS
        h5open(h5path, "cw") do fileIDallData
            #copy cross sections, ELPs and metadata to new H5 file
            HDF5.copy_object(fileIDallData["/crossSections&ELPs"], fileID, "crossSections&ELPs")
            HDF5.copy_object(fileIDallData["/metadata"], fileID, "metadata")
        end
        
        #Create data group an add scaled and integrated bases
        dataGroup = create_group(fileID,"data")
        for amp in amplitudes
            profile = [standProfile(x,amp,1) for x in r_grid]
            profile = normaliseProfile(r_gridStep, profile)
            scaledBasis = applyConcProfile(basis_mat,profile)
            basis_integrated = integrateAxis(scaledBasis)

            dataGroup["basisSet_$amp"] = hcat(basis_integrated...)
        end
        attributes(fileID["metadata"])["gaussianAmplitudes"] = amplitudes
    end
    return nothing
end

main()