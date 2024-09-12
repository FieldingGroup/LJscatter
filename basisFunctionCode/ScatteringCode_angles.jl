#!/usr/bin/env julia

import Pkg
Pkg.activate(".")

using Base.Iterators
using DelimitedFiles
using Distributions
using Interpolations
using LinearAlgebra
using LsqFit
using Rotations
using StaticArrays
using StatsBase
using HDF5
using Roots


# const cs_file = "cs_improved_1b.csv"
# const r_grid = [x for x in 0.1:0.1:50]
# const depth_step = 0.1 #r_grid step size
# const e0 = 0.2
# const num_electrons = 10000 # Number of random walks to perform
# const jet_radius = 10000 # 10 micron radius liquid microjet
const SolidAngle = false


const loss_channels_normal = Dict("vpp_L" => (0.092, 0.04),
                                  "vp_L"  => (0.061, 0.03),
                                  "vpp_T" => (0.025, 0.024),    
                                  "vp_T"  => (0.01, 0.001), 
                                  "v2"      => (0.204, 0.016),
                                  "v1-3"    => (0.417, 0.05),
                                  "v3"      => (0.460, 0.005),
                                  "v1-3+vL" => (0.500, 0.04),                
                                  "2v1-3"   => (0.835, 0.075),       
                                  "electronic" => (5.0, 0.0), # this used to be 5,5 it is now (5,0), i.e. electrons are discarded
                                  "elastic" => (0.0, 0.0))

# -------- EXTRA DATA ------------
reflected_electrons = 0 #number of electrons reflected at surface. This excludes doubly reflected electrons              
doubly_reflected_electrons = 0 #number of electrons that should have been doubly reflected. These electrons are discarded
#approach_angles = [] #angles at which electrons approach the surface #INCOMPATIBLE WITH TRHEADS
#escape_angles = [] #angles at which electrons escape the surface #INCOMPATIBLE WITH TRHEADS
# -------------------------                   


"""
    load_input(file)

Loads the input file. It should have the following format *variable name*, *value*, *(description)*:

The variable names are 
 *cs_file*,
 *r_grid*,
 *e0*,
 *num_electrons*,
 *jet_radius*
and can be in any order.
"""
function load_input(file)
    #Load input csv
    inputs = readdlm(file,',')
    #Place them in a dictionary
    inputs = Dict(zip(inputs[:,1],inputs[:,2]))
    # Assign each value to a global variable
    global cs_file = inputs["cs_file"]
    global r_grid = [x for x in inputs["r_grid start"]:inputs["r_grid step"]:inputs["r_grid stop"]] #depth grid
    global depth_step = inputs["r_grid step"] #r_grid step size
    global e0 = inputs["e0"] #Escape threshold
    global num_electrons = inputs["num_electrons"] # Number of random walks to perform
    global jet_radius = inputs["jet_radius"] # 10 micron radius liquid microjet
    return nothing
end

"""
    load_cs_csv(file)

The (in)elastic scattering cross sections per-channel as a function of electron are required as input as a
CSV file. 
The cross sections are already inter/extrapolated.
"""
function load_cs_csv(file)
    #Load cross section data
    data, headers = readdlm(file, ',',header=true)
    #Remove first column aka eKEs
    table = data[:,2:end]
    #Make list of channels in the same order as in table
    channels = [loss_channels_normal[channel] for channel in headers[2:end]]
    return table, channels

end


"""
    rand_normal_positive(μ, fwhm)
Generate a number from a normal distribution but force it to be
positive.

This function is used exclusively for the randomised eKE loss per
inelastic collision. The widths of the librational resonance are
given as FWHMs in the paper but we need the standard deviation.
"""
function rand_normal_positive(μ, fwhm)
    while true
        x = rand(Normal(μ, fwhm / 2.35482))
        if x >= 0
            return x
        end
    end
end


"""
    eval_cs_lookup(eKE, table)
Returns the cross sections at a specific eKE"""
function eval_cs_lookup(eKE, table)
    idx = round(Int32, eKE * 100) + 1 #+1 as julia's indexing is one-based 
    return table[idx, :]
end

"""
    lose_eKE(eKE, table::Matrix{Float64}, channels)
Samples the cross-sections at a certain eKE and returns
the energy loss and the step length sampled from and exponential
distribution of the mean free path,
"""
function lose_eKE(eKE, table::Matrix{Float64}, channels)
    w = eval_cs_lookup(eKE, table)

    # Choice of channel is weighted by the relative cross sections.
    channel = sample(channels, weights(w))

    # Compute the (positive) energy loss.
    loss = rand_normal_positive(channel...) #three dots just splits the tuple

    # The total cross section.
    if SolidAngle
        σ = sum(w) * 1e-2 * (1 / (4π))
    else
        σ = sum(w) * 1e-2
    end
    # Number density of scatterers.
    ρ = (1e-21 / 18.015) * 6.022141e23

    # The mean free path.
    Λ = 1 / (ρ * σ)

    # Additionally return the random step length from an exponential
    # distribution.
    return loss, rand(Exponential(Λ))
end

"""
    random_walk_jet(eKE, pos, table, channels)
Recursively runs a single random walk step, by updating the energy and position
of the electron. The energy loss and step length are sampled using `lose_eKE` and 
all scattering events are assumed to be isotropic.
When electrons reach the surface they can escape if they have enough energy or
be reflected. Multiple reflections in
"""
function random_walk_jet(eKE, pos, table, channels; e0=e0)
    if eKE < 0.0
        return NaN
    end

    loss, dr = lose_eKE(eKE, table, channels)

    # Use a random orthonormal rotation matrix to define the isotropic
    # step.
    rotation_matrix = rand(RotMatrix{3})
    dr_vector = SVector(dr, 0.0, 0.0)

    # Rotate the movement vector and update the electron position.
    dpos = rotation_matrix * dr_vector
    new_pos = pos + dpos

    # The jet is considered to be infinite in the z-axis, with slices
    # along the xy-plane forming the circular cross section of the
    # jet. Therefore, we only consider the x- and y-coordinates
    # (pos[1] and pos[2]) when checking whether the electron has left
    # the jet.
    r = norm(new_pos[1:2])
    
    #If the radial coordinate r is greater than the jet radius,
    #the electron may escape or be reflected.
    #An electron is reflected if the its energy perpendicular to the surface
    #is higher than the potential energy barrier given b the electron affinity.
    if (r >= jet_radius)
        #Determine the coordinates of the point at which the electron escapes and
        # the component of the eKE (relative to the bottom of the conduction band)
        #perpendicular to the surface
        norm_projection, escape_point, approach_angle = escape_filter(eKE+e0, new_pos, pos; R = jet_radius)
        if norm_projection >= e0
            return approach_angle
        else
            new_pos = reflect_electron(new_pos, escape_point)
            #If the angle of approach to the surface is very shallow,
            #and the electron does not have enough energy to escape,
            #multiple reflection must occur within one step.
            #This is very rare and hard to deal with,
            #so (for now?) this electrons will be discarded, but recorded.
            if  norm(new_pos[1:2])>=jet_radius
                println("P1,P2,Q: $(new_pos), $(pos), $(escape_point)")
                return -1
            end
        end
    end
    # Determine the eKE loss for the collision and take an
    # isotropic random step.
    # And recursively run this function again, accumulating the
    # eKE loss and the change in position.
    return random_walk_jet(eKE - loss, new_pos, table, channels)
end

"""
    escape_filter(eKE, new_pos, old_pos; r = jet_radius)
This function determines the intersection between the vector going through two points
and the circumference of the jet. It determines the velocity component of the vector
perpendicular to the surface and finds the corresponding kinetic energy.
Input energy must be wrt to the bottom of the conduction band.
The approach angle can also be determined.
"""
function escape_filter(eKE, new_pos, old_pos; R = jet_radius)
    P1 = old_pos
    P2 = new_pos
    
    #Determine parameters of line going through P1 and P2 and
    #find the intersections with the surface of the circle.

    #If points are in a vertical line, y = mx+c cannot be used.
    if P1[1] == P2[1]
        x1 = x2 = P1[1]
        y1 = sqrt(R^2 - x1^2) #Pythagoras
        y2 = -y1
    else
        #y = mx+c
        m = (P2[2] - P1[2])/(P2[1] - P1[1]) #Δy/Δx
        c = P2[2] - m*P2[1] #y-mx

        #Intersection points
        #Solving system of equations y^2 + x^2 = r; y = mx+c
        # which gives (1+m^2)x^2 + 2mcx + c^2-R^2 = 0
        A = 1+m^2
        B = 2*m*c
        C = c^2-R^2

        delta = B^2-(4*A*C)
        # if the line is close to vertical,
        #delta can be negative due to floating point error
        #thus the module Roots should be used. (since it is slower, it is not always used)
        #This is also very rare.
        if delta <0
            #Intersection points between line and right semi-circle
            println("delta<0: P1,P2 :$P1 ,$P2 ")
            f(y) =sqrt(R^2-y^2) - ((y-c)/m)
            ys = find_zeros(f,(-R,R))
            y1 = ys[1]
            x1 = (y1-c)/m
            try
                y2=ys[2]
                x2 = (y2-c)/m
            catch
            end
        else
            #Solve quadratic
            x1 = (-B + sqrt(delta))/(2*A)
            y1 = m*x1+c
            x2 = (-B - sqrt(delta))/(2*A)
            y2 = m*x2+c
        end
    end
    #There are two intersection points per line.
    #Only the point Q between P1 and P2 is relevant
    if P1[1] <= x1 <= P2[1] || P2[1] <= x1 <= P1[1]
        Q = x1,y1
    else
        Q = x2,y2
    end
    #Q in polar coordinates is described by (r,θ). 
    #The slope of the normal to the circle at r,θ is given by tan(θ),
    #which is equal to y/x (cartesian coordinates).
    if Q[1] == 0 #the normal vector is vertical
        normal_vec = [0,1,0]
    else
        normal_slope = Q[2]/Q[1]
        normal_vec = [1,normal_slope,0]
        normal_vec = normal_vec / norm(normal_vec)
    end
    #Projecting P2-P1 onto the normal
    P_vec = (P2-P1)
    P_vec = P_vec / norm(P_vec)
    # Energy projection = E (P.n)^2
    E_proj = eKE * dot(P_vec, normal_vec)^2

    #Determine approach angle, cos ϕ = (P.n)
    ϕ = acos(dot(P_vec, normal_vec))
    return E_proj, Q, ϕ
end

"""
    reflect_electron(new_pos, Q)
Reflects electron on the surface of the jet.
It determines the equation of the tangent at the exit point and of the perpendicular line going through
the unreflected point (P2), it finds their intersection A and reflects the electron through it.
This is to determine the *x* and *y* coordinates of the reflection.
The *z* coordinate is the same for both the transmitted and reflected electron

Does not include double reflection
"""
function reflect_electron(new_pos, Q)
    P2 = new_pos
    #Finding the point of reflection A(x0,y0)
    #Checking for special cases (ie when slopes would be Inf)
    if Q[1] == 0 # normal is vertical and A's coordinates are (P2x,±r)
        x0, y0 = P2[1], Q[2]
    elseif Q[2] == 0 # tangent is vertical and A's coordinates are (±r,P2y)
        x0, y0 = Q[1], P2[2]
    else #normal and tangent are not horizontal/vertical
        #Parameters of the normal going through P2
        normal_slope = Q[2]/Q[1]
        normal_intercept = P2[2] - normal_slope * P2[1]
        #Parameters of tangent
        tangent_slope = -1/normal_slope
        tangent_intercept = Q[2] - tangent_slope * Q[1]
        #Find intersection point A
        x0 = (tangent_intercept - normal_intercept) / (normal_slope - tangent_slope)
        y0 = normal_slope * x0 + normal_intercept
    end
    #The reflected position is obtained by subtracting the difference between P2 and
    #A from P2 twice. z stays the same.
    refl_pos = P2 - 2*[P2[1]-x0, P2[2]-y0, 0]
    if norm(refl_pos[1:2]) >= jet_radius
        global doubly_reflected_electrons
        doubly_reflected_electrons += 1
        # println("double refl needed")
        # throw(ErrorException("The electron ended up outside the jet at $(refl_pos).\nDouble reflection is not supported yet.\nP1: $P1 \nP2: $P2 \nQ: $Q"))
    else
        global reflected_electrons
        reflected_electrons +=1
    end
    return refl_pos
end


"""
    run_jet_walk(starting_eKEs, depth, table, channels; e0=e0)
This function runs `random_walk_jet` for all starting eKEs at a certain depth.
"""
function run_jet_walk(starting_eKEs, depth, table, channels; e0=e0)
    pos = SVector(jet_radius - depth, 0.0, 0.0)
    out = similar(starting_eKEs)
    for i in eachindex(out)
        out[i] = random_walk_jet(starting_eKEs[i], pos, table, channels; e0=e0)
    end
    return out
end


function integrate_over_r(starting_eKEs, table, channels;
                          e0=e0, r_grid=1:0.25:1001)
    outputs = similar(r_grid,Vector{Float64})
    Threads.@threads for i in eachindex(r_grid)
        r = r_grid[i]
        o = run_jet_walk(starting_eKEs, r, table, channels; e0=e0)::Vector{Float64}

        # Remove electrons with eKE <= 0 eV
        o = filter((x) -> x != NaN, o)

        # Bin the data onto the 0 : 0.01 : 5 eV grid and get the counts per bin.
        d = bin_data(o)
        y = d[:, 2]

        outputs[i] =  y .* (2 * pi * (jet_radius - r) * depth_step)

        
    end

    return outputs
end

function apply_pdf(basis, dist; r_grid=1:0.25:1001)
    reduced = similar(basis)
    # At each depth, we need to apply the PDF to the basis functions.
    #
    # In the simplest case, the PDF comes from a distribution provided by
    # Distributions.jl (e.g. Normal() or Exponential()) and these are evaluated
    # simply using pdf(dist, x) for each x in r_grid.
    #
    # For a custom profile, you would need an array with the probability at each
    # point in r_grid. Presumably, this array would need to be normalised to sum
    # to 1.
    p = [ pdf(dist, x) for x in r_grid ]
    for (i, b) in enumerate(basis)
        basis_w_pdf = b .* p
        reduced[i] = [vec(reduce(+, hcat(basis_w_pdf...)', dims=1))]
    end

    reduced = map(x -> x[1], reduced)
    return reduced
end



# Take a set of data points and transfer it onto a binned grid. Typically,
# we employ a grid spacing of 0.01 eV between 0.01 and 5.00 eV. The default
# grid starts one point earlier due to a quirk in Histogram.
function bin_data(data; bins=0.00:0.01:2π)
    f = fit(Histogram, data, bins,closed=:right)

    # The last bin is removed due to the previously mentioned quirk that makes
    # this list one item too long.
    x = collect(f.edges[1])[1:end-1]

    # The weights are the number of counts in each bin.
    y = f.weights
    return hcat(x, y)
end

function main()
    load_input("input_scattering.csv")
    table,channels = load_cs_csv(cs_file)

    basis_grid = [ ones(num_electrons) * i for i in 0.01:0.01:5.0 ]
    basis_10 = map(x -> integrate_over_r(x, table, channels; e0=e0, r_grid=r_grid), basis_grid)
    
    basis_uni = apply_pdf(basis_10, Uniform(first(r_grid),last(r_grid)), r_grid=r_grid)

    writedlm("bases/basis_angles.txt", hcat(basis_uni...))

    basis_mat = cat(map(x ->hcat(x...) ,basis_10)...,dims=3)
    fid = h5open("bases/allAngles.h5","cw") #this creates the file if it doesn't exist
    close(fid)
    fid = h5open("bases/allAngles.h5","w") #this overwrites the file
    fid["bases$(e0)"] = basis_mat
    fid["reflected_electrons"] = reflected_electrons
    fid["doubly_reflected_electrons"] = doubly_reflected_electrons
    close(fid)

end
main()
