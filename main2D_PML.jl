# 2D FDTD simulation with PML absorbing boundary conditions
# source of the method : Perfectly Matched Layer for the Absorption of Electromagnetic Waves by Jean-Pierre Bérenger, 1993

using PyPlot
using PyCall
tqdm = pyimport("tqdm")

const c  = 299792458.0       # speed of light           [m/s] = 1 / sqrt(mu0*eps0)
const λ  = 5.0e-2            # wavelength               [m]
const μ0 = 1.256e-6          # mangetic permeability    [N/A^2]
const E0 = 30.0              # electric field amplitude [V/m]
const Z0 = μ0 * c            # impedance of free space  [V/A]
const ϵ0 = 1 / μ0 / c^2      # electric permittivity    [F/m]
const Q_Hz0 = 1.0 * E0 / λ   # Hz source amplitude    [V/m^2]
const τ     = 1.0 * λ / c    # Hz source duration     [s]
const npml = 100             # number of PML layers

struct Field # 2D field
    Ex::Array{Float64, 2}
    Ey::Array{Float64, 2}
    Hz::Array{Float64, 2}
    Hzx::Array{Float64, 2}
    Hzy::Array{Float64, 2}
end

struct Material # 2D material
    ϵr::Array{Float64, 2}
    μr::Array{Float64, 2}
    σ::Array{Float64, 2}
end

struct Grid # 2D grid
    nx::Int
    ny::Int
    dx::Float64
    dy::Float64
    dt::Float64
    t0::Float64
    w0::Float64
    x::Array{Float64, 1}
    y::Array{Float64, 1}
end

function init_field(grid::Grid)
    Ex = zeros(grid.nx, grid.ny)
    Ey = zeros(grid.nx, grid.ny)
    Hz = zeros(grid.nx, grid.ny)
    Hzx = zeros(grid.nx, grid.ny)
    Hzy = zeros(grid.nx, grid.ny)
    return Field(Ex, Ey, Hz, Hzx, Hzy)
end

function init_material(grid::Grid)
    ϵr = ones(grid.nx, grid.ny)
    μr = ones(grid.nx, grid.ny)
    σ  = zeros(grid.nx, grid.ny)

    # Add some ice
    ϵr[:, 250:end] .= 3.2
    σ[:, 250:end] .= 5.0e-8

    # Add some rocks
    ϵr[:, 600:end] .= 5.0
    σ[:, 600:end] .= 1.0e-2


    # Here define reflective surface
    # σ[ceil(Int, grid.nx/2 - 100):ceil(Int, grid.nx/2 + 100), 150] .= 3.5e7

    return Material(ϵr, μr, σ)
end

function init_grid_and_pml()
    # Grid
    xsize = 5.0 # [m]
    ysize = 5.0 # [m]
    dx = 0.005 # grid spacing
    dy = 0.005 # grid spacing
    nx = ceil(Int, xsize / dx) + 1 # number of grid points
    ny = ceil(Int, ysize / dy) + 1 # number of grid points
    dt = min(dx, dy) / (2.1 * c) # time step
    t0 = 6.0 * τ # time delay
    w0 = 2 * π * c / λ # angular frequency
    
    # Grid coordinates
    x = 0:dx:(nx-1)*dx 
    y = 0:dx:(ny-1)*dx

    # Grid object
    grid = Grid(nx, ny, dx, dy, dt, t0, w0, x, y)

    # PML
    pmlfac = 2.0e2/μ0 # scaling factor
    pmlexp = 1.0   # scaling exponent

    # Initialize damping arrays
    qx = zeros(Float64, nx+1, ny+1)
    qy = zeros(Float64, nx+1, ny+1)

    # Calculate damping arrays
    for a in 1:npml
        qx[a,:] .= pmlfac * (npml - a)^pmlexp
        qx[nx-a+1,:] .= pmlfac * (npml - a)^pmlexp
        qy[:,a] .= pmlfac * (npml - a)^pmlexp
        qy[:,ny-a+1] .= pmlfac * (npml - a)^pmlexp
    end

    return grid, qx, qy, npml
end

function update_field!(field::Field, mat::Material, grid::Grid, qx::Array{Float64, 2}, qy::Array{Float64, 2}, t::Real)

    # Source position
    sidx = div(grid.nx, 2)
    sidy = 155

    # Source on Hz
    field.Hzx[sidx-25:sidx+25,sidy] .+= grid.dt / μ0 * gaussian_source(t - grid.t0, Q_Hz0)
    field.Hzy[sidx-25:sidx+25,sidy] .+= grid.dt / μ0 * gaussian_source(t - grid.t0, Q_Hz0)

    # Solve Hz field
    for i in 2:grid.nx, j in 2:grid.ny
        field.Hzx[i, j] = field.Hzx[i, j]*(1 - grid.dt * qy[i, j]) + (grid.dt / μ0) * ((field.Ex[i, j] - field.Ex[i, j-1]) / grid.dy)
        field.Hzy[i, j] = field.Hzy[i, j]*(1 - grid.dt * qx[i, j]) - (grid.dt / μ0) * ((field.Ey[i, j] - field.Ey[i-1, j]) / grid.dx)
    end

    # Combine Hzx and Hzy to get the gradient of Hz
    @. field.Hz = field.Hzx + field.Hzy
    
    # Solve Ex and Ey fields
    for i in 1:grid.nx-1, j in 1:grid.ny-1
        # Update Ex with conductivity
        field.Ex[i, j] = (field.Ex[i, j] * (1 - grid.dt * qy[i, j]) 
            + (grid.dt / (ϵ0 * mat.ϵr[i, j])) * (field.Hz[i, j+1] - field.Hz[i, j]) / grid.dy) / (1 + grid.dt * mat.σ[i, j] / (2 * ϵ0 * mat.ϵr[i, j]))
        # Update Ey with conductivity
        field.Ey[i, j] = (field.Ey[i, j] * (1 - grid.dt * qx[i, j]) 
            - (grid.dt / (ϵ0 * mat.ϵr[i, j])) * (field.Hz[i+1, j] - field.Hz[i, j]) / grid.dx) / (1 + grid.dt * mat.σ[i, j] / (2 * ϵ0 * mat.ϵr[i, j]))
    end
end

function gaussian_source(Δt::Real, Q0::Real)
    return Q0 * exp(- Δt^2 / τ^2)
end

function plot_loop_field(field::Field, grid::Grid, it::Int, mod::Int)
    if it % mod == 0
        figure(1, figsize=(5, 12))
        # imshow(field.Ex[npml+1:end-npml, npml+1:end-npml])
        imshow(field.Ex')
        title("Ex field")
        plot([ceil(Int, grid.nx/2-100), ceil(Int, grid.nx/2+100)], [150, 150], "r")
        clim(-0.5, 0.5)

        if it == 1
            colorbar()
        end
        pause(0.000001)
    end
end

function plot_trace(trace::Array{Float64, 1}, tvec::Array{Float64, 1})
    figure(2)
    plot(tvec, trace)
    title("Trace")
    xlabel("Time [s]")
    ylabel("Amplitude [V/m]")
    show()
end

function main()
    
    # Initialize grid and PML
    grid, qx, qy = init_grid_and_pml()
    field = init_field(grid)
    mat = init_material(grid)
    
    nt = 4000 # number of time steps
    t = 0.0 # time
    
    # Recorded field
    trace_rec = zeros(nt)

    # time vector for comparison
    tvec = 0:grid.dt:(nt-1)*grid.dt

    println("Start simulation...")

    # Set up tqdm progress bar
    pbar = tqdm.tqdm(1:nt, desc="Simulating", ncols=100)

    # Main loop
    for it in pbar
        # Update & plot field
        update_field!(field, mat, grid, qx, qy, t)
        plot_loop_field(field, grid, it, 50)

        # Record trace
        trace_rec[it] = field.Ex[div(grid.nx, 2), 150]

        # Increment time
        t += grid.dt
    end

    # Plot trace
    plot_trace(trace_rec, tvec)
end


main()