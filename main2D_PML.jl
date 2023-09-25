# 2D FDTD simulation with PML absorbing boundary conditions
# source of the method : Perfectly Matched Layer for the Absorption of Electromagnetic Waves by Jean-Pierre Bérenger, 1993

using PyPlot

const c    = 299792458.0    # speed of light           [m/s] = 1 / sqrt(mu0*eps0)
const lam  = 12.0           # wavelength               [m]
const mu0  = 1.256e-6       # mangetic permeability    [N/A^2]
const E0   = 30.0           # electric field amplitude [V/m]

const Z0   = mu0 * c        # impedance of free space  [V/A]
const eps0 = 1 / mu0 / c^2  # electric permittivity    [F/m]

const Q_Hz0  = 1.0 * E0 / lam # Hz source amplitude    [V/m^2]
const tau    = 1 * lam / c  # Hz source duration     [s]

const npml   = 60          # number of PML layers

struct Field # 2D field
    Ex::Array{Float64, 2}
    Ey::Array{Float64, 2}
    Hz::Array{Float64, 2}
    Hzx::Array{Float64, 2}
    Hzy::Array{Float64, 2}
end

struct Material # 2D material
    epr::Array{Float64, 2}
    mur::Array{Float64, 2}
end

struct Grid # 2D grid
    nx::Int
    ny::Int
    dx::Float64
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
    epr = ones(grid.nx, grid.ny)
    mur = ones(grid.nx, grid.ny)

    # Here define material and geometry
    epr[:, 250:400] .= 3.0
    epr[:, 400:-1] .= 5.0

    return Material(epr, mur)
end

function init_grid_and_pml()
    # Grid
    nx = 500 # number of grid points
    ny = 500 # number of grid points
    dx = lam / 20.0 # grid spacing
    dt = dx / (2.1 * c) # time step
    t0 = 6.0 * tau # time delay
    w0 = 2 * π * c / lam # angular frequency
    
    # Grid coordinates
    x = 0:dx:(nx-1)*dx 
    y = 0:dx:(ny-1)*dx

    # Grid object
    grid = Grid(nx, ny, dx, dt, t0, w0, x, y)

    # PML
    pmlfac = 1/mu0 # scaling factor
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
    sidy = div(grid.ny, 3)

    # Source on Hz
    field.Hzx[sidx,sidy] += grid.dt / mu0 * gaussian_source(t - grid.t0, Q_Hz0)
    field.Hzy[sidx-40:sidx+40,sidy] .+= grid.dt / mu0 * gaussian_source(t - grid.t0, Q_Hz0)

    # Solve Hz field
    for i in 2:grid.nx, j in 2:grid.ny
        field.Hzx[i, j] = field.Hzx[i, j]*(1 - grid.dt * qy[i, j]) + (grid.dt / mu0) * ((field.Ex[i, j] - field.Ex[i, j-1]) / grid.dx)
        field.Hzy[i, j] = field.Hzy[i, j]*(1 - grid.dt * qx[i, j]) - (grid.dt / mu0) * ((field.Ey[i, j] - field.Ey[i-1, j]) / grid.dx)
    end

    # Combine Hzx and Hzy to get the gradient of Hz
    @. field.Hz = field.Hzx + field.Hzy
    
    # Solve Ex and Ey fields
    for i in 1:grid.nx-1, j in 1:grid.ny-1
        field.Ex[i, j] = field.Ex[i, j]*(1 - grid.dt * qy[i, j]) + (grid.dt / (eps0 * mat.epr[i, j])) * (field.Hz[i, j+1] - field.Hz[i, j]) / grid.dx
        field.Ey[i, j] = field.Ey[i, j]*(1 - grid.dt * qx[i, j]) - (grid.dt / (eps0 * mat.epr[i, j])) * (field.Hz[i+1, j] - field.Hz[i, j]) / grid.dx
    end
end

# Source functions
function gaussian_source(Δt::Real, Q0::Real)
    return Q0 * exp(- Δt^2 / tau^2)
end

function plot_loop_field(field::Field, grid::Grid, it::Int, mod::Int)
    if it % mod == 0
        figure(1)
        imshow(field.Ex[npml+1:end-npml, npml+1:end-npml])
        clim(-0.1, 0.1)
        if it == mod
            colorbar()
        end
        pause(0.000001)
    end
end

function main()
    
    # Initialize grid and PML
    grid, qx, qy = init_grid_and_pml()
    field = init_field(grid)
    mat = init_material(grid)
    
    nt = 2000 # number of time steps
    t = 0.0 # time
    
    # Recorded field
    trace_rec = zeros(nt)

    # time vector for comparison
    tvec = 0:grid.dt:(nt-1)*grid.dt

    # Main loop
    for it in 1:nt

        # Update & plot field
        update_field!(field, mat, grid, qx, qy, t)
        plot_loop_field(field, grid, it, 40)

        # Record trace
        trace_rec[it] = field.Ex[div(grid.nx, 2)-5, div(grid.ny, 3)-5]

        # Increment time
        t += grid.dt
    end

    # Plot recorded trace
    figure(2)
    plot(tvec, trace_rec)
    title("Recorded trace")
    xlabel("Time [s]")
    ylabel("Electric field [V/m]")
    savefig("trace_rec.png")
    show()
end

main()