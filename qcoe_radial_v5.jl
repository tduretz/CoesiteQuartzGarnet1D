# Coesite-Quartz transformation with garnet host + SiO2 matrix 
using CairoMakie
using SpecialFunctions
using Printf
using MAT
import LinearAlgebra:norm
kyr     = 1000*365.25*24*3600

function main()

    BC = :ConstantPressureRate
    # BC = :ConstantDilationRate

    # Spatial domain:
    rmin    = 0.0        # m
    rmax    = 0.005       # m
    ε̇bg     = 2e-15
    dPdtbg  = -2/(1e-9*(1e3*kyr))
    # material properties (see Table 1):
    ρ0i     = 2600       # kg.m-3
    ρ0f     = 2900       # kg.m-3
    ρ0g     = 3200       # kg.m-3
    βg      = 1.0/(170e9) # kg.m-3
    Pr      = 3.0e9      # Pa
    dPr     = 1e8        # Pa
    Pbg     = 4.5e9      # Pa
    K       = 80e9       # Pa
    βq      = 1.0/K        # Pa-1  
    τk      = 0*3.1558e9   # s
    ηi      = 1e19       # Pa.s
    ηg      = 1e23       # Pa.s
    ηq      = 1e19       # Pa.s
    G       = 40e9       # Pa
    ϕ       = 20.0*π/180.
    C       = 1e7
    # Numerical parameters
    nt      = 10000    # nb of steps
    Δt      = 20*0.5*2.5*1e10#τk/10      # t step
    ηe      = G*Δt       # ηe definition
    ncr     = 50         # number of cells 
    niter   = 3e5        # nb of iterations
    # Pseudo-transient parameters
    αV         = 1.0
    αP         = 0.25
    θ          = 25*(rmax-rmin)/ncr

    nout       = 5000
    rel_tol_V  = 5e-6
    rel_tol_P  = 5e-6
    # Define grid with variable spacing
    rv = LinRange(-rmin, rmax, ncr+1)
    rc = 0.5.*(rv[2:end] .+ rv[1:end-1])
    Δrc = rc[2] - rc[1]
    Δri = Δrc
    # Alocate arrays
    Rv, ∂V∂τ         = zeros(size(rv)), zeros(size(rv))
    Rp, ∂P∂τ         = zeros(size(rc)), zeros(size(rc))
    Vr, Vrc          = zeros(size(rv)),  zeros(size(rc))
    ∇V               = zeros(size(rc))
    P, P0            = zeros(size(rc)), zeros(size(rc))
    β                = zeros(size(rc))
    ρ, ρ0,  ρr       = zeros(size(rc)), zeros(size(rc)), zeros(size(rc))
    τrr, τrr0        = zeros(size(rc)), zeros(size(rc))
    τϕϕ, τϕϕ0, τII   = zeros(size(rc)), zeros(size(rc)), zeros(size(rc))
    σrr, σrri        = zeros(size(rc)), zeros(size(rc,1)-1)
    σϕϕ, σϕϕi, σΘΘ   = zeros(size(rc)), zeros(size(rc,1)-1), zeros(size(rc))
    ε̇rr, ε̇ϕϕ, εrr    = zeros(size(rc)), zeros(size(rc)), zeros(size(rc))
    Ėrr, Ėϕϕ, ĖII_ve = zeros(size(rc)), zeros(size(rc)), zeros(size(rc))
    phc              = zeros(size(rc))
    ηc, ηvec         = zeros(size(rc)), zeros(size(rc))
    X, X0            = zeros(size(rc)), zeros(size(rc))
    F, λ̇             = zeros(size(rc)), zeros(size(rc))
    # Initialisation
    Vr              .= ε̇bg.*(rv.-0*(rmax-rmin)/2)
    P               .= Pbg
    ρ               .= ρ0i .* exp.(β.*P)
    # phc[rc .> 0.003  .&& rc .< 0.0045] .= 1
    # phc[rc .> 0.0055 .&& rc .< 0.007 ] .= 1
    # # phc[rc .< 0.002 .|| rc .> 0.008] .= 2
    phc[rc .> 0.00  .&& rc .< 0.001] .= 1
    phc[rc .> 0.002 .&& rc .< 0.003 ] .= 1
    # phc[rc .< 0.002 .|| rc .> 0.008] .= 2
    β               .= βq
    β[phc.==1]      .= βg 
    ηc              .= ηq  
    ηc[phc.==1]     .= ηg
    ηc[phc.==0]     .= ηi
    t                = 0.0
    rdPdt            = 0.
    # for plots:
    p_inc_vec  = zeros(nt)
    p_mat_vec  = zeros(nt)
    div_vec    = zeros(nt)
    ispl_vec   = zeros(nt)
    t_vec      = zeros(nt)
    maxτII_vec = zeros(nt)
    X_inc_vec  = zeros(nt)
    # # Equation of state
    P_eq = LinRange(0.1e9,5e9,200)
    ρ_eq = zeros(size(P_eq))
    ρ_eq .= ρ0i .+ (ρ0f-ρ0i) .* ( 1.0 .- (0.5 .* erfc.((P_eq .- Pr) ./dPr) ))  
    # p = plot(P_eq/1e9, ρ_eq.*exp.(βq.*P_eq), marker = 2, 
    #         title = "ρ", 
    #         xlabel = "P [GPa]",
    #         ylabel = "ρ [kg.m-3]")
    # display(p)

    # # Equilibrium x
    # X_eq = 1.0 .- 0.5 .* erfc.((P_eq .- Pr) ./ dPr)
    # p = plot(P_eq, X_eq, marker = 2, 
    # title = "X", 
    # xlabel = "P [GPa]",
    # ylabel = "X [-]")
    # display(p)

    # # Kinetics
    # P    = 2e9
    # X_t0 = 0.0
    # X_t   = zeros(nt)
    # for it=1:nt
    #     X_t[it] = (X_t0 .* τk + (X_eq .* Δt) ./ (τk + Δt)
    # end
    # p = plot(1:(nt*Δt), X_t)
    # display(p)

    # Initialise density
    ρr          .= ρ0i .+ (ρ0f-ρ0i) .* ( 1.0 .- (0.5 .* erfc.((P .- Pr) ./dPr) )) 
    ρr[phc.==1] .= ρ0g
    ρr[phc.==2] .= ρ0i
    ρ           .= ρr .*exp.( β .*P)

    # TIME LOOP
    for it = 1:nt
        @show it
        # Initialise t step
        t    += Δt
        # Store old value
        τrr0 .= τrr
        τϕϕ0 .= τϕϕ
        ρ0   .= ρ
        X0   .= X
        P0   .= P
        # For error monitoring
        rr0, rp0, rdPdt0 = 0., 0., 0.
        # Pseudo transient iterations
        for iter = 1:niter
            # Kinematics:
            Vrc   .= 0.5 * (Vr[1:end-1] + Vr[2:end])
            ∇V    .= 1.0 ./rc .* diff(rv.*Vr, dims=1) ./ Δrc      
            ε̇rr   .= diff(Vr,dims=1) ./ Δrc - 1/3 .* ∇V
            ε̇ϕϕ   .= Vrc./rc  - 1/3 .* ∇V
            Ėϕϕ   .= ε̇ϕϕ .+ τϕϕ0./2.0./ηe
            Ėrr   .= ε̇rr .+ τrr0./2.0./ηe    
            # Rheology:
            ηvec  .= (1.0/ηe .+ 1.0./ηc).^(-1)
            τrr   .= 2.0*ηvec.*Ėrr
            τϕϕ   .= 2.0*ηvec.*Ėϕϕ
            σrr   .= -P .+ τrr
            σϕϕ   .= -P .+ τϕϕ
            # Plasticity 
            τII        .= sqrt.(0.5*(τrr.^2 .+ τϕϕ.^2 .+ (-τrr-τϕϕ).^2))
            F          .= τII .- C.*cosd(ϕ) .- P.*sin(ϕ)
            ĖII_ve     .= sqrt.(0.5*(Ėrr.^2 .+ Ėϕϕ.^2 .+ (-Ėrr-Ėϕϕ).^2))
            sumF = sum(F.>0)

            λ̇          .= 0.
            λ̇[F.>0]    .= F[F.>0]./ηvec[F.>0]
            
            ηvec[F.>0] .= (τII[F.>0] .- λ̇[F.>0].*ηvec[F.>0])./2.0./ĖII_ve[F.>0]
            τrr   .= 2.0*ηvec.*Ėrr
            τϕϕ   .= 2.0*ηvec.*Ėϕϕ
            σrr   .= -P .+ τrr
            σϕϕ   .= -P .+ τϕϕ
            σΘΘ   .= -P .+ (-τrr.-τϕϕ)
            τII   .= sqrt.(0.5*(τrr.^2 .+ τϕϕ.^2 .+ (.-τrr.-τϕϕ).^2))
            F     .= τII .- C.*cosd(ϕ) .- P.*sin(ϕ)
            # end
            # Non-linearity
            # X     .= (X0 .* τk + (1.0 .- 0.5 .* erfc.((P .- Pr) ./ dPr)) .* Δt) ./ (τk + Δt)
            # X[phc.==1] .= 0 
            # ρr    .= X.*ρ0f .+ (1.0 .-X).*ρ0i
            ρr    .= ρ0i .+ (ρ0f-ρ0i) .* ( 1.0 .- (0.5 .* erfc.((P .- Pr) ./dPr) )) 
            ρr[phc.==1] .= ρ0g
            ρr[phc.==2] .= ρ0i
            ρ     .= ρr .*exp.( β .*P)
            # Interpolation
            σrri  .= 0.5 .* (σrr[1:end-1] .+ σrr[2:end])
            σϕϕi  .= 0.5 .* (σϕϕ[1:end-1] .+ σϕϕ[2:end])
            # Residuals
            Rv[2:end-1] .=  diff(σrr,dims=1)./Δri .+ 1.0 ./ rv[2:end-1] .* (σrri.-σϕϕi)
            Rp          .= - ( (log.(ρ) .- log.(ρ0))./Δt .+ ∇V )
            
            if BC==:ConstantPressureRate
                dPdtE = ( (P[end]-P0[end])/Δt )  
                RvBC = dPdtE - dPdtbg
            elseif BC==:ConstantDilationRate
                RvBC = 0.
            end
            
            # Checks
            if iter==1 || mod(iter,nout) == 0
                # Errors
                rr     = norm(Rv)/length(Rv)
                rp     = norm(Rp)/length(Rp)
                if isnan(rr) error("NaN au fromage!") end
                if iter==1 rr0 = rr; rp0 = rp; rdPdt0 = rdPdt; end
                @show rdPdt
                @printf("Iteration %05d --- max F = %2.2e\n", iter, maximum(F))
                @printf("||dPdt|| = %2.10e %2.10e\n", rdPdt, rdPdt/rdPdt0)
                @printf("||fr||   = %2.10e %2.10e\n", rr, rr/rr0)
                @printf("||fp||   = %2.10e %2.10e\n", rp, rp/rp0)
                if ((rr/rr0<rel_tol_V || rr<rel_tol_V) && (rp/rp0<rel_tol_P || rp<rel_tol_P) && (rdPdt/rdPdt0<rel_tol_P || rdPdt<rel_tol_P)) break; end
            end
            # PT steps
            ΔtV    = αV*minimum(Δrc)^2/maximum(ηvec)/2.1
            ΔtP    = αP*minimum(ηvec)/ncr/2.1
            # Rate updates
            ∂V∂τ   .= Rv .+ (1.0-θ) .* ∂V∂τ 
            rdPdt   = RvBC .+ (1.0-1.0) * rdPdt 
            ∂P∂τ   .= Rp
            # Solution updates
            Vr     .+= ΔtV.*∂V∂τ
            P      .+= ΔtP.*∂P∂τ  
            # P .= .-(σrr .+ σϕϕ .+ σΘΘ)./3.0
            Vr[end] += Δt*rdPdt*1e-3/1e21

        end
        @show Vr[end]
        @show ( (P[end]-P0[end])/Δt ) * (1e-9*(1e3*kyr))
        # Compute deviators
        τII .= sqrt.(0.5 .* ( τrr.^2 .+ (-τϕϕ.-τrr).^2 .+ τϕϕ.^2))
        εrr .+= ε̇rr.*Δt 
        # Store solutions
        ispl_vec[it]   = sum(F.>0.) > 0.
        p_inc_vec[it]  = P[1]
        p_mat_vec[it]  = P[end]
        div_vec[it]    = (Vr[end]-Vr[1])/(rmax-rmin)
        t_vec[it]      = t
        maxτII_vec[it] = maximum(τII)
        X_inc_vec[it]  = X[1]
        # Visualisation
        if mod(it, 10)==0 || it==1
            f = Figure(resolution = (1200, 1000), fontsize=22)
            # -------------- #
            ax1 = Axis(f[1, 1], title = "Density - pressure", xlabel = "P [GPa]", ylabel = "ρ [kg.m-3]")
            lines!(ax1, P_eq/1e9, ρ_eq.*exp.(βq.*P_eq))
            scatter!(ax1, P[:]./1e9, ρ[:])
            xlims!(ax1, 0, 5)
            # -------------- #
            ax6 = Axis(f[1, 2], title = "Yield", ylabel = "τII [MPa]", xlabel = "P [GPa]")
            lines!(ax6, P_eq/1e9, (C*cos(ϕ) .+ P_eq.*sin(ϕ))./1e6) 
            scatter!(ax6, P[:]./1e9, τII./1e6)
            xlims!(ax6, 0, 5)
            # -------------- #
            ax2 = Axis(f[2, 1], title = "ρ", xlabel = "r [m]", ylabel = "ρ [kg.m-3]")
            lines!(ax2, rc, ρ, label="ρ [kg.m-3]") 
            # ax2 = Axis(f[2, 1], title = "εrr", xlabel = "r [m]", ylabel = "εrr [%]")
            # lines!(ax2, rc, εrr.*100., label="ρ [kg.m-3]")
            # -------------- #
            # ax3 = Axis(f[2, 2], title = "τII", xlabel = "r [m]", ylabel = "τII [MPa]")
            # lines!(ax3, rc, τII./1e6, label="τII [MPa]") 
            ax3 = Axis(f[2, 2], title = "P", xlabel = "r [m]", ylabel = "P [GPa]")
            lines!(ax3, rc, P./1e9, label="P [GPa]") 
            lines!(ax3, rc, .-(σrr .+ σϕϕ .+ σΘΘ)./3.0./1e9, label="P [GPa]")
            # -------------- #
            # ax4 = Axis(f[4, 1], title = "Vr", xlabel = "r [m]", ylabel = "Vr [m/s]")
            # lines!(ax4, rv, Vr, label="Vr [m/s]") 
            ax4 = Axis(f[3, 1], title = "Divergence evolution", xlabel = "t [ky]", ylabel = "∇vbg [1/s]")
            div1= copy(div_vec[1:it])
            pl1 = copy(ispl_vec[1:it])
            t1  = copy(t_vec[1:it])
            lines!(ax4, t1./kyr, log10.(abs.(div1)), label="log₁₀ ∇vbg [1/s]") 
            scatter!(ax4, t1[pl1.>1e-10]./kyr,  log10.(abs.(div1[pl1.>1e-10])), color=:red)   
            # -------------- #
            ax5 = Axis(f[3, 2], title = "Pressure evolution", xlabel = "t [ky]", ylabel = "P [GPa]")
            p1  = copy(p_mat_vec[1:it])
            pl1 = copy(ispl_vec[1:it])
            t1  = copy(t_vec[1:it])
            lines!(ax5, t1./kyr, p1./1e9, color=:blue, label="" , marker=0, linewidth=1.0)   
            scatter!(ax5, t1[pl1.>1e-10]./kyr, p1[pl1.>1e-10]./1e9, color=:red)   
            # -------------- #
            display(f)

            if BC==:ConstantPressureRate
                name = "$(BC)$(@sprintf("%04d", Int64(dPdtbg*(1e-9*(1e3*kyr)))))"
            else
                name = "$(BC)$(@sprintf("%1.0e", ((ε̇bg))))"
            end

            save("img/"*name*"_ts$(@sprintf("%04d", it)).png", f, px_per_unit = 4)
        
            file = matopen("out/"*name*".mat", "w")
            write(file, "p_mat_vec", p_mat_vec)
            write(file, "ispl_vec", ispl_vec)
            write(file, "div_vec", div_vec)
            write(file, "t_vec", t_vec)
            close(file)

            # Stop run
            if maximum(P)<1e9
                break
            end
        end
    end
end # End function

main()