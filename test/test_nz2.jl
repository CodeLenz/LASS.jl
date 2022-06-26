
@testset " Lass / dLass (Nz=2)" begin

    # Number of realizations
    n = 1_000_000

    # Divisions in each dimension
    Nb = [50;50]

    # Load problem
    include("probnz2.jl")

    # Densities
    vp = []
    push!(vp, LogNormal(log(3.0), 0.1) )
    push!(vp, LogNormal(log(10.0),0.1) )

    # Realizations
    distrib = zeros(2,n)
    for i=1:2
        distrib[i,:] .= rand(vp[i],n)
    end

    # Make bins
    bins = Generate_bins(distrib, Nb)

    # Evaluate E[f] and Var[f]
    E, Var = Lass(bins,f)

    # Evaluate dE[f] and dVar[f]
    dE, dVar = dLass(bins,f, dfx,1)

    # Reference values
    refE,refVar,refdE,refdVar = Referencias()

    # Tests
    @test isapprox(E,refE,atol=1E0)

    @test isapprox(dE[1],refdE,atol=1E0)

    @test isapprox(Var,refVar,atol=1E0)

    @test isapprox(dVar[1],refdVar,atol=1E0)


end