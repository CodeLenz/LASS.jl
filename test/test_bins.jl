
@testset "Test bins (nz=2)" begin

    # Divisions in each dimension
    Nb = [3;3]
    #                   |      4
    #             1     |      
    #                 2 |  
    #          ---------|---------------->
    #                   |  3
    #                   |
    #                   | 
    # Non random points in R²
    distrib = [-1.0  -0.2   0.3  0.6  -1.0  0.6 ;
                1.0  -0.2  -0.1  1.2  -0.1  1.2]
  
    # As the limits are
    # x in [-1 0.6]
    # y in [-0.2 1.2]          
    # 
    # and we are using a 3x3 grid, each division is
    #  
    # dx = 1.6/3 = 0.533333333333
    # dy = 1.4/3 = 0.466666666666 
    #
    # Such that we have 9 bins
    #
    # First "row"
    # b1 => x=[-1:-0.466666]       and y=[-0.2:0.266666]
    # b2 => x=[-0.466666:0.066666] and y=[-0.2:0.266666]
    # b3 => x=[0.066666:0.6]       and y=[-0.2:0.266666]
    # 
    # Second "row"
    # b4 => x=[-1:-0.466666]       and y=[0.266666:0.733333]
    # b5 => x=[-0.466666:0.066666] and y=[0.266666:0.733333]
    # b6 => x=[0.066666:0.6]       and y=[0.266666:0.733333]
    #
    # Third "row"
    # b7 => x=[-1:-0.466666]       and y=[0.733333:1.2]
    # b8 => x=[-0.466666:0.066666] and y=[0.733333:1.2]
    # b9 => x=[0.066666:0.6]       and y=[0.733333:1.2]
    #
    #
    # Thus, the map from points to bins should be
    #
    # P1 -> b7
    # P2 -> b2
    # P3 -> b3
    # P4 -> b9
    # P5 -> b1
    # P6 -> b9
    #
    #
    # Dict{Int64, NBin_data} with 5 entries:
    #  5 => NBin_data{Float64}(9, [0.6, 1.2], 0.333333)
    #  4 => NBin_data{Float64}(7, [-1.0, 1.0], 0.166667)
    #  2 => NBin_data{Float64}(2, [-0.2, -0.2], 0.166667)
    #  3 => NBin_data{Float64}(3, [0.3, -0.1], 0.166667)
    #  1 => NBin_data{Float64}(1, [-1.0, -0.1], 0.166667)
    
    # Make bins
    bins = Generate_bins(distrib, Nb)
  
    # First test
    @test length(bins)==5

    # Each bin
    for i=1:5
       
       ID = bins[i].ID
       μ  = bins[i].μ
       w  = bins[i].w

       # Tests
       if ID==1
          @test all(μ.==vec(distrib[:,5]))
          @test isapprox(w,1/6)
       elseif ID==2
            @test all(μ.==vec(distrib[:,2]))
            @test isapprox(w,1/6)
       elseif ID==3
            @test all(μ.==vec(distrib[:,3]))
            @test isapprox(w,1/6)
       elseif ID==7
            @test all(μ.==vec(distrib[:,1]))
            @test isapprox(w,1/6)
       elseif ID==9
            pos = [4;6]
            @test all(μ.==mean(distrib[:,pos],dims=2))
            @test isapprox(w,2/6)
       end    
             

    end

end