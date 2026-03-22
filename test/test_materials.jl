@testset "Materials & Constants" begin
    # Magnetic length at B = 1 T should be ≈ 25.66 nm
    ℓ_B = magnetic_length(1.0)
    @test isapprox(ℓ_B * 1e9, 25.66, atol=0.1)

    # Cyclotron energy at B = 1 T for GaAs: ħω_c/k_B ≈ 20 K
    ħω_c = cyclotron_energy(1.0)
    T_cyc = ħω_c / CONSTANTS.k_B
    @test isapprox(T_cyc, 20.0, atol=2.0)

    # Filling factor: ν = 1/3 at B ≈ 12.4 T for n_e = 1.0e15 m⁻² (= 1.0×10¹¹ cm⁻²)
    n_e = 1.0e15
    ν = filling_factor(12.4, n_e)
    @test isapprox(ν, 1/3, atol=0.02)

    # Inverse: B at ν = 1/3
    B13 = B_at_filling(1/3, n_e)
    @test isapprox(B13, 12.4, atol=0.5)

    # Integer QHE: R_xy at ν = 1 = h/e² ≈ 25812.807 Ω
    R1 = iqhe_Rxy(1)
    @test isapprox(R1, 25812.807, atol=0.01)

    # Coulomb energy scale at B = 10 T: E_C ≈ 160 K (e²/(4πε₀ε_r ℓ_B))
    E_C = coulomb_energy(10.0)
    E_C_K = E_C / CONSTANTS.k_B
    @test isapprox(E_C_K, 160.0, atol=10.0)
end
