module subst

using Calculus
# define types first
struct substance
    name::String
    density::Function
    enthalpy::Function
    heat_capacity::Function
    viscosity::Function
    thermal_conductivity::Function
    alt_heat_capacity::Function
    visc_kin::Function
    exp_coeff::Function
    therm_diff::Function
    Pr::Function
end

function substance(name, x::Vector{Function})
    dens, enth, Cp, visc, thermCond = x
    function alt_heat_capacity(p, T)
        return derivative(T_ -> enth(p, T_), T)
    end
    function visc_kin(p, T)
        return visc(p, T)/dens(p, T)
    end
    function exp_coeff(p, T)
        dRdT = derivative(T_ -> dens(p, T_), T)
        return (-1/dens(p, T))*dRdT
    end
    function thermDiff(p, T)
        if Cp == nothing
            Cp_ = alt_heat_capacity(p, T)
        else
            Cp_ = Cp(p, T)
        end
        return thermCond(p, T) / (dens(p,T)*Cp_)
    end
    function Pr(p, T)
        if Cp == nothing
            Cp_ = alt_heat_capacity(p, T)
        else
            Cp_ = Cp(p, T)
        end
        return Cp_ * visc(p, T) / thermCond(p, T)
    end
    return substance(name, dens, enth, Cp, visc, thermCond, alt_heat_capacity, visc_kin, exp_coeff, thermDiff, Pr)
end

function na_enthalpy(p, T)
    """This function takes the temperature of liquid sodium in K and returns the enthalpy in kJ/kg.
    The reference state is taken as solid sodium at 298 K. The correlation takes this form:
    H(l, T) - H(s, 298.15) = 365.77 + 1.6582 T - 4.2395 × 10^-4 * T^2 + 1.4847 × 10^-7 * T^3 + 2992.6 * T^-1
    This correlation is recommended by a report prepared by Argonne National Lab that can be found here:
    http://www.ne.anl.gov/eda/ANL-RE-95-2.pdf
    The correlation is valid between 371 and 2000 K."""
    a = -365.77
    b = 1.6582
    c = 4.2395E-4
    d = 1.4847E-7
    e = 2992.6
    return (a + b * T - c * T ^ 2 + d * T ^ 3 + e / T) * 1000
end

function na_Cp(p, T)
    """This function takes the temperature of liquid sodium and returns the heat
    capacity at constant pressure in kJ/kg*K. The Argonne report gives a very
    complicated correlation for Cp but notes that it's almost equal to the
    derivative of enthalpy with temperature up to about 1800 K. So I used the derivative here."""
    a = -365.77
    b = 1.6582
    c = 4.2395E-4
    d = 1.4847E-7
    e = 2992.6
    return (b - 2*c * T + 3 * d * T ^ 2 - e / (T^2)) * 1000
    #return derivative(enthalpy, T, dx = 1E-6)
end

function na_density(p, T)
    """This function takes the temperature of liquid sodium in K and returns the density in kg/m^3.
    The correlation takes this form:
    rho = rho_c + f * (1 - Tr) + g * (1 - Tr)^h
    where rho_c is the critical density of sodium, Tr is the reduced temperature T/T_c,
    and f g and h are constants given in the function.
    This correlation is recommended by a report prepared by Argonne National Lab that can be found here:
    http://www.ne.anl.gov/eda/ANL-RE-95-2.pdf
    The correlation is valid between 371 and ~2000 K."""
    rho_c = 219
    f = 275.32
    g = 511.58
    h = 0.5
    Tc = 2503.7
    Tr = T / Tc
    return rho_c + f * (1 - Tr) + g * (1 - Tr) ^ h
end

function na_thermal_conductivity(p, T)
    """This function takes the temperature of liquid sodium in K and returns the thermal conductivity in W/m*K.
    The reference state is taken as solid sodium at 298 K. The correlation takes this form:
    k = 124.67 - 0.11381 * T + 5.5226 × 10^-5 * T^2 - 1.1842 × 10^-8 * T^3
    This correlation is recommended by a report prepared by Argonne National Lab that can be found here:
    http://www.ne.anl.gov/eda/ANL-RE-95-2.pdf
    The correlation is valid between 371 and 2000 K."""
    a = 124.67
    b = 0.11381
    c = 5.5226E-5
    d = 1.1842E-8
    return a - b * T + c * T ^ 2 - d * T ^ 3
end

function na_dynamic_viscosity(p, T)
    """This function takes the temperature of liquid sodium in K and returns the dynamic viscosity in Pa * s.
    The reference state is taken as solid sodium at 298 K. The correlation takes this form:
    ln eta = - 6.4406 - 0.3958 ln (T) + 556.835 / T
    This correlation is recommended by a report prepared by Argonne National Lab that can be found here:
    http://www.ne.anl.gov/eda/ANL-RE-95-2.pdf
    The correlation is valid between 371 and 2500 K."""
    a = -6.4406
    b = -0.3958
    c = 556.835
    return exp(a + b * log(T) + c / T)
end

sodium_subst = substance("sodium", [na_density, na_enthalpy, na_Cp, na_dynamic_viscosity, na_thermal_conductivity])



"""
Struct material
Used for solid materials in pipes. So far includes just name, thermal
conductivity, and surface roughness.
"""
struct material
    name::String
    thermal_conductivity::Function
    roughness::Real
end

function stainless_steel_thermal_conductivity(T)
    # Source: http://www.mace.manchester.ac.uk/project/research/structures/strucfire/materialInFire/Steel/StainlessSteel/thermalProperties.htm
    return 14.6 + 1.27e-2 *(T-273.15)
end

function copper_thermal_conductivity(T)
    # Source: http://www-ferp.ucsd.edu/LIB/PROPS/PANOS/cu.html
    return 14.6 + 1.27e-2 *(T-273.15)
end

ss_roughness_by_finish = Dict("2D" => 1E-6,
                              "2B" => 0.5E-6,
                              "2R" => 0.2E-6,
                              "BA" => 0.2E-6,
                              "2BB" => 0.1E-6)
# These surface roughness values are in m and were obtained from the following
# corporate site: http://www.outokumpu.com/en/products-properties/
# more-stainless/stainless-steel-surface-finishes/cold-rolled-finishes/Pages/default.aspx
# Note that the given values are the maximum roughness values in the range on
# the website, except for finish 2BB, which didn't have a roughness value
# so I took the lowest value for finish 2B since the compnay says this:
# "Due to the fact that surface roughness of the finish 2BB is lower than that
# of 2B, some modifications on lubrication during forming might be needed."

ss_mat = material("Stainless Steel", stainless_steel_thermal_conductivity,
                  ss_roughness_by_finish["2B"])

cu_roughness = 0.03e-3
# Copper roughness.
# Source: http://www.pressure-drop.com/Online-Calculator/rauh.html

cu_mat = material("Copper", copper_thermal_conductivity, cu_roughness)

end # ends module subst

module pr

using Polynomials.Poly, Polynomials.roots
using Roots.fzero

function R()
  return 8.31446
end

struct pengRobinson{T<:Real}
    pc::T
    Tc::T
    a::T
    b::T
    ω::T
    κ::T
    α::Function
    A::Function
    B::Function
    c1::Function
    c2::Function
    c3::Function
end

function pengRobinson(pc, Tc, ω)
  a = 0.45724 * R() ^2 * Tc^2 /pc
  b = 0.0778*R() * Tc /pc
  κ = 0.37464 + 1.54226*ω - 0.26992*ω^2
  function α(T)
    Tr = T/Tc
    return (1+κ*(1-√(Tr)))^2
  end
  function A(p, T)
    return a * α(T) * p / (R() ^2 * T^2)
  end
  function B(p, T)
    return b*p/(T*R())
  end
  function c1(p, T)
    return (B(p,T) - 1)
  end
  function c2(p, T)
    return A(p, T) - 2B(p,T) - 3B(p,T)^2
  end
  function c3(p,T)
    return B(p,T)^3 + B(p,T)^2 - A(p,T)*B(p,T)
  end
  return pengRobinson(pc, Tc, a,b,ω,κ,α,A,B,c1,c2,c3)
end

function get_reals(rvec)
  ans = []
  for num in rvec
    if imag(num) == 0
      ans = push(ans, num)
    end
  end
  return ans
end

function z_roots(p_::pengRobinson, p, T)
  pvec = Poly([p_.c3(p,T),p_.c2(p,T),p_.c1(p,T),1])
  return get_reals(roots(pvec))
end

function v_roots(p_::pengRobinson, p, T)
  return z_roots(p_, p, T) * R() * T / p
end

function pressure(p_::pengRobinson, V, T)
  return R()*T / (V-p_.b) - p_.a*p_.α(T) / (V^2 + 2*V*p_.b - p_.b^2)
end

function vv(p_::pengRobinson, p, T)
  vg = R()*T/p
  zer = fzero(V_ -> pressure(p_, V_, T) - p, vg)
  return zer
end

function vl(p_::pengRobinson, p, T)
  vg = 1.1*p_.b
  zer = fzero(V_ -> pressure(p_, V_, T) - p, vg)
  return zer
end

end #end module pr

import pr
function __main__()
  n2 = pr.pengRobinson(34.0e5, 126.2, 0.0377215)
  @printf "P calc: %f" pr.pressure(n2, 0.024778451307952636, 298.15)
  @printf "N₂ a: %f" n2.a
  @printf "N₂ κ: %f" n2.κ
  @printf "N₂ b: %f" n2.b
  @printf "N₂ a: %f" n2.a
  @printf "N₂ a: %f" n2.a
  @printf "P calc: %f" pr.pressure(n2, 0.024778451307952636, 298.15)
  @printf "P calc: %f" pr.pressure(n2, 0.024778451307952636, 298.15)
  @printf "P calc: %f" pr.pressure(n2, 0.024778451307952636, 298.15)
  @printf "Vapor volume: %f" pr.vv(n2, 1.0e5, 298.15)

end
#__main__()
