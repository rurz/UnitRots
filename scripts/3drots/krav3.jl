# Preamble ########################################################################################
using Pkg
Pkg.activate(Base.current_project())

#cd("..") # Go to the parent directory # Activate when using on terminal bash
###################################################################################################
# Packages needed #################################################################################
using UnitRots
using DelimitedFiles

using PyPlot
ion()
pygui(true)
###################################################################################################

k2(j, n₁, n₂) = kmatrix(j)[:, n₁] * transpose(kmatrix(j)[:, n₂])

function k3(j, n₁, n₂, n₃)
    m = zeros(2*j+1, 2*j+1, 2*j+1)
    for q₁ in -j:j
        for q₂ in -j:j
            for q₃ in -j:j
                m[q₁ + j + 1, q₂ + j + 1, q₃ + j + 1] = kmatrix(j)[q₁ + j + 1, n₁ + 1] * kmatrix(j)[q₂ + j + 1, n₂ + 1] * kmatrix(j)[q₃ + j + 1, n₃ + 1]
            end
        end
    end
    return m
end