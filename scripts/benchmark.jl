# Preamble ########################################################################################
using Pkg
Pkg.activate(Base.current_project())

cd("..") # Go to the parent directory
###################################################################################################
# Packages needed #################################################################################
using UnitRots
using DelimitedFiles
using ImageCore
using ImageIO
using ImageMagick
using FileIO
using PyPlot
using BenchmarkTools
###################################################################################################
# Angle of rotation ###############################################################################

const θ = π/4

###################################################################################################
# Retrieve of data from image #####################################################################

@info "Retrieving the image and extracting the two-dimensional field"
input = Int64.(channelview(Gray.(load("data/benchmark/delta31.png"))))

###################################################################################################
# Identification of the grayscale field dimensions ################################################

j = (size(input)[1] - 1)/2;

N = Integer(2*j + 1);

@info "The size of the field is $N"

###################################################################################################
# Rotation of the image ###########################################################################

@info "Performing rotation at θ = $θ"

output = real(Emodes(input, θ))

###################################################################################################
begin
	fig, (ax1, ax2) = subplots(1, 2, figsize = (20,10), constrained_layout = true)

	ax1.imshow(real(input), cmap = "gray")
    ax1.set_xticks([-0.5:1:N;])
    ax1.set_yticks([-0.5:1:N;])
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.tick_params(axis = "both", length = 0)
    ax1.grid(which = "major", color = "gray", linestyle = "-", linewidth = 1.5, alpha = 0.5)

    ax2.imshow(nnorm(output), cmap = "gray")
    ax2.set_xticks([-0.5:1:N;])
    ax2.set_yticks([-0.5:1:N;])
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.tick_params(axis = "both", length = 0)
    ax2.grid(which = "major", color = "gray", linestyle = "-", linewidth = 1.5, alpha = 0.5)

    savefig("figures/benchmark.png", dpi = 300)
end
