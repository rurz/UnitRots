
export clip, nnorm, adiag, minv, maxv

function clip(M::Matrix)
    v = zeros(Float64, size(M))
    for i in 1:size(M)[1]
        for j in 1:size(M)[1]
            if M[i,j] < 0
                v[i, j] = 0
            elseif M[i,j] > 1
                v[i, j] = 1
            else
                v[i, j] = M[i, j]
            end
        end
    end
    return v
end

function nnorm(M::Matrix)
    v = zeros(Float64, size(M))
    mv = minimum(M)
    Mv = maximum(M)
    for i in 1:size(M)[1]
        for j in 1:size(M)[1]
            v[i, j] = (M[i, j] - mv) * ((1 - 0)/(Mv - mv)) + 0
        end
    end
    return v
end

function adiag(M::Matrix)
	v = zeros(size(M)[1])
	for i in 1:(size(M)[1])
		v[i] = M[(size(M)[1] + 1) - i , i]
	end
	return v
end

xax(M::Matrix) = range(1, size(M)[1], step = 1)

minv(M::Matrix) = minimum(adiag(M))
maxv(M::Matrix) = maximum(adiag(M))

