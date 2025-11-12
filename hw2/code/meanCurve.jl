using Arpack
using Printf
using SparseArrays

include("utils.jl")

filename = "meshes/moomoo.off"
# filename = "meshes/166.off"
X, T = readoff(filename)
nv = size(X, 1)
nt = size(T, 1)

# ADD CODE TO COMPUTE SURFACE AREA HERE #############
function surfaceArea(X, T)
    A = 0
    for i in axes(T,1)
        a = X[T[i,2],:] - X[T[i,1],:]
        b = X[T[i,3],:] - X[T[i,1],:]
        A += norm( cross( a,b ) ) / 2
    end
    return A
end
# END HOMEWORK ASSIGNMENT #########

@printf("The surface area of %s is %f\n", filename, surfaceArea(X, T))

# ADD CODE TO COMPUTE COTANGENT LAPLACIAN HERE ##########
function cotangent( u :: Vector{Float64}, v :: Vector{Float64} ) :: Float64
    return dot( u,v ) / norm( cross( u,v ) )
end

function cotLaplacian(X, T)
    L = spzeros(nv,nv)
    for t in axes(T,1)
        for e in 1:3
            vert = circshift( 1:3, e )
            i = T[t,vert[1]]
            j = T[t,vert[2]]
            k = T[t,vert[3]]
            w = cotangent( X[i,:] - X[k,:], X[j,:] - X[k,:] )
            L[i,j] -= w
            L[j,i] -= w
            L[i,i] += w
            L[j,j] += w
        end
    end
    return L
end
# END HOMEWORK ASSIGNMENT #######

# Sanity checks: Laplacian is symmetric and positive definite
L = cotLaplacian(X, T)
eigenvals, _ = eigs(L, nev=10, which=:SM)
println(eigenvals)
println(norm(L - L'))

# ADD CODE FOR DIVIDED DIFFERENCES HERE ######
function ∇A( X, T, d )
    h = 0.001
    d = vec(d)
    res = zeros( size( X,1 ) )
    for i in axes( X, 1 )
        Xᵢ = @view X[i,:]
        Xᵢ .+= h * d
        res[i] += surfaceArea( X, T ) / 2h
        Xᵢ .-= 2h * d
        res[i] -= surfaceArea( X, T ) / 2h
        Xᵢ .+= h * d
    end
    return res
end

function dividedDifferences(X, T)
    gradApprox = zeros(nv, 3);
    gradApprox[:,1] = ∇A( X,T, [1 0 0] )
    gradApprox[:,2] = ∇A( X,T, [0 1 0] )
    gradApprox[:,3] = ∇A( X,T, [0 0 1] )
    return gradApprox
end
# END HOMEWORK ASSIGNMENT #########

# Check that gradApprox and .5*L*X are similar
#println(norm(.5 .* (L*X) - dividedDifferences(X,T)))

# ADD CODE FOR COMPUTING THE BARYCENTRIC AREA VECTOR HERE ######
function barycentricArea(X, T)
    A = zeros(nv)
    for t in eachrow( T )
        u = X[t[2],:] - X[t[1],:]
        v = X[t[3],:] - X[t[1],:]
        a = norm( cross( u,v ) ) / 2
        for i in t
            A[i] += a / 3
        end
    end
    return A
end
# END HOMEWORK ASSIGNMENT ########

# ADD CODE FOR COMPUTING POINTWISE MEAN CURVATURE HERE #####
function meanCurvature(X, T)
    A = barycentricArea( X,T )
    L = cotLaplacian(X,T)
    H = norm.( eachrow( 0.5 .* L * X ./ A ) )
    return H
end
# END HOMEWORK ASSIGNMENT #######

H = meanCurvature(X, T)
scene = showdescriptor(X, T, H)


## Mean curvature flow ##
function curvatureFlowEuler(X, T)
    Xt = copy(X)
    maxiters = 10
    # ADD CODE FOR THE EXPLICIT INTEGRATOR HERE ####
    τ = 0.5
    for t=1:maxiters
        M⁻ = spdiagm( 1 ./ barycentricArea( Xt, T ) )
        L  = cotLaplacian( X, T ) 
        Xt -= τ .* M⁻ * L * Xt
    end
    # END HOMEWORK ASSIGNMENT #####
    H = meanCurvature(Xt, T)
    # Uncomment to show mean curvature at the end
    scene = showdescriptor(X, T, H)
end

function curvatureFlowImplicit(X, T)
    Xt = copy(X)
    maxiters = 100
    # ADD CODE FOR SEMI-IMPLICIT INTEGRATOR HERE ####
    τ = 0.05
    spId = sparse( I, nv,nv )
    for t=1:maxiters
        M⁻ = spdiagm( 1 ./ barycentricArea( Xt, T ) )
        L  = cotLaplacian( X, T ) 
        Xt = ( spId + τ .* M⁻ * L ) \ Xt
    end
    # END HOMEWORK ASSIGNMENT
    H = meanCurvature(Xt, T)
    # Uncomment to show mean curvature at the end
    scene = showdescriptor(X, T, H)
end

# Uncomment either the implicit or explicit flow to see your results
curvatureFlowEuler(X,T)
curvatureFlowImplicit(X,T)
