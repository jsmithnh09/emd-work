
## Bad input from doppler (end not-a-knot conditions, i.e. n ≤ 3)
X = [-493, 1417, 4887]
V = [0.0072, 0.0072, 0.0072]
Xq = collect(1:4096)




coeffs = vcat(V[1], dvdx)
coeffs[3] = (diff(dvdx) / dx2)[1]
coeffs[2] -= coeffs[3]*dx[1]


M = length(Xq)
dx, dv = diff(X), diff(V)
dx2 = X[3] - X[1] # skipping second "tau" point
dvdx = dv ./ dx
coeffs = zeros(3) # ax² + bx + c
coeffs[1] = (diff(dvdx) / dx2)[1] # d²v / dx²
coeffs[2] = dvdx[1] - coeffs[1]*dx[1]
coeffs[3] = V[1]
Vq = fill(coeffs[1], M)
xs = collect(Xq) .- fill(X[1], M)
for i = 2:3
    Vq = Vq .* xs + fill(coeffs[i], M)
end





N = length(X)
M = length(t)
dx, dv = diff(x), diff(v)
dx2 = x[3] - x[1]
dvdx = dv ./ dx
Vc = vcat(V[1], dvdx)
Vc[3] = (diff(dvdx) / dx2)[1] # d²v / dx²
Vc[2] -= Vc[3]*dx[1]

bounds = (x[1], x[3])
coeffs = Vc[end:-1:1]

## Perform the piecewise interpolation for k=3
xs = collect(t) .- fill(bounds[1], M)
vf = fill(coeffs[1], M)
for i = 2:3
    vf = vf .* xs + fill(coeffs[i], M)
end