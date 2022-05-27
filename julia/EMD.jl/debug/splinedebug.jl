
## Bad input from doppler (end not-a-knot conditions, i.e. n ≤ 3)
x = [-493, 1417, 4887]
v = [0.0072, 0.0072, 0.0072]
t = collect(1:4096)

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