# Global model parameters


N = 100
D = 2
R = 50


# Species parameters

# u_c : matrix of dimension D X R with niche optima for each species and environmental variable for colonization

u_c = matrix(nr = D, nc = R)
u_c[1,] = seq(0.1,0.9,length=R)
u_c[2,] = rep(0.5, R)

# u_e : matrix of dimension D X R with niche optima for each species and environmental variable for extinction

u_e = matrix(nr = D, nc = R)
u_e[1,] = seq(0.1,0.9,length=R)
u_e[2,] = rep(0.5, R)

# s_c : matrix of dimension D X R with niche breadth for each species and environmental variable

s_c = matrix(0.5, nr = D, nc = R)

# s_e : matrix of dimension D X R with niche breadth for each species and environmental variable

s_e = matrix(0.5, nr = D, nc = R)

# A: matrix of dimension R X R with the effect of each species (column) on other species (rows)

A = matrix(-1, nr = R, nc = R)
diag(A) = 0

pars = list(u_c = u_c, 
            u_e = u_e, 
            s_c = s_c, 
            s_e = s_e, 
            alpha = 0.3, # alpha : scalar, average dispersal distance
            m = 0.001, # m : scalar, immigration probability
            
            c_0 = rep(1, R), 
            e_0 = rep(0.05, R), # e_0 : vector of dimension R with extinction probability at 0 interaction
            c_max = rep(0.8, R), # c_max: vector of dimension R maximal colonization probability at infinite positive interactions
            e_min = rep(0.05, R), # e_min: vector of dimension R maximal extinction probability at infinite positive interactions
            
            d_c = 0.5, # d_c: scalar, strength of the effect of interactions on the colonization probability
            d_e = 0.5, # d_e: scalar, strength of the effect of interactions on the extinction probability
            A = A, 
            SA = NULL # SA: spatial autocorrelation - Null defaults to no spatial autocorrelation. If a numeric value is provided, it will create an Environment with SA of the given length. 
)








