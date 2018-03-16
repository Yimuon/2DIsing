#2D Classical Ising Model
#loop summation with modified random walk
#Yimu Bao

#transfer matrix from (x_start, y_start) to (x_end, y_end) in l steps (<x1,y1|W^*(l)|x2,y2>)
#x_start, y_start, x_end, y_end are the coordinates of starting and ending sites.
#mu_start, mu_end denote the direction of next move (mu = 1,2,3,4)
#1: right; 2: up; 3: left; 4: down.
#l is the step that connects two sites
#theta*pi/2 is the local change of angle by the tangent vector. (theta = 0,1,2,3)

function trsfmat(x_start, y_start, mu_start, x_end, y_end, mu_end, l, theta)
    #invalid input
    if x_start <= 0 || x_start > Nx || y_start <= 0 || y_start > Ny ||
        x_end <= 0 || x_end > Nx || y_end <= 0 || y_end > Ny
        return 0
    elseif mu_start <=0 || mu_start > 4 || mu_end <= 0 || mu_end > 4 ||
        l < 1 || theta < 0 || theta >= 8
        println("Invalid input")
        return 0
    end

    #Return
    if abs(x_end-x_start)+abs(y_end-y_start) > l
        return 0
    end

    #Next step
    if mu_start == 1
        x_next = x_start + 1
        y_next = y_start
    elseif mu_start == 2
        x_next = x_start
        y_next = y_start + 1
    elseif mu_start == 3
        x_next = x_start - 1
        y_next = y_start
    elseif mu_start == 4
        x_next = x_start
        y_next = y_start - 1
    else
        println("Error mu_start!!!")
        return 0
    end

    Nloop = 0
    if l == 1
        if x_next == x_end && y_next == y_end
            if mod(mu_end-mu_start, 4) == 2
                return 0
            else
                if mod(theta+(mod(mu_end-mu_start+2, 4)-2), 8) == 0
                    return -0.5 #avoid double counting of directed loops
                elseif mod(theta+(mod(mu_end-mu_start+2, 4)-2), 8) == 4
                    return 0.5
                else
                    println("Not a closed loop")
                    return 0
                end
            end
        else
            #Random walk terminated
            return 0
        end
    else
        for mu_next = 1:4
            #No U-turn condition
            if mod(mu_next - mu_start, 4) != 2
                theta_next = mod(theta+(mod(mu_next-mu_start+2, 4)-2), 8)
                Nloop += trsfmat(x_next, y_next, mu_next, x_end, y_end, mu_end, l-1, theta_next)
            end
        end
        return Nloop
    end
end

#Free energy per site: betaf = betaF/N = ln(Z)/N
#K: beta*t
#lmax: maximal length of loop
function freeenergy(K, lmax)
    t = tanh(K)
    x0 = floor(Nx/2)
    y0 = floor(Ny/2)
    betaf = 0#log(2cosh(K)^2)+1
    for mu = 1:4
        for l = 2:2:lmax
            betaf += t^l/l*trsfmat(x0,y0,mu,x0,y0,mu,l,0)
            # println(l)
            # println(betaf)
        end
    end
    return betaf
end

#Correlation function
#additional minus sign is due to the convention of positive loop
#when computing free energy, simple loop contribute +1,
#which is to the opposite when computing correlation function
function corrfunc(x1, y1, x2, y2, K, lmax)
    t = tanh(K)
    g = 0
    if mod(x2-x1+y2-y1, 2) == 0
        for l = 2:2:lmax
            for mu = 1:4
                g += -t^l*trsfmat(x1,y1,mu,x2,y2,mu,l,0)
            end
        end
    elseif mod(x2-x1+y2-y1, 2) == 1
        for l = 1:2:lmax
            for mu = 1:4
                g += -t^l*trsfmat(x1,y1,mu,x2,y2,mu,l,0)
            end
        end
    else
        println("Error input coordinates")
    end
    return g
end

#Const
Nx = 41
Ny = 41
K = log(sqrt(2)+1)/2 #beta*t = J/kBT
# t = tanh(K)
# tc = sqrt(2)-1

# #Plot correlation function
# r = 1:10
# gcorr = zeros(10)
# for i = 1:10
#     gcorr[i] = corrfunc(15,21,15+i,21,K,lmax)
# end
# using PyPlot
# plot(collect(r),gcorr)

#Extrapolation method
n = 4:2:20
#Coeff of expansion: f = \sum a_n t^n
a = zeros(size(n))
for i = 1:size(n,1)
    a[i] = trsfmat(21,21,1,21,21,1,n[i],0)*4/n[i]
end
#a_n/a_{n-1}
y = a[2:end]./a[1:end-1]
#1/n
x = 1./n[2:end]
using PyPlot
scatter(x,y)
plot(x,y,color = "orange")

#Critical tc & Kc by extrapolation
tc_inv_sq = x[end]/(x[end-1]-x[end])*(y[end]-y[end-1])+y[end]
tc = 1/sqrt(tc_inv_sq)
println("Critical value of t: tc = $(tc)")
println("Theoretical value of tc: tc_thy = $(sqrt(2)-1)")
Kc = atanh(tc)
println("Critical value of K: Kc = $(Kc)")

#a_n: 4:2:22
#     1.0
#     2.0
#     4.5
#    12.0
#    37.3333
#   130.0
#   490.25
#  1958.67
#  8174.2
# 35302.0

#Critical exponent
#alpha
