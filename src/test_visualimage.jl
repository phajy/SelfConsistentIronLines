using Gradus, Plots

m=KerrMetric(M=1.0,a=0.9)


d1 = ShakuraSunyaev(m, eddington_ratio = 0.1)
d2 = ShakuraSunyaev(m, eddington_ratio = 0.2)
d3 = ShakuraSunyaev(m, eddington_ratio = 0.3)


function Gradus.cross_section(d, ρ)
    y=Gradus.isco(m)
    if ρ < y
        return -one(typeof(ρ))
    end
    H = (3 / 2) * inv(η) * (d.Ṁ / d.Ṁedd) * (1 - sqrt(y / ρ))
end

r_start=-50
r_end=50
n=r_end-r_start
r=LinRange(r_start,r_end,n)


function all_height(disk)
    height=zeros(r_end)
    height_b=zeros(r_end)
    for i in 1:r_end
        height[i]=cross_section(disk,i)
        height_b[i]=0-cross_section(disk,i)
    end
    a=reverse(height)
    b=reverse(height_b)
    top=vcat(a,height)
    bottom=vcat(b,height_b) 
    return top, bottom
end

d_h1,d_h1b=all_height(d1)
d_h2,d_h2b=all_height(d2)
d_h3,d_h3b=all_height(d3)

plot(r,[d_h1],fillrange=d_h1b,alpha=0.2,label="Edd=0.1",title="Disk Thickness")
plot!(r,[d_h2],fillrange=d_h2b,alpha=0.2,label="Edd=0.2")
plot!(r,[d_h3],fillrange=d_h3b,alpha=0.2,label="Edd=0.3")





