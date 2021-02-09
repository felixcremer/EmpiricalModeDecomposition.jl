using EmpiricalModeDecomposition
using Makie
using Statistics
using ColorBrewer
const EMD = EmpiricalModeDecomposition

maxcol, mincol, meancol = palette("Set2",3)

ts, x , data = EmpiricalModeDecomposition.maketestdata(123)
ts, x, data = EMD.fosso()
x = imfs[2]


maxes = Int[]
mins = Int[]
yvec = Node(x)
t = Node(1.0)
lw = @lift(3/$t)
mw = @lift(5/$t)
EMD.localmaxmin!(x, maxes, mins)
y_max = Node(x[maxes])
x_max = Node(ts[maxes])
scene = lines(ts, yvec, limits = FRect(0, minimum(x), maximum(ts),maximum(x)-minimum(x)))
st = Stepper(scene, "/home/crem_fe/Documents/presentations/MPI_20201125/sift_animation/")
step!(st)
smin = EMD.get_edgepoint(x, ts, mins, first, isless)
smax = EMD.get_edgepoint(x, ts, maxes, first, !isless)
emin = EMD.get_edgepoint(x, ts, mins, last, isless)
emax = EMD.get_edgepoint(x, ts, maxes, last, !isless)
 
maxTS = EMD.interpolate([first(ts); ts[maxes]; last(ts)],[smax; x[maxes]; emax],ts,EMD.DataInterp())

minTS = EMD.interpolate([first(ts); ts[mins]; last(ts)],[smin; x[mins]; emin],ts,EMD.DataInterp())

#meanTS = mean.(maxTS, minTS)

maxts = Node(maxTS)
mints = Node(minTS)
scatter!(x_max,y_max, color=maxcol, markersize=mw)
step!(st)

lines!(ts, maxts, color=maxcol, linewidth=lw)
step!(st)
y_min = Node(x[mins])
x_min = Node(ts[mins])
scatter!(x_min,y_min, color=mincol, markersize=mw)
step!(st)
lines!(ts, mints, color=mincol, linewidth=lw)
step!(st)

meants = @lift(($mints .+ $maxts) ./ 2)
lines!(ts, meants, color=meancol, linewidth=lw)
step!(st)
ϵ = var(x) * 0.01
imfs = EMD.eemd(x, ts, 100,10)
stop(step) = step.s <= ϵ
#lineplot= scene[end]
display(scene)
framerate = 5
record(scene, "emd_sifting_$(framerate).mp4", framerate=framerate) do io
    for (i, step) in enumerate(EMD.halt(SiftIterable(x, ts, 4), stop))
        @show sum(abs, step.yvec)
        @show i
        #t[] = min(3,i)

        #@show step
        imf= step.yvec
        yvec[] = imf
        x_max.val = Int[]
        y_max[] = Int[]
        x_min.val = Int[]
        y_min[] = Int[]
        maxts[] = fill(NaN, length(ts))
        mints[] = fill(NaN, length(ts))
        sleep(1/10)
        recordframe!(io)
        #@show step.maxes


        if i>1
        smin = EMD.get_edgepoint(step.yvec, step.xvec, step.mins, first, isless)
        smax = EMD.get_edgepoint(step.yvec, step.xvec, step.maxes, first, !isless)
        emin = EMD.get_edgepoint(step.yvec, step.xvec, step.mins, last, isless)
        emax = EMD.get_edgepoint(step.yvec, step.xvec, step.maxes, last, !isless)
        x_max.val = ts[step.maxes]
        y_max[] = imf[step.maxes]
        sleep(1/10)
        recordframe!(io)
        step!(st)

        maxts[] = EMD.interpolate([first(step.xvec); step.xvec[step.maxes]; last(step.xvec)],[smax; step.yvec[step.maxes]; emax],step.xvec,EMD.DataInterp())

        sleep(1/10)
        recordframe!(io)
        step!(st)

        x_min.val = ts[step.mins]
        y_min[] = imf[step.mins]
        sleep(1/10)
        recordframe!(io)
        step!(st)
        mints[] = EMD.interpolate([first(step.xvec); step.xvec[step.mins]; last(step.xvec)],[smin; step.yvec[step.mins]; emin],step.xvec,EMD.DataInterp())

        sleep(1/10)
        recordframe!(io)
        step!(st)

        #lines!(scene, ts, maxTS)
        end
        sleep(1/10)
        recordframe!(io)
        step!(st)

    end
end

using AbstractPlotting.MakieLayout
imfs = EMD.emd(x,ts)
sc, layout = layoutscene(outer_padding, resolution=(1200,700), backgroundcolor= RGBf0(0.98,0.98,0.98))
display(sc)
ax1 = layout[1,1] = LAxis(sc, title="Original")
ax2 = layout[1,2] = LAxis(sc, title="IMF 1")
ax3 = layout[1,3] = LAxis(sc, title="IMF 2")
ax4 = layout[2,1] = LAxis(sc, title="IMF 3")
ax5 = layout[2,2] = LAxis(sc, title="Residual")
#ax6 = layout[2,3] = LAxis(sc, title="Data Components")

lines!(ax1, ts, x)
lines!(ax2, ts, imfs[1])
lines!(ax3, ts, imfs[2])
lines!(ax4, ts, imfs[3])
lines!(ax5, ts, imfs[4])
#lines!(ax6, ts, data[1])
#lines!(ax6, ts, data[2])
#lines!(ax6, ts, data[3])
#lines!(ax6, ts, data[4])
hidexdecorations!(ax1, grid=false)
hidexdecorations!(ax2, grid=false)
hidexdecorations!(ax3, grid=false)
trim!(layout)

linkxaxes!(ax1, ax4, ax2, ax3, ax5)
linkxaxes!(ax2, ax5)
linkxaxes!(ax3, ax6)

imfpath = "/home/crem_fe/Documents/presentations/MPI_20201125/imf_testdata.png"
Makie.save(imfpath, sc)
using AbstractPlotting
