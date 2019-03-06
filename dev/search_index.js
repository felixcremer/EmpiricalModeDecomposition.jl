var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#EmpiricalModeDecomposition.jl-1",
    "page": "Home",
    "title": "EmpiricalModeDecomposition.jl",
    "category": "section",
    "text": "This package implements the Empirical Mode Decomposition (EMD), which was developed by [Huang et al., 1998], and derivatives thereof.To install the package, use] install https://github.com/felixcremer/EmpiricalModeDecomposition.jl.git"
},

{
    "location": "emd/#",
    "page": "Empirical Mode Decomposition",
    "title": "Empirical Mode Decomposition",
    "category": "page",
    "text": ""
},

{
    "location": "emd/#Empirical-Mode-Decomposition-1",
    "page": "Empirical Mode Decomposition",
    "title": "Empirical Mode Decomposition",
    "category": "section",
    "text": "The empirical mode decomposition was developed by Huang in [Huang 1988]. It decomposes the signal into several Intrinsic Mode Functions (IMF). But unlike the Fourier or wavelet transformation the IMFs are computed adaptively from the data. The algorithms works as follows:Detect the local extrema of the original dataset x(t)\nCompute a spline s_mins_max through the local minima and maxima respectively\nCompute the mean of s_min and s_max and subtract it from x(t).  \nmath  h_{1,k}(t)=x(t)-\\frac{s_{min}-s_{max}}{2}  \nRepeat the previous steps with h_1k as x(t) until h_1k suffice a stopping criteria. h_1k is the first IMF c_1(t) of x(t)\nRepeat the steps 1 till 4 with the difference between x(t) and c_1(t) until there is no IMF to extractThis procedure decomposes the dataset x(t) into N IMFs and one residual:x(t)=sum_i=1^N c_i(t)+r(t)"
},

{
    "location": "api/#",
    "page": "Library",
    "title": "Library",
    "category": "page",
    "text": ""
},

{
    "location": "api/#EmpiricalModeDecomposition.ceemd-Tuple{Any,Any}",
    "page": "Library",
    "title": "EmpiricalModeDecomposition.ceemd",
    "category": "method",
    "text": "ceemd(measurements, xvec; num_imfs=6)\n\nCompute the Complete Empirical Mode Decomposition of the time series with values measurements and time steps xvec. numimfs is the Number of Intrinsic Mode Functions. Returns a list of numimfs + 1 Vectors of the same size as measurements.\n\n\n\n\n\n"
},

{
    "location": "api/#EmpiricalModeDecomposition.eemd",
    "page": "Library",
    "title": "EmpiricalModeDecomposition.eemd",
    "category": "function",
    "text": "eemd(measurements, xvec, numtrails=100)\n\nReturn the Intrinsic Mode Functions and the residual of the ensemble Empirical Mode Decomposition of the measurements given on time steps xvec.\n\n\n\n\n\n"
},

{
    "location": "api/#EmpiricalModeDecomposition.emd",
    "page": "Library",
    "title": "EmpiricalModeDecomposition.emd",
    "category": "function",
    "text": "emd(measurements, xvec)\n\nReturn the Intrinsic Mode Functions and the residual of the Empirical Mode Decomposition of the measurements given on time steps given in xvec.\n\n\n\n\n\n"
},

{
    "location": "api/#EmpiricalModeDecomposition.maketestdata-Tuple{Any}",
    "page": "Library",
    "title": "EmpiricalModeDecomposition.maketestdata",
    "category": "method",
    "text": "maketestdata(seed)\n\nReturn a simple example time series with the composing parts\n\n\n\n\n\n"
},

{
    "location": "api/#EmpiricalModeDecomposition.EMDIterable",
    "page": "Library",
    "title": "EmpiricalModeDecomposition.EMDIterable",
    "category": "type",
    "text": "Iterator for the Empirical Mode Decomposition. The time series values are an AbstractVector of type T, and the time positions are an AbstractVector of type U.\n\n\n\n\n\n"
},

{
    "location": "api/#EmpiricalModeDecomposition.SiftIterable",
    "page": "Library",
    "title": "EmpiricalModeDecomposition.SiftIterable",
    "category": "type",
    "text": "SiftIterable{T, U}\n\nIterator for the sifting algorithm. The time series values are an AbstractVector of type T, and the time positions are an AbstractVector of type U. Fields:\n\n\n\n\n\n"
},

{
    "location": "api/#Public-API-1",
    "page": "Library",
    "title": "Public API",
    "category": "section",
    "text": "DocTestSetup= quote\nusing EmpiricalModeDecomposition\nendModules = [EmpiricalModeDecomposition]\nPrivate = false\nOrder = [:function, :type]"
},

{
    "location": "api/#EmpiricalModeDecomposition.get_edgepoint-NTuple{5,Any}",
    "page": "Library",
    "title": "EmpiricalModeDecomposition.get_edgepoint",
    "category": "method",
    "text": "get_edgepoint(y, xvec, extremas, pos, comp)\n\nCompute the edgepoint which should be used as the extrema on the edge for the spline computation.\n\n\n\n\n\n"
},

{
    "location": "api/#EmpiricalModeDecomposition.ismonotonic-Union{Tuple{AbstractArray{T,1}}, Tuple{T}} where T",
    "page": "Library",
    "title": "EmpiricalModeDecomposition.ismonotonic",
    "category": "method",
    "text": "ismonotonic(x::AbstractVector)\n\nCheck wether x is monotonic. This means, every value is either larger or smaller than the preceding value.\n\n\n\n\n\n"
},

{
    "location": "api/#EmpiricalModeDecomposition.localmaxmin!-Tuple{Any,Array{Int64,1},Array{Int64,1}}",
    "page": "Library",
    "title": "EmpiricalModeDecomposition.localmaxmin!",
    "category": "method",
    "text": "localmaxmin!(x, maxes, mins)\n\nDetect the local extrema of x. Push the maxima into maxes and the minima into mins.\n\n\n\n\n\n"
},

{
    "location": "api/#EmpiricalModeDecomposition.sift",
    "page": "Library",
    "title": "EmpiricalModeDecomposition.sift",
    "category": "function",
    "text": "sift(y, xvec)\n\nSift the vector y whose points have x coordinates given by xvec.\n\n\n\n\n\n"
},

{
    "location": "api/#EmpiricalModeDecomposition.zerocrossing!-Tuple{Any,Any}",
    "page": "Library",
    "title": "EmpiricalModeDecomposition.zerocrossing!",
    "category": "method",
    "text": "zerocrossing!(y, crosses)\n\nCompute the indices of zerocrossings of a vector It searches for elements which are either zero or near a signflip and pushes the indices into crosses.\n\n\n\n\n\n"
},

{
    "location": "api/#EmpiricalModeDecomposition.CEEMDIterable",
    "page": "Library",
    "title": "EmpiricalModeDecomposition.CEEMDIterable",
    "category": "type",
    "text": "Iterator for the Complete Empirical Mode Decomposition. The time series values are an AbstractVector of type T, and the time positions are an AbstractVector of type U.\n\n\n\n\n\n"
},

{
    "location": "api/#EmpiricalModeDecomposition.CEEMDState",
    "page": "Library",
    "title": "EmpiricalModeDecomposition.CEEMDState",
    "category": "type",
    "text": "Intermediate results of the CEEMD Iteration\n\n\n\n\n\n"
},

{
    "location": "api/#EmpiricalModeDecomposition.SiftState",
    "page": "Library",
    "title": "EmpiricalModeDecomposition.SiftState",
    "category": "type",
    "text": "SiftState\n\nHandle the  intermediate results of the sifting.\n\n\n\n\n\n"
},

{
    "location": "api/#Internal-API-1",
    "page": "Library",
    "title": "Internal API",
    "category": "section",
    "text": "Modules = [EmpiricalModeDecomposition]\nPublic = false\nOrder = [:function, :type]"
},

]}
