TT-MRF
======

TT-MRF is a library for Markov random field inference with Tensor Train decomposition approach.

Install
=======

Install the [TT-Toolbox](https://github.com/oseledets/TT-Toolbox) (just download it and run `setup.m` to add everything important to MATLAB path).

Run `setup.m` to add required folders to MATLAB path.

Optional dependencies
---------------------

If you want to compare with state-of-the-art techniques for computing partition function and marginal distributions, install [LibDAI](http://staff.science.uva.nl/~jmooij1/libDAI/) library with the MATLAB support.

If you need access to state-of-the-art MAP-inference routines OR want to be able to load models in OpenGM format, install [OpenGM](hci.iwr.uni-heidelberg.de/opengm2/) with MATLAB, HDF5, TRW-S and Maxflow support. Example bash input for compiling OpenGM:
``` bash
cd opengm_folder
cmake . -DWITH_MATLAB=ON -DWITH_HDF5=ON -DBUILD_MATLAB_WRAPPER=ON -DWITH_TRWS=ON -DWITH_MAXFLOW=ON
make
```


Model format
==============

We use custom graphical model instance format. You can build problems like this:
``` MATLAB
% Build 5x4 grid spin glass model with temperature = 2.
Model = generate_spin_glass_model(5, 4, 2);

% Load model in OpenGM format.
Model = load_opengm_model('matching/matching0.h5');

% Load model in uai format.
Model = load_uai_model('RBM.uai');
```

Format details:
``` MATLAB
Model.libdaiFactors     [Cell array] factors of the model in the LibDAI format
Model.numNodes          [Number] number of Model variables
Model.modeSizes         [Vector 1 x d] sizes of variables (e.g. x_1 is from {1, ..., modeSizes(1)})
Model.description       [String] text description
Model.type              [String] Type: 'Spin glass', 'OpenGM' or 'UAI'

% Problem specific, spin glass
Model.grid_n            [Number] vertical size of spin glass model grid
Model.grid_m            [Number] horizontal size of spin glass model grid
Model.temperature       [Number]
Model.unaryWeights      [Matrix n x m]
Model.unaryType         [String] 'number' if all unary weights equals to one number;
                            'matrix' matrix with unary weights was specified during model generation;
                            'rand' if weights were generated from uniform distribution
Model.unaryDistr        [Vector 1 x 2] unary wights uniform distribution support
                                          (e.g. [-1, 1] means that weights are from U(-1, 1))
Model.edgeWeights       [Vector numEdges x 1] all pairwise weights
Model.edgeType          [String] 'number' or 'rand', see details in unryType description
Model.edgesDistr        [Vector 1 x 2] pairwise weight uniform distribution support
```


Example code
==============
``` MATLAB
% Build 5x4 grid spin glass model with temperature = 2
% and pairwise weights generated from uniform distribution on [0, 1].
Model = generate_spin_glass_model(5, 4, 2, 'J', 'rand', 'J_distr', [0, 1]);

% Compute logarithm of the partition function using Tensor Train approach
% with rounding precision equals 1e-6.
logZ = compute_log_z(Model, 1e-6);

% Compute logarithm of the partition function using junction tree method from the libDAI library.
[logZ_JT, ~, ~] = dai_jtree(Model.libdaiFactors, { 1 }, '[updates=HUGIN]');

relError = (logZ_JT - logZ) / logZ_JT;
disp(['Computed logarithm of the partition function with relative error ', num2str(relError)])
```





