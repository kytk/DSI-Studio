#ifndef QBI_VOXEL_H
#define QBI_VOXEL_H

const unsigned int equator_sample_count = 40;
typedef boost::mpl::vector<
ReadDWIData,
CorrectB0,
QBIReconstruction<equator_sample_count>,
ODFDeconvolusion,
ODFDecomposition,
DetermineFiberDirections,
QAScaling,
OutputODF,
SaveFA,
SaveDirIndex
> qbi_process;



// for ODF deconvolution
typedef boost::mpl::vector<
ReadDWIData,
CorrectB0,
QBIReconstruction<equator_sample_count>,
DetermineFiberDirections,
EstimateResponseFunction
> qbi_estimate_response_function;

// SH
typedef boost::mpl::vector<
ReadDWIData,
CorrectB0,
SHDecomposition<8>,
ODFDeconvolusion,
ODFDecomposition,
DetermineFiberDirections,
QAScaling,
OutputODF,
SaveFA,
SaveDirIndex
> qbi_sh_process;


// for ODF deconvolution
typedef boost::mpl::vector<
ReadDWIData,
CorrectB0,
SHDecomposition<8>,
DetermineFiberDirections,
EstimateResponseFunction
> qbi_sh_estimate_response_function;

#endif//QBI_VOXEL_H
