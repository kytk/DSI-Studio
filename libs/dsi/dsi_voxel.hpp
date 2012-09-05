#ifndef DSI_VOXEL_H
#define DSI_VOXEL_H


typedef boost::mpl::vector<
ReadDWIData,
QSpace2Pdf,
Pdf2Odf,
QAScaling,
ODFDeconvolusion,
ODFDecomposition,
DetermineFiberDirections,
OutputODF,
SaveFA,
SaveDirIndex
> dsi_process;

// for ODF deconvolution
typedef boost::mpl::vector<
ReadDWIData,
QSpace2Pdf,
Pdf2Odf,
DetermineFiberDirections,
EstimateResponseFunction
> dsi_estimate_response_function;


#endif
