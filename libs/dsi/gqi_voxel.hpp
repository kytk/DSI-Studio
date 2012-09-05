#ifndef DDI_VOXEL_H
#define DDI_VOXEL_H
#include "odf_process.hpp"
typedef boost::mpl::vector<
ReadDWIData,
CorrectB0,
QSpace2Odf,
QAScaling,
ODFDeconvolusion,
ODFDecomposition,
DetermineFiberDirections,
OutputODF,
SaveFA,
SaveDirIndex
> gqi_process;


typedef boost::mpl::vector<
ReadDWIData,
CorrectB0,
QSpace2Odf,
QAScaling,
DetermineFiberDirections,
AnaQSpace2Odf,
OutputODF,
SaveFA,
SaveDir
> gqi22_process;


// for ODF deconvolution
typedef boost::mpl::vector<
ReadDWIData,
CorrectB0,
QSpace2Odf,
DetermineFiberDirections,
RecordQA,
EstimateResponseFunction
> gqi_estimate_response_function;


typedef boost::mpl::vector<
ReadDWIData,
CorrectB0,
QSpace2OdfR2,
ODFDeconvolusion,
ODFDecomposition,
DetermineFiberDirections,
QAScaling,
OutputODF,
SaveFA,
SaveDirIndex
> gqi_r2_process;

// for ODF deconvolution
typedef boost::mpl::vector<
ReadDWIData,
CorrectB0,
QSpace2OdfR2,
EstimateResponseFunction
> gqi_r2_estimate_response_function;


typedef boost::mpl::vector<
ReadDWIData,
QSpaceAdoptor<40>
> gqi_adaptor_process;

typedef boost::mpl::vector<
                        GQI_MNI,
                        //GQI_phantom,
                        ODFDecomposition,
                        DetermineFiberDirections,
                        OutputODF,
                        SaveFA,
                        SaveDirIndex
				> gqi_mni_process;

typedef boost::mpl::vector<
ReadDWIData,
ODFLoader,
DetermineFiberDirections,
SaveFA,
SaveDirIndex
> reprocess_odf;

#endif//DDI_VOXEL_H
