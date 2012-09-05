#ifndef DTI_VOXEL_H
#define DTI_VOXEL_H

typedef boost::mpl::vector<
ReadDWIData,
ADCProfile,
Dwi2Tensor,
TensorEigenAnalysis
//OutputODF
> dti_transformation_process;

#endif//DTI_VOXEL_H
