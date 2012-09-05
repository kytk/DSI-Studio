#include "stdafx.h"
#include "odf_fem_interface_static_link.hpp"
#include "prog_interface_static_link.h"
#include "mat_file.hpp"
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/math/distributions/students_t.hpp>



#include "basic_voxel.hpp"
#include "image_model.hpp"
#include "odf_decomposition.hpp"
#include "odf_process.hpp"

#include "sample_model.hpp"
#include "space_mapping.hpp"


#include "dti_process.hpp"
#include "dsi_process.hpp"
#include "qbi_process.hpp"
#include "sh_process.hpp"
#include "gqi_process.hpp"
#include "gqi_mni_reconstruction.hpp"

#include "odf_deconvolusion.hpp"

#include "dti_voxel.hpp"
#include "dsi_voxel.hpp"
#include "qbi_voxel.hpp"
#include "gqi_voxel.hpp"
#include "image_model.hpp"



#include "odf_decomposition.hpp"
#include "mix_gaussian_model.hpp"
#include "racian_noise.hpp"

#include "layout.hpp"


boost::mt19937 RacianNoise::generator(static_cast<unsigned> (std::time(0)));
boost::normal_distribution<float> RacianNoise::normal;
boost::uniform_real<float> RacianNoise::uniform(0.0,1.0);
boost::variate_generator<boost::mt19937&,
boost::normal_distribution<float> > RacianNoise::gen_normal(RacianNoise::generator,RacianNoise::normal);
boost::variate_generator<boost::mt19937&,
boost::uniform_real<float> > RacianNoise::gen_uniform(RacianNoise::generator,RacianNoise::uniform);

const char* process_position = 0;


std::string error_msg;


extern "C"
    void* init_reconstruction(const char* file_name)
{
    std::auto_ptr<ImageModel> image(new ImageModel);
    if (!image->load_from_file(file_name))
        return 0;
    return image.release();
}
extern "C"
    void free_reconstruction(ImageModel* image_model)
{
    delete image_model;
}


extern "C"
    const float* get_b_table(ImageModel* image_model,unsigned int& b_number)
{
    unsigned int row;
    const float* table;
    image_model->mat_reader->get_matrix("b_table",row,b_number,table);
    return table;
}


extern "C"
    const unsigned short* get_dimension(ImageModel* image_model)
{
    static unsigned short dim[3];
    dim[0] = image_model->voxel.matrix_width;
    dim[1] = image_model->voxel.matrix_height;
    dim[2] = image_model->voxel.slice_number;
    return dim;
}

extern "C"
    const float* get_voxel_size(ImageModel* image_model)
{
    return image_model->voxel.voxel_size;
}

extern "C"
    unsigned char* get_mask_image(ImageModel* image_model)
{
    return &*image_model->mask.begin();
}

extern "C"
    char* check_reconstruction(ImageModel* image_model)
{
    static char ava[13];
    ava[0] = image_model->avaliable<CheckDSI>();
    ava[1] = image_model->avaliable<CheckDTI>();
    ava[2] = ava[3] = image_model->avaliable<CheckHARDI>();
    ava[4] = 1;
    ava[5] = 1;
    ava[6] = 1;
    ava[7] = 1;
    return ava;
}



extern "C"
    const char* reconstruction(ImageModel* image_model,unsigned int* options,const float* param_values)
{
    unsigned int method = options[0];
    unsigned int threshold_count = options[1];
    unsigned int odf_fold = options[2];
    unsigned int need_odf = options[3];
    unsigned int gfa_as_threshold = options[4];
    unsigned int de_odf = options[5];
    unsigned int decompose = options[6];
    unsigned int half_sphere = options[7];
    unsigned int odf_deconvolusion_iterative = options[8];
    unsigned int max_fiber_number = options[9];
    static std::string output_name;
    try
    {
        process_position = __FUNCTION__;
        ti_initialize(odf_fold);
        image_model->thread_count = threshold_count;
        image_model->voxel.param = param_values;
        image_model->voxel.full_odf_dimension = ti_vertices_count();
        image_model->voxel.need_odf = need_odf;
        image_model->voxel.gfa_as_threshold = gfa_as_threshold;
        image_model->voxel.odf_deconvolusion = de_odf;
        image_model->voxel.odf_decomposition = decompose;
        image_model->voxel.half_sphere = half_sphere;
        image_model->voxel.odf_deconvolusion_iterative = odf_deconvolusion_iterative;
        image_model->voxel.max_fiber_number = max_fiber_number;
        std::ostringstream out;
        out << ".odf" << odf_fold;
        out << ".f" << max_fiber_number;
        if (need_odf)
            out << "rec";
        if (gfa_as_threshold)
            out << ".gfa_mask";
        if (half_sphere)
            out << ".hs";
        if (de_odf)
            out << ".de" << param_values[2];
        if (odf_deconvolusion_iterative)
            out << ".iter" << (int) std::floor(param_values[3]+0.5);
        if (decompose)
            out << ".dec";
        switch (method)
        {
        case 0: //DSI local max
            if (de_odf)
            {
                if (!image_model->reconstruct<dsi_estimate_response_function>())
                    return "reconstruction calceled";
                begin_prog("preparing deconvolution");
            }
            out << ".dsi."<< (int)param_values[0] << ".fib";
            if (!image_model->reconstruct<dsi_process>(out.str()))
                return "reconstruction canceled";
            break;
        case 1://DTI
            out << ".dti.fib";
            image_model->voxel.max_fiber_number = 1;
            if (!image_model->reconstruct<dti_transformation_process>(out.str()))
                return "reconstruction canceled";
            break;

        case 2://QBI
            if (de_odf)
            {
                if (!image_model->reconstruct<qbi_estimate_response_function>())
                    return "reconstruction calceled";
                begin_prog("preparing deconvolution");
            }
            out << ".qbi."<< param_values[0] << "_" << param_values[1] << ".fib";
            if (!image_model->reconstruct<qbi_process>(out.str()))
                return "reconstruction canceled";
            break;
        case 3://QBI
            if (de_odf)
            {
                if (!image_model->reconstruct<qbi_sh_estimate_response_function>())
                    return "reconstruction calceled";

                begin_prog("preparing deconvolution");
            }
            out << ".qbi.sh."<< param_values[0] << ".fib";
            if (!image_model->reconstruct<qbi_sh_process>(out.str()))
                return "reconstruction canceled";
            break;

        case 4://GQI
            if (de_odf)
            {
                if (!image_model->reconstruct<gqi_estimate_response_function>())
                    return "reconstruction calceled";
                begin_prog("preparing deconvolution");
            }
            out << ".gqi."<< param_values[0] << ".fib";
            if (!image_model->reconstruct<gqi_process>(out.str()))
                return "reconstruction canceled";
            break;

        case 5://GQI
            if (de_odf)
            {
                if (!image_model->reconstruct<gqi_r2_estimate_response_function>())
                    return "reconstruction calceled";
                begin_prog("preparing deconvolution");
            }
            out << ".gqi.r2."<< param_values[0] << ".fib";
            if (!image_model->reconstruct<gqi_r2_process>(out.str()))
                return "reconstruction canceled";
            break;
        case 6:
            out << ".gqi.hardi"<< param_values[0] << ".src";
            if (!image_model->reconstruct<gqi_adaptor_process>(out.str()))
                return "reconstruction canceled";
            break;
        case 7:
            // run gqi to get the spin quantity
            std::fill(image_model->mask.begin(),image_model->mask.end(),1.0);
            if (!image_model->reconstruct<gqi_estimate_response_function>())
                return "reconstruction calceled";
            out << ".mni.gqi.";
            if(param_values[2] > 0)
                out << ".r2.";
            out << param_values[0] << "vs" << param_values[1] << ".fib";
            begin_prog("preparing deformation");
            if (!image_model->reconstruct<gqi_mni_process>(out.str()))
                return "reconstruction canceled";
            break;
        }
        output_name = image_model->file_name + out.str() + ".gz";
    }
    catch (std::exception& exp)
    {
        error_msg = "Reconstruction error. Please contact the program developer with the following information:";
        error_msg += exp.what();
        error_msg += " happened";
        if (process_position)
        {
            error_msg += " in ";
            error_msg += process_position;
        }
        return error_msg.c_str();
    }
    catch (...)
    {
        error_msg = "Reconstruction error. Please contact the program developer with the following information:";
        error_msg += "unexpected exception happened in ";
        if (process_position)
        {
            error_msg += " in ";
            error_msg += process_position;
        }
        return error_msg.c_str();
    }
    return output_name.c_str();
}


bool output_odfs(const image::basic_image<unsigned char,3>& mni_mask,
                 const char* out_name,
                 const char* ext,
                 const std::vector<std::vector<float> >& odfs,
                 bool record_odf = true)
{
    begin_prog("output");
    ImageModel image_model;
    image_model.set_dimension(mni_width,mni_height,mni_depth);
    image_model.voxel.full_odf_dimension = ti_vertices_count();
    image_model.voxel.q_count = 0;
    image_model.voxel.need_odf = true;
    image_model.voxel.odf_decomposition = false;
    image_model.voxel.odf_deconvolusion = false;
    image_model.voxel.gfa_as_threshold = false;
    image_model.voxel.half_sphere = false;
    image_model.voxel.max_fiber_number = 3;
    image_model.voxel.qa_scaling = 1;
    image_model.thread_count = 1;
    std::copy(mni_mask.begin(),mni_mask.end(),image_model.mask.begin());
    std::fill(image_model.voxel.voxel_size,image_model.voxel.voxel_size+3,2.0);
    image_model.file_name = out_name;
    image_model.voxel.param = (const float*)&odfs;
    image_model.voxel.need_odf = record_odf;
    if (prog_aborted() || !image_model.reconstruct<reprocess_odf>(ext))
        return false;
    return true;
}


extern "C"
    bool odf_average(const char* out_name,
                     const char* const * file_names,
                     unsigned int num_files)
{
    image::basic_image<unsigned char,3> mask;
    std::vector<std::vector<float> > odfs;
    begin_prog("averaging");
    can_cancel(true);
    unsigned int half_vertex_count = 0;
    for (unsigned int index = 0;check_prog(index,num_files);++index)
    {
        const char* file_name = file_names[index];
        MatFile reader;
        if(!reader.load_from_file(file_name))
        {
            std::cout << "Cannot open file" << std::endl;
            std::cout << file_name << std::endl;
            return false;
        }
        if(index == 0)
        {
            unsigned int row,col;
            const float* odf_buffer;
            const float* fa0;
            const short* face_buffer;
            const unsigned short* dimension;
            reader.get_matrix("dimension",row,col,dimension);
            mask.resize(image::geometry<3>(dimension));
            reader.get_matrix("fa0",row,col,fa0);
            for(unsigned int index = 0;index < mask.size();++index)
                mask[index] = fa0[index] == 0 ? 0:1;

            reader.get_matrix("odf_vertices",row,col,odf_buffer);
            if (!odf_buffer)
            {
                std::cout << "Cannot find odf vertices" << std::endl;
                std::cout << file_name << std::endl;
                return false;
            }
            half_vertex_count = col >> 1;
            ti_load_odf_vertices(col,odf_buffer);
            reader.get_matrix("odf_faces",row,col,face_buffer);
            if (!face_buffer)
            {
                std::cout << "Cannot find odf faces" << std::endl;
                std::cout << file_name << std::endl;
                return false;
            }
            ti_load_odf_faces(col,face_buffer);
        }
        else
        // check odf consistency
        {
            unsigned int row,col;
            const float* odf_buffer;
            reader.get_matrix("odf_vertices",row,col,odf_buffer);
            if (!odf_buffer)
                return false;
            if(col != ti_vertices_count())
                return false;
            for (unsigned int index = 0;index < col;++index,odf_buffer += 3)
            {
                if(ti_vertices(index)[0] != odf_buffer[0] ||
                   ti_vertices(index)[1] != odf_buffer[1] ||
                   ti_vertices(index)[2] != odf_buffer[2])
                {
                    std::cout << "Inconsistent ODF" << std::endl;
                    std::cout << file_name << std::endl;
                    return false;
                }
            }
        }
        std::vector<const float*> odf_bufs;
        std::vector<unsigned int> odf_bufs_size;
        //get_odf_bufs(reader,odf_bufs,odf_bufs_size);
        {
            odf_bufs.clear();
            odf_bufs_size.clear();
            for (unsigned int odf_index = 0;1;++odf_index)
            {
                unsigned int row,col;
                std::ostringstream out;
                out << "odf" << odf_index;
                const float* odf_buf = 0;
                reader.get_matrix(out.str().c_str(),row,col,odf_buf);
                if (!odf_buf)
                    break;
                odf_bufs.push_back(odf_buf);
                odf_bufs_size.push_back(row*col);
            }
        }
        odfs.resize(odf_bufs.size());
        for(unsigned int i = 0;i < odf_bufs.size();++i)
        {
            odfs[i].resize(odf_bufs_size[i]);
            image::add(odfs[i].begin(),odfs[i].end(),odf_bufs[i]);
            for(unsigned int j = 0;j < odf_bufs_size[i];)
            {
                unsigned int next_j = j + half_vertex_count;
                image::minus_constant(odfs[i].begin() + j,
                                      odfs[i].begin() + next_j,
                    *std::min_element(odf_bufs[i]+j,odf_bufs[i]+next_j));
                j = next_j;
            }
        }
    }
    if (prog_aborted())
        return false;

    for (unsigned int odf_index = 0;odf_index < odfs.size();++odf_index)
        for (unsigned int j = 0;j < odfs[odf_index].size();++j)
            odfs[odf_index][j] /= (double)num_files;

    output_odfs(mask,out_name,".mean.odf.fib",odfs);
    output_odfs(mask,out_name,".mean.fib",odfs,false);
}

extern "C"
    bool generate_simulation(
        const char* bvec_file_name,unsigned char s0_snr,float mean_dif,unsigned char odf_fold,
        const char* fa_iteration,
        const char* crossing_angle_iteration,
        unsigned char repeat_num)
{
    ti_initialize(odf_fold);
    Layout layout(s0_snr,mean_dif);
    if (!layout.load_b_table(bvec_file_name))
        return false;
    std::vector<float> fa;
    std::vector<float> angle;
    {
        std::string fa_iteration_str(fa_iteration);
        std::istringstream tmp(fa_iteration_str);
        std::copy(std::istream_iterator<float>(tmp),
                  std::istream_iterator<float>(),std::back_inserter(fa));
    }
    {
        std::string crossing_angle_iteration_str(crossing_angle_iteration);
        std::istringstream tmp(crossing_angle_iteration_str);
        std::copy(std::istream_iterator<float>(tmp),
                  std::istream_iterator<float>(),std::back_inserter(angle));
    }

    layout.createLayout(fa,angle,repeat_num);
    std::ostringstream out;
    out << bvec_file_name << "_snr" << (int)s0_snr << "_dif" << mean_dif << "_odf" << (int)odf_fold << "_n" << (int)repeat_num << ".src";
    layout.generate(out.str().c_str());
    return true;
}








