#include <iostream>
#include <iterator>
#include <string>
#include "image/image.hpp"
#include "boost/program_options.hpp"
#include "dsi_interface_static_link.h"
#include "mapping/fa_template.hpp"
extern fa_template fa_template_imp;
namespace po = boost::program_options;


void init_mask(const image::basic_image<unsigned char, 3>& image,
               image::basic_image<unsigned char, 3>& mask);

/**
 perform reconstruction
 */
int rec(int ac, char *av[])
{
    po::options_description rec_desc("reconstruction options");
    rec_desc.add_options()
    ("help", "help message")
    ("action", po::value<std::string>(), "rec:diffusion reconstruction trk:fiber tracking")
    ("source", po::value<std::string>(), "assign the .src file name")
    ("mask", po::value<std::string>(), "assign the mask file")
    ("method", po::value<int>(), "reconstruction methods (0:dsi, 1:dti, 2:qbi_frt, 3:qbi_sh, 4:gqi)")
    ("odf_order", po::value<int>(), "set odf dimensions (4:162 direcitons, 5:252 directions, 6:362 directions, 8:642 directions)")
    ("record_odf", po::value<int>()->default_value(0), "set record odf to on using --record_odf=1")
    ("thread", po::value<int>()->default_value(2), "set the multi-thread count --thread=2")
    ("num_fiber", po::value<int>(), "maximum fibers resolved per voxel, default=3")
    ("decomposition", po::value<int>(), "set the decomposition")
    ("param0", po::value<float>(), "set parameters")
    ("param1", po::value<float>(), "set parameters")
    ("param2", po::value<float>(), "set parameters")
    ("param3", po::value<float>(), "set parameters")
    ;

    if(!ac)
    {
        std::cout << rec_desc << std::endl;
        return 1;
    }
    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).options(rec_desc).allow_unregistered().run(), vm);
    po::notify(vm);

    std::ofstream out("log.txt");

    std::string file_name = vm["source"].as<std::string>();
    ImageModel *handle = 0;
    out << "loading..." <<std::endl;
    if (!(handle = (ImageModel *)init_reconstruction(file_name.c_str())))
    {
        out << "Invalid src file format. " << std::endl;
        return 1;
    }
    out << "src read" <<std::endl;

    float param[4] = {0,0,0,0};
    int method_index = 0;
    int odf_order = 8; // for 362 directions
    int num_fiber = 3;
    if (vm.count("odf_order"))
        odf_order = vm["odf_order"].as<int>();
    if (vm.count("num_fiber"))
        num_fiber = vm["num_fiber"].as<int>();


    method_index = vm["method"].as<int>();
    if(method_index == 0)
        param[0] = 17.0;
    if(method_index == 2)
    {
        param[0] = 5;
        param[1] = 15;
    }
    if(method_index == 3)
        param[0] = 0.006;
    if(method_index == 4)
        param[0] = 1.2;
    if(method_index == 7)
    {
        if(!fa_template_imp.load_from_file("FMRIB58_FA_1mm.nii"))
        {
            out << "Cannot find FA template" << std::endl;
            free_reconstruction(handle);
            return -1;
        }
        param[0] = 1.2;
        param[1] = 2.0;
    }
    param[3] = 0.0002;

    if (vm.count("param0"))
        param[0] = vm["param0"].as<float>();
    if (vm.count("param1"))
        param[1] = vm["param1"].as<float>();
    if (vm.count("param2"))
        param[2] = vm["param2"].as<float>();
    if (vm.count("param3"))
        param[3] = vm["param3"].as<float>();

    unsigned int options[10];
    options[0] = method_index;
    options[1] = vm["thread"].as<int>(); //thread count
    options[2] = odf_order;
    options[3] = vm["record_odf"].as<int>();
    options[4] = 0;//vm.count("gfa_mask") ? 1 : 0;
    options[5] = 0;//vm.count("deconvolution") ? 1:0;
    options[6] = vm.count("decomposition") ? 1:0;
    options[7] = vm.count("half_sphere") ? 1:0;
    options[8] = 0;//vm.count("recursive_sharpening") ? 1:0;
    options[9] = num_fiber;//# of resolved fibers
    {
        out << "method=" << method_index << std::endl;
        out << "odf_order=" << odf_order << std::endl;
        out << "num_fiber=" << num_fiber << std::endl;
        out << "record_odf=" << vm.count("record_odf") << std::endl;
    }

    {
        const unsigned short* dim = get_dimension(handle);
        image::geometry<3> geometry(dim[0],dim[1],dim[2]);
        if(vm.count("mask"))
        {
            std::string mask_file = vm["mask"].as<std::string>();
            out << "reading mask..." << mask_file << std::endl;
            image::io::nifti header;
            if(!header.load_from_file(mask_file.c_str()))
            {
                out << "fail reading the mask...using default mask" << std::endl;
                goto no_mask;
            }
            image::basic_image<unsigned char,3> mask;
            header >> mask;
            if(header.nif_header.srow_x[0] < 0 || !header.is_nii)
                image::flip_y(mask);
            else
                image::flip_xy(mask);

            if(mask.geometry() != geometry)
            {
                out << "In consistent the mask dimension...using default mask" << std::endl;
                goto no_mask;
            }
            std::copy(mask.begin(),mask.end(),get_mask_image(handle));
        }
        else
            if(method_index != 7) //QSDR does not need mask
        {
            no_mask:
            out << "init mask..." <<std::endl;
            image::basic_image<unsigned char,3> image(geometry);
            image::basic_image<unsigned char,3> mask(geometry);
            std::copy(get_mask_image(handle),get_mask_image(handle)+image.size(),image.begin());
            init_mask(image,mask);
            std::copy(mask.begin(),mask.end(),get_mask_image(handle));
        }
        }
    out << "start reconstruction..." <<std::endl;
    const char* msg = reconstruction(handle,options,param);
    if (!msg)
        out << "Reconstruction finished:" << msg << std::endl;
    free_reconstruction(handle);
    return 0;
}
