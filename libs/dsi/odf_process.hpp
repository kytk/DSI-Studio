#ifndef ODF_TRANSFORMATION_PROCESS_HPP
#define ODF_TRANSFORMATION_PROCESS_HPP
#include <boost/lambda/lambda.hpp>
#include "basic_process.hpp"
#include "basic_voxel.hpp"

class ReadDWIData : public BaseProcess{
public:
    virtual void init(Voxel&) {}
    virtual void run(Voxel& voxel, VoxelData& voxel_data)
    {
        for (unsigned int index = 0; index < voxel.q_count; ++index)
            voxel_data.space[index] = voxel.image_model->dwi_data[index][voxel_data.voxel_index];
    }
    virtual void end(Voxel&,MatFile& mat_writer) {}
};
struct Accumulator
{
    float operator()(const std::vector<float>& odf)
    {
        float sum = std::accumulate(odf.begin(),odf.end(),0.0);
        sum /= ((float)odf.size());
        return sum;
    }
};


struct GeneralizedFA
{
    float operator()(const std::vector<float>& odf)
    {
        float m1 = 0.0;
        float m2 = 0.0;
        std::vector<float>::const_iterator iter = odf.begin();
        std::vector<float>::const_iterator end = odf.end();
        for (;iter != end; ++iter)
        {
            float t = *iter;
            m1 += t;
            m2 += t*t;
        }
        m1 *= m1;
        m1 /= ((float)odf.size());
        if (m2 == 0.0)
            return 0.0;
        return std::sqrt(((float)odf.size())/((float)odf.size()-1.0)*(m2-m1)/m2);
    }
};

struct RemoveIsotropicPart
{
    float operator()(std::vector<float>& odf)
    {
        float min_odf = *std::min_element(odf.begin(),odf.end());
        std::for_each(odf.begin(),odf.end(),boost::lambda::_1 -= min_odf);
        return min_odf;
    }
};

const unsigned int odf_block_size = 20000;
struct OutputODF : public BaseProcess
{
    std::vector<std::vector<float> > odf_data;
    std::vector<unsigned int> odf_index_map;
    unsigned int block_index;
    unsigned int odf_index;
public:
    virtual void init(Voxel& voxel)
    {
        process_position = __FUNCTION__;
        odf_data.clear();
        if (voxel.need_odf)
        {
            unsigned int total_count = 0;
            odf_index_map.resize(voxel.image_model->mask.size());
            for (unsigned int index = 0;index < voxel.image_model->mask.size();++index)
                if (voxel.image_model->mask[index])
                {
                    odf_index_map[index] = total_count;
                    ++total_count;
                }
            try
            {
                std::vector<unsigned int> size_list;
                while (1)
                {

                    if (total_count > odf_block_size)
                    {
                        size_list.push_back(odf_block_size);
                        total_count -= odf_block_size;
                    }
                    else
                    {
                        size_list.push_back(total_count);
                        break;
                    }
                }
                odf_data.resize(size_list.size());
                for (unsigned int index = 0;index < odf_data.size();++index)
                    odf_data[index].resize(size_list[index]*(voxel.full_odf_dimension >> 1));
            }
            catch (...)
            {
                odf_data.clear();
                voxel.need_odf = false;
            }

            block_index = 0;
            odf_index = 0;
        }

    }
    virtual void run(Voxel& voxel,VoxelData& data)
    {

        process_position = __FUNCTION__;
        if (voxel.need_odf)
        {
            unsigned int odf_index = odf_index_map[data.voxel_index];
            std::copy(data.odf.begin(),data.odf.end(),
                      odf_data[odf_index/odf_block_size].begin() + (odf_index%odf_block_size)*(voxel.full_odf_dimension >> 1));
        }

    }
    virtual void end(Voxel& voxel,MatFile& mat_writer)
    {

        process_position = __FUNCTION__;
        if (!voxel.need_odf)
            return;
        {
            set_title("output odfs");
            for (unsigned int index = 0;index < odf_data.size();++index)
            {
                if (!voxel.odf_deconvolusion)
                    std::for_each(odf_data[index].begin(),odf_data[index].end(),boost::lambda::_1 /= voxel.qa_scaling);
                std::ostringstream out;
                out << "odf" << index;
                mat_writer.add_matrix(out.str().c_str(),&*odf_data[index].begin(),
                                      voxel.full_odf_dimension >> 1,
                                      odf_data[index].size()/(voxel.full_odf_dimension >> 1));
            }
            odf_data.clear();
        }

    }
};

struct ODFLoader : public BaseProcess
{
    const float* odf_data;
    unsigned int cur_pos1;
        unsigned int cur_pos2;
	std::vector<std::vector<float> >* odfs;
public:
    virtual void init(Voxel& voxel)
    {
        process_position = __FUNCTION__;
        odfs = (std::vector<std::vector<float> >*)voxel.param;
        cur_pos1 = 0;
        cur_pos2 = 0;
    }
    virtual void run(Voxel& voxel, VoxelData& data)
    {
        process_position = __FUNCTION__;
        std::copy((*odfs)[cur_pos1].begin() + cur_pos2,(*odfs)[cur_pos1].begin() + cur_pos2 + data.odf.size(),data.odf.begin());
        cur_pos2 += data.odf.size();
		if(cur_pos2 >= (*odfs)[cur_pos1].size())
		{
			cur_pos2 = 0;
			++cur_pos1;
		}
    }
    virtual void end(Voxel& voxel,MatFile& mat_writer)
    {
        process_position = __FUNCTION__;
        if (voxel.need_odf)
        {
            set_title("output odfs");
            for (unsigned int index = 0;index < (*odfs).size();++index)
            {
                std::ostringstream out;
                out << "odf" << index;
                mat_writer.add_matrix(out.str().c_str(),&*(*odfs)[index].begin(),
                                      voxel.full_odf_dimension >> 1,
                                      (*odfs)[index].size()/(voxel.full_odf_dimension >> 1));
            }
        }
    }
};

struct QAScaling : public BaseProcess
{
protected:
    std::vector<float> max_odf;
    float iso_qa;
    image::vector<3,int> max_odf_pos;
    std::vector<float> iso;
protected:
    boost::mutex mutex;
public:
    virtual void init(Voxel& voxel)
    {
        process_position = __FUNCTION__;
        iso.clear();
        iso.resize(voxel.total_size);
        iso_qa = 0;
        voxel.qa_scaling = 1.0;
    }
    virtual void run(Voxel& voxel, VoxelData& data)
    {
        process_position = __FUNCTION__;
        if (voxel.odf_deconvolusion)
            return;
        iso[data.voxel_index] =
           *std::min_element(data.odf.begin(),data.odf.end());

        float sd = image::standard_deviation(data.odf.begin(),data.odf.end());
        // get the isotropic part and find the free water diffusion

        boost::mutex::scoped_lock lock(mutex);
        if (iso_qa < iso[data.voxel_index]-3.0*sd)
        {
            iso_qa = iso[data.voxel_index]-3.0*sd;
            max_odf = data.odf;
            int p_index = data.voxel_index;
            int x = p_index % voxel.matrix_width;
            p_index -= x;
            p_index /= voxel.matrix_width;
            int y = p_index % voxel.matrix_height;
            p_index -= y;
            p_index /= voxel.matrix_height;
            max_odf_pos[0] = x;
            max_odf_pos[1] = y;
            max_odf_pos[2] = p_index;
        }

    }
    virtual void end(Voxel& voxel,MatFile& mat_writer)
    {
        if (voxel.odf_deconvolusion)
            return;
        voxel.qa_scaling = std::accumulate(max_odf.begin(),max_odf.end(),0.0)
                            /(max_odf.size());

        // scaled to 1mm cubic
        if(voxel.voxel_size[0] == 0.0 ||
           voxel.voxel_size[1] == 0.0 ||
           voxel.voxel_size[2] == 0.0)
            throw std::runtime_error("No spatial information found in src file. Recreate src file or contact developer for assistance");
        voxel.qa_scaling /= voxel.voxel_size[0];
        voxel.qa_scaling /= voxel.voxel_size[1];
        voxel.qa_scaling /= voxel.voxel_size[2];

        mat_writer.add_matrix("qa_scaling",&voxel.qa_scaling,1,1);
        mat_writer.add_matrix("max_odf_pos",max_odf_pos.begin(),3,1);
        std::for_each(max_odf.begin(),max_odf.end(),boost::lambda::_1 /= voxel.qa_scaling);
        mat_writer.add_matrix("max_odf",&*max_odf.begin(),1,max_odf.size());
        mat_writer.add_matrix("iso",&*iso.begin(),1,iso.size());


    }
};


struct SaveFA : public BaseProcess
{
protected:
    std::vector<float> gfa;
    std::vector<std::vector<float> > fa;
public:
    virtual void init(Voxel& voxel)
    {
        process_position = __FUNCTION__;
        fa.resize(voxel.max_fiber_number);
        for (unsigned int index = 0;index < voxel.max_fiber_number;++index)
            fa[index].resize(voxel.total_size);
        gfa.clear();
        gfa.resize(voxel.total_size);
    }
    virtual void run(Voxel& voxel, VoxelData& data)
    {
        process_position = __FUNCTION__;
        gfa[data.voxel_index] = GeneralizedFA()(data.odf);
        for (unsigned int index = 0;index < voxel.max_fiber_number;++index)
            fa[index][data.voxel_index] = data.fa[index];
    }
    virtual void end(Voxel& voxel,MatFile& mat_writer)
    {
        set_title("output gfa");
        mat_writer.add_matrix("gfa",&*gfa.begin(),1,gfa.size());
        if (!voxel.odf_deconvolusion && voxel.qa_scaling != 0.0)
            for (unsigned int i = 0;i < voxel.max_fiber_number;++i)
                std::for_each(fa[i].begin(),fa[i].end(),boost::lambda::_1 /= voxel.qa_scaling);

        for (unsigned int index = 0;index < voxel.max_fiber_number;++index)
        {
            std::ostringstream out;
            out << index;
            std::string num = out.str();
            std::string fa_str = "fa";
            fa_str += num;
            set_title(fa_str.c_str());
            if (voxel.gfa_as_threshold)
                mat_writer.add_matrix(fa_str.c_str(),&*gfa.begin(),1,gfa.size());
            else
                mat_writer.add_matrix(fa_str.c_str(),&*fa[index].begin(),1,fa[index].size());
        }
    }
};



struct SaveDirIndex : public BaseProcess
{
protected:
    std::vector<std::vector<short> > findex;
public:
    virtual void init(Voxel& voxel)
    {
        process_position = __FUNCTION__;
        findex.resize(voxel.max_fiber_number);
        for (unsigned int index = 0;index < voxel.max_fiber_number;++index)
            findex[index].resize(voxel.total_size);
    }
    virtual void run(Voxel& voxel, VoxelData& data)
    {
        process_position = __FUNCTION__;
        for (unsigned int index = 0;index < voxel.max_fiber_number;++index)
            findex[index][data.voxel_index] = data.dir_index[index];
    }
    virtual void end(Voxel& voxel,MatFile& mat_writer)
    {
        for (unsigned int index = 0;index < voxel.max_fiber_number;++index)
        {
            std::ostringstream out;
            out << index;
            std::string num = out.str();
            std::string index_str = "index";
            index_str += num;
            set_title(index_str.c_str());
            mat_writer.add_matrix(index_str.c_str(),&*findex[index].begin(),1,findex[index].size());
        }
    }
};


struct SaveDir : public BaseProcess
{
protected:
    std::vector<std::vector<float> > dir;
public:
    virtual void init(Voxel& voxel)
    {
        process_position = __FUNCTION__;
        dir.resize(voxel.max_fiber_number);
        for (unsigned int index = 0;index < voxel.max_fiber_number;++index)
            dir[index].resize(voxel.total_size*3);
    }
    virtual void run(Voxel& voxel, VoxelData& data)
    {
        process_position = __FUNCTION__;
        unsigned int dir_index = data.voxel_index;
        for (unsigned int index = 0;index < voxel.max_fiber_number;++index)
            std::copy(data.dir[index].begin(),data.dir[index].end(),dir[index].begin() + dir_index + dir_index + dir_index);
    }
    virtual void end(Voxel& voxel,MatFile& mat_writer)
    {
        for (unsigned int index = 0;index < voxel.max_fiber_number;++index)
        {
            std::ostringstream out;
            out << index;
            std::string num = out.str();
            std::string index_str = "dir";
            index_str += num;
            set_title(index_str.c_str());
            mat_writer.add_matrix(index_str.c_str(),&*dir[index].begin(),1,dir[index].size());
        }
    }
};




struct SearchLocalMaximum
{
    std::vector<std::vector<unsigned short> > neighbor;
    std::map<float,unsigned short,std::greater<float> > max_table;
    void init(void)
    {
        process_position = __FUNCTION__;
        unsigned int half_odf_size = ti_vertices_count() >> 1;
        unsigned int faces_count = ti_faces_count();
        neighbor.resize(ti_vertices_count() >> 1);
        for (unsigned int index = 0;index < faces_count;++index)
        {
            short* i = ti_faces(index);
            short i1 = i[0];
            short i2 = i[1];
            short i3 = i[2];
            if (i1 >= half_odf_size)
                i1 -= half_odf_size;
            if (i2 >= half_odf_size)
                i2 -= half_odf_size;
            if (i3 >= half_odf_size)
                i3 -= half_odf_size;
            neighbor[i1].push_back(i2);
            neighbor[i1].push_back(i3);
            neighbor[i2].push_back(i1);
            neighbor[i2].push_back(i3);
            neighbor[i3].push_back(i1);
            neighbor[i3].push_back(i2);
        }
    }
    void search(const std::vector<float>& old_odf)
    {
        max_table.clear();
        for (unsigned int index = 0;index < neighbor.size();++index)
        {
            float value = old_odf[index];
            bool is_max = true;
            std::vector<unsigned short>& nei = neighbor[index];
            for (unsigned int j = 0;j < nei.size();++j)
            {
                if (value < old_odf[nei[j]])
                {
                    is_max = false;
                    break;
                }
            }
            if (is_max)
                max_table[value] = (unsigned short)index;
        }
    }
};


struct ODFDecomposition : public BaseProcess
{
protected:
    ConvolutedOdfComponent odf_component;
    boost::mutex mutex;
public:
    virtual void init(Voxel& voxel)
    {
        process_position = __FUNCTION__;
        if (voxel.odf_decomposition)
        {
            odf_component.icosa_components.resize(ti_vertices_count() >> 1);
            for (unsigned int index = 0;index < ti_vertices_count() >> 1;++index)
                odf_component.icosa_components[index].initialize(index);
        }
    }
    virtual void run(Voxel& voxel, VoxelData& data)
    {
        process_position = __FUNCTION__;
        if (voxel.odf_decomposition)
        {
            float min_odf = 0;
            boost::mutex::scoped_lock lock(mutex);
            std::vector<float> dodf(data.odf.size());
            min_odf = RemoveIsotropicPart()(data.odf);
            for (unsigned int index = 0;index < data.odf.size()/2;++index)
            {
                odf_component.decomposeODF(data.odf);
                unsigned short mv = odf_component.getMainVector();
                dodf[mv] += odf_component.getExtractedODF()[mv];
            }
            data.odf.swap(dodf);
        }
    }
};


struct DetermineFiberDirections : public BaseProcess
{
    SearchLocalMaximum lm;
    boost::mutex mutex;
public:
    virtual void init(Voxel& voxel)
    {
        process_position = __FUNCTION__;
        lm.init();
    }

    virtual void run(Voxel& voxel,VoxelData& data)
    {
        process_position = __FUNCTION__;
        boost::mutex::scoped_lock lock(mutex);
        data.min_odf = *std::min_element(data.odf.begin(),data.odf.end());
        lm.search(data.odf);
        std::map<float,unsigned short,std::greater<float> >::const_iterator iter = lm.max_table.begin();
        std::map<float,unsigned short,std::greater<float> >::const_iterator end = lm.max_table.end();
        for (unsigned int index = 0;iter != end && index < voxel.max_fiber_number;++index,++iter)
        {
            data.dir_index[index] = iter->second;
            data.fa[index] = iter->first - data.min_odf;
        }
    }
};


class AnaQSpace2Odf : public QSpace2Odf
{
public:
    float max_distribution;
public:
    std::vector<image::vector<3,float> > b_vec;
    struct gradient_fun
    {
        const std::vector<image::vector<3,float> >& b_vec;
        const std::vector<float>& signal;
        gradient_fun(const std::vector<image::vector<3,float> >& b_vec_,
                     const std::vector<float>& signal_):b_vec(b_vec_),signal(signal_) {}

        image::vector<3,float> operator()(const image::vector<3,float>& x0) const
        {
            image::vector<3,float> nx(x0);
            nx.normalize();
            image::vector<3,float> sum;
            for (unsigned int j = 0;j < b_vec.size();++j)
            {
                float prod = b_vec[j]*nx;
                float temp = std::cos(prod)-boost::math::sinc_pi(prod);
                if (std::abs(prod) < 0.0000001)
                    prod = prod > 0.0 ? 0.0000001:-0.0000001;
                prod = temp/prod;
                prod *= signal[j];

                image::vector<3,float> b_vec_ortho(b_vec[j]);
                image::vector<3,float> b_proj(nx);
                b_proj *= (b_vec_ortho*nx);
                b_vec_ortho -= b_proj;
                b_vec_ortho *= prod;
                sum += b_vec_ortho;
            }
            sum /= b_vec.size();
            return sum;
        }
    };

    struct fun
    {
        const std::vector<image::vector<3,float> >& b_vec;
        const std::vector<float>& signal;
        fun(const std::vector<image::vector<3,float> >& b_vec_,
            const std::vector<float>& signal_):b_vec(b_vec_),signal(signal_) {}

        float operator()(const image::vector<3,float>& x0) const
        {
            float sum = 0.0;
            image::vector<3,float> nx(x0);
            nx.normalize();
            for (unsigned int j = 0;j < b_vec.size();++j)
                sum += signal[j] * boost::math::sinc_pi(b_vec[j]*nx);
            return -sum;// to maximize the sum
        }
    };

public:
    virtual void init(Voxel& voxel)
    {
        QSpace2Odf::init(voxel);
        process_position = __FUNCTION__;

        b_vec.resize(voxel.bvalues.size());
        float sigma = voxel.param[0]; //optimal 1.24
        for (unsigned int index = 0; index < voxel.bvalues.size(); ++index)
        {
            b_vec[index] = voxel.bvectors[index];
            b_vec[index] *= std::sqrt(voxel.bvalues[index]*0.01506)*sigma; // £^G£_
        }
    }
    virtual void run(Voxel& voxel, VoxelData& data)
    {
        process_position = __FUNCTION__;

        gradient_fun g(b_vec,data.space);
        fun f(b_vec,data.space);
        std::map<float,image::vector<3,float>,std::greater<float> > result;
        for (unsigned int index = 0;index < voxel.max_fiber_number;++index)
        {
            if (data.fa[index] <= 0.0)
                break;

            image::vector<3,float> initial_dir(ti_vertices(data.dir_index[index]));
            image::optimization::gradient_descent(initial_dir,f,g,(double)0.01);
            double fa = -data.min_odf-f(initial_dir);
            initial_dir.normalize();
            if (std::abs(image::vector<3,float>(ti_vertices(data.dir_index[index])) * initial_dir) < 0.939692623) // the angle greater than 20 degrees
                result[data.fa[index]] = image::vector<3,float>(ti_vertices(data.dir_index[index]));
            else
                result[fa] = initial_dir;
        }

        std::map<float,image::vector<3,float>,std::greater<float> >::const_iterator iter = result.begin();
        std::map<float,image::vector<3,float>,std::greater<float> >::const_iterator end = result.end();
        for (unsigned int index = 0;iter != end;++index,++iter)
        {
            data.fa[index] = iter->first;
            data.dir[index] = iter->second;
        }

    }
};



#endif//ODF_TRANSFORMATION_PROCESS_HPP
