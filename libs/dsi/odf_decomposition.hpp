#ifndef ODF_DECOMPOSITION_HPP
#define ODF_DECOMPOSITION_HPP
#include <boost/lambda/lambda.hpp>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <map>
#include "common.hpp"

/*
class BasicSamplingModel
{
public:
    static const unsigned int layer_count = 20;
    static const unsigned int main_equater_resolution = layer_count*4;
private:
    std::vector<image::vector<3,float> > base_sampling_points;
    std::vector<unsigned int> layers_index;
public:
    BasicSamplingModel(void)
    {
        base_sampling_points.clear();
        base_sampling_points.push_back(image::vector<3,float>(0.0,0.0,1.0));
        for (unsigned int layer_index = 0;layer_index < layer_count;++layer_index)
        {
            float theta = M_PI/2.0/((float)layer_count);
            theta *= layer_index+1;
            float z = std::cos(theta);
            float r = std::sin(theta);
            unsigned int equater_resolution = ((float)main_equater_resolution)*r;
            if (equater_resolution < 2)
                equater_resolution = 2;

            for (unsigned int index = 0;index < equater_resolution;++index)
            {
                float phi = M_PI*2.0/((float)equater_resolution);
                phi *= index;
                base_sampling_points.push_back(image::vector<3,float>(std::cos(phi)*r,std::sin(phi)*r,z));
            }
            layers_index.push_back(base_sampling_points.size());
        }
    }
    unsigned int size(void) const
    {
        return base_sampling_points.size();
    }
    const image::vector<3,float>& operator[](unsigned int index) const
    {
        return base_sampling_points[index];
    }
    unsigned int layer_at(unsigned int index) const
    {
        return layers_index[index];
    }
};
*/


class FiberDistribution
{
public:
    static const unsigned int discrete_count = 45;
    // the resolution of the response function from 0 degree to 90 degree
    std::vector<unsigned int> odf_label;
    std::vector<float> odf_position;
public:
    void initialize(unsigned int axis_index)
    {
        odf_label.resize(ti_vertices_count() >> 1);
        odf_position.resize(ti_vertices_count() >> 1);
        for(unsigned int index = 0; index < ti_vertices_count() >> 1; ++index)
        {
            float cos_value = std::abs(ti_vertices_cos(index,axis_index));
            if(cos_value > 1.0)
                cos_value = 1.0;
            odf_position[index] = ((float)(discrete_count))*std::acos(cos_value)*2.0/M_PI;
            odf_label[index] = (unsigned int)std::floor(odf_position[index]+0.5);
        }
    }
    void get_response_function(unsigned int axis_index,const std::vector<float>& odf,std::vector<float>& response_function) const
    {
        response_function.resize(discrete_count+2); // The response function of a fiber

        std::fill(response_function.begin(),response_function.end(),odf[axis_index]/2.0);
        //std::fill(response_function.begin(),response_function.end(),0);
        // axial maximum
        for(unsigned int index = 0; index < odf_label.size(); ++index)
            if(response_function[odf_label[index]] > odf[index]/2.0)
                response_function[odf_label[index]] = odf[index]/2.0;
        // centrifugal decreasing
        float dis = response_function.front();
        for(unsigned int index = 1; index < response_function.size(); ++index)
            if(dis < response_function[index])
                response_function[index] = dis;
            else
                dis = response_function[index];
    }

    void decompose(unsigned int axis_index,std::vector<float>& odf,std::vector<float>& decomposed_odf) const
    {
        std::vector<float> response_function; // The response function of a fiber
        get_response_function(axis_index,odf,response_function);

        for(unsigned int index = 0; index < odf.size(); ++index)
        {
            float dis = estimate_dis(response_function,index);
            if(dis > odf[index])
                dis = odf[index];
            decomposed_odf[index] = dis;
            odf[index] -= dis;
        }
    }
    float estimate_dis(const std::vector<float>& response_function,unsigned int odf_index) const
    {
        float position = odf_position[odf_index];
        float floor = std::floor(position);
        float r = position-floor;
        unsigned int dis_index = floor;
        return response_function[dis_index]*(1.0-r) + response_function[dis_index+1]*r;
    }
};


class ConvolutedOdfComponent
{
public:
    std::vector<FiberDistribution> icosa_components;
private:
    unsigned short main_vector_index;
    std::vector<float> extracted_odf;
public:
    void decomposeODF(std::vector<float>& odf,unsigned int main_index_)
    {
        main_vector_index = main_index_;
        extracted_odf.resize(odf.size());
        icosa_components[main_vector_index].decompose(main_vector_index,odf,extracted_odf);
    }
    void decomposeODF(std::vector<float>& odf)
    {
        decomposeODF(odf,std::max_element(odf.begin(),odf.end())-odf.begin());
    }
    void get_response_function(const std::vector<float>& odf,std::vector<float>& function)
    {
        main_vector_index = std::max_element(odf.begin(),odf.end())-odf.begin();
        icosa_components[main_vector_index].get_response_function(main_vector_index,odf,function);
    }
    const std::vector<float>& getExtractedODF(void) const
    {
        return extracted_odf;
    }
    unsigned short getMainVector(void) const
    {
        return main_vector_index;
    }

};

#endif//ODF_DECOMPOSITION_HPP
