#ifndef INDEX_ALGORITHM_HPP
#define INDEX_ALGORITHM_HPP
#include <vector>
#include <deque>
#include <map>
#include <algorithm>
#include "image/utility/pixel_index.hpp"
#include "image/utility/basic_image.hpp"

namespace image
{
/**
    connected neighbors
    0 1 0
    1 x 1
    0 1 0
*/
inline void get_connected_neighbors(const pixel_index<2>& index,const geometry<2>& geo,
                                        std::vector<pixel_index<2> >& iterations)
{
    iterations.clear();
    iterations.reserve(4);
    if (index.x() >= 1)
        iterations.push_back(pixel_index<2>(index.x()-1,index.y(),index.index()-1));

    if (index.x()+1 < geo.width())
        iterations.push_back(pixel_index<2>(index.x()+1,index.y(),index.index()+1));

    if (index.y() >= 1)
        iterations.push_back(pixel_index<2>(index.x(),index.y()-1,index.index()-geo.width()));

    if (index.y()+1 < geo.height())
        iterations.push_back(pixel_index<2>(index.x(),index.y()+1,index.index()+geo.width()));
}

/**
    connected neighbors
    0 1 0
    1 x 1
    0 1 0
*/
inline void get_connected_neighbors(const pixel_index<3>& index,const geometry<3>& geo,
                                        std::vector<pixel_index<3> >& iterations)
{
    iterations.clear();
    iterations.reserve(6);

    if (index.x() >= 1)
        iterations.push_back(pixel_index<3>(index.x()-1,index.y(),index.z(),index.index()-1));

    if (index.x()+1 < geo.width())
        iterations.push_back(pixel_index<3>(index.x()+1,index.y(),index.z(),index.index()+1));

    if (index.y() >= 1)
        iterations.push_back(pixel_index<3>(index.x(),index.y()-1,index.z(),index.index()-geo.width()));

    if (index.y()+1 < geo.height())
        iterations.push_back(pixel_index<3>(index.x(),index.y()+1,index.z(),index.index()+geo.width()));

    if (index.z() >= 1)
        iterations.push_back(pixel_index<3>(index.x(),index.y(),index.z()-1,index.index()-geo.plane_size()));

    if (index.z()+1 < geo.depth())
        iterations.push_back(pixel_index<3>(index.x(),index.y(),index.z()+1,index.index()+geo.plane_size()));
}

/**
    1 1 1
    1 x 1
    1 1 1
*/
inline void get_neighbors(const pixel_index<2>& index,const geometry<2>& geo,
                                        std::vector<pixel_index<2> >& iterations)
{
    iterations.clear();
    iterations.reserve(8);
    bool has_left = index.x() >= 1;
    bool has_right = index.x()+1 < geo.width();
    unsigned int x_left,x_right;
    if(has_left)
        x_left = index.x()-1;
    if(has_right)
        x_right = index.x()+1;
    if (index.y() >= 1)
    {
        unsigned int y_top = index.y()-1;
        unsigned int base_index = index.index()-geo.width();
        if (has_left)
            iterations.push_back(pixel_index<2>(x_left,y_top,base_index-1));

        iterations.push_back(pixel_index<2>(index.x()  ,y_top,base_index));
        if (has_right)
            iterations.push_back(pixel_index<2>(x_right,y_top,base_index+1));
    }
    {
        if (has_left)
            iterations.push_back(pixel_index<2>(x_left,index.y(),index.index()-1));

        //iterations.push_back(pixel_index<2>(index.x()  ,index.y(),index.index()));
        if (has_right)
            iterations.push_back(pixel_index<2>(x_right,index.y(),index.index()+1));
    }
    if (index.y()+1 < geo.height())
    {
        unsigned int y_bottom = index.y()+1;
        unsigned int base_index = index.index()+geo.width();
        if (has_left)
            iterations.push_back(pixel_index<2>(x_left,y_bottom,base_index-1));

        iterations.push_back(pixel_index<2>(index.x()  ,y_bottom,base_index));
        if (has_right)
            iterations.push_back(pixel_index<2>(x_right,y_bottom,base_index+1));
    }

}


inline void get_neighbors(const pixel_index<3>& index,const geometry<3>& geo,
                                        std::vector<pixel_index<3> >& iterations)
{
    iterations.clear();
    iterations.reserve(26);
    unsigned int z_offset = geo.plane_size();
    unsigned int y_offset = geo.width();
    bool has_left = index.x() >= 1;
    bool has_right = index.x()+1 < geo.width();
    bool has_top = index.y() >= 1;
    bool has_bottom = index.y()+1 < geo.height();
    unsigned int x_left,x_right,y_top,y_bottom;
    if(has_left)
        x_left = index.x()-1;
    if(has_right)
        x_right = index.x()+1;
    if(has_top)
        y_top = index.y()-1;
    if(has_bottom)
        y_bottom = index.y()+1;
    if (index.z() >= 1)
    {
        unsigned int z =  index.z()-1;
        unsigned int base_index = index.index()-z_offset;
        if (has_top)
        {
            unsigned int base_index2 = base_index - y_offset;
            if (has_left)
                iterations.push_back(pixel_index<3>(x_left,y_top,z,base_index2-1));

            iterations.push_back(pixel_index<3>(index.x()  ,y_top,z,base_index2));

            if (has_right)
                iterations.push_back(pixel_index<3>(x_right,y_top,z,base_index2+1));
        }
        {
            if (has_left)
                iterations.push_back(pixel_index<3>(x_left,index.y(),z,base_index-1));

            iterations.push_back(pixel_index<3>(index.x()  ,index.y(),z,base_index));

            if (has_right)
                iterations.push_back(pixel_index<3>(x_right,index.y(),z,base_index+1));
        }
        if (has_bottom)
        {
            unsigned int base_index2 = base_index + y_offset;
            if (has_left)
                iterations.push_back(pixel_index<3>(x_left,y_bottom,z,base_index2-1));

            iterations.push_back(pixel_index<3>(index.x()  ,y_bottom,z,base_index2));

            if (has_right)
                iterations.push_back(pixel_index<3>(x_right,y_bottom,z,base_index2+1));
        }
    }

    {
        if (has_top)
        {
            unsigned int base_index2 = index.index() - y_offset;
            if (has_left)
                iterations.push_back(pixel_index<3>(x_left,y_top,index.z(),base_index2-1));

            iterations.push_back(pixel_index<3>(index.x()  ,y_top,index.z(),base_index2));

            if (has_right)
                iterations.push_back(pixel_index<3>(x_right,y_top,index.z(),base_index2+1));
        }
        {
            if (has_left)
                iterations.push_back(pixel_index<3>(x_left,index.y(),index.z(),index.index()-1));

            //iterations.push_back(pixel_index<3>(index.x()  ,index.y(),index.z(),index.index()  ));

            if (has_right)
                iterations.push_back(pixel_index<3>(x_right,index.y(),index.z(),index.index()+1));
        }
        if (has_bottom)
        {
            unsigned int base_index2 = index.index() + y_offset;
            if (has_left)
                iterations.push_back(pixel_index<3>(x_left,y_bottom,index.z(),base_index2-1));

            iterations.push_back(pixel_index<3>(index.x()  ,y_bottom,index.z(),base_index2));

            if (has_right)
                iterations.push_back(pixel_index<3>(x_right,y_bottom,index.z(),base_index2+1));
        }

    }
    if (index.z()+1 < geo.depth())
    {
        unsigned int z = index.z()+1;
        unsigned int base_index = index.index()+z_offset;
        if (has_top)
        {
            unsigned int base_index2 = base_index - y_offset;
            if (has_left)
                iterations.push_back(pixel_index<3>(x_left,y_top,z,base_index2-1));

            iterations.push_back(pixel_index<3>(index.x()  ,y_top,z,base_index2));

            if (has_right)
                iterations.push_back(pixel_index<3>(x_right,y_top,z,base_index2+1));
        }
        {
            if (has_left)
                iterations.push_back(pixel_index<3>(x_left,index.y(),z,base_index-1));

            iterations.push_back(pixel_index<3>(index.x()  ,index.y(),z,base_index));

            if (has_right)
                iterations.push_back(pixel_index<3>(x_right,index.y(),z,base_index+1));
        }
        if (has_bottom)
        {
            unsigned int base_index2 = base_index + y_offset;
            if (has_left)
                iterations.push_back(pixel_index<3>(x_left,y_bottom,z,base_index2-1));

            iterations.push_back(pixel_index<3>(index.x()  ,y_bottom,z,base_index2));

            if (has_right)
                iterations.push_back(pixel_index<3>(x_right,y_bottom,z,base_index2+1));
        }
    }

}


template<int Dim>
inline void get_neighbors(const pixel_index<Dim>& index,const geometry<Dim>& geo,unsigned int range,std::vector<pixel_index<Dim> >& iterations)
{
    iterations.clear();
    iterations.reserve(9);
    throw;
}

inline void get_neighbors(const pixel_index<2>& index,const geometry<2>& geo,unsigned int range,std::vector<pixel_index<2> >& iterations)
{
    iterations.clear();
    iterations.reserve(9);
    unsigned int fx = (index.x() > range) ? index.x() - range:0;
    unsigned int fy = (index.y() > range) ? index.y() - range:0;
    unsigned int tx = std::min<unsigned int>(index.x() + range,geo.width()-1);
    unsigned int ty = std::min<unsigned int>(index.y() + range,geo.height()-1);
    unsigned int y_index = fy*geo.width()+fx;
    for (unsigned int y = fy;y <= ty;++y,y_index += geo.width())
    {
        unsigned int x_index = y_index;
        for (unsigned int x = fx;x <= tx;++x,++x_index)
            iterations.push_back(pixel_index<2>(x,y,x_index));
    }
}

inline void get_neighbors(const pixel_index<3>& index,const geometry<3>& geo,unsigned int range,std::vector<pixel_index<3> >& iterations)
{
    iterations.clear();
    iterations.reserve(26);
    unsigned int wh = geo.plane_size();
    unsigned int fx = (index.x() > range) ? index.x() - range:0;
    unsigned int fy = (index.y() > range) ? index.y() - range:0;
    unsigned int fz = (index.z() > range) ? index.z() - range:0;
    unsigned int tx = std::min<unsigned int>(index.x() + range,geo.width()-1);
    unsigned int ty = std::min<unsigned int>(index.y() + range,geo.height()-1);
    unsigned int tz = std::min<unsigned int>(index.z() + range,geo.depth()-1);
    unsigned int z_index = (fz*geo.height()+fy)*geo.width()+fx;
    for (unsigned int z = fz;z <= tz;++z,z_index += wh)
    {
        unsigned int y_index = z_index;
        for (unsigned int y = fy;y <= ty;++y,y_index += geo.width())
        {
            unsigned int x_index = y_index;
            for (unsigned int x = fx;x <= tx;++x,++x_index)
                iterations.push_back(pixel_index<3>(x,y,z,x_index));
        }
    }
}


template<unsigned int dim>
class neighbor_index_shift;

template<>
class neighbor_index_shift<2>
{
public:
    std::vector<int> index_shift;
public:
    neighbor_index_shift(const geometry<2>& geo)
    {
        int w = geo.width();
            for (int y = -1;y <= 1; ++y)
            {
                int yw = y*w;
                for (int x = -1;x <= 1; ++x)
                    index_shift.push_back(x + yw);
            }
    }
};


template<>
class neighbor_index_shift<3>
{
public:
    std::vector<int> index_shift;
public:
    neighbor_index_shift(const geometry<3>& geo)
    {
        int wh = geo.plane_size();
        int w = geo.width();
        for (int z = -1;z <= 1; ++z)
        {
            int zwh = z*wh;
            for (int y = -1;y <= 1; ++y)
            {
                int yw = y*w;
                for (int x = -1;x <= 1; ++x)
                    index_shift.push_back(x + yw + zwh);
            }
        }
    }
};


template<unsigned int dim>
class neighbor_index_shift_narrow;

template<>
class neighbor_index_shift_narrow<2>
{
public:
    std::vector<int> index_shift;
public:
    neighbor_index_shift_narrow(const geometry<2>& geo)
    {
        index_shift.push_back(-geo.width());
        index_shift.push_back(-1);
        index_shift.push_back(0);
        index_shift.push_back(1);
        index_shift.push_back(geo.width());
    }
};


template<>
class neighbor_index_shift_narrow<3>
{
public:
    std::vector<int> index_shift;
public:
    neighbor_index_shift_narrow(const geometry<3>& geo)
    {
        index_shift.push_back(-geo.plane_size());
        index_shift.push_back(-geo.width());
        index_shift.push_back(-1);
        index_shift.push_back(0);
        index_shift.push_back(1);
        index_shift.push_back(geo.width());
        index_shift.push_back(geo.plane_size());
    }
};



}
#endif
