#ifndef TESSELLATED_ICOSAHEDRON_HPP
#define TESSELLATED_ICOSAHEDRON_HPP
#include <boost/lambda/lambda.hpp>
#include <image/image.hpp>
#include <cmath>
#include <vector>
#include <functional>
#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif


template<unsigned int fold>
class TessellatedIcosahedron
{
public:
    unsigned int vertex_count;
    unsigned int half_vertex_count;
private:
    float edge_length,face_dis,angle_res;
    unsigned short cur_vertex;
    void add_vertex(const image::vector<3,float>& vertex)
    {
        vertices[cur_vertex] = vertex;
        vertices[cur_vertex + half_vertex_count] = -vertex;
		++cur_vertex;
    }
        void get_mid_vertex(const image::vector<3,float>& v1,
                                            const image::vector<3,float>& v2,
                                                image::vector<3,float>& v3) const
	{
		v3 = v1;
		v3 += v2;
		v3.normalize();
	}
        void get_mid_vertex(const image::vector<3,float>& v1,
                                            const image::vector<3,float>& v2,
                                                const image::vector<3,float>& v3,
                                                image::vector<3,float>& v4) const
	{
		v4 = v1;
		v4 += v2;
		v4 += v3;
		v4.normalize();
	}
public:
    std::vector<image::vector<3,float> > vertices;
    std::vector<image::vector<3,short> > faces;
    // obtain the segmentation from vector u to v
    void get_edge_segmentation(unsigned short from_vertex,unsigned short to_vertex,
                               std::vector<unsigned short>& edge,unsigned int order)
    {
		edge.push_back(from_vertex);
        if (order == 2)
        {
                        image::vector<3,float> uv;
			get_mid_vertex(vertices[to_vertex],vertices[from_vertex],uv);
			edge.push_back(cur_vertex);
            add_vertex(uv);
        }
		else
		if(order == 4)
		{
                        image::vector<3,float> mid;
			get_mid_vertex(vertices[to_vertex],vertices[from_vertex],mid);
                        image::vector<3,float> u = mid;
			get_mid_vertex(mid,vertices[from_vertex],u);
			edge.push_back(cur_vertex);
            add_vertex(u);
			edge.push_back(cur_vertex);
            add_vertex(mid);
			get_mid_vertex(mid,vertices[to_vertex],u);
			edge.push_back(cur_vertex);
            add_vertex(u);
		}
        else
        {
            image::vector<3,float> uv = vertices[to_vertex];
            uv -= vertices[from_vertex];
            uv.normalize();
            for (unsigned int index = 1;index < order;++index)
            {
                float phi = angle_res*index;
                phi /= (float)order;
                float y = face_dis/std::cos(phi-angle_res*0.5);

                                image::vector<3,float> vec = uv;
                vec *= std::sqrt(1.0 + y*(y-2.0*std::cos(phi)));
                vec += vertices[from_vertex];
                vec.normalize();

                edge.push_back(cur_vertex);
                add_vertex(vec);
            }
        }
		edge.push_back(to_vertex);
    }
	unsigned short opposite(unsigned short v1) const
	{
		return (v1 < half_vertex_count) ? v1 + half_vertex_count:v1 - half_vertex_count;
	}
    void add_face(short v1,short v2,short v3)
    {
                faces.push_back(image::vector<3,short>(v1,v2,v3));
                faces.push_back(image::vector<3,short>(opposite(v1),opposite(v3),opposite(v2)));
    }
	template<typename input_type1,typename input_type2>
        void add_faces_in_line(input_type1 up,input_type2 down,unsigned int up_count)
	{
                for(unsigned int index = 0;index < up_count;++index)
		{
			add_face(up[index],down[index],down[index+1]);
			if(index + 1 < up_count)
			    add_face(up[index+1],up[index],down[index+1]);
		}
	}
	template<typename input_type1,typename input_type2,typename input_type3>
	void add_faces_in_triangle(input_type1 edge0,
							   input_type2 edge1,
                                                           input_type3 edge2,unsigned int folding)
    {
		if(folding <= 1)
			return;
		if(folding == 2)
		{
			add_face(edge2[1],edge0[0],edge0[1]);
			add_face(edge0[1],edge1[0],edge1[1]);
			add_face(edge1[1],edge2[0],edge2[1]);
			add_face(edge0[1],edge1[1],edge2[1]);
			return;
		}
		std::vector<unsigned short> mid_line(folding);
		mid_line.front() = edge2[folding-1];
		mid_line.back() = edge1[1];
                unsigned int old_cur_vertex = cur_vertex;
                for(unsigned int index = 1;index < mid_line.size()-1;++index,++cur_vertex)
		    mid_line[index] = cur_vertex;
		add_faces_in_triangle(mid_line,edge1+1,edge2,folding-1);
		add_faces_in_line(mid_line,edge0,folding);
		cur_vertex = old_cur_vertex;
	}
	// This function obtain the tessellated points and faces from a icosaheron triangle
    template<typename input_type1,typename input_type2,typename input_type3>
        void build_faces(input_type1 edge0,input_type2 edge1,input_type3 edge2,unsigned int folding)
    {
		add_faces_in_triangle(edge0,edge1,edge2,folding);
		if(folding == 3)
		{
                        image::vector<3,float> mid;
			get_mid_vertex(vertices[edge0[0]],vertices[edge1[0]],vertices[edge2[0]],mid);
			add_vertex(mid);
			return;
		}
		if(folding == 4)
		{
                        image::vector<3,float> u;
			get_mid_vertex(vertices[edge0[2]],vertices[edge2[2]],u);
			add_vertex(u);
			get_mid_vertex(vertices[edge0[2]],vertices[edge1[2]],u);
			add_vertex(u);
			get_mid_vertex(vertices[edge1[2]],vertices[edge2[2]],u);
			add_vertex(u);
			return;
		}
		if(folding == 5)
		{
                        image::vector<3,float> u1,u2,u3,u4,u5,u6;
			get_mid_vertex(vertices[edge0[0]],vertices[edge0[3]],vertices[edge2[2]],u1);
			get_mid_vertex(vertices[edge1[0]],vertices[edge1[3]],vertices[edge0[2]],u3);
			get_mid_vertex(vertices[edge2[0]],vertices[edge2[3]],vertices[edge1[2]],u6);
			get_mid_vertex(u1,u3,u2);
			get_mid_vertex(u3,u6,u5);
			get_mid_vertex(u6,u1,u4);
			add_vertex(u1);
			add_vertex(u2);
			add_vertex(u3);
			add_vertex(u4);
			add_vertex(u5);
			add_vertex(u6);
			return;
		}
		if(folding == 6)
		{
                        image::vector<3,float> u1,u2,u3,u4,u5,u6,u7,u8,u9,u10;// (6-1)*(6-2)/2
			get_mid_vertex(vertices[edge0[0]],vertices[edge1[0]],vertices[edge2[0]],u6);

			get_mid_vertex(vertices[edge0[0]],vertices[edge0[3]],vertices[edge2[3]],u1);
			get_mid_vertex(vertices[edge1[0]],vertices[edge1[3]],vertices[edge0[3]],u4);
			get_mid_vertex(vertices[edge2[0]],vertices[edge2[3]],vertices[edge1[3]],u10);

			get_mid_vertex(u6,vertices[edge0[2]],u2);
			get_mid_vertex(u6,vertices[edge0[4]],u3);

			get_mid_vertex(u6,vertices[edge2[4]],u5);
			get_mid_vertex(u6,vertices[edge1[2]],u7);

			get_mid_vertex(u6,vertices[edge2[2]],u8);
			get_mid_vertex(u6,vertices[edge1[4]],u9);

			add_vertex(u1);
			add_vertex(u2);
			add_vertex(u3);
			add_vertex(u4);
			add_vertex(u5);
			add_vertex(u6);
			add_vertex(u7);
			add_vertex(u8);
			add_vertex(u9);
			add_vertex(u10);
			return;
		}
		if(folding == 8)
		{
                        image::vector<3,float> u[22]; // (8-1)*(8-2)/2=21, add one for index started from1
			get_mid_vertex(vertices[edge0[4]],vertices[edge2[4]],u[8]);
			get_mid_vertex(vertices[edge1[4]],vertices[edge2[4]],u[17]);
			get_mid_vertex(vertices[edge0[4]],vertices[edge1[4]],u[10]);

			get_mid_vertex(vertices[edge0[2]],vertices[edge2[6]],u[1]);
			get_mid_vertex(vertices[edge0[2]],u[8],u[2]);
			get_mid_vertex(vertices[edge0[4]],u[8],u[3]);
			get_mid_vertex(vertices[edge0[4]],u[10],u[4]);
			get_mid_vertex(vertices[edge0[6]],u[10],u[5]);
			get_mid_vertex(vertices[edge0[6]],vertices[edge1[2]],u[6]);

			get_mid_vertex(vertices[edge2[6]],u[8],u[7]);
			get_mid_vertex(u[10],u[8],u[9]);
			get_mid_vertex(vertices[edge1[2]],u[10],u[11]);

			get_mid_vertex(vertices[edge2[4]],u[8],u[12]);
			get_mid_vertex(u[17],u[8],u[13]);
			get_mid_vertex(u[17],u[10],u[14]);
			get_mid_vertex(vertices[edge1[4]],u[10],u[15]);

			get_mid_vertex(vertices[edge2[4]],u[17],u[16]);
			get_mid_vertex(vertices[edge1[4]],u[17],u[18]);

			get_mid_vertex(vertices[edge2[2]],u[17],u[19]);
			get_mid_vertex(vertices[edge1[6]],u[17],u[20]);
			get_mid_vertex(vertices[edge2[2]],vertices[edge1[6]],u[21]);
                        for(unsigned int index = 1;index <= 21;++index)
				add_vertex(u[index]);
			return;
		}
    }

    void build_icosahedron(void)
    {
        // the top vertex
        add_vertex(image::vector<3,float>(0,0,1.0));
        //central vertices around the upper staggered circles
        float sqrt5 = std::sqrt(5.0);
        float height = 1.0/sqrt5;
        float radius = 2.0/sqrt5;
        for (unsigned int index = 0;index < 5;++index)
            add_vertex(image::vector<3,float>(
                                   std::cos(((float)index)*M_PI*0.4)*radius,
                                   std::sin(((float)index)*M_PI*0.4)*radius,height));

        edge_length = std::sqrt(2.0-2.0/sqrt5);
        face_dis = std::sqrt(0.5+0.5/sqrt5);
        angle_res = std::acos(1.0/sqrt5);

        std::vector<std::vector<unsigned short> > edges(15);
        // top hat
        get_edge_segmentation(0,1,edges[0],fold);
        get_edge_segmentation(0,2,edges[1],fold);
        get_edge_segmentation(0,3,edges[2],fold);
        get_edge_segmentation(0,4,edges[3],fold);
        get_edge_segmentation(0,5,edges[4],fold);
        // the edge of hat
        get_edge_segmentation(1,2,edges[5],fold);
        get_edge_segmentation(2,3,edges[6],fold);
        get_edge_segmentation(3,4,edges[7],fold);
        get_edge_segmentation(4,5,edges[8],fold);
        get_edge_segmentation(5,1,edges[9],fold);
        // skirt
        get_edge_segmentation(1,opposite(4),edges[10],fold);
        get_edge_segmentation(2,opposite(5),edges[11],fold);
        get_edge_segmentation(3,opposite(1),edges[12],fold);
        get_edge_segmentation(4,opposite(2),edges[13],fold);
        get_edge_segmentation(5,opposite(3),edges[14],fold);

		std::vector<std::vector<unsigned short> > redges(15);
                for(unsigned int index = 0;index < redges.size();++index)
		    std::transform(edges[index].begin(),
						   edges[index].end(),
						   std::back_inserter(redges[index]),
						   (boost::lambda::_1 + half_vertex_count) % vertex_count);

		// hat faces
		build_faces(edges[0].begin(),edges[5].begin(),edges[1].rbegin(),fold);
		build_faces(edges[1].begin(),edges[6].begin(),edges[2].rbegin(),fold);
		build_faces(edges[2].begin(),edges[7].begin(),edges[3].rbegin(),fold);
		build_faces(edges[3].begin(),edges[8].begin(),edges[4].rbegin(),fold);
		build_faces(edges[4].begin(),edges[9].begin(),edges[0].rbegin(),fold);
		// skirt faces
		build_faces(edges[10].begin(),redges[13].begin(),edges[5].rbegin(),fold);
		build_faces(edges[11].begin(),redges[14].begin(),edges[6].rbegin(),fold);
		build_faces(edges[12].begin(),redges[10].begin(),edges[7].rbegin(),fold);
		build_faces(edges[13].begin(),redges[11].begin(),edges[8].rbegin(),fold);
		build_faces(edges[14].begin(),redges[12].begin(),edges[9].rbegin(),fold);

    }
	void check_vertex(void)
	{

		std::vector<float> min_cos(vertices.size());
                for(unsigned int i = 0;i < vertex_count;++i)
		{
			float value = 0.0;
                        for(unsigned int j = 0;j < vertex_count;++j)
			if(j != i && j != opposite(i) && std::abs(vertices[i]*vertices[j]) > value)
				value = std::abs(vertices[i]*vertices[j]);
			min_cos[i] = value;
		}

	}
	void check_face(void)
	{
                std::vector<unsigned int> count(vertex_count);
		std::vector<float> dis;
                for(unsigned int index = 0;index < faces.size();++index)
		{
		    ++count[faces[index][0]];
		    ++count[faces[index][1]];
		    ++count[faces[index][2]];
			dis.push_back(vertices[faces[index][0]]*vertices[faces[index][1]]);
			dis.push_back(vertices[faces[index][1]]*vertices[faces[index][2]]);
			dis.push_back(vertices[faces[index][2]]*vertices[faces[index][0]]);
		}
	}
	void sort_vertices(void)
	{
                std::vector<image::vector<3,float> > sorted_vertices(vertices.begin(),vertices.end());
                std::sort(sorted_vertices.begin(),sorted_vertices.end(),std::greater<image::vector<3,float> >());
                for(unsigned int index = 0;index < half_vertex_count;++index)
			sorted_vertices[index+half_vertex_count] = -sorted_vertices[index];
		std::vector<unsigned short> index_map(vertex_count);
                for(unsigned int i = 0;i < vertex_count;++i)
		{
                        for(unsigned int j = 0;j < vertex_count;++j)
				if(vertices[i] == sorted_vertices[j])
				{
					index_map[i] = j;
					break;
				}
		}
                for(unsigned int index = 0;index < faces.size();++index)
		{
			faces[index][0] = index_map[faces[index][0]];
			faces[index][1] = index_map[faces[index][1]];
			faces[index][2] = index_map[faces[index][2]];
		}
		sorted_vertices.swap(vertices);
		std::sort(faces.begin(),faces.end());
	}
public:
    TessellatedIcosahedron(void):
		vertex_count(fold*fold*10+2),
		half_vertex_count(vertex_count >> 1),
		cur_vertex(0),vertices(vertex_count)
    {
        build_icosahedron();
		sort_vertices();
		#ifdef _DEBUG
		check_vertex();
		check_face();
		#endif

    }
};

#endif//TESSELLATED_ICOSAHEDRON_HPP
