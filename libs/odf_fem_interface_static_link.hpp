#ifndef ODF_FEM_INTERFACE_STATIC_LINK_H
#define ODF_FEM_INTERFACE_STATIC_LINK_H


extern "C"{
void ti_load_odf_vertices(unsigned int vertices_count_,const float* odf_buffer);
void ti_load_odf_faces(unsigned int faces_count_,const short* face_buffer);
void ti_initialize(unsigned int fold);
unsigned int ti_vertices_count(void);
unsigned int ti_faces_count(void);
float* ti_vertices(unsigned int index);
short ti_vertices_pair(short v1);
float ti_vertices_cos(unsigned int v1,unsigned int v2);
short* ti_faces(unsigned int index);
short ti_discretize(float vx,float vy,float vz);

}

#endif
