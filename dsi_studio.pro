QT += core \
    gui \
    opengl \
    charts
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
CONFIG += c++11
#CONFIG += console
TARGET = dsi_studio
TEMPLATE = app
INCLUDEPATH += ./plot



win32* {
# GPU computation
# LIBS += -L"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v8.0\lib\x64" -lcudart_static -lcublas
# INCLUDEPATH += "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v8.0\include"

INCLUDEPATH += ../include
QMAKE_CXXFLAGS += -wd4244 -wd4267 -wd4018
LIBS += -lOpenGL32 -lGlu32
RC_FILE = dsi_studio.rc
}

linux* {
QMAKE_CXXFLAGS += -fpermissive
LIBS += lGLU -lz
}

mac{

INCLUDEPATH += /Users/admin/include
LIBS += -lz
ICON = dsi_studio.icns
}


INCLUDEPATH += libs \
    libs/dsi \
    libs/tracking \
    libs/mapping
HEADERS += mainwindow.h \
    dicom/dicom_parser.h \
    dicom/dwi_header.hpp \
    libs/dsi/tessellated_icosahedron.hpp \
    libs/dsi/space_mapping.hpp \
    libs/dsi/sh_process.hpp \
    libs/dsi/sample_model.hpp \
    libs/dsi/racian_noise.hpp \
    libs/dsi/qbi_process.hpp \
    libs/dsi/odf_process.hpp \
    libs/dsi/odf_deconvolusion.hpp \
    libs/dsi/odf_decomposition.hpp \
    libs/dsi/mix_gaussian_model.hpp \
    libs/dsi/layout.hpp \
    libs/dsi/image_model.hpp \
    libs/dsi/gqi_process.hpp \
    libs/dsi/gqi_mni_reconstruction.hpp \
    libs/dsi/dti_process.hpp \
    libs/dsi/dsi_process.hpp \
    libs/dsi/basic_voxel.hpp \
    SliceModel.h \
    tracking/tracking_window.h \
    reconstruction/reconstruction_window.h \
    tracking/slice_view_scene.h \
    opengl/glwidget.h \
    libs/tracking/tracking_method.hpp \
    libs/tracking/roi.hpp \
    libs/tracking/interpolation_process.hpp \
    libs/tracking/fib_data.hpp \
    libs/tracking/basic_process.hpp \
    libs/tracking/tract_cluster.hpp \
    tracking/region/regiontablewidget.h \
    tracking/region/Regions.h \
    tracking/region/RegionModel.h \
    libs/tracking/tract_model.hpp \
    tracking/tract/tracttablewidget.h \
    opengl/renderingtablewidget.h \
    qcolorcombobox.h \
    libs/tracking/tracking_thread.hpp \
    libs/prog_interface_static_link.h \
    simulation.h \
    libs/mapping/atlas.hpp \
    view_image.h \
    libs/gzip_interface.hpp \
    libs/dsi/racian_noise.hpp \
    libs/dsi/mix_gaussian_model.hpp \
    libs/dsi/layout.hpp \
    manual_alignment.h \
    tracking/tract_report.hpp \
    tracking/color_bar_dialog.hpp \
    tracking/connectivity_matrix_dialog.h \
    tracking/atlasdialog.h \
    filebrowser.h \
    program_option.hpp \
    qcompletelineedit.h \
    libs/mapping/connectometry_db.hpp \
    connectometry/createdbdialog.h \
    libs/dsi/qbi_voxel.hpp \
    connectometry/individual_connectometry.hpp \
    connectometry/match_db.h \
    connectometry/db_window.h \
    connectometry/group_connectometry.hpp \
    connectometry/group_connectometry_analysis.h \
    regtoolbox.h \
    connectometry/nn_connectometry.h \
    connectometry/nn_connectometry_analysis.h

FORMS += mainwindow.ui \
    tracking/tracking_window.ui \
    reconstruction/reconstruction_window.ui \
    dicom/dicom_parser.ui \
    simulation.ui \
    view_image.ui \
    manual_alignment.ui \
    tracking/tract_report.ui \
    tracking/color_bar_dialog.ui \
    tracking/connectivity_matrix_dialog.ui \
    tracking/atlasdialog.ui \
    filebrowser.ui \
    connectometry/createdbdialog.ui \
    connectometry/individual_connectometry.ui \
    connectometry/match_db.ui \
    connectometry/db_window.ui \
    connectometry/group_connectometry.ui \
    regtoolbox.ui \
    connectometry/nn_connectometry.ui
RESOURCES += \
    icons.qrc
SOURCES += main.cpp \
    mainwindow.cpp \
    dicom/dicom_parser.cpp \
    dicom/dwi_header.cpp \
    libs/utility/prog_interface.cpp \
    libs/dsi/sample_model.cpp \
    libs/dsi/dsi_interface_imp.cpp \
    libs/tracking/interpolation_process.cpp \
    libs/tracking/tract_cluster.cpp \
    SliceModel.cpp \
    tracking/tracking_window.cpp \
    reconstruction/reconstruction_window.cpp \
    tracking/slice_view_scene.cpp \
    opengl/glwidget.cpp \
    tracking/region/regiontablewidget.cpp \
    tracking/region/Regions.cpp \
    tracking/region/RegionModel.cpp \
    libs/tracking/tract_model.cpp \
    tracking/tract/tracttablewidget.cpp \
    opengl/renderingtablewidget.cpp \
    qcolorcombobox.cpp \
    cmd/trk.cpp \
    cmd/rec.cpp \
    simulation.cpp \
    cmd/src.cpp \
    libs/mapping/atlas.cpp \
    cmd/ana.cpp \
    view_image.cpp \
    manual_alignment.cpp \
    tracking/tract_report.cpp \
    tracking/color_bar_dialog.cpp \
    cmd/exp.cpp \
    tracking/connectivity_matrix_dialog.cpp \
    libs/dsi/tessellated_icosahedron.cpp \
    cmd/atl.cpp \
    tracking/atlasdialog.cpp \
    cmd/cnt.cpp \
    cmd/vis.cpp \
    filebrowser.cpp \
    qcompletelineedit.cpp \
    libs/tracking/fib_data.cpp \
    libs/tracking/tracking_thread.cpp \
    cmd/ren.cpp \
    libs/mapping/connectometry_db.cpp \
    connectometry/createdbdialog.cpp \
    connectometry/individual_connectometry.cpp \
    connectometry/match_db.cpp \
    connectometry/db_window.cpp \
    connectometry/group_connectometry.cpp \
    connectometry/group_connectometry_analysis.cpp \
    regtoolbox.cpp \
    cmd/cnn.cpp \
    cmd/qc.cpp \
    libs/dsi/basic_voxel.cpp \
    libs/dsi/image_model.cpp \
    connectometry/nn_connectometry.cpp \
    connectometry/nn_connectometry_analysis.cpp

OTHER_FILES += \
    options.txt \
    dicom_tag.txt \
    FreeSurferColorLUT.txt \
    shader_fragment.txt \
    shader_vertex.txt

DISTFILES += \
    shader_fragment2.txt \
    shader_vertex2.txt \

