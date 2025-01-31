#include <iostream>
#include <iterator>
#include <string>
#include <cstdio>
#include <QApplication>
#include <QMessageBox>
#include <QStyleFactory>
#include <QDir>
#include "mainwindow.h"
#include "tipl/tipl.hpp"
#include "mapping/atlas.hpp"
#include <iostream>
#include <iterator>
#include "program_option.hpp"
#include "cmd/cnt.cpp" // Qt project cannot build cnt.cpp without adding this.

std::string arg_file_name;
std::string
        fib_template_file_name_2mm,
        device_content_file;
std::vector<std::string> fa_template_list,
                         iso_template_list,
                         track_atlas_file_list;
std::vector<std::vector<std::string> > template_atlas_list;


int rec(void);
int trk(void);
int src(void);
int ana(void);
int exp(void);
int atl(void);
int cnt(void);
int cnt_ind(void);
int vis(void);
int ren(void);
int cnn(void);
int qc(void);
int reg(void);
int atk(void);


size_t match_template(float volume)
{
    float min_dif = std::numeric_limits<float>::max();
    size_t matched_index = 0;
    for(size_t i = 0;i < fa_template_list.size();++i)
    {
        gz_nifti read;
        if(!read.load_from_file(fa_template_list[i].c_str()))
            continue;
        float v = float(read.nif_header2.dim[1]*read.nif_header2.dim[2]*read.nif_header2.dim[3])*
                float(read.nif_header2.pixdim[1]*read.nif_header2.pixdim[2]*read.nif_header2.pixdim[3]);
        v = std::fabs(v-volume);
        if(v < min_dif)
        {
            min_dif = v;
            matched_index = i;
        }
    }
    return matched_index;
}

QStringList search_files(QString dir,QString filter)
{
    QStringList dir_list,src_list;
    dir_list << dir;
    for(unsigned int i = 0;i < dir_list.size();++i)
    {
        QDir cur_dir = dir_list[i];
        QStringList new_list = cur_dir.entryList(QStringList(""),QDir::AllDirs|QDir::NoDotAndDotDot);
        for(unsigned int index = 0;index < new_list.size();++index)
            dir_list << cur_dir.absolutePath() + "/" + new_list[index];
        QStringList file_list = cur_dir.entryList(QStringList(filter),QDir::Files|QDir::NoSymLinks);
        for (unsigned int index = 0;index < file_list.size();++index)
            src_list << dir_list[i] + "/" + file_list[index];
    }
    return src_list;
}

std::string find_full_path(QString name)
{
    QString filename = QCoreApplication::applicationDirPath() + name;
    if(QFileInfo(filename).exists())
        return filename.toStdString();
    filename = QDir::currentPath() + name;
    if(QFileInfo(filename).exists())
        return filename.toStdString();
    return std::string();
}

bool load_file_name(void)
{
    fib_template_file_name_2mm = find_full_path("/hcp1065_2mm.fib.gz");
    device_content_file = find_full_path("/device.txt");

    // search for all template
    {
        QDir dir = QCoreApplication::applicationDirPath()+ "/atlas";
        if(!dir.exists())
            dir = QDir::currentPath()+ "/atlas";

        QStringList name_list = dir.entryList(QStringList("*"),QDir::Dirs|QDir::NoSymLinks);

        // Make ICBM152 the default
        for(int i = 0;i < name_list.size();++i)
        {
            if(name_list[i].contains("ICBM"))
            {
                QString item_to_move = name_list[i];
                name_list.erase(name_list.begin()+i);
                name_list.insert(name_list.begin(),item_to_move);
            }
        }
        for(int i = 0;i < name_list.size();++i)
        {
            QDir template_dir = dir.absolutePath() + "/" + name_list[i];
            QString qa_file_path = template_dir.absolutePath() + "/" + name_list[i] + ".QA.nii.gz";
            QString iso_file_path = template_dir.absolutePath() + "/" + name_list[i] + ".ISO.nii.gz";
            QString tt_file_path = template_dir.absolutePath() + "/" + name_list[i] + ".tt.gz";
            if(!QFileInfo(qa_file_path).exists())
                continue;
            // setup QA and ISO template
            fa_template_list.push_back(qa_file_path.toStdString());
            if(QFileInfo(iso_file_path).exists())
                iso_template_list.push_back(iso_file_path.toStdString());
            else
                iso_template_list.push_back(std::string());

            if(QFileInfo(iso_file_path).exists())
                track_atlas_file_list.push_back(tt_file_path.toStdString());
            else
                track_atlas_file_list.push_back(std::string());
            // find related atlases
            {
                QStringList atlas_list = template_dir.entryList(QStringList("*.nii"),QDir::Files|QDir::NoSymLinks);
                atlas_list << template_dir.entryList(QStringList("*.nii.gz"),QDir::Files|QDir::NoSymLinks);
                atlas_list.sort();
                std::vector<std::string> atlas_file_list;
                for(int index = 0;index < atlas_list.size();++index)
                    if(QFileInfo(atlas_list[index]).baseName() != name_list[i])
                        atlas_file_list.push_back((template_dir.absolutePath() + "/" + atlas_list[index]).toStdString());
                template_atlas_list.push_back(std::move(atlas_file_list));
            }
        }
        if(fa_template_list.empty())
            return false;
    }
    return true;
}

void init_application(void)
{
    QApplication::setOrganizationName("LabSolver");
    QApplication::setApplicationName("DSI Studio");

    #ifdef __APPLE__
    QFont font;
    font.setFamily(QString::fromUtf8("Arial"));
    QApplication::setFont(font);
    #endif
    QSettings settings;
    QString style = settings.value("styles","Fusion").toString();
    if(style != "default" && !style.isEmpty())
        QApplication::setStyle(style);

    if(!load_file_name())
        QMessageBox::information(nullptr,"Error",
        "Cannot find FA template in the template folder. Please download dsi_studio_other_files.zip from DSI Studio website and place them with the DSI Studio executives.");
}

program_option po;
int run_action(std::shared_ptr<QApplication> gui)
{
    std::string action = po.get("action");
    if(action == std::string("rec"))
        return rec();
    if(action == std::string("trk"))
        return trk();
    if(action == std::string("atk"))
        return atk();
    if(action == std::string("src"))
        return src();
    if(action == std::string("ana"))
        return ana();
    if(action == std::string("exp"))
        return exp();
    if(action == std::string("atl"))
        return atl();
    if(action == std::string("cnt"))
        return cnt();
    if(action == std::string("cnt_ind"))
        return cnt_ind();
    if(action == std::string("ren"))
        return ren();
    if(action == std::string("cnn"))
        return cnn();
    if(action == std::string("qc"))
        return qc();
    if(action == std::string("reg"))
        return reg();
    if(action == std::string("vis"))
    {
        vis();
        if(po.get("stay_open") == std::string("1"))
            gui->exec();
        return 0;
    }
    std::cout << "Unknown action:" << action << std::endl;
    return 1;
}
void get_filenames_from(const std::string param,std::vector<std::string>& filenames);
int run_cmd(int ac, char *av[])
{
    try
    {
        std::cout << "DSI Studio " << __DATE__ << ", Fang-Cheng Yeh" << std::endl;
        if(!po.parse(ac,av))
        {
            std::cout << po.error_msg << std::endl;
            return 1;
        }
        std::shared_ptr<QApplication> gui;
        std::shared_ptr<QCoreApplication> cmd;
        for (int i = 1; i < ac; ++i)
            if ((std::string(av[i]) == std::string("--action=cnt") && po.get("no_tractogram",1) != 1) ||
                std::string(av[i]) == std::string("--action=vis"))
            {
                gui.reset(new QApplication(ac, av));
                init_application();
                std::cout << "Starting GUI-based command line interface." << std::endl;
                break;
            }

        if(!gui.get())
        {
            cmd.reset(new QCoreApplication(ac, av));
            if(!load_file_name())
            {
                std::cout << "Cannot find FA template in the template folder. Please download dsi_studio_other_files.zip from DSI Studio website and place them with the DSI Studio executives." << std::endl;
                return 1;
            }
            cmd->setOrganizationName("LabSolver");
            cmd->setApplicationName("DSI Studio");
        }
        if (!po.has("action"))
        {
            std::cout << "invalid command, use --help for more detail" << std::endl;
            return 1;
        }
        std::string source = po.get("source");
        if(po.get("action") != std::string("atk") && // atk handle * by itself
           (source.find('*') != std::string::npos ||
            source.find(',') != std::string::npos))
        {
            std::vector<std::string> source_files;
            get_filenames_from("source",source_files);
            for (size_t i = 0;i < source_files.size();++i)
            {
                std::cout << "Process file:" << source_files[i] << std::endl;
                po.set("source",source_files[i]);
                po.set_used(0);
                if(run_action(gui) == 1)
                {
                    std::cout << "Terminated due to error." << std::endl;
                    return 1;
                }
            }
        }
        else
            return run_action(gui);
    }
    catch(const std::exception& e ) {
        std::cout << e.what() << std::endl;
    }
    catch(...)
    {
        std::cout << "unknown error occured" << std::endl;
    }
    return 0;
}



int main(int ac, char *av[])
{
    if(ac > 2)
        return run_cmd(ac,av);
    if(QString(av[1]).endsWith(".txt") || QString(av[1]).endsWith(".log"))
        return run_cmd(ac,av);
    if(ac == 2)
        arg_file_name = av[1];


    QApplication a(ac,av);
    init_application();
    MainWindow w;

    has_gui = true;

    // presentation mode
    QStringList fib_list = QDir(QCoreApplication::applicationDirPath()+ "/presentation").
                            entryList(QStringList("*fib.gz") << QString("*_qa.nii.gz"),QDir::Files|QDir::NoSymLinks);
    if(fib_list.size())
    {
        w.hide();
        w.loadFib(QCoreApplication::applicationDirPath() + "/presentation/" + fib_list[0],true);
    }
    else
    {
        w.show();
        w.setWindowTitle(QString("DSI Studio ") + __DATE__ + " build");
    }
    return a.exec();
}
