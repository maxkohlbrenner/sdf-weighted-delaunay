#include <fmt/core.h>
#include "igl/read_triangle_mesh.h"
#include "igl/signed_distance.h"
#include "igl/marching_tets.h"

#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"

#include "src/IRT.hpp"

#include <filesystem>
namespace fs = std::filesystem;

void sample_sdf(int res, Eigen::MatrixXd &G, std::function<double(double,double,double)> sdf){
    G.resize(res*res*res, 4);
    int nx=res;
    int ny=res;
    int nz=res;
    for (int z=0 ; z<nz ; z++) {
        for (int y=0 ; y<ny ; y++) {
            for (int x=0 ; x<nx ; x++) {
                double x_ = (1./(nx-1) * x - 0.5)*2;
                double y_ = (1./(ny-1) * y - 0.5)*2;
                double z_ = (1./(nz-1) * z - 0.5)*2;
                G.row(z*nx*ny+y*nx+x) = Eigen::RowVector4d(x_,y_,z_,sdf(x_,y_,z_));
            }
        }
    }
}

int main(int argc, char *argv[])
{

    std::string outpath;
    int N = 10;
    int max_refinement = 1000;
    int show_steps = 100;
    bool display=true;
    bool write=false;
    bool is_delaunay=false;

    if(argc<2) {
        fmt::print("usage: {} path_to_obj (N) (max_refinement) (show_steps) (outpath)\n", argv[0]);
        return 0;
    }
    if (argc >=3) {
        N = atoi(argv[2]);
    }
    if (argc >= 4) {
        max_refinement = atoi(argv[3]);
    }
    if (argc >= 5) {
        show_steps = atoi(argv[4]);
    }

    if (argc >= 6) {
        is_delaunay = atoi(argv[5]);
    }

    if (argc >= 7) {
        outpath = argv[6];
        write=true;
        display=false;
    }
    fmt::print("Start Incremental\n");
    fmt::print("    N              = {}\n",N);
    fmt::print("    max_refinement = {}\n",max_refinement);
    fmt::print("    show_steps     = {}\n",show_steps);
    if (write) {
    fmt::print("    outpath      = {}\n",outpath);
    }

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    if(igl::readOBJ(argv[1], V,F)) {

        // in the beginning for them to be available later
        Eigen::MatrixXd P;
        igl::AABB<Eigen::MatrixXd,3> tree;
        igl::FastWindingNumberBVH fwn_bvh;

        V.rowwise() -= V.colwise().minCoeff();
        V /= V.maxCoeff();
        V.rowwise() -= 0.5 * V.colwise().maxCoeff();
        V *= 0.75/0.5;

        fmt::print("V.shape: ({},{})\n",V.rows(),V.cols());
        fmt::print("V min/max: ({},{})\n",V.minCoeff(),V.maxCoeff());

        tree.init(V,F);
        igl::fast_winding_number(V.cast<float>(), F, 2, fwn_bvh);
        auto sdf = [&](double x, double y, double z) {
            Eigen::VectorXd S(1);
            P = Eigen::RowVector3d(x,y,z);
            igl::signed_distance_fast_winding_number(P,V,F,tree,fwn_bvh,S);
            return S(0);
        };
        
        Eigen::MatrixXd S;
        sample_sdf(N,S,sdf);
        fmt::print("S.shape: ({},{})\n",S.rows(),S.cols());

        IncrementalContouring IC(S,is_delaunay);
      
        std::vector<int> Rs;
        std::vector<Eigen::MatrixXd> Vs_mt;
        std::vector<Eigen::MatrixXi> Fs_mt;

        std::vector<Eigen::MatrixXd> Vp_I;
        std::vector<Eigen::VectorXd> Sp_I;

        auto log_recs = [&](){

            // extract current triangulation and SDF values
            Eigen::MatrixXd Vp;
            Eigen::MatrixXi Tp;
            Eigen::VectorXd Sp;
            IC.to_ifs(Vp,Tp,Sp);

            // log positions and sdf values 
            Vp_I.push_back(Vp);
            Sp_I.push_back(Sp);

            // perform marching tets
            Eigen::MatrixXd V_mt;
            Eigen::MatrixXi F_mt;
            igl::marching_tets(Vp,Tp, Sp, 0.0, V_mt, F_mt);

            // log MT results and resolution
            Vs_mt.push_back(V_mt);
            Fs_mt.push_back(F_mt);
            Rs.push_back(IC.m_maxind);
        };

        log_recs();

        fmt::print("Start {} Refinement\n", (IC.m_is_delaunay)? "DELAUNAY": "REGULAR");

        int rsteps = 0;
        while (rsteps < max_refinement) {
            int nr = IC.refine(show_steps,sdf);
            fmt::print("performed {} refinement steps\n", nr);
            log_recs();

            rsteps += nr;
            if (nr < show_steps) break;
        }

        if (display) {
            polyscope::init();
            polyscope::registerSurfaceMesh("Ground Truth", V, F);
            auto s = polyscope::registerPointCloud("Samples", S.block(0,0,S.rows(),3));
            Eigen::VectorXd R = S.col(3);
            Eigen::VectorXd signs(R.size());
            for (int i=0; i<R.size(); i++){
                signs(i) = (R(i)>=0.)?1. : -1.;
                R(i) = fabs(R(i));
            }
            s->addScalarQuantity("radii", R);
            s->setPointRadiusQuantity("radii",false);
            s->addScalarQuantity("Sign", signs)->setEnabled(true);
            s->setEnabled(false);

            for (int i=0; i<Vs_mt.size(); i++) {
                auto rec = polyscope::registerSurfaceMesh(fmt::format("MT Rec {} ({} samples)",i,Rs[i]), Vs_mt[i], Fs_mt[i]);
                if (i!= Vs_mt.size()-1) rec->setEnabled(false); 
            }

            polyscope::SlicePlane* psPlane = polyscope::addSceneSlicePlane();
            psPlane->setDrawPlane(false);
            psPlane->setDrawWidget(true);
            polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
            polyscope::show();
        } 
        if (write) {
            fmt::print("Write files to folder {}\n", outpath);
            std::string filename = std::string(fs::path(argv[1]).filename());
            fmt::print("filename is : {}\n", filename);

            // (D)elaunay (Refinement) / (R)egular (R)efinement
            std::string method_string((is_delaunay)? "DTMTR":"RegMTR");

            const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
            for (int i=0; i<Vs_mt.size(); i++) {

                std::ofstream vertfile(fmt::format("{}/{}_{}_{}_{}_{}_V.csv",outpath,filename, method_string,N,i,Rs[i]));
                vertfile << Vs_mt[i].format(CSVFormat);
                vertfile.close();

                std::ofstream trifile(fmt::format("{}/{}_{}_{}_{}_{}_F.csv",outpath,filename,method_string,N,i,Rs[i]));
                trifile << Fs_mt[i].format(CSVFormat);
                trifile.close();

                // log sample positions
                std::ofstream posfile(fmt::format("{}/{}_{}_{}_{}_{}_P.csv",outpath,filename,method_string,N,i,Rs[i]));
                posfile << Vp_I[i].format(CSVFormat);
                posfile.close();

                // log sample positions
                std::ofstream valfile(fmt::format("{}/{}_{}_{}_{}_{}_S.csv",outpath,filename,method_string,N,i,Rs[i]));
                valfile << Sp_I[i].format(CSVFormat);
                valfile.close();

            }

            for (int i=0; i<Rs.size(); i++) {

                // write out a regular sdf grid in an appropriate resolution
                Eigen::MatrixXd S;
                int N = std::ceil(std::pow(Rs[i],1./3.));
                sample_sdf(N,S,sdf);
                std::ofstream sdfgridfile(fmt::format("{}/{}_{}_sdfgrid.csv",outpath,filename,N));
                sdfgridfile <<S.format(CSVFormat);
                sdfgridfile.close();
            }
        }

    } 

    return 1;
}
