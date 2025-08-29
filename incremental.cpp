#include <fmt/core.h>
#include "igl/read_triangle_mesh.h"
#include "igl/writeOBJ.h"
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
        fmt::print("usage: {} filepath (N) (max_refinement) (delaunay) (outpath)\n", argv[0]);
        return 0;
    }
    if (argc >=3) {
        N = atoi(argv[2]);
    }
    if (argc >= 4) {
        max_refinement = atoi(argv[3]);
    }

    if (argc >= 5) {
        is_delaunay = atoi(argv[4]);
    }

    if (argc >= 6) {
        outpath = argv[5];
        write=true;
        display=false;
    }
    fmt::print("Start Incremental\n");
    fmt::print("    N              = {}\n",N);
    fmt::print("    max_refinement = {}\n",max_refinement);
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

        int nr=0;
        if (max_refinement > 0) {
            fmt::print("Start {} Refinement\n", (IC.m_is_delaunay)? "DELAUNAY": "REGULAR");
            nr = IC.refine(max_refinement,sdf);
        } else {
            fmt::print("No Refinement");
        }

        // Primal Contouring:
        // extract current triangulation and SDF values
        Eigen::MatrixXd Vp;
        Eigen::MatrixXi Tp;
        Eigen::VectorXd Sp;
        IC.to_ifs(Vp,Tp,Sp);

        // perform marching tets
        Eigen::MatrixXd V_mt;
        Eigen::MatrixXi F_mt;
        igl::marching_tets(Vp,Tp, Sp, 0.0, V_mt, F_mt);

        // Dual Contouring:
        Eigen::MatrixXd Vd;
        Eigen::MatrixXi Fd;
        IC.extract_dual_contour_ifs_triangulated(Vd, Fd);

        if (display) {
            polyscope::init();

            auto s = polyscope::registerPointCloud("Samples", Vp);
            Eigen::VectorXd R = Sp;
            Eigen::VectorXd signs(R.size());
            for (int i=0; i<R.size(); i++){
                signs(i) = (R(i)>=0.)?1. : -1.;
                R(i) = fabs(R(i));
            }
            s->addScalarQuantity("radii", R);
            s->setPointRadiusQuantity("radii",false);
            s->addScalarQuantity("Sign", signs)->setColorMap("coolwarm")->setEnabled(true);
            s->setEnabled(false);
            polyscope::SlicePlane* psPlane = polyscope::addSceneSlicePlane();
            psPlane->setDrawPlane(false);
            psPlane->setDrawWidget(true);

            auto mpc_ = polyscope::registerSurfaceMesh("Primal Contouring", V_mt, F_mt);
            auto mdc_ = polyscope::registerSurfaceMesh("Dual   Contouring", Vd, Fd);
            mpc_->setIgnoreSlicePlane(psPlane->name, true);
            mdc_->setIgnoreSlicePlane(psPlane->name, true);

            polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
            polyscope::show();
        }
        if (write) {
            fmt::print("Write files to folder {}\n", outpath);
            // fmt::print("filename is : {}\n", filename);
            // (D)elaunay (Refinement) / (R)egular (R)efinement
            std::string filename = std::string(fs::path(argv[1]).filename());
            std::string method_string((is_delaunay)? "DTMTR":"RegMTR");
            std::string filepath = fmt::format("{}/{}_{}_{}_{}.obj",outpath,filename,method_string,N,nr);
            igl::writeOBJ(filepath, V_mt, F_mt);
        }
    }

    return 1;
}
