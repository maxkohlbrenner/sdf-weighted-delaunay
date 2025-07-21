#include <fmt/core.h>
#include "igl/read_triangle_mesh.h"
#include "igl/signed_distance.h"
#include "igl/marching_tets.h"

#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::FT                                                       Weight;
typedef K::Point_3                                                  Point;
typedef K::Weighted_point_3                                         Weighted_point;
typedef CGAL::Regular_triangulation_vertex_base_3<K>                Vb0;
typedef CGAL::Triangulation_vertex_base_with_info_3<int, K, Vb0>    Vb;
typedef CGAL::Regular_triangulation_cell_base_3<K>                  Cb;
typedef CGAL::Triangulation_data_structure_3<Vb,Cb>                 Tds;
typedef CGAL::Regular_triangulation_3<K, Tds>                       Rt;
typedef Tds::Vertex_handle                                          Vertex_handle;
typedef Tds::Cell_handle                                            Cell_handle;

#include <filesystem>
namespace fs = std::filesystem;

// ----------------------- REFINEMENT ------------------------------
class Candidate {
    public:
    Vertex_handle m_v0;
    Vertex_handle m_v1;
    Vertex_handle m_v2;
    Vertex_handle m_v3;
    double m_score;
    Point m_point;

	Candidate(  Vertex_handle v0,
                Vertex_handle v1,
                Vertex_handle v2,
                Vertex_handle v3,
				Point p,
                double score ) : 
                m_v0{v0},m_v1{v1},m_v2{v2},m_v3{v3},m_point{p},m_score{score}
    {}
  // ~CCandidate() {} 
};

struct less_refinement
{ 
  bool operator()(const Candidate& c1, 
                  const Candidate& c2) const
  {
    return (c1.m_score < c2.m_score);
  }
};
typedef typename std::priority_queue<Candidate, std::vector<Candidate>,less_refinement > PQueue;

// ------------------------------------------------------------------

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


void get_RT_faces(const Rt &T, Eigen::MatrixXi &F){
    F.resize(T.number_of_finite_cells(),4);
    int i=0;
    for (auto ch: T.finite_cell_handles()){
        for (int vi=0; vi<4; vi++) {
            F(i,vi) = ch->vertex(vi)->info();
        }    
        i++;
    }
}


class IncrementalContouring {       
    public:             
    Rt m_T;
    PQueue m_refinement_candidates;
    bool m_is_delaunay;

    // used for IFS later
    int m_maxind=0;
    std::vector<Point>  m_positions;
    std::vector<double> m_sdf_values;

    IncrementalContouring(const Eigen::MatrixXd &G, bool is_delaunay=false):m_is_delaunay(is_delaunay){

        fill_regular_tri(G,m_is_delaunay);
        if (m_is_delaunay) {
            calculate_dual_vertices();
        } else {
            add_delaunay_refinement_candidates();
        }
    }
   
    void fill_regular_tri(const Eigen::MatrixXd &G, bool delaunay=false,bool pbyp=false){

        // Create Regular Tri
        m_T = Rt();
        m_maxind=0;

        std::vector<Weighted_point> P;
        int number_of_points = G.rows();
        for (int i=0; i<number_of_points; i++){
            Point p(G(i,0), G(i,1), G(i,2));
            double sdfv = G(i,3); 
            Weight w = (delaunay)? 0. :  (Weight) sdfv*sdfv;
            P.push_back(Weighted_point(p, w));

            m_positions.push_back(p);
            m_sdf_values.push_back(sdfv);
        }

        if (pbyp) {
            // insert vertices one after the other in order to add the info() with the index
            for (Weighted_point p : P){
                Rt::Vertex_handle vh = m_T.insert(p);
                vh->info() = m_maxind++;
            }
        } else {
            // create list with info and insert at once (CGAL can sort and optimize)
            std::vector<std::pair<Weighted_point, int>> P_with_info; 
            for (Weighted_point p : P) P_with_info.push_back(std::pair(p,m_maxind++));
            
            m_T.insert(P_with_info.begin(),P_with_info.end());
        }

        assert( m_T.is_valid() );
        assert( m_T.dimension() == 3 );
    }

    Weighted_point dual_weighted_point(Cell_handle ch){
            Point ps = m_T.dual(ch);

            // calculate weight of the dual tet as the average dist. to the vertices
            double pd = 0.;
            for (int vi=0; vi<4; vi++){
                pd += (ps - ch->vertex(vi)->point().point()).squared_length() - ch->vertex(vi)->point().weight();
            }
            pd /= 4.;

            return Weighted_point(ps,pd);
    }

    bool tet_contains_sign_flip(Cell_handle fh){
        int vis[4] = {fh->vertex(0)->info(),fh->vertex(1)->info(),fh->vertex(2)->info(),fh->vertex(3)->info()};
        bool first_sign = m_sdf_values[vis[0]] >= 0.;
        for (int i=1; i<4; i++){
            if (first_sign != (m_sdf_values[vis[i]] >= 0.) ){
                return true;
            }
        }
        return false;
    }

    void add_delaunay_refinement_candidates(){
        // add circumcenter of all tets that contain contour
        for (auto fh:m_T.finite_cell_handles()){
            if (tet_contains_sign_flip(fh)) {
                Weighted_point pw(m_T.dual(fh), 0.);
                m_refinement_candidates.push(Candidate(fh->vertex(0),fh->vertex(1),fh->vertex(2),fh->vertex(3),pw.point(),pw.weight()));
            }
        }
    }

    void calculate_dual_vertices(){
        for (auto fh:m_T.finite_cell_handles()){
            Weighted_point pw = dual_weighted_point(fh);
            m_refinement_candidates.push(Candidate(fh->vertex(0),fh->vertex(1),fh->vertex(2),fh->vertex(3),pw.point(),pw.weight()));
        }
    }

    int refine(int steps, std::function<double(double,double,double)> sdf, bool debug=false){

        int ninserts=0;
        while (!m_refinement_candidates.empty()){
            Candidate c = m_refinement_candidates.top();
            m_refinement_candidates.pop();

            if(debug) {
                fmt::print("    point: ({},{},{})\n", c.m_point.x(),c.m_point.y(),c.m_point.z());
                fmt::print("    weight: {}\n", c.m_score);
                fmt::print("    dual of cell ({},{},{},{})\n", c.m_v0->info(),c.m_v1->info(),c.m_v2->info(),c.m_v3->info());
            }

            bool in_range = (c.m_point-CGAL::ORIGIN).squared_length() < 2.;
            if (!in_range) {
                if (debug)
                    fmt::print("    ... SKIP! Point out of unit cube\n");
                continue;
            }

            Cell_handle cell = NULL;
            bool is_cell = m_T.is_cell(c.m_v0,c.m_v1,c.m_v2,c.m_v3,cell);
            if(is_cell)
            {
                // fmt::print("--> is cell\n");
                Point point         = c.m_point;
                double sdfvalue = sdf(point.x(),point.y(),point.z());
                Weighted_point pw(point, (m_is_delaunay)? 0. : sdfvalue*sdfvalue);
                Vertex_handle v = m_T.insert(pw,cell); 
                // fmt::print("    inserted!\n");
                v->info() = m_maxind++;
                m_positions.push_back(point);
                m_sdf_values.push_back(sdfvalue);

                // fmt::print("    info is: {}\n", v->info());
                if(++ninserts >= steps) return ninserts;

                std::list<Cell_handle> cells;
                m_T.incident_cells(v,std::back_inserter(cells));

                typename std::list<Cell_handle>::iterator it;
                for(it = cells.begin();
                        it != cells.end();
                        it++)
                {
                    Cell_handle c = *it;
                    if(m_T.is_infinite(c))
                        continue;
                    if (!m_is_delaunay || tet_contains_sign_flip(c)){
                        Weighted_point pw = (m_is_delaunay)? Weighted_point(m_T.dual(c), 0.) : dual_weighted_point(c);
                        m_refinement_candidates.push(Candidate(c->vertex(0),c->vertex(1),c->vertex(2),c->vertex(3),pw.point(),pw.weight()));
                    }
                }

            } else {
                if(debug) {
                    fmt::print("    (SKIP, cell no longer exists)\n");
                }
            }

            // fmt::print("...next one\n\n");

        }

        return ninserts;
    }

    void face_set(Eigen::MatrixXi &F){
        F.resize(m_T.number_of_finite_cells(),4);
        int ci=0;
        for (auto ch : m_T.finite_cell_handles()){
            F.row(ci++) = Eigen::RowVector4i(ch->vertex(0)->info(),
                                                ch->vertex(1)->info(),
                                                ch->vertex(2)->info(),
                                                ch->vertex(3)->info()
            );
        }
    }

    void to_ifs(Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::VectorXd &sdfvalues){
        face_set(T);
        V.resize(m_maxind,3);
        sdfvalues.resize(m_maxind);
        for (int i=0; i<m_maxind; i++) {
            V.row(i)    = Eigen::RowVector3d(m_positions[i].x(),m_positions[i].y(),m_positions[i].z());
            sdfvalues(i) = m_sdf_values[i];
        }
    }

};

int main(int argc, char *argv[])
{

    std::string outpath;
    int N = 10;
    int max_refinement = 1000;
    int show_steps = 100;
    bool render=true;
    bool write=false;
    bool is_delaunay=false;

    if(argc<2) {
        fmt::print("usage: {} (N) (max_refinement) (show_steps) (outpath)\n", argv[0]);
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
        render=false;
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

        std::vector<Point> dual_points;
        while (!IC.m_refinement_candidates.empty()){
            Candidate c = IC.m_refinement_candidates.top();
            IC.m_refinement_candidates.pop();
            dual_points.push_back(c.m_point);
        }

        if (render) {

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

            // polyscope::registerPointCloud("Dual Points", dual_points)->setEnabled(false);
            // auto pm = polyscope::registerTetMesh("Primal Mesh", IC.m_positions,Tp);
            // pm->setEnabled(false);

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
        // if (false) { // disable writing for timings
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
