#pragma once
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
   
    void fill_regular_tri(const Eigen::MatrixXd &G, bool delaunay=false){

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

        // create list with info and insert at once (CGAL can sort and optimize)
        std::vector<std::pair<Weighted_point, int>> P_with_info; 
        for (Weighted_point p : P) P_with_info.push_back(std::pair(p,m_maxind++));
        
        m_T.insert(P_with_info.begin(),P_with_info.end());

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


    void triangulate_polygon_ifs(const std::vector<Point> &dual_vertices, const std::vector<std::vector<int>> &dual_faces, Eigen::MatrixXd &V, Eigen::MatrixXi &F){
        
        V.resize(dual_vertices.size()+dual_faces.size(),3);
        int n_dual_tris = 0;
        for (auto f: dual_faces) n_dual_tris += f.size();
        F.resize(n_dual_tris,3);

        int i=0;
        for (; i<dual_vertices.size(); i++) {
            V.row(i) = Eigen::RowVector3d (dual_vertices[i].x(), dual_vertices[i].y(), dual_vertices[i].z());
        }
        
        int fi=0;
        for (auto f: dual_faces){
            int vc=0;
            int vn;
            Eigen::RowVector3d center_vertex = Eigen::RowVector3d::Zero();
            for (; vc<f.size(); vc++){
                center_vertex += V.row(f[vc]);
                vn = (vc+1)%f.size(); 
                F.row(fi++) = Eigen::RowVector3i(i,f[vc],f[vn]);
            }
            center_vertex /= f.size();
            V.row(i) = center_vertex;
            i++;
        }
    }

    void extract_dual_contour_ifs(std::vector<Point> &dual_vertices, std::vector<std::vector<int>> &dual_faces){

        typedef Tds::Cell_handle    Cell_handle;
        typedef Tds::Vertex_handle  Vertex_handle;

        std::vector<std::vector<Cell_handle>> df_cells;

        // add contour cell handles of an edge
        auto add_edge = [&](Tds::Edge e, bool reverse){

                std::vector<Cell_handle> facet;
                bool valid = true;
                // traverse adjacent cells and add the vertices to the facet (if they are finite)
                Tds::Cell_circulator cc = m_T.incident_cells(e);
                Tds::Cell_circulator strt = cc;

                std::function<Tds::Cell_circulator()> inc = [&](){return ++cc;};
                if (reverse) {
                    inc = [&](){return --cc;};
                }

                do {
                    if (m_T.is_infinite(cc)){
                        valid=false;
                        break;
                    } 
                    Cell_handle ch = cc;
                    facet.push_back(ch);
                    /*
                    Point p = T.dual(cc);
                    facet.push_back(p);
                    */
                } while (inc() != strt);
                
                df_cells.push_back(facet);
        };
        
        // add contour cell handle for all edges
        for (auto e: m_T.finite_edges()){
            Tds::Cell_handle ch = e.first;
            int vi_cell = e.second;
            int vj_cell = e.third;
            int vi = ch->vertex(vi_cell)->info();
            int vj = ch->vertex(vj_cell)->info();
           
            if ((m_sdf_values[vi] < 0.) && (m_sdf_values[vj] >= 0.)){
                add_edge(e, false);
            } else if ((m_sdf_values[vi] >= 0.) && (m_sdf_values[vj] <  0.)) {
                add_edge(e, true);
            }

        }

        // transform contour cell handles to IFS
        std::unordered_map<Cell_handle, int> cell_map;
        dual_vertices.clear();
        dual_faces.clear();

        for (auto fc_handles : df_cells){
            std::vector<int> facet;
            for (auto ch :fc_handles) {
                if (cell_map.find(ch) == cell_map.end()){
                    // create dual vertex if it doesnt exist yer
                    dual_vertices.push_back(m_T.dual(ch));
                    cell_map[ch] = dual_vertices.size()-1;
                }
                facet.push_back(cell_map[ch]);
            }
            dual_faces.push_back(facet);
        }
    }

    void extract_dual_contour_ifs_triangulated(Eigen::MatrixXd &V, Eigen::MatrixXi &F){

        std::vector<Point>              dual_vertices;
        std::vector<std::vector<int>>   dual_faces;
        fmt::print("call contour ifs\n");
        extract_dual_contour_ifs(dual_vertices,dual_faces);
        fmt::print("call triangulate polygon ifs\n");
        triangulate_polygon_ifs(dual_vertices,dual_faces, V, F);
    }

};
