#include "IRT.hpp"
#include <fmt/core.h>

IncrementalContouring::IncrementalContouring(const Eigen::MatrixXd &G, bool is_delaunay):m_is_delaunay(is_delaunay){
    fill_regular_tri(G,m_is_delaunay);
    if (m_is_delaunay) {
        calculate_dual_vertices();
    } else {
        add_delaunay_refinement_candidates();
    }
}
   
void IncrementalContouring::fill_regular_tri(const Eigen::MatrixXd &G, bool delaunay){
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

Weighted_point IncrementalContouring::dual_weighted_point(Cell_handle ch){
        Point ps = m_T.dual(ch);

        // calculate weight of the dual tet as the average dist. to the vertices
        double pd = 0.;
        for (int vi=0; vi<4; vi++){
            pd += (ps - ch->vertex(vi)->point().point()).squared_length() - ch->vertex(vi)->point().weight();
        }
        pd /= 4.;

        return Weighted_point(ps,pd);
}

bool IncrementalContouring::tet_contains_sign_flip(Cell_handle fh){
    int vis[4] = {fh->vertex(0)->info(),fh->vertex(1)->info(),fh->vertex(2)->info(),fh->vertex(3)->info()};
    bool first_sign = m_sdf_values[vis[0]] >= 0.;
    for (int i=1; i<4; i++){
        if (first_sign != (m_sdf_values[vis[i]] >= 0.) ){
            return true;
        }
    }
    return false;
}

void IncrementalContouring::add_delaunay_refinement_candidates(){
    // add circumcenter of all tets that contain contour
    for (auto fh:m_T.finite_cell_handles()){
        if (tet_contains_sign_flip(fh)) {
            Weighted_point pw(m_T.dual(fh), 0.);
            m_refinement_candidates.push(Candidate(fh->vertex(0),fh->vertex(1),fh->vertex(2),fh->vertex(3),pw.point(),pw.weight()));
        }
    }
}

void IncrementalContouring::calculate_dual_vertices(){
    for (auto fh:m_T.finite_cell_handles()){
        Weighted_point pw = dual_weighted_point(fh);
        m_refinement_candidates.push(Candidate(fh->vertex(0),fh->vertex(1),fh->vertex(2),fh->vertex(3),pw.point(),pw.weight()));
    }
}

int IncrementalContouring::refine(int steps, std::function<double(double,double,double)> sdf, bool debug){

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

void IncrementalContouring::face_set(Eigen::MatrixXi &F){
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

void IncrementalContouring::to_ifs(Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::VectorXd &sdfvalues){
    face_set(T);
    V.resize(m_maxind,3);
    sdfvalues.resize(m_maxind);
    for (int i=0; i<m_maxind; i++) {
        V.row(i)    = Eigen::RowVector3d(m_positions[i].x(),m_positions[i].y(),m_positions[i].z());
        sdfvalues(i) = m_sdf_values[i];
    }
}
