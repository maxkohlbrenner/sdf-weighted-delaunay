#pragma once

#include <queue>

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

    // functions
    public:
    IncrementalContouring(const Eigen::MatrixXd &G, bool is_delaunay=false);
    int refine(int steps, std::function<double(double,double,double)> sdf, bool debug=false);
    void to_ifs(Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::VectorXd &sdfvalues);

    private:
    void fill_regular_tri(const Eigen::MatrixXd &G, bool delaunay=false);
    Weighted_point dual_weighted_point(Cell_handle ch);
    bool tet_contains_sign_flip(Cell_handle fh);
    void add_delaunay_refinement_candidates();
    void calculate_dual_vertices();
    void face_set(Eigen::MatrixXi &F);
};
