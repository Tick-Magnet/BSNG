#ifndef __KDTREE2_HPP
#define __KDTREE2_HPP

#include "utils.h"
#include <vector>

typedef struct {
    float lower, upper;
} interval;

class kdtree2_node;
class searchrecord;

struct kdtree2_result {
    float dis;  // Square Euclidean distance
    int idx;    // Index of the neighbor
};

class kdtree2_result_vector : public std::vector<kdtree2_result> {
public:
    void push_element_and_heapify(kdtree2_result&);
    float replace_maxpri_elt_return_new_maxpri(kdtree2_result&);
    float max_value();
};

class kdtree2 {
public:
    const array2dfloat& the_data;
    const int N;
    int dim;
    bool sort_results;
    const bool rearrange;

    kdtree2(array2dfloat& data_in, bool rearrange_in = true, int dim_in = -1);
    ~kdtree2();

    void n_nearest_brute_force(std::vector<float>& qv, int nn, kdtree2_result_vector& result);
    void n_nearest(std::vector<float>& qv, int nn, kdtree2_result_vector& result);
    void n_nearest_around_point(int idxin, int correltime, int nn, kdtree2_result_vector& result);
    void r_nearest(std::vector<float>& qv, float r2, kdtree2_result_vector& result);
    void r_nearest_around_point(int idxin, int correltime, float r2, kdtree2_result_vector& result);
    int r_count(std::vector<float>& qv, float r2);
    int r_count_around_point(int idxin, int correltime, float r2);
    std::vector<int>* getIndex() { return &ind; }

    friend class kdtree2_node;
    friend class searchrecord;
private:
    kdtree2_node* root;
    const array2dfloat* data;
    std::vector<int> ind;
    array2dfloat rearranged_data;
    static const int bucketsize = 12;

    void set_data(array2dfloat& din);
    void build_tree();
    kdtree2_node* build_tree_for_range(int l, int u, kdtree2_node* parent);
    void select_on_coordinate(int c, int k, int l, int u);
    int select_on_coordinate_value(int c, float alpha, int l, int u);
    void spread_in_coordinate(int c, int l, int u, interval& interv);
};

class kdtree2_node {
public:
    kdtree2_node(int dim);
    ~kdtree2_node();

private:
    friend class kdtree2;

    int cut_dim;
    float cut_val, cut_val_left, cut_val_right;
    int l, u;
    std::vector<interval> box;
    kdtree2_node* left;
    kdtree2_node* right;

    void search(searchrecord& sr);
    bool box_in_search_range(searchrecord& sr);
    void check_query_in_bound(searchrecord& sr);
    void process_terminal_node(searchrecord& sr);
    void process_terminal_node_fixedball(searchrecord& sr);
};

#endif
