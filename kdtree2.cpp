#include "kdtree2.hpp"

// Utility function to compute the square of a float.
inline float squared(const float x) {
  return (x * x);
}

// Utility function to swap two integers.
inline void swap(int& a, int& b) {
  int tmp;
  tmp = a;
  a = b;
  b = tmp;
}

// Utility function to swap two floats.
inline void swap(float& a, float& b) {
  float tmp;
  tmp = a;
  a = b;
  b = tmp;
}

// Implementation of operator< for kdtree2_result comparison.
inline bool operator<(const kdtree2_result& e1, const kdtree2_result& e2) {
  return (e1.dis < e2.dis);
}

// Implementation of KDTREE2_RESULT_VECTOR member functions.
float kdtree2_result_vector::max_value() {
  return (*begin()).dis; // Return the distance of the very first element.
}

void kdtree2_result_vector::push_element_and_heapify(kdtree2_result& e) {
  push_back(e);              // Add the element to the vector.
  push_heap(begin(), end()); // Heapify the vector with the new element.
}

float kdtree2_result_vector::replace_maxpri_elt_return_new_maxpri(kdtree2_result& e) {
  // Remove the maximum priority element from the queue, replace it with 'e', and return its priority.
  pop_heap(begin(), end());
  pop_back();
  push_back(e);              // Insert the new element.
  push_heap(begin(), end()); // Heapify the vector.
  return (*this)[0].dis;     // Return the new maximum priority element's distance.
}

// Implementation of KDTREE2 constructor.
kdtree2::kdtree2(array2dfloat& data_in, bool rearrange_in, int dim_in)
    : the_data(data_in),
      N(data_in.size()),
      dim(data_in[0].size()),
      sort_results(false),
      rearrange(rearrange_in),
      root(NULL),
      data(NULL),
      ind(N) {
  #ifdef _DEBUG
  not_regular_median = 0;
  regular_median = 0;
  #endif

  if (dim_in > 0)
    dim = dim_in;

  build_tree();
  
  if (rearrange) {
    rearranged_data.resize(N);
    for (int ll = 0; ll < N; ll++)
      rearranged_data[ll].resize(dim);

    // Permute the data for rearranged tree.
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < dim; j++) {
        rearranged_data[i][j] = the_data[ind[i]][j];
      }
    }
    data = &rearranged_data;
  } else {
    data = &the_data;
  }
}

// Implementation of KDTREE2 destructor.
kdtree2::~kdtree2() {
  delete root;
}

// Implementation of KDTREE2 building routines.
void kdtree2::build_tree() {
  for (int i = 0; i < N; i++)
    ind[i] = i;

  root = build_tree_for_range(0, N - 1, NULL);
}

kdtree2_node* kdtree2::build_tree_for_range(int l, int u, kdtree2_node* parent) {
  kdtree2_node* node = new kdtree2_node(dim);

  if (u < l) {
    return (NULL);
  }

  if ((u - l) <= bucketsize) {
    // Create a terminal node.
    for (int i = 0; i < dim; i++) {
      spread_in_coordinate(i, l, u, node->box[i]);
    }

    node->cut_dim = 0;
    node->cut_val = 0.0;
    node->l = l;
    node->u = u;
    node->left = node->right = NULL;
  } else {
    int c = -1;
    float maxspread = 0.0;
    int m;

    for (int i = 0; i < dim; i++) {
      if ((parent == NULL) || (parent->cut_dim == i)) {
        spread_in_coordinate(i, l, u, node->box[i]);
      } else {
        node->box[i] = parent->box[i];
      }

      float spread = node->box[i].upper - node->box[i].lower;

      if (spread > maxspread) {
        maxspread = spread;
        c = i;
      }
    }

    if (false) {
      m = (l + u) / 2;
      select_on_coordinate(c, m, l, u);
    } else {
      float sum;
      float average;

      if (true) {
        sum = 0.0;
        for (int k = l; k <= u; k++) {
          sum += the_data[ind[k]][c];
        }
        average = sum / static_cast<float>(u - l + 1);
      } else {
        average = (node->box[c].upper + node->box[c].lower) * 0.5;
      }

      m = select_on_coordinate_value(c, average, l, u);
    }

    if (m <= l || m >= u) {
      m = (l + u) / 2;
    }

    node->cut_dim = c;
    node->l = l;
    node->u = u;
    node->left = build_tree_for_range(l, m, node);
    node->right = build_tree_for_range(m + 1, u, node);

    if (node->right == NULL) {
      for (int i = 0; i < dim; i++)
        node->box[i] = node->left->box[i];
      node->cut_val = node->left->box[c].upper;
      node->cut_val_left = node->cut_val_right = node->cut_val;
    } else if (node->left == NULL) {
      for (int i = 0; i < dim; i++)
        node->box[i] = node->right->box[i];
      node->cut_val = node->right->box[c].upper;
      node->cut_val_left = node->cut_val_right = node->cut_val;
    } else {
      node->cut_val_right = node->right->box[c].lower;
      node->cut_val_left = node->left->box[c].upper;
      node->cut_val = (node->cut_val_left + node->cut_val_right) / 2.0;

      for (int i = 0; i < dim; i++) {
        node->box[i].upper = max(node->left->box[i].upper,
                                 node->right->box[i].upper);
        node->box[i].lower = min(node->left->box[i].lower,
                                 node->right->box[i].lower);
      }
    }
  }
  return (node);
}

void kdtree2::spread_in_coordinate(int c, int l, int u, interval& interv) {
  float smin, smax;
  float lmin, lmax;
  int i;

  smin = the_data[ind[l]][c];
  smax = smin;

  for (i = l + 2; i <= u; i += 2) {
    lmin = the_data[ind[i - 1]][c];
    lmax = the_data[ind[i]][c];

    if (lmin > lmax) {
      swap(lmin, lmax);
    }

    if (smin > lmin) smin = lmin;
    if (smax < lmax) smax = lmax;
  }

  if (i == u + 1) {
    float last = the_data[ind[u]][c];
    if (smin > last) smin = last;
    if (smax < last) smax = last;
  }
  interv.lower = smin;
  interv.upper = smax;
}

// Select elements within a specified range based on a coordinate value.
void kdtree2::select_on_coordinate(int c, int k, int l, int u) {
  while (l < u) {
    int t = ind[l];
    int m = l;

    for (int i = l + 1; i <= u; i++) {
      if (the_data[ind[i]][c] < the_data[t][c]) {
        m++;
        swap(ind[i], ind[m]);
      }
    }

    swap(ind[l], ind[m]);

    if (m <= k) l = m + 1;
    if (m >= k) u = m - 1;
  }
}

// Select elements within a specified range based on a coordinate value.
int kdtree2::select_on_coordinate_value(int c, float alpha, int l, int u) {
  int lb = l, ub = u;

  while (lb < ub) {
    if (the_data[ind[lb]][c] <= alpha) {
      lb++; // Good where it is.
    } else {
      swap(ind[lb], ind[ub]);
      ub--;
    }
  }

  if (the_data[ind[lb]][c] <= alpha)
    return (lb);
  else
    return (lb - 1);
}

// Search record class to store search-related information.
class searchrecord {
private:
  friend class kdtree2;
  friend class kdtree2_node;

  vector<float>& qv;
  int dim;
  bool rearrange;
  unsigned int nn;
  float ballsize;
  int centeridx, correltime;

  kdtree2_result_vector& result;
  const array2dfloat* data;
  const vector<int>& ind;

public:
  searchrecord(vector<float>& qv_in, kdtree2& tree_in, kdtree2_result_vector& result_in)
      : qv(qv_in),
        result(result_in),
        data(tree_in.data),
        ind(tree_in.ind) {
    dim = tree_in.dim;
    rearrange = tree_in.rearrange;
    static const float infinity = 1.0e38;
    ballsize = infinity; 
    nn = 0;
  };
};

// Brute-force search for N nearest neighbors.
void kdtree2::n_nearest_brute_force(vector<float>& qv, int nn, kdtree2_result_vector& result) {
  result.clear();

  for (int i = 0; i < N; i++) {
    float dis = 0.0;
    kdtree2_result e;
    for (int j = 0; j < dim; j++) {
      dis += squared(the_data[i][j] - qv[j]);
    }
    e.dis = dis;
    e.idx = i;
    result.push_back(e);
  }
  sort(result.begin(), result.end());
}

// Search for N nearest neighbors.
void kdtree2::n_nearest(vector<float>& qv, int nn, kdtree2_result_vector& result) {
  searchrecord sr(qv, *this, result);
  vector<float> vdiff(dim, 0.0);

  result.clear();

  sr.centeridx = -1;
  sr.correltime = 0;
  sr.nn = nn;

  root->search(sr);

  if (sort_results) sort(result.begin(), result.end());
}

// Search for N nearest neighbors around a point.
void kdtree2::n_nearest_around_point(int idxin, int correltime, int nn, kdtree2_result_vector& result) {
  vector<float> qv(dim);

  result.clear();

  for (int i = 0; i < dim; i++) {
    qv[i] = the_data[idxin][i];
  }

  {
    searchrecord sr(qv, *this, result);
    sr.centeridx = idxin;
    sr.correltime = correltime;
    sr.nn = nn;
    root->search(sr);
  }

  if (sort_results) sort(result.begin(), result.end());
}

// Search for elements within a specified radius.
void kdtree2::r_nearest(vector<float>& qv, float r2, kdtree2_result_vector& result) {
  searchrecord sr(qv, *this, result);
  vector<float> vdiff(dim, 0.0);

  result.clear();

  sr.centeridx = -1;
  sr.correltime = 0;
  sr.nn = 0;
  sr.ballsize = r2;

  root->search(sr);

  if (sort_results) sort(result.begin(), result.end());
}

// Count elements within a specified radius.
int kdtree2::r_count(vector<float>& qv, float r2) {
  {
    kdtree2_result_vector result;
    searchrecord sr(qv, *this, result);

    sr.centeridx = -1;
    sr.correltime = 0;
    sr.nn = 0;
    sr.ballsize = r2;

    root->search(sr);
    return (result.size());
  }
}

// Search for elements within a specified radius around a point.
void kdtree2::r_nearest_around_point(int idxin, int correltime, float r2, kdtree2_result_vector& result) {
  vector<float> qv(dim);

  result.clear();

  for (int i = 0; i < dim; i++) {
    qv[i] = the_data[idxin][i];
  }

  {
    searchrecord sr(qv, *this, result);
    sr.centeridx = idxin;
    sr.correltime = correltime;
    sr.ballsize = r2;
    sr.nn = 0;
    root->search(sr);
  }

  if (sort_results) sort(result.begin(), result.end());
}

// Count elements within a specified radius around a point.
int kdtree2::r_count_around_point(int idxin, int correltime, float r2) {
  vector<float> qv(dim);

  for (int i = 0; i < dim; i++) {
    qv[i] = the_data[idxin][i];
  }

  {
    kdtree2_result_vector result;
    searchrecord sr(qv, *this, result);

    sr.centeridx = idxin;
    sr.correltime = correltime;
    sr.ballsize = r2;
    sr.nn = 0;
    root->search(sr);
    return (result.size());
  }
}

// KDTREE2_NODE constructor.
kdtree2_node::kdtree2_node(int dim) : box(dim) {
  left = right = NULL;
}

// KDTREE2_NODE destructor.
kdtree2_node::~kdtree2_node() {
  if (left != NULL) delete left;
  if (right != NULL) delete right;
}

// Perform a search in the k-d tree.
void kdtree2_node::search(searchrecord& sr) {
  // The core search routine.
  // This uses true distance to bounding box as the
  // criterion to search the secondary node. 
  //
  // This results in somewhat fewer searches of the secondary nodes
  // than 'search', which uses the vdiff vector,  but as this
  // takes more computational time, the overall performance may not
  // be improved in actual run time. 
  //

  if ((left == NULL) && (right == NULL)) {
    // We are on a terminal node.
    if (sr.nn == 0) {
      process_terminal_node_fixedball(sr);
    } else {
      process_terminal_node(sr);
    }
  } else {
    kdtree2_node* ncloser, *nfarther;

    float extra;
    float qval = sr.qv[cut_dim]; 
    // Value of the wall boundary on the cut dimension. 
    if (qval < cut_val) {
      ncloser = left;
      nfarther = right;
      extra = cut_val_right - qval;
    } else {
      ncloser = right;
      nfarther = left;
      extra = qval - cut_val_left; 
    }

    if (ncloser != NULL) ncloser->search(sr);

    if ((nfarther != NULL) && (squared(extra) < sr.ballsize)) {
      // First cut
      if (nfarther->box_in_search_range(sr)) {
        nfarther->search(sr); 
      }      
    }
  }
}

// Calculate distance from a point to a boundary.
inline float dis_from_bnd(float x, float amin, float amax) {
  if (x > amax) {
    return (x - amax); 
  } else if (x < amin) {
    return (amin - x);
  } else {
    return 0.0;
  }
}

// Check if the bounding box is within the search range.
inline bool kdtree2_node::box_in_search_range(searchrecord& sr) {
  //
  // Does the bounding box, represented by minbox[*],maxbox[*]
  // have any point which is within 'sr.ballsize' to 'sr.qv'??
  //
 
  int dim = sr.dim;
  float dis2 = 0.0; 
  float ballsize = sr.ballsize; 
  for (int i = 0; i < dim; i++) {
    dis2 += squared(dis_from_bnd(sr.qv[i], box[i].lower, box[i].upper));
    if (dis2 > ballsize) {
      return false;
    }
  }
  return true;
}

// Process a terminal node in the k-d tree with variable-sized ball.
void kdtree2_node::process_terminal_node(searchrecord& sr) {
  int centeridx = sr.centeridx;
  int correltime = sr.correltime;
  unsigned int nn = sr.nn; 
  int dim = sr.dim;
  float ballsize = sr.ballsize;
  //
  bool rearrange = sr.rearrange; 
  const array2dfloat& data = *sr.data;

  const bool debug = false;

  if (debug) {
    printf("Processing terminal node %d, %d\n", l, u);
    cout << "Query vector = [";
    for (int i = 0; i < dim; i++) cout << sr.qv[i] << ','; 
    cout << "]\n";
    cout << "nn = " << nn << '\n';
    check_query_in_bound(sr);
  }

  for (int i = l; i <= u; i++) {
    int indexofi;  // sr.ind[i]; 
    float dis;
    bool early_exit; 

    if (rearrange) {
      early_exit = false;
      dis = 0.0;
      for (int k = 0; k < dim; k++) {
        dis += squared(data[i][k] - sr.qv[k]);
        if (dis > ballsize) {
          early_exit = true; 
          break;
        }
      }
      if (early_exit) continue; // Next iteration of the main loop.
      // Why do we do things like this? Because if we take an early
      // exit (due to distance being too large) which is common, then
      // we need not read in the actual point index, thus saving main
      // memory bandwidth. If the distance to the point is less than the
      // ballsize, though, then we need the index.
      //
      indexofi = sr.ind[i];
    } else {
      // 
      // But if we are not using the rearranged data, then
      // we must always 
      indexofi = sr.ind[i];
      early_exit = false;
      dis = 0.0;
      for (int k = 0; k < dim; k++) {
        dis += squared(data[indexofi][k] - sr.qv[k]);
        if (dis > ballsize) {
          early_exit = true; 
          break;
        }
      }
      if (early_exit) continue; // Next iteration of the main loop.
    } // End if rearrange. 
    
    if (centeridx > 0) {
      // We are doing decorrelation interval.
      if (abs(indexofi - centeridx) < correltime) continue; // Skip this point. 
    }

    // Here the point must be added to the list.
    //
    // Two choices for any point. The list so far is either
    // undersized, or it is not.
    //
    if (sr.result.size() < nn) {
      kdtree2_result e;
      e.idx = indexofi;
      e.dis = dis;
      sr.result.push_element_and_heapify(e); 
      if (debug) cout << "Unilaterally pushed dis=" << dis;
      if (sr.result.size() == nn) ballsize = sr.result.max_value();
      // Set the ball radius to the largest on the list (maximum priority).
      if (debug) {
        cout << " ballsize = " << ballsize << "\n"; 
        cout << "sr.result.size() = "  << sr.result.size() << '\n';
      }
    } else {
      //
      // If we get here then the current node has a squared 
      // distance smaller
      // than the last on the list, and belongs on the list.
      // 
      kdtree2_result e;
      e.idx = indexofi;
      e.dis = dis;
      ballsize = sr.result.replace_maxpri_elt_return_new_maxpri(e); 
      if (debug) {
        cout << "Replaced maximum dis with dis=" << dis << 
          " new ballsize =" << ballsize << '\n';
      }
    }
  } // Main loop
  sr.ballsize = ballsize;
}

// Process a terminal node in the k-d tree with a fixed-sized ball.
void kdtree2_node::process_terminal_node_fixedball(searchrecord& sr) {
  int centeridx = sr.centeridx;
  int correltime = sr.correltime;
  int dim = sr.dim;
  float ballsize = sr.ballsize;
  //
  bool rearrange = sr.rearrange; 
  const array2dfloat& data = *sr.data;

  for (int i = l; i <= u; i++) {
    int indexofi = sr.ind[i]; 
    float dis;
    bool early_exit; 

    if (rearrange) {
      early_exit = false;
      dis = 0.0;
      for (int k = 0; k < dim; k++) {
        dis += squared(data[i][k] - sr.qv[k]);
        if (dis > ballsize) {
          early_exit = true; 
          break;
        }
      }
      if (early_exit) continue; // Next iteration of main loop.
      // Why do we do things like this? Because if we take an early
      // exit (due to distance being too large) which is common, then
      // we need not read in the actual point index, thus saving main
      // memory bandwidth. If the distance to the point is less than the
      // ballsize, though, then we need the index.
      //
      indexofi = sr.ind[i];
    } else {
      // 
      // But if we are not using the rearranged data, then
      // we must always 
      indexofi = sr.ind[i];
      early_exit = false;
      dis = 0.0;
      for (int k = 0; k < dim; k++) {
        dis += squared(data[indexofi][k] - sr.qv[k]);
        if (dis > ballsize) {
          early_exit = true; 
          break;
        }
      }
      if (early_exit) continue; // Next iteration of main loop.
    } // End if rearrange. 
    
    if (centeridx > 0) {
      // We are doing decorrelation interval.
      if (abs(indexofi - centeridx) < correltime) continue; // Skip this point. 
    }

    {
      kdtree2_result e;
      e.idx = indexofi;
      e.dis = dis;
      sr.result.push_back(e);
    }
  }
}