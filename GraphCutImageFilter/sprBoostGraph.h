#ifndef BOOSTGRAPH_H
#define BOOSTGRAPH_H

#include <boost/config.hpp>
#include <iostream>
#include <string>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/read_dimacs.hpp>
#include <boost/graph/graph_utility.hpp>

//#include <boost\config.hpp>
//#include <iostream>
//#include <string>
//#include <boost\graph\adjacency_list.hpp>
//#include <boost\graph\boykov_kolmogorov_max_flow.hpp>
//#include <boost\graph\read_dimacs.hpp>
//#include <boost\graph\graph_utility.hpp>

using namespace boost;

//typedef adjacency_list_traits < vecS, vecS, directedS > Traits;
typedef adjacency_list_traits < setS, vecS, directedS > Traits;
//typedef adjacency_list < vecS, vecS, directedS,
typedef adjacency_list < setS, vecS, directedS,
  property < vertex_name_t, std::string,
  property < vertex_index_t, long,
  property < vertex_color_t, boost::default_color_type,
  property < vertex_distance_t, long,
  property < vertex_predecessor_t, Traits::edge_descriptor > > > > >,

  property < edge_capacity_t, float,
  property < edge_residual_capacity_t, float,
  property < edge_reverse_t, Traits::edge_descriptor > > > > Graph;

//typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
//typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
//typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

typedef  graph_traits<Graph>::vertices_size_type vertices_size_type;
typedef  graph_traits<Graph>::vertex_descriptor vertex_descriptor;
typedef  graph_traits<Graph>::edge_descriptor edge_descriptor;


#endif // BOOSTGRAPH_H
