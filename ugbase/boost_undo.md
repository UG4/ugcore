Dateien in `ugcore/ugbase`, die <boost/function.hpp> einbinden

- ugcore/ugbase/bindings/vrl/user_data.cpp
- ugcore/ugbase/lib_grid/grid/grid.h
- ugcore/ugbase/lib_grid/parallelization/parallelization_util.h
- ugcore/ugbase/lib_grid/algorithms/extrusion/expand_layers_arte3D.cpp
- ugcore/ugbase/lib_grid/algorithms/extrusion/expand_layers.cpp
- ugcore/ugbase/lib_grid/algorithms/extrusion/expand_layers_arte.cpp
- ugcore/ugbase/lib_grid/algorithms/extrusion/ArteExpandFracs3D.h
- ugcore/ugbase/lib_grid/algorithms/extrusion/ArteExpandFracs3D.cpp
- ugcore/ugbase/registry/registry.h
- ugcore/ugbase/lib_disc/function_spaces/grid_function_util.h
- ugcore/ugbase/lib_disc/function_spaces/integrate.h
- ugcore/ugbase/lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h
- ugcore/ugbase/lib_disc/spatial_disc/elem_disc/dirac_source/lagrange_dirac_source.h
- ugcore/ugbase/common/util/message_hub.h
- ugcore/ugbase/pcl/pcl_layout_tests.h

Weitere Dateiabhängigkeiten von <boost/*>

- ugcore/ugbase/bindings/vrl/user_data.cpp: <boost/function.hpp>
- ugcore/ugbase/bridge/disc_bridges/user_data_bridge.cpp: <boost/function.hpp>
- ugcore/ugbase/bridge/disc_bridges/user_data_bridge.cpp.orig: <boost/function.hpp>
- ugcore/ugbase/bridge/misc_bridges/test_bridge.cpp: <boost/bind.hpp>
- ugcore/ugbase/common/util/message_hub.h: <boost/function.hpp>, <boost/function_equal.hpp>, <boost/bind.hpp>
- ugcore/ugbase/common/util/factory.h: <boost/static_assert.hpp>, <boost/mpl/for_each.hpp>, <boost/mpl/vector.hpp>, <boost/type_traits/is_base_of.hpp>
- ugcore/ugbase/common/util/archivar.h: <boost/static_assert.hpp>, <boost/mpl/for_each.hpp>, <boost/type_traits/is_base_of.hpp>
- ugcore/ugbase/common/util/field.h: <boost/serialization/split_member.hpp>
- ugcore/ugbase/common/util/base64_file_writer.cpp: <boost/archive/iterators/transform_width.hpp>, <boost/archive/iterators/base64_from_binary.hpp>, <boost/archive/iterators/ostream_iterator.hpp>
- ugcore/ugbase/common/util/smart_pointer.h: <boost/pointee.hpp>
- ugcore/ugbase/common/util/bucket_sorter.hpp: (auskommentiert) <boost/limits.hpp>
- ugcore/ugbase/common/util/end_boost_list.h: <boost/mpl/list.hpp>
- ugcore/ugbase/common/boost_serialization.h: <boost/serialization/access.hpp>, <boost/serialization/export.hpp>, <boost/serialization/level.hpp>, <boost/serialization/nvp.hpp>, <boost/serialization/version.hpp>
- ugcore/ugbase/lib_grid/grid/grid.h: <boost/function.hpp>
- ugcore/ugbase/lib_grid/grid_objects/grid_dim_traits.h: <boost/mpl/list.hpp>
- ugcore/ugbase/lib_grid/parallelization/parallelization_util.h: <boost/function.hpp>
- ugcore/ugbase/lib_grid/refinement/projectors/projectors.h: <boost/mpl/pair.hpp>, <boost/mpl/string.hpp>, <boost/mpl/vector.hpp>
- ugcore/ugbase/lib_grid/refinement/projectors/neurite_projector.h: <boost/serialization/split_member.hpp>
- ugcore/ugbase/lib_grid/refinement/projectors/neurite_projector.cpp: <boost/lexical_cast.hpp>
- ugcore/ugbase/lib_grid/file_io/file_io_ugx.cpp: <boost/archive/text_oarchive.hpp>, <boost/archive/text_iarchive.hpp>
- ugcore/ugbase/lib_grid/file_io/file_io_lgb.cpp: <boost/archive/text_oarchive.hpp>, <boost/archive/text_iarchive.hpp>
- ugcore/ugbase/lib_grid/file_io/file_io_swc.cpp: <boost/lexical_cast.hpp>
- ugcore/ugbase/lib_grid/algorithms/extrusion/expand_layers.cpp: <boost/function.hpp>
- ugcore/ugbase/lib_grid/algorithms/extrusion/expand_layers_arte.cpp: <boost/function.hpp>
- ugcore/ugbase/lib_grid/algorithms/extrusion/expand_layers_arte3D.cpp: <boost/function.hpp>
- ugcore/ugbase/lib_grid/algorithms/extrusion/ArteExpandFracs3D.h/.cpp: <boost/function.hpp>
- ugcore/ugbase/lib_grid/algorithms/serialization.cpp: (auskommentiert) <boost/archive/text_oarchive.hpp>, <boost/archive/text_iarchive.hpp>
- ugcore/ugbase/lib_grid/tools/periodic_boundary_manager_impl.hpp: <boost/mpl/map.hpp>, <boost/mpl/at.hpp>
- ugcore/ugbase/lib_algebra/ordering_strategies/algorithms/topological_ordering.h: <boost/graph/adjacency_list.hpp>, <boost/graph/graph_traits.hpp>, <boost/graph/properties.hpp>
- ugcore/ugbase/lib_algebra/ordering_strategies/algorithms/SCC_ordering.h: <boost/graph/adjacency_list.hpp>, <boost/graph/graph_traits.hpp>, <boost/graph/properties.hpp>, <boost/graph/graph_utility.hpp>, <boost/graph/strong_components.hpp>
- ugcore/ugbase/lib_algebra/ordering_strategies/algorithms/boost_minimum_degree_ordering.h: <boost/graph/adjacency_list.hpp>, <boost/graph/graph_traits.hpp>, <boost/graph/properties.hpp>, <boost/graph/minimum_degree_ordering.hpp>
- ugcore/ugbase/lib_algebra/ordering_strategies/algorithms/boost_cuthill_mckee_ordering.h: <boost/graph/adjacency_list.hpp>, <boost/graph/graph_traits.hpp>, <boost/graph/properties.hpp>, <boost/graph/cuthill_mckee_ordering.hpp>
- ugcore/ugbase/lib_algebra/ordering_strategies/algorithms/util.h: <boost/graph/adjacency_list.hpp>, <boost/graph/graph_traits.hpp>
- ugcore/ugbase/lib_algebra/graph_interface/undirected_boost.h: <boost/iterator/counting_iterator.hpp>
- ugcore/ugbase/lib_algebra/graph_interface/boost_util.h: <boost/iterator/filter_iterator.hpp>
- ugcore/ugbase/lib_algebra/graph_interface/parallel_matrix.h: <boost/iterator/filter_iterator.hpp>
- ugcore/ugbase/lib_algebra/graph_interface/sparsematrix_boost.h: <boost/graph/properties.hpp>, <boost/iterator/counting_iterator.hpp>
- ugcore/ugbase/lib_algebra/graph_interface/undirected.h: <boost/geometry/iterators/concatenate_iterator.hpp>, <boost/graph/adjacency_list.hpp>
- ugcore/ugbase/lib_algebra/operator/interface/constrained_linear_iterator.h: <boost/core/enable_if.hpp>, <boost/type_traits/is_base_of.hpp>
- ugcore/ugbase/lib_algebra/common/matrixio/matrix_io_mtx.h: <boost/algorithm/string.hpp>, <boost/lexical_cast.hpp>
- ugcore/ugbase/registry/class.h: <boost/type_traits.hpp>, <boost/optional.hpp>
- ugcore/ugbase/registry/registry.h: <boost/function.hpp>, <boost/type_traits.hpp>
- ugcore/ugbase/lib_disc/function_spaces/grid_function_util.h: <boost/function.hpp>
- ugcore/ugbase/lib_disc/function_spaces/integrate.h: <boost/function.hpp>
- ugcore/ugbase/lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h: <boost/function.hpp>
- ugcore/ugbase/lib_disc/spatial_disc/elem_disc/dirac_source/lagrange_dirac_source.h: <boost/function.hpp>
- ugcore/ugbase/lib_disc/spatial_disc/elem_disc/err_est_data.h: <boost/mpl/for_each.hpp>
- ugcore/ugbase/lib_disc/spatial_disc/disc_util/conv_shape.h: <boost/mpl/range_c.hpp>, <boost/mpl/for_each.hpp>
- ugcore/ugbase/lib_disc/domain_impl.h: <boost/archive/text_oarchive.hpp>, <boost/archive/text_iarchive.hpp>
- ugcore/ugbase/lib_disc/reference_element/element_list_traits.h: <boost/mpl/transform_view.hpp>, <boost/mpl/fold.hpp>, <boost/mpl/min_max.hpp>
- ugcore/ugbase/bridge/util_algebra_dependent.h: <boost/mpl/if.hpp>, <boost/mpl/list.hpp>, <boost/mpl/empty.hpp>, <boost/mpl/front.hpp>, <boost/mpl/pop_front.hpp>
- ugcore/ugbase/bridge/util_domain_dependent.h: <boost/mpl/if.hpp>, <boost/mpl/list.hpp>, <boost/mpl/empty.hpp>, <boost/mpl/front.hpp>, <boost/mpl/pop_front.hpp>
- ugcore/ugbase/lib_disc/ordering_strategies/algorithms/riverorder.h: <boost/graph/adjacency_list.hpp>, <boost/graph/graph_traits.hpp>, <boost/graph/properties.hpp>
- ugcore/ugbase/lib_disc/ordering_strategies/algorithms/directional_ordering.cpp: <boost/graph/adjacency_list.hpp>, <boost/graph/graph_traits.hpp>, <boost/graph/properties.hpp>, <boost/graph/cuthill_mckee_ordering.hpp>

Hinweis: Die Liste basiert auf einer Quelltextsuche nach "#include <boost/...>" unter `ugcore/ugbase`. Manche Einträge sind auskommentiert oder mehrfach vorhanden (z.B. Header + generierte Build-Dependency-Dateien). Wenn du möchtest, kann ich die Liste deduplizieren, vollständig als CSV/JSON exportieren oder fehlende Treffer nachprüfen.

Dateien in `ugcore/ugbase`, die <boost/function.hpp> einbinden

- ugcore/ugbase/bindings/vrl/user_data.cpp
- ugcore/ugbase/lib_grid/grid/grid.h
- ugcore/ugbase/lib_grid/parallelization/parallelization_util.h
- ugcore/ugbase/lib_grid/algorithms/extrusion/expand_layers_arte3D.cpp
- ugcore/ugbase/lib_grid/algorithms/extrusion/expand_layers.cpp
- ugcore/ugbase/lib_grid/algorithms/extrusion/expand_layers_arte.cpp
- ugcore/ugbase/lib_grid/algorithms/extrusion/ArteExpandFracs3D.h
- ugcore/ugbase/lib_grid/algorithms/extrusion/ArteExpandFracs3D.cpp
- ugcore/ugbase/registry/registry.h
- ugcore/ugbase/lib_disc/function_spaces/grid_function_util.h
- ugcore/ugbase/lib_disc/function_spaces/integrate.h
- ugcore/ugbase/lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h
- ugcore/ugbase/lib_disc/spatial_disc/elem_disc/dirac_source/lagrange_dirac_source.h
- ugcore/ugbase/common/util/message_hub.h
- ugcore/ugbase/pcl/pcl_layout_tests.h
