# Übersicht
In folgenden Dateien in `ugcore/ugbase`, wurden die Abhängigkeiten von <boost/function.hpp> entfernt:

## Dateiliste
o ugcore/ugbase/bindings/vrl/user_data.cpp
f ugcore/ugbase/lib_grid/grid/grid.h
f ugcore/ugbase/lib_grid/parallelization/parallelization_util.h
o ugcore/ugbase/lib_grid/algorithms/extrusion/expand_layers_arte3D.cpp
o ugcore/ugbase/lib_grid/algorithms/extrusion/expand_layers.cpp
o ugcore/ugbase/lib_grid/algorithms/extrusion/expand_layers_arte.cpp
o ugcore/ugbase/lib_grid/algorithms/extrusion/ArteExpandFracs3D.h
 ougcore/ugbase/lib_grid/algorithms/extrusion/ArteExpandFracs3D.cpp
f ugcore/ugbase/registry/registry.h
o ugcore/ugbase/lib_disc/function_spaces/grid_function_util.h
o ugcore/ugbase/lib_disc/function_spaces/integrate.h
m ugcore/ugbase/lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h
o ugcore/ugbase/lib_disc/spatial_disc/elem_disc/dirac_source/lagrange_dirac_source.h
f ugcore/ugbase/common/util/message_hub.h
f ugcore/ugbase/pcl/pcl_layout_tests.h

## Legende
f: Done (fixed)
o: Done (obsolet)
m: Neue Header (MPL).




Weitere Dateiabhängigkeiten von <boost/*>

# Resolved
## Type traits
- ugcore/ugbase/lib_algebra/operator/interface/constrained_linear_iterator.h: 
f: <boost/core/enable_if.hpp>, 
f: <boost/type_traits/is_base_of.hpp>
- ugcore/ugbase/registry/class.h: 
f: <boost/type_traits.hpp>, 
<boost/optional.hpp> => CXX17
- ugcore/ugbase/registry/registry.h:  
f: <boost/type_traits.hpp>
- ugcore/ugbase/common/util/bucket_sorter.hpp: (auskommentiert) 
f: <boost/limits.hpp>


## Bind
- ugcore/ugbase/bridge/misc_bridges/test_bridge.cpp: 
o: <boost/bind.hpp>
- ugcore/ugbase/common/util/message_hub.h: 
f: <boost/bind.hpp>

## Pointee
f: ugcore/ugbase/common/util/smart_pointer.h: <boost/pointee.hpp>

# Unresolved
## Lexical cast
- ugcore/ugbase/lib_grid/refinement/projectors/neurite_projector.cpp: 
<boost/lexical_cast.hpp>
- ugcore/ugbase/lib_grid/file_io/file_io_swc.cpp: 
<boost/lexical_cast.hpp>
- ugcore/ugbase/lib_algebra/common/matrixio/matrix_io_mtx.h: 
<boost/lexical_cast.hpp>,
<boost/algorithm/string.hpp>


## Iterator
- ugcore/ugbase/lib_algebra/graph_interface/undirected_boost.h: 
<boost/iterator/counting_iterator.hpp>
- ugcore/ugbase/lib_algebra/graph_interface/boost_util.h:
 <boost/iterator/filter_iterator.hpp>
- ugcore/ugbase/lib_algebra/graph_interface/parallel_matrix.h: 
<boost/iterator/filter_iterator.hpp>
- ugcore/ugbase/lib_algebra/graph_interface/sparsematrix_boost.h: 
<boost/iterator/counting_iterator.hpp>
<boost/graph/properties.hpp> 
- ugcore/ugbase/lib_algebra/graph_interface/undirected.h: 
<boost/geometry/iterators/concatenate_iterator.hpp>,
<boost/graph/adjacency_list.hpp>


## Archive
- ugcore/ugbase/common/util/base64_file_writer.cpp: 
 <boost/archive/iterators/transform_width.hpp>,
 <boost/archive/iterators/base64_from_binary.hpp>,
 <boost/archive/iterators/ostream_iterator.hpp>
- ugcore/ugbase/lib_disc/domain_impl.h: 
<boost/archive/text_oarchive.hpp>, 
<boost/archive/text_iarchive.hpp>
- ugcore/ugbase/lib_grid/file_io/file_io_ugx.cpp:
 <boost/archive/text_oarchive.hpp>,
  <boost/archive/text_iarchive.hpp>
- ugcore/ugbase/lib_grid/file_io/file_io_lgb.cpp:
 <boost/archive/text_oarchive.hpp>,
  <boost/archive/text_iarchive.hpp>
- ugcore/ugbase/lib_grid/algorithms/serialization.cpp: (auskommentiert) 
<boost/archive/text_oarchive.hpp>, 
<boost/archive/text_iarchive.hpp>

## MPL
- ugcore/ugbase/common/util/factory.h: 
 f: <boost/static_assert.hpp>
 f: <boost/type_traits/is_base_of.hpp>
    <boost/mpl/for_each.hpp>, 
    <boost/mpl/vector.hpp>,
- ugcore/ugbase/common/util/archivar.h: 
 f: <boost/static_assert.hpp>,
 f: <boost/type_traits/is_base_of.hpp>
 <boost/mpl/for_each.hpp>, 
- ugcore/ugbase/common/util/end_boost_list.h: 
<boost/mpl/list.hpp>
- ugcore/ugbase/lib_disc/spatial_disc/elem_disc/err_est_data.h: 
<boost/mpl/for_each.hpp>
- ugcore/ugbase/lib_disc/spatial_disc/disc_util/conv_shape.h: 
<boost/mpl/range_c.hpp>, 
<boost/mpl/for_each.hpp>
- ugcore/ugbase/lib_grid/tools/periodic_boundary_manager_impl.hpp: 
<boost/mpl/map.hpp>, 
<boost/mpl/at.hpp>
- ugcore/ugbase/lib_disc/reference_element/element_list_traits.h: 
<boost/mpl/transform_view.hpp>, 
<boost/mpl/fold.hpp>,
 <boost/mpl/min_max.hpp>
- ugcore/ugbase/bridge/util_algebra_dependent.h: 
<boost/mpl/if.hpp>, 
<boost/mpl/list.hpp>, 
<boost/mpl/empty.hpp>, 
<boost/mpl/front.hpp>, 
<boost/mpl/pop_front.hpp>
- ugcore/ugbase/bridge/util_domain_dependent.h: 
<boost/mpl/if.hpp>, 
<boost/mpl/list.hpp>, 
<boost/mpl/empty.hpp>, 
<boost/mpl/front.hpp>, 
<boost/mpl/pop_front.hpp>
- ugcore/ugbase/lib_grid/grid_objects/grid_dim_traits.h:
 <boost/mpl/list.hpp>
- ugcore/ugbase/lib_grid/refinement/projectors/projectors.h: 
<boost/mpl/pair.hpp>, 
<boost/mpl/string.hpp>, 
<boost/mpl/vector.hpp>

## Graph
- ugcore/ugbase/lib_disc/ordering_strategies/algorithms/riverorder.h: 
<boost/graph/adjacency_list.hpp>, 
<boost/graph/graph_traits.hpp>, 
<boost/graph/properties.hpp>
- ugcore/ugbase/lib_disc/ordering_strategies/algorithms/directional_ordering.cpp:
 <boost/graph/adjacency_list.hpp>, 
 <boost/graph/graph_traits.hpp>, 
 <boost/graph/properties.hpp>,
 <boost/graph/cuthill_mckee_ordering.hpp>
 - ugcore/ugbase/lib_algebra/ordering_strategies/algorithms/topological_ordering.h: 
<boost/graph/adjacency_list.hpp>, 
<boost/graph/graph_traits.hpp>, 
<boost/graph/properties.hpp>
- ugcore/ugbase/lib_algebra/ordering_strategies/algorithms/SCC_ordering.h: 
<boost/graph/adjacency_list.hpp>, 
<boost/graph/graph_traits.hpp>, 
<boost/graph/properties.hpp>, 
<boost/graph/graph_utility.hpp>, 
<boost/graph/strong_components.hpp>
- ugcore/ugbase/lib_algebra/ordering_strategies/algorithms/boost_minimum_degree_ordering.h: 
<boost/graph/adjacency_list.hpp>, 
<boost/graph/graph_traits.hpp>, 
<boost/graph/properties.hpp>, 
<boost/graph/minimum_degree_ordering.hpp>
- ugcore/ugbase/lib_algebra/ordering_strategies/algorithms/boost_cuthill_mckee_ordering.h: 
<boost/graph/adjacency_list.hpp>,
 <boost/graph/graph_traits.hpp>, 
 <boost/graph/properties.hpp>, 
 <boost/graph/cuthill_mckee_ordering.hpp>
- ugcore/ugbase/lib_algebra/ordering_strategies/algorithms/util.h: 
<boost/graph/adjacency_list.hpp>, 
<boost/graph/graph_traits.hpp>

## Serialization
- ugcore/ugbase/common/util/field.h: 
<boost/serialization/split_member.hpp>
- ugcore/ugbase/common/boost_serialization.h: 
<boost/serialization/access.hpp>,
 <boost/serialization/export.hpp>,
<boost/serialization/level.hpp>, 
<boost/serialization/nvp.hpp>, 
<boost/serialization/version.hpp>
- ugcore/ugbase/lib_grid/refinement/projectors/neurite_projector.h: 
<boost/serialization/split_member.hpp>


Hinweis: Die Liste basiert auf einer Quelltextsuche nach "#include <boost/...>" unter `ugcore/ugbase`. Manche Einträge sind auskommentiert oder mehrfach vorhanden (z.B. Header + generierte Build-Dependency-Dateien). 


## Migrationsvorschlag — Reihenfolge (priorisiert, kurz)

1. Funktionale-APIs (niedriges Risiko): Ersetze boost::function / boost::bind durch std::function, std::bind oder besser: C++11-Lambdas. Dateien zuerst prüfen: alle Einträge in boost_undo.md, z.B. user_data.cpp, ugcore/ugbase/lib_grid/.... Aufwand: niedrig. Test: kompilieren + Unittests.

2. Smart-Pointer (niedriges Risiko): Ersetze boost::shared_ptr, boost::weak_ptr, boost::make_shared durch std::shared_ptr, std::weak_ptr, std::make_shared. Aufwand: niedrig–mittel (API-Suche/typedefs). Test: kompilieren.

3. Type Traits & Metaprogramming (niedriges bis mittleres Risiko): Migriere boost::is_base_of, boost::enable_if, boost::type_traits → <type_traits> (std::is_base_of, std::enable_if_t/std::enable_if). Für MPL- (Boost.MPL) heavy uses: erst analysieren; viele MPL-Patterns bleiben (siehe bridge/*, lib_disc/*). Aufwand: mittel. Test: kompilieren.

4. String/Conversion (niedriges Risiko): Ersetze boost::lexical_cast durch std::to_string, std::stoi/std::stod oder std::stringstream wo passend. Aufwand: niedrig.

5. Threading (mittleres Risiko): Ersetze boost::thread, boost::mutex, boost::condition_variable durch std::thread, std::mutex, std::condition_variable. Achte auf thread-API-Differenzen (interrupts, join/detach). Aufwand: mittel. Test: Laufzeit- und race-tests.

6. Iterators / Algorithmic helpers (mittel): Manche Boost-Iteratoren/Adaptoren (z.B. filter_iterator, counting_iterator) fehlen in STL — entweder implementieren kleine wrappers oder verwenden einfache loops/lambdas. Aufwand: mittel.

7. Small utilities (niedrig): Ersetze boost::pointee, boost::function_equal etc. durch kleine helper-Funktionen / Standard-Pattern. Aufwand: niedrig.

8. I/O / Archiving / Serialization / Graph / MPL (hohes Risiko, defer):

Boost.Serialization / Boost.Archive (boost/archive/...) hat keine direkte C++11-Ersatzlösung — entweder behalten oder auf Alternativen ( cereal, protobuf ) migrieren (großer Aufwand).
Boost.Graph, Boost.MPL, Boost.Geometry sind größere Subsysteme; nur migrieren wenn zwingend (hoher Aufwand / Design-Entscheidung).
Empfehlung: diese Schritte ans Ende setzen und erst planen, wenn kleinere Ersetzungen abgeschlossen und Tests stabil sind.
9. Build & CI (parallel, kontinuierlich): Nach jeder Gruppen-Migration: Branch erstellen, CMake anpassen (schrittweise Boost-Componenten entfernen), komplette Build + Tests laufen lassen. Dokumentation der Änderungen in boost_undo.md/CHANGES.

10. Final: Aufräumen & Deduplizieren: Entferne jetzt unnötige #include <boost/...>-Zeilen, dedupliziere boost_undo.md-Fundstellen, committen und ggf. PR.

## Kurz-Checklist pro Migrationsschritt:

a. Grep alle Vorkommen (Scope: ugbase / Dateien in boost_undo.md).
b. Ändern mit kleinen, isolierten Kommits.
c. Kompilieren + Tests ausführen.
d. CI grün → weitermachen.
