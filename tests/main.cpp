#include <cstdlib>
#include <iostream>

void test_surface_mesh_basic();
void test_surface_mesh_split_edge_and_split_face();
void test_manifold_mesh_single_triangle_boundary_loop();
void test_manifold_mesh_tetrahedron_closed();
void test_manifold_vertex_halfedges();
void test_property_edge_length_updates();
void test_split_face_into_triangles();

int main() {
  test_surface_mesh_basic();
  test_surface_mesh_split_edge_and_split_face();
  test_manifold_mesh_single_triangle_boundary_loop();
  test_manifold_mesh_tetrahedron_closed();
  test_manifold_vertex_halfedges();
  test_property_edge_length_updates();
  test_split_face_into_triangles();

  std::cout << "mesh1_tests: OK\n";
  return EXIT_SUCCESS;
}
