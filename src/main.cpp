#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<HalfedgeMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

// Some algorithm parameters
float tCoef = 1.0;
std::unique_ptr<HeatMethodDistanceSolver> solver;
int vertexInd = 0;

// Example computation function -- this one computes and registers a scalar
// quantity
void doWork() {

  if (solver == nullptr) {
    solver.reset(new HeatMethodDistanceSolver(*geometry, tCoef));
  }

  /*
  std::vector<int> inds{1631};
  std::vector<SurfacePoint> points;
  for (int ind : inds) {
    Face f = mesh->face(ind);
    points.emplace_back(f, Vector3::constant(1. / 3.));
  }
  VertexData<double> distance = solver->computeDistance(points);
  */

  VertexData<double> distance =
      solver->computeDistance(mesh->vertex(vertexInd));

  psMesh->addVertexDistanceQuantity("distance", distance);
}

void myCallback() {

  if (ImGui::Button("make it so")) {
    doWork();
  }

  if (ImGui::InputFloat("tCoef", &tCoef)) {
    solver.reset();
  }

  ImGui::InputInt("vertex index", &vertexInd);
}

int main(int argc, char **argv) {

  // Configure the argument parser
  args::ArgumentParser parser("A demo of the Vector Heat Method");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  if (!inputFilename) {
    std::cerr << "Please specify a mesh file as argument" << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize polyscope
  polyscope::init();

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // Load mesh
  std::tie(mesh, geometry) = loadMesh(args::get(inputFilename));

  // Register the mesh with polyscope
  psMesh = polyscope::registerSurfaceMesh(
      polyscope::guessNiceNameFromPath(args::get(inputFilename)),
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
