#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include <sstream>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<HalfedgeMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

// Some algorithm parameters
float tCoef = 1.0;
std::unique_ptr<VectorHeatMethodSolver> solver;
int vertexInd = 0;

// Manage a list of sources
struct SourceVert {
  Vertex vertex;
  float scalarVal = 1.;
  float vectorMag = 1.;
  float vectorAngleRad = 0.;
};
std::vector<SourceVert> sourcePoints;

bool vizFirstRun = true;
void updateSourceSetViz() {

  // Scalar balls around sources
  std::vector<std::pair<size_t, double>> sourcePairs;
  for (SourceVert& s : sourcePoints) {
    size_t ind = geometry->vertexIndices[s.vertex];
    sourcePairs.emplace_back(ind, s.scalarVal);
  }
  auto scalarQ = polyscope::getSurfaceMesh()->addVertexIsolatedScalarQuantity("source scalars", sourcePairs);
  scalarQ->pointRadius *= 2.;
  scalarQ->cMap = polyscope::gl::ColorMapID::REDS;
  if (vizFirstRun) {
    scalarQ->setEnabled(true);
  }

  // Vectors at sources
  VertexData<Vector2> sourceVectors(*mesh, Vector2::zero());
  for (SourceVert& s : sourcePoints) {
    sourceVectors[s.vertex] = Vector2::fromAngle(s.vectorAngleRad) * s.vectorMag;
  }
  auto vectorQ = polyscope::getSurfaceMesh()->addVertexIntrinsicVectorQuantity("source vectors", sourceVectors);
  vectorQ->lengthMult *= 2.;
  vectorQ->radiusMult *= 4.;
  vectorQ->vectorColor = glm::vec3{227 / 255., 52 / 255., 28 / 255.};
  if (vizFirstRun) {
    vectorQ->setEnabled(true);
  }

  vizFirstRun = false;
}

void addVertexSource(size_t ind) {
  Vertex v = mesh->vertex(ind);

  // Make sure not already used
  for (SourceVert& s : sourcePoints) {
    if (s.vertex == v) {
      std::stringstream ss;
      ss << "Vertex " << v;
      std::string vStr = ss.str();
      polyscope::warning("Vertex " + vStr + " is already a source");
      return;
    }
  }

  SourceVert newV;
  newV.vertex = v;
  sourcePoints.push_back(newV);
  updateSourceSetViz();
}

void scalarExtension() {

  if (solver == nullptr) {
    solver.reset(new VectorHeatMethodSolver(*geometry, tCoef));
  }

  if (sourcePoints.size() == 0) {
    polyscope::warning("no source points set");
    return;
  }

  std::vector<std::tuple<SurfacePoint, double>> points;
  for (SourceVert& s : sourcePoints) {
    points.emplace_back(s.vertex, s.scalarVal);
  }

  VertexData<double> scalarExtension = solver->extendScalar(points);

  psMesh->addVertexScalarQuantity("scalar extension", scalarExtension);
}

void vectorTransport() {

  if (solver == nullptr) {
    solver.reset(new VectorHeatMethodSolver(*geometry, tCoef));
  }

  if (sourcePoints.size() == 0) {
    polyscope::warning("no source points set");
    return;
  }

  std::vector<std::tuple<SurfacePoint, Vector2>> points;
  for (SourceVert& s : sourcePoints) {
    points.emplace_back(s.vertex, Vector2::fromAngle(s.vectorAngleRad) * s.vectorMag);
  }
  VertexData<Vector2> vectorExtension = solver->transportTangentVectors(points);

  psMesh->addVertexIntrinsicVectorQuantity("vector extension", vectorExtension);
}

void buildPointsMenu() {

  bool anyChanged = false;

  ImGui::PushItemWidth(200);

  int id = 0;
  int eraseInd = -1;
  for (SourceVert& s : sourcePoints) {
    std::stringstream ss;
    ss << "Vertex " << s.vertex;
    std::string vStr = ss.str();
    ImGui::PushID(vStr.c_str());

    ImGui::TextUnformatted(vStr.c_str());

    ImGui::SameLine();
    if (ImGui::Button("delete")) {
      eraseInd = id;
      anyChanged = true;
    }
    ImGui::Indent();

    if (ImGui::InputFloat("scalar value", &s.scalarVal)) anyChanged = true;
    if (ImGui::InputFloat("vector mag", &s.vectorMag)) anyChanged = true;
    if (ImGui::SliderAngle("vector angle", &s.vectorAngleRad)) anyChanged = true;

    ImGui::Unindent();
    ImGui::PopID();
  }
  ImGui::PopItemWidth();

  // actually do erase, if requested
  if (eraseInd != -1) {
    sourcePoints.erase(sourcePoints.begin() + eraseInd);
  }

  if (ImGui::Button("add point")) {
    long long int pickVert = polyscope::getSurfaceMesh()->selectVertex();
    if (pickVert >= 0) {
      addVertexSource(pickVert);
      anyChanged = true;
    }
  }

  if (anyChanged) {
    updateSourceSetViz();
  }
}

void myCallback() {

  ImGui::TextUnformatted("Algorithm options:");
  ImGui::PushItemWidth(300);
  if (ImGui::InputFloat("tCoef", &tCoef)) {
    solver.reset();
  }
  ImGui::PopItemWidth();

  // Build the list of source points
  if (ImGui::TreeNode("source points")) {
    buildPointsMenu();
    ImGui::TreePop();
  }

  if (ImGui::Button("run scalar extension")) {
    scalarExtension();
  }

  if (ImGui::Button("run vector transport")) {
    vectorTransport();
  }
}

int main(int argc, char** argv) {

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
  psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(args::get(inputFilename)),
                                          geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                          polyscopePermutations(*mesh));


  // Set vertex tangent spaces
  geometry->requireVertexTangentBasis();
  VertexData<Vector3> vBasisX(*mesh);
  for(Vertex v : mesh->vertices()) {
    vBasisX[v] = geometry->vertexTangentBasis[v][0];
  }
  polyscope::getSurfaceMesh()->setVertexTangentBasisX(vBasisX);
  
  // Set face tangent spaces
  geometry->requireFaceTangentBasis();
  FaceData<Vector3> fBasisX(*mesh);
  for(Face f : mesh->faces()) {
    fBasisX[f] = geometry->faceTangentBasis[f][0];
  }
  polyscope::getSurfaceMesh()->setFaceTangentBasisX(fBasisX);

  // To start, pick two vertices as sources
  geometry->requireVertexIndices();
  addVertexSource(0);
  addVertexSource(mesh->nVertices() / 2);

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
