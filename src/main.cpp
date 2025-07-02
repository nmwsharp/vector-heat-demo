#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/surface_centers.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/utilities/vector3.h"
#include "polyscope/options.h"
#include "polyscope/pick.h"
#include "polyscope/point_cloud.h"
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
int pCenter = 2;
bool useIntrinsicTriangulation = true;
std::unique_ptr<SignpostIntrinsicTriangulation> signpostTri;

// other UI other
std::string UImode = ""; // "extension", "logmap", or "centers"

// Manage a list of sources for extension
struct SourceVert {
  Vertex vertex;
  float scalarVal = 1.;
  float vectorMag = 1.;
  float vectorAngleRad = 0.;
};
std::vector<SourceVert> extensionSourcePoints;

// Manage a single source vertex for logmap
Vertex logmapSourceVertex;
SurfacePoint logmapSourcePoint;
bool logMapContinuous = false;
Vector3 logmapContinuousLastXDir = Vector3{1., 0., 0.};
Vector3 logmapContinuousLastNormal = Vector3{0., 0., 1.};
VertexData<Vector2> lastLogmap;
LogMapStrategy logMapStrategy = LogMapStrategy::AffineLocal;
float logMapRadius = -1.;
float logMapRadiusUpper = -1.;
bool limitLogMapRadius = false;

// Manage a list of sites for centers/averages
struct SiteVert {
  Vertex vertex;
  float weight = 1.;
};
std::vector<SiteVert> centerSiteVerts;

std::tuple<VertexData<Vector3>, VertexData<Vector3>> getTangentVectors() {
  VertexData<Vector3> basisX(*mesh);
  VertexData<Vector3> basisY(*mesh);
  for (Vertex v : mesh->vertices()) {
    basisX[v] = geometry->vertexTangentBasis[v][0];
    basisY[v] = geometry->vertexTangentBasis[v][1];
  }
  return std::make_tuple(basisX, basisY);
}

bool vizFirstRun = true;
void updateSourceSetViz() {

  // Scalar balls around sources
  std::vector<Vector3> sourcePositions;
  std::vector<double> sourceValues;
  for (SourceVert& s : extensionSourcePoints) {
    size_t ind = geometry->vertexIndices[s.vertex];
    sourcePositions.push_back(geometry->inputVertexPositions[s.vertex]);
    sourceValues.push_back(s.scalarVal);
  }
  auto pointQ = polyscope::registerPointCloud("source points", sourcePositions);
  auto scalarQ = pointQ->addScalarQuantity("source scalars", sourceValues);
  pointQ->setPointRadius(0.015);
  scalarQ->setColorMap("reds");
  if (vizFirstRun) {
    scalarQ->setEnabled(true);
  }

  // Vectors at sources
  VertexData<Vector3> basisX, basisY;
  std::tie(basisX, basisY) = getTangentVectors();
  std::vector<Vector3> sourceVectors;
  for (SourceVert& s : extensionSourcePoints) {
    Vector2 vec = Vector2::fromAngle(s.vectorAngleRad) * s.vectorMag;
    Vector3 vec3D = basisX[s.vertex] * vec.x + basisY[s.vertex] * vec.y;
    sourceVectors.push_back(vec3D);
  }
  auto vectorQ = pointQ->addVectorQuantity("source vectors", sourceVectors);
  vectorQ->setVectorLengthScale(.05);
  vectorQ->setVectorRadius(.005);
  vectorQ->setVectorColor(glm::vec3{227 / 255., 52 / 255., 28 / 255.});
  if (vizFirstRun) {
    vectorQ->setEnabled(true);
  }

  vizFirstRun = false;
}

bool vizFirstRunSite = true;
void updateSiteSetViz() {

  // Scalar balls around sources
  std::vector<Vector3> sitePositions;
  std::vector<double> siteValues;
  for (SiteVert& s : centerSiteVerts) {
    size_t ind = geometry->vertexIndices[s.vertex];
    sitePositions.push_back(geometry->inputVertexPositions[s.vertex]);
    siteValues.push_back(s.weight);
  }
  auto pointQ = polyscope::registerPointCloud("site points", sitePositions);
  auto scalarQ = pointQ->addScalarQuantity("source scalars", siteValues);
  pointQ->setPointRadius(0.015);
  scalarQ->setColorMap("blues");
  if (vizFirstRunSite) {
    scalarQ->setEnabled(true);
  }

  vizFirstRunSite = false;
}
void addVertexSource(size_t ind) {
  Vertex v = mesh->vertex(ind);

  // Make sure not already used
  for (SourceVert& s : extensionSourcePoints) {
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
  extensionSourcePoints.push_back(newV);
  updateSourceSetViz();
}

void updateLogmapSourceViz() {
  auto pointQ = polyscope::registerPointCloud("logmap source",
                                              std::vector<Vector3>{geometry->inputVertexPositions[logmapSourceVertex]});
  pointQ->setPointRadius(0.015);
}


void addVertexSite(size_t ind) {
  Vertex v = mesh->vertex(ind);

  // Make sure not already used
  for (SiteVert& s : centerSiteVerts) {
    if (s.vertex == v) {
      std::stringstream ss;
      ss << "Vertex " << v;
      std::string vStr = ss.str();
      polyscope::warning("Vertex " + vStr + " is already a site");
      return;
    }
  }

  SiteVert newV;
  newV.vertex = v;
  newV.weight = 1.0;
  centerSiteVerts.push_back(newV);
  updateSiteSetViz();
}

void clearCachedSolvers() { solver.reset(); }

void ensureHaveIntrinsicTriangulation() {
  if (signpostTri != nullptr) {
    return;
  }

  signpostTri.reset(new SignpostIntrinsicTriangulation(*mesh, *geometry));
  signpostTri->flipToDelaunay();
  // signpostTri->delaunayRefine(20);
}

void ensureHaveSolver() {
  if (solver != nullptr) {
    return;
  }

  if (useIntrinsicTriangulation) {
    ensureHaveIntrinsicTriangulation();
    solver.reset(new VectorHeatMethodSolver(*signpostTri, tCoef));
  } else {
    solver.reset(new VectorHeatMethodSolver(*geometry, tCoef));
  }
}

void vectorExtension() {
  ensureHaveSolver();

  if (extensionSourcePoints.size() == 0) {
    polyscope::warning("no source points set");
    return;
  }

  // Prep the data, remap to intrinsic triangulation if using
  std::vector<std::tuple<SurfacePoint, double>> points;
  for (SourceVert& s : extensionSourcePoints) {
    if (useIntrinsicTriangulation) {
      points.emplace_back(signpostTri->equivalentPointOnIntrinsic(SurfacePoint(s.vertex)).vertex, s.scalarVal);
    } else {
      points.emplace_back(s.vertex, s.scalarVal);
    }
  }

  // Run the algorithm
  VertexData<double> scalarExtension = solver->extendScalar(points);

  // Copy back to original mesh, if using intrinsic triangulation
  if (useIntrinsicTriangulation) {
    VertexData<double> scalarExtensionOnInput = signpostTri->restrictToInput(scalarExtension);
    scalarExtension = scalarExtensionOnInput;
  }

  auto psScalar = psMesh->addVertexScalarQuantity("scalar extension", scalarExtension);
  psScalar->setEnabled(true);
}

void vectorTransport() {
  ensureHaveSolver();

  if (extensionSourcePoints.size() == 0) {
    polyscope::warning("no source points set");
    return;
  }

  // Prep the data, remap to intrinsic triangulation if using
  std::vector<std::tuple<SurfacePoint, Vector2>> points;
  for (SourceVert& s : extensionSourcePoints) {
    Vector2 vec = Vector2::fromAngle(s.vectorAngleRad) * s.vectorMag;
    if (useIntrinsicTriangulation) {
      // tangent spaces are aligned by construction
      points.emplace_back(signpostTri->equivalentPointOnIntrinsic(SurfacePoint(s.vertex)).vertex, vec);
    } else {
      points.emplace_back(s.vertex, vec);
    }
  }

  // Run the algorithm
  VertexData<Vector2> vectorExtension = solver->transportTangentVectors(points);

  // Copy back to original mesh, if using intrinsic triangulation
  if (useIntrinsicTriangulation) {
    VertexData<Vector2> vectorExtensionOnInput = signpostTri->restrictToInput(vectorExtension);
    vectorExtension = vectorExtensionOnInput;
  }

  VertexData<Vector3> basisX, basisY;
  std::tie(basisX, basisY) = getTangentVectors();
  auto vectorQ = psMesh->addVertexTangentVectorQuantity("extended vectors", vectorExtension, basisX, basisY);
  vectorQ->setEnabled(true);
}

void clipLogMap(VertexData<Vector2>& logmap) {
  for (Vertex v : mesh->vertices()) {
    Vector2 logmapVal = logmap[v];
    float dist = norm(logmapVal);
    if (dist > logMapRadius) {
      // logmap[v] = logmapVal * (logMapRadius / dist);
      logmap[v] = Vector2{std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity()};
      // logmap[v] = Vector2{0., 0.};
    }
  }
}


void computeLogMapFromVertex() {
  ensureHaveSolver();

  // Remap the vertex to the intrinsic triangulation, if using
  Vertex sourceVert = useIntrinsicTriangulation
                          ? signpostTri->equivalentPointOnIntrinsic(SurfacePoint(logmapSourceVertex)).vertex
                          : logmapSourceVertex;

  // Run the algorithm
  VertexData<Vector2> logmap = solver->computeLogMap(sourceVert, logMapStrategy);

  // Copy back to original mesh, if using intrinsic triangulation
  if (useIntrinsicTriangulation) {
    VertexData<Vector2> logmapOnInput = signpostTri->restrictToInput(logmap);
    logmap = logmapOnInput;
  }

  if (limitLogMapRadius) {
    clipLogMap(logmap);
  }

  auto psLogmap = psMesh->addLocalParameterizationQuantity("logmap", logmap);
  psLogmap->setEnabled(true);
}

void computeLogMapFromSurfacePoint() {
  if (logmapSourcePoint == SurfacePoint()) {
    return;
  }
  if (logmapSourcePoint.type != SurfacePointType::Face) {
    throw std::runtime_error("this UI assumes logmap source point is a face point, something must be wrong");
  }

  ensureHaveSolver();

  // Remap the source to the intrinsic triangulation, if using
  SurfacePoint sourcePoint =
      useIntrinsicTriangulation ? signpostTri->equivalentPointOnIntrinsic(logmapSourcePoint) : logmapSourcePoint;

  // Run the algorithm
  VertexData<Vector2> logmap = solver->computeLogMap(sourcePoint, logMapStrategy);

  // Copy back to original mesh, if using intrinsic triangulation
  if (useIntrinsicTriangulation) {
    VertexData<Vector2> logmapOnInput = signpostTri->restrictToInput(logmap);
    logmap = logmapOnInput;
  }

  psMesh->addLocalParameterizationQuantity("logmap pre-rot", logmap);

  if (lastLogmap.size() > 0) { // Rotate the logmap to match the last x-direction

    // We could also read this off directly via a change-of-basis, but this is tricky (though possible) for the
    // intrinsic triangulation

    Vector2 rot{0., 0.};
    for (Vertex v : mesh->vertices()) {
      Vector2 newVal = logmap[v];
      Vector2 oldVal = lastLogmap[v];
      double weight = 1.0; // NOTE:  not obvious what the best weight is
      // double weight = norm(newVal);
      // double weight = 1. / norm(newVal);
      rot += weight * unit(newVal / oldVal);
    }
    rot = unit(rot);

    Vector2 invRot = rot.inv();
    for (Vertex v : mesh->vertices()) {
      logmap[v] = invRot * logmap[v];
    }
  }

  lastLogmap = logmap;

  if (limitLogMapRadius) {
    clipLogMap(logmap);
  }

  auto psLogmap = psMesh->addLocalParameterizationQuantity("logmap", logmap);
  psLogmap->setEnabled(true);
}

void computeCenter() {
  if (!(pCenter == 1 || pCenter == 2)) {
    polyscope::warning("p must be 1 or 2");
    return;
  }
  if (centerSiteVerts.size() == 0) {
    polyscope::warning("must select at least one site");
    return;
  }

  ensureHaveSolver();

  SurfacePoint center;
  if (useIntrinsicTriangulation) {

    // Remap the sites to the intrinsic triangulation
    VertexData<double> dist = VertexData<double>(*signpostTri->intrinsicMesh, 0.);
    for (SiteVert& s : centerSiteVerts) {
      dist[signpostTri->equivalentPointOnIntrinsic(SurfacePoint(s.vertex)).vertex] += s.weight;
    }

    // Run the algorithm
    center = findCenter(*signpostTri->intrinsicMesh, *signpostTri, *solver, dist, pCenter, logMapStrategy);

    // Transform back to original mesh
    center = signpostTri->equivalentPointOnInput(center);

  } else {

    VertexData<double> dist = VertexData<double>(*mesh, 0.);
    for (SiteVert& s : centerSiteVerts) {
      dist[s.vertex] += s.weight;
    }

    // Run the algorithm
    center = findCenter(*mesh, *geometry, *solver, dist, pCenter, logMapStrategy);
  }


  // Visualize
  Vector3 centerPos = center.interpolate(geometry->inputVertexPositions);
  std::vector<Vector3> centerPosCloud{centerPos};
  auto pointQ = polyscope::registerPointCloud("center", centerPosCloud);
  pointQ->setPointRadius(0.02);
}

void buildPointsMenu() {

  bool anyChanged = false;

  ImGui::PushItemWidth(200);

  int id = 0;
  int eraseInd = -1;
  for (SourceVert& s : extensionSourcePoints) {
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
    extensionSourcePoints.erase(extensionSourcePoints.begin() + eraseInd);
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

void buildSitesMenu() {

  bool anyChanged = false;

  ImGui::PushItemWidth(200);

  int id = 0;
  int eraseInd = -1;
  for (SiteVert& s : centerSiteVerts) {
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

    if (ImGui::InputFloat("weight", &s.weight)) anyChanged = true;

    ImGui::Unindent();
    ImGui::PopID();
  }
  ImGui::PopItemWidth();

  // actually do erase, if requested
  if (eraseInd != -1) {
    centerSiteVerts.erase(centerSiteVerts.begin() + eraseInd);
  }

  if (ImGui::Button("add site")) {
    long long int pickVert = polyscope::getSurfaceMesh()->selectVertex();
    if (pickVert >= 0) {
      addVertexSite(pickVert);
      anyChanged = true;
    }
  }

  if (anyChanged) {
    updateSiteSetViz();
  }
}


void buildLogmapStrategyMenu() {
  // Combobox to select which strategy to use
  std::map<LogMapStrategy, std::string> logMapStrategyNames = {{LogMapStrategy::AffineLocal, "AffineLocal"},
                                                               {LogMapStrategy::AffineAdaptive, "AffineAdaptive"},
                                                               {LogMapStrategy::VectorHeat, "VectorHeat"}};
  if (ImGui::BeginCombo("logmap strategy", logMapStrategyNames[logMapStrategy].c_str())) {
    if (ImGui::Selectable("VectorHeat", logMapStrategy == LogMapStrategy::VectorHeat)) {
      logMapStrategy = LogMapStrategy::VectorHeat;
    }
    if (ImGui::Selectable("AffineLocal", logMapStrategy == LogMapStrategy::AffineLocal)) {
      logMapStrategy = LogMapStrategy::AffineLocal;
    }
    if (ImGui::Selectable("AffineAdaptive", logMapStrategy == LogMapStrategy::AffineAdaptive)) {
      logMapStrategy = LogMapStrategy::AffineAdaptive;
    }
    ImGui::EndCombo();
  }
}

void myCallback() {

  std::string newUImode = UImode;

  ImGuiTabBarFlags tab_bar_flags = ImGuiTabBarFlags_None;
  if (ImGui::BeginTabBar("MyTabBar", tab_bar_flags)) {
    if (ImGui::BeginTabItem("Extension")) {
      newUImode = "extension";
      psMesh->setSelectionMode(polyscope::MeshSelectionMode::Auto);

      ImGui::TextWrapped(
          "Use the Vector Heat Method to extend scalars/vectors from isolated points across the surface, according "
          "to "
          "parallel transport along shortest geodesics. Use the menu below to "
          "select one or more source points, and set scalar/vector values per-point. The points are visualized as "
          "balls at the surface, and the menu pane on the left can toggle viewing scalar/vector values.");

      ImGui::NewLine();

      ImGui::TextUnformatted("Algorithm options:");
      ImGui::PushItemWidth(100);
      if (ImGui::InputFloat("tCoef", &tCoef)) {
        solver.reset();
      }
      ImGui::PopItemWidth();

      if (ImGui::Checkbox("use intrinsic delaunay triangulation", &useIntrinsicTriangulation)) {
        clearCachedSolvers();
      }

      // Build the list of source points
      if (ImGui::TreeNode("select source points")) {
        buildPointsMenu();
        ImGui::TreePop();
      }

      if (ImGui::Button("run scalar extension")) {
        vectorExtension();
      }

      if (ImGui::Button("run vector transport")) {
        vectorTransport();
      }


      ImGui::EndTabItem();
    }
    if (ImGui::BeginTabItem("Log map")) {
      newUImode = "logmap";

      ImGui::TextWrapped("Compute the logarithmic map, which is a local 2d coordinate system along the surface with "
                         "its origin at a given source point.");

      ImGui::Separator();

      ImGui::TextUnformatted("Several variants of the algorithm are available:");

      ImGui::BulletText("VectorHeat");
      ImGui::Indent();
      ImGui::BulletText("the original algorithm from 'The Vector Heat Method'");
      ImGui::BulletText("fast, but may have some distortion artifacts");
      ImGui::Unindent();

      ImGui::BulletText("AffineLocal");
      ImGui::Indent();
      ImGui::BulletText("the fast local algorithm from 'The Affine Heat Method'");
      ImGui::BulletText("fast, highly accurate near source");
      ImGui::Unindent();
      ImGui::BulletText("AffineAdaptive");
      ImGui::Indent();
      ImGui::BulletText("the global algorithm from 'The Affine Heat Method'");
      ImGui::BulletText("highest quality, but slower for repeated solves");
      ImGui::Unindent();

      buildLogmapStrategyMenu();

      if (ImGui::Checkbox("use intrinsic delaunay triangulation", &useIntrinsicTriangulation)) {
        clearCachedSolvers();
      }

      ImGui::Separator();

      if (logMapContinuous) {
        psMesh->setSelectionMode(polyscope::MeshSelectionMode::FacesOnly);

        ImGui::TextWrapped(
            "Hold ctrl (cmd on macOS) and mouse-over the mesh to contniuously recompute the logmap each frame");

        // Get the point on the mesh under the mouse, if there is one
        ImGuiIO& io = ImGui::GetIO();
        if (io.KeyCtrl) {

          if (polyscope::hasPointCloud("logmap source")) {
            // hide it, so we can click-past it when picking below
            polyscope::getPointCloud("logmap source")->setEnabled(false);
          }

          glm::vec2 screenCoords{io.MousePos.x, io.MousePos.y};
          polyscope::PickResult pickResult = polyscope::pickAtScreenCoords(screenCoords);
          if (pickResult.structure == psMesh) {
            polyscope::SurfaceMeshPickResult surfacePickResult = psMesh->interpretPickResult(pickResult);

            if (surfacePickResult.elementType == polyscope::MeshElement::FACE) {
              // this should be the only possibility if there's a hit, we set selection mode to be faces-only

              logmapSourcePoint = SurfacePoint(mesh->face(surfacePickResult.index),
                                               Vector3{surfacePickResult.baryCoords.x, surfacePickResult.baryCoords.y,
                                                       surfacePickResult.baryCoords.z});

              // Visualize the source point
              Vector3 logmapSourcePointPos = logmapSourcePoint.interpolate(geometry->inputVertexPositions);
              auto pointQ = polyscope::registerPointCloud("logmap source", std::vector<Vector3>{logmapSourcePointPos});
              pointQ->setPointRadius(0.005);
            }
          }

          psMesh->setSelectionMode(polyscope::MeshSelectionMode::Auto);

          // un-hide, if we hid it
          if (polyscope::hasPointCloud("logmap source")) {
            // hide it, so we can click-past it when picking below
            polyscope::getPointCloud("logmap source")->setEnabled(true);
          }

          computeLogMapFromSurfacePoint();
        }

      } else {
        psMesh->setSelectionMode(polyscope::MeshSelectionMode::Auto);

        if (ImGui::Button("select logmap source")) {
          long long int pickVert = polyscope::getSurfaceMesh()->selectVertex();
          if (pickVert >= 0) {
            logmapSourceVertex = mesh->vertex(pickVert);
            updateLogmapSourceViz();
          }
        }

        if (ImGui::Button("compute log map")) {
          computeLogMapFromVertex();
        }
      }

      // Checkbox to toggle continuous updates
      ImGui::Checkbox("continuous updates", &logMapContinuous);

      ImGui::Checkbox("limit log map radius", &limitLogMapRadius);
      if (limitLogMapRadius) {
        ImGui::SliderFloat("log map radius", &logMapRadius, 0.0, logMapRadiusUpper);
      }

      ImGui::EndTabItem();
    }
    if (ImGui::BeginTabItem("Centers")) {
      newUImode = "centers";
      psMesh->setSelectionMode(polyscope::MeshSelectionMode::Auto);

      ImGui::TextWrapped("Compute 'center' of a set of points on a surface, also known as a Karcher or Frechet mean. "
                         "The approach is an iterative algorithm built on top of the logarithmic map.");

      buildLogmapStrategyMenu();

      if (ImGui::Checkbox("use intrinsic delaunay triangulation", &useIntrinsicTriangulation)) {
        clearCachedSolvers();
      }

      if (ImGui::TreeNode("select sites to compute center of")) {
        buildSitesMenu();
        ImGui::TreePop();
      }


      ImGui::PushItemWidth(200);
      ImGui::InputInt("p norm", &pCenter);
      ImGui::PopItemWidth();

      if (ImGui::Button("find center")) {
        computeCenter();
      }

      ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
  }

  if (newUImode != UImode) {
    UImode = newUImode;

    // clear out old viz
    polyscope::removeStructure("source points", false);
    polyscope::removeStructure("logmap source", false);
    polyscope::removeStructure("site points", false);
    polyscope::removeStructure("center", false);
    psMesh->removeAllQuantities();

    if (UImode == "extension") {
      updateSourceSetViz();
    } else if (UImode == "logmap") {
      updateLogmapSourceViz();
    } else if (UImode == "centers") {
      updateSiteSetViz();
    }
  }
}

int main(int argc, char** argv) {

  // Configure the argument parser
  args::ArgumentParser parser("A demo of the Vector Heat Method");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help&) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError& e) {
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

  // Some options
  polyscope::options::giveFocusOnShow = true;
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::TileReflection;

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // Load mesh
  std::tie(mesh, geometry) = loadMesh(args::get(inputFilename));

  // Register the mesh with polyscope
  psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(args::get(inputFilename)),
                                          geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                          polyscopePermutations(*mesh));


  geometry->requireVertexIndices();
  geometry->requireVertexTangentBasis();

  // To start, pick two vertices as sources
  addVertexSource(0);
  addVertexSource(mesh->nVertices() / 2);
  extensionSourcePoints[1].scalarVal = 3.0;


  addVertexSite(0);
  addVertexSite(2 * mesh->nVertices() / 3);
  addVertexSite(mesh->nVertices() / 3);

  logmapSourceVertex = mesh->vertex(0);

  // Diameter estimate
  geometry->requireFaceAreas();
  geometry->requireFaceTangentBasis();
  double surfaceArea = 0.;
  for (Face f : mesh->faces()) {
    surfaceArea += geometry->faceAreas[f];
  }
  logMapRadiusUpper = std::sqrt(surfaceArea);
  logMapRadius = logMapRadiusUpper / 4.;

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
