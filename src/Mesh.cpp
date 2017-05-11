#include "Mesh.hpp"
#include "MeshVertex.hpp"
#include "MeshEdge.hpp"
#include "MeshFace.hpp"
#include "DGP/FilePath.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <unordered_map>
#include <limits>

void
Mesh::draw(Graphics::RenderSystem & render_system, bool draw_edges, bool use_vertex_data, bool send_colors) const
{
  // Three separate passes over the faces is probably faster than using Primitive::POLYGON for each face

  if (draw_edges)
  {
    render_system.pushShapeFlags();
    render_system.setPolygonOffset(true, 1);
  }

  // First try to render as much stuff using triangles as possible
  render_system.beginPrimitive(Graphics::RenderSystem::Primitive::TRIANGLES);
    for (FaceConstIterator fi = facesBegin(); fi != facesEnd(); ++fi)
      if (fi->isTriangle()) drawFace(*fi, render_system, use_vertex_data, send_colors);
  render_system.endPrimitive();

  // Now render all quads
  render_system.beginPrimitive(Graphics::RenderSystem::Primitive::QUADS);
    for (FaceConstIterator fi = facesBegin(); fi != facesEnd(); ++fi)
      if (fi->isQuad()) drawFace(*fi, render_system, use_vertex_data, send_colors);
  render_system.endPrimitive();

  // Finish off with all larger polygons
  for (FaceConstIterator fi = facesBegin(); fi != facesEnd(); ++fi)
    if (fi->numEdges() > 4)
    {
      render_system.beginPrimitive(Graphics::RenderSystem::Primitive::POLYGON);
        drawFace(*fi, render_system, use_vertex_data, send_colors);
      render_system.endPrimitive();
    }

  if (draw_edges)
    render_system.popShapeFlags();

  if (draw_edges)
  {
    render_system.pushShader();
    render_system.pushColorFlags();

      render_system.setShader(NULL);

      render_system.beginPrimitive(Graphics::RenderSystem::Primitive::LINES);
        for (EdgeConstIterator ei = edgesBegin(); ei != edgesEnd(); ++ei)
        {
	       double p = (ei->getEndpoint(0)->getPhi() + ei->getEndpoint(1)->getPhi())/2.0;

          //render_system.setColor(ColorRGBA(0.2, 0.3, 0.7, 1));  // set default edge color
          render_system.setColor(ColorRGBA(p, 0, 1-p, 1));
          render_system.sendVertex(ei->getEndpoint(0)->getPosition());
          render_system.sendVertex(ei->getEndpoint(1)->getPosition());
        }
      render_system.endPrimitive();

    render_system.popColorFlags();
    render_system.popShader();
  }
}

bool
Mesh::loadOFF(std::string const & path)
{
  std::ifstream in(path.c_str());
  if (!in)
  {
    DGP_ERROR << "Could not open '" << path << "' for reading";
    return false;
  }

  clear();

  std::string magic;
  if (!(in >> magic) || magic != "OFF")
  {
    DGP_ERROR << "Header string OFF not found at beginning of file '" << path << '\'';
    return false;
  }

  long nv, nf, ne;
  if (!(in >> nv >> nf >> ne))
  {
    DGP_ERROR << "Could not read element counts from OFF file '" << path << '\'';
    return false;
  }

  if (nv < 0 || nf < 0 || ne < 0)
  {
    DGP_ERROR << "Negative element count in OFF file '" << path << '\'';
    return false;
  }

  std::vector<Vertex *> indexed_vertices;
  Vector3 p;
  for (long i = 0; i < nv; ++i)
  {
    if (!(in >> p[0] >> p[1] >> p[2]))
    {
      DGP_ERROR << "Could not read vertex " << indexed_vertices.size() << " from '" << path << '\'';
      return false;
    }

    Vertex * v = addVertex(p);
    if (!v)
      return false;

    indexed_vertices.push_back(v);
  }

  std::vector<Vertex *> face_vertices;
  long num_face_vertices, vertex_index;
  for (long i = 0; i < nf; ++i)
  {
    if (!(in >> num_face_vertices) || num_face_vertices < 0)
    {
      DGP_ERROR << "Could not read valid vertex count of face " << faces.size() << " from '" << path << '\'';
      return false;
    }

    face_vertices.resize(num_face_vertices);
    for (size_t j = 0; j < face_vertices.size(); ++j)
    {
      if (!(in >> vertex_index))
      {
        DGP_ERROR << "Could not read vertex " << j << " of face " << faces.size() << " from '" << path << '\'';
        return false;
      }

      if (vertex_index < 0 || vertex_index >= (long)vertices.size())
      {
        DGP_ERROR << "Out-of-bounds index " << vertex_index << " of vertex " << j << " of face " << faces.size() << " from '"
                  << path << '\'';
        return false;
      }

      face_vertices[j] = indexed_vertices[(size_t)vertex_index];
    }

    addFace(face_vertices.begin(), face_vertices.end());  // ok if this fails, just skip the face with a warning
  }

  setName(FilePath::objectName(path));

  return true;
}

bool
Mesh::saveOFF(std::string const & path) const
{
  std::ofstream out(path.c_str(), std::ios::binary);
  if (!out)
  {
    DGP_ERROR << "Could not open '" << path << "' for writing";
    return false;
  }

  out << "OFF\n";
  out << numVertices() << ' ' << numFaces() << " 0\n";

  typedef std::unordered_map<Vertex const *, long> VertexIndexMap;
  VertexIndexMap vertex_indices;
  long index = 0;
  for (VertexConstIterator vi = vertices.begin(); vi != vertices.end(); ++vi, ++index)
  {
    Vector3 const & p = vi->getPosition();
    out << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';

    vertex_indices[&(*vi)] = index;
  }

  for (FaceConstIterator fi = faces.begin(); fi != faces.end(); ++fi)
  {
    out << fi->numVertices();

    for (Face::VertexConstIterator vi = fi->verticesBegin(); vi != fi->verticesEnd(); ++vi)
    {
      VertexIndexMap::const_iterator existing = vertex_indices.find(*vi);
      if (existing == vertex_indices.end())
      {
        DGP_ERROR << "Face references vertex absent from mesh '" << path << '\'';
        return false;
      }

      out << ' ' << existing->second;
    }

    out << '\n';
  }

  return true;
}

bool
Mesh::load(std::string const & path)
{
  std::string path_lc = toLower(path);
  bool status = false;
  if (endsWith(path_lc, ".off"))
    status = loadOFF(path);
  else
  {
    DGP_ERROR << "Unsupported mesh format: " << path;
  }
  
  return status;
}

bool
Mesh::save(std::string const & path) const
{
  std::string path_lc = toLower(path);
  if (endsWith(path_lc, ".off"))
    return saveOFF(path);

  DGP_ERROR << "Unsupported mesh format: " << path;
  return false;
}

// ------------------------------------------------------
// --- Geodesic Heat Method methods ---
void 
Mesh::buildCotanOperator(Eigen::SparseMatrix<double>& Lc) const
{
  std::vector<Eigen::Triplet<double>> LcTriplet;
    
  for (VertexConstIterator vci = verticesBegin(); vci != verticesEnd(); vci++) 
  {
  	double sumCoefficients = 0.0;
    for (MeshVertex::EdgeConstIterator eci = vci->edgesBegin(); eci != vci->edgesEnd(); ++eci)
    {
      // (cotA + cotB)
      double coefficient = (*eci)->sumCotans();
      sumCoefficients += coefficient;
            
      Vertex const * otherVertex = (*eci)->getOtherEndpoint(&(*vci));
      LcTriplet.push_back(Eigen::Triplet<double>(vci->getIndex(), otherVertex->getIndex(), coefficient));
    }      	
    
    LcTriplet.push_back(Eigen::Triplet<double>(vci->getIndex(), vci->getIndex(), -sumCoefficients));	    
  
  }        

  Lc.setFromTriplets(LcTriplet.begin(), LcTriplet.end());

}

void 
Mesh::buildAreaMatrix(Eigen::SparseMatrix<double>& A) const
{
  std::vector<Eigen::Triplet<double>> ATriplet;
    
  for (VertexConstIterator vci = verticesBegin(); vci != verticesEnd(); vci++) 
  {
	double areas = 0.0;
	for (MeshVertex::EdgeConstIterator veci = vci->edgesBegin(); veci != vci->edgesEnd(); ++veci)
	{
	  // multipling with 0.5, because same area is calculated twice
	  double area = 0.5 * (*veci)->sumAreas();
	  areas += area;
	}
		
	// one third of the neighbouring areas
	areas *= (1.0/3.0);
		
    ATriplet.push_back(Eigen::Triplet<double>(vci->getIndex(), vci->getIndex(), areas));
  }
    
  A.setFromTriplets(ATriplet.begin(), ATriplet.end());
}

double 
Mesh::computeTimeStep() const
{
  // t = avg edge length ^ 2
  double avgLength = 0.0;
  for (EdgeConstIterator eci = edgesBegin(); eci != edgesEnd(); eci++) 
  {
    avgLength += eci->getLength();
  }    
  avgLength /= (double)edges.size();  
  return avgLength * avgLength;
}

void 
Mesh::computeFaceGradients(Eigen::MatrixXd& gradients, const Eigen::VectorXd& u) const
{
  for (FaceConstIterator fci = faces.begin(); fci != faces.end(); ++fci) 
  {
    Vector3 gradient(0,0,0);
    Vector3 normal = fci->getNormal();
    normal /= normal.length();   // assuring length 1
    
    for (MeshFace::VertexConstIterator vci = fci->verticesBegin(); vci != fci->verticesEnd(); ++vci)
    { // vci get hold ui
      
      for (MeshFace::EdgeConstIterator eci = fci->edgesBegin(); eci != fci->edgesEnd(); ++eci)
      { // finding opposite edge, making sure correct direction
        if (!(*vci)->hasIncidentEdge((*eci)))
        {
          double ui = u((*vci)->getIndex());
  		  Vector3 p0 = (*eci)->getEndpoint(0)->getPosition();
		  Vector3 p1 = (*eci)->getEndpoint(1)->getPosition();
          Vector3 ed = p0 - p1;
          
          if ((*eci)->getEndpoint(1) == fci->getPredecessor(vci))  
            ed = -ed;        // maintain counter-clockwise direction
    	  
    	  gradient += ui * normal.cross(ed);
        }
      }
      
  	}
  	
  	gradient /= (2.0 * fci->area());
  	gradient /= gradient.length();
  	Eigen::Vector3d grad(gradient.x(),gradient.y(),gradient.z());
  	gradients.row(fci->getIndex()) = -grad;
  
  }
}

void 
Mesh::computeIntegratedDivergence(Eigen::VectorXd& integratedDivs,
                            const Eigen::MatrixXd& gradients) const
{
  for (VertexConstIterator vci = vertices.begin(); vci != vertices.end(); ++vci) 
  {
    double integratedDiv = 0.0;
    Vector3 p = vci->getPosition();
 		
 	for (MeshVertex::FaceConstIterator fci = vci->facesBegin(); fci != vci->facesEnd(); ++fci)
    {
      Eigen::Vector3d gradient = gradients.row((*fci)->getIndex());
            
      for(MeshFace::VertexConstIterator fvi = (*fci)->verticesBegin(); fvi != (*fci)->verticesEnd(); ++fvi)
      {
 		if(*fvi == &(*vci))
        {
          ++fvi;
          if(fvi == (*fci)->verticesEnd()) 
            fvi = (*fci)->verticesBegin();
          Vector3 p1 = (*fvi)->getPosition();
                   
		  ++fvi;
          if(fvi == (*fci)->verticesEnd()) 
            fvi = (*fci)->verticesBegin(); 
          Vector3 p2 = (*fvi)->getPosition(); 
         
          Vector3 e1 = p1 - p;
          Vector3 e2 = p2 - p;
          Vector3 ei = p2 - p1;
                  
		  double theta1 = acos((-e2).dot(-ei) / (e2.length() * ei.length()));
          double cot1 = 1.0 / tan(theta1);
                  
          double theta2 = acos((-e1).dot(ei) / (e1.length() * ei.length()));
          double cot2 = 1.0 / tan(theta2);
		  
		  Vector3 grad(gradient[0],gradient[1],gradient[2]);
                  
		  integratedDiv += e1.dot(grad) * cot1 + e2.dot(grad) * cot2;
		  break;
        }
      }        
    } 
    integratedDivs[vci->getIndex()] = 0.5 * integratedDiv;
  }
}

void 
Mesh::setup()
{
  // set id to vertices
  long idv = 0;
  for (VertexIterator vi = verticesBegin(); vi != verticesEnd(); ++vi)
  {
	vi->setIndex(idv);
	++idv;
  }
	
  // set id to faces
  long idf = 0;
  for (FaceIterator fi = facesBegin(); fi != facesEnd(); ++fi)
  {
	fi->setIndex(idf);
	++idf;
  }
	
  long nv = vertices.size();
    
  // build Cotan Operator
  Eigen::SparseMatrix<double> Lc(nv, nv);
  buildCotanOperator(Lc);
  poissonSolver.compute(Lc);
    
  // build diagonal area matrix
  Eigen::SparseMatrix<double> A(nv, nv);
  buildAreaMatrix(A);
    
  // compute time step
  double t = computeTimeStep();
    
  // 1. compute heat flow for time t
  heatSolver.compute(A - t*Lc);
}

void
Mesh::computeGeodesics()
{
  int nv = this->vertices.size();		// number of vertices, nv
  std::srand(time(NULL));
  long npti = rand() % nv;				// needle point index, npti for heat
  needlePointIndex = npti;	
	
  // Step 1: Set energy u for all vertices.
  Eigen::VectorXd u(nv);
  u.setZero();					
  u(npti) = 1;						// Heated needle is placed
	
  u = heatSolver.solve(u);			// Heat is spread

  // 2. evaluate face gradients
  Eigen::MatrixXd gradients((int)faces.size(), 3);
  computeFaceGradients(gradients, u);

  // 3. solve poisson equation
  Eigen::VectorXd integratedDivs(nv);
  computeIntegratedDivergence(integratedDivs, gradients);        
     
  Eigen::VectorXd phis = poissonSolver.solve(integratedDivs);

  // compute max and min phis
  double minPhi = INFINITY;
  double maxPhi = -INFINITY;
  for (VertexIterator vi = verticesBegin(); vi != verticesEnd(); vi++)  
  {
    if (minPhi > phis(vi->getIndex())) minPhi = phis(vi->getIndex());
    if (maxPhi < phis(vi->getIndex())) maxPhi = phis(vi->getIndex());
  }
  double range = maxPhi - minPhi;  
   	
  // set phis in range [0,1]
  for (VertexIterator vi = verticesBegin(); vi != verticesEnd(); vi++) 
    vi->setPhi((1 - ((phis(vi->getIndex()) - minPhi) / range)));                                                           
}