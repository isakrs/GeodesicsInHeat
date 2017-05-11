#include "MeshEdge.hpp"
#include "MeshFace.hpp"
#include "MeshVertex.hpp"
#include "Mesh.hpp"
#include "DGP/Vector4.hpp"

# include <cmath>

MeshEdge *
MeshEdge::nextAroundEndpoint(int i)
{
  debugAssertM(i == 0 || i == 1, "MeshEdge: Invalid endpoint index");

  if (numFaces() > 2)  // non-manifold
    return NULL;

  // Find which incident face has this endpoint as the origin of the edge when stepping round the face. The required edge
  // is then the predecessor of this edge around the face.
  for (FaceIterator fi = facesBegin(); fi != facesEnd(); ++fi)
  {
    Face * face = *fi;
    MeshEdge * prev = face->getPredecessor(this);
    if (prev->hasEndpoint(endpoints[i]))  // found it!
      return prev;
  }

  return NULL;
}

double
MeshEdge::sumCotans() const
{
	double cotans = 0.0;
	Vector3 p0 = this->getEndpoint(0)->getPosition();
	Vector3 p1 = this->getEndpoint(1)->getPosition();
	for (FaceConstIterator fci = facesBegin(); fci != facesEnd(); ++fci)
	{
		for (Face::VertexConstIterator fvci = (*fci)->verticesBegin(); fvci != (*fci)->verticesEnd(); ++fvci)
		{
			Vector3 p = (*fvci)->getPosition();
			if (p != p0 && p != p1)
			{	// p is p2, last corner of triangulated face.
				// making right angle triangle.
				// a is projection of h on p to p1, correctly scaled.
				Vector3 v1 = p1 - p;
				Vector3 v2 = p0 - p;
				double cotan = v1.dot(v2) / v1.cross(v2).length();
				cotans += cotan;
			}
		}
	}
	return cotans;
}

double
MeshEdge::sumAreas() const
{
	double areas = 0.0;
	Vector3 p0 = this->getEndpoint(0)->getPosition();
	Vector3 p1 = this->getEndpoint(1)->getPosition();
	for (FaceConstIterator fci = facesBegin(); fci != facesEnd(); ++fci)
	{
		for (Face::VertexConstIterator fvci = (*fci)->verticesBegin(); fvci != (*fci)->verticesEnd(); ++fvci)
		{
			Vector3 p = (*fvci)->getPosition();
			if (p != p0 && p != p1)
			{	// p is p2, last corner of triangulated face.
				Vector3 v1 = p1 - p;
				Vector3 v2 = p0 - p;
				double area = 0.5 * v1.cross(v2).length();
				areas += area;
			}
		}
	}
	return areas;
}

double 
MeshEdge::getLength() const 
{ 
	Vector3 p0 = this->getEndpoint(0)->getPosition();
	Vector3 p1 = this->getEndpoint(1)->getPosition();
    return (p1-p0).length(); 
}
