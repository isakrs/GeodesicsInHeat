#include "Mesh.hpp"
#include "Viewer.hpp"
#include <cstdlib>

int
usage(int argc, char * argv[])
{
  DGP_CONSOLE << "";
  DGP_CONSOLE << "Usage: " << argv[0] << " <mesh-in> [<mesh-out>]";
  DGP_CONSOLE << "";

  return -1;
}

int
main(int argc, char * argv[])
{
  if (argc < 2)
    return usage(argc, argv);

  std::string in_path = argv[1];

  std::string out_path;
  if (argc >= 3)
  {
      out_path = argv[2];
  }

  Mesh mesh;
  if (!mesh.load(in_path))
    return -1;

  DGP_CONSOLE << "Read mesh '" << mesh.getName() << "' with " << mesh.numVertices() << " vertices, " << mesh.numEdges()
              << " edges and " << mesh.numFaces() << " faces from " << in_path;

  if (!out_path.empty())
  {
    if (!mesh.save(out_path))
      return -1;

    DGP_CONSOLE << "Saved mesh to " << out_path;
  }
	
  mesh.setup();
    
  std::cout << "Press: g in window to place heat point, i.e. point to find distances from." << std::endl;
  
  Viewer viewer;
  viewer.setObject(&mesh);
  viewer.launch(argc, argv);
    
  return 0;
}
