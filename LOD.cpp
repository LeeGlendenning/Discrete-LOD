#include <iostream>
#include <vector>
#include <math.h>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

using namespace std;

struct MyTraits : public OpenMesh::DefaultTraits
{
	VertexAttributes(OpenMesh::Attributes::Status);
	FaceAttributes(OpenMesh::Attributes::Status);
	EdgeAttributes(OpenMesh::Attributes::Status);
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  MyMesh;

struct Error {
	double error;
	MyMesh::VertexHandle vertex;
	int index;
	vector <MyMesh::VertexHandle> adjacent_verts;
};

struct Node {
	MyMesh::VertexHandle vertex;
	vector <Node> children;
};

// Shallow tree consisting of a root and a set of children that are leaves
struct LODTree {
	Node root;
	vector <MyMesh::FaceHandle> temp_triangulation;
};

// Each level contains a shallow tree containing removed vertices and their adjacent vertices (children)
struct LODForest {
	vector <vector <LODTree>> levels;
	int cur_level = -1; // Incremented every time reduce_LOD is called
};

MyMesh  mesh;
LODForest forest;

// Counts and returns total number of vertices in mesh
int countVertices() {

	int totalVerts = 0;
	// Loop through and count all vertices in mesh *v_it
	for (MyMesh::VertexIter vertex_iter = mesh.vertices_begin(); vertex_iter != mesh.vertices_end(); ++vertex_iter) {
		totalVerts++;
	}
	
	return totalVerts;
}


vector<MyMesh::VertexHandle> get_adjacent_vertices(MyMesh::VertexHandle v) {
	vector <MyMesh::VertexHandle> adjacent_verts;
	// Store all adjacent vertices
	for (MyMesh::VertexVertexIter vv_it = mesh.vv_iter(v); vv_it.is_valid(); ++vv_it)
	{
		adjacent_verts.push_back(*vv_it);
	}

	return adjacent_verts;
}

OpenMesh::Vec3f findBarycenter(vector <MyMesh::VertexHandle> points) {
	//cout << "Finding barycenter of " << points[0] << ", ..." << endl;
	OpenMesh::Vec3f barycenter = { 0, 0, 0 };

	int i;
	for (i = 0; i < points.size(); i++)
	{
		barycenter[0] += mesh.point(points[i])[0];
		barycenter[1] += mesh.point(points[i])[1];
		barycenter[2] += mesh.point(points[i])[2];
	}

	barycenter[0] /= points.size();
	barycenter[1] /= points.size();
	barycenter[2] /= points.size();

	return barycenter;
}

OpenMesh::Vec3f vecSubtract(OpenMesh::Vec3f v1, OpenMesh::Vec3f v2){
	OpenMesh::Vec3f v;

	v[0] = v1[0] - v2[0]; // x
	v[1] = v1[1] - v2[1]; // y
	v[2] = v1[2] - v2[2]; // z

	return v;
}

OpenMesh::Vec3f findCrossProd(MyMesh::Point p1, MyMesh::Point p2, MyMesh::Point p3) {
	// e1 = (p1-p3), e2 = (p2-p3)
	OpenMesh::Vec3f e1 = vecSubtract(p1, p3), e2 = vecSubtract(p2, p3), cross_prod;

	// Find cross product: e1 x e2
	cross_prod[0] = e1[1] * e2[2] - e1[2] * e2[1]; // x
	cross_prod[1] = e1[2] * e2[0] - e1[0] * e2[2]; // y
	cross_prod[2] = e1[0] * e2[1] - e1[1] * e2[0]; // z

	//cout << "e1: " << e1[0] << ", " << e1[1] << ", " << e1[2] << endl;
	//cout << "e2: " << e2[0] << ", " << e2[1] << ", " << e2[2] << endl;
	//cout << "cross_prod: " << cross_prod << endl;

	return cross_prod;
}

double findDotProd(MyMesh::Point p1, MyMesh::Point p2) {
	//cout << "dotProd = " << (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) << endl;
	return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
}

// V = |(p1-p4) <dot> ((p2-p4) <cross> (p3-p4))| / 6
double findTetrahedronVolume(MyMesh::Point p1, MyMesh::Point p2, MyMesh::Point p3, MyMesh::Point p4) {
	return abs(findDotProd(vecSubtract(p1, p4), findCrossProd(p2, p3, p4))) / 6;
}

// Given a vertex, find volume lost if it were to be removed - called the error
double findError(MyMesh::VertexHandle v) {
	double error = 0.0;
	vector <MyMesh::VertexHandle> adjacent_verts = get_adjacent_vertices(v);
	MyMesh::Point ref_point;

	// Projection of *vertex_it onto mesh if *vertex_it were to be removed
	ref_point = findBarycenter(adjacent_verts);

	// Calculate error contribution for every adjacent vertex
	MyMesh::Point prev = mesh.point(adjacent_verts[adjacent_verts.size() - 1]); // Last adjacent vertex
	MyMesh::Point curr = mesh.point(adjacent_verts[0]); // First adjacent vertex
	
	error += findTetrahedronVolume(prev, curr, ref_point, mesh.point(v));

	int i;
	for (i = 1; i < adjacent_verts.size(); i++) {
		prev = curr;
		curr = mesh.point(adjacent_verts[i]);
		error += findTetrahedronVolume(prev, curr, ref_point, mesh.point(v));
	}
	
	return error;
}

Error calc_new_worst_error(Error* errors, int size_errors, Error worst_error) {
	worst_error.error = -1; // Clear old error value
	int i;
	for (i = 0; i < size_errors; i++) {
		if (errors[i].error > worst_error.error) {
			worst_error.error = errors[i].error;
			worst_error.index = errors[i].index;
			worst_error.vertex = errors[i].vertex;
		}
	}
	
	return worst_error;
}

// Given a percentage of vertices to remove, find the ones with minimum errors
Error* findMinErrors(int size_errors, int totalVerts) {
	
	Error* errors = new Error[size_errors];
	Error worst_error; // Error stored in list having the highest value. Used to avoid comparisons
	worst_error.error = -1;
	cout << "# vertices to remove = " << size_errors << endl;

	int count = 0;

	// Loop through and count all vertices in mesh, find minimum errors
	for (MyMesh::VertexIter vertex_iter = mesh.vertices_begin(); vertex_iter != mesh.vertices_end(); ++vertex_iter) {

		// Do not remove boundaries
		if (!mesh.is_boundary(*vertex_iter)) {

			// Determine if new error is less than errors in list
			if (count >= size_errors)
			{
				double temp_error = findError(*vertex_iter);
				
				if (temp_error < worst_error.error)
				{
					// Create new Error structure and replace worst_error with it in list
					Error new_error;
					new_error.vertex = *vertex_iter;
					new_error.error = temp_error;
					new_error.index = worst_error.index;
					errors[new_error.index] = new_error;

					worst_error = calc_new_worst_error(errors, size_errors, worst_error);
				}
			}
			else // Just put error in list
			{
				// Create a new Error structure
				Error new_error;
				new_error.vertex = *vertex_iter;
				new_error.error = findError(*vertex_iter);
				new_error.index = count;

				// Update worst_error
				if (new_error.error > worst_error.error) {
					worst_error.error = new_error.error;
					worst_error.index = count;
					worst_error.vertex = *vertex_iter;
				}
				errors[new_error.index] = new_error;
			}

			count++;
		}
	}
	
	return errors;
}

LODTree remove_vertex(Error e) {
	LODTree tree;
	
	//Create LODTree
	Node root;
	root.vertex = e.vertex;

	int i;
	// Add children (adjacent vertices) to root
	for (i = 0; i < e.adjacent_verts.size(); i++) {
		Node child;
		child.vertex = e.adjacent_verts[i];
		root.children.push_back(child);
	}

	tree.root = root;
	forest.levels[forest.cur_level].push_back(tree);// Add LODTree to forest

	// Actually remove vertex from mesh
	mesh.delete_vertex(e.vertex, false);

	return tree;
}

// Fill empty space created when removing vertex by arbitrarily triangulating it
LODTree retriangulate(LODTree tree) {

	// Route stderr to a file because OpenMesh prints an error every time an invalid face is created
	// but does not provide a way to check if a face is valid before trying to create it so to re-triangulate
	// we must create the face and then check if it is valid or not.
	// This error message is printed to the screen by default and happens a lot of times for a complex mesh.
	std::ofstream error("error.txt");
	std::streambuf *errbuf = std::cerr.rdbuf(error.rdbuf());

	int i;
	// Here we use a naive retriangulation method. Since faces require a counter clockwise ordering of vertices,
	// and we do not know the order of the tree's children and the set may not be convex, we must try all 
	// combinations of faces (all subsets of size 3) to triangulate the set
	for (i = 0; i < tree.root.children.size()-2; i++) {
		if (tree.root.children.size() > 1) {
			
			MyMesh::FaceHandle f;
			int j;
			for (j = i + 1; j < tree.root.children.size() - 1; j++) 
			{
				int k;
				for (k = j + 1; k < tree.root.children.size(); k++) 
				{
					f = mesh.add_face({ tree.root.children[i].vertex, tree.root.children[k].vertex, tree.root.children[j].vertex });
					if (f.is_valid()) {
						//cout << "Success" << endl;
						tree.temp_triangulation.push_back(f); // Store the triangulation so it can be removed when increasing LOD
					}
				}
			}

		}
		else { //else it should have been caught as a boundary vertex but didn't. Don't think this enters anymore
			break;
		}
		
	}

	std::cerr.rdbuf(errbuf); // Re-route stderr to the screen

	return tree;
}

// Add new level to the LODForest
void add_forest_level() {
	vector <LODTree> level;
	forest.levels.push_back(level);
	forest.cur_level++;
}

void remove_temp_triangulation(LODTree tree) {
	int i;
	// Remove old triangulation of face
	for (i = 0; i < tree.temp_triangulation.size(); i++) {
		if (!mesh.status(tree.temp_triangulation[i]).deleted())
			mesh.delete_face(tree.temp_triangulation[i], false);
	}
}


void reduce_LOD(double percent_reduction) {
	add_forest_level();
	
	int totalVerts = countVertices();

	int size_errors = percent_reduction * totalVerts;

	Error* min_errors = findMinErrors(size_errors, totalVerts);

	int i;
	// Reduce LOD by removing the "percent_reduction" amount of vertices found to have the lowest error
	for (i = 0; i < size_errors; i++) {

		// Store all adjacent vertices of vertices to be deleted
		min_errors[i].adjacent_verts = get_adjacent_vertices(min_errors[i].vertex);
			
		LODTree tree = remove_vertex(min_errors[i]); // Maintains the LODForest
		
		forest.levels[forest.cur_level][i] = retriangulate(tree); // Fix the hole created by removing vertex
		
		// Invalidates all handles except those passed to the function
		// First, this is very slow to garbage collect after each face triangulation
		// Second, this won't work because some faces depend on vertices that may have already been deleted.
		//mesh.garbage_collection(&forest.levels[forest.cur_level][i].temp_triangulation);
	}
	
	mesh.garbage_collection();
}

// Connect all children of vert to it
void add_tree_to_mesh(LODTree tree) {
	int i;
	// Connect new vertex to its children
	for (i = 0; i < tree.root.children.size()-1; i++) {
		mesh.add_face({ tree.root.vertex, tree.root.children[i + 1].vertex, tree.root.children[i].vertex });
	}
	
}

// Increases LOD to the next level (however much % was specified when reducing)
void increase_LOD() {
	// If not already at highest LOD
	if (forest.cur_level >= 0) {
		int i;
		
		cout << "# trees: " << forest.levels[forest.cur_level].size() << endl;
		// Remove the temporary triangulation that replaced the removed vertex
		// This will leave a hole
		for (i = 0; i < forest.levels[forest.cur_level].size(); i++) {
			if (i < 10)
				remove_temp_triangulation(forest.levels[forest.cur_level][i]);
		}
		mesh.garbage_collection();

		// Add all vertices (defined as LODTrees) that have been removed during previous LOD reduction
		// This will fill the hole
		for (i = 0; i < forest.levels[forest.cur_level].size(); i++) {
			// Add previously removed vertex to mesh
			mesh.add_vertex(mesh.point(forest.levels[forest.cur_level][i].root.vertex));
			// Connect vertex by edges (and remove old triangulation of the space)
			add_tree_to_mesh(forest.levels[forest.cur_level][i]);
		}

		forest.cur_level--;
	}
}

int main(int argc, char **argv)
{
	// Check command line options
	if (argc != 3)
	{
		std::cerr << "Usage:  " << argv[0] << " <input_file> <output_file>\n";
		return 1;
	}
	// Read mesh from input_file argv[1]
	if (!OpenMesh::IO::read_mesh(mesh, argv[1]))
	{
		std::cerr << "Error: Cannot read mesh from " << argv[1] << std::endl;
		return 1;
	}

	// Reduce LOD by 5%
	reduce_LOD(0.10);
	// Reduce LOD by 5%
	reduce_LOD(0.15);
	// Reduce LOD by 5%
	reduce_LOD(0.20);
	// Reduce LOD by 5%
	reduce_LOD(0.25);
	// Reduce LOD by 5%
	reduce_LOD(0.30);
	// Reduce LOD by 5%
	reduce_LOD(0.35);


	// Write mesh after LOD is decreased
	if (!OpenMesh::IO::write_mesh(mesh, "output_decrease.obj")) {
		std::cerr << "Error: cannot write mesh to " << "output_decrease.obj" << std::endl;
		return 1;
	}

	/*increase_LOD();

	// Write mesh after LOD is increased
	if (!OpenMesh::IO::write_mesh(mesh, "output_increase.obj")) {
		std::cerr << "Error: cannot write mesh to " << "output_increase.obj" << std::endl;
		return 1;
	}*/

	return 0;
}
