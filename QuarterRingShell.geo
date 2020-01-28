// Gmsh project created on Wed Feb 13 16:44:14 2019
Point(1) = {0.010, 0, 0, 1.0};
Point(2) = {0.010, 0.002, 0, 1.0};
Line(1) = {1, 2};
// Circular extrusion of line 1
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Line{1};Layers{12};Recombine;
}
//Axial divisions and Recombine
Transfinite Line {1,2} = 2 Using Progression 1;
// Symmetry Z
Physical Line(11) = {1};
// Symmetry X
Physical Line(12) = {2};
// Bottom
Physical Line(13) = {3};
// Top
Physical Line(14) = {4};
// Wall (Inner)
Physical Surface(15) = {5};
// Generate 3D mesh
Mesh 3;
// Remove duplicate entities
Coherence Mesh;
