// Gmsh project created on Wed Feb 13 16:44:14 2019
//Point(1) = {0.009295, 0, 0, 1.0};
//Point(2) = {0.010705, 0, 0, 1.0};
Point(1) = {0.010, 0, 0, 1.0};
Point(2) = {0.01141, 0, 0, 1.0};
Line(1) = {1, 2};
Extrude {0, 0.001, 0} {
  Line{1};
}
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Surface{5};Layers{12};Recombine;
}
//Thickness divisions
Transfinite Line {1,2,7,9} = 2 Using Progression 1;
//Axial divisions and Recombine
Transfinite Line {3,4,8,10} = 2 Using Progression 1;
Transfinite Surface {5,27};
Recombine Surface {5,27};
// Symmetry Z
Physical Surface(98) = {5};
// Symmetry X
Physical Surface(99) = {27};
// Top
Physical Surface(93) = {22};
// Bottom
Physical Surface(94) = {14};
// Inner wall
Physical Surface(95) = {26};
// Outer wall
Physical Surface(96) = {18};
// Volume of Cylinder
Physical Volume(97) = {1};
// Generate 3D mesh
Mesh 3;
// Remove duplicate entities
Coherence Mesh;
