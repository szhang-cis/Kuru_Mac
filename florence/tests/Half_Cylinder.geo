// Gmsh project created on Wed Feb 13 16:44:14 2019
Point(1) = {9.295, 0, 0, 1.0};
Point(2) = {10.705, 0, 0, 1.0};
Line(1) = {1, 2};
Extrude {0, 100.0, 0} {
  Line{1};
}
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Surface{5};Layers{13};Recombine;
}
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Surface{27};Layers{13};Recombine;
}
//Thickness divisions
Transfinite Line {1,7,29,2,9,31} = 2 Using Progression 1;
//Axial divisions and Recombine
Transfinite Line {3,4,8,10,30,32} = 30 Using Progression 1;
Transfinite Surface {5,27,49};
Recombine Surface {5,27,49};
// Symmetry
Physical Surface(98) = {5, 49};
// Top
Physical Surface(93) = {44, 22};
// Bottom
Physical Surface(94) = {36, 14};
// Inner wall
Physical Surface(95) = {26, 48};
// Outer wall
Physical Surface(96) = {40, 18};
// Volume of Cylinder
Physical Volume(97) = {2, 1};
// Generate 3D mesh
Mesh 3;
// Remove duplicate entities
Coherence Mesh;
