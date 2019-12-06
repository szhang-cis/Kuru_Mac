// Gmsh project created on Tue Jan 22 18:33:14 2019
Point(1) = {9.295, 0, 0, 1.0};
Point(2) = {10.705, 0, 0, 1.0};
Line(1) = {1, 2};
Extrude {0, 90.0, 0} {
  Line{1};Layers{18};Recombine;
}
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Surface{5};Layers{10};Recombine;
}
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Surface{27};Layers{10};Recombine;
}
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Surface{49};Layers{10};Recombine;
}
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Surface{71};Layers{10};Recombine;
}
// Top
Physical Surface(93) = {66, 88, 22, 44};
// Bottom
Physical Surface(94) = {58, 80, 14, 36};
// Inner wall
Physical Surface(95) = {70, 48, 26, 92};
// Outer wall
Physical Surface(96) = {84, 62, 40, 18};
// Volume of Cylinder
Physical Volume(97) = {3, 4, 1, 2};
// Generate 3D mesh
Mesh 3;
// Remove duplicate entities
Coherence Mesh;
