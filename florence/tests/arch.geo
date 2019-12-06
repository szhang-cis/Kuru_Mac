// Gmsh project created on Tue Jan 22 18:33:14 2019
Point(1) = {54.295, 0, 0, 1};
Point(2) = {55.705, 0, 0, 1};
Line(1) = {1, 2};
// Base circle for the arch
Extrude {{0, 0, 1}, {45.0, 0, 0}, Pi/2} {
  Line{1};Layers{10};Recombine;
}
Extrude {{0, 0, 1}, {45.0, 0, 0}, Pi/2} {
  Line{2};Layers{10};Recombine;
}
Extrude {{0, 0, 1}, {45.0, 0, 0}, Pi/2} {
  Line{6};Layers{10};Recombine;
}
Extrude {{0, 0, 1}, {45.0, 0, 0}, Pi/2} {
  Line{10};Layers{10};Recombine;
}
// Arch
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Surface{9, 5, 17, 13};Layers{15};Recombine;
}
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Surface{61, 83, 39, 105};Layers{15};Recombine;
}
// X left circular surface
Physical Surface(194) = {127, 149, 193, 171};
// X right circular surface
Physical Surface(195) = {13, 9, 5, 17};
// Inner surface
Physical Surface(196) = {104, 82, 60, 38, 126, 148, 192, 170};
// Outer surface
Physical Surface(197) = {118, 140, 162, 184, 52, 74, 96, 30};
// Volume of Toroid
Physical Volume(198) = {4, 1, 2, 3, 8, 6, 7, 5};
// Generate 3D Mesh
Mesh 3;
// Remove duplicate entities
Coherence Mesh;
