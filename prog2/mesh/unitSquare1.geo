cl1 = 1;

Point(1) = {0, 0, 0, cl1};
Point(2) = {1, 0, 0, cl1};
Point(3) = {1, 1, 0, cl1};
Point(4) = {0, 1, 0, cl1};
Line(1) = {2, 3};
Line(2) = {4, 3};
Line(3) = {1, 4};
Line(4) = {2, 1};
Line Loop(6) = {1, -2, -3, -4};
Plane Surface(6) = {6};
Physical Line(2) = {2, 3, 4};
Physical Surface(1) = {6};
Physical Line(3) = {1};
