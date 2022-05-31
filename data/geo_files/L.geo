/*--------------
    File L.geo
/*-------------- */

// Parameters
If (!Exists(h)) h=0.05; EndIf
If (!Exists(L)) L=3.; EndIf // Length
If (!Exists(l)) l=1.; EndIf // Length

// Points

Point(101) = {0,0,0,h};
Point(102) = {L,0,0,h};
Point(103) = {L,l/2.,0,h};
Point(104) = {L/2.,l/2.,0,h};
Point(105) = {L/2.,l,0,h};
Point(106) = {0,l,0,h};

// Curves

Line(11) = {101,102};
Line(12) = {102,103};
Line(13) = {103,104};
Line(14) = {104,105};
Line(15) = {105,106};
Line(16) = {106,101};

Line Loop(21) = {11,12,13,14,15,16};

Plane Surface(1) = {21};

// Physical
Physical Line("Gamma") = {11,12,13,14,15,16};
Physical Surface("Omega") = {1};
