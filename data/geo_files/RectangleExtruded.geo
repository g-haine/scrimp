/*--------------------------
    File Rectangle.geo
/*-------------------------- */

// Omega
/*
/           L
/    _______________
/   |               |
/   |               |
/   |               |
/   |               | l
/   |               |  
/   |               |
/   |_______________|
/ (0,0)          
/
*/

// Parameters
If (!Exists(h)) h=0.05; EndIf
If (!Exists(layer)) layer=1./3.; EndIf
If (!Exists(L)) L=3.; EndIf // Length
If (!Exists(l)) l=1.; EndIf // Length

// Points

Point(101) = {0,0,0,h};
Point(102) = {L,0,0,h};
Point(103) = {L,l,0,h};
Point(104) = {0,l,0,h};

// Curves

Line(11) = {101,102};
Line(12) = {102,103};
Line(13) = {103,104};
Line(14) = {104,101};

Line Loop(21) = {11,12,13,14};

Point(105) = {layer*L,layer*l,0,h};
Point(106) = {(1-layer)*L,layer*l,0,h};
Point(107) = {(1-layer)*L,(1-layer)*l,0,h};
Point(108) = {layer*L,(1-layer)*l,0,h};

Line(15) = {105,106};
Line(16) = {106,107};
Line(17) = {107,108};
Line(18) = {108,105};

Line Loop(22) = {-15,-16,-17,-18};

Plane Surface(1) = {21,22};

// Physical
Physical Line("Gamma") = {11,12,13,14,15,16,17,18};
Physical Surface("Omega") = {1};

