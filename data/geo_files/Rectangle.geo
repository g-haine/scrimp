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
If (!Exists(h)) h=0.08; EndIf
If (!Exists(layer)) layer=0.35; EndIf
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

Plane Surface(1) = {21};


// Physical
Physical Line("Gamma") = {11,12,13,14};
Physical Surface("Omega") = {1};

