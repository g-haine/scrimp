/*--------------------------
    File Rectangle.geo
/*-------------------------- */

// Omega
/*
/           L
/    _______________
/   |               |
/   |   _________   |
/   |  |         |  |
/   |  |         |  | l
/   |  |_________|  |  
/   |               |
/   |_______________|
/ (0,0)          
/
*/

// Parameters
If (!Exists(h)) h=0.05; EndIf
If (!Exists(hmin)) hmin=0.2*h; EndIf
If (!Exists(layer)) layer=0.1; EndIf
If (!Exists(L)) L=3.; EndIf // Length
If (!Exists(l)) l=1.; EndIf // Length

// Points
Point(101) = {0,0,0,hmin};
Point(102) = {L,0,0,hmin/2.};
Point(103) = {L,l,0,hmin/5.};
Point(104) = {0,l,0,hmin/3.};

// Curves
Line(11) = {101,102};
Line(12) = {102,103};
Line(13) = {103,104};
Line(14) = {104,101};

Line Loop(21) = {11,12,13,14};
Plane Surface(1) = {21};

// Inside lines for boundary layer refinement
Point(105) = {layer*L,layer*l,0,h};
Point(106) = {(1-layer)*L,layer*l,0,h/2.};
Point(107) = {(1-layer)*L,(1-layer)*l,0,h/5.};
Point(108) = {layer*L,(1-layer)*l,0,h/3.};

Line(15) = {105,106};
Line(16) = {106,107};
Line(17) = {107,108};
Line(18) = {108,105};
Line{15,16,17,18} In Surface{1};

// Physical
Physical Line("Bottom") = {11};
Physical Line("Right") = {12};
Physical Line("Top") = {13};
Physical Line("Left") = {14};

Physical Surface("Omega") = {1};

