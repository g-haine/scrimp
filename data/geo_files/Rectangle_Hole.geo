/*--------------------------
    File Rectangle_Hole.geo
/*-------------------------- */

// Omega
/*
/                     L
/    __________________________________
/   |                                  |
/   |           ___                    |
/   |          /   \                   |
/   |         |  _r_|------------------|  l
/   |          \___/                   |  
/   |                                  |
/   |__________________________________|
/ (0,0)          
/
*/

// Parameters
If (!Exists(h)) h=0.1; EndIf // Rectangle
If (!Exists(hmin)) hmin=0.05; EndIf // Hole

If (!Exists(L)) L=3.; EndIf // Length
If (!Exists(l)) l=1.; EndIf // Section

If (!Exists(r)) r=0.25; EndIf // Hole diameter
If (!Exists(cx)) cx=0.5; EndIf // x-coordinate of its center
If (!Exists(cy)) cy=0.5; EndIf // y-coordinate of its center

// Points

Point(100) = {0,0,0,h};
Point(101) = {L,0,0,h};
Point(102) = {L,l/2.,0,hmin};
Point(103) = {L,l,0,h};
Point(104) = {0,l,0,h};

Point(105) = {cx,cy,0,hmin};
Point(106) = {cx+r,cy,0,hmin/3.};
Point(107) = {cx,cy+r,0,hmin/2.};
Point(108) = {cx-r,cy,0,hmin};
Point(109) = {cx,cy-r,0,hmin/2.};

// Curves

Line(10) = {100,101};
Line(11) = {101,102};
Line(12) = {102,103};
Line(13) = {103,104};
Line(14) = {104,100};

Line Loop(20) = {10,11,12,13,14};

Circle(30) = {106,105,107};
Circle(31) = {107,105,108};
Circle(32) = {108,105,109};
Circle(33) = {109,105,106};

Line Loop(40) = {30,31,32,33};
Plane Surface(1) = {20,40};

Line(15) = {106,102};
Line{15} In Surface{1};

// Physical
Physical Line("Inflow") = {14};
Physical Line("Wall") = {10,13};

Physical Line("Hole") = {30,31,32,33};
Physical Line("Outflow") = {11,12};
Physical Surface("Omega") = {1};

