/+dub.sdl:
dependency "mir" version="~>3.2.3"
dependency "lubeck" version="~>1.5.1"
+/
import std.stdio;
import std.string;
import std.conv;
import std.exception : assertThrown;
import str;
import std.math;
import mir.ndslice.slice: sliced;
import kaleidic.lubeck: mldivide;
// parallel stuff
import core.thread: Thread;
import std.range;
import std.parallelism:parallel;
import std.parallelism:taskPool;
import std.parallelism:totalCPUs;

// parallel loop

// create tmp structure 
// it contains the Sphere dor & bool type
// whether it belongs to the surface (true)
// or it overlaps with the other radii

struct DT {
  bool INC;
  Sphere Sp;
}

DT GETDOT (ulong gr, ulong j, Atom [] mol, ulong i) {
  // generate dot on each spehere for each atom
  double phi = pi*(3.0 - sqrt(5.0));
  DT dt ;
  double y = 1 - (j/(to!double(gr)-1))*2;
  double radius = sqrt(1-y*y);
  double theta = phi*to!double(j);
  double x = cos(theta)*radius;
  double z = sin(theta)*radius;
  Sphere d1;
  d1.El = "X";
  d1.x = x*mol[i].R1;
  d1.y = y*mol[i].R1;
  d1.z = z*mol[i].R1;
  d1.x = d1.x+mol[i].x;
  d1.y = d1.y+mol[i].y;
  d1.z = d1.z+mol[i].z;
  // check overlapping
  dt.INC = true;
  for (auto m=0;m<mol.length;m++) {
    if (m != i) {
      if (distDA(d1,mol[m]) < mol[m].R1) {
        dt.INC = false;
        break;
      }
    }  
  }
  // if belong to the Sphere - calculate potential
  if (dt.INC == true) {
    double Vi=0.0;
    for (auto j1=0; j1<mol.length;j1++) {
      Vi = Vi + mol[j1].Z/distDA(d1,mol[j1]);
    }
    d1.V = Vi;
  }
  dt.Sp = d1;
  return dt;
}

// this is the function to calculate Aij cells
// in parallel

double GETAij (ulong i, ulong j, Sphere dot_i, Sphere dot_j) {
  double Aij;
  if (i == j) {
    Aij = 1.0694*sqrt((4*pi/dot_i.S));
    //Aij = 1.0694*pow((4*pi/dot_i.S),0.5);
  }
  else {
    Aij = 1/(distDD(dot_i,dot_j));
  }
  return Aij;  
}



// this is the function to calculate electrostatic

double EL (Atom [] mol, ulong gr, Solv solv) {
  shared Sphere [] dot;
  double N1=0;
  for (auto i=0;i<mol.length;i++) {
    // Area of atom sphere
    // https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere/44164075#44164075
    double S = 4*pi*mol[i].R1*mol[i].R1;
    double Si = S/to!double(gr);
    auto I = std.range.iota(0,gr);
    Sphere [][] cache;
    cache.length = totalCPUs;
//    writeln ("NCPUs = ",totalCPUs);
    foreach (j; parallel(I)) {
      DT dt = GETDOT (gr, j, mol, i);
      if (dt.INC == true) {
        dt.Sp.S = Si;
        auto m1 = taskPool.workerIndex;
        cache[m1] ~=dt.Sp;
      } 
    }
//    writeln("Parallel finished");
    for (auto j=0; j<cache.length;j++){
      for (auto m=0;m<cache[j].length;m++) {
        dot ~= cache[j][m];
      }
    }
  }
//  writeln("Calculation of V and Aij ", dot.length);
  // eq. 5
  double [] A;
  A.length = dot.length*dot.length;
  double [] B;
  for (auto i=0; i<dot.length;i++){
    B = B ~ (-1*dot[i].V);
  }
  // all combinations
  // i x j
  shared double [][] Aij;
  Aij.length = dot.length;
  for (auto i=0; i<dot.length;i++) {
    Aij[i].length = dot.length;
  }
  auto I = std.range.iota(0,dot.length*dot.length);
  foreach (j;parallel(I)) {
    // convert 1d to 2d index
    // m, n
    // https://softwareengineering.stackexchange.com/questions/212808/treating-a-1d-data-structure-as-2d-grid
    auto m = j % dot.length;
    auto n = j / dot.length;
    Aij[m][n] = GETAij (m, n, dot[m], dot[n]); 
  }
//  writeln("foreach done");
  for (auto i=0; i<dot.length;i++){
    for (auto j=0;j<dot.length;j++) {
      auto m = i + dot.length*j;
      A[m] = Aij[i][j];
      //A = A ~ Aij[i][j];
    }
  }
//  writeln("Population V, Aij done");
  //writeln (A);
  // solving the system of linear equations
  // see here
  // https://github.com/kaleidicassociates/lubeck/blob/master/example/source/app.d
  auto a = A.sliced(dot.length,dot.length);
  auto b = B.sliced;
  auto x = a.mldivide(b);
//  writeln(a);
  //writeln(x);
  double Z = 0.0;
  double qiVi = 0.0;
  for (auto i=0;i<x.length;i++) {
    dot[i].Z = x[i];
    qiVi = qiVi + dot[i].Z*dot[i].V;
    Z=Z+x[i];
}
  double eps = solv.eps;
  double xe;
  // get total charge
  double Zm = 0;
  for (auto i=0;i<mol.length;i++){
    Zm = Zm + mol[i].Z;
  }
  double offset = 0.1;
  // see COSMO wiki
  if (abs(Zm) > abs(offset)) {
    xe = 0.0;
  }
  else {
    xe = 0.5;
  }
  double E = 0.5*((eps-1)/(eps+xe))*qiVi*Eel*Eel*(1/1e-10)*6.02e23*k/(4.184*1000);
  //writeln(Z); // total charge on cavity
  return E;
}
