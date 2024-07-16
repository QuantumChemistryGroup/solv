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


// function to calculate the dispersion

struct DT {
  bool INC; // yes if does not belong to molecule spheres
  double ED; // dispersion
  double ER; // repulsion
  bool INCV; // wheteher the dot is included for mol volume
}

DT GETDOT (Atom com, Atom [] mol, ulong n, Solv solv, double w, double st, ulong num) {
  DT dt;
  dt.INC = false;
  // convert 1d to 3d
  long V = 2*num+1;
  long i = n % V;
  long j = (n / V)%V;
  long k = n / (V*V);
  // get i,j,k signed
  // as it should be -num;+num;
  i=i-num;
  j=j-num;
  k=k-num;
  // writeln(i, " ", j, " ", k, " ");
  // generate dot
  Atom dot;
  dot.x = com.x + st*to!double(i);
  dot.y = com.y + st*to!double(j);
  dot.z = com.z + st*to!double(k);
  // check if overlaps with the spheres of atoms  
  dt.INC = false;
  dt.ED=0;
  dt.ER=0;
  // ietarate over elements in Solv
  for (auto l=0; l<solv.EL.length; l++) {
    // iterate over molecule
    // check if should be included
    bool INC = true;
    for (auto m=0;m<mol.length;m++) {
      double d1 = distAA(dot, mol[m]);
      // get sum of radii solv + solute
      double d2 = mol[m].R2 + solv.R3[l]; // important radii vdW and cav are different.
      if (d1 < d2) {
        INC = false;
        break;
      }
    }
    if (INC == true) {
      dt.INC = true;
      for (auto m=0;m<mol.length;m++) {
        double d1 = distAA(dot, mol[m]);
        // get sum of radii solv + solute
        double d2 = mol[m].R3 + solv.R3[l];
        if (d1 > d2) {
          dt.INC = true;
          // calculate all
          double Rij0 = sqrt(2*mol[m].R3*2*solv.R3[l]);
          double z = d1/Rij0;
          double A = 0.214; // in kcal/mol
          double Eds = -1*mol[m].k*solv.k[l]*(A/(pow(z,6)));
          double C = 47.0e3;
          double a = 12.35;
          double Erp = mol[m].k*solv.k[l]*C*exp(-1*a*z);
          // account for element index
          Eds = Eds*to!double(solv.N[l])*w*(solv.p/1000);
          Erp = Erp*to!double(solv.N[l])*w*(solv.p/1000);
          dt.ED=dt.ED+Eds;
          dt.ER=dt.ER+Erp;
          //writeln(i, " ", j, " ", k, " ", l, " ", m, " ", dt.ED, " ", dt.ER);
        }
      }
    }
  }
  // check wheteher the dot is in molecular cavity
  for (auto m=0;m<mol.length;m++) {
    dt.INCV = false;
    double d1 = distAA(dot, mol[m]);
    // check whether the dot is in atomic radius (for cavitation)
    if (d1 < mol[m].R2) {
      dt.INCV = true;
      break;
    }
  } 
  return dt;
}


double [3] DISPREP (Atom [] mol, Solv solv, ulong num) {
  // first we need to calculate COM
  double M = 0;
  for (auto i=0;i<mol.length;i++) {
    M = M + mol[i].m;
  }
  Atom com; com.x=0;com.y=0;com.z=0;
  for (auto i=0;i<mol.length;i++) {
    double [3] tmp = [0,0,0];
    tmp[0] = (1/M)*mol[i].m*mol[i].x;
    tmp[1] = (1/M)*mol[i].m*mol[i].y;
    tmp[2] = (1/M)*mol[i].m*mol[i].z;
    com.x = com.x + tmp[0];
    com.y = com.y + tmp[1];
    com.z = com.z + tmp[2];
  }
  // then we need to calculate the longest distance to com
  double dist = 0.0;
  for (auto i=0;i<mol.length;i++) {
    double d1 = distAA(mol[i],com);
    if (dist < d1) {
      dist = d1;
    }
  }
  // longest axis value is dist + 15 A
  double X = dist + 23.0;
  // get the weight for each mesh dot
  // should give density upon integration
  double nm = to!double(num);
  double w = pow((2*X),3)/(pow(2*nm+1,3)); 
  double st = X/nm;

  // i,j,k - coordinates of all cubic dots
  // dimension for i = 2num + 1; j = 2num + 1; k = 2num + 1;
  ulong dm = pow((2*num+1),3);
  // https://softwareengineering.stackexchange.com/questions/212808/treating-a-1d-data-structure-as-2d-grid
  // 3d to 1d
  auto N = std.range.iota(0,dm);
  // array for each thread
  double [][] cache;
  cache.length = totalCPUs;
  for (auto i=0;i<totalCPUs;i++) {
    cache[i].length = 3;
    // disp
    cache[i][0] = 0;
    // rep
    cache[i][1] = 0;
    // number of dots in the sphere
    // needed to calculate Volume
    cache[i][2] = 0;
  }
  ulong A = 2*num+1;
  foreach (n; parallel(N)) {
    DT dt = GETDOT (com, mol, n, solv, w, st, num);
    auto m1 = taskPool.workerIndex;
    if (dt.INC == true) {
      cache[m1][0]=cache[m1][0] + dt.ED;
      cache[m1][1]=cache[m1][1] + dt.ER;
    }
    if (dt.INCV == true) {
      cache[m1][2] = cache[m1][2] + 1;
    }      
  }
  // sum over all threads
  double [3] EdrV = [0.0,0.0,0.0]; 
  double Nv = 0; // number of dots in the cavity
  for (auto i=0; i<cache.length;i++) {
    EdrV[0] = EdrV[0] + cache[i][0];
    EdrV[1] = EdrV[1] + cache[i][1];
    Nv = Nv + cache[i][2];
  }
  // calculate the molecular volume;
  double Vcub = pow(2*X,3); // volume of the cube
  double Vd = Vcub/dm; // volume per dot
  double Vm = Nv * Vd;
  EdrV[2] = Vm;
  return EdrV;
}









