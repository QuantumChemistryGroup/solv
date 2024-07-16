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
  d1.x = x*mol[i].R2; // another set of radii is used.
  d1.y = y*mol[i].R2;
  d1.z = z*mol[i].R2;
  d1.x = d1.x+mol[i].x;
  d1.y = d1.y+mol[i].y;
  d1.z = d1.z+mol[i].z;
  // check overlapping
  dt.INC = true;
  for (auto m=0;m<mol.length;m++) {
    if (m != i) {
      if (distDA(d1,mol[m]) < mol[m].R2) {
        dt.INC = false;
        break;
      }
    }  
  } 
  return dt;
}


// this is the function to calculate cavitation energy

double [3] ECAV (Atom [] mol, ulong gr, Solv solv, double Vm) {
  double Sexp=0; // exposed surface for the whole molecule
  double [3] Gcav = [0.0,0.0,0.0]; // cavitation for all molecule
  // PC, SPT-S, SPT-V
  for (auto i=0;i<mol.length;i++) {
    shared Sphere [] dot;
    // Area of atom sphere
    // https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere/44164075#44164075
    double S = 4*pi*mol[i].R2*mol[i].R2;
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
  // here we calculate the cavitation for individual atom
  ulong N = dot.length;
  double A = (to!double(N)/to!double(gr));
  Sexp = Sexp + A*S;
  // constants
  double T = 298.00;
  double P = 101325;
  double Rg = kb*Na;
  double rho = solv.den*1e3*Na/solv.M;
  double sigma1 = 2*solv.rs*solv.param[0]; // scaling factor for PC
  double sigma2 = 2*mol[i].R2;
  double R = sigma2/sigma1;
  double y = (pi*rho*((sigma1*1e-10)*(sigma1*1e-10)*(sigma1*1e-10))/6.0);
  double G1 = -1*log(1-y);
  double G2 = (3*y/(1-y))*R;
  double G31 = (3*y/(1-y)+4.5*(pow(y/(1-y),2)));
  double G32 = G31*pow(R,2);
  double G3 = y*P/(rho*Rg*T)*pow(R,3);
  double G = Rg*T*(G1 + G2 + G32 + G3)/4184;
  double Rv = 1e10*pow(solv.M*3/(solv.den*Na*1e3*4*pi),1/3.);
//  writeln("Vol: ", Rv); // radii from molecular volume
//  writeln("Einstein - Stocks ", solv.rs); // radii from JCC
  // for water and for octanol, Rv / Rs = 1.40
  // may be use it later if we have no diff and vis.
  double Gcavi = A*G;
//  writeln("Atom: ", mol[i].El, i+1, " S0: ", S, " A: ", A, " A*S: ", A*S, " N*Si: ", N*Si, " N: ", N, " Si: ", Si, " G: ", G, " A*G: ", A*G);
  Gcav[0] = Gcav[0] + Gcavi;
  }
  // constants as above, the only difference is radius of cavity
  double T = 298.00;
  double P = 101325;
  double Rg = kb*Na;
  double rho = solv.den*1e3*Na/solv.M;
  // related to cavity radius R
  double [2] RAD = [0.0,0.0]; // SPT-S, SPT-V
  RAD[0] = pow((Sexp/(4*pi)),0.5); // SPT-S
  RAD[1] = pow((3*Vm/(4*pi)),(1/3.0));
  writeln ("Rmol (SPT-S) = ", RAD[0]); // SPT-S radius
  writeln ("Rmol (SPT-V) = ", RAD[1]); // SPT-V radius
//  writeln(RAD, " ", Sexp);
  for (auto i=0; i<RAD.length; i++) {
    double sigma1 = 0;
    if (i == 0) {
      sigma1 = 2*solv.rs*solv.param[1]; // SPT-S
    }
    else {
      sigma1 = 2*solv.rs*solv.param[2]; // SPT-V
    }
    double y = (pi*rho*((sigma1*1e-10)*(sigma1*1e-10)*(sigma1*1e-10))/6.0);
    double G1 = -1*log(1-y);
    double G31 = (3*y/(1-y)+4.5*(pow(y/(1-y),2)));
    double R0 = RAD[i];
    double sigma2 = 2*R0;
    double R = sigma2/sigma1;
    double G2 = (3*y/(1-y))*R;
    double G32 = G31*pow(R,2);
    double G3 = y*P/(rho*Rg*T)*pow(R,3);
    double G = Rg*T*(G1 + G2 + G32 + G3)/4184;
    Gcav[i+1] = G;
  }
  
  return Gcav;
}
