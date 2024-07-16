import std.stdio;
import std.string;
import std.conv;
import std.exception : assertThrown;
import el;
import std.typecons;
import std.math;


struct Atom {
  string El;
  double x;
  double y;
  double z;
  double R1; // radii for electrostatic
  double R2; // radii for non-electrostatic cavitation
  double R3; // radii for non-electrostatic dispersion
  double k; // coefficient to calculate the dispersion / repulsion
  double Z; // atomic charge
  double m;
}

struct Solv {
  string Nm;
  double den;
  double p; // den*6.02/M
  double rs; // solvent radius
  double M; // molar mass
  string [] EL;
  ulong [] N;
  double [] k;
  double [] R3;
  double eps;
  double D; // diffusion, to calculate Rsolv
  double n; // viscosity
  double [] param; // parameters to scale Rsolv for PC, SPT-S, SPT-V 
}

struct Sphere {
  string El;
  double x;
  double y;
  double z;
  double S;
  double Z;
  double V;
}

double pi = 3.141592653589793;
double kb = 1.380649E-23;
double Na = 6.02e23;
// elementary charge
double Eel = 1.60217663e-19;
// electric constant
double k = 1/(4*3.141592653589793*8.85e-12);

double distDA (Sphere A, Atom B) {
  double dx2 = (A.x-B.x)*(A.x-B.x);
  double dy2 = (A.y-B.y)*(A.y-B.y);
  double dz2 = (A.z-B.z)*(A.z-B.z);
  double d = pow((dx2 + dy2 + dz2),0.5);
  return d;
}

double distAA (Atom A, Atom B) {
  double dx2 = (A.x-B.x)*(A.x-B.x);
  double dy2 = (A.y-B.y)*(A.y-B.y);
  double dz2 = (A.z-B.z)*(A.z-B.z);
  double d = pow((dx2 + dy2 + dz2),0.5);
  return d;
}


double distDD (Sphere A, Sphere B) {
  double dx2 = (A.x-B.x)*(A.x-B.x);
  double dy2 = (A.y-B.y)*(A.y-B.y);
  double dz2 = (A.z-B.z)*(A.z-B.z);
  double d = pow((dx2 + dy2 + dz2),0.5);
  return d;
}

// function to round the double numbers
double round(double x, uint places)
{
   double pwr = pow(10.0,places);
   return std.math.round(x * pwr) / pwr;
}  
