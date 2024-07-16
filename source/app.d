/+dub.sdl:
dependency "mir" version="~>3.2.3"
+/
import std.stdio;
import std.string;
import std.conv;
import std.exception : assertThrown;
import str;
import el;
import std.typecons;
import std.math;
import cav;
import core.stdc.stdlib : exit;
import disp;
import atom;
import std.file: thisExePath;

// this is the D code to calculate the solvation
// the first argument is the file with coordinates and charges

void main (string[] argv)
{
try {
  // get the path of executable
  auto path = thisExePath();
  auto path0 = path.split("/");
  auto lena = path0[path0.length-1].length;
  path = path[0..path.length-lena];
//  writeln(path);
  // solvent name is the second argument
  string solvname = toLower(argv[2]); 
  // grid is the third argument
  ulong gr = to!ulong(argv[3]); 
  // num is another grid to build the mesh
  ulong num = to!ulong(argv[4]);
  // read each line in file til end of file
  // reading the file which is the first argument
  File f = File (argv[1], "r");
  Atom [] mol;
  while (!f.eof()) {
    string line = chomp(f.readln());
    auto d = line.split;
    // read the file by line
    // only string with length of 5 are read
    if (d.length > 3) {
      Atom A;
      try {
        A.El = d[0];
        A.x = to!double(d[1]);
        A.y = to!double(d[2]);
        A.z = to!double(d[3]);
        A.Z = to!double(d[4]);
      }
      catch (Throwable) {
        writeln("Check your mol.xyz.");
        writeln(d);
        writeln("Should be: El  x  y  z  Charge");
        exit(0);
      }
      mol = mol ~ A;
    }
  }
  f.close();

  // read electrostatic radii
  double [string] radel;
  f = File (path~"radii_el.txt", "r");
  while (!f.eof()) {
    string line = chomp(f.readln());
    auto d = line.split;
    if (d.length == 2) {
        radel[toLower(d[0])] = to!double(d[1]);
    }
  }
  // read non-electrostatic cavitation radii
  double [string] radne;
  f = File (path~"radii_cav.txt", "r");
  while (!f.eof()) {
    string line = chomp(f.readln());
    auto d = line.split;
    if (d.length == 2) {
        radne[toLower(d[0])] = to!double(d[1]);
    }
  }

  // read non-electrostatic dispersion radii
  double [string] raddisp;
  f = File (path~"radii_disp.txt", "r");
  while (!f.eof()) {
    string line = chomp(f.readln());
    auto d = line.split;
    if (d.length == 2) {
        raddisp[toLower(d[0])] = to!double(d[1]);
    }
  }


  // read dispersion k
  double [string] kdisp;
  f = File (path~"ki.txt", "r");
  while (!f.eof()) {
    string line = chomp(f.readln());
    auto d = line.split;
    if (d.length == 2) {
        kdisp[toLower(d[0])] = to!double(d[1]);
    }
  }

  // reading the file with the solvent parameters
  f = File (path~"solv.txt", "r");
  // read each line in file til end of file
  Solv solv;
  bool Found = false;
  // check if Rsolv is specified manually
  bool Rsolv = false;
  bool PARAM = false;
  while (!f.eof()) {
    string line = chomp(f.readln());
    auto d = line.split;
    // read the file by line
    // only string with length of 5 are read
    if ( d.length > 1 && toLower(d[0]) == "name") {  
      for (auto i=0;i<d.length;i++) {
        if (toLower(d[i]) == solvname) {
          Found = true;
          solv.Nm = solvname;
        }
      }
    }
    if (Found == true) {
      do {
       line = chomp(f.readln());
       d = line.split;
       if (d.length > 1) {
         if (toLower(d[0]) == "den") {
           // get density
           try {
             solv.den = to!double(d[1]);
           }
           catch (Throwable) {
             writeln("Check solv.txt for ", solvname, " parameters.");
             exit(0);
           }
         }
         else if (toLower(d[0]) == "elements") {
           for (auto i=1;i<d.length;i++) {
             solv.EL ~= toLower(d[i]);
           }
         }
         else if (toLower(d[0]) == "index") {
           for (auto i=1;i<d.length;i++) {
             solv.N ~= to!ulong(d[i]);
           }
         }
         else if (toLower(d[0]) == "param") {
           auto j = 0;
           for (auto i=1;i<d.length;i++) {
             solv.param ~= to!double(d[i]);
             j=j+1;
           }
           if (j != 3) {
                writeln ("Check param for ", solv.Nm);
                writeln ("Will be overwritten with 1 1 1");
             }
           else {
             PARAM = true;
           }  
         }
         // instead I am now calculating Rsolv from 
         // diffusion and viscosity if rsolv is not given
       else if (toLower(d[0]) == "rsolv") {
             solv.rs = to!double(d[1]);
             Rsolv = true;
         }      
         else if (toLower(d[0]) == "eps") {
             try {
               solv.eps = to!double(d[1]);
             }
             catch (Throwable) {
               writeln("Check solv.txt for ", solvname, " parameters.");
               exit(0);
             }
         }      
         else if (toLower(d[0]) == "diff") {
             try {
               solv.D = to!double(d[1]);
             }
             catch (Throwable) {
               writeln("Check solv.txt for ", solvname, " parameters.");
               exit(0);
             }
         }      
         else if (toLower(d[0]) == "vis") {
             try {
               solv.n = to!double(d[1]);
             }
             catch (Throwable) {
               writeln("Check solv.txt for ", solvname, " parameters.");
               exit(0);
             }
         }      
       }
     } while (d.length  > 1 && toLower(d[1]) != "name");
    break;
    }   
  }
  f.close();

  if (Found == false) {
    writeln ("No solvent");
    exit(0);
  }

  // calculate Rsolv
  // see https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.26139
  if (Rsolv != true) {
    solv.rs = 1.27*kb*298.15/(6*pi*solv.n*(1e-3)*(1e-10)*solv.D);
  }
  writeln ("Solvent: ",solv.Nm);
  writeln ("Rsolv = ",solv.rs," (unscaled)");

  // parameters to scale Rsolv
  if (PARAM == false) {
     solv.param.length = 3;
     solv.param = [1,1,1]; 
  }
  writeln ("param: ",solv.param);

  // assign mass to solvent
  double M = 0;
  for (auto i=0;i<solv.EL.length;i++) {
    double m=0;
    try {
       m = mass(toLower(solv.EL[i]))*solv.N[i];
       M = M + m;
    }
    catch (Throwable) {
       writeln("No mass for element ",mol[i].El);
       exit(0);
    }
  }
  solv.M = M;
  writeln("Solvent mass = ", solv.M);
  writeln("Solvent den = ", solv.den);
  // this is the correction for 1atm -> x=1
  auto VM = solv.M/solv.den;
  auto RG = kb*Na;
  auto PN = 101325.0;
  auto TN = 298.15;
  auto Gtr = RG*TN*log(RG*TN/((PN/1000.0)*VM))/4184.0;
  auto V1 = VM/Na; // volume of one molecule
//  auto RV = (V1*3/(4*pi))^(1.0/3);// radius from volume
  auto RV = cbrt((V1*3/(4*pi*1000)))*1e10;// radius from volume
  writeln("Rsolv(vol) = ", RV);
  writeln ("Gtr (1atm->x=1) = ", Gtr);
  // assign p - num den for solvent
  solv.p = solv.den*Na*1e-23/(10*solv.M);

  // assign radii

  for (auto i=0;i<mol.length;i++) {
    try {
       mol[i].R1 = radel[toLower(mol[i].El)];
    }
    catch (Throwable) {
       writeln ("Warning: no R1 for: ", mol[i].El, ", assign: X");
       mol[i].R1 = radel[toLower("X")];
    }
    try {
       mol[i].R2 = radne[toLower(mol[i].El)];
    }
    catch (Throwable) {
       writeln ("Warning: no R2 for: ", mol[i].El, ", assign: X");
       mol[i].R2 = radne[toLower("X")];
    }
    try {
       mol[i].R3 = raddisp[toLower(mol[i].El)];
    }
    catch (Throwable) {
       writeln ("Warning: no R3 for: ", mol[i].El, ", assign: X");
       mol[i].R3 = raddisp[toLower("X")];
    }
    try {
       mol[i].k = kdisp[toLower(mol[i].El)];
    }
    catch (Throwable) {
       writeln ("Warning: no k for: ", mol[i].El, ", assign: X");
       mol[i].k = kdisp[toLower("X")];
    }
  //  writeln(mol[i]);
  }

  // assign mass

  for (auto i=0;i<mol.length;i++) {
    try {
       mol[i].m = mass(toLower(mol[i].El));
    }
    catch (Throwable) {
       writeln("No mass for element ",mol[i].El);
       exit(0);
    }
  }
  // assign k and R to solv elements
  for (auto i=0;i<solv.EL.length;i++) {
    try {
       solv.k ~= kdisp[toLower(solv.EL[i])];
    }
    catch (Throwable) {
       solv.k ~= kdisp[toLower("X")];
    }
    try {
       solv.R3 ~= raddisp[toLower(solv.EL[i])];
    }
    catch (Throwable) {
       solv.R3 ~= raddisp[toLower("X")];
    }
  }
  //writeln(mol);
  //writeln (solv);
  double E = EL (mol, gr, solv);
  writeln("Gel = ", round(E,2));
  double [3] EdrV = [0.0,0.0,0.0]; // disp, rep, vol
  EdrV = DISPREP (mol, solv, num);
  writeln("Gdisp = ", round(EdrV[0],2));
  writeln("Grep = ", round(EdrV[1],2));
  double [3] Ecav = [0.0,0.0,0.0]; // PC, SPT-S, SPT-V
  Ecav = ECAV (mol, gr*100, solv, EdrV[2]); // the last term is the volume
  writeln("Gcav = ", round(Ecav[0],2), " ", round(Ecav[1],2), " ", round(Ecav[2],2));
  double [3] Gsolv = [0.0,0.0,0.0]; // PC, SPT-S, SPT-V
  Gsolv[0] = E + EdrV[0] + EdrV[1] + Ecav[0];
  Gsolv[1] = E + EdrV[0] + EdrV[1] + Ecav[1];
  Gsolv[2] = E + EdrV[0] + EdrV[1] + Ecav[2];
  writeln("Gsolv = ", round(Gsolv[0],2), " ", round(Gsolv[1],2), " ", round(Gsolv[2],2));
  }
  catch (Throwable) {
    writeln ("Too few arguments. ", argv.length-1, ". Exit.");
    writeln ("Arg 1: mol.xyz with charges.");
    writeln ("Arg 2: solvent as in solv.txt");
    writeln ("Arg 3: int (grid) for Electrostatics/Cavitation");
    writeln ("Arg 4: int (grid) for Disp/Rep");
    exit(0);
  }

}

