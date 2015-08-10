{ mkDerivation, base, feast, gfortran, global-lock, hmatrix, hspec
, mtl, openblas, parallel-io-simple, primitive, QuickCheck
, sparse-linear, stdenv, suitesparse, transformers, vector
}:
mkDerivation {
  pname = "feast";
  version = "0.1.0.0";
  src = ./.;
  buildDepends = [
    base global-lock hmatrix mtl parallel-io-simple primitive
    sparse-linear suitesparse transformers vector
  ];
  testDepends = [
    base hspec QuickCheck sparse-linear suitesparse vector
  ];
  extraLibraries = [ feast gfortran openblas ];
  description = "Haskell bindings to the FEAST eigensolver library";
  license = stdenv.lib.licenses.gpl2;
}
