{ pkgs ? (import <nixpkgs> {}) }:

with pkgs;

haskellPackages.callPackage ./. {
  suitesparse = suitesparse_4_4_1;
  globalLock = haskellPackages.callPackage ./global-lock.nix {};
  loops = haskellPackages.callPackage ../loops {};
}
