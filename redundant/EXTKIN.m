(* ::Package:: *)
(* EXTKIN.m -- Exact Rational Kinematics Package for Spinor-Helicity Variables *)
(*
   Generates random rational kinematic points in 4D spinor-helicity formalism
   for massless n-particle scattering processes. Supports generic kinematics,
   kinematics on specified poles, and comprehensive validation.

   All quantities are exact rationals -- no floating point anywhere.

   Usage:
     Get["path/to/EXTKIN/EXTKIN.m"]
     kin = RandomRationalKinematics[6]
     ValidateKinematics[kin]
*)

BeginPackage["EXTKIN`"];

(* ---- Public symbol declarations ---- *)

(* Spinors.m *)
AngleBracket::usage =
  "AngleBracket[lambdaList, i, j] computes the angle bracket <ij> = \
eps_{ab} lambda_i^a lambda_j^b. Returns an exact rational.";

SquareBracket::usage =
  "SquareBracket[lambdaTildeList, i, j] computes the square bracket [ij] = \
eps_{ab} lambdaTilde_i^a lambdaTilde_j^b. Returns an exact rational.";

MomentumFromSpinors::usage =
  "MomentumFromSpinors[lambda, lambdaTilde] returns the 2x2 momentum matrix \
p^{a adot} = lambda^a lambdaTilde^{adot}. Accepts single particle (two \
2-vectors) or all particles (two lists of 2-vectors).";

RandomRationalSpinors::usage =
  "RandomRationalSpinors[n] generates n pairs of rational 2-component \
spinors {lambda_i, lambdaTilde_i} satisfying exact momentum conservation. \
Options: \"NonDegenerate\" (True), \"Range\" (9), \"MaxAttempts\" (1000).";

(* Kinematics.m *)
MandelstamFromSpinors::usage =
  "MandelstamFromSpinors[lambdaList, lambdaTildeList] computes all \
Mandelstam invariants s_I for subsets I with 2 <= |I| <= n-2. Returns \
an Association with sorted subset keys.";

RandomRationalKinematics::usage =
  "RandomRationalKinematics[n] generates a complete rational kinematic \
point: spinors, momenta, Mandelstam invariants, and bracket products. \
Returns an Association with keys \"n\", \"lambda\", \"lambdaTilde\", \
\"momenta\", \"mandelstams\", \"spinorProducts\".";

(* Poles.m *)
RandomRationalKinematicsOnPole::usage =
  "RandomRationalKinematicsOnPole[n, poles] generates rational kinematics \
where specified Mandelstam invariants vanish exactly. poles is a list of \
subsets, e.g. {{1,2}, {3,4,5}}. Options: \"TwoParticlePoleType\" \
(\"Random\"/\"Angle\"/\"Square\" or per-pole rules), \"NonDegenerate\" \
(True, resamples if any non-pole Mandelstam accidentally vanishes), \
\"Range\" (9), \"MaxAttempts\" (1000).";

(* Validation.m *)
ValidateKinematics::usage =
  "ValidateKinematics[kin] checks momentum conservation, on-shell \
conditions, bracket consistency, Mandelstam identities, and exact \
rationality. Returns True if all pass, or {False, issues} with \
diagnostic messages.";

(* ---- Error messages ---- *)
RandomRationalSpinors::nodegen =
  "Failed to generate non-degenerate spinors after `1` attempts.";

RandomRationalKinematicsOnPole::failed =
  "Failed after `1` attempts: `2`";

(* ---- Load implementations ---- *)

Begin["`Private`"];

$EXTKINDir = DirectoryName[$InputFileName];

Get[FileNameJoin[{$EXTKINDir, "Spinors.m"}]];
Get[FileNameJoin[{$EXTKINDir, "Kinematics.m"}]];
Get[FileNameJoin[{$EXTKINDir, "Poles.m"}]];
Get[FileNameJoin[{$EXTKINDir, "Validation.m"}]];

End[]; (* `Private` *)

EndPackage[];
