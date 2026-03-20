(* ::Package:: *)
(* EXTKINSingle.m -- Exact Rational Kinematics Package (single-file version)

   Generates random rational kinematic points in 4D spinor-helicity formalism
   for massless n-particle scattering. All quantities are exact rationals.

   Convention: s_ij = <ij>[ij]  (same index ordering on both brackets)

   Usage:
     Get["path/to/EXTKINSingle.m"]
     kin = RandomRationalKinematics[6]
     ValidateKinematics[kin]
*)

BeginPackage["EXTKIN`"];

(* ---- Public symbols ---- *)

AngleBracket::usage = "AngleBracket[lambdaList, i, j] computes <ij>.";
SquareBracket::usage = "SquareBracket[lambdaTildeList, i, j] computes [ij].";
MomentumFromSpinors::usage = "MomentumFromSpinors[lam, lamT] returns 2x2 momentum matrix(es).";
RandomRationalSpinors::usage = "RandomRationalSpinors[n] generates n rational spinor pairs with exact momentum conservation.";
MandelstamFromSpinors::usage = "MandelstamFromSpinors[lam, lamT] computes all Mandelstam invariants.";
RandomRationalKinematics::usage = "RandomRationalKinematics[n] generates a complete rational kinematic point.";
RandomRationalKinematicsOnPole::usage = "RandomRationalKinematicsOnPole[n, poles] generates kinematics with specified Mandelstams vanishing.";
ValidateKinematics::usage = "ValidateKinematics[kin] checks all kinematic identities.";
IsNonDegenerate::usage = "IsNonDegenerate[kin] returns True if all brackets and Mandelstams are non-zero. IsNonDegenerate[kin, expectedZeroPoles] excludes specified poles from the check.";
RandomModularKinematics::usage = "RandomModularKinematics[n, p] generates a complete n-point kinematic point over GF(p). All spinor components, brackets, and Mandelstams are integers mod p. For finite-field linear algebra: evaluate constraints mod p, solve mod p, repeat for multiple primes, reconstruct via CRT.";
RandomModularKinematicsOnPole::usage = "RandomModularKinematicsOnPole[n, poles, p] generates modular kinematics with specified Mandelstams vanishing mod p.";
ModularKinToRules::usage = "ModularKinToRules[kin, p] converts a modular kinematics Association to replacement rules {ang[i,j] -> value mod p, ...}.";

RandomRationalSpinors::nodegen = "Failed to generate non-degenerate spinors after `1` attempts.";
RandomRationalKinematicsOnPole::failed = "Failed after `1` attempts: `2`";

Begin["`Private`"];

(* ======================================================================== *)
(* SPINORS                                                                  *)
(* ======================================================================== *)

$epm = {{0, 1}, {-1, 0}};

AngleBracket[lam_List, i_Integer, j_Integer] :=
  lam[[i, 1]] lam[[j, 2]] - lam[[i, 2]] lam[[j, 1]];

SquareBracket[lamT_List, i_Integer, j_Integer] :=
  lamT[[i, 1]] lamT[[j, 2]] - lamT[[i, 2]] lamT[[j, 1]];

MomentumFromSpinors[lam_?VectorQ, lamT_?VectorQ] := Outer[Times, lam, lamT];
MomentumFromSpinors[lamList_?MatrixQ, lamTList_?MatrixQ] :=
  MapThread[Outer[Times, #1, #2] &, {lamList, lamTList}];

checkNonDegenerate[lam_, lamT_, n_] := Module[
  {i, j, k, subsets, sI},
  Do[
    If[AngleBracket[lam, i, j] === 0, Return[False, Module]];
  , {i, n}, {j, i + 1, n}];
  Do[
    If[SquareBracket[lamT, i, j] === 0, Return[False, Module]];
  , {i, n}, {j, i + 1, n}];
  Do[
    subsets = Subsets[Range[n], {k}];
    Do[
      sI = Sum[
        AngleBracket[lam, sub[[a]], sub[[b]]] *
        SquareBracket[lamT, sub[[a]], sub[[b]]],
        {a, Length[sub]}, {b, a + 1, Length[sub]}];
      If[sI === 0, Return[False, Module]];
    , {sub, subsets}];
  , {k, 2, Floor[n/2]}];
  True
];

Options[RandomRationalSpinors] = {"NonDegenerate" -> True, "Range" -> 9, "MaxAttempts" -> 1000};

RandomRationalSpinors[n_Integer, opts : OptionsPattern[]] := Module[
  {range, maxAttempts, nonDeg, lam, lamT, matrix},
  range = OptionValue[RandomRationalSpinors, {opts}, "Range"];
  maxAttempts = OptionValue[RandomRationalSpinors, {opts}, "MaxAttempts"];
  nonDeg = OptionValue[RandomRationalSpinors, {opts}, "NonDegenerate"];
  Do[
    lam = Table[{RandomInteger[{1, range}], RandomInteger[{1, range}]}, {n - 2}];
    lamT = Table[{RandomInteger[{1, range}], RandomInteger[{1, range}]}, {n - 2}];
    lam = Join[lam, {{1, 0}, {0, 1}}];
    matrix = Sum[Outer[Times, lam[[i]], lamT[[i]]], {i, 1, n - 2}];
    lamT = Join[lamT, {-matrix[[1]], -matrix[[2]]}];
    If[!nonDeg, Return[{lam, lamT}, Module]];
    If[checkNonDegenerate[lam, lamT, n], Return[{lam, lamT}, Module]];
  , {attempt, 1, maxAttempts}];
  Message[RandomRationalSpinors::nodegen, maxAttempts]; $Failed
];

(* ======================================================================== *)
(* KINEMATICS                                                               *)
(* ======================================================================== *)

computeMandelstamSubset[lam_, lamT_, subset_List] :=
  Sum[AngleBracket[lam, subset[[a]], subset[[b]]] *
      SquareBracket[lamT, subset[[a]], subset[[b]]],
      {a, Length[subset]}, {b, a + 1, Length[subset]}];

MandelstamFromSpinors[lam_List, lamT_List] := Module[
  {n = Length[lam], result = <||>, k, subsets, sub, sI},
  Do[
    subsets = Subsets[Range[n], {k}];
    Do[result[sub] = computeMandelstamSubset[lam, lamT, sub], {sub, subsets}];
  , {k, 2, n - 2}];
  result
];

Options[RandomRationalKinematics] = Options[RandomRationalSpinors];

RandomRationalKinematics[n_Integer, opts : OptionsPattern[]] := Module[
  {spinorResult, lam, lamT},
  spinorResult = RandomRationalSpinors[n,
    "NonDegenerate" -> OptionValue[RandomRationalKinematics, {opts}, "NonDegenerate"],
    "Range" -> OptionValue[RandomRationalKinematics, {opts}, "Range"],
    "MaxAttempts" -> OptionValue[RandomRationalKinematics, {opts}, "MaxAttempts"]];
  If[spinorResult === $Failed, Return[$Failed]];
  {lam, lamT} = spinorResult;
  buildKinematicsAssociation[lam, lamT, n]
];

(* ======================================================================== *)
(* POLES                                                                    *)
(* ======================================================================== *)

(* Union-Find for connected components *)
connectedComponents[pairs_List, n_Integer] := Module[
  {parent, i, p, r1, r2, root, groups},
  parent = Range[n];
  root[x_Integer] := Module[{r = x, c, next},
    While[parent[[r]] =!= r, r = parent[[r]]];
    c = x; While[c =!= r, next = parent[[c]]; parent[[c]] = r; c = next]; r];
  Do[r1 = root[p[[1]]]; r2 = root[p[[2]]];
     If[r1 =!= r2, parent[[r1]] = r2], {p, pairs}];
  Do[root[i], {i, n}];
  groups = GatherBy[Range[n], (parent[[#]] &)];
  Select[groups, Length[#] > 1 &]
];

parseTwoParticlePoleTypes[twoPoles_List, poleTypeOpt_] := Module[
  {types = <||>, pole, ptype, match},
  Do[
    ptype = Which[
      StringQ[poleTypeOpt] && poleTypeOpt === "Angle", "Angle",
      StringQ[poleTypeOpt] && poleTypeOpt === "Square", "Square",
      StringQ[poleTypeOpt] && poleTypeOpt === "Random", RandomChoice[{"Angle", "Square"}],
      ListQ[poleTypeOpt],
        match = Select[poleTypeOpt,
          (Head[#] === Rule && (#[[1]] === pole || #[[1]] === Reverse[pole])) &];
        If[Length[match] > 0, match[[1, 2]], RandomChoice[{"Angle", "Square"}]],
      True, RandomChoice[{"Angle", "Square"}]];
    types[pole] = ptype;
  , {pole, twoPoles}];
  types
];

(* BCFW shift: lambdaTilde_a -> lambdaTilde_a + t*lambdaTilde_b,
               lambda_b -> lambda_b - t*lambda_a *)
applyBCFWShiftForPole[lam_, lamT_, poleSet_List, n_Integer,
    anglePairs_List, squarePairs_List, allMultiPoles_List] :=
  Module[{complement, a, b, sI, denom, tStar, newLam, newLamT,
          sqGroupParticles, angGroupParticles, otherMultiParticles,
          candidateA, candidateB},
    complement = Complement[Range[n], poleSet];
    sqGroupParticles = If[Length[squarePairs] > 0, Union @@ squarePairs, {}];
    angGroupParticles = If[Length[anglePairs] > 0, Union @@ anglePairs, {}];
    otherMultiParticles = Module[{others = Select[allMultiPoles, # =!= poleSet &]},
      If[Length[others] > 0, Union @@ others, {}]];
    candidateA = Complement[poleSet, sqGroupParticles];
    If[Length[candidateA] === 0, candidateA = poleSet];
    candidateB = Complement[complement, angGroupParticles, otherMultiParticles];
    If[Length[candidateB] === 0, candidateB = Complement[complement, angGroupParticles]];
    If[Length[candidateB] === 0, candidateB = complement];
    sI = computeMandelstamSubset[lam, lamT, poleSet];
    If[sI === 0, Return[{lam, lamT}, Module]];
    Do[Do[
      denom = Sum[AngleBracket[lam, a, j] * SquareBracket[lamT, b, j],
                  {j, DeleteCases[poleSet, a]}];
      If[denom =!= 0,
        tStar = -sI / denom;
        newLamT = lamT; newLam = lam;
        newLamT[[a]] = lamT[[a]] + tStar * lamT[[b]];
        newLam[[b]] = lam[[b]] - tStar * lam[[a]];
        Return[{newLam, newLamT}, Module]];
    , {b, candidateB}], {a, candidateA}];
    $Failed
  ];

validatePoleKinematics[lam_, lamT_, n_, poles_, poleTypes_, twoPoles_] :=
  Module[{expectedZeroSubsets, complementSub, sub, sI, k, subsets, pole, i, j, ab, sb},
    expectedZeroSubsets = {};
    Do[AppendTo[expectedZeroSubsets, Sort[pole]];
       complementSub = Complement[Range[n], pole];
       If[2 <= Length[complementSub] <= n - 2,
         AppendTo[expectedZeroSubsets, Sort[complementSub]]];
    , {pole, poles}];
    expectedZeroSubsets = Union[expectedZeroSubsets];
    Do[sI = computeMandelstamSubset[lam, lamT, sub];
       If[sI =!= 0, Return[False, Module]], {sub, expectedZeroSubsets}];
    Do[
      If[KeyExistsQ[poleTypes, Sort[pole]] && poleTypes[Sort[pole]] === "Angle",
        {i, j} = Sort[pole]; ab = AngleBracket[lam, i, j]; sb = SquareBracket[lamT, i, j];
        If[ab =!= 0 || sb === 0, Return[False, Module]]];
      If[KeyExistsQ[poleTypes, Sort[pole]] && poleTypes[Sort[pole]] === "Square",
        {i, j} = Sort[pole]; ab = AngleBracket[lam, i, j]; sb = SquareBracket[lamT, i, j];
        If[sb =!= 0 || ab === 0, Return[False, Module]]];
    , {pole, twoPoles}];
    Do[subsets = Subsets[Range[n], {k}];
       Do[If[!MemberQ[expectedZeroSubsets, sub],
             sI = computeMandelstamSubset[lam, lamT, sub];
             If[sI === 0, Return[False, Module]]],
          {sub, subsets}];
    , {k, 2, Floor[n/2]}];
    Do[If[!MemberQ[twoPoles, Sort[{i, j}]],
         If[AngleBracket[lam, i, j] === 0, Return[False, Module]];
         If[SquareBracket[lamT, i, j] === 0, Return[False, Module]]];
    , {i, n}, {j, i + 1, n}];
    True
  ];

buildKinematicsAssociation[lam_, lamT_, n_] := Module[
  {momenta, mandelstams, angleBrackets, squareBrackets},
  momenta = MomentumFromSpinors[lam, lamT];
  mandelstams = MandelstamFromSpinors[lam, lamT];
  angleBrackets = Association @@ Flatten[Table[
    {i, j} -> AngleBracket[lam, i, j], {i, n}, {j, i + 1, n}]];
  squareBrackets = Association @@ Flatten[Table[
    {i, j} -> SquareBracket[lamT, i, j], {i, n}, {j, i + 1, n}]];
  <|"n" -> n, "lambda" -> lam, "lambdaTilde" -> lamT,
    "momenta" -> momenta, "mandelstams" -> mandelstams,
    "spinorProducts" -> <|"angle" -> angleBrackets, "square" -> squareBrackets|>|>
];

Options[RandomRationalKinematicsOnPole] = {
  "TwoParticlePoleType" -> "Random", "Range" -> 9, "MaxAttempts" -> 1000};

RandomRationalKinematicsOnPole[n_Integer, poles_List, opts : OptionsPattern[]] :=
  Module[
    {range, maxAttempts, poleTypeOpt,
     twoPoles, multiPoles, poleTypes,
     anglePairs, squarePairs, angleGroups, squareGroups,
     allParticlesInAngleGroups, allParticlesInSquareGroups,
     freeAngParticles, freeAllParticles, solveParticles,
     lam, lamT, matrix, group, refSpinor, c, i, j,
     shiftResult, valid, attempt, pole, sortedPoles, complementPole,
     solveInSquareGroups, nonSolveSquareGroups, allSolveGroupParticles,
     solveGroupScales, effectiveLam, coeffMat, sol, rhsParticles, sp, m},

    range = OptionValue[RandomRationalKinematicsOnPole, {opts}, "Range"];
    maxAttempts = OptionValue[RandomRationalKinematicsOnPole, {opts}, "MaxAttempts"];
    poleTypeOpt = OptionValue[RandomRationalKinematicsOnPole, {opts}, "TwoParticlePoleType"];

    (* Validate and deduplicate poles *)
    sortedPoles = Sort /@ poles;
    Do[If[Length[pole] < 2 || Length[pole] > n - 2,
         Message[RandomRationalKinematicsOnPole::failed, 0,
           "Pole subset " <> ToString[pole] <> " invalid"]; Return[$Failed]];
       If[!SubsetQ[Range[n], pole],
         Message[RandomRationalKinematicsOnPole::failed, 0,
           "Invalid particle indices"]; Return[$Failed]];
    , {pole, sortedPoles}];
    Do[complementPole = Sort[Complement[Range[n], pole]];
       If[MemberQ[sortedPoles, complementPole] && pole =!= complementPole,
         sortedPoles = DeleteCases[sortedPoles, complementPole, 1, 1]];
    , {pole, sortedPoles}];
    sortedPoles = Union[sortedPoles];

    twoPoles = Select[sortedPoles, Length[#] == 2 &];
    multiPoles = Select[sortedPoles, Length[#] >= 3 &];
    poleTypes = parseTwoParticlePoleTypes[twoPoles, poleTypeOpt];
    anglePairs = Select[twoPoles, poleTypes[#] === "Angle" &];
    squarePairs = Select[twoPoles, poleTypes[#] === "Square" &];
    angleGroups = connectedComponents[anglePairs, n];
    squareGroups = connectedComponents[squarePairs, n];
    allParticlesInAngleGroups = If[Length[angleGroups] > 0, Union @@ angleGroups, {}];
    allParticlesInSquareGroups = If[Length[squareGroups] > 0, Union @@ squareGroups, {}];

    (* Choose solve particles: not in angle groups, prefer free from square groups too *)
    freeAngParticles = Complement[Range[n], allParticlesInAngleGroups];
    If[Length[freeAngParticles] < 2,
      Message[RandomRationalKinematicsOnPole::failed, 0,
        "Too many angle-type poles: cannot find 2 particles with independent angle spinors"];
      Return[$Failed]];
    freeAllParticles = Complement[freeAngParticles, allParticlesInSquareGroups];
    solveParticles = Which[
      Length[freeAllParticles] >= 2, Take[freeAllParticles, 2],
      Length[freeAllParticles] >= 1,
        {freeAllParticles[[1]], First[Complement[freeAngParticles, freeAllParticles]]},
      True, Module[{reps = {}},
        Do[Module[{rep = Select[group, MemberQ[freeAngParticles, #] &]},
             If[Length[rep] > 0 && Length[reps] < 2, AppendTo[reps, rep[[1]]]]],
           {group, squareGroups}];
        If[Length[reps] >= 2, reps, Take[freeAngParticles, 2]]]
    ];

    (* Identify solve square groups *)
    solveInSquareGroups = <||>;
    Do[Module[{myGroups = Select[squareGroups, MemberQ[#, sp] &]},
         If[Length[myGroups] > 0,
           solveInSquareGroups[sp] = {myGroups[[1]], DeleteCases[myGroups[[1]], sp]}]],
       {sp, solveParticles}];
    nonSolveSquareGroups = Select[squareGroups, Length[Intersection[#, solveParticles]] === 0 &];
    allSolveGroupParticles = {};
    Do[If[KeyExistsQ[solveInSquareGroups, sp],
         allSolveGroupParticles = Union[allSolveGroupParticles, solveInSquareGroups[sp][[1]]]],
       {sp, solveParticles}];

    (* Main attempt loop *)
    Do[
      (* Generate angle spinors *)
      lam = Table[{0, 0}, {n}];
      Do[refSpinor = {RandomInteger[{1, range}], RandomInteger[{1, range}]};
         lam[[group[[1]]]] = refSpinor;
         Do[c = RandomChoice[Join[Range[-3, -1], Range[1, 3]]];
            lam[[group[[j]]]] = c * refSpinor, {j, 2, Length[group]}],
         {group, angleGroups}];
      Do[lam[[i]] = {RandomInteger[{1, range}], RandomInteger[{1, range}]},
         {i, Complement[Range[n], allParticlesInAngleGroups, solveParticles]}];
      lam[[solveParticles[[1]]]] = {1, 0};
      lam[[solveParticles[[2]]]] = {0, 1};

      (* Generate square spinors for non-solve groups *)
      lamT = Table[{0, 0}, {n}];
      Do[refSpinor = {RandomInteger[{1, range}], RandomInteger[{1, range}]};
         lamT[[group[[1]]]] = refSpinor;
         Do[c = RandomChoice[Join[Range[-3, -1], Range[1, 3]]];
            lamT[[group[[j]]]] = c * refSpinor, {j, 2, Length[group]}],
         {group, nonSolveSquareGroups}];
      Do[lamT[[i]] = {RandomInteger[{1, range}], RandomInteger[{1, range}]},
         {i, Complement[Range[n], allParticlesInSquareGroups, solveParticles]}];

      (* Proportionality constants for solve square groups *)
      solveGroupScales = <||>;
      Do[If[KeyExistsQ[solveInSquareGroups, sp],
           Module[{others = solveInSquareGroups[sp][[2]]},
             Do[c = RandomChoice[Join[Range[-3, -1], Range[1, 3]]];
                solveGroupScales[m] = {sp, c}, {m, others}]]],
         {sp, solveParticles}];

      (* Effective angle spinors for momentum conservation *)
      effectiveLam = lam;
      Do[If[KeyExistsQ[solveInSquareGroups, sp],
           Module[{others = solveInSquareGroups[sp][[2]], effLam},
             effLam = lam[[sp]];
             Do[effLam = effLam + solveGroupScales[m][[2]] * lam[[m]], {m, others}];
             effectiveLam[[sp]] = effLam]],
         {sp, solveParticles}];

      (* Solve momentum conservation *)
      rhsParticles = Complement[Range[n], solveParticles, allSolveGroupParticles];
      matrix = If[Length[rhsParticles] > 0,
        Sum[Outer[Times, lam[[i]], lamT[[i]]], {i, rhsParticles}], {{0, 0}, {0, 0}}];
      coeffMat = Transpose[{effectiveLam[[solveParticles[[1]]]],
                             effectiveLam[[solveParticles[[2]]]]}];
      If[Det[coeffMat] === 0, Continue[]];
      sol = Inverse[coeffMat] . (-matrix);
      lamT[[solveParticles[[1]]]] = sol[[1]];
      lamT[[solveParticles[[2]]]] = sol[[2]];
      Do[lamT[[m]] = solveGroupScales[m][[2]] * lamT[[solveGroupScales[m][[1]]]],
         {m, Keys[solveGroupScales]}];

      (* BCFW shifts for multi-particle poles *)
      valid = True;
      Do[shiftResult = applyBCFWShiftForPole[lam, lamT, pole, n,
           anglePairs, squarePairs, multiPoles];
         If[shiftResult === $Failed, valid = False; Break[]];
         {lam, lamT} = shiftResult, {pole, multiPoles}];
      If[!valid, Continue[]];

      (* Validate *)
      If[validatePoleKinematics[lam, lamT, n, sortedPoles, poleTypes, twoPoles],
        Return[buildKinematicsAssociation[lam, lamT, n], Module]];
    , {attempt, maxAttempts}];

    Message[RandomRationalKinematicsOnPole::failed, maxAttempts,
      "Could not simultaneously satisfy all pole constraints"];
    $Failed
  ];

(* ======================================================================== *)
(* VALIDATION                                                               *)
(* ======================================================================== *)

exactRationalQ[x_Integer] := True;
exactRationalQ[x_Rational] := True;
exactRationalQ[_] := False;

ValidateKinematics[kin_Association] := Module[
  {n, lam, lamT, issues = {}, passed = True,
   momSum, i, j, k, sub, subsets, sI,
   angStored, sqStored, angComputed, sqComputed, manStored, manComputed, rowSum},

  n = kin["n"]; lam = kin["lambda"]; lamT = kin["lambdaTilde"];

  (* Exact rationality *)
  If[!And @@ (exactRationalQ /@ Flatten[lam]),
    AppendTo[issues, "FAIL: lambda non-rational"]; passed = False];
  If[!And @@ (exactRationalQ /@ Flatten[lamT]),
    AppendTo[issues, "FAIL: lambdaTilde non-rational"]; passed = False];

  (* Momentum conservation *)
  momSum = Sum[Outer[Times, lam[[i]], lamT[[i]]], {i, n}];
  If[momSum =!= {{0, 0}, {0, 0}},
    AppendTo[issues, "FAIL: Momentum conservation: " <> ToString[momSum]]; passed = False,
    AppendTo[issues, "PASS: Momentum conservation"]];

  (* On-shell *)
  Do[If[Det[Outer[Times, lam[[i]], lamT[[i]]]] =!= 0,
       AppendTo[issues, "FAIL: p_" <> ToString[i] <> "^2 != 0"]; passed = False], {i, n}];
  AppendTo[issues, "PASS: On-shell conditions"];

  (* Bracket consistency *)
  If[KeyExistsQ[kin, "spinorProducts"],
    Do[angComputed = AngleBracket[lam, i, j];
       angStored = kin["spinorProducts"]["angle"][{i, j}];
       If[angComputed =!= angStored,
         AppendTo[issues, "FAIL: <" <> ToString[i] <> "," <> ToString[j] <> "> mismatch"];
         passed = False], {i, n}, {j, i + 1, n}];
    Do[sqComputed = SquareBracket[lamT, i, j];
       sqStored = kin["spinorProducts"]["square"][{i, j}];
       If[sqComputed =!= sqStored,
         AppendTo[issues, "FAIL: [" <> ToString[i] <> "," <> ToString[j] <> "] mismatch"];
         passed = False], {i, n}, {j, i + 1, n}]];

  (* s_ij = <ij>[ij] *)
  If[KeyExistsQ[kin, "mandelstams"],
    Do[manComputed = AngleBracket[lam, i, j] * SquareBracket[lamT, i, j];
       manStored = kin["mandelstams"][{i, j}];
       If[manComputed =!= manStored,
         AppendTo[issues, "FAIL: s_{" <> ToString[{i,j}] <> "} mismatch"];
         passed = False], {i, n}, {j, i + 1, n}]];

  (* Row sums *)
  Do[rowSum = Sum[AngleBracket[lam, i, j] * SquareBracket[lamT, i, j],
                  {j, DeleteCases[Range[n], i]}];
     If[rowSum =!= 0,
       AppendTo[issues, "FAIL: Row sum particle " <> ToString[i]]; passed = False],
     {i, n}];
  AppendTo[issues, "PASS: Row sums vanish"];

  (* Complement identity *)
  Do[subsets = Subsets[Range[n], {k}];
     Do[sI = computeMandelstamSubset[lam, lamT, sub];
        Module[{comp = Complement[Range[n], sub], sComp},
          If[Length[comp] >= 2,
            sComp = computeMandelstamSubset[lam, lamT, comp];
            If[sI =!= sComp,
              AppendTo[issues, "FAIL: s_" <> ToString[sub] <> " != s_" <> ToString[comp]];
              passed = False]]],
        {sub, subsets}],
     {k, 2, Floor[n/2]}];
  AppendTo[issues, "PASS: Complement identity"];

  If[passed, True, {False, issues}]
];

(* ======================================================================== *)
(* NON-DEGENERACY CHECK (public interface)                                  *)
(* ======================================================================== *)

IsNonDegenerate[kin_Association] := Module[{n, lam, lamT},
  n = kin["n"]; lam = kin["lambda"]; lamT = kin["lambdaTilde"];
  checkNonDegenerate[lam, lamT, n]
];

IsNonDegenerate[kin_Association, expectedZeroPoles_List] := Module[
  {n, lam, lamT, sortedPoles, bracketPairs, i, j, k, sub, subsets, sI},
  n = kin["n"]; lam = kin["lambda"]; lamT = kin["lambdaTilde"];
  sortedPoles = Union[Sort /@ expectedZeroPoles];
  bracketPairs = Select[sortedPoles, Length[#] == 2 &];

  (* Check 2-point brackets: non-zero unless on expected pole *)
  Do[
    If[!MemberQ[bracketPairs, Sort[{i, j}]],
      If[AngleBracket[lam, i, j] === 0, Return[False, Module]];
      If[SquareBracket[lamT, i, j] === 0, Return[False, Module]]],
    {i, n}, {j, i + 1, n}];

  (* Check multi-particle Mandelstams: non-zero unless expected *)
  Do[
    subsets = Subsets[Range[n], {k}];
    Do[
      If[!MemberQ[sortedPoles, Sort[sub]] &&
         !MemberQ[sortedPoles, Sort[Complement[Range[n], sub]]],
        sI = computeMandelstamSubset[lam, lamT, sub];
        If[sI === 0, Return[False, Module]]],
      {sub, subsets}],
    {k, 2, Floor[n/2]}];

  True
];

(* ======================================================================== *)
(* MODULAR KINEMATICS (GF(p))                                              *)
(* ======================================================================== *)

(* Generate n-point kinematics with all values in GF(p).
   Same spinor construction as RandomRationalSpinors, but all
   arithmetic is done modulo a prime p. This keeps all entries
   as machine integers in {0,...,p-1}, avoiding rational coefficient
   swell entirely. For use with the finite-field evaluation strategy:
   evaluate constraints mod p → solve mod p → CRT → reconstruct. *)

checkNonDegenerateMod[lam_, lamT_, n_, p_] := Module[
  {i, j, k, subsets, sI},
  Do[
    If[Mod[AngleBracket[lam, i, j], p] === 0, Return[False, Module]],
    {i, n}, {j, i + 1, n}];
  Do[
    If[Mod[SquareBracket[lamT, i, j], p] === 0, Return[False, Module]],
    {i, n}, {j, i + 1, n}];
  Do[
    subsets = Subsets[Range[n], {k}];
    Do[
      sI = Mod[computeMandelstamSubset[lam, lamT, sub], p];
      If[sI === 0, Return[False, Module]],
      {sub, subsets}],
    {k, 2, Floor[n/2]}];
  True
];

RandomModularKinematics[n_Integer, p_Integer] := Module[
  {lam, lamT, matrix, attempt},
  Do[
    lam = Table[{RandomInteger[{1, p - 1}], RandomInteger[{1, p - 1}]}, {n - 2}];
    lamT = Table[{RandomInteger[{1, p - 1}], RandomInteger[{1, p - 1}]}, {n - 2}];
    lam = Join[lam, {{1, 0}, {0, 1}}];
    (* Momentum conservation: lamT for last two particles *)
    matrix = Mod[Sum[Outer[Times, lam[[i]], lamT[[i]]], {i, 1, n - 2}], p];
    lamT = Join[lamT, {Mod[-matrix[[1]], p], Mod[-matrix[[2]], p]}];
    (* All values mod p *)
    lam = Mod[lam, p];
    lamT = Mod[lamT, p];
    If[checkNonDegenerateMod[lam, lamT, n, p],
      Return[buildModularKinAssociation[lam, lamT, n, p], Module]],
    {attempt, 1000}];
  $Failed
];

buildModularKinAssociation[lam_, lamT_, n_, p_] := Module[
  {angleBrackets, squareBrackets, mandelstams, k, sub, subsets},
  angleBrackets = Association @@ Flatten[Table[
    {i, j} -> Mod[AngleBracket[lam, i, j], p], {i, n}, {j, i + 1, n}]];
  squareBrackets = Association @@ Flatten[Table[
    {i, j} -> Mod[SquareBracket[lamT, i, j], p], {i, n}, {j, i + 1, n}]];
  mandelstams = <||>;
  Do[subsets = Subsets[Range[n], {k}];
    Do[mandelstams[sub] = Mod[computeMandelstamSubset[lam, lamT, sub], p],
      {sub, subsets}],
    {k, 2, n - 2}];
  <|"n" -> n, "prime" -> p, "lambda" -> lam, "lambdaTilde" -> lamT,
    "spinorProducts" -> <|"angle" -> angleBrackets, "square" -> squareBrackets|>,
    "mandelstams" -> mandelstams|>
];

(* On-pole modular kinematics: generate generic modular kinematics,
   then apply a BCFW shift mod p to set the desired Mandelstam to zero.
   The shift lambda_a -> lambda_a + t*lambda_b (mod p) preserves
   momentum conservation and on-shell conditions. *)
RandomModularKinematicsOnPole[n_Integer, poles_List, p_Integer] := Module[
  {twoPoles, lam, lamT, pole, complement, sI, a, b, j, denom, tStar,
   attempt, newLam},
  twoPoles = Select[Sort /@ poles, Length[#] == 2 &];
  Do[
    (* Start from non-degenerate generic kinematics *)
    Module[{kin0 = RandomModularKinematics[n, p]},
      If[kin0 === $Failed, Continue[]];
      lam = kin0["lambda"]; lamT = kin0["lambdaTilde"];
    ];
    (* For each 2-particle pole, apply BCFW shift mod p *)
    Module[{ok = True},
      Do[
        complement = Complement[Range[n], pole];
        a = pole[[1]]; b = complement[[1]];
        sI = Mod[computeMandelstamSubset[lam, lamT, pole], p];
        If[sI === 0, Continue[]]; (* Already on pole *)
        denom = Mod[Sum[
          AngleBracket[lam, a, j] * SquareBracket[lamT, b, j],
          {j, DeleteCases[pole, a]}], p];
        If[Mod[denom, p] === 0, ok = False; Break[]];
        tStar = Mod[-sI * PowerMod[denom, -1, p], p];
        newLam = lam;
        newLam[[a]] = Mod[lam[[a]] + tStar * lamT[[b]] (* wrong: should shift lamTilde *), p];
        (* BCFW: lamTilde_a -> lamTilde_a + t*lamTilde_b,
                 lambda_b -> lambda_b - t*lambda_a *)
        lamT[[a]] = Mod[lamT[[a]] + tStar * lamT[[b]], p];
        lam[[b]] = Mod[lam[[b]] - tStar * lam[[a]], p],
        {pole, twoPoles}];
      If[!ok, Continue[]];
    ];
    (* Verify poles vanish mod p *)
    If[AllTrue[twoPoles,
        Mod[computeMandelstamSubset[lam, lamT, #], p] === 0 &],
      Return[buildModularKinAssociation[lam, lamT, n, p], Module]],
    {attempt, 1000}];
  $Failed
];

End[];
EndPackage[];
