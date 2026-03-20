(* Spinors.m -- Spinor-helicity variable generation and bracket computations *)
(* Part of the EXTKIN package *)

(* Levi-Civita tensor in 2D spinor space *)
$epm = {{0, 1}, {-1, 0}};

(* ---- Bracket computations ---- *)

AngleBracket[lam_List, i_Integer, j_Integer] :=
  lam[[i, 1]] lam[[j, 2]] - lam[[i, 2]] lam[[j, 1]];

SquareBracket[lamT_List, i_Integer, j_Integer] :=
  lamT[[i, 1]] lamT[[j, 2]] - lamT[[i, 2]] lamT[[j, 1]];

(* ---- Momentum construction ---- *)

(* Single particle: two 2-component spinors -> 2x2 momentum matrix *)
MomentumFromSpinors[lam_?VectorQ, lamT_?VectorQ] :=
  Outer[Times, lam, lamT];

(* All particles: list of angle spinors, list of square spinors -> list of 2x2 matrices *)
MomentumFromSpinors[lamList_?MatrixQ, lamTList_?MatrixQ] :=
  MapThread[Outer[Times, #1, #2] &, {lamList, lamTList}];

(* ---- Non-degeneracy check ---- *)

checkNonDegenerate[lam_, lamT_, n_] := Module[
  {i, j, k, subsets, sI},

  (* Check all angle brackets nonzero *)
  Do[
    If[AngleBracket[lam, i, j] === 0, Return[False, Module]];
  , {i, n}, {j, i + 1, n}];

  (* Check all square brackets nonzero *)
  Do[
    If[SquareBracket[lamT, i, j] === 0, Return[False, Module]];
  , {i, n}, {j, i + 1, n}];

  (* Check multi-particle Mandelstams for subsets of size 2..Floor[n/2] *)
  Do[
    subsets = Subsets[Range[n], {k}];
    Do[
      sI = Sum[
        AngleBracket[lam, sub[[a]], sub[[b]]] *
        SquareBracket[lamT, sub[[a]], sub[[b]]],
        {a, Length[sub]}, {b, a + 1, Length[sub]}
      ];
      If[sI === 0, Return[False, Module]];
    , {sub, subsets}];
  , {k, 2, Floor[n/2]}];

  True
];

(* ---- Random rational spinor generation ---- *)

Options[RandomRationalSpinors] = {
  "NonDegenerate" -> True,
  "Range" -> 9,
  "MaxAttempts" -> 1000
};

RandomRationalSpinors[n_Integer, opts : OptionsPattern[]] := Module[
  {range, maxAttempts, nonDeg, lam, lamT, matrix, attempt},

  range = OptionValue[RandomRationalSpinors, {opts}, "Range"];
  maxAttempts = OptionValue[RandomRationalSpinors, {opts}, "MaxAttempts"];
  nonDeg = OptionValue[RandomRationalSpinors, {opts}, "NonDegenerate"];

  Do[
    (* Generate n-2 free spinor pairs with positive integer entries *)
    lam = Table[{RandomInteger[{1, range}], RandomInteger[{1, range}]}, {n - 2}];
    lamT = Table[{RandomInteger[{1, range}], RandomInteger[{1, range}]}, {n - 2}];

    (* Last two angle spinors form a basis *)
    lam = Join[lam, {{1, 0}, {0, 1}}];

    (* Solve for last two square spinors from momentum conservation:
       Sum_i lambda_i (x) lambdaTilde_i = 0
       With lambda_{n-1} = {1,0}, lambda_n = {0,1}, the solution is:
       lambdaTilde_{n-1} = -M[[1]], lambdaTilde_n = -M[[2]]
       where M = Sum_{i=1}^{n-2} lambda_i (x) lambdaTilde_i *)
    matrix = Sum[Outer[Times, lam[[i]], lamT[[i]]], {i, 1, n - 2}];
    lamT = Join[lamT, {-matrix[[1]], -matrix[[2]]}];

    (* If non-degeneracy not required, return immediately *)
    If[!nonDeg, Return[{lam, lamT}, Module]];

    (* Check non-degeneracy *)
    If[checkNonDegenerate[lam, lamT, n],
      Return[{lam, lamT}, Module]
    ];
  , {attempt, 1, maxAttempts}];

  Message[RandomRationalSpinors::nodegen, maxAttempts];
  $Failed
];
