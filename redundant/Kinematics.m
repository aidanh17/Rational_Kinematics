(* Kinematics.m -- Mandelstam invariants and kinematic point generation *)
(* Part of the EXTKIN package *)

(* ---- Mandelstam computation from spinors ---- *)

(* Compute the Mandelstam invariant for a single subset I:
   s_I = Sum_{i<j in I} <ij>[ij] *)
computeMandelstamSubset[lam_, lamT_, subset_List] :=
  Sum[
    AngleBracket[lam, subset[[a]], subset[[b]]] *
    SquareBracket[lamT, subset[[a]], subset[[b]]],
    {a, Length[subset]}, {b, a + 1, Length[subset]}
  ];

(* Compute all independent Mandelstam invariants.
   Returns an Association with sorted subset keys -> rational values.
   Computes s_I for all subsets I with 2 <= |I| <= n-2.
   Note: s_I = s_{complement(I)} by momentum conservation. *)
MandelstamFromSpinors[lam_List, lamT_List] := Module[
  {n = Length[lam], result = <||>, k, subsets, sub, sI},

  Do[
    subsets = Subsets[Range[n], {k}];
    Do[
      sI = computeMandelstamSubset[lam, lamT, sub];
      result[sub] = sI;
    , {sub, subsets}];
  , {k, 2, n - 2}];

  result
];

(* ---- Full kinematic point generation ---- *)

Options[RandomRationalKinematics] = Options[RandomRationalSpinors];

RandomRationalKinematics[n_Integer, opts : OptionsPattern[]] := Module[
  {spinorResult, lam, lamT, momenta, mandelstams, angleBrackets, squareBrackets},

  spinorResult = RandomRationalSpinors[n,
    "NonDegenerate" -> OptionValue[RandomRationalKinematics, {opts}, "NonDegenerate"],
    "Range" -> OptionValue[RandomRationalKinematics, {opts}, "Range"],
    "MaxAttempts" -> OptionValue[RandomRationalKinematics, {opts}, "MaxAttempts"]
  ];
  If[spinorResult === $Failed, Return[$Failed]];

  {lam, lamT} = spinorResult;

  momenta = MomentumFromSpinors[lam, lamT];
  mandelstams = MandelstamFromSpinors[lam, lamT];

  angleBrackets = Association @@ Flatten[Table[
    {i, j} -> AngleBracket[lam, i, j],
    {i, n}, {j, i + 1, n}
  ]];

  squareBrackets = Association @@ Flatten[Table[
    {i, j} -> SquareBracket[lamT, i, j],
    {i, n}, {j, i + 1, n}
  ]];

  <|
    "n" -> n,
    "lambda" -> lam,
    "lambdaTilde" -> lamT,
    "momenta" -> momenta,
    "mandelstams" -> mandelstams,
    "spinorProducts" -> <|"angle" -> angleBrackets, "square" -> squareBrackets|>
  |>
];
