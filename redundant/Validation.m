(* Validation.m -- Kinematic consistency checks *)
(* Part of the EXTKIN package *)

(* Check if a value is an exact rational or integer (no floats) *)
exactRationalQ[x_Integer] := True;
exactRationalQ[x_Rational] := True;
exactRationalQ[_] := False;

ValidateKinematics[kin_Association] := Module[
  {n, lam, lamT, issues = {}, passed = True,
   momSum, i, j, k, sub, subsets,
   angStored, sqStored, angComputed, sqComputed,
   manStored, manComputed, rowSum, sI},

  n = kin["n"];
  lam = kin["lambda"];
  lamT = kin["lambdaTilde"];

  (* ---- Check exact rationality ---- *)
  If[!And @@ (exactRationalQ /@ Flatten[lam]),
    AppendTo[issues, "FAIL: lambda contains non-rational entries"];
    passed = False;
  ];
  If[!And @@ (exactRationalQ /@ Flatten[lamT]),
    AppendTo[issues, "FAIL: lambdaTilde contains non-rational entries"];
    passed = False;
  ];

  (* ---- Check momentum conservation ---- *)
  momSum = Sum[Outer[Times, lam[[i]], lamT[[i]]], {i, n}];
  If[momSum =!= {{0, 0}, {0, 0}},
    AppendTo[issues, "FAIL: Momentum conservation violated, sum = " <> ToString[momSum]];
    passed = False;
  ,
    AppendTo[issues, "PASS: Momentum conservation"];
  ];

  (* ---- Check on-shell conditions (p_i^2 = 0) ---- *)
  (* For rank-1 matrices p_i = lambda_i (x) lambdaTilde_i, det = 0 always.
     We verify this explicitly. *)
  Do[
    If[Det[Outer[Times, lam[[i]], lamT[[i]]]] =!= 0,
      AppendTo[issues, "FAIL: On-shell condition p_" <> ToString[i] <> "^2 != 0"];
      passed = False;
    ];
  , {i, n}];
  If[passed || !MemberQ[issues, _String?(StringMatchQ[#, "FAIL: On-shell*"] &)],
    AppendTo[issues, "PASS: On-shell conditions (all p_i^2 = 0)"];
  ];

  (* ---- Check bracket consistency ---- *)
  If[KeyExistsQ[kin, "spinorProducts"],
    Do[
      angComputed = AngleBracket[lam, i, j];
      angStored = kin["spinorProducts"]["angle"][{i, j}];
      If[angComputed =!= angStored,
        AppendTo[issues, "FAIL: Angle bracket mismatch <" <> ToString[i] <> "," <> ToString[j] <> ">"];
        passed = False;
      ];
    , {i, n}, {j, i + 1, n}];

    Do[
      sqComputed = SquareBracket[lamT, i, j];
      sqStored = kin["spinorProducts"]["square"][{i, j}];
      If[sqComputed =!= sqStored,
        AppendTo[issues, "FAIL: Square bracket mismatch [" <> ToString[i] <> "," <> ToString[j] <> "]"];
        passed = False;
      ];
    , {i, n}, {j, i + 1, n}];
  ];

  (* ---- Check Mandelstam consistency: s_ij = <ij>[ij] ---- *)
  If[KeyExistsQ[kin, "mandelstams"],
    Do[
      manComputed = AngleBracket[lam, i, j] * SquareBracket[lamT, i, j];
      manStored = kin["mandelstams"][{i, j}];
      If[manComputed =!= manStored,
        AppendTo[issues, "FAIL: Mandelstam s_{" <> ToString[i] <> ToString[j] <> "} mismatch"];
        passed = False;
      ];
    , {i, n}, {j, i + 1, n}];
  ];

  (* ---- Check row sums: Sum_{j != i} s_ij = 0 ---- *)
  Do[
    rowSum = Sum[
      AngleBracket[lam, i, j] * SquareBracket[lamT, i, j],
      {j, DeleteCases[Range[n], i]}
    ];
    If[rowSum =!= 0,
      AppendTo[issues, "FAIL: Row sum for particle " <> ToString[i] <> " = " <> ToString[rowSum] <> " != 0"];
      passed = False;
    ];
  , {i, n}];
  If[passed || !MemberQ[issues, _String?(StringMatchQ[#, "FAIL: Row sum*"] &)],
    AppendTo[issues, "PASS: Row sums vanish (Mandelstam identities)"];
  ];

  (* ---- Check complement identity: s_I = s_{complement(I)} ---- *)
  Do[
    subsets = Subsets[Range[n], {k}];
    Do[
      sI = computeMandelstamSubset[lam, lamT, sub];
      Module[{comp = Complement[Range[n], sub], sComp},
        If[Length[comp] >= 2,
          sComp = computeMandelstamSubset[lam, lamT, comp];
          If[sI =!= sComp,
            AppendTo[issues, "FAIL: s_" <> ToString[sub] <> " != s_" <> ToString[comp]];
            passed = False;
          ];
        ];
      ];
    , {sub, subsets}];
  , {k, 2, Floor[n/2]}];
  If[passed || !MemberQ[issues, _String?(StringMatchQ[#, "FAIL: s_*"] &)],
    AppendTo[issues, "PASS: Complement identity s_I = s_{I^c}"];
  ];

  If[passed,
    True,
    {False, issues}
  ]
];
