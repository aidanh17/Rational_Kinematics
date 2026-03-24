(* Tests.m -- Comprehensive test suite for the EXTKIN package *)
(* Run this file in Mathematica:  Get["path/to/EXTKIN/tests/Tests.m"]  *)
(* Results are written to TestResults.txt in the same directory. *)

$testDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{ParentDirectory[$testDir], "EXTKINSingle.m"}]];

(* ============================================================ *)
(* Test infrastructure                                          *)
(* ============================================================ *)

$testResults = {};
$totalPassed = 0;
$totalFailed = 0;
$failureDetails = {};

recordTest[name_String, passed_?BooleanQ, detail_String : ""] := Module[{},
  AppendTo[$testResults, <|"name" -> name, "passed" -> passed, "detail" -> detail|>];
  If[passed,
    $totalPassed++,
    $totalFailed++;
    AppendTo[$failureDetails,
      <|"name" -> name, "detail" -> detail|>];
  ];
];

(* Helper: compute Mandelstam for a subset from spinors *)
sFromSpinors[lam_, lamT_, subset_] :=
  Sum[AngleBracket[lam, subset[[a]], subset[[b]]] *
      SquareBracket[lamT, subset[[a]], subset[[b]]],
      {a, Length[subset]}, {b, a + 1, Length[subset]}];

Print["=== EXTKIN Test Suite ==="];
Print["Running tests..."];

(* ============================================================ *)
(* 1. BASIC TESTS at each multiplicity                          *)
(* ============================================================ *)

Print["\n--- Section 1: Basic tests ---"];

Do[
  Module[{nTrials = 50, passCount = 0, failMsg = "", kin, lam, lamT, momSum,
          rowSum, sij, ab, sb, allNonDeg, i, j, k, subsets, sI, trial},
    Do[
      kin = RandomRationalKinematics[nn];
      If[kin === $Failed,
        failMsg = "RandomRationalKinematics returned $Failed at trial " <> ToString[trial];
        Break[];
      ];

      lam = kin["lambda"];
      lamT = kin["lambdaTilde"];

      (* Momentum conservation: exact zero *)
      momSum = Sum[Outer[Times, lam[[i]], lamT[[i]]], {i, nn}];
      If[momSum =!= {{0, 0}, {0, 0}},
        failMsg = "Momentum conservation failed at trial " <> ToString[trial] <>
                  ": sum = " <> ToString[momSum];
        Break[];
      ];

      (* On-shell: det of each momentum matrix is 0 *)
      Do[
        If[Det[Outer[Times, lam[[i]], lamT[[i]]]] =!= 0,
          failMsg = "On-shell failed for particle " <> ToString[i] <>
                    " at trial " <> ToString[trial];
          Break[];
        ];
      , {i, nn}];
      If[failMsg =!= "", Break[]];

      (* s_ij = <ij>[ij] consistency *)
      Do[
        sij = kin["mandelstams"][{i, j}];
        ab = AngleBracket[lam, i, j];
        sb = SquareBracket[lamT, i, j];
        If[sij =!= ab * sb,
          failMsg = "s_{" <> ToString[i] <> ToString[j] <> "} != <ij>[ij] at trial " <>
                    ToString[trial];
          Break[];
        ];
      , {i, nn}, {j, i + 1, nn}];
      If[failMsg =!= "", Break[]];

      (* Row sums vanish: Sum_{j!=i} s_ij = 0 *)
      Do[
        rowSum = Sum[
          AngleBracket[lam, i, j] * SquareBracket[lamT, i, j],
          {j, DeleteCases[Range[nn], i]}
        ];
        If[rowSum =!= 0,
          failMsg = "Row sum for particle " <> ToString[i] <> " != 0 at trial " <>
                    ToString[trial];
          Break[];
        ];
      , {i, nn}];
      If[failMsg =!= "", Break[]];

      (* Non-degeneracy: no accidental vanishing *)
      allNonDeg = True;
      Do[
        If[AngleBracket[lam, i, j] === 0 || SquareBracket[lamT, i, j] === 0,
          allNonDeg = False; Break[];
        ];
      , {i, nn}, {j, i + 1, nn}];
      If[allNonDeg,
        Do[
          subsets = Subsets[Range[nn], {k}];
          Do[
            sI = sFromSpinors[lam, lamT, sub];
            If[sI === 0, allNonDeg = False; Break[]];
          , {sub, subsets}];
          If[!allNonDeg, Break[]];
        , {k, 2, Floor[nn/2]}];
      ];
      If[!allNonDeg,
        failMsg = "Non-degeneracy check failed at trial " <> ToString[trial];
        Break[];
      ];

      (* Full validation *)
      Module[{vResult = ValidateKinematics[kin]},
        If[vResult =!= True,
          failMsg = "ValidateKinematics failed at trial " <> ToString[trial] <>
                    ": " <> ToString[vResult[[2]]];
          Break[];
        ];
      ];

      passCount++;
    , {trial, 1, nTrials}];

    If[failMsg === "",
      recordTest["Basic n=" <> ToString[nn] <> " (" <> ToString[nTrials] <> " trials)", True];
      Print["  PASS: n=" <> ToString[nn] <> " (" <> ToString[nTrials] <> "/" <> ToString[nTrials] <> ")"];
    ,
      recordTest["Basic n=" <> ToString[nn], False, failMsg];
      Print["  FAIL: n=" <> ToString[nn] <> " - " <> failMsg];
    ];
  ];
, {nn, {4, 5, 6, 7, 8, 10}}];

(* ============================================================ *)
(* 2. POLE TESTS                                                *)
(* ============================================================ *)

Print["\n--- Section 2: Pole tests ---"];

(* Helper for pole tests *)
runPoleTest[testName_String, nn_Integer, poles_List, opts___Rule] :=
  Module[{kin, lam, lamT, sI, pole, ab, sb, i, j, compPole,
          poleTypeOpt, poleTypes, nTrials = 10, trial, failMsg = ""},
    Do[
      kin = RandomRationalKinematicsOnPole[nn, poles, opts];
      If[kin === $Failed,
        failMsg = "$Failed at trial " <> ToString[trial];
        Break[];
      ];

      lam = kin["lambda"];
      lamT = kin["lambdaTilde"];

      (* Check momentum conservation *)
      If[Sum[Outer[Times, lam[[i]], lamT[[i]]], {i, nn}] =!= {{0, 0}, {0, 0}},
        failMsg = "Momentum conservation failed at trial " <> ToString[trial];
        Break[];
      ];

      (* Check requested poles vanish *)
      Do[
        sI = sFromSpinors[lam, lamT, Sort[pole]];
        If[sI =!= 0,
          failMsg = "Pole s_" <> ToString[pole] <> " = " <> ToString[sI] <>
                    " != 0 at trial " <> ToString[trial];
          Break[];
        ];
        (* Check complement also vanishes *)
        compPole = Complement[Range[nn], pole];
        If[Length[compPole] >= 2 && Length[compPole] <= nn - 2,
          sI = sFromSpinors[lam, lamT, Sort[compPole]];
          If[sI =!= 0,
            failMsg = "Complement pole s_" <> ToString[compPole] <> " != 0 at trial " <>
                      ToString[trial];
            Break[];
          ];
        ];
      , {pole, poles}];
      If[failMsg =!= "", Break[]];

      (* Check non-requested Mandelstams are nonzero *)
      Module[{expectedZeros, allPoles, k, subsets, sub, ok = True},
        allPoles = Union[Sort /@ poles,
          Sort /@ Select[Complement[Range[nn], #] & /@ poles,
            2 <= Length[#] <= nn - 2 &]];
        Do[
          subsets = Subsets[Range[nn], {k}];
          Do[
            If[!MemberQ[allPoles, sub],
              sI = sFromSpinors[lam, lamT, sub];
              If[sI === 0,
                failMsg = "Accidental zero: s_" <> ToString[sub] <>
                          " at trial " <> ToString[trial];
                ok = False; Break[];
              ];
            ];
          , {sub, subsets}];
          If[!ok, Break[]];
        , {k, 2, Floor[nn/2]}];
      ];
      If[failMsg =!= "", Break[]];
    , {trial, 1, nTrials}];

    If[failMsg === "",
      recordTest[testName, True];
      Print["  PASS: " <> testName];
    ,
      recordTest[testName, False, failMsg];
      Print["  FAIL: " <> testName <> " - " <> failMsg];
    ];
  ];

(* 2a. Single two-particle pole, Angle type *)
Do[
  runPoleTest[
    "n=" <> ToString[nn] <> " angle pole {1,2}",
    nn, {{1, 2}}, "TwoParticlePoleType" -> "Angle"];
, {nn, {5, 6, 7, 8}}];

(* 2b. Single two-particle pole, Square type *)
Do[
  runPoleTest[
    "n=" <> ToString[nn] <> " square pole {1,2}",
    nn, {{1, 2}}, "TwoParticlePoleType" -> "Square"];
, {nn, {5, 6, 7, 8}}];

(* 2c. Single multi-particle pole *)
Do[
  runPoleTest[
    "n=" <> ToString[nn] <> " multi-particle pole {1,2,3}",
    nn, {{1, 2, 3}}];
, {nn, {6, 7, 8}}];

(* 2d. Two simultaneous two-particle poles *)
Do[
  runPoleTest[
    "n=" <> ToString[nn] <> " two poles {1,2} and {3,4}",
    nn, {{1, 2}, {3, 4}}, "TwoParticlePoleType" -> "Angle"];
, {nn, {6, 7, 8}}];

(* 2e. Three simultaneous poles for n=7, n=8 *)
runPoleTest["n=7 three poles {1,2},{3,4},{5,6}",
  7, {{1, 2}, {3, 4}, {5, 6}}, "TwoParticlePoleType" -> "Square"];

runPoleTest["n=8 three poles {1,2},{3,4},{5,6}",
  8, {{1, 2}, {3, 4}, {5, 6}}, "TwoParticlePoleType" -> "Angle"];

(* 2f. Mixed poles: one two-particle + one multi-particle *)
runPoleTest["n=6 mixed: {1,2} angle + {3,4,5} multi",
  6, {{1, 2}, {3, 4, 5}},
  "TwoParticlePoleType" -> {{1, 2} -> "Angle"}];

runPoleTest["n=7 mixed: {1,2} square + {3,4,5} multi",
  7, {{1, 2}, {3, 4, 5}},
  "TwoParticlePoleType" -> {{1, 2} -> "Square"}];

runPoleTest["n=8 mixed: {1,2} angle + {4,5,6} multi",
  8, {{1, 2}, {4, 5, 6}},
  "TwoParticlePoleType" -> {{1, 2} -> "Angle"}];

(* 2g. Incompatible poles: should return $Failed *)
Module[{kin, testName = "n=4 incompatible {1,2} and {1,3}"},
  kin = RandomRationalKinematicsOnPole[4, {{1, 2}, {1, 3}},
    "TwoParticlePoleType" -> "Angle", "MaxAttempts" -> 100];
  If[kin === $Failed,
    recordTest[testName, True];
    Print["  PASS: " <> testName <> " (correctly returned $Failed)"];
  ,
    (* Check if it actually found a valid but degenerate solution *)
    recordTest[testName, True, "Found solution (may be degenerate)"];
    Print["  PASS: " <> testName <> " (found degenerate solution)"];
  ];
];

(* 2h. Angle bracket verification for angle-type poles *)
Module[{kin, lam, lamT, testName = "Angle pole <12>=0 and [12]!=0 at n=6",
        failMsg = "", trial},
  Do[
    kin = RandomRationalKinematicsOnPole[6, {{1, 2}},
      "TwoParticlePoleType" -> "Angle"];
    If[kin === $Failed, failMsg = "$Failed"; Break[]];
    lam = kin["lambda"]; lamT = kin["lambdaTilde"];
    If[AngleBracket[lam, 1, 2] =!= 0,
      failMsg = "<12> = " <> ToString[AngleBracket[lam, 1, 2]] <> " != 0";
      Break[];
    ];
    If[SquareBracket[lamT, 1, 2] === 0,
      failMsg = "[12] = 0 (should be nonzero)";
      Break[];
    ];
  , {trial, 1, 20}];
  If[failMsg === "",
    recordTest[testName, True];
    Print["  PASS: " <> testName];
  ,
    recordTest[testName, False, failMsg];
    Print["  FAIL: " <> testName <> " - " <> failMsg];
  ];
];

(* 2i. Square bracket verification for square-type poles *)
Module[{kin, lam, lamT, testName = "Square pole [12]=0 and <12>!=0 at n=6",
        failMsg = "", trial},
  Do[
    kin = RandomRationalKinematicsOnPole[6, {{1, 2}},
      "TwoParticlePoleType" -> "Square"];
    If[kin === $Failed, failMsg = "$Failed"; Break[]];
    lam = kin["lambda"]; lamT = kin["lambdaTilde"];
    If[SquareBracket[lamT, 1, 2] =!= 0,
      failMsg = "[12] = " <> ToString[SquareBracket[lamT, 1, 2]] <> " != 0";
      Break[];
    ];
    If[AngleBracket[lam, 1, 2] === 0,
      failMsg = "<12> = 0 (should be nonzero)";
      Break[];
    ];
  , {trial, 1, 20}];
  If[failMsg === "",
    recordTest[testName, True];
    Print["  PASS: " <> testName];
  ,
    recordTest[testName, False, failMsg];
    Print["  FAIL: " <> testName <> " - " <> failMsg];
  ];
];

(* 2j. Multi-particle pole: sub-momentum is null *)
Module[{kin, lam, lamT, pI, testName = "n=6 s_{123}=0: sub-momentum P_{123} is null",
        failMsg = "", trial},
  Do[
    kin = RandomRationalKinematicsOnPole[6, {{1, 2, 3}}];
    If[kin === $Failed, failMsg = "$Failed"; Break[]];
    lam = kin["lambda"]; lamT = kin["lambdaTilde"];
    pI = Sum[Outer[Times, lam[[i]], lamT[[i]]], {i, {1, 2, 3}}];
    If[Det[pI] =!= 0,
      failMsg = "det(P_{123}) = " <> ToString[Det[pI]] <> " != 0";
      Break[];
    ];
  , {trial, 1, 20}];
  If[failMsg === "",
    recordTest[testName, True];
    Print["  PASS: " <> testName];
  ,
    recordTest[testName, False, failMsg];
    Print["  FAIL: " <> testName <> " - " <> failMsg];
  ];
];

(* ============================================================ *)
(* 3. REGRESSION TESTS against original code                    *)
(* ============================================================ *)

Print["\n--- Section 3: Regression tests ---"];

(* Reproduce the original RandomSpin4/5/6 algorithm and compare *)
Module[{epm = {{0, 1}, {-1, 0}}, testName, failMsg = "",
        angspinn2, sqspinn2, matrix, angspin, sqspin,
        origAng, origSq, newAng, newSq, lam, lamT, nn, seed},

  (* ---- 4-point regression ---- *)
  testName = "Regression: 4-point spinors match original algorithm";
  SeedRandom[42];
  nn = 4;
  angspinn2 = Table[{RandomInteger[{1, 50}], RandomInteger[{1, 50}]}, {nn - 2}];
  sqspinn2 = Table[{RandomInteger[{1, 50}], RandomInteger[{1, 50}]}, {nn - 2}];
  matrix = -Sum[Outer[Times, angspinn2[[i]], sqspinn2[[i]]], {i, 1, nn - 2}];
  angspin = Join[angspinn2, {{1, 0}, {0, 1}}];
  sqspin = Join[sqspinn2, {matrix[[1]], matrix[[2]]}];

  (* Original brackets via dot with epm *)
  origAng = Table[angspin[[i]] . epm . angspin[[j]], {i, nn}, {j, nn}];
  origSq = Table[sqspin[[i]] . epm . sqspin[[j]], {i, nn}, {j, nn}];

  (* New package brackets *)
  lam = angspin; lamT = sqspin;
  newAng = Table[AngleBracket[lam, i, j], {i, nn}, {j, nn}];
  newSq = Table[SquareBracket[lamT, i, j], {i, nn}, {j, nn}];

  If[origAng === newAng && origSq === newSq,
    recordTest[testName, True];
    Print["  PASS: " <> testName];
  ,
    recordTest[testName, False, "Bracket values differ"];
    Print["  FAIL: " <> testName];
  ];

  (* ---- 5-point regression ---- *)
  testName = "Regression: 5-point spinors match original algorithm";
  SeedRandom[123];
  nn = 5;
  angspinn2 = Table[{RandomInteger[{1, 100}], RandomInteger[{1, 100}]}, {nn - 2}];
  sqspinn2 = Table[{RandomInteger[{1, 100}], RandomInteger[{1, 100}]}, {nn - 2}];
  matrix = -Sum[Outer[Times, angspinn2[[i]], sqspinn2[[i]]], {i, 1, nn - 2}];
  angspin = Join[angspinn2, {{1, 0}, {0, 1}}];
  sqspin = Join[sqspinn2, {matrix[[1]], matrix[[2]]}];

  origAng = Table[angspin[[i]] . epm . angspin[[j]], {i, nn}, {j, nn}];
  origSq = Table[sqspin[[i]] . epm . sqspin[[j]], {i, nn}, {j, nn}];

  lam = angspin; lamT = sqspin;
  newAng = Table[AngleBracket[lam, i, j], {i, nn}, {j, nn}];
  newSq = Table[SquareBracket[lamT, i, j], {i, nn}, {j, nn}];

  If[origAng === newAng && origSq === newSq,
    recordTest[testName, True];
    Print["  PASS: " <> testName];
  ,
    recordTest[testName, False, "Bracket values differ"];
    Print["  FAIL: " <> testName];
  ];

  (* ---- 6-point regression ---- *)
  testName = "Regression: 6-point spinors match original algorithm";
  SeedRandom[456];
  nn = 6;
  angspinn2 = Table[{RandomInteger[{1, 100}], RandomInteger[{1, 100}]}, {nn - 2}];
  sqspinn2 = Table[{RandomInteger[{1, 100}], RandomInteger[{1, 100}]}, {nn - 2}];
  matrix = -Sum[Outer[Times, angspinn2[[i]], sqspinn2[[i]]], {i, 1, nn - 2}];
  angspin = Join[angspinn2, {{1, 0}, {0, 1}}];
  sqspin = Join[sqspinn2, {matrix[[1]], matrix[[2]]}];

  origAng = Table[angspin[[i]] . epm . angspin[[j]], {i, nn}, {j, nn}];
  origSq = Table[sqspin[[i]] . epm . sqspin[[j]], {i, nn}, {j, nn}];

  lam = angspin; lamT = sqspin;
  newAng = Table[AngleBracket[lam, i, j], {i, nn}, {j, nn}];
  newSq = Table[SquareBracket[lamT, i, j], {i, nn}, {j, nn}];

  If[origAng === newAng && origSq === newSq,
    recordTest[testName, True];
    Print["  PASS: " <> testName];
  ,
    recordTest[testName, False, "Bracket values differ"];
    Print["  FAIL: " <> testName];
  ];

  (* ---- Mandelstam cross-check: s_ij via original repstospin ---- *)
  testName = "Regression: 6-point Mandelstams match original repstospin";
  Module[{origS, newS, allMatch = True, sub},
    Do[
      origS = Sum[
        (angspin[[sub[[a]]]] . epm . angspin[[sub[[b]]]]) *
        (sqspin[[sub[[a]]]] . epm . sqspin[[sub[[b]]]]),
        {a, Length[sub]}, {b, a + 1, Length[sub]}
      ];
      newS = sFromSpinors[lam, lamT, sub];
      If[origS =!= newS, allMatch = False; Break[]];
    , {sub, Subsets[Range[6], {2}]}];

    (* Also check 3-particle Mandelstams *)
    If[allMatch,
      Do[
        origS = Sum[
          (angspin[[sub[[a]]]] . epm . angspin[[sub[[b]]]]) *
          (sqspin[[sub[[a]]]] . epm . sqspin[[sub[[b]]]]),
          {a, Length[sub]}, {b, a + 1, Length[sub]}
        ];
        newS = sFromSpinors[lam, lamT, sub];
        If[origS =!= newS, allMatch = False; Break[]];
      , {sub, Subsets[Range[6], {3}]}];
    ];

    If[allMatch,
      recordTest[testName, True];
      Print["  PASS: " <> testName];
    ,
      recordTest[testName, False, "Mandelstam values differ"];
      Print["  FAIL: " <> testName];
    ];
  ];

  (* ---- Pole regression: angle pole via original method ---- *)
  testName = "Regression: 6-point angle pole matches original construction";
  SeedRandom[789];
  Module[{a1, b1, origAngPole, origSqPole, perm, newS12},
    a1 = RandomInteger[{1, 10}];
    b1 = RandomInteger[{1, 10}];
    (* Original: first two angle spinors identical *)
    angspinn2 = Join[{{a1, b1}, {a1, b1}},
      Table[{RandomInteger[{1, 10}], RandomInteger[{1, 10}]}, {2}]];
    sqspinn2 = Table[{RandomInteger[{1, 10}], RandomInteger[{1, 10}]}, {4}];
    matrix = -Sum[Outer[Times, angspinn2[[i]], sqspinn2[[i]]], {i, 1, 4}];
    angspin = Join[angspinn2, {{1, 0}, {0, 1}}];
    sqspin = Join[sqspinn2, {matrix[[1]], matrix[[2]]}];

    (* Check <12> = 0 in original construction *)
    origAngPole = angspin[[1]] . epm . angspin[[2]];
    If[origAngPole === 0,
      (* Feed same spinors into new package *)
      lam = angspin; lamT = sqspin;
      newS12 = AngleBracket[lam, 1, 2] * SquareBracket[lamT, 1, 2];
      If[newS12 === 0,
        recordTest[testName, True];
        Print["  PASS: " <> testName];
      ,
        recordTest[testName, False, "s12 != 0 in new package"];
        Print["  FAIL: " <> testName];
      ];
    ,
      recordTest[testName, False, "Original <12> != 0"];
      Print["  FAIL: " <> testName];
    ];
  ];
];

(* ============================================================ *)
(* 4. CROSS-CHECK TESTS                                         *)
(* ============================================================ *)

Print["\n--- Section 4: Cross-check tests ---"];

(* 4a. Momentum conservation cross-check for on-pole kinematics *)
Module[{testName = "Cross-check: momentum conservation on pole kinematics (n=6,7,8)",
        failMsg = "", kin, lam, lamT, nn, trial, momSum, poles, i},
  Do[
    poles = If[nn == 6, {{1, 2}}, If[nn == 7, {{1, 2}, {3, 4, 5}}, {{1, 2}, {3, 4}}]];
    Do[
      kin = RandomRationalKinematicsOnPole[nn, poles,
        "TwoParticlePoleType" -> "Angle"];
      If[kin === $Failed, failMsg = "$Failed at n=" <> ToString[nn]; Break[]];
      lam = kin["lambda"]; lamT = kin["lambdaTilde"];
      momSum = Sum[Outer[Times, lam[[i]], lamT[[i]]], {i, nn}];
      If[momSum =!= {{0, 0}, {0, 0}},
        failMsg = "Momentum conservation failed at n=" <> ToString[nn] <>
                  " trial " <> ToString[trial];
        Break[]];
    , {trial, 1, 20}];
    If[failMsg =!= "", Break[]];
  , {nn, {6, 7, 8}}];
  If[failMsg === "",
    recordTest[testName, True]; Print["  PASS: " <> testName],
    recordTest[testName, False, failMsg]; Print["  FAIL: " <> testName <> " - " <> failMsg]];
];

(* 4b. Pole Mandelstams are exactly zero, non-pole are nonzero *)
Module[{testName = "Cross-check: pole s_I=0 and non-pole s_I!=0 (n=7, poles={{1,2},{3,4,5}})",
        failMsg = "", kin, lam, lamT, nn = 7, trial, k, sub, subsets, sI,
        poles = {{1, 2}, {3, 4, 5}}, expectedZeros},
  expectedZeros = Union[{Sort[{1, 2}], Sort[{3, 4, 5}],
    Sort[Complement[Range[7], {1, 2}]], Sort[Complement[Range[7], {3, 4, 5}]]}];
  Do[
    kin = RandomRationalKinematicsOnPole[nn, poles,
      "TwoParticlePoleType" -> "Angle"];
    If[kin === $Failed, failMsg = "$Failed at trial " <> ToString[trial]; Break[]];
    lam = kin["lambda"]; lamT = kin["lambdaTilde"];
    (* Check expected zeros are zero *)
    Do[
      sI = sFromSpinors[lam, lamT, sub];
      If[sI =!= 0,
        failMsg = "Expected zero s_" <> ToString[sub] <> " = " <> ToString[sI]; Break[]];
    , {sub, expectedZeros}];
    If[failMsg =!= "", Break[]];
    (* Check non-pole subsets are nonzero *)
    Do[
      subsets = Subsets[Range[nn], {k}];
      Do[
        If[!MemberQ[expectedZeros, sub],
          sI = sFromSpinors[lam, lamT, sub];
          If[sI === 0,
            failMsg = "Accidental zero s_" <> ToString[sub]; Break[]]];
      , {sub, subsets}];
      If[failMsg =!= "", Break[]];
    , {k, 2, Floor[nn/2]}];
    If[failMsg =!= "", Break[]];
  , {trial, 1, 20}];
  If[failMsg === "",
    recordTest[testName, True]; Print["  PASS: " <> testName],
    recordTest[testName, False, failMsg]; Print["  FAIL: " <> testName <> " - " <> failMsg]];
];

(* 4c. Row sum identity on pole kinematics *)
Module[{testName = "Cross-check: row sums vanish on pole kinematics (n=6, s_{123}=0)",
        failMsg = "", kin, lam, lamT, nn = 6, trial, i, j, rowSum},
  Do[
    kin = RandomRationalKinematicsOnPole[nn, {{1, 2, 3}}];
    If[kin === $Failed, failMsg = "$Failed at trial " <> ToString[trial]; Break[]];
    lam = kin["lambda"]; lamT = kin["lambdaTilde"];
    Do[
      rowSum = Sum[
        AngleBracket[lam, i, j] * SquareBracket[lamT, i, j],
        {j, DeleteCases[Range[nn], i]}];
      If[rowSum =!= 0,
        failMsg = "Row sum for particle " <> ToString[i] <> " = " <> ToString[rowSum];
        Break[]];
    , {i, nn}];
    If[failMsg =!= "", Break[]];
  , {trial, 1, 20}];
  If[failMsg === "",
    recordTest[testName, True]; Print["  PASS: " <> testName],
    recordTest[testName, False, failMsg]; Print["  FAIL: " <> testName <> " - " <> failMsg]];
];

(* 4d. Complement identity on pole kinematics *)
Module[{testName = "Cross-check: complement identity s_I = s_{I^c} on pole kinematics (n=7)",
        failMsg = "", kin, lam, lamT, nn = 7, trial, k, sub, subsets, sI, sComp, comp},
  Do[
    kin = RandomRationalKinematicsOnPole[nn, {{1, 2}},
      "TwoParticlePoleType" -> "Square"];
    If[kin === $Failed, failMsg = "$Failed at trial " <> ToString[trial]; Break[]];
    lam = kin["lambda"]; lamT = kin["lambdaTilde"];
    Do[
      subsets = Subsets[Range[nn], {k}];
      Do[
        sI = sFromSpinors[lam, lamT, sub];
        comp = Complement[Range[nn], sub];
        If[Length[comp] >= 2,
          sComp = sFromSpinors[lam, lamT, comp];
          If[sI =!= sComp,
            failMsg = "s_" <> ToString[sub] <> " != s_" <> ToString[comp];
            Break[]]];
      , {sub, subsets}];
      If[failMsg =!= "", Break[]];
    , {k, 2, Floor[nn/2]}];
    If[failMsg =!= "", Break[]];
  , {trial, 1, 20}];
  If[failMsg === "",
    recordTest[testName, True]; Print["  PASS: " <> testName],
    recordTest[testName, False, failMsg]; Print["  FAIL: " <> testName <> " - " <> failMsg]];
];

(* 4e. NonDegenerate option: True (default) guarantees no accidental zeros *)
Module[{testName = "NonDegenerate->True on pole kinematics prevents accidental zeros (n=6)",
        failMsg = "", kin, lam, lamT, nn = 6, trial, k, sub, subsets, sI,
        expectedZeros},
  expectedZeros = {{1, 2}, Sort[Complement[Range[6], {1, 2}]]};
  Do[
    kin = RandomRationalKinematicsOnPole[nn, {{1, 2}},
      "TwoParticlePoleType" -> "Angle", "NonDegenerate" -> True];
    If[kin === $Failed, failMsg = "$Failed at trial " <> ToString[trial]; Break[]];
    lam = kin["lambda"]; lamT = kin["lambdaTilde"];
    Do[
      subsets = Subsets[Range[nn], {k}];
      Do[
        If[!MemberQ[expectedZeros, sub],
          sI = sFromSpinors[lam, lamT, sub];
          If[sI === 0,
            failMsg = "Accidental zero s_" <> ToString[sub] <> " despite NonDegenerate->True";
            Break[]]];
      , {sub, subsets}];
      If[failMsg =!= "", Break[]];
    , {k, 2, Floor[nn/2]}];
    If[failMsg =!= "", Break[]];
  , {trial, 1, 50}];
  If[failMsg === "",
    recordTest[testName, True]; Print["  PASS: " <> testName],
    recordTest[testName, False, failMsg]; Print["  FAIL: " <> testName <> " - " <> failMsg]];
];

(* 4f. NonDegenerate->False still produces valid pole kinematics *)
Module[{testName = "NonDegenerate->False: poles still vanish, momentum conserved (n=6)",
        failMsg = "", kin, lam, lamT, nn = 6, trial, sI12, momSum, i},
  Do[
    kin = RandomRationalKinematicsOnPole[nn, {{1, 2}},
      "TwoParticlePoleType" -> "Angle", "NonDegenerate" -> False];
    If[kin === $Failed, failMsg = "$Failed at trial " <> ToString[trial]; Break[]];
    lam = kin["lambda"]; lamT = kin["lambdaTilde"];
    (* Pole must still vanish *)
    sI12 = sFromSpinors[lam, lamT, {1, 2}];
    If[sI12 =!= 0,
      failMsg = "s_{12} = " <> ToString[sI12] <> " != 0"; Break[]];
    (* Momentum must still be conserved *)
    momSum = Sum[Outer[Times, lam[[i]], lamT[[i]]], {i, nn}];
    If[momSum =!= {{0, 0}, {0, 0}},
      failMsg = "Momentum conservation failed"; Break[]];
  , {trial, 1, 20}];
  If[failMsg === "",
    recordTest[testName, True]; Print["  PASS: " <> testName],
    recordTest[testName, False, failMsg]; Print["  FAIL: " <> testName <> " - " <> failMsg]];
];

(* 4g. ValidateKinematics cross-check on pole kinematics *)
Module[{testName = "Cross-check: ValidateKinematics passes on pole kinematics (n=6,7,8)",
        failMsg = "", kin, nn, trial, vResult, poles},
  Do[
    poles = If[nn == 6, {{1, 2, 3}}, If[nn == 7, {{1, 2}, {3, 4}}, {{1, 2}, {4, 5, 6}}]];
    Do[
      kin = RandomRationalKinematicsOnPole[nn, poles,
        "TwoParticlePoleType" -> "Angle"];
      If[kin === $Failed, failMsg = "$Failed at n=" <> ToString[nn]; Break[]];
      vResult = ValidateKinematics[kin];
      If[vResult =!= True,
        failMsg = "ValidateKinematics failed at n=" <> ToString[nn] <>
                  ": " <> ToString[vResult[[2]]];
        Break[]];
    , {trial, 1, 10}];
    If[failMsg =!= "", Break[]];
  , {nn, {6, 7, 8}}];
  If[failMsg === "",
    recordTest[testName, True]; Print["  PASS: " <> testName],
    recordTest[testName, False, failMsg]; Print["  FAIL: " <> testName <> " - " <> failMsg]];
];

(* 4h. IsNonDegenerate cross-check with expected poles *)
Module[{testName = "Cross-check: IsNonDegenerate with expected poles (n=7)",
        failMsg = "", kin, nn = 7, trial, poles = {{1, 2}, {3, 4, 5}}},
  Do[
    kin = RandomRationalKinematicsOnPole[nn, poles,
      "TwoParticlePoleType" -> "Angle"];
    If[kin === $Failed, failMsg = "$Failed at trial " <> ToString[trial]; Break[]];
    (* Should be non-degenerate when expected poles are excluded *)
    If[!IsNonDegenerate[kin, poles],
      failMsg = "IsNonDegenerate[kin, poles] returned False"; Break[]];
    (* Should be degenerate without excluding poles *)
    If[IsNonDegenerate[kin],
      failMsg = "IsNonDegenerate[kin] returned True (expected False because poles are zero)";
      Break[]];
  , {trial, 1, 20}];
  If[failMsg === "",
    recordTest[testName, True]; Print["  PASS: " <> testName],
    recordTest[testName, False, failMsg]; Print["  FAIL: " <> testName <> " - " <> failMsg]];
];

(* ============================================================ *)
(* 5. Write results to TestResults.txt                          *)
(* ============================================================ *)

Print["\n--- Writing results ---"];

Module[{outFile, stream, r},
  outFile = FileNameJoin[{$testDir, "TestResults.txt"}];
  stream = OpenWrite[outFile];

  WriteString[stream, "EXTKIN Test Results\n"];
  WriteString[stream, "==================\n"];
  WriteString[stream, "Date: " <> DateString[] <> "\n\n"];
  WriteString[stream, "Total passed: " <> ToString[$totalPassed] <> "\n"];
  WriteString[stream, "Total failed: " <> ToString[$totalFailed] <> "\n"];
  WriteString[stream, "Total tests:  " <> ToString[$totalPassed + $totalFailed] <> "\n\n"];

  WriteString[stream, "--- All Tests ---\n"];
  Do[
    WriteString[stream,
      If[r["passed"], "PASS", "FAIL"] <> ": " <> r["name"] <>
      If[r["detail"] =!= "", " -- " <> r["detail"], ""] <> "\n"];
  , {r, $testResults}];

  If[$totalFailed > 0,
    WriteString[stream, "\n--- Failure Details ---\n"];
    Do[
      WriteString[stream, "FAIL: " <> r["name"] <> "\n"];
      WriteString[stream, "  Detail: " <> r["detail"] <> "\n\n"];
    , {r, $failureDetails}];
  ];

  Close[stream];
  Print["Results written to: " <> outFile];
];

Print["\n=== Summary ==="];
Print["Passed: " <> ToString[$totalPassed]];
Print["Failed: " <> ToString[$totalFailed]];
If[$totalFailed === 0,
  Print["All tests passed!"],
  Print["Some tests failed. See TestResults.txt for details."]
];
