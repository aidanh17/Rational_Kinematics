(* Poles.m -- Generation of rational kinematics on poles *)
(* Part of the EXTKIN package *)

(* ---- Union-Find for connected components ---- *)

connectedComponents[pairs_List, n_Integer] := Module[
  {parent, i, p, r1, r2, root, groups},

  parent = Range[n];

  (* Find root with path compression *)
  root[x_Integer] := Module[{r = x, c, next},
    While[parent[[r]] =!= r, r = parent[[r]]];
    c = x;
    While[c =!= r, next = parent[[c]]; parent[[c]] = r; c = next];
    r
  ];

  (* Union step *)
  Do[
    r1 = root[p[[1]]];
    r2 = root[p[[2]]];
    If[r1 =!= r2, parent[[r1]] = r2];
  , {p, pairs}];

  Do[root[i], {i, n}];

  groups = GatherBy[Range[n], (parent[[#]] &)];
  Select[groups, Length[#] > 1 &]
];

(* ---- Parse two-particle pole type specifications ---- *)

parseTwoParticlePoleTypes[twoPoles_List, poleTypeOpt_] := Module[
  {types = <||>, pole, ptype, match},

  Do[
    ptype = Which[
      StringQ[poleTypeOpt] && poleTypeOpt === "Angle", "Angle",
      StringQ[poleTypeOpt] && poleTypeOpt === "Square", "Square",
      StringQ[poleTypeOpt] && poleTypeOpt === "Random",
        RandomChoice[{"Angle", "Square"}],

      ListQ[poleTypeOpt],
        match = Select[poleTypeOpt,
          (Head[#] === Rule &&
           (#[[1]] === pole || #[[1]] === Reverse[pole])) &];
        If[Length[match] > 0,
          match[[1, 2]],
          RandomChoice[{"Angle", "Square"}]
        ],

      True, RandomChoice[{"Angle", "Square"}]
    ];
    types[pole] = ptype;
  , {pole, twoPoles}];

  types
];

(* ---- BCFW shift to set a multi-particle Mandelstam to zero ---- *)

applyBCFWShiftForPole[lam_, lamT_, poleSet_List, n_Integer,
    anglePairs_List, squarePairs_List, allMultiPoles_List] :=
  Module[
    {complement, a, b, sI, denom, tStar, newLam, newLamT,
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
    If[Length[candidateB] === 0,
      candidateB = Complement[complement, angGroupParticles]];
    If[Length[candidateB] === 0,
      candidateB = complement];

    sI = computeMandelstamSubset[lam, lamT, poleSet];
    If[sI === 0, Return[{lam, lamT}, Module]];

    Do[
      Do[
        denom = Sum[
          AngleBracket[lam, a, j] * SquareBracket[lamT, b, j],
          {j, DeleteCases[poleSet, a]}
        ];

        If[denom =!= 0,
          tStar = -sI / denom;
          newLamT = lamT;
          newLam = lam;
          newLamT[[a]] = lamT[[a]] + tStar * lamT[[b]];
          newLam[[b]] = lam[[b]] - tStar * lam[[a]];
          Return[{newLam, newLamT}, Module];
        ];
      , {b, candidateB}];
    , {a, candidateA}];

    $Failed
  ];

(* ---- Validation for pole kinematics ---- *)

validatePoleKinematics[lam_, lamT_, n_, poles_, poleTypes_, twoPoles_] :=
  Module[
    {expectedZeroSubsets, complementSub, sub, sI, k, subsets,
     pole, i, j, ab, sb},

    expectedZeroSubsets = {};
    Do[
      AppendTo[expectedZeroSubsets, Sort[pole]];
      complementSub = Complement[Range[n], pole];
      If[Length[complementSub] >= 2 && Length[complementSub] <= n - 2,
        AppendTo[expectedZeroSubsets, Sort[complementSub]];
      ];
    , {pole, poles}];
    expectedZeroSubsets = Union[expectedZeroSubsets];

    Do[
      sI = computeMandelstamSubset[lam, lamT, sub];
      If[sI =!= 0, Return[False, Module]];
    , {sub, expectedZeroSubsets}];

    Do[
      If[KeyExistsQ[poleTypes, Sort[pole]] && poleTypes[Sort[pole]] === "Angle",
        {i, j} = Sort[pole];
        ab = AngleBracket[lam, i, j];
        sb = SquareBracket[lamT, i, j];
        If[ab =!= 0 || sb === 0, Return[False, Module]];
      ];
      If[KeyExistsQ[poleTypes, Sort[pole]] && poleTypes[Sort[pole]] === "Square",
        {i, j} = Sort[pole];
        ab = AngleBracket[lam, i, j];
        sb = SquareBracket[lamT, i, j];
        If[sb =!= 0 || ab === 0, Return[False, Module]];
      ];
    , {pole, twoPoles}];

    Do[
      subsets = Subsets[Range[n], {k}];
      Do[
        If[!MemberQ[expectedZeroSubsets, sub],
          sI = computeMandelstamSubset[lam, lamT, sub];
          If[sI === 0, Return[False, Module]];
        ];
      , {sub, subsets}];
    , {k, 2, Floor[n/2]}];

    Do[
      If[!MemberQ[twoPoles, Sort[{i, j}]],
        If[AngleBracket[lam, i, j] === 0, Return[False, Module]];
        If[SquareBracket[lamT, i, j] === 0, Return[False, Module]];
      ];
    , {i, n}, {j, i + 1, n}];

    True
  ];

(* ---- Build kinematics association ---- *)

buildKinematicsAssociation[lam_, lamT_, n_] := Module[
  {momenta, mandelstams, angleBrackets, squareBrackets},

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

(* ---- Main pole generation function ---- *)

Options[RandomRationalKinematicsOnPole] = {
  "TwoParticlePoleType" -> "Random",
  "Range" -> 9,
  "MaxAttempts" -> 1000
};

RandomRationalKinematicsOnPole[n_Integer, poles_List, opts : OptionsPattern[]] :=
  Module[
    {range, maxAttempts, poleTypeOpt,
     twoPoles, multiPoles, poleTypes,
     anglePairs, squarePairs, angleGroups, squareGroups,
     allParticlesInAngleGroups, allParticlesInSquareGroups,
     freeAngParticles, freeAllParticles, solveParticles,
     lam, lamT, matrix, group, refSpinor, c, i, j,
     shiftResult, valid, attempt, pole,
     sortedPoles, complementPole,
     (* Effective-spinor solve variables *)
     solveInSquareGroups, nonSolveSquareGroups, allSolveGroupParticles,
     solveGroupScales, effectiveLam, coeffMat, sol, rhsParticles,
     sp, m},

    (* Parse options *)
    range = OptionValue[RandomRationalKinematicsOnPole, {opts}, "Range"];
    maxAttempts = OptionValue[RandomRationalKinematicsOnPole, {opts}, "MaxAttempts"];
    poleTypeOpt = OptionValue[RandomRationalKinematicsOnPole, {opts}, "TwoParticlePoleType"];

    (* Sort and validate poles *)
    sortedPoles = Sort /@ poles;
    Do[
      If[Length[pole] < 2 || Length[pole] > n - 2,
        Message[RandomRationalKinematicsOnPole::failed, 0,
          "Pole subset " <> ToString[pole] <> " must have 2 to " <> ToString[n-2] <> " elements"];
        Return[$Failed]
      ];
      If[!SubsetQ[Range[n], pole],
        Message[RandomRationalKinematicsOnPole::failed, 0,
          "Pole subset " <> ToString[pole] <> " contains invalid particle indices"];
        Return[$Failed]
      ];
    , {pole, sortedPoles}];

    (* Remove redundant complement poles *)
    Do[
      complementPole = Sort[Complement[Range[n], pole]];
      If[MemberQ[sortedPoles, complementPole] && pole =!= complementPole,
        sortedPoles = DeleteCases[sortedPoles, complementPole, 1, 1]
      ];
    , {pole, sortedPoles}];
    sortedPoles = Union[sortedPoles];

    (* Classify poles *)
    twoPoles = Select[sortedPoles, Length[#] == 2 &];
    multiPoles = Select[sortedPoles, Length[#] >= 3 &];

    (* Determine pole types for two-particle poles *)
    poleTypes = parseTwoParticlePoleTypes[twoPoles, poleTypeOpt];

    anglePairs = Select[twoPoles, poleTypes[#] === "Angle" &];
    squarePairs = Select[twoPoles, poleTypes[#] === "Square" &];

    (* Build proportionality groups *)
    angleGroups = connectedComponents[anglePairs, n];
    squareGroups = connectedComponents[squarePairs, n];

    allParticlesInAngleGroups = If[Length[angleGroups] > 0, Union @@ angleGroups, {}];
    allParticlesInSquareGroups = If[Length[squareGroups] > 0, Union @@ squareGroups, {}];

    (* ---- Choose solve particles ---- *)
    (* Must not be in any angle group (need independent angle spinors).
       Prefer free from all groups. If not enough, pick from different square groups. *)
    freeAngParticles = Complement[Range[n], allParticlesInAngleGroups];
    If[Length[freeAngParticles] < 2,
      Message[RandomRationalKinematicsOnPole::failed, 0,
        "Too many angle-type poles: cannot find 2 particles with independent angle spinors"];
      Return[$Failed]
    ];

    freeAllParticles = Complement[freeAngParticles, allParticlesInSquareGroups];

    solveParticles = Which[
      (* Best: both free from all groups *)
      Length[freeAllParticles] >= 2,
        Take[freeAllParticles, 2],

      (* One free, one from a square group *)
      Length[freeAllParticles] >= 1,
        {freeAllParticles[[1]],
         First[Complement[freeAngParticles, freeAllParticles]]},

      (* Both from square groups -- pick from different groups *)
      True,
        Module[{reps = {}},
          Do[
            Module[{rep = Select[group, MemberQ[freeAngParticles, #] &]},
              If[Length[rep] > 0 && Length[reps] < 2,
                AppendTo[reps, rep[[1]]]]];
          , {group, squareGroups}];
          If[Length[reps] >= 2, reps, Take[freeAngParticles, 2]]
        ]
    ];

    (* ---- Identify solve square groups (groups containing a solve particle) ---- *)
    solveInSquareGroups = <||>; (* sp -> {group, otherMembers} *)
    Do[
      Module[{myGroups = Select[squareGroups, MemberQ[#, sp] &]},
        If[Length[myGroups] > 0,
          solveInSquareGroups[sp] = {myGroups[[1]], DeleteCases[myGroups[[1]], sp]};
        ];
      ];
    , {sp, solveParticles}];

    nonSolveSquareGroups = Select[squareGroups,
      Length[Intersection[#, solveParticles]] === 0 &];

    allSolveGroupParticles = {};
    Do[
      If[KeyExistsQ[solveInSquareGroups, sp],
        allSolveGroupParticles = Union[allSolveGroupParticles, solveInSquareGroups[sp][[1]]];
      ];
    , {sp, solveParticles}];

    (* ---- Main attempt loop ---- *)
    Do[
      (* ---- Generate angle spinors ---- *)
      lam = Table[{0, 0}, {n}];

      Do[
        refSpinor = {RandomInteger[{1, range}], RandomInteger[{1, range}]};
        lam[[group[[1]]]] = refSpinor;
        Do[
          c = RandomChoice[Join[Range[-3, -1], Range[1, 3]]];
          lam[[group[[j]]]] = c * refSpinor;
        , {j, 2, Length[group]}];
      , {group, angleGroups}];

      Do[
        lam[[i]] = {RandomInteger[{1, range}], RandomInteger[{1, range}]};
      , {i, Complement[Range[n], allParticlesInAngleGroups, solveParticles]}];

      lam[[solveParticles[[1]]]] = {1, 0};
      lam[[solveParticles[[2]]]] = {0, 1};

      (* ---- Generate square spinors ---- *)
      lamT = Table[{0, 0}, {n}];

      (* Non-solve square groups: generate proportional spinors as usual *)
      Do[
        refSpinor = {RandomInteger[{1, range}], RandomInteger[{1, range}]};
        lamT[[group[[1]]]] = refSpinor;
        Do[
          c = RandomChoice[Join[Range[-3, -1], Range[1, 3]]];
          lamT[[group[[j]]]] = c * refSpinor;
        , {j, 2, Length[group]}];
      , {group, nonSolveSquareGroups}];

      (* Free non-solve particles not in any square group *)
      Do[
        lamT[[i]] = {RandomInteger[{1, range}], RandomInteger[{1, range}]};
      , {i, Complement[Range[n], allParticlesInSquareGroups, solveParticles]}];

      (* NOTE: Solve particles and their square-group members have lamT = {0,0}
         for now. They will be set by the effective-spinor momentum conservation solve. *)

      (* ---- Generate proportionality constants for solve square groups ---- *)
      solveGroupScales = <||>; (* member -> {solveParticle, scale} *)
      Do[
        If[KeyExistsQ[solveInSquareGroups, sp],
          Module[{others = solveInSquareGroups[sp][[2]]},
            Do[
              c = RandomChoice[Join[Range[-3, -1], Range[1, 3]]];
              solveGroupScales[m] = {sp, c};
            , {m, others}];
          ];
        ];
      , {sp, solveParticles}];

      (* ---- Compute effective angle spinors for solve particles ---- *)
      (* If solve particle sp is in square group G, the effective spinor is:
         A_sp = lambda_sp + Sum_{m in G\{sp}} c_m * lambda_m
         This accounts for the fact that lambdaTilde_m = c_m * lambdaTilde_sp *)
      effectiveLam = lam;
      Do[
        If[KeyExistsQ[solveInSquareGroups, sp],
          Module[{others = solveInSquareGroups[sp][[2]], effLam},
            effLam = lam[[sp]];
            Do[
              effLam = effLam + solveGroupScales[m][[2]] * lam[[m]];
            , {m, others}];
            effectiveLam[[sp]] = effLam;
          ];
        ];
      , {sp, solveParticles}];

      (* ---- Solve momentum conservation using effective spinors ---- *)
      (* Equation: A_{s1} (x) lamT_{s1} + A_{s2} (x) lamT_{s2} = -RHS
         where RHS = Sum over particles NOT in any solve group *)
      rhsParticles = Complement[Range[n], solveParticles, allSolveGroupParticles];
      matrix = If[Length[rhsParticles] > 0,
        Sum[Outer[Times, lam[[i]], lamT[[i]]], {i, rhsParticles}],
        {{0, 0}, {0, 0}}
      ];

      coeffMat = Transpose[{effectiveLam[[solveParticles[[1]]]],
                             effectiveLam[[solveParticles[[2]]]]}];

      If[Det[coeffMat] === 0, Continue[]]; (* Degenerate effective spinors, retry *)

      sol = Inverse[coeffMat] . (-matrix);
      lamT[[solveParticles[[1]]]] = sol[[1]];
      lamT[[solveParticles[[2]]]] = sol[[2]];

      (* Set dependent square spinors for solve square groups *)
      Do[
        lamT[[m]] = solveGroupScales[m][[2]] * lamT[[solveGroupScales[m][[1]]]];
      , {m, Keys[solveGroupScales]}];

      (* ---- Apply BCFW shifts for multi-particle poles ---- *)
      valid = True;
      Do[
        shiftResult = applyBCFWShiftForPole[lam, lamT, pole, n,
          anglePairs, squarePairs, multiPoles];
        If[shiftResult === $Failed,
          valid = False;
          Break[];
        ];
        {lam, lamT} = shiftResult;
      , {pole, multiPoles}];

      If[!valid, Continue[]];

      (* ---- Validate all constraints ---- *)
      If[validatePoleKinematics[lam, lamT, n, sortedPoles, poleTypes, twoPoles],
        Return[
          buildKinematicsAssociation[lam, lamT, n],
          Module
        ]
      ];
    , {attempt, maxAttempts}];

    Message[RandomRationalKinematicsOnPole::failed, maxAttempts,
      "Could not simultaneously satisfy all pole constraints"];
    $Failed
  ];
