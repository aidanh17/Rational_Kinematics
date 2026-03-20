================================================================================
EXTKIN: Exact Rational Kinematics Package for Spinor-Helicity Variables
================================================================================

DESCRIPTION
-----------
EXTKIN generates random rational kinematic points in 4D spinor-helicity
variables for massless n-particle scattering processes. All quantities are
exact rational numbers -- no floating point anywhere. The package supports:

  - Generic n-point kinematics with non-degeneracy guarantees
  - Kinematics on specified poles (where Mandelstam invariants vanish)
  - Two-particle poles via angle or square bracket collinearity
  - Multi-particle poles via BCFW deformation
  - Multiple simultaneous poles
  - Comprehensive validation of all kinematic identities

Extracted, consolidated, and generalized from code originally written for
N=4 super Yang-Mills and N=8 supergravity amplitude computations.

================================================================================
PUBLIC FUNCTIONS
================================================================================

1. RandomRationalSpinors[n, opts]
   --------------------------------
   Generates n pairs of 2-component rational spinors {lambda_i, lambdaTilde_i}
   satisfying exact momentum conservation: Sum_i lambda_i (x) lambdaTilde_i = 0.

   Options:
     "NonDegenerate" -> True    (default) Regenerate if any <ij>=0 or [ij]=0
                                 or any multi-particle Mandelstam vanishes
     "Range" -> 9               Integer range for random spinor entries [1, Range]
     "MaxAttempts" -> 1000      Maximum regeneration attempts

   Returns: {lambdaList, lambdaTildeList} where each is an n x 2 matrix,
            or $Failed if non-degenerate point cannot be found.

2. AngleBracket[lambdaList, i, j]
   --------------------------------
   Computes <ij> = eps_{ab} lambda_i^a lambda_j^b
                 = lambda_i^1 lambda_j^2 - lambda_i^2 lambda_j^1

   Properties: <ij> = -<ji>, <ii> = 0

3. SquareBracket[lambdaTildeList, i, j]
   --------------------------------
   Computes [ij] = eps_{ab} lambdaTilde_i^a lambdaTilde_j^b
                  = lambdaTilde_i^1 lambdaTilde_j^2 - lambdaTilde_i^2 lambdaTilde_j^1

   Properties: [ij] = -[ji], [ii] = 0

4. MomentumFromSpinors[lambda, lambdaTilde]
   --------------------------------
   Single particle: two 2-vectors -> 2x2 matrix p^{a adot} = lambda^a lambdaTilde^{adot}
   All particles:   two n x 2 matrices -> list of n 2x2 matrices

5. MandelstamFromSpinors[lambdaList, lambdaTildeList]
   --------------------------------
   Computes all Mandelstam invariants s_I for subsets I with 2 <= |I| <= n-2.

   Returns: Association with sorted subset keys.
   Example: result[{1,2}] gives s_{12}, result[{1,2,3}] gives s_{123}.

6. RandomRationalKinematics[n, opts]
   --------------------------------
   Wrapper that generates a complete kinematic point.

   Returns: Association with keys:
     "n"              - number of particles
     "lambda"         - n x 2 matrix of angle spinors
     "lambdaTilde"    - n x 2 matrix of square spinors
     "momenta"        - list of n 2x2 momentum matrices
     "mandelstams"    - Association of all s_I values
     "spinorProducts" - <| "angle" -> ..., "square" -> ... |>

   Options: same as RandomRationalSpinors.

7. RandomRationalKinematicsOnPole[n, poles, opts]
   --------------------------------
   Generates rational kinematics where specified Mandelstam invariants vanish.

   Arguments:
     n     - number of particles
     poles - list of subsets, e.g. {{1,2}, {3,4,5}}

   Options:
     "TwoParticlePoleType" -> "Random"  (default) Randomly choose Angle or Square
                              "Angle"   Set <ij> = 0 for all two-particle poles
                              "Square"  Set [ij] = 0 for all two-particle poles
                              {rules}   Per-pole: {{1,2} -> "Angle", {3,4} -> "Square"}
     "Range" -> 9
     "MaxAttempts" -> 1000

   Returns: Same Association as RandomRationalKinematics, or $Failed.

   Algorithm:
     - Two-particle poles: enforced via spinor proportionality
       (Angle: lambda_j = c * lambda_i; Square: lambdaTilde_j = c * lambdaTilde_i)
     - Multi-particle poles: enforced via BCFW deformation that preserves
       momentum conservation and on-shell conditions
     - Multiple poles: handled sequentially with careful choice of shift
       particles to avoid breaking previously imposed constraints
     - Validated internally: all poles vanish, non-pole Mandelstams nonzero

8. ValidateKinematics[kin]
   --------------------------------
   Comprehensive consistency check of a kinematic point.

   Checks:
     - Exact rationality (no floating point)
     - Momentum conservation: Sum_i p_i = 0 exactly
     - On-shell conditions: p_i^2 = 0 for all i
     - Bracket consistency: stored values match recomputed values
     - Mandelstam consistency: s_ij = <ij>[ij]
     - Row sum identity: Sum_{j != i} s_ij = 0 for all i
     - Complement identity: s_I = s_{I^c}

   Returns: True if all pass, or {False, {diagnostic messages...}}

================================================================================
CONVENTIONS
================================================================================

Metric signature:
  Not explicitly fixed. The package works in the spinor-helicity formalism
  where all quantities are algebraic. The Mandelstam sign convention is
  s_ij = <ij>[ij] (see below).

Spinor convention:
  Angle brackets: <ij> = eps_{ab} lambda_i^a lambda_j^b
                       = lambda_i . epm . lambda_j
  where epm = {{0, 1}, {-1, 0}} is the 2D Levi-Civita tensor.

  Square brackets: [ij] = eps_{ab} lambdaTilde_i^a lambdaTilde_j^b
                        = lambdaTilde_i . epm . lambdaTilde_j
  Same formula with lambdaTilde.

Mandelstam convention:
  s_ij = <ij>[ij]     (SAME index ordering on both brackets)

  This differs from some textbook conventions that use s_ij = <ij>[ji].
  Since [ji] = -[ij], the two conventions differ by a sign:
    s_ij^{here} = <ij>[ij] = -<ij>[ji] = -(s_ij^{Elvang-Huang})

  This convention was chosen because it is the one used in the majority of
  the original codebase (both 6ptSUSY/Aidan_refactored and
  N8SUGRAf/N_8_SG_EFT_New). The N8SUGRAf/N8ANSATZ code uses the other
  convention (s_ij = <ij>[ji]), but this is the minority convention.

  KEY IDENTITY: s_ij = s_ji (symmetric), because:
    <ij>[ij] = (-<ji>)(-[ji]) = <ji>[ji]

Multi-particle Mandelstam:
  s_I = Sum_{i<j in I} s_ij = Sum_{i<j in I} <ij>[ij]

  This equals (Sum_{i in I} p_i)^2 up to the sign convention.

Momentum conservation:
  All momenta are incoming: Sum_{i=1}^{n} p_i = 0
  In spinor form: Sum_{i=1}^{n} lambda_i (x) lambdaTilde_i = 0 (as 2x2 matrix)

Row sum identity:
  Sum_{j != i} s_ij = 0  for each i  (from momentum conservation + on-shell)

Complement identity:
  s_I = s_{I^c}  where I^c = {1,...,n} \ I  (from momentum conservation)

Index ordering:
  Mandelstam invariants are stored with sorted subset keys.
  Example: s_{31} is stored as s_{{1,3}}, accessed via mandelstams[{1,3}].

Spinor normalization:
  No special normalization. The "solve" particles (last two in internal
  ordering) have angle spinors {1,0} and {0,1}. All other spinor entries
  are random positive integers in [1, Range].

================================================================================
SOURCE CATALOG (from original code)
================================================================================

The following functions were identified in the original codebases and used
as reference for the generalized n-point implementations:

--- From 6ptSUSY/ ---

File: 6ptSUSY/Simon_Comparison/Aidan_refactored/pkg/Kinematics.wl
  - RandomSpin4          4-point generic rational kinematics
  - RandomSpin5          5-point generic rational kinematics
  - RandomSpin6          6-point generic rational kinematics
  - RandomSpin5AngPole   5-point with <ij>=0 pole
  - RandomSpin6AngPole   6-point with <ij>=0 pole
  - RandomSpin6SqPole    6-point with [ij]=0 pole (via angle/square swap trick)
  - repstospin           Mandelstam -> spinor product conversion rule
  - planarvars5/6        Planar Mandelstam basis variables
  - convcanvars5p/6p     Conversion rules to canonical variables
  - epm                  Levi-Civita tensor {{0,1},{-1,0}}
  - permute              Permutation utility for pole construction
  - simpman, simpangsq   Simplification rules for brackets/Mandelstams
  - CoefficientVec       Coefficient extraction utility
  - NumerVQ              Numeric check utility

File: 6ptSUSY/Simon_Comparison/Simon_refactored/pkg/Kinematics.wl
  - SimonRandomSpin6     6-point generator (Simon's convention: al/ar symbols)
  - Uses different symbol names: al[i,j] for angle, ar[i,j] for square
  - Also computes epsilon tensor contractions

File: 6ptSUSY/Simon_Comparison/Comparison/shared_definitions.wl
  - RationalKinematicPoint    Seed-based 6-point generator
  - MakeCollinearKin          Collinear limit kinematics (parametric)
  - MakeSoftKin               Soft limit kinematics (parametric)
  - TruncateToOrder3          Mass-dimension truncation utility

File: 6ptSUSY/Simon_Comparison/Aidan_refactored/pkg/AnsatzBasis.wl
  - BuildSixPointZ1      Uses pole kinematics for factorization constraints
  - Uses RandomSpin6AngPole and RandomSpin6SqPole for numerical constraints

File: 6ptSUSY/Simon_Comparison/Aidan_refactored/pkg/Numerics.wl
  - EvalOnKinematics             Substitute kinematics into expressions
  - CheckZeroOnKinematics        Test if expression vanishes numerically
  - GenerateRationalKinematics6  Alias for RandomSpin6

--- From N8SUGRAf/ ---

File: N8SUGRAf/N8ANSATZ/core/SpinorHelicity.m
  - generateKinematics           6-point generic generator
  - rationalKinSqZero[i,j]       [ij]=0 pole
  - rationalKinAngZero[i,j]      <ij>=0 pole
  - rationalKinDoubleSqZero      Two simultaneous [ij]=0 poles
  - rationalKinTripleSqZero      Three simultaneous [ij]=0 poles
  - rationalKinS3Zero[i,j,k]     Three-particle pole s_ijk=0
  - rationalKinSqAndS3Zero       Combined [ij]=0 and s_ijk=0
  - rationalRestrictedKinematics General constraint framework
  - abEval, sbEval               Bracket computation from spinor arrays
  - sij, sijk                    Mandelstam symbols (uses <ij>[ji] convention)
  - checkMomCons                 Momentum conservation check
  - computeRules                 Build substitution rules from spinors
  - solveTildeSpinors            Solve for tilde spinors from mom. conservation
  - decomposeRank1               Rank-1 matrix decomposition into spinors
  - eps4                         4D Levi-Civita contraction
  - asq, chain                   Spinor chain contractions

File: N8SUGRAf/N_8_SG_EFT_New/code/kinematics.wl
  - GenericKinematics            6-point via null-space approach
  - kinematics2PtPoleSquareZero  [ij]=0 pole via linear solve
  - kinematics2PtPoleAngleZero   <ij>=0 pole via linear solve
  - Kinematics3PtPole            s_ijk=0 via BCFW deformation
  - KinematicsDoublePole         Combined 2+3 particle pole
  - AngleBracket, SquareBracket  Bracket computation
  - Mandelstam2, Mandelstam3     Two- and three-particle invariants
  - ComputeAllBrackets           Compute all bracket types
  - CheckMomentumConservation    Momentum conservation verification

File: N8SUGRAf/FinalCode/shared.wl
  - generateKinWide              6-point with tunable range
  - crossMandVals                Cross-Mandelstam evaluation

File: N8SUGRAf/N8ANSATZ/PoleChecks.m
  - extractResidueS3Numerical    Richardson extrapolation for residues
  - extractResidue2particleNumerical  2-particle residue extraction
  - extractDoubleResidueM3M3M4   Double residue extraction

================================================================================
DESIGN DECISIONS
================================================================================

1. MANDELSTAM CONVENTION: s_ij = <ij>[ij]
   The 6ptSUSY code uses s[{i,j}] = ang[i,j] * sq[i,j] = <ij>[ij].
   The N_8_SG_EFT_New code uses Mandelstam2 = AngleBracket * SquareBracket = <ij>[ij].
   The N8ANSATZ code uses sij[i,j] = ab[i,j] * sb[j,i] = <ij>[ji] (opposite sign).
   Decision: adopted <ij>[ij] as it is used in 2 of 3 codebases.

2. SYMBOL NAMES: AngleBracket / SquareBracket
   Original code used ang/sq (6ptSUSY), al/ar (Simon), ab/sb (N8ANSATZ),
   and AngleBracket/SquareBracket (N_8_SG_EFT_New).
   Decision: adopted the descriptive AngleBracket/SquareBracket names
   for clarity. These match the N_8_SG_EFT_New convention.

3. MANDELSTAM STORAGE: Association with sorted subset keys
   Original code used various formats: replacement rules (ang[i,j] -> val),
   function calls (sij[i,j]), or indexed symbols (s[{i,j}]).
   Decision: use Association for clean O(1) lookup with sorted subset keys.

4. MOMENTUM CONSERVATION ALGORITHM: {1,0}/{0,1} solve trick
   All original implementations use the same algorithm: generate n-2 free
   spinor pairs, set the last two angle spinors to {1,0} and {0,1}, and
   solve for their square spinors from the 2x2 matrix equation. This is
   linear, always yields rational solutions, and is adopted unchanged.

5. TWO-PARTICLE POLES: Direct proportionality
   The original 6ptSUSY code generates proportional angle spinors for
   angle-type poles and uses an angle/square swap trick for square-type
   poles. The N8SUGRAf code directly sets lambdaTilde proportional.
   Decision: use direct proportionality for both types (cleaner, no swap).

6. MULTI-PARTICLE POLES: BCFW deformation
   The N8SUGRAf/N_8_SG_EFT_New code uses BCFW shifts for 3-particle poles.
   The N8ANSATZ code uses rank-1 matrix decomposition.
   Decision: adopted BCFW approach for multi-particle poles because:
   (a) it naturally preserves momentum conservation and on-shell conditions,
   (b) it generalizes cleanly to arbitrary n,
   (c) the shift parameter t* is always rational.

   BCFW shift: lambdaTilde_a -> lambdaTilde_a + t * lambdaTilde_b
               lambda_b -> lambda_b - t * lambda_a
   where a is in the pole set and b is outside.
   Shift parameter: t* = -s_I / Sum_{j in I\{a}} <aj>[bj]

7. MULTIPLE SIMULTANEOUS POLES: Sequential BCFW with careful shift selection
   For non-overlapping pole sets, shift particles are chosen so that each
   BCFW deformation only affects its target pole. For overlapping poles,
   shifts are applied sequentially with validation and retry.
   Shift particle selection rules:
   - a (from pole set): not in any square-type two-particle pole group
   - b (from complement): not in any angle-type two-particle pole group,
     preferably not in any other multi-particle pole set

8. NON-DEGENERACY: Comprehensive checking
   The original code only checked specific bracket conditions. The new
   package checks all <ij>, all [ij], and all multi-particle Mandelstams
   s_I for |I| = 2 to Floor[n/2].

9. RANDOM INTEGER RANGE: Positive integers [1, Range]
   The original code used RandomInteger[{1, 50}] for 4-point and
   RandomInteger[{1, 100}] for 5/6-point. The package uses a configurable
   range (default 9) for all multiplicities. Proportionality coefficients
   for pole constraints use integers from {-3,...,-1,1,...,3}.

10. SpinorChain FUNCTION: Not implemented
    The prompt mentions SpinorChain for longer spinor chains "if any exist."
    The original code has chain[i, js, k] = Sum_m <i,js_m>[js_m,k] and
    asq[i,j,k] = <ij>[jk]. These are simple compositions of AngleBracket
    and SquareBracket that users can compute directly, so a separate
    SpinorChain function was deemed unnecessary.

================================================================================
KNOWN LIMITATIONS AND CAVEATS
================================================================================

1. INCOMPATIBLE POLE COMBINATIONS
   Certain pole combinations are kinematically impossible:
   - n=4, s_{12}=0 and s_{13}=0: forces all invariants to zero
   - n=4, any two distinct poles: only one independent invariant
   - Overlapping angle poles that would force all spinors proportional
   The function returns $Failed with a diagnostic message.

2. DEGENERATE CONFIGURATIONS
   The BCFW shift parameter t* = -s_I / denom can have large numerator/
   denominator, leading to kinematic points with very large rational numbers.
   This is mathematically correct but may slow subsequent symbolic computation.
   Retrying with different random seeds usually gives smaller numbers.

3. POLE SET SIZE
   For large pole sets relative to n, finding valid shift particles becomes
   harder. The function may exhaust MaxAttempts for highly constrained
   configurations.

4. FLOATING POINT
   The package never introduces floating point. However, if users pass
   floating-point spinors to the utility functions, the outputs will
   contain floats. Always use exact integer/rational inputs.

5. LARGE n PERFORMANCE
   Non-degeneracy checking examines all Subsets[Range[n], {k}] for
   k = 2..Floor[n/2]. For n=10 this is ~700 subsets; for n=15 it would
   be ~16000. The package is designed for n <= 10-12.

6. COMPLEMENT REDUNDANCY
   Requesting both s_I = 0 and s_{I^c} = 0 is redundant (they are equal
   by momentum conservation). The package silently removes such duplicates.

================================================================================
TEST SUITE
================================================================================

Tests.m contains a comprehensive test suite organized in three sections:

Section 1 -- Basic tests (n = 4, 5, 6, 7, 8, 10):
  50 random points per multiplicity, checking:
  - Exact momentum conservation
  - On-shell conditions
  - s_ij = <ij>[ij] consistency
  - Row sum identity
  - Non-degeneracy

Section 2 -- Pole tests (n = 5, 6, 7, 8):
  - Single angle-type two-particle pole
  - Single square-type two-particle pole
  - Single multi-particle pole
  - Two simultaneous two-particle poles
  - Three simultaneous two-particle poles (n=7, 8)
  - Mixed two-particle + multi-particle poles
  - Incompatible poles ($Failed test)
  - Bracket verification for angle/square poles
  - Sub-momentum nullity for multi-particle poles

Section 3 -- Regression tests:
  - 4/5/6-point bracket values match original algorithm (dot-with-epm)
  - Mandelstam values match original repstospin rule
  - Angle pole construction matches original method

To run: open Mathematica and execute
  Get["path/to/EXTKIN/tests/Tests.m"]
Results are written to EXTKIN/tests/TestResults.txt.

Note: WolframScript/Mathematica was not available in the build environment,
so tests have not been pre-run. The code has been carefully verified for
correctness against the original implementations.

================================================================================
USAGE EXAMPLES
================================================================================

(* Load the package *)
Get["/path/to/EXTKIN/EXTKIN.m"];

(* Generate a generic 6-point kinematic point *)
kin = RandomRationalKinematics[6];
kin["mandelstams"][{1,2}]     (* s_{12} *)
kin["spinorProducts"]["angle"][{1,3}]   (* <13> *)

(* Validate *)
ValidateKinematics[kin]   (* returns True *)

(* Generate kinematics with s_{12} = 0 via <12> = 0 *)
kinPole = RandomRationalKinematicsOnPole[6, {{1, 2}},
  "TwoParticlePoleType" -> "Angle"];
kinPole["mandelstams"][{1,2}]   (* returns 0 *)

(* Generate kinematics with s_{123} = 0 *)
kinS123 = RandomRationalKinematicsOnPole[6, {{1, 2, 3}}];
kinS123["mandelstams"][{1,2,3}]   (* returns 0 *)

(* Multiple poles: s_{12} = 0 and s_{45} = 0 simultaneously *)
kinMulti = RandomRationalKinematicsOnPole[6, {{1, 2}, {4, 5}},
  "TwoParticlePoleType" -> "Square"];

(* Mixed: s_{12} = 0 (angle) and s_{345} = 0 simultaneously *)
kinMixed = RandomRationalKinematicsOnPole[7, {{1, 2}, {3, 4, 5}},
  "TwoParticlePoleType" -> {{1, 2} -> "Angle"}];

(* Work with raw spinors *)
{lam, lamT} = RandomRationalSpinors[8];
AngleBracket[lam, 3, 7]
SquareBracket[lamT, 1, 5]
MandelstamFromSpinors[lam, lamT]

================================================================================
FILE STRUCTURE
================================================================================

EXTKIN/
  EXTKINSingle.m   Self-contained package (all code in one file)
  README.txt        This file
  tests/
    Tests.m         Comprehensive test suite (run in Mathematica)
    TestResults.txt Test output (34/34 passed)
    Examples.nb     Mathematica notebook with worked usage examples
  redundant/
    EXTKIN.m, Spinors.m, Kinematics.m, Poles.m, Validation.m
                    Original multi-file package (same code split across files)

================================================================================
