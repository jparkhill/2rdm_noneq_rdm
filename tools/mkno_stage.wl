(* Build ONE stage of the full-cumulant 1-2-2 derivation, under hard memory and
   time guards, and checkpoint it to results/stage_<name>.m.

   Usage:  MKNO_STAGE=rv1v2 tools/wl -file tools/mkno_stage.wl

   Stages: r1v1 v1r1 rv1v2 v1rv2 v2v1r v2rv1 combine
   Each runs in its own kernel with the whole machine to itself, so one
   pathological stage cannot cost the others.  MemoryConstrained makes an
   over-large stage report its size instead of being OOM-killed silently.
*)
$HistoryLength = 0;
$RecursionLimit = 16384;
Off[ReplaceRepeated::rrlim]; Off[General::stop]; Off[Join::heads];

res  = "/root/workspace/results/";
want = Environment["MKNO_STAGE"];
(* Default 12 GB.  The container's cgroup hard-caps RAM at 8 GB (memory.max,
   read-only, remount denied) but allows up to 8 GB of swap, and the host has
   ~7 GB free, so the effective ceiling before an OOM kill is ~15 GB.  12 GB
   leaves headroom and lets a big stage spill to swap -- slow, but it finishes
   instead of dying.  /proc/meminfo reports the HOST's 32 GB and 32 CPUs; it is
   not namespaced, so ignore it and trust /sys/fs/cgroup/memory.max. *)
cap  = With[{e = ToExpression@Environment["MKNO_MEMCAP"]},
            If[IntegerQ[e], e, 12]] * 1024^3;          (* bytes *)
tcap = With[{e = ToExpression@Environment["MKNO_TIMECAP"]},
            If[IntegerQ[e], e, 1500]];                (* seconds *)

inputs = Get["/root/workspace/docs/MKNO_inputs.m"];
Do[Module[{i = row[[1]], s = row[[2]]},
   If[i <= 1600 && StringContainsQ[s, ":="],
     Quiet@Check[TimeConstrained[ReleaseHold[ToExpression[s]], 60], Null]]],
  {row, inputs}];
Print["definitions loaded; stage = ", want, "; memcap = ",
      Round[cap/1024^3], " GB; timecap = ", tcap, " s"];

lam[l_] := Count[l, x_ /; ! FreeQ[x, \[Lambda]]];
load[n_] := Module[{f = res <> "stage_" <> n <> ".m"},
   If[! FileExistsQ[f], Print["MISSING PREREQ: ", n]; Abort[]]; Get[f]];

build[n_, expr_] := Module[{f = res <> "stage_" <> n <> ".m", t0, v},
  If[FileExistsQ[f],
    v = Get[f];
    Print["[cached] ", n, ": ", Length[v], " terms, ", lam[v], " with lambda"];
    Return[v]];
  t0 = AbsoluteTime[];
  v = MemoryConstrained[TimeConstrained[expr, tcap, $TimedOut], cap, $OutOfMem];
  Which[
    v === $TimedOut, Print["TIMEOUT building ", n, " after ", tcap, "s"]; Abort[],
    v === $OutOfMem, Print["OUT OF MEMORY building ", n, " at ",
                           Round[cap/1024^3], " GB"]; Abort[]];
  Put[v, f];
  Print["[built ] ", n, ": ", Length[v], " terms, ", lam[v], " with lambda   (",
        Round[AbsoluteTime[] - t0], "s, peak ", Round[MaxMemoryUsed[]/1024^2], " MB)"];
  v];

SubstituteDeltas = False;
fac1 = NRta1[{"p"}, {"q"}];
fac2 = NRta2[{"r", "t"}, {"s", "u"}];
fac3 = NRta2[{"v", "x"}, {"w", "y"}];

Switch[want,
 "r1v1",  build["r1v1",  AllStepsOfProductRule[fac1, fac2]],
 "v1r1",  build["v1r1",  AllStepsOfProductRule[fac2, fac1]],
 "rv1v2", build["rv1v2", AllStepsOfProductRule[load["r1v1"], fac3]],
 "v1rv2", build["v1rv2", AllStepsOfProductRule[load["v1r1"], fac3]],
 "v2v1r", build["v2v1r", AllStepsOfProductRule[fac3, load["v1r1"]]],
 "v2rv1", build["v2rv1", AllStepsOfProductRule[fac3, load["r1v1"]]],
 "combine",
   Module[{a = load["rv1v2"], b = load["v2v1r"], c = load["v1rv2"], d = load["v2rv1"], all},
    all = Join[a, b, MultiplyScalar[c, -1], MultiplyScalar[d, -1]];
    Print["AllTerms (raw): ", Length[all], " terms, ", lam[all], " with lambda"];
    sAll = build["simplified_full",  AssociativeSimplify[all]];
    sNo  = build["simplified_nocum", AssociativeSimplify[all //. NoCumulants]];
    Print["DELTA (terms the old derivations deleted): ", Length[sAll] - Length[sNo]];
    Put[sAll, res <> "mkno_full_cumulant.m"]; Put[sNo, res <> "mkno_full_nocum.m"];
    census = Tally[Cases[sAll, t_List :> {Count[t, {{\[Lambda], _}, __}],
       Sort[Cases[t, {{\[Lambda], k_}, __} :> k]]}, {1}]];
    Print["lambda census {nLambda, ranks} -> count:"];
    Do[Print["   ", c], {c, SortBy[census, -Last[#] &]}];
    WriteString[res <> "mkno_full_cumulant_stats.txt",
      "full  = " <> ToString[Length[sAll]] <> " terms, " <> ToString[lam[sAll]] <>
      " carry lambda\nnocum = " <> ToString[Length[sNo]] <> " terms\ncensus: " <>
      ToString[census] <> "\n"]],
 _, Print["set MKNO_STAGE to one of: r1v1 v1r1 rv1v2 v1rv2 v2v1r v2rv1 combine"]];
Print["STAGE DONE"];
