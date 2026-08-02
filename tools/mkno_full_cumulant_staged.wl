(* Staged/checkpointed version of tools/mkno_full_cumulant.wl.

   Same derivation -- full 1-2-2 triple commutator, cumulants RETAINED -- but
   each of the four operator products is written to results/ as soon as it is
   built, and any stage already on disk is reloaded instead of recomputed.
   So an OOM, a timeout, or a manual stop costs one stage, not the run.

   Stage files:  results/stage_<name>.m
   Run repeatedly until it prints DONE.
*)
$HistoryLength = 0;
$RecursionLimit = 16384;
Off[ReplaceRepeated::rrlim]; Off[General::stop]; Off[Join::heads];

res = "/root/workspace/results/";
inputs = Get["/root/workspace/docs/MKNO_inputs.m"];
Do[Module[{i = row[[1]], s = row[[2]]},
   If[i <= 1600 && StringContainsQ[s, ":="],
     Quiet@Check[TimeConstrained[ReleaseHold[ToExpression[s]], 60], Null]]],
  {row, inputs}];
Print["definitions loaded"];

lamCount[l_] := Count[l, x_ /; ! FreeQ[x, \[Lambda]]];

stage[name_String, expr_] := Module[{f = res <> "stage_" <> name <> ".m", v},
  If[FileExistsQ[f],
    v = Get[f];
    Print["  [cached] ", name, ": ", Length[v], " terms, ", lamCount[v], " with lambda"],
    Module[{t0 = AbsoluteTime[]},
      v = expr;
      Put[v, f];
      Print["  [built ] ", name, ": ", Length[v], " terms, ", lamCount[v],
            " with lambda   (", Round[AbsoluteTime[] - t0], "s, ",
            Round[MemoryInUse[]/1024^2], " MB)"]]];
  v];

SubstituteDeltas = False;
fac1 = NRta1[{"p"}, {"q"}];
fac2 = NRta2[{"r", "t"}, {"s", "u"}];
fac3 = NRta2[{"v", "x"}, {"w", "y"}];

r1v1  = stage["r1v1",  AllStepsOfProductRule[fac1, fac2]];
v1r1  = stage["v1r1",  AllStepsOfProductRule[fac2, fac1]];
rv1v2 = stage["rv1v2", AllStepsOfProductRule[r1v1, fac3]];
v1rv2 = stage["v1rv2", AllStepsOfProductRule[v1r1, fac3]];
v2v1r = stage["v2v1r", AllStepsOfProductRule[fac3, v1r1]];
v2rv1 = stage["v2rv1", AllStepsOfProductRule[fac3, r1v1]];

AllTerms = Join[rv1v2, v2v1r, MultiplyScalar[v1rv2, -1], MultiplyScalar[v2rv1, -1]];
Print["AllTerms (raw): ", Length[AllTerms], " terms, ", lamCount[AllTerms], " with lambda"];

sAll = stage["simplified_full", AssociativeSimplify[AllTerms]];
sNo  = stage["simplified_nocum", AssociativeSimplify[AllTerms //. NoCumulants]];
Print["DELTA (terms the old derivations deleted): ", Length[sAll] - Length[sNo]];

Put[sAll, res <> "mkno_full_cumulant.m"];
Put[sNo,  res <> "mkno_full_nocum.m"];
census = Tally[Cases[sAll, t_List :> {Count[t, {{\[Lambda], _}, __}],
    Sort[Cases[t, {{\[Lambda], k_}, __} :> k]]}, {1}]];
Print["lambda census {nLambda, ranks} -> count:"];
Do[Print["   ", c], {c, SortBy[census, -Last[#] &]}];
WriteString[res <> "mkno_full_cumulant_stats.txt",
  "full  = " <> ToString[Length[sAll]] <> " terms, " <> ToString[lamCount[sAll]] <>
  " carry lambda\nnocum = " <> ToString[Length[sNo]] <> " terms\ncensus: " <>
  ToString[census] <> "\n"];
Print["DONE"];
