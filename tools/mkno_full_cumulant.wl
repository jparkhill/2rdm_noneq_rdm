(* FULL-SIZE derivation of the 1-2-2 triple commutator WITH cumulants retained.

   This is the notebook's cell 1611/1614/1618/1619/1620/1621/1622 chain, with
   two deliberate changes:

     1. AllStepsOfProductRule1PV is undefined (it lived in an unsaved kernel
        session), so we use the general AllStepsOfProductRule everywhere -- the
        same substitution the notebook itself makes in the 2-1-1 chain at cells
        1630-1638.
     2. NoCumulants is NOT applied.  Every lambda term is kept.

   Output: results/mkno_full_cumulant.m  (the raw term list)
           results/mkno_full_cumulant_stats.txt
*)
$HistoryLength = 0;
$RecursionLimit = 16384;
Off[ReplaceRepeated::rrlim]; Off[General::stop]; Off[Join::heads];

inputs = Get["/root/workspace/docs/MKNO_inputs.m"];
nd = 0;
Do[Module[{i = row[[1]], s = row[[2]]},
   If[i <= 1600 && StringContainsQ[s, ":="],
     Quiet@Check[TimeConstrained[ReleaseHold[ToExpression[s]], 60], Null]; nd++]],
  {row, inputs}];
Print["definitions loaded: ", nd];

lamCount[l_] := Count[l, x_ /; ! FreeQ[x, \[Lambda]]];
report[name_, l_] := Print["  ", name, ": ", Length[l], " terms, ",
   lamCount[l], " with lambda"];

SubstituteDeltas = False;
fac1 = NRta1[{"p"}, {"q"}];
fac2 = NRta2[{"r", "t"}, {"s", "u"}];
fac3 = NRta2[{"v", "x"}, {"w", "y"}];
Print["factors: ", Length /@ {fac1, fac2, fac3}];

t0 = AbsoluteTime[];
r1v1 = AllStepsOfProductRule[fac1, fac2];   report["r*v1  ", r1v1];
v1r1 = AllStepsOfProductRule[fac2, fac1];   report["v1*r  ", v1r1];
Print["  (", Round[AbsoluteTime[] - t0], "s)"];

t0 = AbsoluteTime[];
rv1v2 = AllStepsOfProductRule[r1v1, fac3];  report["r v1 v2", rv1v2];
Print["  (", Round[AbsoluteTime[] - t0], "s)"];
t0 = AbsoluteTime[];
v1rv2 = AllStepsOfProductRule[v1r1, fac3];  report["v1 r v2", v1rv2];
v2v1r = AllStepsOfProductRule[fac3, v1r1];  report["v2 v1 r", v2v1r];
v2rv1 = AllStepsOfProductRule[fac3, r1v1];  report["v2 r v1", v2rv1];
Print["  (", Round[AbsoluteTime[] - t0], "s)"];

AllTerms = Join[rv1v2, v2v1r, MultiplyScalar[v1rv2, -1], MultiplyScalar[v2rv1, -1]];
report["AllTerms (raw commutator)", AllTerms];

t0 = AbsoluteTime[];
sAll = AssociativeSimplify[AllTerms];
Print["AssociativeSimplify -> ", Length[sAll], " terms, ",
      lamCount[sAll], " with lambda   (", Round[AbsoluteTime[] - t0], "s)"];

(* the same object with the notebook's NoCumulants applied, for a clean diff *)
sNo = AssociativeSimplify[AllTerms //. NoCumulants];
Print["with NoCumulants applied -> ", Length[sNo], " terms, ",
      lamCount[sNo], " with lambda"];
Print["DELTA (terms the old derivations deleted): ", Length[sAll] - Length[sNo]];

Put[sAll, "/root/workspace/results/mkno_full_cumulant.m"];
Put[sNo,  "/root/workspace/results/mkno_full_nocum.m"];

(* term-shape census: how many lambda factors, and of what rank *)
census = Tally[Cases[sAll, t_List :> {Count[t, {{\[Lambda], _}, __}],
    Sort[Cases[t, {{\[Lambda], k_}, __} :> k]]}, {1}]];
Print["lambda census {nLambda, ranks} -> count:"];
Do[Print["   ", c], {c, SortBy[census, -Last[#] &]}];

WriteString["/root/workspace/results/mkno_full_cumulant_stats.txt",
  "full  = " <> ToString[Length[sAll]] <> " terms, " <> ToString[lamCount[sAll]] <>
  " carry lambda\nnocum = " <> ToString[Length[sNo]] <> " terms\n" <>
  "census: " <> ToString[census] <> "\n"];
Print["DONE"];
