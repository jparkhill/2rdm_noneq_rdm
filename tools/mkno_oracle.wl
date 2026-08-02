(* Targeted oracle run: load ONLY the engine definitions from MKNO.nb, then
   re-derive the 1-2-2 triple commutator (notebook cells 1611-1624) from
   scratch and compare with the 96-term list scraped from the stored output.

   A straight top-to-bottom replay dies around cell 1549 on a RecursionLimit
   in unrelated scratch work, so we skip everything that is not a definition. *)
$HistoryLength = 0;
$RecursionLimit = 8192;
inputs = Get["/root/workspace/docs/MKNO_inputs.m"];
defQ[src_] := StringContainsQ[src, ":="];
nd = 0;
Do[Module[{idx = row[[1]], src = row[[2]]},
   If[idx <= 1600 && defQ[src],
     Quiet@Check[TimeConstrained[ReleaseHold[ToExpression[src]], 60], Null];
     nd++]], {row, inputs}];
Print["definition cells evaluated: ", nd];
Print["AllStepsOfProductRule1PV defined: ",
      ValueQ[AllStepsOfProductRule1PV] || DownValues[AllStepsOfProductRule1PV] =!= {}];

(* --- re-derive, exactly as cells 1611/1614/1618/1619/1620/1621/1622 do --- *)
SubstituteDeltas = False;
fac1 = NRta1[{"p"}, {"q"}];
fac2 = NRta2[{"r", "t"}, {"s", "u"}];
fac3 = NRta2[{"v", "x"}, {"w", "y"}];
Print["fac1/2/3 built: ", {Length[fac1], Length[fac2], Length[fac3]}];

t0 = AbsoluteTime[];
temp1 = AllStepsOfProductRule[fac1, fac2];
Print["temp1 terms: ", Length[temp1], "   (", Round[AbsoluteTime[] - t0], "s)"];

rv1v2 = AllStepsOfProductRule1PV[temp1, fac3];
v1rv2 = AllStepsOfProductRule1PV[AllStepsOfProductRule[fac2, fac1], fac3];
v2v1r = AllStepsOfProductRule1PV[fac3, AllStepsOfProductRule[fac2, fac1]];
v2rv1 = AllStepsOfProductRule1PV[fac3, AllStepsOfProductRule[fac1, fac2]];
Print["four products: ", Length /@ {rv1v2, v1rv2, v2v1r, v2rv1}];

AllTerms = Join[rv1v2, v2v1r, MultiplyScalar[v1rv2, -1], MultiplyScalar[v2rv1, -1]];
Print["AllTerms: ", Length[AllTerms]];
sAll = AssociativeSimplify[AllTerms];
Print["AssociativeSimplify[AllTerms] length = ", Length[sAll],
      "    (notebook cell 1624 reported 96)"];
Export["/root/workspace/results/mkno_oracle_1624.m", sAll];
Print["wrote results/mkno_oracle_1624.m"];
