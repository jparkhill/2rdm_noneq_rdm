(* Build small product-rule test cases so the defect can be localised by rank. *)
$HistoryLength = 0; $RecursionLimit = 16384;
Off[ReplaceRepeated::rrlim]; Off[General::stop]; Off[Join::heads];
res = "/root/workspace/results/";
inputs = Get["/root/workspace/docs/MKNO_inputs.m"];
Do[Module[{i = row[[1]], s = row[[2]]},
   If[i <= 1600 && StringContainsQ[s, ":="],
     Quiet@Check[TimeConstrained[ReleaseHold[ToExpression[s]], 60], Null]]],
  {row, inputs}];
lam[l_] := Count[l, x_ /; ! FreeQ[x, \[Lambda]]];
SubstituteDeltas = False;
b = NRta2[{"r","t"},{"s","u"}];
c = NRta2[{"v","x"},{"w","y"}];
Print["--- NRta5 availability probe ---"];
probe = NRta5[{"a","b","c","d","e"},{"f","g","h","i","j"}];
Print["NRta5 head: ", Head[probe]];
Print["--- rank2 x rank2 ---"];
p22 = AllStepsOfProductRule[b, c];
Print["p22: ", Length[p22], " terms, ", lam[p22], " with lambda"];
Put[p22, res <> "stage_p22.m"];
Print["max residual rank: ", Max[Cases[p22, {ta, r_} :> r, Infinity]]];
Print["STAGE DONE"];
