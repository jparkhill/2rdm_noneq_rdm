(* Minimal reproducer: does the product rule stay exact when the residual
   operator reaches rank 5?

     T = {a^{p,r,t}_{q,s,u}}   (rank 3)
     C = {a^{v,x}_{w,y}}       (rank 2)
     T . C  ->  residual operators up to rank 5.

   Also traps any call of Signature on a non-list, printing the caller's stack,
   since exactly one such error appeared in the rv1v2 build. *)
$HistoryLength = 0; $RecursionLimit = 16384;
Off[ReplaceRepeated::rrlim]; Off[General::stop]; Off[Join::heads];
res = "/root/workspace/results/";
inputs = Get["/root/workspace/docs/MKNO_inputs.m"];
Do[Module[{i = row[[1]], s = row[[2]]},
   If[i <= 1600 && StringContainsQ[s, ":="],
     Quiet@Check[TimeConstrained[ReleaseHold[ToExpression[s]], 60], Null]]],
  {row, inputs}];

(* --- trap --- *)
Unprotect[Signature];
Signature[x_ /; ! ListQ[x]] := (
   Print["SIGTRAP arg=", HoldForm[x], "  head=", Head[x]];
   Print["SIGTRAP stack=", Take[Stack[_], UpTo[8]]];
   1);
Protect[Signature];

lam[l_] := Count[l, x_ /; ! FreeQ[x, \[Lambda]]];
SubstituteDeltas = False;
t3 = NRta3[{"p","r","t"}, {"q","s","u"}];
c2 = NRta2[{"v","x"}, {"w","y"}];
Print["NRta3 terms: ", Length[t3], "   NRta2 terms: ", Length[c2]];
prod = AllStepsOfProductRule[t3, c2];
Print["p32: ", Length[prod], " terms, ", lam[prod], " with lambda"];
Print["max residual rank: ", Max[Cases[prod, {ta, r_} :> r, Infinity]]];
Put[prod, res <> "stage_p32.m"];
Print["STAGE DONE"];
