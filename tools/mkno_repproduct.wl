(* Multiply ORBIT REPRESENTATIVES instead of the full term list.

   If the product rule commutes with the external index symmetries of the
   factors -- and it should, since those indices are only permuted, never
   contracted against each other -- then

       (antisymmetrized reps) * fac3  ==  antisymmetrize( reps * fac3 )

   so the engine never has to build the redundant images at all.  This script
   builds reps*fac3; tools/orbit_reduce.py then checks the expansion against
   the full stage_rv1v2.m built the ordinary way. *)
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
fac3 = NRta2[{"v", "x"}, {"w", "y"}];
reps = Get[res <> "stage_r1v1_reps.m"];
Print["representatives in: ", Length[reps]];
t0 = AbsoluteTime[];
prod = MemoryConstrained[TimeConstrained[
   AllStepsOfProductRule[reps, fac3], 2400, $TimedOut], 5*1024^3, $OutOfMem];
If[prod === $TimedOut || prod === $OutOfMem,
  Print["FAILED: ", prod]; Abort[]];
Print["[built ] rv1v2 from reps: ", Length[prod], " terms, ", lam[prod],
      " with lambda   (", Round[AbsoluteTime[] - t0], "s, peak ",
      Round[MaxMemoryUsed[]/1024^2], " MB)"];
Put[prod, res <> "stage_rv1v2_fromreps.m"];
Print["STAGE DONE"];
