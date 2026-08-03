$HistoryLength = 0; $RecursionLimit = 16384;
Off[ReplaceRepeated::rrlim]; Off[General::stop]; Off[Join::heads];
inputs = Get["/root/workspace/docs/MKNO_inputs.m"];
Do[Module[{i = row[[1]], s = row[[2]]},
   If[i <= 1600 && StringContainsQ[s, ":="],
     Quiet@Check[TimeConstrained[ReleaseHold[ToExpression[s]], 60], Null]]],
  {row, inputs}];
Print["definitions loaded"];
r4 = TildeSubstitution[{p1,p2,p3,p4},{q1,q2,q3,q4}];
r5 = Quiet@TildeSubstitution[{p1,p2,p3,p4,p5},{q1,q2,q3,q4,q5}];
Print["rank4 head=", Head[r4], "  length=", If[ListQ[r4], Length[r4], "n/a"]];
Print["rank5 head=", Head[r5], "  isNull=", r5 === Null];
n5 = Quiet@NRTildeSubstitution[{p1,p2,p3,p4,p5},{q1,q2,q3,q4,q5}];
Print["NR rank5 isNull=", n5 === Null];
Print["ta5 defined=", ValueQ[ta5[{p1,p2,p3,p4,p5},{q1,q2,q3,q4,q5}]] || Head[ta5[{p1,p2,p3,p4,p5},{q1,q2,q3,q4,q5}]] === List];
