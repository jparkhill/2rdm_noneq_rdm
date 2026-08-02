(* Replay the MKNO notebook's input cells in order, under the local Engine.
   Purpose: turn the notebook into a re-runnable derivation oracle, and find
   out which cells are self-contained definitions vs. broken/interactive.
   Logs one line per cell to results/mkno_replay.log. *)
$HistoryLength = 0;
inputs = Get["/root/workspace/docs/MKNO_inputs.m"];
limit  = ToExpression@Environment["MKNO_LIMIT"];
If[!IntegerQ[limit], limit = 1700];
skip   = {1594, 1601};                    (* megabyte paste-cells *)
log    = OpenWrite["/root/workspace/results/mkno_replay.log"];
Do[
  Module[{idx = row[[1]], src = row[[2]], held, t, res, status},
   If[idx > limit || MemberQ[skip, idx], Continue[]];
   held = Quiet@Check[ToExpression[src], $Failed];
   {t, res} = Quiet@Check[
      AbsoluteTiming@TimeConstrained[ReleaseHold[held], 90, $TimedOut],
      {0., $Failed}];
   status = Which[res === $TimedOut, "TIMEOUT", res === $Failed, "ERROR", True, "ok"];
   WriteString[log, StringJoin[ToString[idx], "\t", status, "\t",
       ToString[NumberForm[t, {5, 2}]], "\t",
       StringTake[ToString[Head[res]], UpTo[24]], "\t",
       StringTake[StringReplace[ToString[res, InputForm], {"\n" -> " "}], UpTo[110]], "\n"]];
   ], {row, inputs}];
Close[log];
Print["replay complete through cell ", limit];
