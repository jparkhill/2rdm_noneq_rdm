(* Extend the cell-62 dispatchers to rank 5.

   ta5/Na5/NRta5 are defined in cell 125 and cell 146 already extends
   InversionOfVacNormal and InversionOfGamma to rank 5, but only cell 62
   dispatches tilde substitution and it stops at rank 4. Past that it called
   Print, which returns Null, and DistributeOperatorProduct multiplied that
   Null into the product -- so a 1-2-2 product, which reaches rank 5 the
   moment cumulants are retained, was silently wrong.

   Throw replaces Print so the next missing rank fails loudly instead. *)

Clear[TildeSubstitution, ASubstitution, NRTildeSubstitution];

TildeSubstitution[ui_, li_] := Switch[Length[ui],
  1, ta1[ui, li], 2, ta2[ui, li], 3, ta3[ui, li], 4, ta4[ui, li], 5, ta5[ui, li],
  _, Throw[{"TildeSubstitution: unimplemented rank", Length[ui]}, "MKNORank"]];

ASubstitution[ui_, li_] := Switch[Length[ui],
  1, Na1[ui, li], 2, Na2[ui, li], 3, Na3[ui, li], 4, Na4[ui, li], 5, Na5[ui, li],
  _, Throw[{"ASubstitution: unimplemented rank", Length[ui]}, "MKNORank"]];

NRTildeSubstitution[ui_, li_] := Switch[Length[ui],
  1, NRta1[ui, li], 2, NRta2[ui, li], 3, NRta3[ui, li], 4, NRta4[ui, li],
  5, NRta5[ui, li],
  _, Throw[{"NRTildeSubstitution: unimplemented rank", Length[ui]}, "MKNORank"]];
