

## [1] Implementation Of Generalized Normal Ordering based on ' brute force' permutation of density matrices.
**Input [3]**
```
CellGroupData[{ Cell["Explanation. ", "Subsection"
```


### [4] Explanation.
> <p and q are actually exceptions. They are reserved for external indices. 
> B,M, and L are reserved for RI integral indices. >
**Input [7]**
```
[""], "Input"
```


### [9] Step 1, Automatic definition of an MK tilde - a operator
**Input [12]**
```
[ RowBox[{"an", "[", "1", "]"}]], "Input"
```
**Input [15]**
```
[ RowBox[{"ta1", "[", {b},{c}, "]"}]], "Input"
```
**Input [19]**
```
[ RowBox[{"an", "[", "2", "]"}]], "Input"
```
**Input [22]**
```
[{ RowBox[{"ta2", "[", {\<b\>\<d\>},{\<c\>\<e\>}, "]"}], " ", RowBox[{"RecognizableForm", "[", "%", "]"}]}], "Input"
```
> <This is MK' s correct 126 (Mol. Phys.), Signs correct, everything.>
**Input [28]**
```
[ RowBox[{"an", "[", "3", "]"}]], "Input"
```
**Input [31]**
```
[{ RowBox[{ RowBox[{"ta3", "[", {\<b\>\<d\>\<f\>},{\<c\>\<e\>\<g\>}, "]"}], ";"}], " ", RowBox[{"RecognizableForm", "[", "%", "]"}], " ", RowBox[{"Length", "[", "%", "]"}]}], "Input"
```
**Input [37]**
```
[ RowBox[{"an", "[", "4", "]"}]], "Input"
```
**Input [40]**
```
[{ RowBox[{ RowBox[{"tmp", "=", RowBox[{"ta4", "[", {\<b\>\<d\>\<f\>\<h\>},{\<c\>\<e\>\<g\>\<i\>}, "]"}]}], ";"}], " ", RowBox[{"Length", "[", "tmp", "]"}], " ", RowBox[{"RecognizableForm", "[", "tmp", "]"}], " ", RowBox[{"Length", "[", RowBox[{"RecognizableForm", "[", "tmp", "]"}], "]"}]}], "Input"
```
> <It' s also useful to have versions of these which are completely free of any reference to the tilde operators so that the a - components of a - tilde can be completely normal ordered. Otherwise you can lose a degree of simplification in later steps, especially for high rank expressions. Also must be careful to only sort factors when they are scalars, and remember mathematica assumes multiplication is associative! I do this after defining multiplication... >


### [47] Expressions for Normal operators in terms of Mk - Normal Operators
> <These are basically the MK normal order definitions, but the signs are flipped. They are important because they are the only terms which MK give the signs for every contribution, explicitly, so they are an important check of the work so far. Note that the normal operators are at most LINEAR in the MK normal operators.  >
**Input [49]**
```
[{ RowBox[{ RowBox[{ RowBox[{"Na1", "[", ui_,li_, "]"}], ":=", RowBox[{"Join", "[", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", {ta1},ui,li, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", "ui", ",", "li"}], "}"}]}], "}"}]}], "}"}], "]"}]}], ";"}], " ", RowBox[{ RowBox[{"Na2", "[", uis_,lis_, "]"}], ":=", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", {ta2},uis,lis, "}"}]}], "}"}], "}"}], ",", RowBox[{"Grassmann", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", RowBox[{"Take", "[", uis,{1}, "]"}], ",", RowBox[{"Take", "[", lis,{1}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {ta1}, ",", RowBox[{"Take", "[", uis,{2}, "]"}], ",", RowBox[{"Take", "[", lis,{2}, "]"}]}], "}"}], ",", "4"}], "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "2"}], "}"}], ",", "uis", ",", "lis"}], "}"}]}], "}"}], "}"}]}], "]"}]}], " ", RowBox[{ RowBox[{"Na3", "[", uis_,lis_, "]"}], ":=", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", {ta3},uis,lis, "}"}]}], "}"}], "}"}], ",", " ", RowBox[{"Grassmann", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", RowBox[{"Take", "[", uis,{1}, "]"}], ",", RowBox[{"Take", "[", lis,{1}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {ta2}, ",", RowBox[{"Take", "[", uis,{23}, "]"}], ",", RowBox[{"Take", "[", lis,{23}, "]"}]}], "}"}], ",", "9"}], "]"}], ",", " ", RowBox[{"Grassmann", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "2"}], "}"}], ",", RowBox[{"Take", "[", uis,{12}, "]"}], ",", RowBox[{"Take", "[", lis,{12}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {ta1}, ",", RowBox[{"Take", "[", uis,{3}, "]"}], ",", RowBox[{"Take", "[", lis,{3}, "]"}]}], "}"}], ",", "9"}], "]"}], ",", " ", RowBox[{"{", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "3"}], "}"}], ",", "uis", ",", "lis"}], "}"}]}], "}"}], "}"}]}], "]"}]}], " ", RowBox[{ RowBox[{ RowBox[{"Na4", "[", uis_,lis_, "]"}], ":=", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", {ta4},uis,lis, "}"}]}], "}"}], "}"}], ",", " ", RowBox[{"Grassmann", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", RowBox[{"Take", "[", uis,{1}, "]"}], ",", RowBox[{"Take", "[", lis,{1}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {ta3}, ",", RowBox[{"Take", "[", uis,{24}, "]"}], ",", RowBox[{"Take", "[", lis,{24}, "]"}]}], "}"}], ",", "16"}], "]"}], ",", RowBox[{"Grassmann", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "2"}], "}"}], ",", RowBox[{"Take", "[", uis,{12}, "]"}], ",", RowBox[{"Take", "[", lis,{12}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {ta2}, ",", RowBox[{"Take", "[", uis,{34}, "]"}], ",", RowBox[{"Take", "[", lis,{34}, "]"}]}], "}"}], ",", "36"}], "]"}], ",", " ", RowBox[{"Grassmann", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "3"}], "}"}], ",", RowBox[{"Take", "[", uis,{13}, "]"}], ",", RowBox[{"Take", "[", lis,{13}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {ta1}, ",", RowBox[{"Take", "[", uis,{4}, "]"}], ",", RowBox[{"Take", "[", lis,{4}, "]"}]}], "}"}], ",", "16"}], "]"}], ",", " ", RowBox[{"{", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "4"}], "}"}], ",", "uis", ",", "lis"}], "}"}]}], "}"}], "}"}]}], "]"}]}], ";"}]}], "Input"
```
**Input [51]**
```
[{ RowBox[{"an", "[", "1", "]"}], " ", RowBox[{"an", "[", "2", "]"}], " ", RowBox[{"an", "[", "3", "]"}], " ", RowBox[{"an", "[", "4", "]"}]}], "Input"
```
**Input [57]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{"Na1", "[", {b},{c}, "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"Na2", "[", {\<b\>\<d\>},{\<c\>\<e\>}, "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"Na3", "[", {\<b\>\<d\>\<f\>},{\<c\>\<e\>\<g\>}, "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"Na4", "[", {\<b\>\<d\>\<f\>\<h\>},{\<c\>\<e\>\<g\>\<i\>}, "]"}], "]"}]}], "Input"
```
> The above is PRECISELY 27 b.
> <So let' s say step 1 is complete and we can now express the MK normal ordered operators in terms of the lower. We are now ready to attempt to derive Expressions for products of Normal ordered operators. >


### [66] Ordinary Normal ordering.
> <These are presicely MK 9 a and 9 b  up to signed permutations of the indices. We need these contractions to perform the products of ordinary Normal Operators that occur in the products of the MK normal operators. >
**Input [74]**
```
CellGroupData[{ Cell["Operator Mulitplication. ", "Subsection"
```


### [75] Operator Mulitplication.
**Input [76]**
```
[{ RowBox[{ RowBox[{ RowBox[{"TermProduct", "[", t1_,t2_, "]"}], ":=", RowBox[{"Module", "[", RowBox[{ {sign1sign2}, ",", RowBox[{"(", " ", RowBox[{ RowBox[{"sign1", "=", RowBox[{"t1", "[", RowBox[{"[", 1,2, "]"}], "]"}]}], ";", " ", RowBox[{"sign2", "=", RowBox[{"t2", "[", RowBox[{"[", 1,2, "]"}], "]"}]}], ";", " ", RowBox[{"Return", "[", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", {\<Sign\>sign1*sign2}, "}"}], ",", RowBox[{"Drop", "[", t1,1, "]"}], ",", RowBox[{"Drop", "[", t2,1, "]"}]}], "]"}], "]"}]}], " ", ")"}]}], "]"}]}], ";"}], " ", RowBox[{ RowBox[{ RowBox[{"DistributeOperatorProduct", "[", o1_,o2_, "]"}], ":=", RowBox[{"(", " ", RowBox[{"Return", "[", " ", RowBox[{ RowBox[{"Reap", "[", " ", RowBox[{"Do", "[", " ", RowBox[{ RowBox[{ RowBox[{"Sow", "[", RowBox[{"TermProduct", "[", RowBox[{ RowBox[{"o1", "[", RowBox[{"[", "i", "]"}], "]"}], ",", RowBox[{"o2", "[", RowBox[{"[", "j", "]"}], "]"}]}], "]"}], "]"}], ";"}], " ", ",", RowBox[{"{", RowBox[{"i", ",", "1", ",", RowBox[{"Length", "[", "o1", "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{"j", ",", "1", ",", RowBox[{"Length", "[", "o2", "]"}]}], "}"}]}], "]"}], " ", "]"}], "[", RowBox[{"[", 2,1, "]"}], "]"}], " ", "]"}], " ", ")"}]}], ";"}], " ", RowBox[{ RowBox[{ RowBox[{"TermMultiplyScalar", "[", t1_,sc_, "]"}], ":=", RowBox[{"Module", "[", RowBox[{ {sign1}, ",", RowBox[{"(", " ", RowBox[{ RowBox[{"sign1", "=", RowBox[{"t1", "[", RowBox[{"[", 1,2, "]"}], "]"}]}], ";", " ", RowBox[{"Return", "[", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", {\<Sign\>sign1*sc}, "}"}], ",", RowBox[{"Drop", "[", t1,1, "]"}]}], "]"}], "]"}]}], " ", ")"}]}], "]"}]}], ";"}], " ", RowBox[{ RowBox[{ RowBox[{"MultiplyScalar", "[", o1_,sc_, "]"}], ":=", RowBox[{"(", " ", RowBox[{"Return", "[", " ", RowBox[{ RowBox[{"Reap", "[", " ", RowBox[{"Do", "[", " ", RowBox[{ RowBox[{ RowBox[{"Sow", "[", RowBox[{"TermMultiplyScalar", "[", RowBox[{ RowBox[{"o1", "[", RowBox[{"[", "i", "]"}], "]"}], ",", "sc"}], "]"}], "]"}], ";"}], " ", ",", RowBox[{"{", RowBox[{"i", ",", "1", ",", RowBox[{"Length", "[", "o1", "]"}]}], "}"}]}], "]"}], " ", "]"}], "[", RowBox[{"[", 2,1, "]"}], "]"}], " ", "]"}], " ", ")"}]}], ";"}], " ", RowBox[{ RowBox[{ RowBox[{"TermUnifyScalar", "[", "t1_", "]"}], ":=", RowBox[{"Module", "[", RowBox[{ {sign1}, ",", RowBox[{"(", " ", RowBox[{ RowBox[{"sign1", "=", RowBox[{"t1", "[", RowBox[{"[", 1,2, "]"}], "]"}]}], ";", " ", RowBox[{"Return", "[", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\"\<Sign\>\"", ",", RowBox[{"Sign", "[", "sign1", "]"}]}], "}"}], "}"}], ",", RowBox[{"Drop", "[", t1,1, "]"}]}], "]"}], "]"}]}], " ", ")"}]}], "]"}]}], ";"}], " ", RowBox[{ RowBox[{ RowBox[{"UnifyScalars", "[", "o1_", "]"}], ":=", RowBox[{"(", " ", RowBox[{"Return", "[", " ", RowBox[{ RowBox[{"Reap", "[", " ", RowBox[{"Do", "[", " ", RowBox[{ RowBox[{ RowBox[{"Sow", "[", RowBox[{"TermUnifyScalar", "[", RowBox[{"o1", "[", RowBox[{"[", "i", "]"}], "]"}], "]"}], "]"}], ";"}], " ", ",", RowBox[{"{", RowBox[{"i", ",", "1", ",", RowBox[{"Length", "[", "o1", "]"}]}], "}"}]}], "]"}], " ", "]"}], "[", RowBox[{"[", 2,1, "]"}], "]"}], " ", "]"}], " ", ")"}]}], ";"}]}], "Input"
```


### [78] NonRecursive definition of MK operators in terms of a operators.
**Input [82]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{"NRta3", "[", {\<b\>\<d\>\<f\>},{\<c\>\<e\>\<g\>}, "]"}], "]"}], " ", RowBox[{ RowBox[{"RecognizableForm", "[", RowBox[{"NRta4", "[", {\<b\>\<d\>\<f\>\<h\>},{\<c\>\<e\>\<g\>\<i\>}, "]"}], "]"}], ";"}]}], "Input"
```


### [85] One - One Products.
> <The above is the first line of MK-JCP 37. 
> We' ll leave this here and take a diversion to calculate the a' s and gammas now after checking the 2-2 product looks right.  >
**Input [96]**
```
[{ RowBox[{"an", "[", "2", "]"}], " ", RowBox[{"an", "[", 2,5, "]"}]}], "Input"
```
**Input [100]**
```
[{ RowBox[{"fac1", "=", RowBox[{"ta2", "[", {\<b\>\<d\>},{\<c\>\<e\>}, "]"}]}], " ", RowBox[{"fac2", "=", RowBox[{"ta2", "[", {\<f\>\<h\>},{\<g\>\<i\>}, "]"}]}], " ", RowBox[{ RowBox[{"prod", "=", RowBox[{"DistributeOperatorProduct", "[", fac1,fac2, "]"}]}], ";"}], " ", RowBox[{ RowBox[{"nrmprod", "=", RowBox[{"ContractVacuumNormalProducts", "[", "prod", "]"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", "nrmprod", "]"}], " ", RowBox[{"Clear", "[", fac1,fac2,prod,nrmprod, "]"}]}], "Input",\ CellChangeTimes->{{3.625939461268899*^9, 3.625939502585144*^9}, { 3.625939556460792*^9, 3.625939559812894*^9}, {3.625943768973413*^9, 3.625943769142695*^9}}
```
> <Next Contracted a' s must be replaced with the equivalent tilde a expression determined by the inversion worked out in the previous section. The density matrices also need to be replaced by their cumulant decomposition. >


### [106] Cumulant Decomposition of Gamma up to rank 4.
> <This Requires Products of Grassman Products which are a whole new clusterfuck. I[CloseCurlyQuote]ll do them up to four for now. 
> NOTE: These routines assume the arguement tensors have normal antisymmetry. If you ran them on delta tensors for example the results would be gibberish. >
**Input [111]**
```
[{ RowBox[{ RowBox[{"decompgamma1", "=", RowBox[{"Gd1", "[", {b},{c}, "]"}]}], ";"}], " ", RowBox[{ RowBox[{"decompgamma2", "=", RowBox[{"Gd2", "[", {\<b\>\<d\>},{\<c\>\<e\>}, "]"}]}], ";"}], " ", RowBox[{ RowBox[{"decompgamma3", "=", RowBox[{"Gd3", "[", {\<b\>\<d\>\<f\>},{\<c\>\<e\>\<g\>}, "]"}]}], ";"}], " ", RowBox[{ RowBox[{ RowBox[{"decompgamma4", "=", RowBox[{"Gd4", "[", {\<b\>\<d\>\<f\>\<h\>},{\<c\>\<e\>\<g\>\<i\>}, "]"}]}], ";"}], " "}], " ", RowBox[{"RecognizableForm", "[", "decompgamma1", "]"}], " ", RowBox[{"RecognizableForm", "[", "decompgamma2", "]"}], " ", RowBox[{"RecognizableForm", "[", "decompgamma3", "]"}], " ", RowBox[{"RecognizableForm", "[", "decompgamma4", "]"}], " ", RowBox[{"Length", "[", RowBox[{"RecognizableForm", "[", "decompgamma1", "]"}], "]"}], " ", RowBox[{"Length", "[", RowBox[{"RecognizableForm", "[", "decompgamma2", "]"}], "]"}], " ", RowBox[{"Length", "[", RowBox[{"RecognizableForm", "[", "decompgamma3", "]"}], "]"}], " ", RowBox[{"Length", "[", RowBox[{"RecognizableForm", "[", "decompgamma4", "]"}], "]"}]}], "Input"
```
**Input [121]**
```
[ 1+4+6+3+1+12+6+24+12+6+24+12+3+8+3+6], "Input"
```
> <That[CloseCurlyQuote]s exactly the same number of terms for [Gamma]4 counted by MK. We are now ready to complete the normal ordering. :) >


### [125] Finally Finishing The Product Rule (Testing both 1-1 and 1-2 products)
> <This is exactly MK' s Expression 144 from the Mol. Phys. Paper. 
> Which is also Eq. 39 of the JCP paper. We will also check the 1-2 product which should match 42a): >
> <This is directly comparable with MK' s 42 a (JCP). Although you could check me. >


### [140] Associative Simplification by Mapping into and out of ordinary expressions
**Output [143]**
```
CellGroupData[{ Cell["Evaluation Of 1-particle Triple Commutator", "Subsection"
```


### [144] Evaluation Of 1-particle Triple Commutator
> <If Fac1 or Fac2 contain Tilde Normal operators they need to be exchanged for a operators. >
**Input [147]**
```
[{ RowBox[{ RowBox[{"fac1", "=", RowBox[{"NRta1", "[", {"\<p\>"},{"\<q\>"}, "]"}]}], ";"}], " ", RowBox[{ RowBox[{"fac2", "=", RowBox[{"NRta1", "[", {"\<r\>"},{"\<s\>"}, "]"}]}], ";"}], " ", RowBox[{ RowBox[{"fac3", "=", RowBox[{"NRta1", "[", {"\<t\>"},{"\<u\>"}, "]"}]}], ";"}]}], "Input"
```
> <Hm ... Derive the Blocking Factor Type Of expression to compare with the two particle case. >
**Input [160]**
```
[ RowBox[{"TCLSimp", "=", RowBox[{"ListForm", "[", RowBox[{"Expand", "[", RowBox[{ RowBox[{ RowBox[{"(", d^{{"r"}}_{{"q"}}, ")"}], " ", RowBox[{"(", RowBox[{ RowBox[{ g^{{"p"}}_{{"u"}}, " ", eta^{{"t"}}_{{"s"}}}], "-", RowBox[{ g^{{"t"}}_{{"s"}}, " ", eta^{{"p"}}_{{"u"}}}]}], ")"}]}], "-", RowBox[{ RowBox[{"(", d^{{"p"}}_{{"s"}}, ")"}], " ", RowBox[{"(", RowBox[{ RowBox[{ g^{{"r"}}_{{"u"}}, " ", eta^{{"t"}}_{{"q"}}}], "-", RowBox[{ g^{{"t"}}_{{"q"}}, " ", eta^{{"r"}}_{{"u"}}}]}], ")"}]}]}], "]"}], "]"}]}]], "Input"
```
**Output [162]**
```
CellGroupData[{ Cell[[{ RowBox[{ RowBox[{ RowBox[{"STwo", "[", p_, q_, r_, s_, "]"}], " ", ":=", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {r}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {s}}], "}"}]}], "}"}], ",", " ", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {r}}], "}"}]}], "}"}], ",", " ", RowBox[{"{", RowBox[{ {\<Sign\>-1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {r}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {r}}], "}"}]}], "}"}]}], "}"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"TCLSimp2", "=", RowBox[{"MultiplyTermsByTensor", "[", RowBox[{ RowBox[{"MultiplyTermsByTensor", "[", RowBox[{"TCLSimp", ",", RowBox[{"{", RowBox[{"{", {\<VTs\>2},{"\<r\>"},{"\<s\>"}, "}"}], "}"}]}], "]"}], ",", RowBox[{"{", RowBox[{"{", {\<VTt\>2},{"\<t\>"},{"\<u\>"}, "}"}], "}"}]}], "]"}]}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"TCLSimp3", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"DistributeOperatorProduct", "[", RowBox[{"TCLSimp2", ",", RowBox[{"STwo", "[", "\<r\>","\<s\>","\<t\>","\<u\>", "]"}]}], "]"}]}], "]"}]}], "]"}]}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"TCLSimp4", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", DiagonalGamma1,TCLSimp3, "]"}]}], "]"}]}], "]"}]}], "]"}], " "}], "Input"
```
**Output [165]**
```
[ RowBox[{"Total", "[", RowBox[{"MultiplyTermsByTensor", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"\<\"Sign\"\>", " ", TemplateBox[{{VTs2},RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], ",", TemplateBox[{{VTs2},RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"g", " ", TemplateBox[{{VTs2},RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], ",", TemplateBox[{{VTs2},RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"p\"\>", " ", TemplateBox[{{VTs2},RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"u\"\>", " ", TemplateBox[{{VTs2},RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"d", " ", TemplateBox[{{VTs2},RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], ",", TemplateBox[{{VTs2},RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"r\"\>", " ", TemplateBox[{{VTs2},RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"q\"\>", " ", TemplateBox[{{VTs2},RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"eta", " ", TemplateBox[{{VTs2},RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], ",", TemplateBox[{{VTs2},RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"t\"\>", " ", TemplateBox[{{VTs2},RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"s\"\>", " ", TemplateBox[{{VTs2},RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", {VTs2},{"r"},{"s"}, "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], "}"}]}], "}"}]}], "}"}], ",", {\<"t"\>}}], "]"}], "]"}]], "Output"
```
**Output [170]**
```
[ RowBox[{ RowBox[{ RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", RowBox[{ {-\<"Sign"\>-1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"-", "g"}], ",", -1}], "}"}], ",", {-\<"p"\>}, ",", {-\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"-", "d"}], ",", -1}], "}"}], ",", {-\<"r"\>}, ",", {-\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"-", "eta"}], ",", -1}], "}"}], ",", {-\<"t"\>}, ",", {-\<"s"\>}}], "}"}]}], "}"}]}], "}"}], "}"}], ",", RowBox[{"MultiplyTermsByTensor", "[", RowBox[{"{", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>},{\<"s"\>}, "}"}], "}"}], "]"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"t"\>}}], "}"}]}], "}"}]}], "]"}], "]"}], "]"}], "[[", 1,2, "]]"}], "+", RowBox[{ RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}]}], "}"}], "}"}], ",", RowBox[{"MultiplyTermsByTensor", "[", RowBox[{"{", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>},{\<"s"\>}, "}"}], "}"}], "]"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"u"\>}}], "}"}]}], "}"}]}], "]"}], "]"}], "]"}], "[[", 1,2, "]]"}], "+", RowBox[{ RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}]}], "}"}], "}"}], ",", RowBox[{"MultiplyTermsByTensor", "[", RowBox[{"{", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>},{\<"s"\>}, "}"}], "}"}], "]"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"t"\>}}], "}"}]}], "}"}]}], "]"}], "]"}], "]"}], "[[", 1,2, "]]"}], "+", RowBox[{ RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"{", RowBox[{ RowBox[{"{", \<"Sign"\>,{\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"u"\>}}], "}"}]}], "}"}], "]"}], "]"}], "[[", 1,2, "]]"}], "+", RowBox[{ RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"{", RowBox[{ RowBox[{"{", \<"Sign"\>,{\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"t"\>}}], "}"}]}], "}"}], "]"}], "]"}], "[[", 1,2, "]]"}], "+", RowBox[{ RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"{", RowBox[{ RowBox[{"{", \<"Sign"\>,{-\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"t"\>}}], "}"}]}], "}"}], "]"}], "]"}], "[[", 1,2, "]]"}]}]], "Output"
```
**Output [171]**
```
[ RowBox[{ RowBox[{ RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"DiagonalGamma1", "[", RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", RowBox[{ {-\<"Sign"\>-1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"-", "g"}], ",", -1}], "}"}], ",", {-\<"p"\>}, ",", {-\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"-", "d"}], ",", -1}], "}"}], ",", {-\<"r"\>}, ",", {-\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"-", "eta"}], ",", -1}], "}"}], ",", {-\<"t"\>}, ",", {-\<"s"\>}}], "}"}]}], "}"}]}], "}"}], "}"}], ",", RowBox[{"MultiplyTermsByTensor", "[", RowBox[{"{", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>},{\<"s"\>}, "}"}], "}"}], "]"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"t"\>}}], "}"}]}], "}"}]}], "]"}], "]"}], "]"}], "]"}], "]"}], "]"}], "[[", 1,2, "]]"}], "+", RowBox[{ RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"DiagonalGamma1", "[", RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}]}], "}"}], "}"}], ",", RowBox[{"MultiplyTermsByTensor", "[", RowBox[{"{", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>},{\<"s"\>}, "}"}], "}"}], "]"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"u"\>}}], "}"}]}], "}"}]}], "]"}], "]"}], "]"}], "]"}], "]"}], "]"}], "[[", 1,2, "]]"}], "+", RowBox[{ RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"DiagonalGamma1", "[", RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}]}], "}"}], "}"}], ",", RowBox[{"MultiplyTermsByTensor", "[", RowBox[{"{", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>},{\<"s"\>}, "}"}], "}"}], "]"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"t"\>}}], "}"}]}], "}"}]}], "]"}], "]"}], "]"}], "]"}], "]"}], "]"}], "[[", 1,2, "]]"}], "+", RowBox[{ RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"DiagonalGamma1", "[", RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"{", RowBox[{ RowBox[{"{", \<"Sign"\>,{\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"u"\>}}], "}"}]}], "}"}], "]"}], "]"}], "]"}], "]"}], "]"}], "[[", 1,2, "]]"}], "+", RowBox[{ RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"DiagonalGamma1", "[", RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"{", RowBox[{ RowBox[{"{", \<"Sign"\>,{\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"t"\>}}], "}"}]}], "}"}], "]"}], "]"}], "]"}], "]"}], "]"}], "[[", 1,2, "]]"}], "+", RowBox[{ RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"DiagonalGamma1", "[", RowBox[{"CollapseMultiDeltas", "[", RowBox[{"CollapseDeltas", "[", RowBox[{"{", RowBox[{ RowBox[{"{", \<"Sign"\>,{-\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"t"\>}}], "}"}]}], "}"}], "]"}], "]"}], "]"}], "]"}], "]"}], "[[", 1,2, "]]"}]}]], "Output"
```


### [173] Accelerated Versions To treat more than 4 Particles
**Input [176]**
```
[{ RowBox[{ RowBox[{ RowBox[{"ta5", "[", uis_,lis_, "]"}], ":=", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", {a5},uis,lis, "}"}]}], "}"}], "}"}], ",", " ", RowBox[{"AssociativeSimplify", "[", RowBox[{"Grassmann", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", RowBox[{"Take", "[", uis,{1}, "]"}], ",", RowBox[{"Take", "[", lis,{1}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {ta4}, ",", RowBox[{"Take", "[", uis,{25}, "]"}], ",", RowBox[{"Take", "[", lis,{25}, "]"}]}], "}"}], ",", -25}], "]"}], "]"}], ",", " ", RowBox[{"AssociativeSimplify", "[", RowBox[{"Grassmann", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "2"}], "}"}], ",", RowBox[{"Take", "[", RowBox[{"uis", ",", RowBox[{"Range", "[", "2", "]"}]}], "]"}], ",", RowBox[{"Take", "[", RowBox[{"lis", ",", RowBox[{"Range", "[", "2", "]"}]}], "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {ta3}, ",", RowBox[{"Take", "[", uis,{35}, "]"}], ",", RowBox[{"Take", "[", lis,{35}, "]"}]}], "}"}], ",", -100}], "]"}], "]"}], ",", " ", RowBox[{"AssociativeSimplify", "[", RowBox[{"Grassmann", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "3"}], "}"}], ",", RowBox[{"Take", "[", uis,{13}, "]"}], ",", RowBox[{"Take", "[", lis,{13}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {ta2}, ",", RowBox[{"Take", "[", uis,{45}, "]"}], ",", RowBox[{"Take", "[", lis,{45}, "]"}]}], "}"}], ",", -100}], "]"}], "]"}], ",", " ", RowBox[{"AssociativeSimplify", "[", RowBox[{"Grassmann", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "4"}], "}"}], ",", RowBox[{"Take", "[", uis,{14}, "]"}], ",", RowBox[{"Take", "[", lis,{14}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {ta1}, ",", RowBox[{"Take", "[", uis,{5}, "]"}], ",", RowBox[{"Take", "[", lis,{5}, "]"}]}], "}"}], ",", -25}], "]"}], "]"}], ",", " ", RowBox[{"{", RowBox[{"{", RowBox[{ {\<Sign\>-1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "5"}], "}"}], ",", "uis", ",", "lis"}], "}"}]}], "}"}], "}"}]}], "]"}]}], ";"}], " ", RowBox[{ RowBox[{ RowBox[{"Na5", "[", uis_,lis_, "]"}], ":=", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", {ta5},uis,lis, "}"}]}], "}"}], "}"}], ",", " ", RowBox[{"AssociativeSimplify", "[", RowBox[{"Grassmann", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", RowBox[{"Take", "[", uis,{1}, "]"}], ",", RowBox[{"Take", "[", lis,{1}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {ta4}, ",", RowBox[{"Take", "[", uis,{25}, "]"}], ",", RowBox[{"Take", "[", lis,{25}, "]"}]}], "}"}], ",", "25"}], "]"}], "]"}], ",", " ", RowBox[{"AssociativeSimplify", "[", RowBox[{"Grassmann", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "2"}], "}"}], ",", RowBox[{"Take", "[", RowBox[{"uis", ",", RowBox[{"Range", "[", "2", "]"}]}], "]"}], ",", RowBox[{"Take", "[", RowBox[{"lis", ",", RowBox[{"Range", "[", "2", "]"}]}], "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {ta3}, ",", RowBox[{"Take", "[", uis,{35}, "]"}], ",", RowBox[{"Take", "[", lis,{35}, "]"}]}], "}"}], ",", "100"}], "]"}], "]"}], ",", " ", RowBox[{"AssociativeSimplify", "[", RowBox[{"Grassmann", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "3"}], "}"}], ",", RowBox[{"Take", "[", uis,{13}, "]"}], ",", RowBox[{"Take", "[", lis,{13}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {ta2}, ",", RowBox[{"Take", "[", uis,{45}, "]"}], ",", RowBox[{"Take", "[", lis,{45}, "]"}]}], "}"}], ",", "100"}], "]"}], "]"}], ",", " ", RowBox[{"AssociativeSimplify", "[", RowBox[{"Grassmann", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "4"}], "}"}], ",", RowBox[{"Take", "[", uis,{14}, "]"}], ",", RowBox[{"Take", "[", lis,{14}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {ta1}, ",", RowBox[{"Take", "[", uis,{5}, "]"}], ",", RowBox[{"Take", "[", lis,{5}, "]"}]}], "}"}], ",", "25"}], "]"}], "]"}], ",", " ", RowBox[{"{", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "5"}], "}"}], ",", "uis", ",", "lis"}], "}"}]}], "}"}], "}"}]}], "]"}]}], ";"}], " ", RowBox[{ RowBox[{ RowBox[{"NRta5", "[", uis_,lis_, "]"}], ":=", RowBox[{"Module", "[", RowBox[{ {tmptmp2tmp3}, ",", RowBox[{"(", " ", RowBox[{ RowBox[{"tmp", "=", " ", RowBox[{"AssociativeSimplify", "[", RowBox[{"AllTildesToNormals", "[", RowBox[{"ta5", "[", uis,lis, "]"}], "]"}], "]"}]}], ";", " ", RowBox[{"tmp2", "=", RowBox[{"AssociativeSimplify", "[", RowBox[{"AllTildesToNormals", "[", "tmp", "]"}], "]"}]}], ";", " ", RowBox[{"tmp3", "=", RowBox[{"AssociativeSimplify", "[", RowBox[{"AllTildesToNormals", "[", "tmp2", "]"}], "]"}]}], ";", " ", RowBox[{"Return", "[", RowBox[{"AllTildesToNormals", "[", "tmp3", "]"}], "]"}], ";"}], " ", ")"}]}], "]"}]}], ";"}], " ", RowBox[{ RowBox[{ RowBox[{"Gd5", "[", uis_,lis_, "]"}], ":=", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"lam", ",", "5"}], "}"}], ",", "uis", ",", "lis"}], "}"}]}], "}"}], "}"}], ",", " ", RowBox[{"Grassmann", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", RowBox[{"Take", "[", uis,{1}, "]"}], ",", RowBox[{"Take", "[", lis,{1}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"lam", ",", "4"}], "}"}], ",", RowBox[{"Take", "[", uis,{25}, "]"}], ",", RowBox[{"Take", "[", lis,{25}, "]"}]}], "}"}], ",", "25"}], "]"}], ",", " ", RowBox[{"Grassmann", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"lam", ",", "2"}], "}"}], ",", RowBox[{"Take", "[", uis,{12}, "]"}], ",", RowBox[{"Take", "[", lis,{12}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"lam", ",", "3"}], "}"}], ",", RowBox[{"Take", "[", uis,{35}, "]"}], ",", RowBox[{"Take", "[", lis,{35}, "]"}]}], "}"}], ",", "100"}], "]"}], ",", " ", RowBox[{"Grassmann3", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"} ...[TRUNC 8103]
```
**Input [178]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"Gd3", "[", {\<b\>\<d\>\<f\>},{\<c\>\<e\>\<g\>}, "]"}], "]"}]], "Input"
```
> <This Takes 10 minutes on my desktop. I chucked a timelimit in the associative simplify in Grassmann4 to let this end.>
> <We still need The gamma1^5 grassman product (actually that' s most important.) Somehow this goes in just 1.1 seconds on my desktop I have no idea why it[CloseCurlyQuote]s so much faster. >
**Input [192]**
```
[ RowBox[{ RowBox[{"Timing", "[", RowBox[{"RecognizableForm", "[", RowBox[{"Grassmann5", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", RowBox[{"Take", "[", uis,{1}, "]"}], ",", RowBox[{"Take", "[", lis,{1}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", RowBox[{"Take", "[", uis,{2}, "]"}], ",", RowBox[{"Take", "[", lis,{2}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", RowBox[{"Take", "[", uis,{3}, "]"}], ",", RowBox[{"Take", "[", lis,{3}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", RowBox[{"Take", "[", uis,{4}, "]"}], ",", RowBox[{"Take", "[", lis,{4}, "]"}]}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", RowBox[{"Take", "[", uis,{5}, "]"}], ",", RowBox[{"Take", "[", lis,{5}, "]"}]}], "}"}], ",", -120}], "]"}], "]"}], "]"}], "[", RowBox[{"[", "1", "]"}], "]"}]], "Input"
```
> So Gd5 should take about 20 seconds ...
**Input [199]**
```
[ RowBox[{ RowBox[{"Timing", "[", RowBox[{"RecognizableForm", "[", RowBox[{"Gd5", "[", {\<b\>\<d\>\<f\>\<h\>\<j\>},{\<c\>\<e\>\<g\>\<i\>\<k\>}, "]"}], "]"}], "]"}], "[", RowBox[{"[", "1", "]"}], "]"}]], "Input"
```
**Input [206]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{ RowBox[{"Gd5", "[", {\<b\>\<d\>\<f\>\<h\>\<j\>},{\<c\>\<e\>\<g\>\<i\>\<k\>}, "]"}], "//.", {NoCumulantsNormalVaccuumRule}}], "]"}]], "Input"
```
> <Which it indeed does ... Finally we are ready to complete the decomposition of a5!. >


### [215] LatexForm
**Input [223]**
```
[ RowBox[{"Timing", "[", RowBox[{ RowBox[{"temp2", "=", RowBox[{"AllStepsOfProductRule", "[", fac2,fac1, "]"}]}], ";"}], "]"}]], "Input"
```
**Input [226]**
```
[{ RowBox[{ RowBox[{"L1Terms", "=", RowBox[{"Join", "[", RowBox[{"temp1", ",", RowBox[{"MultiplyScalar", "[", temp2,-1, "]"}]}], "]"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", "L1Terms", "]"}]}], "Input"
```
**Input [229]**
```
[ RowBox[{"LatexRecognizableForm", "[", "L1Terms", "]"}]], "Input"
```
**Input [232]**
```
[{ RowBox[{"SimplifiedL1Terms", "=", RowBox[{ RowBox[{"RecognizableForm", "[", L1Terms//.{NoCumulantsNormalVaccuumRule}, "]"}], "//", "FullSimplify"}]}], " ", RowBox[{ RowBox[{"RecognizableForm", "[", temp1//.{NoCumulantsNormalVaccuumRule}, "]"}], "//", "FullSimplify"}]}], "Input"
```
**Input [236]**
```
[ RowBox[{"Length", "[", "SimplifiedL1Terms", "]"}]], "Input"
```


### [239] Triple Commutator of 1-2-2 particle operators required for Correlated dynamics Theory with Projectors.


#### [249] Everything Together
> This whole computation takes about 20 mins.
**Input [253]**
```
[{ RowBox[{ RowBox[{"com21", "=", RowBox[{ RowBox[{"AllStepsOfProductRulePProjector", "[", com1,v2, "]"}], "//.", {NoCumulantsNormalVaccuumRule}}]}], ";"}], " ", RowBox[{ RowBox[{"com22", "=", RowBox[{ RowBox[{"AllStepsOfProductRulePProjector", "[", v2,com1, "]"}], "//.", {NoCumulantsNormalVaccuumRule}}]}], ";"}], " ", RowBox[{ RowBox[{"com2", "=", RowBox[{"Join", "[", RowBox[{"com21", ",", RowBox[{"MultiplyScalar", "[", com22,-1, "]"}]}], "]"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", "com2", "]"}]}], "Input"
```
> Check the Hermiticity of Com2.
**Input [270]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"com2", "-", RowBox[{"(", RowBox[{"com2", "/.", RowBox[{"{", RowBox[{ RowBox[{"\"\<p\>\"", "->", "\"\<q\>\""}], ",", RowBox[{"\"\<q\>\"", "->", "\"\<p\>\""}]}], "}"}]}], ")"}]}], "]"}]], "Input"
```
> This uses the results from the previous section ...
**Input [276]**
```
[{ RowBox[{ RowBox[{"v1", "=", RowBox[{"NRta2", "[", {\<r\>\<s\>},{\<t\>\<u\>}, "]"}]}], ";"}], " ", RowBox[{ RowBox[{"v2", "=", RowBox[{"NRta2", "[", {\<v\>\<w\>},{\<x\>\<y\>}, "]"}]}], ";"}]}], "Input"
```
**Input [279]**
```
[ RowBox[{"RecognizableForm", "[", "WithV1andV2", "]"}]], "Input"
```
> <One way to reduce the number of terms that must be computed would be to assume a diagonal representation of the density (especially since the equations are cubic in the density). >
**Input [288]**
```
[{ RowBox[{ RowBox[{"TCLNoDelta1s", "=", RowBox[{"AssociativeSimplify", "[", RowBox[{"Map", "[", CollapseDeltas,WithV1andV2, "]"}], "]"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", "%", "]"}], " ", RowBox[{"Length", "[", "TCLNoDelta1s", "]"}]}], "Input"
```
**Input [296]**
```
[{ RowBox[{ RowBox[{"ThreeGam", "=", RowBox[{"ListForm", "[", RowBox[{ RowBox[{"2", " ", "VTs"^{{ry}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{"2", " ", "VTs"^{{rt}}_{{wx}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{"2", " ", "VTs"^{{vy}}_{{su}}, " ", "VTt"^{{rt}}_{{wx}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "+", RowBox[{"2", " ", "VTs"^{{tv}}_{{sy}}, " ", "VTt"^{{rx}}_{{uw}}, " ", g^{{"y"}}_{{"q"}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"x"}}}]}], "]"}]}], ";"}], " ", RowBox[{ RowBox[{"ThreeEta", "=", RowBox[{"ListForm", "[", RowBox[{ RowBox[{ -2, " ", "VTs"^{{rt}}_{{wx}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{"2", " ", "VTs"^{{vx}}_{{su}}, " ", "VTt"^{{rt}}_{{wy}}, " ", g^{{"y"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}, " ", eta^{{"p"}}_{{"x"}}}], "+", RowBox[{"2", " ", "VTs"^{{ry}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"y"}}}], "+", RowBox[{"2", " ", "VTs"^{{tv}}_{{sx}}, " ", "VTt"^{{ry}}_{{uw}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"y"}}}]}], "]"}]}], ";"}]}], "Input"
```
**Input [297]**
```
[ The combinations are (13), (32),(21),(44)], "Input"
```
**Input [299]**
```
[ RowBox[{ RowBox[{"2", " ", "VTs"^{{ry}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{"2", " ", "VTs"^{{ry}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"y"}}}]}]], "Input"
```
**Input [302]**
```
[ RowBox[{ RowBox[{ -2, " ", "VTs"^{{vy}}_{{su}}, " ", "VTt"^{{rt}}_{{wx}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{"2", " ", "VTs"^{{vx}}_{{su}}, " ", "VTt"^{{rt}}_{{wy}}, " ", g^{{"y"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}, " ", eta^{{"p"}}_{{"x"}}}]}]], "Input"
```
**Input [305]**
```
[ RowBox[{ RowBox[{ -2, " ", "VTs"^{{rt}}_{{wx}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{"2", " ", "VTs"^{{rt}}_{{wx}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}]}]], "Input"
```
**Input [308]**
```
[ RowBox[{ RowBox[{ +2, " ", "VTs"^{{tv}}_{{sy}}, " ", "VTt"^{{rx}}_{{uw}}, " ", g^{{"y"}}_{{"q"}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"x"}}}], "+", RowBox[{"2", " ", "VTs"^{{tv}}_{{sx}}, " ", "VTt"^{{ry}}_{{uw}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"y"}}}]}]], "Input"
```
**Input [314]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{"tmp123", "=", RowBox[{"AssociativeSimplify", "[", RowBox[{"Map", "[", CollapseDeltas,FromCode//.EtaToDelta, "]"}], "]"}]}], "]"}], " ", RowBox[{"Length", "[", "tmp123", "]"}]}], "Input"
```
**Input [318]**
```
[{ RowBox[{"RecognizableForm", "[", "FromCode", "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"FromCodeND", "=", RowBox[{"Map", "[", CollapseDeltas,FromCode, "]"}]}], "]"}], " ", RowBox[{"LatexRecognizableForm", "[", "FromCode", "]"}], " ", RowBox[{"LatexRecognizableForm", "[", "FromCodeND", "]"}]}], "Input"
```


### [324] Secular With Diagonal.
**Output [325]**
```
CellGroupData[{ Cell[[{ RowBox[{"RecognizableForm", "[", "FromCodeND", "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", CollapseDeltas,FromCodeND//.EtaToDelta, "]"}]}], "]"}], "]"}]}], "Input"
```
**Input [326]**
```
[{ RowBox[{"RecognizableForm", "[", "FromCodeND", "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", CollapseDeltas,FromCodeND//.EtaToDelta, "]"}]}], "]"}], "]"}]}], "Input"
```
**Output [327]**
```
[ RowBox[{ RowBox[{ -(1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pv}}_{{su}}, " ", "VTt"^{{rt}}_{{wx}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{qw}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{tv}}_{{qs}}, " ", "VTt"^{{rx}}_{{uw}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"x"}}}]}]], "Output"
```
**Output [328]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uv}}, " ", "VTt"^{{uv}}_{{st}}, " ", g^{{"t"}}_{{"q"}}, " ", g^{{"s"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{uv}}_{{qs}}, " ", "VTt"^{{rt}}_{{uv}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"p"}}_{{"t"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pv}}_{{su}}, " ", "VTt"^{{rt}}_{{qv}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"u"}}_{{"t"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{qv}}, " ", "VTt"^{{pv}}_{{su}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"u"}}_{{"t"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tw}}_{{sv}}, " ", g^{{"v"}}_{{"q"}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"u"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pw}}_{{su}}, " ", "VTt"^{{rt}}_{{vw}}, " ", g^{{"v"}}_{{"q"}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"u"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{qw}}, " ", "VTt"^{{vw}}_{{su}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"p"}}_{{"v"}}}], "-", RowBox[{ "VTs"^{{tw}}_{{qs}}, " ", "VTt"^{{rv}}_{{uw}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"p"}}_{{"v"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pv}}_{{su}}, " ", "VTt"^{{rt}}_{{qw}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{qw}}, " ", "VTt"^{{pv}}_{{su}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pv}}_{{su}}, " ", "VTt"^{{rt}}_{{wx}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{qw}}, " ", "VTt"^{{vx}}_{{su}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"x"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{tv}}_{{qs}}, " ", "VTt"^{{rx}}_{{uw}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"x"}}}]}]], "Output"
```
**Input [329]**
```
[{ RowBox[{ RowBox[{ RowBox[{"STwo", "[", p_, q_, r_, s_, "]"}], " ", ":=", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {r}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {s}}], "}"}]}], "}"}], ",", " ", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {r}}], "}"}]}], "}"}], ",", " ", RowBox[{"{", RowBox[{ {\<Sign\>-1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {r}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {r}}], "}"}]}], "}"}]}], "}"}]}], ";"}], " ", RowBox[{ RowBox[{ RowBox[{"SymSecATerm", "[", "ATerm_", "]"}], ":=", RowBox[{"(", " ", RowBox[{"DistributeOperatorProduct", "[", RowBox[{ RowBox[{"DistributeOperatorProduct", "[", RowBox[{ {ATerm}, ",", RowBox[{"STwo", "[", RowBox[{ RowBox[{"ATerm", "[", RowBox[{"[", 2,2,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 2,2,2, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,3,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,3,2, "]"}], "]"}]}], "]"}]}], "]"}], ",", RowBox[{"STwo", "[", RowBox[{ RowBox[{"ATerm", "[", RowBox[{"[", 3,2,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,2,2, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 2,3,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 2,3,2, "]"}], "]"}]}], "]"}]}], "]"}], " ", ")"}]}], ";"}], " ", RowBox[{ RowBox[{ RowBox[{"SymSecATerm2", "[", "ATerm_", "]"}], ":=", RowBox[{"(", " ", RowBox[{"DistributeOperatorProduct", "[", RowBox[{ RowBox[{"DistributeOperatorProduct", "[", RowBox[{ {ATerm}, ",", RowBox[{"STwo", "[", RowBox[{ RowBox[{"ATerm", "[", RowBox[{"[", 2,2,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 2,2,2, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,3,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,3,2, "]"}], "]"}]}], "]"}]}], "]"}], ",", RowBox[{"STwo", "[", RowBox[{ RowBox[{"ATerm", "[", RowBox[{"[", 3,2,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,2,2, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 2,3,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 2,3,2, "]"}], "]"}]}], "]"}]}], "]"}], " ", ")"}]}], ";"}]}], "Input"
```
**Output [330]**
```
CellGroupData[{ Cell[[ RowBox[{"FromCodeND", "[", RowBox[{"[", "1", "]"}], "]"}]], "Input"], Cell[[ RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"r"\>},{\<"u"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"t"\>\<"v"\>},{\<"s"\>\<"x"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"x"\>}, ",", {\<"q"\>}}], "}"}]}], "}"}]], "Output"
```
**Input [331]**
```
[ RowBox[{"FromCodeND", "[", RowBox[{"[", "1", "]"}], "]"}]], "Input"
```
**Output [332]**
```
[ RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"r"\>},{\<"u"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"t"\>\<"v"\>},{\<"s"\>\<"x"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"x"\>}, ",", {\<"q"\>}}], "}"}]}], "}"}]], "Output"
```
**Output [333]**
```
CellGroupData[{ Cell[[{ RowBox[{"RecognizableForm", "[", RowBox[{"STwo", "[", a,b,c,d, "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", 1;;1, "]"}], "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"SymSecATerm", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", "1", "]"}], "]"}], "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"SymSecATerm2", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", "1", "]"}], "]"}], "]"}], "]"}]}], "Input"
```
**Input [334]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{"STwo", "[", a,b,c,d, "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", 1;;1, "]"}], "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"SymSecATerm", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", "1", "]"}], "]"}], "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"SymSecATerm2", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", "1", "]"}], "]"}], "]"}], "]"}]}], "Input"
```
**Output [335]**
```
[ RowBox[{ RowBox[{ d^{{b}}_{{c}}, " ", d^{{a}}_{{d}}}], "+", RowBox[{ d^{{a}}_{{c}}, " ", d^{{b}}_{{d}}}], "-", RowBox[{ d^{{a}}_{{c}}, " ", d^{{b}}_{{c}}, " ", d^{{a}}_{{d}}, " ", d^{{b}}_{{d}}}]}]], "Output"
```
**Output [336]**
```
[ RowBox[{ -(1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}]], "Output"
```
**Output [337]**
```
[ RowBox[{ RowBox[{ -(1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"r"}}_{{"s"}}, " ", d^{{"v"}}_{{"u"}}, " ", d^{{"t"}}_{{"w"}}, " ", d^{{"p"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"r"}}_{{"s"}}, " ", d^{{"t"}}_{{"u"}}, " ", d^{{"v"}}_{{"w"}}, " ", d^{{"p"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"r"}}_{{"s"}}, " ", d^{{"t"}}_{{"u"}}, " ", d^{{"v"}}_{{"u"}}, " ", d^{{"t"}}_{{"w"}}, " ", d^{{"v"}}_{{"w"}}, " ", d^{{"p"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"s"}}, " ", d^{{"v"}}_{{"u"}}, " ", d^{{"t"}}_{{"w"}}, " ", d^{{"r"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"s"}}, " ", d^{{"t"}}_{{"u"}}, " ", d^{{"v"}}_{{"w"}}, " ", d^{{"r"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"s"}}, " ", d^{{"t"}}_{{"u"}}, " ", d^{{"v"}}_{{"u"}}, " ", d^{{"t"}}_{{"w"}}, " ", d^{{"v"}}_{{"w"}}, " ", d^{{"r"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"s"}}, " ", d^{{"r"}}_{{"s"}}, " ", d^{{"v"}}_{{"u"}}, " ", d^{{"t"}}_{{"w"}}, " ", d^{{"p"}}_{{"x"}}, " ", d^{{"r"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"s"}}, " ", d^{{"r"}}_{{"s"}}, " ", d^{{"t"}}_{{"u"}}, " ", d^{{"v"}}_{{"w"}}, " ", d^{{"p"}}_{{"x"}}, " ", d^{{"r"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"s"}}, " ", d^{{"r"}}_{{"s"}}, " ", d^{{"t"}}_{{"u"}}, " ", d^{{"v"}}_{{"u"}}, " ", d^{{"t"}}_{{"w"}}, " ", d^{{"v"}}_{{"w"}}, " ", d^{{"p"}}_{{"x"}}, " ", d^{{"r"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}]}]], "Output"
```
**Output [338]**
```
[ RowBox[{ RowBox[{ -(1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"r"}}_{{"s"}}, " ", d^{{"v"}}_{{"u"}}, " ", d^{{"t"}}_{{"w"}}, " ", d^{{"p"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"r"}}_{{"s"}}, " ", d^{{"t"}}_{{"u"}}, " ", d^{{"v"}}_{{"w"}}, " ", d^{{"p"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"r"}}_{{"s"}}, " ", d^{{"t"}}_{{"u"}}, " ", d^{{"v"}}_{{"u"}}, " ", d^{{"t"}}_{{"w"}}, " ", d^{{"v"}}_{{"w"}}, " ", d^{{"p"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"s"}}, " ", d^{{"v"}}_{{"u"}}, " ", d^{{"t"}}_{{"w"}}, " ", d^{{"r"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"s"}}, " ", d^{{"t"}}_{{"u"}}, " ", d^{{"v"}}_{{"w"}}, " ", d^{{"r"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"s"}}, " ", d^{{"t"}}_{{"u"}}, " ", d^{{"v"}}_{{"u"}}, " ", d^{{"t"}}_{{"w"}}, " ", d^{{"v"}}_{{"w"}}, " ", d^{{"r"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"s"}}, " ", d^{{"r"}}_{{"s"}}, " ", d^{{"v"}}_{{"u"}}, " ", d^{{"t"}}_{{"w"}}, " ", d^{{"p"}}_{{"x"}}, " ", d^{{"r"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"s"}}, " ", d^{{"r"}}_{{"s"}}, " ", d^{{"t"}}_{{"u"}}, " ", d^{{"v"}}_{{"w"}}, " ", d^{{"p"}}_{{"x"}}, " ", d^{{"r"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"s"}}, " ", d^{{"r"}}_{{"s"}}, " ", d^{{"t"}}_{{"u"}}, " ", d^{{"v"}}_{{"u"}}, " ", d^{{"t"}}_{{"w"}}, " ", d^{{"v"}}_{{"w"}}, " ", d^{{"p"}}_{{"x"}}, " ", d^{{"r"}}_{{"x"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}]}]], "Output"
```
**Output [342]**
```
CellGroupData[{ Cell[[ RowBox[{"RecognizableForm", "[", RowBox[{"Sec8FromCodeND", "=", RowBox[{"AssociativeSimplify", "[", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SecularizeATerm,FromCodeND, "]"}], ",", "1"}], "]"}]}], "]"}]}], "]"}], "]"}]}], "]"}]], "Input"
```
**Input [343]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"Sec8FromCodeND", "=", RowBox[{"AssociativeSimplify", "[", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SecularizeATerm,FromCodeND, "]"}], ",", "1"}], "]"}]}], "]"}]}], "]"}], "]"}]}], "]"}]], "Input"
```
**Output [344]**
```
[ RowBox[{ RowBox[{ RowBox[{"-", "VTs"^{{pr}}_{{pr}}}], " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"p"}}_{{"r"}}, " ", eta^{{"r"}}_{{"p"}}, " ", eta^{{"p"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"r"}}_{{"p"}}, " ", eta^{{"p"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{ps}}_{{pr}}, " ", "VTt"^{{pr}}_{{ps}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"r"}}_{{"p"}}, " ", eta^{{"p"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"r"}}_{{"p"}}, " ", g^{{"p"}}_{{"r"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"q"}}_{{"r"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"q"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{ps}}, " ", "VTt"^{{ps}}_{{pr}}, " ", g^{{"s"}}_{{"p"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qs}}_{{qr}}, " ", "VTt"^{{qr}}_{{qs}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"q"}}_{{"s"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{rs}}, " ", "VTt"^{{rs}}_{{pr}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{rs}}_{{qr}}, " ", "VTt"^{{qr}}_{{rs}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"p"}}_{{"r"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"q"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"q"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{qr}}_{{qs}}, " ", "VTt"^{{qs}}_{{qr}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"q"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"p"}}_{{"r"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{ps}}_{{pr}}, " ", "VTt"^{{pr}}_{{ps}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{ps}}_{{pr}}, " ", "VTt"^{{pr}}_{{ps}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"r"}}_{{"p"}}, " ", g^{{"p"}}_{{"r"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"q"}}_{{"r"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"q"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{ps}}, " ", "VTt"^{{ps}}_{{pr}}, " ", g^{{"s"}}_{{"p"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{qs}}_{{qr}}, " ", "VTt"^{{qr}}_{{qs}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"q"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{rs}}, " ", "VTt"^{{rs}}_{{pr}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{rs}}_{{qr}}, " ", "VTt"^{{qr}}_{{rs}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{ps}}, " ", "VTt"^{{ps}}_{{pr}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qs}}_{{qr}}, " ", "VTt"^{{qr}}_{{qs}}, " ", g^{{"q"}}_{{"q"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{rs}}, " ", "VTt"^{{rs}}_{{pr}}, " ", g^{{"r"}}_{{"r"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{rs}}_{{qr}}, " ", "VTt"^{{qr}}_{{rs}}, " ", g^{{"r"}}_{{"r"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"s"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{qr}}, " ", "VTt"^{{qr}}_{{st}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"s"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{qr}}, " ", "VTt"^{{qr}}_{{st}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"p"}}_{{"r"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qr}}_{{qs}}, " ", "VTt"^{{qs}}_{{ ...[TRUNC 24143]
```
> Check that' s still got the same diagonal EOM
**Output [346]**
```
CellGroupData[{ Cell[[ RowBox[{"RecognizableForm", "[", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", DiagonalGamma1,SecFromCodeND, "]"}]}], "]"}]}], "]"}], "]"}]], "Input"
```
**Input [347]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", DiagonalGamma1,SecFromCodeND, "]"}]}], "]"}]}], "]"}], "]"}]], "Input"
```
**Output [348]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{pr}}, " ", "VTt"^{{pr}}_{{st}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pt}}_{{rs}}, " ", "VTt"^{{rs}}_{{pt}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rs}}_{{pt}}, " ", "VTt"^{{pt}}_{{rs}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}]}]], "Output"
```


#### [350] Simple Code for all these contractions
**Input [357]**
```
[ RowBox[{"CodeTerm", "[", RowBox[{"spinsec8", "[", RowBox[{"[", "1", "]"}], "]"}], "]"}]], "Input"
```
**Input [361]**
```
[{ RowBox[{"Length", "[", "HasQ", "]"}], " ", RowBox[{"Length", "[", "HasNoQ", "]"}]}], "Input"
```
**Input [365]**
```
[{ RowBox[{"Length", "[", "HasQT", "]"}], " ", RowBox[{"Length", "[", "HasQNoT", "]"}], " ", RowBox[{"Length", "[", "HasTNoQ", "]"}]}], "Input"
```
**Output [373]**
```
CellGroupData[{ Cell[[ RowBox[{"HasQTandQinVs", "[", RowBox[{"[", "1", "]"}], "]"}]], "Input"
```
**Input [374]**
```
[ RowBox[{"HasQTandQinVs", "[", RowBox[{"[", "1", "]"}], "]"}]], "Input"
```
**Output [375]**
```
[ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"r"\>},{\<"s"\>\<"t"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"s"\>\<"t"\>},{\<"p"\>\<"r"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"s"\>\<"t"\>},{\<"q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"r"\>},{\<"s"\>\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"r"\>},{\<"s"\>\<"t"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"s"\>\<"t"\>},{\<"p"\>\<"r"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"s"\>\<"t"\>},{\<"q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"r"\>},{\<"s"\>\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"t"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"r"\>\<"s"\>},{\<"p"\>\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"s"\>},{\<"q"\>\<"t"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"t"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"t"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"r"\>\<"s"\>},{\<"p"\>\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"s"\>},{\<"q"\>\<"t"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"t"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"t"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"r"\>\<"s"\>},{\<"p"\>\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)} ...[TRUNC 10589]
```


#### [384] Secular8 Version.
**Input [391]**
```
[ RowBox[{"CodeTerm", "[", RowBox[{"spinsec8", "[", RowBox[{"[", "1", "]"}], "]"}], "]"}]], "Input"
```
**Input [395]**
```
[{ RowBox[{"Length", "[", "HasQ", "]"}], " ", RowBox[{"Length", "[", "HasNoQ", "]"}]}], "Input"
```
**Input [399]**
```
[{ RowBox[{"Length", "[", "HasQT", "]"}], " ", RowBox[{"Length", "[", "HasQNoT", "]"}], " ", RowBox[{"Length", "[", "HasTNoQ", "]"}]}], "Input"
```
**Output [407]**
```
CellGroupData[{ Cell[[ RowBox[{"HasQTandQinVs", "[", RowBox[{"[", "1", "]"}], "]"}]], "Input"
```
**Input [408]**
```
[ RowBox[{"HasQTandQinVs", "[", RowBox[{"[", "1", "]"}], "]"}]], "Input"
```
**Output [409]**
```
[ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"s"\>},{\<"p"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"p"\>\<"r"\>},{\<"p"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"p"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"r"\>},{\<"p"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"p"\>\<"s"\>},{\<"p"\>\<"r"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"p"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"s"\>},{\<"q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"q"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"r"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"r"\>\<"s"\>},{\<"p"\>\<"r"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"s"\>},{\<"q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"r"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"s"\>},{\<"q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"q"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"s"\>},{\<"p"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"p"\>\<"r"\>},{\<"p"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"p"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"s"\>},{\<"p"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"p"\>\<"r"\>},{\<"p"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"p"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"r"\>},{\<"p"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"p"\>\<"s"\>},{\<"p"\>\<"r"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"p"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"s"\>},{\<"q"\>\ ...[TRUNC 58403]
```
**Output [417]**
```
CellGroupData[{ Cell["SecularTwo Secular approximation Used for off diagonal. ", "Subsection"
```


### [418] SecularTwo Secular approximation Used for off diagonal.
**Output [419]**
```
CellGroupData[{ Cell[[ RowBox[{"RecognizableForm", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", 2;;2, "]"}], "]"}], "]"}]], "Input"
```
**Input [420]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", 2;;2, "]"}], "]"}], "]"}]], "Input"
```
**Output [421]**
```
[ RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{qw}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}]], "Output"
```
**Input [422]**
```
[{ RowBox[{ RowBox[{ RowBox[{"STwo", "[", p_, q_, r_, s_, "]"}], " ", ":=", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {r}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {s}}], "}"}]}], "}"}], ",", " ", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {r}}], "}"}]}], "}"}], ",", " ", RowBox[{"{", RowBox[{ {\<Sign\>-1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {r}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {r}}], "}"}]}], "}"}]}], "}"}]}], ";"}], " ", RowBox[{ RowBox[{"SymSecATerm", "[", "ATerm_", "]"}], ":=", RowBox[{"(", " ", RowBox[{"DistributeOperatorProduct", "[", RowBox[{ RowBox[{"DistributeOperatorProduct", "[", RowBox[{ {ATerm}, ",", RowBox[{"STwo", "[", RowBox[{ RowBox[{"ATerm", "[", RowBox[{"[", 2,2,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,2,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 2,3,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,3,1, "]"}], "]"}]}], "]"}]}], "]"}], ",", RowBox[{"STwo", "[", RowBox[{ RowBox[{"ATerm", "[", RowBox[{"[", 2,2,2, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,2,2, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 2,3,2, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,3,2, "]"}], "]"}]}], "]"}]}], "]"}], " ", ")"}]}], " ", RowBox[{ RowBox[{"SymSecATermTwo", "[", "ATerm_", "]"}], ":=", RowBox[{"(", " ", RowBox[{"DistributeOperatorProduct", "[", RowBox[{ RowBox[{"DistributeOperatorProduct", "[", RowBox[{ {ATerm}, ",", RowBox[{"STwo", "[", RowBox[{ RowBox[{"ATerm", "[", RowBox[{"[", 2,2,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,2,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 2,3,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,3,1, "]"}], "]"}]}], "]"}]}], "]"}], ",", RowBox[{"STwo", "[", RowBox[{ RowBox[{"ATerm", "[", RowBox[{"[", 2,2,2, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,2,2, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 2,3,2, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,3,2, "]"}], "]"}]}], "]"}]}], "]"}], " ", ")"}]}]}], "Input"
```
**Output [423]**
```
CellGroupData[{ Cell[[ RowBox[{"RecognizableForm", "[", RowBox[{"symsectcl2", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"SymSecATerm", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", "2", "]"}], "]"}], "]"}]}], "]"}]}], "]"}]}], "]"}]], "Input"
```
**Input [424]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"symsectcl2", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"SymSecATerm", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", "2", "]"}], "]"}], "]"}]}], "]"}]}], "]"}]}], "]"}]], "Input"
```
**Output [425]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"p"}}_{{"r"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{qr}}_{{qs}}, " ", "VTt"^{{qs}}_{{qr}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{rs}}_{{rs}}, " ", g^{{"p"}}_{{"r"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"s"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{qr}}_{{qt}}, " ", "VTt"^{{st}}_{{rs}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"p"}}_{{"t"}}, " ", eta^{{"s"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qs}}_{{qs}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{st}}_{{st}}, " ", g^{{"p"}}_{{"s"}}, " ", g^{{"r"}}_{{"t"}}, " ", eta^{{"t"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rs}}_{{qr}}, " ", "VTt"^{{qr}}_{{rs}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"p"}}_{{"r"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{rs}}_{{qt}}, " ", "VTt"^{{qt}}_{{rs}}, " ", g^{{"t"}}_{{"q"}}, " ", g^{{"p"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{qr}}, " ", "VTt"^{{qs}}_{{st}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"t"}}_{{"t"}}}]}]], "Output"
```
**Output [426]**
```
CellGroupData[{ Cell[[ RowBox[{"RecognizableForm", "[", RowBox[{"AssociativeSimplify", "[", RowBox[{"SpinLabelTerms", "[", "symsectcl2", "]"}], "]"}], "]"}]], "Input",\ CellChangeTimes->{{3.670235699649184*^9, 3.670235721829818*^9}}], Cell[[ InterpretationBox[ RowBox[{"\<\"Unable to Parse Expression back into List...\"\>", "", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"\<\"VTs\"\>", ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"r"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"s"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"s"\>},{\<"q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"s"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"s"\>},{\<"q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"s"\>},{\<"s"\>\<"t"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"t"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"t"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"s"\>},{\<"q"\>\<"t"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"r"\>\<"s"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"s"\>\<"t"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"t"\>}, "}"}]}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"2", ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"r"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"s"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"s"\>},{\<"q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"s"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"s"\>},{\<"q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"s"\>},{\<"s"\>\<"t"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"t"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"t"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"s"\>},{\<"q"\>\<"t"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"r"\>\<"s"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"s"\>\<"t"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"t"\>}, "}"}]}], "}"}]}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"\<\"q\"\>", ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"r"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"s"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"s"\>},{\<"q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"s"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"s"\>},{\<"q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"s"\>},{\<"s"\>\<"t"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"t"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"t"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"s"\>},{\<"q"\>\<"t"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"r"\>\<"s"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"s"\>\<"t"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"t"\>}, "}"}]}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"r\"\>", ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"r"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"s"\>},{\<"q"\>\<"r"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"s"\>},{\<"q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"q"\>\<"r"\>},{\<"q"\>\<"s"\>}, "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"s"\>},{\<"q"\>\<"s"\> ...[TRUNC 14563]
```
**Input [427]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"AssociativeSimplify", "[", RowBox[{"SpinLabelTerms", "[", "symsectcl2", "]"}], "]"}], "]"}]], "Input",\ CellChangeTimes->{{3.670235699649184*^9, 3.670235721829818*^9}}
```
**Output [429]**
```
[ InterpretationBox[ TagBox[ FrameBox[GridBox[{ { ItemBox[ TagBox[ RowBox[{"{", RowBox[{ RowBox[{ RowBox[{"\<\"Error\"\>", "[[", 1,2, "]]"}], "+", RowBox[{"2", " ", SuperscriptBox[ RowBox[{"(", SubsuperscriptBox["\<\"q\"\>", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], "}"}], ",", RowBox[{"{", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], "}"}]}], "}"}], ",", "3"}], "}"}], RowBox[{"{", RowBox[{"\<\"r\"\>", ",", TemplateBox[{"7"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], "}"}], ",", RowBox[{"{", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], "}"}]}], "}"}]}], "}"}]], ")"}], "2"]}]}], ",", TemplateBox[{"7"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{"\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], ",", RowBox[{"{", RowBox[{ RowBox[{ RowBox[{"\<\"Error\"\>", "[[", 1,2, "]]"}], "+", RowBox[{"\<\"Sign\"\>", " ", SuperscriptBox[ RowBox[{"(", SubsuperscriptBox["\<\"q\"\>", RowBox[{"{", RowBox[{ RowBox[{"{", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], "}"}], ",", "3"}], "}"}], RowBox[{"{", RowBox[{"\<\"r\"\>", ",", TemplateBox[{"7"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], ",", RowBox[{"{", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], "}"}]}], "}"}]], ")"}], "2"]}]}], ",", RowBox[{"{", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{ RowBox[{"\<\"Error\"\>", "[[", 1,2, "]]"}], "+", RowBox[{"\<\"VTs\"\>", " ", SuperscriptBox[ RowBox[{"(", SubsuperscriptBox["\<\"q\"\>", RowBox[{"{", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], "}"}], RowBox[{"{", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], "}"}]], ")"}], "2"]}]}], ",", RowBox[{ RowBox[{"\<\"Error\"\>", "[[", 1,2, "]]"}], "+", RowBox[{"2", " ", SuperscriptBox[ RowBox[{"(", SubsuperscriptBox["\<\"q\"\>", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor ->  ...[TRUNC 13034]
```
**Output [430]**
```
CellGroupData[{ Cell["The diagonal/secular approximation I like best. ", "Subsection"
```


### [431] The diagonal/secular approximation I like best.
**Output [437]**
```
CellGroupData[{ Cell[[{ RowBox[{"RecognizableForm", "[", RowBox[{"dsec", "=", RowBox[{"Map", "[", RowBox[{"RelabelTerm", ",", RowBox[{"ListForm", "[", RowBox[{ RowBox[{ -(1)/(2), " ", "VTs"^{{pr}}_{{uv}}, " ", "VTt"^{{uv}}_{{rp}}, " ", g^{{"u"}}_{{"u"}}, " ", g^{{"v"}}_{{"v"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pv}}_{{rt}}, " ", "VTt"^{{rt}}_{{vp}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"v"}}_{{"v"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"t"}}_{{"t"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{qv}}, " ", "VTt"^{{vq}}_{{rt}}, " ", g^{{"v"}}_{{"v"}}, " ", g^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"t"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{tv}}_{{qr}}, " ", "VTt"^{{rq}}_{{tv}}, " ", g^{{"t"}}_{{"t"}}, " ", g^{{"v"}}_{{"v"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"p"}}_{{"q"}}}]}], "]"}]}], "]"}]}], "]"}], " ", RowBox[{"Map", "[", NumIndicesInTerm,dsec, "]"}], " ", RowBox[{"LatexRecognizableForm", "[", "dsec", "]"}]}], "Input"
```
**Input [438]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{"dsec", "=", RowBox[{"Map", "[", RowBox[{"RelabelTerm", ",", RowBox[{"ListForm", "[", RowBox[{ RowBox[{ -(1)/(2), " ", "VTs"^{{pr}}_{{uv}}, " ", "VTt"^{{uv}}_{{rp}}, " ", g^{{"u"}}_{{"u"}}, " ", g^{{"v"}}_{{"v"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pv}}_{{rt}}, " ", "VTt"^{{rt}}_{{vp}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"v"}}_{{"v"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"t"}}_{{"t"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{qv}}, " ", "VTt"^{{vq}}_{{rt}}, " ", g^{{"v"}}_{{"v"}}, " ", g^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"t"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{tv}}_{{qr}}, " ", "VTt"^{{rq}}_{{tv}}, " ", g^{{"t"}}_{{"t"}}, " ", g^{{"v"}}_{{"v"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"p"}}_{{"q"}}}]}], "]"}]}], "]"}]}], "]"}], " ", RowBox[{"Map", "[", NumIndicesInTerm,dsec, "]"}], " ", RowBox[{"LatexRecognizableForm", "[", "dsec", "]"}]}], "Input"
```
**Output [439]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{qr}}, " ", "VTt"^{{qr}}_{{st}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pt}}_{{rs}}, " ", "VTt"^{{rs}}_{{pt}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rs}}_{{qt}}, " ", "VTt"^{{qt}}_{{rs}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}]}]], "Output"
```
**Output [440]**
```
[ RowBox[{"{", RowBox[{ RowBox[{"NumIndicesInTerm", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"s"\>\<"t"\>},{\<"q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"r"\>},{\<"s"\>\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], "]"}], ",", RowBox[{"NumIndicesInTerm", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"r"\>},{\<"s"\>\<"t"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"s"\>\<"t"\>},{\<"p"\>\<"r"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], "]"}], ",", RowBox[{"NumIndicesInTerm", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"s"\>},{\<"q"\>\<"t"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"q"\>\<"t"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], "]"}], ",", RowBox[{"NumIndicesInTerm", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"t"\>},{\<"r"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"r"\>\<"s"\>},{\<"p"\>\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], "]"}]}], "}"}]], "Output"
```
**Output [441]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "\<\"\\\\gamma_{s}^{s}\"\>", " ", "\<\"\\\\gamma_{t}^{t}\"\>", " ", "\<\"VTs_{s t}^{p r}\"\>", " ", "\<\"VTt_{p r}^{s t}\"\>", " ", "\<\"eta_{q}^{p}\"\>", " ", "\<\"eta_{r}^{r}\"\>"}], "+", RowBox[{ (1)/(2), " ", "\<\"\\\\gamma_{s}^{s}\"\>", " ", "\<\"\\\\gamma_{t}^{t}\"\>", " ", "\<\"VTs_{q r}^{s t}\"\>", " ", "\<\"VTt_{s t}^{q r}\"\>", " ", "\<\"eta_{q}^{p}\"\>", " ", "\<\"eta_{r}^{r}\"\>"}], "-", RowBox[{ (1)/(2), " ", "\<\"\\\\gamma_{q}^{p}\"\>", " ", "\<\"\\\\gamma_{t}^{t}\"\>", " ", "\<\"VTs_{r s}^{p t}\"\>", " ", "\<\"VTt_{p t}^{r s}\"\>", " ", "\<\"eta_{r}^{r}\"\>", " ", "\<\"eta_{s}^{s}\"\>"}], "-", RowBox[{ (1)/(2), " ", "\<\"\\\\gamma_{q}^{p}\"\>", " ", "\<\"\\\\gamma_{t}^{t}\"\>", " ", "\<\"VTs_{q t}^{r s}\"\>", " ", "\<\"VTt_{r s}^{q t}\"\>", " ", "\<\"eta_{r}^{r}\"\>", " ", "\<\"eta_{s}^{s}\"\>"}]}]], "Output"
```
> <Only the last term is the only nonzero if pp is empty ... should be like a relaxation process into p,r. Filled levels s and t move into p and r. so this term should be strictly positive. For that to be true, p,r must be higher in energy than s,t. actually if the secularizer enforced zero denominator; that would work. >
**Output [448]**
```
CellGroupData[{ Cell[[ RowBox[{"RecognizableForm", "[", RowBox[{"AssociativeSimplify", "[", RowBox[{"(", RowBox[{"SpinLabelTerms", "[", "ddsec", "]"}], ")"}], "]"}], "]"}]], "Input"
```
**Input [449]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"AssociativeSimplify", "[", RowBox[{"(", RowBox[{"SpinLabelTerms", "[", "ddsec", "]"}], ")"}], "]"}], "]"}]], "Input"
```
**Output [451]**
```
[ InterpretationBox[ TagBox[ FrameBox[GridBox[{ { ItemBox[ TagBox[ RowBox[{"{", RowBox[{ RowBox[{ RowBox[{"\<\"Error\"\>", "[[", 1,2, "]]"}], "+", RowBox[{"2", " ", SubsuperscriptBox["\<\"p\"\>", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], "}"}], ",", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )]}], "}"}], ",", "3"}], "}"}], RowBox[{"{", RowBox[{"\<\"r\"\>", ",", RowBox[{"{", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], "}"}], ",", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], ",", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )]}], "}"}]], " ", SubsuperscriptBox["\<\"s\"\>", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"\<\"s\"\>", ",", RowBox[{"{", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], "}"}], ",", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"t"\>},{\<"r"\>\<"s"\>}, "}"}]}], "}"}]}], "}"}], ",", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )]}], "}"}], ",", "3"}], "}"}], RowBox[{"{", RowBox[{"\<\"t\"\>", ",", RowBox[{"{", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], "}"}], ",", RowBox[{"{", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTt"\>2},{\<"r"\>\<"s"\>},{\<"p"\>\<"t"\>}, "}"}], ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"t"\>},{\<"r"\>\<"s"\>}, "}"}]}], "}"}]}], "}"}]]}]}], ",", TemplateBox[{"3"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{"\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )]}], "}"}], Short[#, 5]& ], BaseStyle->{Deployed -> False}, StripOnInput->False]}, {GridBox[{ { TagBox[ TooltipBox[ StyleBox[ StyleBox[ DynamicBox[ToBoxes[ FEPrivate`FrontEndResource[ "FEStrings", "sizeBriefExplanation"], StandardForm], ImageSizeCache->{58., {2., 8.}}], StripOnInput->False, DynamicUpdating->True], "SuggestionsBarCategory", StripOnInput->False], StyleBox[ DynamicBox[ ToBoxes[ FEPrivate`FrontEndResource["FEStrings", "sizeExplanation"], St ...[TRUNC 10137]
```


### [453] Spin Adaptation
**Input [458]**
```
[ RowBox[{"Length", "[", "AllTermsSpinLabeled", "]"}]], "Input"
```
**Output [466]**
```
CellGroupData[{ Cell[[{ RowBox[{"RecognizableForm", "[", RowBox[{ RowBox[{"SpinLabelTerm", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", "1", "]"}], "]"}], "]"}], "//.", {SingleBarRule1SingleBarRule2SingleBarRule3SingleBarRule4}}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{ RowBox[{"SpinLabelTerm", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", "2", "]"}], "]"}], "]"}], "//.", {SingleBarRule1SingleBarRule2SingleBarRule3SingleBarRule4}}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{ RowBox[{"SpinLabelTerm", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", "3", "]"}], "]"}], "]"}], "//.", {SingleBarRule1SingleBarRule2SingleBarRule3SingleBarRule4}}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{ RowBox[{"SpinLabelTerm", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", "4", "]"}], "]"}], "]"}], "//.", {SingleBarRule1SingleBarRule2SingleBarRule3SingleBarRule4}}], "]"}]}], "Input"], Cell[[ RowBox[{ RowBox[{ "vTs"^{{rP}}_{{wU}}, " ", "vTt"^{{vT}}_{{sX}}, " ", g^{{"U"}}_{{"T"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"X"}}_{{"Q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ "vTs"^{{rP}}_{{uW}}, " ", "vTt"^{{tV}}_{{sX}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"X"}}_{{"Q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ "vTs"^{{RP}}_{{UW}}, " ", "vTt"^{{VT}}_{{SX}}, " ", g^{{"U"}}_{{"T"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"X"}}_{{"Q"}}, " ", eta^{{"S"}}_{{"R"}}}], "+", RowBox[{ "vTs"^{{RP}}_{{WU}}, " ", "vTt"^{{VT}}_{{SX}}, " ", g^{{"U"}}_{{"T"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"X"}}_{{"Q"}}, " ", eta^{{"S"}}_{{"R"}}}], "+", RowBox[{ "vTs"^{{RP}}_{{UW}}, " ", "vTt"^{{VT}}_{{XS}}, " ", g^{{"U"}}_{{"T"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"X"}}_{{"Q"}}, " ", eta^{{"S"}}_{{"R"}}}], "-", RowBox[{ "vTs"^{{RP}}_{{WU}}, " ", "vTt"^{{VT}}_{{XS}}, " ", g^{{"U"}}_{{"T"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"X"}}_{{"Q"}}, " ", eta^{{"S"}}_{{"R"}}}]}]], "Output"
```
**Input [467]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{ RowBox[{"SpinLabelTerm", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", "1", "]"}], "]"}], "]"}], "//.", {SingleBarRule1SingleBarRule2SingleBarRule3SingleBarRule4}}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{ RowBox[{"SpinLabelTerm", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", "2", "]"}], "]"}], "]"}], "//.", {SingleBarRule1SingleBarRule2SingleBarRule3SingleBarRule4}}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{ RowBox[{"SpinLabelTerm", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", "3", "]"}], "]"}], "]"}], "//.", {SingleBarRule1SingleBarRule2SingleBarRule3SingleBarRule4}}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{ RowBox[{"SpinLabelTerm", "[", RowBox[{"FromCodeND", "[", RowBox[{"[", "4", "]"}], "]"}], "]"}], "//.", {SingleBarRule1SingleBarRule2SingleBarRule3SingleBarRule4}}], "]"}]}], "Input"
```
**Output [468]**
```
[ RowBox[{ RowBox[{ "vTs"^{{rP}}_{{wU}}, " ", "vTt"^{{vT}}_{{sX}}, " ", g^{{"U"}}_{{"T"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"X"}}_{{"Q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ "vTs"^{{rP}}_{{uW}}, " ", "vTt"^{{tV}}_{{sX}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"X"}}_{{"Q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ "vTs"^{{RP}}_{{UW}}, " ", "vTt"^{{VT}}_{{SX}}, " ", g^{{"U"}}_{{"T"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"X"}}_{{"Q"}}, " ", eta^{{"S"}}_{{"R"}}}], "+", RowBox[{ "vTs"^{{RP}}_{{WU}}, " ", "vTt"^{{VT}}_{{SX}}, " ", g^{{"U"}}_{{"T"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"X"}}_{{"Q"}}, " ", eta^{{"S"}}_{{"R"}}}], "+", RowBox[{ "vTs"^{{RP}}_{{UW}}, " ", "vTt"^{{VT}}_{{XS}}, " ", g^{{"U"}}_{{"T"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"X"}}_{{"Q"}}, " ", eta^{{"S"}}_{{"R"}}}], "-", RowBox[{ "vTs"^{{RP}}_{{WU}}, " ", "vTt"^{{VT}}_{{XS}}, " ", g^{{"U"}}_{{"T"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"X"}}_{{"Q"}}, " ", eta^{{"S"}}_{{"R"}}}]}]], "Output"
```
**Output [469]**
```
[ RowBox[{ RowBox[{ RowBox[{"-", "vTs"^{{tR}}_{{wQ}}}], " ", "vTt"^{{vY}}_{{uS}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"P"}}_{{"Y"}}, " ", eta^{{"S"}}_{{"R"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{ "vTs"^{{rT}}_{{wQ}}, " ", "vTt"^{{vY}}_{{sU}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"P"}}_{{"Y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"U"}}_{{"T"}}}], "+", RowBox[{ "vTs"^{{TR}}_{{QW}}, " ", "vTt"^{{YV}}_{{SU}}, " ", g^{{"W"}}_{{"V"}}, " ", g^{{"P"}}_{{"Y"}}, " ", eta^{{"S"}}_{{"R"}}, " ", eta^{{"U"}}_{{"T"}}}], "-", RowBox[{ "vTs"^{{TR}}_{{WQ}}, " ", "vTt"^{{YV}}_{{SU}}, " ", g^{{"W"}}_{{"V"}}, " ", g^{{"P"}}_{{"Y"}}, " ", eta^{{"S"}}_{{"R"}}, " ", eta^{{"U"}}_{{"T"}}}], "-", RowBox[{ "vTs"^{{TR}}_{{QW}}, " ", "vTt"^{{YV}}_{{US}}, " ", g^{{"W"}}_{{"V"}}, " ", g^{{"P"}}_{{"Y"}}, " ", eta^{{"S"}}_{{"R"}}, " ", eta^{{"U"}}_{{"T"}}}], "+", RowBox[{ "vTs"^{{TR}}_{{WQ}}, " ", "vTt"^{{YV}}_{{US}}, " ", g^{{"W"}}_{{"V"}}, " ", g^{{"P"}}_{{"Y"}}, " ", eta^{{"S"}}_{{"R"}}, " ", eta^{{"U"}}_{{"T"}}}]}]], "Output"
```
**Output [470]**
```
[ RowBox[{ RowBox[{ RowBox[{"-", "vTs"^{{vP}}_{{uS}}}], " ", "vTt"^{{tR}}_{{wX}}, " ", g^{{"X"}}_{{"Q"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"S"}}_{{"R"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{ "vTs"^{{vP}}_{{sU}}, " ", "vTt"^{{rT}}_{{wX}}, " ", g^{{"X"}}_{{"Q"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"U"}}_{{"T"}}}], "+", RowBox[{ "vTs"^{{VP}}_{{SU}}, " ", "vTt"^{{TR}}_{{WX}}, " ", g^{{"X"}}_{{"Q"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"S"}}_{{"R"}}, " ", eta^{{"U"}}_{{"T"}}}], "-", RowBox[{ "vTs"^{{VP}}_{{US}}, " ", "vTt"^{{TR}}_{{WX}}, " ", g^{{"X"}}_{{"Q"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"S"}}_{{"R"}}, " ", eta^{{"U"}}_{{"T"}}}], "-", RowBox[{ "vTs"^{{VP}}_{{SU}}, " ", "vTt"^{{TR}}_{{XW}}, " ", g^{{"X"}}_{{"Q"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"S"}}_{{"R"}}, " ", eta^{{"U"}}_{{"T"}}}], "+", RowBox[{ "vTs"^{{VP}}_{{US}}, " ", "vTt"^{{TR}}_{{XW}}, " ", g^{{"X"}}_{{"Q"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"S"}}_{{"R"}}, " ", eta^{{"U"}}_{{"T"}}}]}]], "Output"
```
**Output [471]**
```
[ RowBox[{ RowBox[{ "vTs"^{{vT}}_{{sQ}}, " ", "vTt"^{{rX}}_{{wU}}, " ", g^{{"U"}}_{{"T"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"P"}}_{{"X"}}}], "+", RowBox[{ "vTs"^{{tV}}_{{sQ}}, " ", "vTt"^{{rX}}_{{uW}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"P"}}_{{"X"}}}], "-", RowBox[{ "vTs"^{{VT}}_{{QS}}, " ", "vTt"^{{XR}}_{{UW}}, " ", g^{{"U"}}_{{"T"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"S"}}_{{"R"}}, " ", eta^{{"P"}}_{{"X"}}}], "+", RowBox[{ "vTs"^{{VT}}_{{SQ}}, " ", "vTt"^{{XR}}_{{UW}}, " ", g^{{"U"}}_{{"T"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"S"}}_{{"R"}}, " ", eta^{{"P"}}_{{"X"}}}], "+", RowBox[{ "vTs"^{{VT}}_{{QS}}, " ", "vTt"^{{XR}}_{{WU}}, " ", g^{{"U"}}_{{"T"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"S"}}_{{"R"}}, " ", eta^{{"P"}}_{{"X"}}}], "-", RowBox[{ "vTs"^{{VT}}_{{SQ}}, " ", "vTt"^{{XR}}_{{WU}}, " ", g^{{"U"}}_{{"T"}}, " ", g^{{"W"}}_{{"V"}}, " ", eta^{{"S"}}_{{"R"}}, " ", eta^{{"P"}}_{{"X"}}}]}]], "Output"
```


#### [473] Code For Diag - Secular Approximation
**Input [476]**
```
[{ RowBox[{"RecognizableForm", "[", "dsec", "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"FromCode5", "=", RowBox[{ RowBox[{"SpinLabelTerms", "[", "dsec", "]"}], "//.", {BSpinsBSpint}}]}], "]"}]}], "Input"
```
> I think these can be separated into spin - cases.
**Output [483]**
```
CellGroupData[{ Cell[[{ RowBox[{"ListForm", "[", RowBox[{ RowBox[{ "VTs"^{{Pr}}_{{St}}, " ", "VTt"^{{St}}_{{Pr}}, " ", g^{{"S"}}_{{"S"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"P"}}_{{"Q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{PR}}_{{ST}}, " ", "VTt"^{{ST}}_{{PR}}, " ", g^{{"S"}}_{{"S"}}, " ", g^{{"T"}}_{{"T"}}, " ", eta^{{"P"}}_{{"Q"}}, " ", eta^{{"R"}}_{{"R"}}}], "+", RowBox[{ "VTs"^{{Pr}}_{{Ts}}, " ", "VTt"^{{Ts}}_{{Pr}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"T"}}_{{"T"}}, " ", eta^{{"P"}}_{{"Q"}}, " ", eta^{{"r"}}_{{"r"}}}]}], "]"}], " ", RowBox[{"ListForm", "[", RowBox[{ RowBox[{ "VTs"^{{St}}_{{Qr}}, " ", "VTt"^{{Qr}}_{{St}}, " ", g^{{"S"}}_{{"S"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"P"}}_{{"Q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{Ts}}_{{Qr}}, " ", "VTt"^{{Qr}}_{{Ts}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"T"}}_{{"T"}}, " ", eta^{{"P"}}_{{"Q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{ST}}_{{QR}}, " ", "VTt"^{{QR}}_{{ST}}, " ", g^{{"S"}}_{{"S"}}, " ", g^{{"T"}}_{{"T"}}, " ", eta^{{"P"}}_{{"Q"}}, " ", eta^{{"R"}}_{{"R"}}}]}], "]"}], " ", RowBox[{"ListForm", "[", RowBox[{ RowBox[{ RowBox[{"-", "VTs"^{{Pt}}_{{Rs}}}], " ", "VTt"^{{Rs}}_{{Pt}}, " ", g^{{"P"}}_{{"Q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"R"}}_{{"R"}}, " ", eta^{{"s"}}_{{"s"}}}], "-", RowBox[{ "VTs"^{{Pt}}_{{Sr}}, " ", "VTt"^{{Sr}}_{{Pt}}, " ", g^{{"P"}}_{{"Q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"S"}}_{{"S"}}}], "-", RowBox[{ "VTs"^{{PT}}_{{RS}}, " ", "VTt"^{{RS}}_{{PT}}, " ", g^{{"P"}}_{{"Q"}}, " ", g^{{"T"}}_{{"T"}}, " ", eta^{{"R"}}_{{"R"}}, " ", eta^{{"S"}}_{{"S"}}}]}], "]"}], " ", RowBox[{"ListForm", "[", RowBox[{ RowBox[{ RowBox[{"-", "VTs"^{{Sr}}_{{Qt}}}], " ", "VTt"^{{Qt}}_{{Sr}}, " ", g^{{"P"}}_{{"Q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"S"}}_{{"S"}}}], "-", RowBox[{ "VTs"^{{Rs}}_{{Qt}}, " ", "VTt"^{{Qt}}_{{Rs}}, " ", g^{{"P"}}_{{"Q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"R"}}_{{"R"}}, " ", eta^{{"s"}}_{{"s"}}}], "-", RowBox[{ "VTs"^{{RS}}_{{QT}}, " ", "VTt"^{{QT}}_{{RS}}, " ", g^{{"P"}}_{{"Q"}}, " ", g^{{"T"}}_{{"T"}}, " ", eta^{{"R"}}_{{"R"}}, " ", eta^{{"S"}}_{{"S"}}}]}], "]"}]}], "Input"
```
**Input [484]**
```
[{ RowBox[{"ListForm", "[", RowBox[{ RowBox[{ "VTs"^{{Pr}}_{{St}}, " ", "VTt"^{{St}}_{{Pr}}, " ", g^{{"S"}}_{{"S"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"P"}}_{{"Q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{PR}}_{{ST}}, " ", "VTt"^{{ST}}_{{PR}}, " ", g^{{"S"}}_{{"S"}}, " ", g^{{"T"}}_{{"T"}}, " ", eta^{{"P"}}_{{"Q"}}, " ", eta^{{"R"}}_{{"R"}}}], "+", RowBox[{ "VTs"^{{Pr}}_{{Ts}}, " ", "VTt"^{{Ts}}_{{Pr}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"T"}}_{{"T"}}, " ", eta^{{"P"}}_{{"Q"}}, " ", eta^{{"r"}}_{{"r"}}}]}], "]"}], " ", RowBox[{"ListForm", "[", RowBox[{ RowBox[{ "VTs"^{{St}}_{{Qr}}, " ", "VTt"^{{Qr}}_{{St}}, " ", g^{{"S"}}_{{"S"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"P"}}_{{"Q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{Ts}}_{{Qr}}, " ", "VTt"^{{Qr}}_{{Ts}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"T"}}_{{"T"}}, " ", eta^{{"P"}}_{{"Q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{ST}}_{{QR}}, " ", "VTt"^{{QR}}_{{ST}}, " ", g^{{"S"}}_{{"S"}}, " ", g^{{"T"}}_{{"T"}}, " ", eta^{{"P"}}_{{"Q"}}, " ", eta^{{"R"}}_{{"R"}}}]}], "]"}], " ", RowBox[{"ListForm", "[", RowBox[{ RowBox[{ RowBox[{"-", "VTs"^{{Pt}}_{{Rs}}}], " ", "VTt"^{{Rs}}_{{Pt}}, " ", g^{{"P"}}_{{"Q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"R"}}_{{"R"}}, " ", eta^{{"s"}}_{{"s"}}}], "-", RowBox[{ "VTs"^{{Pt}}_{{Sr}}, " ", "VTt"^{{Sr}}_{{Pt}}, " ", g^{{"P"}}_{{"Q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"S"}}_{{"S"}}}], "-", RowBox[{ "VTs"^{{PT}}_{{RS}}, " ", "VTt"^{{RS}}_{{PT}}, " ", g^{{"P"}}_{{"Q"}}, " ", g^{{"T"}}_{{"T"}}, " ", eta^{{"R"}}_{{"R"}}, " ", eta^{{"S"}}_{{"S"}}}]}], "]"}], " ", RowBox[{"ListForm", "[", RowBox[{ RowBox[{ RowBox[{"-", "VTs"^{{Sr}}_{{Qt}}}], " ", "VTt"^{{Qt}}_{{Sr}}, " ", g^{{"P"}}_{{"Q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"S"}}_{{"S"}}}], "-", RowBox[{ "VTs"^{{Rs}}_{{Qt}}, " ", "VTt"^{{Qt}}_{{Rs}}, " ", g^{{"P"}}_{{"Q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"R"}}_{{"R"}}, " ", eta^{{"s"}}_{{"s"}}}], "-", RowBox[{ "VTs"^{{RS}}_{{QT}}, " ", "VTt"^{{QT}}_{{RS}}, " ", g^{{"P"}}_{{"Q"}}, " ", g^{{"T"}}_{{"T"}}, " ", eta^{{"R"}}_{{"R"}}, " ", eta^{{"S"}}_{{"S"}}}]}], "]"}]}], "Input"
```
**Output [485]**
```
[ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"P"\>\<"r"\>},{\<"S"\>\<"t"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"S"\>\<"t"\>},{\<"P"\>\<"r"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"S"\>}, ",", {\<"S"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"P"\>\<"r"\>},{\<"T"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"T"\>\<"s"\>},{\<"P"\>\<"r"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"T"\>}, ",", {\<"T"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"P"\>\<"R"\>},{\<"S"\>\<"T"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"S"\>\<"T"\>},{\<"P"\>\<"R"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"S"\>}, ",", {\<"S"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"T"\>}, ",", {\<"T"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"R"\>}, ",", {\<"R"\>}}], "}"}]}], "}"}]}], "}"}]], "Output"
```
**Output [486]**
```
[ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"S"\>\<"t"\>},{\<"Q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"Q"\>\<"r"\>},{\<"S"\>\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"S"\>}, ",", {\<"S"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"T"\>\<"s"\>},{\<"Q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"Q"\>\<"r"\>},{\<"T"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"T"\>}, ",", {\<"T"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"S"\>\<"T"\>},{\<"Q"\>\<"R"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"Q"\>\<"R"\>},{\<"S"\>\<"T"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"S"\>}, ",", {\<"S"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"T"\>}, ",", {\<"T"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"R"\>}, ",", {\<"R"\>}}], "}"}]}], "}"}]}], "}"}]], "Output"
```
**Output [487]**
```
[ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"P"\>\<"t"\>},{\<"R"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"R"\>\<"s"\>},{\<"P"\>\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"R"\>}, ",", {\<"R"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"P"\>\<"t"\>},{\<"S"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"S"\>\<"r"\>},{\<"P"\>\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"S"\>}, ",", {\<"S"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"P"\>\<"T"\>},{\<"R"\>\<"S"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"R"\>\<"S"\>},{\<"P"\>\<"T"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"T"\>}, ",", {\<"T"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"R"\>}, ",", {\<"R"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"S"\>}, ",", {\<"S"\>}}], "}"}]}], "}"}]}], "}"}]], "Output"
```
**Output [488]**
```
[ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"R"\>\<"s"\>},{\<"Q"\>\<"t"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"Q"\>\<"t"\>},{\<"R"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"R"\>}, ",", {\<"R"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"S"\>\<"r"\>},{\<"Q"\>\<"t"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"Q"\>\<"t"\>},{\<"S"\>\<"r"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"S"\>}, ",", {\<"S"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"R"\>\<"S"\>},{\<"Q"\>\<"T"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"Q"\>\<"T"\>},{\<"R"\>\<"S"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"T"\>}, ",", {\<"T"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"R"\>}, ",", {\<"R"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"S"\>}, ",", {\<"S"\>}}], "}"}]}], "}"}]}], "}"}]], "Output"
```
**Output [489]**
```
CellGroupData[{ Cell[[ RowBox[{"dspinsec", "=", RowBox[{"ListForm", "[", RowBox[{ RowBox[{ "VTs"^{{Pr}}_{{St}}, " ", "VTt"^{{St}}_{{Pr}}, " ", g^{{"S"}}_{{"S"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"P"}}_{{"Q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{St}}_{{Qr}}, " ", "VTt"^{{Qr}}_{{St}}, " ", g^{{"S"}}_{{"S"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"P"}}_{{"Q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{Pt}}_{{Rs}}, " ", "VTt"^{{Rs}}_{{Pt}}, " ", g^{{"P"}}_{{"Q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"R"}}_{{"R"}}, " ", eta^{{"s"}}_{{"s"}}}], "-", RowBox[{ "VTs"^{{Rs}}_{{Qt}}, " ", "VTt"^{{Qt}}_{{Rs}}, " ", g^{{"P"}}_{{"Q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"R"}}_{{"R"}}, " ", eta^{{"s"}}_{{"s"}}}]}], "]"}]}]], "Input"
```
**Input [490]**
```
[ RowBox[{"dspinsec", "=", RowBox[{"ListForm", "[", RowBox[{ RowBox[{ "VTs"^{{Pr}}_{{St}}, " ", "VTt"^{{St}}_{{Pr}}, " ", g^{{"S"}}_{{"S"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"P"}}_{{"Q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{St}}_{{Qr}}, " ", "VTt"^{{Qr}}_{{St}}, " ", g^{{"S"}}_{{"S"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"P"}}_{{"Q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{Pt}}_{{Rs}}, " ", "VTt"^{{Rs}}_{{Pt}}, " ", g^{{"P"}}_{{"Q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"R"}}_{{"R"}}, " ", eta^{{"s"}}_{{"s"}}}], "-", RowBox[{ "VTs"^{{Rs}}_{{Qt}}, " ", "VTt"^{{Qt}}_{{Rs}}, " ", g^{{"P"}}_{{"Q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"R"}}_{{"R"}}, " ", eta^{{"s"}}_{{"s"}}}]}], "]"}]}]], "Input"
```
**Output [491]**
```
[ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"P"\>\<"r"\>},{\<"S"\>\<"t"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"S"\>\<"t"\>},{\<"P"\>\<"r"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"S"\>}, ",", {\<"S"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"S"\>\<"t"\>},{\<"Q"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"Q"\>\<"r"\>},{\<"S"\>\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"S"\>}, ",", {\<"S"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"P"\>\<"t"\>},{\<"R"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"R"\>\<"s"\>},{\<"P"\>\<"t"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"R"\>}, ",", {\<"R"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"R"\>\<"s"\>},{\<"Q"\>\<"t"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"Q"\>\<"t"\>},{\<"R"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"t"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"R"\>}, ",", {\<"R"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}]}], "}"}]], "Output"
```
**Input [492]**
```
[ RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", {\<VTs\>2},{\<P\>\<r\>},{\<S\>\<t\>}, "}"}], ",", RowBox[{"{", {\<VTt\>2},{\<S\>\<t\>},{\<P\>\<r\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<S\>"}, ",", {"\<S\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<t\>"}, ",", {"\<t\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"\<P\>"}, ",", {"\<Q\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"\<r\>"}, ",", {"\<r\>"}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", {\<VTs\>2},{\<S\>\<t\>},{\<Q\>\<r\>}, "}"}], ",", RowBox[{"{", {\<VTt\>2},{\<Q\>\<r\>},{\<S\>\<t\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<S\>"}, ",", {"\<S\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<t\>"}, ",", {"\<t\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"\<P\>"}, ",", {"\<Q\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"\<r\>"}, ",", {"\<r\>"}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<Sign\>-1}, ",", RowBox[{"{", {\<VTs\>2},{\<P\>\<t\>},{\<R\>\<s\>}, "}"}], ",", RowBox[{"{", {\<VTt\>2},{\<R\>\<s\>},{\<P\>\<t\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<P\>"}, ",", {"\<Q\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<t\>"}, ",", {"\<t\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"\<R\>"}, ",", {"\<R\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"\<s\>"}, ",", {"\<s\>"}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<Sign\>-1}, ",", RowBox[{"{", {\<VTs\>2},{\<R\>\<s\>},{\<Q\>\<t\>}, "}"}], ",", RowBox[{"{", {\<VTt\>2},{\<Q\>\<t\>},{\<R\>\<s\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<P\>"}, ",", {"\<Q\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<t\>"}, ",", {"\<t\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"\<R\>"}, ",", {"\<R\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"\<s\>"}, ",", {"\<s\>"}}], "}"}]}], "}"}]}], "}"}], ";"}]], "Input"
```
> So every term goes to two equivalent OS, and one SS term.
**Input [495]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{"FromCode6", "=", RowBox[{"Map", "[", RelabelAndSpinLabel,dspinsec, "]"}]}], "]"}], " ", AllVts={};, " ", PrintForLoops=True;, " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", TermCodeSec,FromCode6, "]"}], "]"}], " ", RowBox[{"Union", "[", "AllVts", "]"}]}], "Input"
```


#### [500] Secular-two Approximation.
**Input [501]**
```
[ RowBox[{ RowBox[{"SymSecATerm", "[", "ATerm_", "]"}], ":=", RowBox[{"(", " ", RowBox[{"DistributeOperatorProduct", "[", RowBox[{ RowBox[{"DistributeOperatorProduct", "[", RowBox[{ {ATerm}, ",", RowBox[{"STwo", "[", RowBox[{ RowBox[{"ATerm", "[", RowBox[{"[", 2,2,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,2,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 2,3,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,3,1, "]"}], "]"}]}], "]"}]}], "]"}], ",", RowBox[{"STwo", "[", RowBox[{ RowBox[{"ATerm", "[", RowBox[{"[", 2,2,2, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,2,2, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 2,3,2, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,3,2, "]"}], "]"}]}], "]"}]}], "]"}], " ", ")"}]}]], "Input"
```
**Input [503]**
```
[{ RowBox[{ RowBox[{"RecognizableForm", "[", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SymSecATerm,FromCode, "]"}], ",", "1"}], "]"}], "]"}], ";"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"FromCodeSec", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SymSecATerm,FromCode, "]"}], ",", "1"}], "]"}]}], "]"}]}], "]"}]}], "]"}]}], "Input"
```
**Input [506]**
```
[{ RowBox[{"Length", "[", "FromCodeSec", "]"}], " ", RowBox[{ RowBox[{ RowBox[{"TermIndices", "[", "ATerm_", "]"}], ":=", RowBox[{"(", " ", RowBox[{"Map", "[", RowBox[{"ToLowerCase", ",", RowBox[{ RowBox[{"Reap", "[", RowBox[{"Do", "[", " ", RowBox[{ RowBox[{"If", " ", "[", RowBox[{ RowBox[{ RowBox[{ RowBox[{"ATerm", "[", RowBox[{"[", i,1,1, "]"}], "]"}], "=!=", "\"\<Deltas\>\""}], " ", "&&", " ", RowBox[{ RowBox[{"ATerm", "[", RowBox[{"[", i,1,1, "]"}], "]"}], "=!=", "\"\<Deltat\>\""}]}], ",", RowBox[{"Sow", "[", RowBox[{"SlowComplement", "[", RowBox[{ RowBox[{"Select", "[", RowBox[{ RowBox[{"Flatten", "[", RowBox[{"ATerm", "[", RowBox[{"[", "i", "]"}], "]"}], "]"}], ",", RowBox[{ RowBox[{"!", " ", RowBox[{"NumberQ", "[", "#", "]"}]}], "&"}]}], "]"}], ",", RowBox[{"{", RowBox[{ "\"\<Deltas\>\"", ",", "\"\<Deltat\>\"", ",", "\"\<Vs\>\"", ",", "\"\<Vt\>\"", ",", "\"\<VTs\>\"", ",", "\"\<VTt\>\"", ",", "g", ",", "eta", ",", "d", ",", "\"\<Sign\>\"", ",", "a", ",", "\"\<B\>\"", ",", "\"\<M\>\"", ",", "\"\<L\>\""}], "}"}]}], "]"}], "]"}]}], "]"}], " ", ",", RowBox[{"{", RowBox[{"i", ",", "2", ",", RowBox[{"Length", "[", "ATerm", "]"}]}], "}"}]}], "]"}], "]"}], "[", RowBox[{"[", 2,1, "]"}], "]"}]}], "]"}], " ", ")"}]}], ";"}], " ", RowBox[{ RowBox[{ RowBox[{"NumIndicesInTerm", "[", "ATerm_", "]"}], ":=", RowBox[{"(", RowBox[{"Length", "[", RowBox[{"Union", "[", RowBox[{"Flatten", "[", RowBox[{"TermIndices", "[", "ATerm", "]"}], "]"}], "]"}], "]"}], ")"}]}], ";"}], " ", RowBox[{"Map", "[", NumIndicesInTerm,FromCodeSec, "]"}]}], "Input"
```
**Output [509]**
```
CellGroupData[{ Cell["The Inhomogeneity", "Subsection"
```


### [510] The Inhomogeneity
> <The inhomogeneity is the contribution of 2 - body operators to the one body EOM. How to derive that... so you need Q{rho_0} but rho_0 is not a state-varaible... the only way to really define that would be to know Rho_0 after you take the expectation value. I think we only want the first term. So >
**Output [517]**
```
CellGroupData[{ Cell[[{ RowBox[{ RowBox[{"RecognizableForm", "[", com11=r1v1//.{NoCumulants}, "]"}], ";"}], " ", RowBox[{ RowBox[{"RecognizableForm", "[", com12=v1r1//.{NoCumulants}, "]"}], ";"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"com1", "=", RowBox[{"AssociativeSimplify", "[", RowBox[{"Join", "[", RowBox[{"com11", ",", RowBox[{"MultiplyScalar", "[", com12,-1, "]"}]}], "]"}], "]"}]}], "]"}]}], "Input"
```
**Input [518]**
```
[{ RowBox[{ RowBox[{"RecognizableForm", "[", com11=r1v1//.{NoCumulants}, "]"}], ";"}], " ", RowBox[{ RowBox[{"RecognizableForm", "[", com12=v1r1//.{NoCumulants}, "]"}], ";"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"com1", "=", RowBox[{"AssociativeSimplify", "[", RowBox[{"Join", "[", RowBox[{"com11", ",", RowBox[{"MultiplyScalar", "[", com12,-1, "]"}]}], "]"}], "]"}]}], "]"}]}], "Input"
```
**Output [519]**
```
[ RowBox[{ RowBox[{ ta^{{ps}}_{{tu}}, " ", g^{{"r"}}_{{"q"}}}], "-", RowBox[{ ta^{{pr}}_{{tu}}, " ", g^{{"s"}}_{{"q"}}}], "-", RowBox[{ ta^{{rs}}_{{qu}}, " ", g^{{"p"}}_{{"t"}}}], "+", RowBox[{ ta^{{rs}}_{{qt}}, " ", g^{{"p"}}_{{"u"}}}], "+", RowBox[{ ta^{{ps}}_{{tu}}, " ", eta^{{"r"}}_{{"q"}}}], "+", RowBox[{ ta^{{"s"}}_{{"u"}}, " ", g^{{"p"}}_{{"t"}}, " ", eta^{{"r"}}_{{"q"}}}], "-", RowBox[{ ta^{{"s"}}_{{"t"}}, " ", g^{{"p"}}_{{"u"}}, " ", eta^{{"r"}}_{{"q"}}}], "-", RowBox[{ ta^{{pr}}_{{tu}}, " ", eta^{{"s"}}_{{"q"}}}], "-", RowBox[{ ta^{{"r"}}_{{"u"}}, " ", g^{{"p"}}_{{"t"}}, " ", eta^{{"s"}}_{{"q"}}}], "+", RowBox[{ ta^{{"r"}}_{{"t"}}, " ", g^{{"p"}}_{{"u"}}, " ", eta^{{"s"}}_{{"q"}}}], "-", RowBox[{ ta^{{rs}}_{{qu}}, " ", eta^{{"p"}}_{{"t"}}}], "-", RowBox[{ ta^{{"s"}}_{{"u"}}, " ", g^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"t"}}}], "+", RowBox[{ ta^{{"r"}}_{{"u"}}, " ", g^{{"s"}}_{{"q"}}, " ", eta^{{"p"}}_{{"t"}}}], "+", RowBox[{ ta^{{rs}}_{{qt}}, " ", eta^{{"p"}}_{{"u"}}}], "+", RowBox[{ ta^{{"s"}}_{{"t"}}, " ", g^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"u"}}}], "-", RowBox[{ ta^{{"r"}}_{{"t"}}, " ", g^{{"s"}}_{{"q"}}, " ", eta^{{"p"}}_{{"u"}}}]}]], "Output"
```
> This whole computation takes about 20 mins.
**Output [522]**
```
CellGroupData[{ Cell[[{ RowBox[{ RowBox[{"com21", "=", RowBox[{ RowBox[{"AllStepsOfProductRulePProjector", "[", com1,v2, "]"}], "//.", {NoCumulantsNormalVaccuumRule}}]}], ";"}], " ", RowBox[{ RowBox[{"com22", "=", RowBox[{ RowBox[{"AllStepsOfProductRulePProjector", "[", v2,com1, "]"}], "//.", {NoCumulantsNormalVaccuumRule}}]}], ";"}], " ", RowBox[{ RowBox[{"com2", "=", RowBox[{"Join", "[", RowBox[{"com21", ",", RowBox[{"MultiplyScalar", "[", com22,-1, "]"}]}], "]"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", "com2", "]"}]}], "Input"
```
**Input [523]**
```
[{ RowBox[{ RowBox[{"com21", "=", RowBox[{ RowBox[{"AllStepsOfProductRulePProjector", "[", com1,v2, "]"}], "//.", {NoCumulantsNormalVaccuumRule}}]}], ";"}], " ", RowBox[{ RowBox[{"com22", "=", RowBox[{ RowBox[{"AllStepsOfProductRulePProjector", "[", v2,com1, "]"}], "//.", {NoCumulantsNormalVaccuumRule}}]}], ";"}], " ", RowBox[{ RowBox[{"com2", "=", RowBox[{"Join", "[", RowBox[{"com21", ",", RowBox[{"MultiplyScalar", "[", com22,-1, "]"}]}], "]"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", "com2", "]"}]}], "Input"
```
**Output [537]**
```
[ RowBox[{ RowBox[{ g^{{"p"}}_{{"u"}}, " ", g^{{"s"}}_{{"x"}}, " ", g^{{"r"}}_{{"y"}}, " ", eta^{{"w"}}_{{"q"}}, " ", eta^{{"v"}}_{{"t"}}}], "-", RowBox[{ g^{{"p"}}_{{"u"}}, " ", g^{{"r"}}_{{"x"}}, " ", g^{{"s"}}_{{"y"}}, " ", eta^{{"w"}}_{{"q"}}, " ", eta^{{"v"}}_{{"t"}}}], "-", RowBox[{ g^{{"p"}}_{{"u"}}, " ", g^{{"s"}}_{{"x"}}, " ", g^{{"r"}}_{{"y"}}, " ", eta^{{"v"}}_{{"q"}}, " ", eta^{{"w"}}_{{"t"}}}], "+", RowBox[{ g^{{"p"}}_{{"u"}}, " ", g^{{"r"}}_{{"x"}}, " ", g^{{"s"}}_{{"y"}}, " ", eta^{{"v"}}_{{"q"}}, " ", eta^{{"w"}}_{{"t"}}}], "+", RowBox[{ g^{{"s"}}_{{"x"}}, " ", g^{{"r"}}_{{"y"}}, " ", eta^{{"w"}}_{{"q"}}, " ", eta^{{"v"}}_{{"t"}}, " ", eta^{{"p"}}_{{"u"}}}], "-", RowBox[{ g^{{"r"}}_{{"x"}}, " ", g^{{"s"}}_{{"y"}}, " ", eta^{{"w"}}_{{"q"}}, " ", eta^{{"v"}}_{{"t"}}, " ", eta^{{"p"}}_{{"u"}}}], "-", RowBox[{ g^{{"s"}}_{{"x"}}, " ", g^{{"r"}}_{{"y"}}, " ", eta^{{"v"}}_{{"q"}}, " ", eta^{{"w"}}_{{"t"}}, " ", eta^{{"p"}}_{{"u"}}}], "+", RowBox[{ g^{{"r"}}_{{"x"}}, " ", g^{{"s"}}_{{"y"}}, " ", eta^{{"v"}}_{{"q"}}, " ", eta^{{"w"}}_{{"t"}}, " ", eta^{{"p"}}_{{"u"}}}], "-", RowBox[{ g^{{"p"}}_{{"t"}}, " ", g^{{"s"}}_{{"x"}}, " ", g^{{"r"}}_{{"y"}}, " ", eta^{{"w"}}_{{"q"}}, " ", eta^{{"v"}}_{{"u"}}}], "+", RowBox[{ g^{{"p"}}_{{"t"}}, " ", g^{{"r"}}_{{"x"}}, " ", g^{{"s"}}_{{"y"}}, " ", eta^{{"w"}}_{{"q"}}, " ", eta^{{"v"}}_{{"u"}}}], "-", RowBox[{ g^{{"s"}}_{{"x"}}, " ", g^{{"r"}}_{{"y"}}, " ", eta^{{"w"}}_{{"q"}}, " ", eta^{{"p"}}_{{"t"}}, " ", eta^{{"v"}}_{{"u"}}}], "+", RowBox[{ g^{{"r"}}_{{"x"}}, " ", g^{{"s"}}_{{"y"}}, " ", eta^{{"w"}}_{{"q"}}, " ", eta^{{"p"}}_{{"t"}}, " ", eta^{{"v"}}_{{"u"}}}], "-", RowBox[{ g^{{"s"}}_{{"q"}}, " ", g^{{"r"}}_{{"x"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"w"}}_{{"t"}}, " ", eta^{{"v"}}_{{"u"}}}], "+", RowBox[{ g^{{"r"}}_{{"q"}}, " ", g^{{"s"}}_{{"x"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"w"}}_{{"t"}}, " ", eta^{{"v"}}_{{"u"}}}], "+", RowBox[{ g^{{"s"}}_{{"q"}}, " ", g^{{"p"}}_{{"x"}}, " ", g^{{"r"}}_{{"y"}}, " ", eta^{{"w"}}_{{"t"}}, " ", eta^{{"v"}}_{{"u"}}}], "-", RowBox[{ g^{{"r"}}_{{"q"}}, " ", g^{{"p"}}_{{"x"}}, " ", g^{{"s"}}_{{"y"}}, " ", eta^{{"w"}}_{{"t"}}, " ", eta^{{"v"}}_{{"u"}}}], "+", RowBox[{ g^{{"s"}}_{{"x"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"w"}}_{{"t"}}, " ", eta^{{"v"}}_{{"u"}}}], "-", RowBox[{ g^{{"p"}}_{{"x"}}, " ", g^{{"s"}}_{{"y"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"w"}}_{{"t"}}, " ", eta^{{"v"}}_{{"u"}}}], "-", RowBox[{ g^{{"r"}}_{{"x"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"q"}}, " ", eta^{{"w"}}_{{"t"}}, " ", eta^{{"v"}}_{{"u"}}}], "+", RowBox[{ g^{{"p"}}_{{"x"}}, " ", g^{{"r"}}_{{"y"}}, " ", eta^{{"s"}}_{{"q"}}, " ", eta^{{"w"}}_{{"t"}}, " ", eta^{{"v"}}_{{"u"}}}], "+", RowBox[{ g^{{"p"}}_{{"t"}}, " ", g^{{"s"}}_{{"x"}}, " ", g^{{"r"}}_{{"y"}}, " ", eta^{{"v"}}_{{"q"}}, " ", eta^{{"w"}}_{{"u"}}}], "-", RowBox[{ g^{{"p"}}_{{"t"}}, " ", g^{{"r"}}_{{"x"}}, " ", g^{{"s"}}_{{"y"}}, " ", eta^{{"v"}}_{{"q"}}, " ", eta^{{"w"}}_{{"u"}}}], "+", RowBox[{ g^{{"s"}}_{{"x"}}, " ", g^{{"r"}}_{{"y"}}, " ", eta^{{"v"}}_{{"q"}}, " ", eta^{{"p"}}_{{"t"}}, " ", eta^{{"w"}}_{{"u"}}}], "-", RowBox[{ g^{{"r"}}_{{"x"}}, " ", g^{{"s"}}_{{"y"}}, " ", eta^{{"v"}}_{{"q"}}, " ", eta^{{"p"}}_{{"t"}}, " ", eta^{{"w"}}_{{"u"}}}], "+", RowBox[{ g^{{"s"}}_{{"q"}}, " ", g^{{"r"}}_{{"x"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"v"}}_{{"t"}}, " ", eta^{{"w"}}_{{"u"}}}], "-", RowBox[{ g^{{"r"}}_{{"q"}}, " ", g^{{"s"}}_{{"x"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"v"}}_{{"t"}}, " ", eta^{{"w"}}_{{"u"}}}], "-", RowBox[{ g^{{"s"}}_{{"q"}}, " ", g^{{"p"}}_{{"x"}}, " ", g^{{"r"}}_{{"y"}}, " ", eta^{{"v"}}_{{"t"}}, " ", eta^{{"w"}}_{{"u"}}}], "+", RowBox[{ g^{{"r"}}_{{"q"}}, " ", g^{{"p"}}_{{"x"}}, " ", g^{{"s"}}_{{"y"}}, " ", eta^{{"v"}}_{{"t"}}, " ", eta^{{"w"}}_{{"u"}}}], "-", RowBox[{ g^{{"s"}}_{{"x"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"v"}}_{{"t"}}, " ", eta^{{"w"}}_{{"u"}}}], "+", RowBox[{ g^{{"p"}}_{{"x"}}, " ", g^{{"s"}}_{{"y"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"v"}}_{{"t"}}, " ", eta^{{"w"}}_{{"u"}}}], "+", RowBox[{ g^{{"r"}}_{{"x"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"q"}}, " ", eta^{{"v"}}_{{"t"}}, " ", eta^{{"w"}}_{{"u"}}}], "-", RowBox[{ g^{{"p"}}_{{"x"}}, " ", g^{{"r"}}_{{"y"}}, " ", eta^{{"s"}}_{{"q"}}, " ", eta^{{"v"}}_{{"t"}}, " ", eta^{{"w"}}_{{"u"}}}], "+", RowBox[{ g^{{"s"}}_{{"q"}}, " ", g^{{"w"}}_{{"t"}}, " ", g^{{"v"}}_{{"u"}}, " ", eta^{{"r"}}_{{"x"}}, " ", eta^{{"p"}}_{{"y"}}}], "-", RowBox[{ g^{{"s"}}_{{"q"}}, " ", g^{{"v"}}_{{"t"}}, " ", g^{{"w"}}_{{"u"}}, " ", eta^{{"r"}}_{{"x"}}, " ", eta^{{"p"}}_{{"y"}}}], "+", RowBox[{ g^{{"w"}}_{{"t"}}, " ", g^{{"v"}}_{{"u"}}, " ", eta^{{"s"}}_{{"q"}}, " ", eta^{{"r"}}_{{"x"}}, " ", eta^{{"p"}}_{{"y"}}}], "-", RowBox[{ g^{{"v"}}_{{"t"}}, " ", g^{{"w"}}_{{"u"}}, " ", eta^{{"s"}}_{{"q"}}, " ", eta^{{"r"}}_{{"x"}}, " ", eta^{{"p"}}_{{"y"}}}], "-", RowBox[{ g^{{"r"}}_{{"q"}}, " ", g^{{"w"}}_{{"t"}}, " ", g^{{"v"}}_{{"u"}}, " ", eta^{{"s"}}_{{"x"}}, " ", eta^{{"p"}}_{{"y"}}}], "+", RowBox[{ g^{{"r"}}_{{"q"}}, " ", g^{{"v"}}_{{"t"}}, " ", g^{{"w"}}_{{"u"}}, " ", eta^{{"s"}}_{{"x"}}, " ", eta^{{"p"}}_{{"y"}}}], "-", RowBox[{ g^{{"w"}}_{{"t"}}, " ", g^{{"v"}}_{{"u"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"s"}}_{{"x"}}, " ", eta^{{"p"}}_{{"y"}}}], "+", RowBox[{ g^{{"v"}}_{{"t"}}, " ", g^{{"w"}}_{{"u"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"s"}}_{{"x"}}, " ", eta^{{"p"}}_{{"y"}}}], "-", RowBox[{ g^{{"s"}}_{{"q"}}, " ", g^{{"w"}}_{{"t"}}, " ", g^{{"v"}}_{{"u"}}, " ", eta^{{"p"}}_{{"x"}}, " ", eta^{{"r"}}_{{"y"}}}], "+", RowBox[{ g^{{"s"}}_{{"q"}}, " ", g^{{"v"}}_{{"t"}}, " ", g^{{"w"}}_{{"u"}}, " ", eta^{{"p"}}_{{"x"}}, " ", eta^{{"r"}}_{{"y"}}}], "-", RowBox[{ g^{{"w"}}_{{"t"}}, " ", g^{{"v"}}_{{"u"}}, " ", eta^{{"s"}}_{{"q"}}, " ", eta^{{"p"}}_{{"x"}}, " ", eta^{{"r"}}_{{"y"}}}], "+", RowBox[{ g^{{"v"}}_{{"t"}}, " ", g^{{"w"}}_{{"u"}}, " ", eta^{{"s"}}_{{"q"}}, " ", eta^{{"p"}}_{{"x"}}, " ", eta^{{"r"}}_{{"y"}}}], "-", RowBo ...[TRUNC 8721]
```
**Input [538]**
```
CellGroupData[{ Cell[[ RowBox[{"RecognizableForm", "[", RowBox[{"com2", "-", RowBox[{"(", RowBox[{"com2", "/.", RowBox[{"{", RowBox[{ RowBox[{"\"\<p\>\"", "->", "\"\<q\>\""}], ",", RowBox[{"\"\<q\>\"", "->", "\"\<p\>\""}]}], "}"}]}], ")"}]}], "]"}]], "Input"], Cell[["0"], "Output"
```
**Input [539]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"com2", "-", RowBox[{"(", RowBox[{"com2", "/.", RowBox[{"{", RowBox[{ RowBox[{"\"\<p\>\"", "->", "\"\<q\>\""}], ",", RowBox[{"\"\<q\>\"", "->", "\"\<p\>\""}]}], "}"}]}], ")"}]}], "]"}]], "Input"
```
**Output [540]**
```
["0"], "Output"
```
**Input [544]**
```
[{ RowBox[{ RowBox[{"TCLNoDelta1sI", "=", RowBox[{"AssociativeSimplify", "[", RowBox[{"Map", "[", CollapseDeltas,WithV1andV2I, "]"}], "]"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", "%", "]"}], " ", RowBox[{"Length", "[", "TCLNoDelta1sI", "]"}]}], "Input"
```
**Output [547]**
```
CellGroupData[{ Cell[[ RowBox[{"RecognizableForm", "[", RowBox[{"FinalTCL2", "=", RowBox[{"AssociativeSimplify", "[", RowBox[{"MultiplyScalar", "[", RowBox[{ RowBox[{"Map", "[", RelabelTerm,TCLNoDelta1sI, "]"}], ",", -1/16}], "]"}], "]"}]}], "]"}]], "Input"
```
**Output [549]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{ry}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{wx}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{vy}}_{{su}}, " ", "VTt"^{{rt}}_{{wx}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{wx}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{tv}}_{{sy}}, " ", "VTt"^{{rx}}_{{uw}}, " ", g^{{"y"}}_{{"q"}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"x"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{vx}}_{{su}}, " ", "VTt"^{{rt}}_{{wy}}, " ", g^{{"y"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}, " ", eta^{{"p"}}_{{"x"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{ry}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"y"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{tv}}_{{sx}}, " ", "VTt"^{{ry}}_{{uw}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"y"}}}]}]], "Output"
```
**Input [550]**
```
[""], "Input"
```


### [552] Resolution of the Identities & Armadillo AutoCoding.
**Input [555]**
```
[ RowBox[{ (*Goal here is to get a unique re-labeling of a termby ordering the loops of the contractions in a well-defined way and then relabeling those loops*), " ", RowBox[{ RowBox[{ RowBox[{"GetIndexLoops", "[", "Indices_", "]"}], ":=", RowBox[{"(", " ", RowBox[{ RowBox[{"PPos", "=", RowBox[{ RowBox[{"Position", "[", Indices,"\<p\>", "]"}], "[", RowBox[{"[", 1,1, "]"}], "]"}]}], ";", " ", RowBox[{"Working", "=", RowBox[{"Complement", "[", RowBox[{ RowBox[{"Range", "[", RowBox[{"Length", "[", "Indices", "]"}], "]"}], ",", {PPos}}], "]"}]}], ";", " ", Loop={PPos}, ";", " ", RowBox[{"LoopIndices", "=", RowBox[{"Indices", "[", RowBox[{"[", "PPos", "]"}], "]"}]}], ";", " ", RowBox[{"Return", "[", RowBox[{ RowBox[{"Reap", "[", RowBox[{"Do", "[", RowBox[{ RowBox[{"(", " ", RowBox[{ FoundMatch=False, ";", " ", RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{"Length", "[", "Working", "]"}], "<", "1"}], ",", RowBox[{"Return", "[", "]"}]}], "]"}], ";", " ", RowBox[{"Do", "[", RowBox[{ RowBox[{"(", " ", RowBox[{ RowBox[{"If", "[", RowBox[{ RowBox[{"j", ">", RowBox[{"Length", "[", "Working", "]"}]}], ",", RowBox[{"Continue", "[", "]"}]}], "]"}], ";", " ", RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{"Length", "[", RowBox[{"Intersection", "[", RowBox[{ RowBox[{"Indices", "[", RowBox[{"[", RowBox[{"Working", "[", RowBox[{"[", "j", "]"}], "]"}], "]"}], "]"}], ",", "LoopIndices"}], "]"}], "]"}], ">", "0"}], ",", RowBox[{"(", " ", RowBox[{ RowBox[{"Loop", "=", RowBox[{"Append", "[", RowBox[{"Loop", ",", RowBox[{"Working", "[", RowBox[{"[", "j", "]"}], "]"}]}], "]"}]}], ";", " ", RowBox[{"Working", "=", RowBox[{"Complement", "[", RowBox[{"Working", ",", RowBox[{"{", RowBox[{"Working", "[", RowBox[{"[", "j", "]"}], "]"}], "}"}]}], "]"}]}], ";", " ", RowBox[{"LoopIndices", "=", RowBox[{"Flatten", "[", RowBox[{"Part", "[", Indices,Loop, "]"}], "]"}]}], ";", " ", FoundMatch=True, ";", " ", RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{"Length", "[", "Working", "]"}], "==", "0"}], ",", RowBox[{"Sow", "[", "Loop", "]"}]}], "]"}], ";"}], " ", ")"}]}], "]"}], ";"}], " ", ")"}], ",", RowBox[{"{", RowBox[{"j", ",", "1", ",", RowBox[{"Length", "[", "Working", "]"}]}], "}"}]}], "]"}], ";", " ", RowBox[{"If", "[", RowBox[{ !FoundMatch, ",", RowBox[{"(", " ", (*This loop is done.*), " ", RowBox[{ RowBox[{"Sow", "[", "Loop", "]"}], ";", " ", RowBox[{"Loop", "=", RowBox[{"{", RowBox[{"Working", "[", RowBox[{"[", "1", "]"}], "]"}], "}"}]}], ";", " ", RowBox[{"Working", "=", RowBox[{"Complement", "[", RowBox[{"Working", ",", RowBox[{"{", RowBox[{"Working", "[", RowBox[{"[", "1", "]"}], "]"}], "}"}]}], "]"}]}], ";", " ", RowBox[{"LoopIndices", "=", RowBox[{"Flatten", "[", RowBox[{"Part", "[", Indices,Loop, "]"}], "]"}]}], ";"}], " ", ")"}]}], "]"}], ";"}], " ", ")"}], ",", RowBox[{"{", RowBox[{"i", ",", "1", ",", RowBox[{"Length", "[", "Working", "]"}]}], "}"}]}], "]"}], "]"}], "[", RowBox[{"[", 2,1, "]"}], "]"}], "]"}]}], " ", ")"}]}], ";"}]}]], "Input"
```


### [565] Secular From Code.
**Input [567]**
```
[ RowBox[{ RowBox[{"SymSecATerm", "[", "ATerm_", "]"}], ":=", RowBox[{"(", " ", RowBox[{"DistributeOperatorProduct", "[", RowBox[{ RowBox[{"DistributeOperatorProduct", "[", RowBox[{ {ATerm}, ",", RowBox[{"STwo", "[", RowBox[{ RowBox[{"ATerm", "[", RowBox[{"[", 2,2,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,2,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 2,3,1, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,3,1, "]"}], "]"}]}], "]"}]}], "]"}], ",", RowBox[{"STwo", "[", RowBox[{ RowBox[{"ATerm", "[", RowBox[{"[", 2,2,2, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,2,2, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 2,3,2, "]"}], "]"}], ",", RowBox[{"ATerm", "[", RowBox[{"[", 3,3,2, "]"}], "]"}]}], "]"}]}], "]"}], " ", ")"}]}]], "Input"
```
**Output [568]**
```
CellGroupData[{ Cell[[ RowBox[{"RecognizableForm", "[", "FromCode", "]"}]], "Input"
```
**Input [569]**
```
[ RowBox[{"RecognizableForm", "[", "FromCode", "]"}]], "Input"
```
**Output [570]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{ry}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"y"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{wx}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", d^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{vy}}_{{su}}, " ", "VTt"^{{rt}}_{{wx}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{tv}}_{{sy}}, " ", "VTt"^{{rx}}_{{uw}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"y"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"x"}}}]}]], "Output"
```
**Input [572]**
```
[{ SubstituteDeltas=True;, " ", RowBox[{"RecognizableForm", "[", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SymSecATerm,FromCode, "]"}], ",", "1"}], "]"}], "]"}]}], "Input"
```
**Input [575]**
```
[{ SubstituteDeltas=False;, " ", RowBox[{"RecognizableForm", "[", RowBox[{"FromCodeSec", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SymSecATerm,FromCode, "]"}], ",", "1"}], "]"}]}], "]"}]}], "]"}]}], "]"}]}], "Input"
```
**Input [577]**
```
[{ RowBox[{ RowBox[{"NoDiagGamma1", "=", RowBox[{ RowBox[{"{", RowBox[{"any1___", ",", RowBox[{"{", RowBox[{ {\<Sign\>a_}, ",", "any5___", ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", "p_", ",", "p_"}], "}"}], ",", "any8___"}], "}"}], ",", "any4___"}], "}"}], "->", {any1any4}}]}], ";"}], " ", RowBox[{ RowBox[{"NoDiagGamma11", "=", RowBox[{ RowBox[{"{", RowBox[{"any1___", ",", RowBox[{"{", RowBox[{ {\<Sign\>a_}, ",", "any5___", ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", "p_", ",", "p_"}], "}"}], ",", "any8___"}], "}"}], ",", "any4___"}], "}"}], "->", {any1any4}}]}], ";"}]}], "Input"
```
**Input [579]**
```
[ FromCodeSec//.{NoDiagGamma1NoDiagGamma11}], "Input"
```
**Input [582]**
```
[{ RowBox[{"Length", "[", "FromCodeSec", "]"}], " ", RowBox[{"Map", "[", NumIndicesInTerm,FromCodeSec, "]"}]}], "Input"
```
**Input [586]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"TCLSimp4", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", DiagonalGamma1,FromCode, "]"}]}], "]"}]}], "]"}]}], "]"}]], "Input"
```
**Input [589]**
```
["TCLSimp4"], "Input"
```
**Input [592]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{"FromCode5", "=", RowBox[{ RowBox[{"SpinLabelTerms", "[", "TCLSimp4", "]"}], "//.", {BSpinsBSpint}}]}], "]"}], " ", RowBox[{ RowBox[{"RecognizableForm", "[", RowBox[{"FromCode6", "=", RowBox[{"Map", "[", RelabelAndSpinLabel,FromCode5, "]"}]}], "]"}], ";"}], " ", AllVts={};, " ", PrintForLoops=True;, " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", TermCodeSec,FromCode6, "]"}], "]"}], " ", RowBox[{"Union", "[", "AllVts", "]"}]}], "Input"
```
**Input [595]**
```
[{ RowBox[{ RowBox[{"RecognizableForm", "[", RowBox[{"FromCode5", "=", RowBox[{ RowBox[{"SpinLabelTerms", "[", "FromCodeSec", "]"}], "//.", {BSpinsBSpint}}]}], "]"}], ";"}], " ", RowBox[{ RowBox[{"RecognizableForm", "[", RowBox[{"FromCode6", "=", RowBox[{"Map", "[", RelabelAndSpinLabel,FromCode5, "]"}]}], "]"}], ";", AllVts={}, ";"}], " ", PrintForLoops=False;, " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", TermCodeSec,FromCode6, "]"}], "]"}], " ", RowBox[{"Union", "[", "AllVts", "]"}]}], "Input"
```
**Output [597]**
```
CellGroupData[{ Cell["Correspondence with CIS(D)", "Subsection"
```


### [598] Correspondence with CIS(D)
**Output [599]**
```
CellGroupData[{ Cell[[ RowBox[{"RecognizableForm", "[", "FinalTCL2", "]"}]], "Input"
```
**Input [600]**
```
[ RowBox[{"RecognizableForm", "[", "FinalTCL2", "]"}]], "Input"
```
**Output [601]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{ry}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{wx}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{vy}}_{{su}}, " ", "VTt"^{{rt}}_{{wx}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{wx}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{tv}}_{{sy}}, " ", "VTt"^{{rx}}_{{uw}}, " ", g^{{"y"}}_{{"q"}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"x"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{vx}}_{{su}}, " ", "VTt"^{{rt}}_{{wy}}, " ", g^{{"y"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}, " ", eta^{{"p"}}_{{"x"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{ry}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"y"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{tv}}_{{sx}}, " ", "VTt"^{{ry}}_{{uw}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"y"}}}]}]], "Output"
```
> <The rule that I've worked out on paper shows that the two diagrams correspond to the two terms in CIS (D). The MP2 like term < [Phi] | V | T2U1 > come from the terms where the quasiparticle character of both lines on PQ are the same. The Excited state correlation term < [Phi] | V | U2 > comes about from the diagram where the quasiparticle character on the pq lines are different. The denominators on these terms are both O2V2 in character however the ESC term has an omega in the denominator too. >
**Input [603]**
```
[{ RowBox[{ RowBox[{ RowBox[{"DiagA", "=", RowBox[{"ListForm", "[", RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{ry}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{wx}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{tv}}_{{sy}}, " ", "VTt"^{{rx}}_{{uw}}, " ", g^{{"y"}}_{{"q"}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"x"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{vx}}_{{su}}, " ", "VTt"^{{rt}}_{{wy}}, " ", g^{{"y"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}, " ", eta^{{"p"}}_{{"x"}}}]}], "]"}]}], ";"}], (*These are ESC-Like and should be negative for rho_ai*)}], " ", RowBox[{ RowBox[{ RowBox[{"DiagB", "=", RowBox[{"ListForm", "[", RowBox[{ RowBox[{ -(1)/(2), " ", "VTs"^{{rt}}_{{wx}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{vy}}_{{su}}, " ", "VTt"^{{rt}}_{{wx}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{ry}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"y"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{tv}}_{{sx}}, " ", "VTt"^{{ry}}_{{uw}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"y"}}}]}], "]"}]}], ";"}], (*These are MP2-like*)}]}], "Input"
```
> <To obtain CIS(d) from these diagrams, first let[CloseCurlyQuote]s consider the MP2 like terms, The rule for these terms are that every line is diagonal except for one matrix coming off of pq. This is the same rule obeyed for the ESC terms. >
**Output [605]**
```
CellGroupData[{ Cell[[{ RowBox[{ RowBox[{ RowBox[{"DiagonalDeltasNoP", "[", "aterm_", "]"}], ":=", RowBox[{ RowBox[{"Reap", "[", RowBox[{"Map", "[", RowBox[{ RowBox[{ RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{"(", RowBox[{ RowBox[{ RowBox[{"#", "[", RowBox[{"[", "1", "]"}], "]"}], "===", RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}]}], " ", "||", " ", RowBox[{ RowBox[{"#", "[", RowBox[{"[", "1", "]"}], "]"}], "===", RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}]}]}], ")"}], "&&", RowBox[{"(", RowBox[{"!", RowBox[{"MemberQ", "[", RowBox[{ RowBox[{"Flatten", "[", "#", "]"}], ",", "\"\<p\>\""}], "]"}]}], ")"}]}], ",", RowBox[{"Sow", "[", RowBox[{"ADeltaForFactor", "[", "#", "]"}], "]"}], ","}], "]"}], "&"}], ",", "aterm"}], "]"}], "]"}], "[", RowBox[{"[", 2,1, "]"}], "]"}]}], ";"}], " ", RowBox[{ RowBox[{ RowBox[{"DiagonalGammaNoPq", "[", "aterm_", "]"}], ":=", RowBox[{"Join", "[", RowBox[{"aterm", ",", RowBox[{"DiagonalDeltasNoPq", "[", "aterm", "]"}]}], "]"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"dDiagA", "=", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", DiagonalGammaNoPq,DiagA, "]"}]}], "]"}]}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"dDiagB", "=", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", DiagonalGammaNoPq,DiagB, "]"}]}], "]"}]}], "]"}]}], "Input"
```
**Output [611]**
```
[ RowBox[{ RowBox[{ RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"DiagonalDeltasNoPq", "[", "]"}], ",", RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"t"\>},{\<"w"\>\<"x"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"v"\>\<"y"\>},{\<"s"\>\<"u"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"x"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}]}], "}"}], "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"t"\>},{\<"w"\>\<"x"\>}, "}"}]}], "}"}], "}"}]}], "]"}], "[[", 1,2, "]]"}], " ", TemplateBox[{"\"Sign\"",RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{"{", {VTs2},{rt},{wx}, "}"}]}], "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{"{", {VTs2},{rt},{wx}, "}"}]}], "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"], " ", TemplateBox[{{Sign-(1)/(2)},RowBox[{"{", RowBox[{ RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {Sign-(1)/(2)}, ",", RowBox[{"{", {VTs2},{rt},{wx}, "}"}], ",", RowBox[{"{", {VTt2},{vy},{su}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"w"}, ",", {"v"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"p"}, ",", {"y"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"x"}, ",", {"q"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"s"}, ",", {"r"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"u"}, ",", {"t"}}], "}"}]}], "}"}], "]"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {Sign-(1)/(2)}, ",", RowBox[{"{", {VTs2},{rt},{wx}, "}"}], ",", RowBox[{"{", {VTt2},{vy},{su}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"w"}, ",", {"v"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"p"}, ",", {"y"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"x"}, ",", {"q"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"s"}, ",", {"r"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"u"}, ",", {"t"}}], "}"}]}], "}"}], "]"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], "+", RowBox[{ RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"DiagonalDeltasNoPq", "[", "]"}], ",", RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"v"\>\<"x"\>},{\<"s"\>\<"u"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"r"\>\<"t"\>},{\<"w"\>\<"y"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"y"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"x"\>}}], "}"}]}], "}"}], "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTs"\>2},{\<"v"\>\<"x"\>},{\<"s"\>\<"u"\>}, "}"}]}], "}"}], "}"}]}], "]"}], "[[", 1,2, "]]"}], " ", TemplateBox[{"\"Sign\"",RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{"{", {VTs2},{vx},{su}, "}"}]}], "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{"{", {VTs2},{vx},{su}, "}"}]}], "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"], " ", TemplateBox[{{Sign-(1)/(2)},RowBox[{"{", RowBox[{ RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {Sign-(1)/(2)}, ",", RowBox[{"{", {VTs2},{vx},{su}, "}"}], ",", RowBox[{"{", {VTt2},{rt},{wy}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"y"}, ",", {"q"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"w"}, ",", {"v"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"s"}, ",", {"r"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"u"}, ",", {"t"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"p"}, ",", {"x"}}], "}"}]}], "}"}], "]"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {Sign-(1)/(2)}, ",", RowBox[{"{", {VTs2},{vx},{su}, "}"}], ",", RowBox[{"{", {VTt2},{rt},{wy}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"y"}, ",", {"q"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"w"}, ",", {"v"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"s"}, ",", {"r"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"u"}, ",", {"t"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"p"}, ",", {"x"}}], "}"}]}], "}"}], "]"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], "+", RowBox[{ RowBox[{ RowBox[{"Join", "[", RowBox ...[TRUNC 11885]
```
**Output [616]**
```
[ RowBox[{ RowBox[{ RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"DiagonalDeltasNoPq", "[", "]"}], ",", RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"t"\>},{\<"w"\>\<"x"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"v"\>\<"y"\>},{\<"s"\>\<"u"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"x"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}]}], "}"}], "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"t"\>},{\<"w"\>\<"x"\>}, "}"}]}], "}"}], "}"}]}], "]"}], "[[", 1,2, "]]"}], " ", TemplateBox[{"\"Sign\"",RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{"{", {VTs2},{rt},{wx}, "}"}]}], "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{"{", {VTs2},{rt},{wx}, "}"}]}], "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"], " ", TemplateBox[{{Sign-(1)/(2)},RowBox[{"{", RowBox[{ RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {Sign-(1)/(2)}, ",", RowBox[{"{", {VTs2},{rt},{wx}, "}"}], ",", RowBox[{"{", {VTt2},{vy},{su}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"x"}, ",", {"q"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"w"}, ",", {"v"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"p"}, ",", {"y"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"s"}, ",", {"r"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"u"}, ",", {"t"}}], "}"}]}], "}"}], "]"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {Sign-(1)/(2)}, ",", RowBox[{"{", {VTs2},{rt},{wx}, "}"}], ",", RowBox[{"{", {VTt2},{vy},{su}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"x"}, ",", {"q"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"w"}, ",", {"v"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"p"}, ",", {"y"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"s"}, ",", {"r"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"u"}, ",", {"t"}}], "}"}]}], "}"}], "]"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], "+", RowBox[{ RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"DiagonalDeltasNoPq", "[", "]"}], ",", RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"v"\>\<"y"\>},{\<"s"\>\<"u"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"r"\>\<"t"\>},{\<"w"\>\<"x"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"x"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}]}], "}"}], "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTs"\>2},{\<"v"\>\<"y"\>},{\<"s"\>\<"u"\>}, "}"}]}], "}"}], "}"}]}], "]"}], "[[", 1,2, "]]"}], " ", TemplateBox[{"\"Sign\"",RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{"{", {VTs2},{vy},{su}, "}"}]}], "}"}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{"{", {VTs2},{vy},{su}, "}"}]}], "}"}], "}"}], ",", "2"}], "}"}]}, "Subsuperscript"], " ", TemplateBox[{{Sign-(1)/(2)},RowBox[{"{", RowBox[{ RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {Sign-(1)/(2)}, ",", RowBox[{"{", {VTs2},{vy},{su}, "}"}], ",", RowBox[{"{", {VTt2},{rt},{wx}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"x"}, ",", {"q"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"w"}, ",", {"v"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"p"}, ",", {"y"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"s"}, ",", {"r"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"u"}, ",", {"t"}}], "}"}]}], "}"}], "]"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{ RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {Sign-(1)/(2)}, ",", RowBox[{"{", {VTs2},{vy},{su}, "}"}], ",", RowBox[{"{", {VTt2},{rt},{wx}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"x"}, ",", {"q"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"w"}, ",", {"v"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"p"}, ",", {"y"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"s"}, ",", {"r"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {"u"}, ",", {"t"}}], "}"}]}], "}"}], "]"}], ",", "2"}], "}"}]}, "Subsuperscript"]}], "+", RowBox[{ RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{" ...[TRUNC 11885]
```
**Input [617]**
```
[ RowBox[{ RowBox[{"IntegralsForTerm", "[", "ATerm_", "]"}], ":=", RowBox[{"(", " ", RowBox[{ RowBox[{"VS", "=", RowBox[{"ATerm", "[", RowBox[{"[", "2", "]"}], "]"}]}], ";", " ", RowBox[{"VT", "=", RowBox[{"ATerm", "[", RowBox[{"[", "3", "]"}], "]"}]}], ";", " ", RowBox[{"samespin", "=", RowBox[{ RowBox[{"UpperCaseQ", "[", RowBox[{"VS", "[", RowBox[{"[", 2,1, "]"}], "]"}], "]"}], "&&", RowBox[{"UpperCaseQ", "[", RowBox[{"VS", "[", RowBox[{"[", 2,2, "]"}], "]"}], "]"}]}]}], ";", " ", RowBox[{"Return", "[", RowBox[{"If", "[", RowBox[{"samespin", ",", RowBox[{ RowBox[{"\"\<Vss[n3*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VS", "[", RowBox[{"[", 2,1, "]"}], "]"}], "]"}], "<>", "\"\<+n2*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VS", "[", RowBox[{"[", 2,2, "]"}], "]"}], "]"}], "<>", "\"\<+n*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VS", "[", RowBox[{"[", 3,1, "]"}], "]"}], "]"}], "<>", "\"\<+\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VS", "[", RowBox[{"[", 3,2, "]"}], "]"}], "]"}]}], RowBox[{"\"\<]*Vts[\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VT", "[", RowBox[{"[", 2,1, "]"}], "]"}], "]"}], "<>", "\"\<+n2*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VT", "[", RowBox[{"[", 2,2, "]"}], "]"}], "]"}], "<>", "\"\<+n*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VT", "[", RowBox[{"[", 3,1, "]"}], "]"}], "]"}], "<>", "\"\<+\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VT", "[", RowBox[{"[", 3,2, "]"}], "]"}], "]"}], "<>", "\"\<]\>\""}]}], ",", RowBox[{ RowBox[{"\"\<Vso[n3*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VS", "[", RowBox[{"[", 2,1, "]"}], "]"}], "]"}], "<>", "\"\<+n2*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VS", "[", RowBox[{"[", 2,2, "]"}], "]"}], "]"}], "<>", "\"\<+n*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VS", "[", RowBox[{"[", 3,1, "]"}], "]"}], "]"}], "<>", "\"\<+\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VS", "[", RowBox[{"[", 3,2, "]"}], "]"}], "]"}]}], RowBox[{"\"\<]*Vto[\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VT", "[", RowBox[{"[", 2,1, "]"}], "]"}], "]"}], "<>", "\"\<+n2*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VT", "[", RowBox[{"[", 2,2, "]"}], "]"}], "]"}], "<>", "\"\<+n*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VT", "[", RowBox[{"[", 3,1, "]"}], "]"}], "]"}], "<>", "\"\<+\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VT", "[", RowBox[{"[", 3,2, "]"}], "]"}], "]"}], "<>", "\"\<]\>\""}]}]}], "]"}], "]"}]}], " ", ")"}]}]], "Input"
```
**Output [618]**
```
CellGroupData[{ Cell[[ RowBox[{"Map", "[", SpinLabelTerm,dDiagA, "]"}]], "Input"], Cell[[ RowBox[{ Join::heads, : , "\<\"Heads \\!\\(DiagonalDeltasNoPq\\) and \ \\!\\(List\\) at positions \\!\\(1\\) \ and \\!\\(2\\) are expected to be the same. \ \\!\\(\\*ButtonBox[\\\">>\\\", ButtonStyle->\\\"Link\\\", \ ButtonFrame->None, ButtonData:>\\\"paclet:ref/message/General/heads\\\", \ ButtonNote -> \\\"Join::heads\\\"]\\)\"\>"}]], "Message", "MSG"
```
**Input [619]**
```
[ RowBox[{"Map", "[", SpinLabelTerm,dDiagA, "]"}]], "Input"
```
**Output [632]**
```
[ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"2", " ", RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"DiagonalDeltasNoPq", "[", "]"}], ",", RowBox[{"ToUpperCase", "[", RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"y"\>},{\<"u"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"t"\>\<"v"\>},{\<"s"\>\<"x"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"x"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], "]"}], "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"y"\>},{\<"u"\>\<"w"\>}, "}"}]}], "}"}], "}"}]}], "]"}], "[[", 1,2, "]]"}]}]}], "}"}], "}"}], ",", RowBox[{"ToUpperCase", "[", RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"y"\>},{\<"u"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"t"\>\<"v"\>},{\<"s"\>\<"x"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"x"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], "]"}], "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"y"\>},{\<"u"\>\<"w"\>}, "}"}]}], "}"}], "}"}]}], "]"}], ",", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"2", " ", RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"DiagonalDeltasNoPq", "[", "]"}], ",", RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"y"\>},{\<"u"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"t"\>\<"v"\>},{\<"s"\>\<"x"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"x"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\<\"SIGN\"\>", ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"\<\"VTS\"\>", ",", RowBox[{"ToUpperCase", "[", "2", "]"}]}], "}"}], ",", {\<"R"\>\<"Y"\>}, ",", {\<"U"\>\<"W"\>}}], "}"}]}], "}"}], "}"}]}], "]"}], "[[", 1,2, "]]"}]}]}], "}"}], "}"}], ",", RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"y"\>},{\<"u"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"t"\>\<"v"\>},{\<"s"\>\<"x"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"x"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\<\"SIGN\"\>", ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"\<\"VTS\"\>", ",", RowBox[{"ToUpperCase", "[", "2", "]"}]}], "}"}], ",", {\<"R"\>\<"Y"\>}, ",", {\<"U"\>\<"W"\>}}], "}"}]}], "}"}], "}"}]}], "]"}], ",", RowBox[{"Join", "[", RowBox[{ RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"2", " ", RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"DiagonalDeltasNoPq", "[", "]"}], ",", RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"y"\>},{\<"u"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"t"\>\<"v"\>},{\<"s"\>\<"x"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"x"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"y"\>},{\<"u"\>\<"w"\>}, "}"}]}], "}"}], "}"}]}], "]"}], "[[", 1,2, "]]"}]}]}], "}"}], "}"}], ",", RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"y"\>},{\<"u"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<" ...[TRUNC 73665]
```
**Output [633]**
```
CellGroupData[{ Cell[[ RowBox[{"Map", "[", RowBox[{"IntegralsForTerm", ",", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SpinLabelTerm,dDiagA, "]"}], ",", "1"}], "]"}]}], "]"}]], "Input"
```
**Input [634]**
```
[ RowBox[{"Map", "[", RowBox[{"IntegralsForTerm", ",", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SpinLabelTerm,dDiagA, "]"}], ",", "1"}], "]"}]}], "]"}]], "Input"
```
**Output [647]**
```
[ InterpretationBox[ TagBox[ FrameBox[GridBox[{ { ItemBox[ TagBox[ RowBox[{"{", RowBox[{ RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{"UpperCaseQ", "[", RowBox[{ RowBox[{"ToUpperCase", "[", RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"y"\>},{\<"u"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"t"\>\<"v"\>},{\<"s"\>\<"x"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"x"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], "]"}], "]"}], "[[", 2,1, "]]"}], "]"}], "&&", RowBox[{"UpperCaseQ", "[", RowBox[{ RowBox[{"ToUpperCase", "[", RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"y"\>},{\<"u"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"t"\>\<"v"\>},{\<"s"\>\<"x"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"x"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}], "]"}], "]"}], "[[", 2,2, "]]"}], "]"}]}], ",", RowBox[{ TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], " ", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )]}], ",", RowBox[{ TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], " ", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )]}]}], "]"}], ",", TemplateBox[{"30"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{"\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], ",", RowBox[{"If", "[", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], "]"}]}], "}"}], Short[#, 5]& ], BaseStyle->{Deployed -> False}, StripOnInput->False]}, {GridBox[{ { TagBox[ TooltipBox[ StyleBox[ StyleBox[ DynamicBox[ToBoxes[ FEPrivate`FrontEndResource[ "FEStrings", "sizeBriefExplanation"], StandardForm], ImageSizeCache->{58., {2., 8.}}], StripOnInput->False, DynamicUpdating->True], "SuggestionsBarCategory", StripOnInput->False], StyleBox[ DynamicBox[ ToBoxes[ FEPrivate`FrontEndResource["FEStrings", "sizeExplanation"], StandardForm]], DynamicUpdating -> True, StripOnInput -> False]], Annotation[#, Style[ Dynamic[ FEPrivate`FrontEndResource["FEStrings", "sizeExplanation"]], DynamicUpdating -> True], "Tooltip"]& ], ButtonBox[ PaneSelectorBox[{False-> StyleBox[ StyleBox[ DynamicBox[ToBoxes[ FEPrivate`FrontEndResource["FEStrings", "sizeShowLess"], StandardForm], ImageSizeCache->{49., {0., 8.}}], StripOnInput->False, DynamicUpdating->True], "SuggestionsBarControl", StripOnInput->False], True-> StyleBox[ StyleBox[ DynamicBox[ToBoxes[ FEPrivate`FrontEndResource["FEStrings", "sizeShowLess"], StandardForm]], StripOnInput->False, DynamicUpdating->True], "SuggestionsBarControlActive", StripOnInput->False]}, Dynamic[ CurrentValue["MouseOver"]], Alignment->Center, FrameMargins->0, ImageSize->{Automatic, 25}], Appearance->None, ButtonFunction:>OutputSizeLimit`ButtonFunction[ Identity, 428, 16138051299609149988, 5/2], Enabled->True, Evaluator->Automatic, Method->"Queued"], ButtonBox[ PaneSelectorBox[{False-> StyleBox[ StyleBox[ DynamicBox[ToBoxes[ FEPrivate`FrontEndResource["FEStrings", "sizeShowMore"], StandardForm], ImageSizeCache->{54., {0., 8.}}], StripOnInput->False, DynamicUpdating->True], "SuggestionsBarControl", StripOnInput->False], Tru ...[TRUNC 8901]
```
**Output [648]**
```
CellGroupData[{ Cell[[ RowBox[{"Map", "[", RowBox[{"IntegralsForTerm", ",", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SpinLabelTerm,dDiagB, "]"}], ",", "1"}], "]"}]}], "]"}]], "Input"
```
**Input [649]**
```
[ RowBox[{"Map", "[", RowBox[{"IntegralsForTerm", ",", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SpinLabelTerm,dDiagB, "]"}], ",", "1"}], "]"}]}], "]"}]], "Input"
```
**Output [662]**
```
[ InterpretationBox[ TagBox[ FrameBox[GridBox[{ { ItemBox[ TagBox[ RowBox[{"{", RowBox[{ RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{"UpperCaseQ", "[", RowBox[{ RowBox[{"ToUpperCase", "[", RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"t"\>},{\<"w"\>\<"x"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"v"\>\<"y"\>},{\<"s"\>\<"u"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"x"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}]}], "}"}], "]"}], "]"}], "[[", 2,1, "]]"}], "]"}], "&&", RowBox[{"UpperCaseQ", "[", RowBox[{ RowBox[{"ToUpperCase", "[", RowBox[{"DiagonalDeltasNoPq", "[", RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"t"\>},{\<"w"\>\<"x"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"v"\>\<"y"\>},{\<"s"\>\<"u"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"x"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}]}], "}"}], "]"}], "]"}], "[[", 2,2, "]]"}], "]"}]}], ",", RowBox[{ TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], " ", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )]}], ",", RowBox[{ TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], " ", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )]}]}], "]"}], ",", TemplateBox[{"30"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{"\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], ",", RowBox[{"If", "[", TemplateBox[{"1"}, "OutputSizeLimit`Skeleton", DisplayFunction->(FrameBox[ RowBox[{ "\" ...\"", #, "\"... \""}], Background -> GrayLevel[0.75], BaseStyle -> { "Deploy", FontColor -> GrayLevel[1], FontSize -> Smaller, ShowStringCharacters -> False}, BaselinePosition -> Baseline, ContentPadding -> False, FrameMargins -> 1, FrameStyle -> GrayLevel[0.75], RoundingRadius -> 7]& )], "]"}]}], "}"}], Short[#, 5]& ], BaseStyle->{Deployed -> False}, StripOnInput->False]}, {GridBox[{ { TagBox[ TooltipBox[ StyleBox[ StyleBox[ DynamicBox[ToBoxes[ FEPrivate`FrontEndResource[ "FEStrings", "sizeBriefExplanation"], StandardForm], ImageSizeCache->{58., {2., 8.}}], StripOnInput->False, DynamicUpdating->True], "SuggestionsBarCategory", StripOnInput->False], StyleBox[ DynamicBox[ ToBoxes[ FEPrivate`FrontEndResource["FEStrings", "sizeExplanation"], StandardForm]], DynamicUpdating -> True, StripOnInput -> False]], Annotation[#, Style[ Dynamic[ FEPrivate`FrontEndResource["FEStrings", "sizeExplanation"]], DynamicUpdating -> True], "Tooltip"]& ], ButtonBox[ PaneSelectorBox[{False-> StyleBox[ StyleBox[ DynamicBox[ToBoxes[ FEPrivate`FrontEndResource["FEStrings", "sizeShowLess"], StandardForm], ImageSizeCache->{49., {0., 8.}}], StripOnInput->False, DynamicUpdating->True], "SuggestionsBarControl", StripOnInput->False], True-> StyleBox[ StyleBox[ DynamicBox[ToBoxes[ FEPrivate`FrontEndResource["FEStrings", "sizeShowLess"], StandardForm]], StripOnInput->False, DynamicUpdating->True], "SuggestionsBarControlActive", StripOnInput->False]}, Dynamic[ CurrentValue["MouseOver"]], Alignment->Center, FrameMargins->0, ImageSize->{Automatic, 25}], Appearance->None, ButtonFunction:>OutputSizeLimit`ButtonFunction[ Identity, 429, 16138051299609149988, 5/2], Enabled->True, Evaluator->Automatic, Method->"Queued"], ButtonBox[ PaneSelectorBox[{False-> StyleBox[ StyleBox[ DynamicBox[ToBoxes[ FEPrivate`FrontEndResource["FEStrings", "sizeShowMore"], StandardForm], ImageSizeCache->{54., {0., 8.}}], StripOnInput->False, DynamicUpdating->True], "SuggestionsBarControl", StripOnInput->False], T ...[TRUNC 8903]
```
**Input [664]**
```
[{ RowBox[{ RowBox[{ RowBox[{"DiagonalDeltasNoPq", "[", "aterm_", "]"}], ":=", RowBox[{ RowBox[{"Reap", "[", RowBox[{"Map", "[", RowBox[{ RowBox[{ RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{"(", RowBox[{ RowBox[{ RowBox[{"#", "[", RowBox[{"[", "1", "]"}], "]"}], "===", RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}]}], " ", "||", " ", RowBox[{ RowBox[{"#", "[", RowBox[{"[", "1", "]"}], "]"}], "===", RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}]}]}], ")"}], "&&", RowBox[{"(", RowBox[{ RowBox[{"!", RowBox[{"MemberQ", "[", RowBox[{ RowBox[{"Flatten", "[", "#", "]"}], ",", "\"\<p\>\""}], "]"}]}], "&&", RowBox[{"!", RowBox[{"MemberQ", "[", RowBox[{ RowBox[{"Flatten", "[", "#", "]"}], ",", "\"\<q\>\""}], "]"}]}]}], ")"}]}], ",", RowBox[{"Sow", "[", RowBox[{"ADeltaForFactor", "[", "#", "]"}], "]"}], ","}], "]"}], "&"}], ",", "aterm"}], "]"}], "]"}], "[", RowBox[{"[", 2,1, "]"}], "]"}]}], ";", RowBox[{ RowBox[{"DiagonalGammaNoP", "[", "aterm_", "]"}], ":=", RowBox[{"Join", "[", RowBox[{"aterm", ",", RowBox[{"DiagonalDeltasNoP", "[", "aterm", "]"}]}], "]"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"dDiagA", "=", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", DiagonalGammaNoP,DiagA, "]"}]}], "]"}]}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"dDiagB", "=", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", DiagonalGammaNoP,DiagB, "]"}]}], "]"}]}], "]"}]}], "Input"
```
**Input [667]**
```
[ giving these terms the usual labeling we obtain... (here's MP2 like terms) :], "Input"
```
**Output [668]**
```
CellGroupData[{ Cell[[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{jy}}_{{bc}}, " ", "VTt"^{{bc}}_{{ij}}, " ", g^{{"a"}}_{{"k"}}}], " ", "+", RowBox[{ (1)/(2), " ", "VTs"^{{bc}}_{{ij}}, " ", "VTt"^{{jy}}_{{bc}}, " ", g^{{"a"}}_{{"k"}}}], " ", "-", RowBox[{ (1)/(2), " ", "VTs"^{{sy}}_{{uw}}, " ", "VTt"^{{uw}}_{{qs}}, " ", g^{{"u"}}_{{"u"}}, " ", g^{{"w"}}_{{"w"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"y"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{uw}}_{{qs}}, " ", "VTt"^{{sy}}_{{uw}}, " ", g^{{"u"}}_{{"u"}}, " ", g^{{"w"}}_{{"w"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"y"}}}]}]], "Input"
```
**Input [669]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{jy}}_{{bc}}, " ", "VTt"^{{bc}}_{{ij}}, " ", g^{{"a"}}_{{"k"}}}], " ", "+", RowBox[{ (1)/(2), " ", "VTs"^{{bc}}_{{ij}}, " ", "VTt"^{{jy}}_{{bc}}, " ", g^{{"a"}}_{{"k"}}}], " ", "-", RowBox[{ (1)/(2), " ", "VTs"^{{sy}}_{{uw}}, " ", "VTt"^{{uw}}_{{qs}}, " ", g^{{"u"}}_{{"u"}}, " ", g^{{"w"}}_{{"w"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"y"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{uw}}_{{qs}}, " ", "VTt"^{{sy}}_{{uw}}, " ", g^{{"u"}}_{{"u"}}, " ", g^{{"w"}}_{{"w"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"y"}}}]}]], "Input"
```
**Output [670]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{bc}}_{{ij}}, " ", "VTt"^{{jy}}_{{bc}}, " ", g^{{"a"}}_{{"k"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{jy}}_{{bc}}, " ", "VTt"^{{bc}}_{{ij}}, " ", g^{{"a"}}_{{"k"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{sy}}_{{uw}}, " ", "VTt"^{{uw}}_{{qs}}, " ", g^{{"u"}}_{{"u"}}, " ", g^{{"w"}}_{{"w"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"y"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{uw}}_{{qs}}, " ", "VTt"^{{sy}}_{{uw}}, " ", g^{{"u"}}_{{"u"}}, " ", g^{{"w"}}_{{"w"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"y"}}}]}]], "Output"
```
> Verify that every Term spin adapts the same way.
**Input [672]**
```
CellGroupData[{ Cell[[ RowBox[{"TlS", "[", "\"\<P\>\"", "]"}]], "Input"
```
**Input [673]**
```
[ RowBox[{"TlS", "[", "\"\<P\>\"", "]"}]], "Input"
```
**Output [674]**
```
["\<\"p\"\>"], "Output"
```
**Input [675]**
```
[ RowBox[{ RowBox[{"IntegralsForTerm", "[", "ATerm_", "]"}], ":=", RowBox[{"(", " ", RowBox[{ RowBox[{"VS", "=", RowBox[{"ATerm", "[", RowBox[{"[", "2", "]"}], "]"}]}], ";", " ", RowBox[{"VT", "=", RowBox[{"ATerm", "[", RowBox[{"[", "3", "]"}], "]"}]}], ";", " ", RowBox[{"samespin", "=", RowBox[{ RowBox[{"UpperCaseQ", "[", RowBox[{"VS", "[", RowBox[{"[", 2,1, "]"}], "]"}], "]"}], "&&", RowBox[{"UpperCaseQ", "[", RowBox[{"VS", "[", RowBox[{"[", 2,2, "]"}], "]"}], "]"}]}]}], ";", " ", RowBox[{"Return", "[", RowBox[{"If", "[", RowBox[{"samespin", ",", RowBox[{ RowBox[{"\"\<Vss[n3*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VS", "[", RowBox[{"[", 2,1, "]"}], "]"}], "]"}], "<>", "\"\<+n2*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VS", "[", RowBox[{"[", 2,2, "]"}], "]"}], "]"}], "<>", "\"\<+n*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VS", "[", RowBox[{"[", 3,1, "]"}], "]"}], "]"}], "<>", "\"\<+\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VS", "[", RowBox[{"[", 3,2, "]"}], "]"}], "]"}]}], RowBox[{"\"\<]*Vts[\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VT", "[", RowBox[{"[", 2,1, "]"}], "]"}], "]"}], "<>", "\"\<+n2*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VT", "[", RowBox[{"[", 2,2, "]"}], "]"}], "]"}], "<>", "\"\<+n*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VT", "[", RowBox[{"[", 3,1, "]"}], "]"}], "]"}], "<>", "\"\<+\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VT", "[", RowBox[{"[", 3,2, "]"}], "]"}], "]"}], "<>", "\"\<]\>\""}]}], ",", RowBox[{ RowBox[{"\"\<Vso[n3*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VS", "[", RowBox[{"[", 2,1, "]"}], "]"}], "]"}], "<>", "\"\<+n2*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VS", "[", RowBox[{"[", 2,2, "]"}], "]"}], "]"}], "<>", "\"\<+n*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VS", "[", RowBox[{"[", 3,1, "]"}], "]"}], "]"}], "<>", "\"\<+\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VS", "[", RowBox[{"[", 3,2, "]"}], "]"}], "]"}]}], RowBox[{"\"\<]*Vto[\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VT", "[", RowBox[{"[", 2,1, "]"}], "]"}], "]"}], "<>", "\"\<+n2*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VT", "[", RowBox[{"[", 2,2, "]"}], "]"}], "]"}], "<>", "\"\<+n*\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VT", "[", RowBox[{"[", 3,1, "]"}], "]"}], "]"}], "<>", "\"\<+\>\"", "<>", RowBox[{"TlS", "[", RowBox[{"VT", "[", RowBox[{"[", 3,2, "]"}], "]"}], "]"}], "<>", "\"\<]\>\""}]}]}], "]"}], "]"}]}], " ", ")"}]}]], "Input"
```
**Output [676]**
```
CellGroupData[{ Cell[[ RowBox[{"Map", "[", SpinLabelTerm,dDiagA, "]"}]], "Input"], Cell[[ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"Y"\>\<"s"\>},{\<"U"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"U"\>\<"w"\>},{\<"Q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"U"\>}, ",", {\<"U"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"w"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"Q"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"Y"\>\<"s"\>},{\<"W"\>\<"u"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"W"\>\<"u"\>},{\<"Q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"W"\>}, ",", {\<"W"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"Q"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"S"\>\<"Y"\>},{\<"U"\>\<"W"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"U"\>\<"W"\>},{\<"Q"\>\<"S"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"U"\>}, ",", {\<"U"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"W"\>}, ",", {\<"W"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"Q"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"S"\>}, ",", {\<"S"\>}}], "}"}]}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"S"\>\<"u"\>},{\<"Q"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"Y"\>\<"w"\>},{\<"S"\>\<"u"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"w"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"Q"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"S"\>}, ",", {\<"S"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"u"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"U"\>\<"s"\>},{\<"Q"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"Y"\>\<"w"\>},{\<"U"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"w"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"Q"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"U"\>}, ",", {\<"U"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"S"\>\<"U"\>},{\<"Q"\>\<"W"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"W"\>\<"Y"\>},{\<"S"\>\<"U"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"W"\>}, ",", {\<"W"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"Q"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"S"\>}, ",", {\<"S"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"U"\>}, ",", {\<"U"\>}}], "}"}]}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"U"\>\<"w"\>},{\<"Q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"X"\>\<"s"\>},{\<"U"\>\<"w"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"Q"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"U"\>}, ",", {\<"U"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"w"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"X"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"W"\>\<"u"\>},{\<"Q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"X"\>\<"s"\>},{\<"W"\>\<"u"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"Q"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox ...[TRUNC 9459]
```
**Input [677]**
```
[ RowBox[{"Map", "[", SpinLabelTerm,dDiagA, "]"}]], "Input"
```
**Output [678]**
```
[ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"Y"\>\<"s"\>},{\<"U"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"U"\>\<"w"\>},{\<"Q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"U"\>}, ",", {\<"U"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"w"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"Q"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"Y"\>\<"s"\>},{\<"W"\>\<"u"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"W"\>\<"u"\>},{\<"Q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"W"\>}, ",", {\<"W"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"Q"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"S"\>\<"Y"\>},{\<"U"\>\<"W"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"U"\>\<"W"\>},{\<"Q"\>\<"S"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"U"\>}, ",", {\<"U"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"W"\>}, ",", {\<"W"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"Q"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"S"\>}, ",", {\<"S"\>}}], "}"}]}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"S"\>\<"u"\>},{\<"Q"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"Y"\>\<"w"\>},{\<"S"\>\<"u"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"w"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"Q"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"S"\>}, ",", {\<"S"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"u"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"U"\>\<"s"\>},{\<"Q"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"Y"\>\<"w"\>},{\<"U"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"w"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"Q"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"U"\>}, ",", {\<"U"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"S"\>\<"U"\>},{\<"Q"\>\<"W"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"W"\>\<"Y"\>},{\<"S"\>\<"U"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"Y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"W"\>}, ",", {\<"W"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"Q"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"S"\>}, ",", {\<"S"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"U"\>}, ",", {\<"U"\>}}], "}"}]}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"U"\>\<"w"\>},{\<"Q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"X"\>\<"s"\>},{\<"U"\>\<"w"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"Q"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"U"\>}, ",", {\<"U"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"w"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"P"\>}, ",", {\<"X"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"W"\>\<"u"\>},{\<"Q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"X"\>\<"s"\>},{\<"W"\>\<"u"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"Q"\>}, ",", {\<"Q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"W"\>}, ",", {\<"W"\>}}], "}"}], ...[TRUNC 9371]
```
**Output [679]**
```
CellGroupData[{ Cell[[ RowBox[{"Map", "[", RowBox[{"IntegralsForTerm", ",", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SpinLabelTerm,dDiagA, "]"}], ",", "1"}], "]"}]}], "]"}]], "Input"
```
**Input [680]**
```
[ RowBox[{"Map", "[", RowBox[{"IntegralsForTerm", ",", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SpinLabelTerm,dDiagA, "]"}], ",", "1"}], "]"}]}], "]"}]], "Input"
```
**Output [681]**
```
[ RowBox[{"{", RowBox[{ RowBox[{"\<\"Vso[n3*y+n2*s+n*u+w\"\>", " ", "\<\"]*Vto[u+n2*w+n*q+s]\"\>"}], ",", RowBox[{"\<\"Vso[n3*y+n2*s+n*w+u\"\>", " ", "\<\"]*Vto[w+n2*u+n*q+s]\"\>"}], ",", RowBox[{"\<\"Vss[n3*s+n2*y+n*u+w\"\>", " ", "\<\"]*Vts[u+n2*w+n*q+s]\"\>"}], ",", RowBox[{"\<\"Vso[n3*s+n2*u+n*q+w\"\>", " ", "\<\"]*Vto[y+n2*w+n*s+u]\"\>"}], ",", RowBox[{"\<\"Vso[n3*u+n2*s+n*q+w\"\>", " ", "\<\"]*Vto[y+n2*w+n*u+s]\"\>"}], ",", RowBox[{"\<\"Vss[n3*s+n2*u+n*q+w\"\>", " ", "\<\"]*Vts[w+n2*y+n*s+u]\"\>"}], ",", RowBox[{"\<\"Vso[n3*u+n2*w+n*q+s\"\>", " ", "\<\"]*Vto[x+n2*s+n*u+w]\"\>"}], ",", RowBox[{"\<\"Vso[n3*w+n2*u+n*q+s\"\>", " ", "\<\"]*Vto[x+n2*s+n*w+u]\"\>"}], ",", RowBox[{"\<\"Vss[n3*u+n2*w+n*q+s\"\>", " ", "\<\"]*Vts[s+n2*x+n*u+w]\"\>"}], ",", RowBox[{"\<\"Vso[n3*x+n2*w+n*s+u\"\>", " ", "\<\"]*Vto[s+n2*u+n*q+w]\"\>"}], ",", RowBox[{"\<\"Vso[n3*x+n2*w+n*u+s\"\>", " ", "\<\"]*Vto[u+n2*s+n*q+w]\"\>"}], ",", RowBox[{"\<\"Vss[n3*w+n2*x+n*s+u\"\>", " ", "\<\"]*Vts[s+n2*u+n*q+w]\"\>"}]}], "}"}]], "Output"
```
**Output [682]**
```
CellGroupData[{ Cell[[ RowBox[{"Map", "[", RowBox[{"IntegralsForTerm", ",", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SpinLabelTerm,dDiagB, "]"}], ",", "1"}], "]"}]}], "]"}]], "Input"
```
**Input [683]**
```
[ RowBox[{"Map", "[", RowBox[{"IntegralsForTerm", ",", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SpinLabelTerm,dDiagB, "]"}], ",", "1"}], "]"}]}], "]"}]], "Input"
```
**Output [684]**
```
[ RowBox[{"{", RowBox[{ RowBox[{"\<\"Vso[n3*s+n2*u+n*q+w\"\>", " ", "\<\"]*Vto[y+n2*w+n*s+u]\"\>"}], ",", RowBox[{"\<\"Vso[n3*u+n2*s+n*q+w\"\>", " ", "\<\"]*Vto[y+n2*w+n*u+s]\"\>"}], ",", RowBox[{"\<\"Vss[n3*s+n2*u+n*q+w\"\>", " ", "\<\"]*Vts[w+n2*y+n*s+u]\"\>"}], ",", RowBox[{"\<\"Vso[n3*y+n2*w+n*s+u\"\>", " ", "\<\"]*Vto[s+n2*u+n*q+w]\"\>"}], ",", RowBox[{"\<\"Vso[n3*y+n2*w+n*u+s\"\>", " ", "\<\"]*Vto[u+n2*s+n*q+w]\"\>"}], ",", RowBox[{"\<\"Vss[n3*w+n2*y+n*s+u\"\>", " ", "\<\"]*Vts[s+n2*u+n*q+w]\"\>"}], ",", RowBox[{"\<\"Vso[n3*y+n2*s+n*u+w\"\>", " ", "\<\"]*Vto[u+n2*w+n*q+s]\"\>"}], ",", RowBox[{"\<\"Vso[n3*y+n2*s+n*w+u\"\>", " ", "\<\"]*Vto[w+n2*u+n*q+s]\"\>"}], ",", RowBox[{"\<\"Vss[n3*s+n2*y+n*u+w\"\>", " ", "\<\"]*Vts[u+n2*w+n*q+s]\"\>"}], ",", RowBox[{"\<\"Vso[n3*u+n2*w+n*q+s\"\>", " ", "\<\"]*Vto[y+n2*s+n*u+w]\"\>"}], ",", RowBox[{"\<\"Vso[n3*w+n2*u+n*q+s\"\>", " ", "\<\"]*Vto[y+n2*s+n*w+u]\"\>"}], ",", RowBox[{"\<\"Vss[n3*u+n2*w+n*q+s\"\>", " ", "\<\"]*Vts[s+n2*y+n*u+w]\"\>"}]}], "}"}]], "Output"
```
**Output [685]**
```
CellGroupData[{ Cell[["dDiagA"], "Input"], Cell[[ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"s"\>\<"y"\>},{\<"u"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"u"\>\<"w"\>},{\<"q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"w"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"q"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"s"\>\<"u"\>},{\<"q"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"w"\>\<"y"\>},{\<"s"\>\<"u"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"w"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"q"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"u"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"u"\>\<"w"\>},{\<"q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"s"\>\<"x"\>},{\<"u"\>\<"w"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"q"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"w"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"x"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"w"\>\<"x"\>},{\<"s"\>\<"u"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"s"\>\<"u"\>},{\<"q"\>\<"w"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"q"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"w"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"x"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"u"\>}}], "}"}]}], "}"}]}], "}"}]], "Output"
```
**Input [686]**
```
["dDiagA"], "Input"
```
**Output [687]**
```
[ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"s"\>\<"y"\>},{\<"u"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"u"\>\<"w"\>},{\<"q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"w"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"q"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"s"\>\<"u"\>},{\<"q"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"w"\>\<"y"\>},{\<"s"\>\<"u"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"w"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"q"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"u"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"u"\>\<"w"\>},{\<"q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"s"\>\<"x"\>},{\<"u"\>\<"w"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"q"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"u"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"w"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"x"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"w"\>\<"x"\>},{\<"s"\>\<"u"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"s"\>\<"u"\>},{\<"q"\>\<"w"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"q"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"w"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"x"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"s"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"u"\>}}], "}"}]}], "}"}]}], "}"}]], "Output"
```
**Input [688]**
```
[ Every term has two opposite spin copies V_os*V_os and one same spin. There are no terms which are ss*os.], "Input"
```
**Output [689]**
```
CellGroupData[{ Cell[[{ RowBox[{ -(1)/(2), " ", "VTs"^{{vv}}_{{ox}}, " ", "VTt"^{{oy}}_{{vv}}, " ", g^{{"a"}}_{{"q"}}, " ", g^{{"w"}}_{{"w"}}, " ", g^{{"p"}}_{{"b"}}, " ", eta^{{"s"}}_{{"s"}}, " ", eta^{{"u"}}_{{"u"}}}], " ", RowBox[{ +(1)/(2), " ", "VTs"^{{si}}_{{uw}}, " ", "VTt"^{{uw}}_{{vb}}, " ", g^{{"u"}}_{{"u"}}, " ", g^{{"w"}}_{{"w"}}, " ", eta^{{"b"}}_{{"j"}}, " ", eta^{{"s"}}_{{"s"}}, " ", eta^{{"a"}}_{{"i"}}}]}], "Input"
```
**Output [691]**
```
[ RowBox[{ -(1)/(2), " ", "VTs"^{{vv}}_{{ox}}, " ", "VTt"^{{oy}}_{{vv}}, " ", g^{{"p"}}_{{"b"}}, " ", g^{{"a"}}_{{"q"}}, " ", g^{{"w"}}_{{"w"}}, " ", eta^{{"s"}}_{{"s"}}, " ", eta^{{"u"}}_{{"u"}}}]], "Output"
```
**Output [692]**
```
[ RowBox[{ (1)/(2), " ", "VTs"^{{si}}_{{uw}}, " ", "VTt"^{{uw}}_{{vb}}, " ", g^{{"u"}}_{{"u"}}, " ", g^{{"w"}}_{{"w"}}, " ", eta^{{"a"}}_{{"i"}}, " ", eta^{{"b"}}_{{"j"}}, " ", eta^{{"s"}}_{{"s"}}}]], "Output"
```
**Output [693]**
```
CellGroupData[{ Cell[[{ RowBox[{"Map", "[", NumIndicesInTerm,FinalTCL2, "]"}], " ", RowBox[{"Map", "[", NumIndicesInTerm,dtcl2, "]"}]}], "Input"
```
**Input [694]**
```
[{ RowBox[{"Map", "[", NumIndicesInTerm,FinalTCL2, "]"}], " ", RowBox[{"Map", "[", NumIndicesInTerm,dtcl2, "]"}]}], "Input"
```
**Output [695]**
```
[ {1010101010101010}], "Output"
```
**Output [696]**
```
["dtcl2"], "Output"
```


### [698] Secular 8 version.
**Output [700]**
```
CellGroupData[{ Cell[[{ RowBox[{"RecognizableForm", "[", "FromCode", "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"FinalTCL2", "=", RowBox[{"Map", "[", CollapseDeltas,FromCode, "]"}]}], "]"}]}], "Input"
```
**Input [701]**
```
[{ RowBox[{"RecognizableForm", "[", "FromCode", "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"FinalTCL2", "=", RowBox[{"Map", "[", CollapseDeltas,FromCode, "]"}]}], "]"}]}], "Input"
```
**Output [702]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{ry}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"y"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{wx}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", d^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{vy}}_{{su}}, " ", "VTt"^{{rt}}_{{wx}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{tv}}_{{sy}}, " ", "VTt"^{{rx}}_{{uw}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"y"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"x"}}}]}]], "Output"
```
**Output [703]**
```
[ RowBox[{ RowBox[{ -(1)/(2), " ", "VTs"^{{pr}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pv}}_{{su}}, " ", "VTt"^{{rt}}_{{wx}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{qw}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{tv}}_{{qs}}, " ", "VTt"^{{rx}}_{{uw}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"x"}}}]}]], "Output"
```
**Output [704]**
```
CellGroupData[{ Cell[[{ RowBox[{"FinalTCL2", "[", RowBox[{"[", "1", "]"}], "]"}], " ", RowBox[{"FinalTCL2", "[", RowBox[{"[", "2", "]"}], "]"}], " ", RowBox[{"FinalTCL2", "[", RowBox[{"[", "3", "]"}], "]"}], " ", RowBox[{"FinalTCL2", "[", RowBox[{"[", "4", "]"}], "]"}]}], "Input"
```
**Input [705]**
```
[{ RowBox[{"FinalTCL2", "[", RowBox[{"[", "1", "]"}], "]"}], " ", RowBox[{"FinalTCL2", "[", RowBox[{"[", "2", "]"}], "]"}], " ", RowBox[{"FinalTCL2", "[", RowBox[{"[", "3", "]"}], "]"}], " ", RowBox[{"FinalTCL2", "[", RowBox[{"[", "4", "]"}], "]"}]}], "Input"
```
**Output [706]**
```
[ RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"r"\>},{\<"u"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"t"\>\<"v"\>},{\<"s"\>\<"x"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"x"\>}, ",", {\<"q"\>}}], "}"}]}], "}"}]], "Output"
```
**Output [707]**
```
[ RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"r"\>\<"t"\>},{\<"q"\>\<"w"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"v"\>\<"y"\>},{\<"s"\>\<"u"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"y"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}]}], "}"}]], "Output"
```
**Output [708]**
```
[ RowBox[{"{", RowBox[{ {\<"Sign"\>(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"v"\>},{\<"s"\>\<"u"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"r"\>\<"t"\>},{\<"w"\>\<"x"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"x"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}]}], "}"}]], "Output"
```
**Output [709]**
```
[ RowBox[{"{", RowBox[{ {\<"Sign"\>-(1)/(2)}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"t"\>\<"v"\>},{\<"q"\>\<"s"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"r"\>\<"x"\>},{\<"u"\>\<"w"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"u"\>}, ",", {\<"t"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"w"\>}, ",", {\<"v"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"x"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"s"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}]], "Output"
```
> <Note that the middle two fail to simplify out (first two terms should cancel), because there is no convention on relabeling repeated indices.  I can fix that like this : >
**Output [714]**
```
CellGroupData[{ Cell[[{ RowBox[{"RecognizableForm", "[", RowBox[{"AssociativeSimplify", "[", RowBox[{"Map", "[", RowBox[{"RelabelTerm", ",", RowBox[{"Map", "[", SecIndexRelabel,sec8d, "]"}]}], "]"}], "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"dCC", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", DiagonalGamma1,FromCode, "]"}]}], "]"}]}], "]"}]}], "]"}]}], "Input"
```
**Input [715]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{"AssociativeSimplify", "[", RowBox[{"Map", "[", RowBox[{"RelabelTerm", ",", RowBox[{"Map", "[", SecIndexRelabel,sec8d, "]"}]}], "]"}], "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"dCC", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", DiagonalGamma1,FromCode, "]"}]}], "]"}]}], "]"}]}], "]"}]}], "Input"
```
**Output [716]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{pr}}, " ", "VTt"^{{pr}}_{{st}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pt}}_{{rs}}, " ", "VTt"^{{rs}}_{{pt}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rs}}_{{pt}}, " ", "VTt"^{{pt}}_{{rs}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}]}]], "Output"
```
**Output [717]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{qr}}, " ", "VTt"^{{pr}}_{{st}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{qr}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rs}}_{{qt}}, " ", "VTt"^{{pt}}_{{rs}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pt}}_{{rs}}, " ", "VTt"^{{rs}}_{{qt}}, " ", g^{{"q"}}_{{"q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}]}]], "Output"
```
> Finally attacking the simplification of the main problem.
**Output [719]**
```
CellGroupData[{ Cell[[ RowBox[{"RecognizableForm", "[", RowBox[{"sec8", "=", RowBox[{"AssociativeSimplify", "[", RowBox[{"Map", "[", RowBox[{"RelabelTerm", ",", RowBox[{"Map", "[", RowBox[{"SecIndexRelabel", ",", RowBox[{ RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", RowBox[{"RelabelTerm", ",", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SecularizeATerm,FinalTCL2, "]"}], ",", "1"}], "]"}]}], "]"}]}], "]"}]}], "]"}], "//.", "NoZeros"}]}], "]"}]}], "]"}], "]"}]}], "]"}]], "Input"
```
**Input [720]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"sec8", "=", RowBox[{"AssociativeSimplify", "[", RowBox[{"Map", "[", RowBox[{"RelabelTerm", ",", RowBox[{"Map", "[", RowBox[{"SecIndexRelabel", ",", RowBox[{ RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", RowBox[{"RelabelTerm", ",", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SecularizeATerm,FinalTCL2, "]"}], ",", "1"}], "]"}]}], "]"}]}], "]"}]}], "]"}], "//.", "NoZeros"}]}], "]"}]}], "]"}], "]"}]}], "]"}]], "Input"
```
**Output [721]**
```
[ RowBox[{ RowBox[{ RowBox[{"-", "VTs"^{{pr}}_{{pr}}}], " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"p"}}_{{"r"}}, " ", eta^{{"r"}}_{{"p"}}, " ", eta^{{"p"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"r"}}_{{"p"}}, " ", eta^{{"p"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{ps}}_{{pr}}, " ", "VTt"^{{pr}}_{{ps}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"r"}}_{{"p"}}, " ", eta^{{"p"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{ps}}_{{ps}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"s"}}_{{"p"}}, " ", eta^{{"p"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"r"}}_{{"p"}}, " ", g^{{"p"}}_{{"r"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"q"}}_{{"r"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"q"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{ps}}, " ", "VTt"^{{ps}}_{{pr}}, " ", g^{{"s"}}_{{"p"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qs}}_{{qr}}, " ", "VTt"^{{qr}}_{{qs}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"q"}}_{{"s"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{rs}}, " ", "VTt"^{{rs}}_{{pr}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{rs}}_{{qr}}, " ", "VTt"^{{qr}}_{{rs}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{ps}}_{{ps}}, " ", g^{{"r"}}_{{"p"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"s"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qs}}_{{qs}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"q"}}_{{"s"}}, " ", eta^{{"s"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"p"}}_{{"r"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"q"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"q"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{qr}}_{{qs}}, " ", "VTt"^{{qs}}_{{qr}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"q"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qs}}_{{qs}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"s"}}_{{"q"}}, " ", eta^{{"q"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"p"}}_{{"r"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{ps}}_{{pr}}, " ", "VTt"^{{pr}}_{{ps}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{ps}}_{{pr}}, " ", "VTt"^{{pr}}_{{ps}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"r"}}_{{"p"}}, " ", g^{{"p"}}_{{"r"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"r"}}_{{"q"}}, " ", g^{{"q"}}_{{"r"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"q"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{ps}}, " ", "VTt"^{{ps}}_{{pr}}, " ", g^{{"s"}}_{{"p"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{qs}}_{{qr}}, " ", "VTt"^{{qr}}_{{qs}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"q"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{rs}}, " ", "VTt"^{{rs}}_{{pr}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{rs}}_{{qr}}, " ", "VTt"^{{qr}}_{{rs}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{ps}}, " ", "VTt"^{{ps}}_{{pr}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qs}}_{{qr}}, " ", "VTt"^{{qr}}_{{qs}}, " ", g^{{"q"}}_{{"q"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{rs}}, " ", "VTt"^{{rs}}_{{pr}}, " ", g^{{"r"}}_{{"r"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{rs}}_{{qr}}, " ", "VTt"^{{qr}}_{{rs}}, " ", g^{{"r"}}_{{"r"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"s"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{qr}}, " ", "VTt"^{{qr}}_{{st}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"s"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{ ...[TRUNC 22424]
```
**Input [722]**
```
CellGroupData[{ Cell[[ RowBox[{"Length", "[", "sec8", "]"}]], "Input"
```
**Input [723]**
```
[ RowBox[{"Length", "[", "sec8", "]"}]], "Input"
```
**Output [724]**
```
["136"], "Output"
```
**Output [725]**
```
CellGroupData[{ Cell[[ RowBox[{"sec8", "[", RowBox[{"[", "1", "]"}], "]"}]], "Input"
```
**Input [726]**
```
[ RowBox[{"sec8", "[", RowBox[{"[", "1", "]"}], "]"}]], "Input"
```
**Output [727]**
```
[ RowBox[{"{", RowBox[{ {\<"Sign"\>-1}, ",", RowBox[{"{", {\<"VTs"\>2},{\<"p"\>\<"r"\>},{\<"p"\>\<"r"\>}, "}"}], ",", RowBox[{"{", {\<"VTt"\>2},{\<"p"\>\<"r"\>},{\<"p"\>\<"r"\>}, "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"q"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"r"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"r"\>}, ",", {\<"p"\>}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"eta", ",", "1"}], "}"}], ",", {\<"p"\>}, ",", {\<"r"\>}}], "}"}]}], "}"}]], "Output"
```
**Input [728]**
```
[ RowBox[{ RowBox[{"NoZerosDenoms", "=", RowBox[{"(", RowBox[{ RowBox[{"{", RowBox[{"any1___", ",", RowBox[{"{", RowBox[{ {\<Sign\>g_}, ",", RowBox[{"{", {\<VTs\>2},{h_i_},{h_i_}, "}"}], ",", " ", "any3___"}], "}"}], ",", "any2___"}], "}"}], "->", {any1any2}}], ")"}]}], ";"}]], "Input"
```
**Output [729]**
```
CellGroupData[{ Cell[[ RowBox[{"RecognizableForm", "[", sec8//.NoZerosDenoms, "]"}]], "Input"
```
**Input [730]**
```
[ RowBox[{"RecognizableForm", "[", sec8//.NoZerosDenoms, "]"}]], "Input"
```
**Output [731]**
```
[ RowBox[{ RowBox[{ "VTs"^{{ps}}_{{pr}}, " ", "VTt"^{{pr}}_{{ps}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"r"}}_{{"p"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{ps}}, " ", "VTt"^{{ps}}_{{pr}}, " ", g^{{"s"}}_{{"p"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qs}}_{{qr}}, " ", "VTt"^{{qr}}_{{qs}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"q"}}_{{"s"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{rs}}, " ", "VTt"^{{rs}}_{{pr}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{rs}}_{{qr}}, " ", "VTt"^{{qr}}_{{rs}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"p"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{qr}}_{{qs}}, " ", "VTt"^{{qs}}_{{qr}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"q"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{ps}}_{{pr}}, " ", "VTt"^{{pr}}_{{ps}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{ps}}_{{pr}}, " ", "VTt"^{{pr}}_{{ps}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{ps}}, " ", "VTt"^{{ps}}_{{pr}}, " ", g^{{"s"}}_{{"p"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{qs}}_{{qr}}, " ", "VTt"^{{qr}}_{{qs}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"q"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{rs}}, " ", "VTt"^{{rs}}_{{pr}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{rs}}_{{qr}}, " ", "VTt"^{{qr}}_{{rs}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{ps}}, " ", "VTt"^{{ps}}_{{pr}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qs}}_{{qr}}, " ", "VTt"^{{qr}}_{{qs}}, " ", g^{{"q"}}_{{"q"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{rs}}, " ", "VTt"^{{rs}}_{{pr}}, " ", g^{{"r"}}_{{"r"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{rs}}_{{qr}}, " ", "VTt"^{{qr}}_{{rs}}, " ", g^{{"r"}}_{{"r"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"s"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{qr}}, " ", "VTt"^{{qr}}_{{st}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"s"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{qr}}, " ", "VTt"^{{qr}}_{{st}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qr}}_{{qs}}, " ", "VTt"^{{qs}}_{{qr}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{qr}}_{{qs}}, " ", "VTt"^{{qs}}_{{qr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{pt}}, " ", "VTt"^{{st}}_{{rs}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"p"}}_{{"t"}}, " ", eta^{{"s"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qr}}_{{qt}}, " ", "VTt"^{{st}}_{{rs}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"p"}}_{{"t"}}, " ", eta^{{"s"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{ps}}, " ", "VTt"^{{st}}_{{rt}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"p"}}_{{"t"}}, " ", eta^{{"t"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qr}}_{{qs}}, " ", "VTt"^{{st}}_{{rt}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"p"}}_{{"t"}}, " ", eta^{{"t"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{rt}}, " ", "VTt"^{{st}}_{{ps}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"r"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{st}}_{{qs}}, " ", "VTt"^{{qr}}_{{rt}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"r"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{rt}}, " ", "VTt"^{{st}}_{{ps}}, " ", g^{{"r"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{st}}_{{qs}}, " ", "VTt"^{{qr}}_{{rt}}, " ", g^{{"r"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{pt}}, " ", "VTt"^{{st}}_{{rs}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"p"}}_{{"t"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{qr}}_{{qt}}, " ", "VTt"^{{st}}_{{rs}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"p"}}_{{"t"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ "VTs"^{{pr}}_{{pt}}, " ", "VTt"^{{st}}_{{rs}}, " ", g^{{"p"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{qr}}_{{qt}}, " ", "VTt"^{{st}}_{{rs}}, " ", g^{{"p"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ "VTs"^{{pr}}_{{ps}}, " ", "VTt"^{{ps}}_{{pr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"p"}}, " ", eta^{{"p"}}_{{"s"}}}], "+", RowBox[{ "VTs"^{{ps}}_{{pr}}, " ", "VTt"^{{pr}}_{{ps}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"q"}}, " ", eta^{{"p"}}_{{"s"}}}], "+", RowBox[{ "VTs"^{{qr}}_{{qs}}, " ", "VTt"^{{qs}}_{{qr}}, " ", g^{{ ...[TRUNC 13281]
```


#### [733] Simple Code for all these contractions
**Input [739]**
```
[ RowBox[{"CodeTerm", "[", RowBox[{"spinsec8", "[", RowBox[{"[", "1", "]"}], "]"}], "]"}]], "Input"
```
**Output [742]**
```
CellGroupData[{ Cell[[{ RowBox[{"Length", "[", "HasQ", "]"}], " ", RowBox[{"Length", "[", "HasNoQ", "]"}]}], "Input"
```
**Input [743]**
```
[{ RowBox[{"Length", "[", "HasQ", "]"}], " ", RowBox[{"Length", "[", "HasNoQ", "]"}]}], "Input"
```
**Output [744]**
```
["104"], "Output"
```
**Output [745]**
```
["0"], "Output"
```
**Output [746]**
```
CellGroupData[{ Cell[[{ RowBox[{"Length", "[", "HasQT", "]"}], " ", RowBox[{"Length", "[", "HasQNoT", "]"}], " ", RowBox[{"Length", "[", "HasTNoQ", "]"}]}], "Input"], Cell[["64"], "Output"
```
**Input [747]**
```
[{ RowBox[{"Length", "[", "HasQT", "]"}], " ", RowBox[{"Length", "[", "HasQNoT", "]"}], " ", RowBox[{"Length", "[", "HasTNoQ", "]"}]}], "Input"
```
**Output [748]**
```
["64"], "Output"
```
**Output [749]**
```
["40"], "Output"
```
**Output [750]**
```
["0"], "Output"
```
**Output [755]**
```
CellGroupData[{ Cell[[ RowBox[{" ", RowBox[{ RowBox[{ RowBox[{"Map", "[", CodeTermHasQ,HasQTandPinVs, "]"}], "//", "MatrixForm"}], " ", RowBox[{ RowBox[{"Map", "[", CodeTermHasP,HasQTandQinVs, "]"}], "//", "MatrixForm"}], " ", RowBox[{ RowBox[{"Map", "[", CodeTerm,HasQNoT, "]"}], "//", "MatrixForm"}], " ", RowBox[{ RowBox[{"Map", "[", CodeTerm,HasTNoQ, "]"}], "//", "MatrixForm"}], " ", RowBox[{ RowBox[{"Map", "[", CodeTerm,HasNoTNoQ, "]"}], "//", "MatrixForm"}]}]}]], "Input"
```
**Output [757]**
```
[ TagBox[ RowBox[{"(", "", TagBox[GridBox[{ {"\<\"Oute(p,p) += -1.*Vi_ss[p*n3+r*n2+s*n+t]*g0(d2,ed2,p,r,s,t,params[\ \\\"CohG0\\\"])*Vi_ss[s*n3+t*n2+p*n+r]*Rho(t,s)*Rho(s,t)*eta(r)*;\"\>"}, {"\<\"Oute(p,p) += \ 1.*Vi[p*n3+r*n2+s*n+t]*g0(d2,ed2,p,r,s,t,params[\\\"CohG0\\\"])*Vi[s*n3+t*n2+\ p*n+r]*rho(s)*rho(t)*eta(r)*;\"\>"}, {"\<\"Oute(p,p) += \ 1.*Vi[p*n3+r*n2+t*n+s]*g0(d2,ed2,p,r,t,s,params[\\\"CohG0\\\"])*Vi[t*n3+s*n2+\ p*n+r]*rho(s)*rho(t)*eta(r)*;\"\>"}, {"\<\"Oute(p,p) += \ 1.*Vi_ss[p*n3+r*n2+s*n+t]*g0(d2,ed2,p,r,s,t,params[\\\"CohG0\\\"])*Vi_ss[s*n3+\ t*n2+p*n+r]*rho(s)*rho(t)*eta(r)*;\"\>"}, {"\<\"Oute(p,s) += \ 2.*Vi_ss[p*n3+r*n2+p*n+t]*g0(d2,ed2,p,r,p,t,params[\\\"CohG0\\\"])*Vi_ss[s*n3+\ t*n2+r*n+s]*Rho(t,s)*Rho(p,t)*eta(r)*;\"\>"}, {"\<\"Oute(p,t) += \ 2.*Vi[p*n3+r*n2+p*n+s]*g0(d2,ed2,p,r,p,s,params[\\\"CohG0\\\"])*Vi[t*n3+s*n2+\ t*n+r]*rho(s)*Rho(p,t)*eta(r)*;\"\>"}, {"\<\"Oute(p,t) += \ 2.*Vi_ss[p*n3+r*n2+p*n+s]*g0(d2,ed2,p,r,p,s,params[\\\"CohG0\\\"])*Vi_ss[s*n3+\ t*n2+r*n+t]*rho(s)*Rho(p,t)*eta(r)*;\"\>"}, {"\<\"Oute(p,p) += -2.*Vi_ss[p*n3+r*n2+r*n+t]*g0(d2,ed2,p,r,r,t,params[\ \\\"CohG0\\\"])*Vi_ss[s*n3+t*n2+p*n+s]*Rho(t,s)*Rho(r,t)*Eta(s,r)*;\"\>"}, {"\<\"Oute(p,p) += \ 2.*Vi[p*n3+r*n2+t*n+r]*g0(d2,ed2,p,r,t,r,params[\\\"CohG0\\\"])*Vi[t*n3+s*n2+\ p*n+s]*Rho(r,s)*rho(t)*Eta(s,r)*;\"\>"}, {"\<\"Oute(p,p) += \ 2.*Vi_ss[p*n3+r*n2+r*n+t]*g0(d2,ed2,p,r,r,t,params[\\\"CohG0\\\"])*Vi_ss[s*n3+\ t*n2+p*n+s]*Rho(r,s)*rho(t)*Eta(s,r)*;\"\>"}, {"\<\"Oute(p,r) += -2.*Vi_ss[p*n3+r*n2+p*n+t]*g0(d2,ed2,p,r,p,t,params[\ \\\"CohG0\\\"])*Vi_ss[s*n3+t*n2+r*n+s]*Rho(t,s)*Rho(p,t)*Eta(s,r)*;\"\>"}, {"\<\"Oute(p,r) += \ 2.*Vi_ss[p*n3+r*n2+p*n+t]*g0(d2,ed2,p,r,p,t,params[\\\"CohG0\\\"])*Vi_ss[s*n3+\ t*n2+r*n+s]*Rho(p,s)*rho(t)*Eta(s,r)*;\"\>"}, {"\<\"Outr(p,t) += -2.*Vi_ss[p*n3+t*n2+p*n+s]*g0(d2,ed2,p,t,p,s,params[\ \\\"CohG0\\\"])*Vi_ss[r*n3+s*n2+r*n+t]*Rho(r,t)*Eta(s,r)*Eta(p,s)*;\"\>"}, {"\<\"Outr(p,r) += \ 2.*Vi_ss[p*n3+t*n2+p*n+s]*g0(d2,ed2,p,t,p,s,params[\\\"CohG0\\\"])*Vi_ss[r*n3+\ s*n2+r*n+t]*rho(t)*Eta(s,r)*Eta(p,s)*;\"\>"}, {"\<\"Outr(p,t) += -1.*Vi_ss[p*n3+t*n2+r*n+s]*g0(d2,ed2,p,t,r,s,params[\ \\\"CohG0\\\"])*Vi_ss[r*n3+s*n2+p*n+t]*Rho(p,t)*Eta(s,r)*Eta(r,s)*;\"\>"}, {"\<\"Outr(p,p) += \ 1.*Vi_ss[p*n3+t*n2+r*n+s]*g0(d2,ed2,p,t,r,s,params[\\\"CohG0\\\"])*Vi_ss[r*n3+\ s*n2+p*n+t]*rho(t)*Eta(s,r)*Eta(r,s)*;\"\>"}, {"\<\"Outr(p,t) += \ 1.*Vi_ss[p*n3+t*n2+r*n+s]*g0(d2,ed2,p,t,r,s,params[\\\"CohG0\\\"])*Vi_ss[r*n3+\ s*n2+p*n+t]*Rho(p,t)*eta(r)*eta(s)*;\"\>"}, {"\<\"Outr(p,p) += \ -1.*Vi[p*n3+t*n2+r*n+s]*g0(d2,ed2,p,t,r,s,params[\\\"CohG0\\\"])*Vi[r*n3+s*n2+\ p*n+t]*rho(t)*eta(r)*eta(s)*;\"\>"}, {"\<\"Outr(p,p) += \ -1.*Vi[p*n3+t*n2+s*n+r]*g0(d2,ed2,p,t,s,r,params[\\\"CohG0\\\"])*Vi[s*n3+r*n2+\ p*n+t]*rho(t)*eta(r)*eta(s)*;\"\>"}, {"\<\"Outr(p,p) += -1.*Vi_ss[p*n3+t*n2+r*n+s]*g0(d2,ed2,p,t,r,s,params[\ \\\"CohG0\\\"])*Vi_ss[r*n3+s*n2+p*n+t]*rho(t)*eta(r)*eta(s)*;\"\>"}, {"\<\"Outr(p,s) += \ 2.*Vi_ss[p*n3+t*n2+r*n+t]*g0(d2,ed2,p,t,r,t,params[\\\"CohG0\\\"])*Vi_ss[r*n3+\ s*n2+p*n+s]*Rho(p,t)*eta(r)*Eta(t,s)*;\"\>"}, {"\<\"Outr(p,p) += \ -2.*Vi[p*n3+t*n2+r*n+t]*g0(d2,ed2,p,t,r,t,params[\\\"CohG0\\\"])*Vi[r*n3+s*n2+\ p*n+s]*Rho(s,t)*eta(r)*Eta(t,s)*;\"\>"}, {"\<\"Outr(p,p) += -2.*Vi_ss[p*n3+t*n2+r*n+t]*g0(d2,ed2,p,t,r,t,params[\ \\\"CohG0\\\"])*Vi_ss[r*n3+s*n2+p*n+s]*Rho(s,t)*eta(r)*Eta(t,s)*;\"\>"}, {"\<\"Outr(p,r) += \ 2.*Vi_ss[p*n3+t*n2+s*n+t]*g0(d2,ed2,p,t,s,t,params[\\\"CohG0\\\"])*Vi_ss[r*n3+\ s*n2+p*n+r]*Rho(p,t)*Eta(s,r)*Eta(t,s)*;\"\>"}, {"\<\"Outr(p,p) += -2.*Vi_ss[p*n3+t*n2+s*n+t]*g0(d2,ed2,p,t,s,t,params[\ \\\"CohG0\\\"])*Vi_ss[r*n3+s*n2+p*n+r]*Rho(r,t)*Eta(s,r)*Eta(t,s)*;\"\>"}, {"\<\"Oute(p,r) += -2.*Vi_ss[p*n3+t*n2+s*n+t]*g0(d2,ed2,p,t,s,t,params[\ \\\"CohG0\\\"])*Vi_ss[r*n3+s*n2+p*n+r]*Rho(s,r)*Rho(t,s)*Eta(p,t)*;\"\>"}, {"\<\"Oute(p,s) += -2.*Vi_ss[p*n3+t*n2+r*n+t]*g0(d2,ed2,p,t,r,t,params[\ \\\"CohG0\\\"])*Vi_ss[r*n3+s*n2+p*n+s]*rho(r)*Rho(t,s)*Eta(p,t)*;\"\>"}, {"\<\"Oute(p,t) += \ 1.*Vi_ss[p*n3+t*n2+r*n+s]*g0(d2,ed2,p,t,r,s,params[\\\"CohG0\\\"])*Vi_ss[r*n3+\ s*n2+p*n+t]*Rho(s,r)*Rho(r,s)*Eta(p,t)*;\"\>"}, {"\<\"Oute(p,t) += -1.*Vi_ss[p*n3+t*n2+r*n+s]*g0(d2,ed2,p,t,r,s,params[\ \\\"CohG0\\\"])*Vi_ss[r*n3+s*n2+p*n+t]*rho(r)*rho(s)*Eta(p,t)*;\"\>"}, {"\<\"Outr(p,t) += \ -2.*Vi[p*n3+s*n2+p*n+r]*g0(d2,ed2,p,s,p,r,params[\\\"CohG0\\\"])*Vi[t*n3+r*n2+\ t*n+s]*rho(s)*eta(r)*Eta(p,t)*;\"\>"}, {"\<\"Outr(p,t) += -2.*Vi_ss[p*n3+s*n2+p*n+r]*g0(d2,ed2,p,s,p,r,params[\ \\\"CohG0\\\"])*Vi_ss[r*n3+t*n2+s*n+t]*rho(s)*eta(r)*Eta(p,t)*;\"\>"}, {"\<\"Outr(p,s) += \ 2.*Vi_ss[p*n3+s*n2+p*n+r]*g0(d2,ed2,p,s,p,r,params[\\\"CohG0\\\"])*Vi_ss[r*n3+\ t*n2+s*n+t]*Rho(t,s)*eta(r)*Eta(p,t)*;\"\>"} }, GridBoxAlignment->{ "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}}, GridBoxSpacings->{"Columns" -> { Offset[0.27999999999999997`], { Offset[0.5599999999999999]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> { Offset[0.2], { Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}], Column], "", ")"}], Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]], "Output"
```
**Output [758]**
```
[ TagBox[ RowBox[{"(", "", TagBox[GridBox[{ {"\<\"Oute2(q,q) += \ -1.*Vi_ss[s*n3+t*n2+q*n+r]*g0(d2,ed2,s,t,q,r,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+r*n2+s*n+t]*Rho(t,s)*Rho(s,t)*eta(r)*;\"\>"}, {"\<\"Oute2(q,q) += \ 1.*Vi[s*n3+t*n2+q*n+r]*g0(d2,ed2,s,t,q,r,params[\\\"CohG0\\\"])*Vi[q*n3+r*n2+\ s*n+t]*rho(s)*rho(t)*eta(r)*;\"\>"}, {"\<\"Oute2(q,q) += \ 1.*Vi[t*n3+s*n2+q*n+r]*g0(d2,ed2,t,s,q,r,params[\\\"CohG0\\\"])*Vi[q*n3+r*n2+\ t*n+s]*rho(s)*rho(t)*eta(r)*;\"\>"}, {"\<\"Oute2(q,q) += 1.*Vi_ss[s*n3+t*n2+q*n+r]*g0(d2,ed2,s,t,q,r,params[\ \\\"CohG0\\\"])*Vi_ss[q*n3+r*n2+s*n+t]*rho(s)*rho(t)*eta(r)*;\"\>"}, {"\<\"Outr2(t,q) += \ -2.*Vi_ss[q*n3+r*n2+q*n+t]*g0(d2,ed2,q,r,q,t,params[\\\"CohG0\\\"])*Vi_ss[s*\ n3+t*n2+r*n+s]*Rho(t,s)*Eta(s,q)*eta(r)*;\"\>"}, {"\<\"Outr2(t,q) += -2.*Vi[q*n3+r*n2+q*n+s]*g0(d2,ed2,q,r,q,s,params[\\\ \"CohG0\\\"])*Vi[t*n3+s*n2+t*n+r]*rho(s)*Eta(t,q)*eta(r)*;\"\>"}, {"\<\"Outr2(t,q) += \ -2.*Vi_ss[q*n3+r*n2+q*n+s]*g0(d2,ed2,q,r,q,s,params[\\\"CohG0\\\"])*Vi_ss[s*\ n3+t*n2+r*n+t]*rho(s)*Eta(t,q)*eta(r)*;\"\>"}, {"\<\"Oute2(q,q) += \ -2.*Vi_ss[s*n3+t*n2+q*n+s]*g0(d2,ed2,s,t,q,s,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+r*n2+r*n+t]*Rho(t,s)*Rho(r,t)*Eta(s,r)*;\"\>"}, {"\<\"Oute2(q,q) += \ 2.*Vi[t*n3+s*n2+q*n+s]*g0(d2,ed2,t,s,q,s,params[\\\"CohG0\\\"])*Vi[q*n3+r*n2+\ t*n+r]*Rho(r,s)*rho(t)*Eta(s,r)*;\"\>"}, {"\<\"Oute2(q,q) += 2.*Vi_ss[s*n3+t*n2+q*n+s]*g0(d2,ed2,s,t,q,s,params[\ \\\"CohG0\\\"])*Vi_ss[q*n3+r*n2+r*n+t]*Rho(r,s)*rho(t)*Eta(s,r)*;\"\>"}, {"\<\"Outr2(t,q) += 2.*Vi_ss[q*n3+r*n2+q*n+t]*g0(d2,ed2,q,r,q,t,params[\ \\\"CohG0\\\"])*Vi_ss[s*n3+t*n2+r*n+s]*Rho(t,s)*Eta(r,q)*Eta(s,r)*;\"\>"}, {"\<\"Outr2(s,q) += \ -2.*Vi_ss[q*n3+r*n2+q*n+t]*g0(d2,ed2,q,r,q,t,params[\\\"CohG0\\\"])*Vi_ss[s*\ n3+t*n2+r*n+s]*rho(t)*Eta(r,q)*Eta(s,r)*;\"\>"}, {"\<\"Oute2(s,q) += 2.*Vi_ss[q*n3+t*n2+q*n+s]*g0(d2,ed2,q,t,q,s,params[\ \\\"CohG0\\\"])*Vi_ss[r*n3+s*n2+r*n+t]*Rho(t,q)*Rho(r,t)*Eta(s,r)*;\"\>"}, {"\<\"Oute2(s,q) += \ -2.*Vi_ss[q*n3+t*n2+q*n+s]*g0(d2,ed2,q,t,q,s,params[\\\"CohG0\\\"])*Vi_ss[r*\ n3+s*n2+r*n+t]*Rho(r,q)*rho(t)*Eta(s,r)*;\"\>"}, {"\<\"Outr2(t,q) += \ -1.*Vi_ss[r*n3+s*n2+q*n+t]*g0(d2,ed2,r,s,q,t,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+t*n2+r*n+s]*Rho(t,q)*Eta(s,r)*Eta(r,s)*;\"\>"}, {"\<\"Outr2(q,q) += 1.*Vi_ss[r*n3+s*n2+q*n+t]*g0(d2,ed2,r,s,q,t,params[\ \\\"CohG0\\\"])*Vi_ss[q*n3+t*n2+r*n+s]*rho(t)*Eta(s,r)*Eta(r,s)*;\"\>"}, {"\<\"Outr2(t,q) += 1.*Vi_ss[r*n3+s*n2+q*n+t]*g0(d2,ed2,r,s,q,t,params[\ \\\"CohG0\\\"])*Vi_ss[q*n3+t*n2+r*n+s]*Rho(t,q)*eta(r)*eta(s)*;\"\>"}, {"\<\"Outr2(q,q) += -1.*Vi[r*n3+s*n2+q*n+t]*g0(d2,ed2,r,s,q,t,params[\\\ \"CohG0\\\"])*Vi[q*n3+t*n2+r*n+s]*rho(t)*eta(r)*eta(s)*;\"\>"}, {"\<\"Outr2(q,q) += -1.*Vi[s*n3+r*n2+q*n+t]*g0(d2,ed2,s,r,q,t,params[\\\ \"CohG0\\\"])*Vi[q*n3+t*n2+s*n+r]*rho(t)*eta(r)*eta(s)*;\"\>"}, {"\<\"Outr2(q,q) += \ -1.*Vi_ss[r*n3+s*n2+q*n+t]*g0(d2,ed2,r,s,q,t,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+t*n2+r*n+s]*rho(t)*eta(r)*eta(s)*;\"\>"}, {"\<\"Outr2(t,q) += 2.*Vi_ss[r*n3+s*n2+q*n+s]*g0(d2,ed2,r,s,q,s,params[\ \\\"CohG0\\\"])*Vi_ss[q*n3+t*n2+r*n+t]*Rho(s,q)*eta(r)*Eta(t,s)*;\"\>"}, {"\<\"Outr2(q,q) += -2.*Vi[r*n3+s*n2+q*n+s]*g0(d2,ed2,r,s,q,s,params[\\\ \"CohG0\\\"])*Vi[q*n3+t*n2+r*n+t]*Rho(s,t)*eta(r)*Eta(t,s)*;\"\>"}, {"\<\"Outr2(q,q) += \ -2.*Vi_ss[r*n3+s*n2+q*n+s]*g0(d2,ed2,r,s,q,s,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+t*n2+r*n+t]*Rho(s,t)*eta(r)*Eta(t,s)*;\"\>"}, {"\<\"Outr2(t,q) += 2.*Vi_ss[r*n3+s*n2+q*n+r]*g0(d2,ed2,r,s,q,r,params[\ \\\"CohG0\\\"])*Vi_ss[q*n3+t*n2+s*n+t]*Rho(r,q)*Eta(s,r)*Eta(t,s)*;\"\>"}, {"\<\"Outr2(q,q) += \ -2.*Vi_ss[r*n3+s*n2+q*n+r]*g0(d2,ed2,r,s,q,r,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+t*n2+s*n+t]*Rho(r,t)*Eta(s,r)*Eta(t,s)*;\"\>"}, {"\<\"Oute2(t,q) += \ -2.*Vi_ss[r*n3+s*n2+q*n+r]*g0(d2,ed2,r,s,q,r,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+t*n2+s*n+t]*Rho(s,r)*Rho(t,s)*Eta(r,q)*;\"\>"}, {"\<\"Oute2(t,q) += \ -2.*Vi_ss[r*n3+s*n2+q*n+s]*g0(d2,ed2,r,s,q,s,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+t*n2+r*n+t]*rho(r)*Rho(t,s)*Eta(s,q)*;\"\>"}, {"\<\"Oute2(t,q) += 1.*Vi_ss[r*n3+s*n2+q*n+t]*g0(d2,ed2,r,s,q,t,params[\ \\\"CohG0\\\"])*Vi_ss[q*n3+t*n2+r*n+s]*Rho(s,r)*Rho(r,s)*Eta(t,q)*;\"\>"}, {"\<\"Oute2(t,q) += \ -1.*Vi_ss[r*n3+s*n2+q*n+t]*g0(d2,ed2,r,s,q,t,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+t*n2+r*n+s]*rho(r)*rho(s)*Eta(t,q)*;\"\>"}, {"\<\"Oute2(t,q) += \ 2.*Vi[q*n3+s*n2+q*n+r]*g0(d2,ed2,q,s,q,r,params[\\\"CohG0\\\"])*Vi[t*n3+r*n2+\ t*n+s]*Rho(t,q)*rho(s)*eta(r)*;\"\>"}, {"\<\"Oute2(t,q) += 2.*Vi_ss[q*n3+s*n2+q*n+r]*g0(d2,ed2,q,s,q,r,params[\ \\\"CohG0\\\"])*Vi_ss[r*n3+t*n2+s*n+t]*Rho(t,q)*rho(s)*eta(r)*;\"\>"}, {"\<\"Oute2(t,q) += \ -2.*Vi_ss[q*n3+s*n2+q*n+r]*g0(d2,ed2,q,s,q,r,params[\\\"CohG0\\\"])*Vi_ss[r*\ n3+t*n2+s*n+t]*Rho(s,q)*Rho(t,s)*eta(r)*;\"\>"} }, GridBoxAlignment->{ "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}}, GridBoxSpacings->{"Columns" -> { Offset[0.27999999999999997`], { Offset[0.5599999999999999]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> { Offset[0.2], { Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}], Column], "", ")"}], Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]], "Output"
```
**Output [759]**
```
[ TagBox[ RowBox[{"(", "", TagBox[GridBox[{ {"\<\"Out(p,q) += \ -2.*Vi[p*n3+r*n2+p*n+s]*g0(d2,ed2,p,r,p,s,params[\\\"CohG0\\\"])*Vi[p*n3+s*n2+\ p*n+r]*rho(p)*rho(s)*Eta(p,q)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi[p*n3+r*n2+s*n+r]*g0(d2,ed2,p,r,s,r,params[\\\"CohG0\\\"])*Vi[s*n3+r*n2+\ p*n+r]*rho(r)*rho(s)*Eta(p,q)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi_ss[p*n3+r*n2+p*n+s]*g0(d2,ed2,p,r,p,s,params[\\\"CohG0\\\"])*Vi_ss[p*\ n3+s*n2+p*n+r]*rho(p)*rho(s)*Eta(p,q)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi_ss[p*n3+r*n2+p*n+s]*g0(d2,ed2,p,r,p,s,params[\\\"CohG0\\\"])*Vi_ss[p*\ n3+s*n2+p*n+r]*Rho(p,q)*rho(r)*Eta(s,p)*Eta(p,s)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi_ss[p*n3+r*n2+p*n+s]*g0(d2,ed2,p,r,p,s,params[\\\"CohG0\\\"])*Vi_ss[p*\ n3+s*n2+p*n+r]*Rho(s,p)*Rho(p,s)*Eta(r,q)*Eta(p,r)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi_ss[p*n3+r*n2+r*n+s]*g0(d2,ed2,p,r,r,s,params[\\\"CohG0\\\"])*Vi_ss[r*\ n3+s*n2+p*n+r]*Rho(p,q)*rho(r)*Eta(s,r)*Eta(r,s)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi_ss[p*n3+r*n2+r*n+s]*g0(d2,ed2,p,r,r,s,params[\\\"CohG0\\\"])*Vi_ss[r*\ n3+s*n2+p*n+r]*rho(r)*rho(s)*Eta(p,q)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi_ss[p*n3+r*n2+r*n+s]*g0(d2,ed2,p,r,r,s,params[\\\"CohG0\\\"])*Vi_ss[r*\ n3+s*n2+p*n+r]*Rho(s,r)*Rho(r,s)*Eta(r,q)*Eta(p,r)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi_ss[p*n3+s*n2+p*n+r]*g0(d2,ed2,p,s,p,r,params[\\\"CohG0\\\"])*Vi_ss[p*\ n3+r*n2+p*n+s]*Rho(s,q)*Rho(p,s)*eta(p)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi_ss[p*n3+s*n2+r*n+s]*g0(d2,ed2,p,s,r,s,params[\\\"CohG0\\\"])*Vi_ss[r*\ n3+s*n2+p*n+s]*Rho(s,q)*Rho(p,s)*eta(r)*eta(s)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi_ss[q*n3+r*n2+q*n+s]*g0(d2,ed2,q,r,q,s,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+s*n2+q*n+r]*Rho(s,q)*Rho(p,s)*eta(q)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi[q*n3+s*n2+q*n+r]*g0(d2,ed2,q,s,q,r,params[\\\"CohG0\\\"])*Vi[q*n3+r*n2+\ q*n+s]*rho(q)*rho(s)*Eta(p,q)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi_ss[q*n3+s*n2+q*n+r]*g0(d2,ed2,q,s,q,r,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+r*n2+q*n+s]*Rho(p,q)*rho(r)*Eta(s,q)*Eta(q,s)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi_ss[q*n3+s*n2+q*n+r]*g0(d2,ed2,q,s,q,r,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+r*n2+q*n+s]*rho(q)*rho(s)*Eta(p,q)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi_ss[q*n3+s*n2+q*n+r]*g0(d2,ed2,q,s,q,r,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+r*n2+q*n+s]*Rho(s,q)*Rho(q,s)*Eta(r,q)*Eta(p,r)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi_ss[r*n3+s*n2+q*n+r]*g0(d2,ed2,r,s,q,r,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+r*n2+r*n+s]*Rho(p,q)*rho(r)*Eta(s,r)*Eta(r,s)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi_ss[r*n3+s*n2+q*n+r]*g0(d2,ed2,r,s,q,r,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+r*n2+r*n+s]*rho(r)*rho(s)*Eta(p,q)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi_ss[r*n3+s*n2+q*n+r]*g0(d2,ed2,r,s,q,r,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+r*n2+r*n+s]*Rho(s,r)*Rho(r,s)*Eta(r,q)*Eta(p,r)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi_ss[r*n3+s*n2+q*n+s]*g0(d2,ed2,r,s,q,s,params[\\\"CohG0\\\"])*Vi_ss[q*\ n3+s*n2+r*n+s]*Rho(s,q)*Rho(p,s)*eta(r)*eta(s)*;\"\>"}, {"\<\"Out(p,q) += \ -2.*Vi[s*n3+r*n2+q*n+r]*g0(d2,ed2,s,r,q,r,params[\\\"CohG0\\\"])*Vi[q*n3+r*n2+\ s*n+r]*rho(r)*rho(s)*Eta(p,q)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += 2.*Vi_ss[p*n3+r*n2+p*n+s]*g0(d2,ed2,p,r,p,s,params[\\\ \"CohG0\\\"])*Vi_ss[p*n3+s*n2+p*n+r]*Rho(s,p)*Rho(p,s)*Eta(p,q)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += 2.*Vi_ss[p*n3+r*n2+r*n+s]*g0(d2,ed2,p,r,r,s,params[\\\ \"CohG0\\\"])*Vi_ss[r*n3+s*n2+p*n+r]*Rho(r,q)*Rho(p,r)*Eta(s,r)*Eta(r,s)*;\"\>\ "}, {"\<\"Out(p,q) += 2.*Vi_ss[p*n3+r*n2+r*n+s]*g0(d2,ed2,p,r,r,s,params[\\\ \"CohG0\\\"])*Vi_ss[r*n3+s*n2+p*n+r]*Rho(s,r)*Rho(r,s)*Eta(p,q)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += \ 2.*Vi[p*n3+s*n2+p*n+r]*g0(d2,ed2,p,s,p,r,params[\\\"CohG0\\\"])*Vi[p*n3+r*n2+\ p*n+s]*Rho(p,q)*rho(s)*eta(p)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += \ 2.*Vi[p*n3+s*n2+r*n+s]*g0(d2,ed2,p,s,r,s,params[\\\"CohG0\\\"])*Vi[r*n3+s*n2+\ p*n+s]*Rho(p,q)*rho(s)*eta(r)*eta(s)*;\"\>"}, {"\<\"Out(p,q) += 2.*Vi_ss[p*n3+s*n2+p*n+r]*g0(d2,ed2,p,s,p,r,params[\\\ \"CohG0\\\"])*Vi_ss[p*n3+r*n2+p*n+s]*rho(p)*rho(r)*Eta(s,q)*Eta(p,s)*;\"\>"}, {"\<\"Out(p,q) += 2.*Vi_ss[p*n3+s*n2+p*n+r]*g0(d2,ed2,p,s,p,r,params[\\\ \"CohG0\\\"])*Vi_ss[p*n3+r*n2+p*n+s]*Rho(p,q)*rho(s)*eta(p)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += 2.*Vi_ss[p*n3+s*n2+p*n+r]*g0(d2,ed2,p,s,p,r,params[\\\ \"CohG0\\\"])*Vi_ss[p*n3+r*n2+p*n+s]*Rho(s,q)*Rho(p,s)*Eta(r,p)*Eta(p,r)*;\"\>\ "}, {"\<\"Out(p,q) += 2.*Vi_ss[p*n3+s*n2+r*n+s]*g0(d2,ed2,p,s,r,s,params[\\\ \"CohG0\\\"])*Vi_ss[r*n3+s*n2+p*n+s]*Rho(p,q)*rho(s)*eta(r)*eta(s)*;\"\>"}, {"\<\"Out(p,q) += 2.*Vi_ss[p*n3+s*n2+r*n+s]*g0(d2,ed2,p,s,r,s,params[\\\ \"CohG0\\\"])*Vi_ss[r*n3+s*n2+p*n+s]*rho(r)*rho(s)*Eta(s,q)*Eta(p,s)*;\"\>"}, {"\<\"Out(p,q) += \ 2.*Vi[q*n3+r*n2+q*n+s]*g0(d2,ed2,q,r,q,s,params[\\\"CohG0\\\"])*Vi[q*n3+s*n2+\ q*n+r]*Rho(p,q)*rho(s)*eta(q)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += 2.*Vi_ss[q*n3+r*n2+q*n+s]*g0(d2,ed2,q,r,q,s,params[\\\ \"CohG0\\\"])*Vi_ss[q*n3+s*n2+q*n+r]*Rho(p,q)*rho(s)*eta(q)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += 2.*Vi_ss[q*n3+r*n2+q*n+s]*g0(d2,ed2,q,r,q,s,params[\\\ \"CohG0\\\"])*Vi_ss[q*n3+s*n2+q*n+r]*rho(q)*rho(r)*Eta(s,q)*Eta(p,s)*;\"\>"}, {"\<\"Out(p,q) += 2.*Vi_ss[q*n3+r*n2+q*n+s]*g0(d2,ed2,q,r,q,s,params[\\\ \"CohG0\\\"])*Vi_ss[q*n3+s*n2+q*n+r]*Rho(s,q)*Rho(p,s)*Eta(r,q)*Eta(q,r)*;\"\>\ "}, {"\<\"Out(p,q) += 2.*Vi_ss[q*n3+s*n2+q*n+r]*g0(d2,ed2,q,s,q,r,params[\\\ \"CohG0\\\"])*Vi_ss[q*n3+r*n2+q*n+s]*Rho(s,q)*Rho(q,s)*Eta(p,q)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += \ 2.*Vi[r*n3+s*n2+q*n+s]*g0(d2,ed2,r,s,q,s,params[\\\"CohG0\\\"])*Vi[q*n3+s*n2+\ r*n+s]*Rho(p,q)*rho(s)*eta(r)*eta(s)*;\"\>"}, {"\<\"Out(p,q) += 2.*Vi_ss[r*n3+s*n2+q*n+r]*g0(d2,ed2,r,s,q,r,params[\\\ \"CohG0\\\"])*Vi_ss[q*n3+r*n2+r*n+s]*Rho(r,q)*Rho(p,r)*Eta(s,r)*Eta(r,s)*;\"\>\ "}, {"\<\"Out(p,q) += 2.*Vi_ss[r*n3+s*n2+q*n+r]*g0(d2,ed2,r,s,q,r,params[\\\ \"CohG0\\\"])*Vi_ss[q*n3+r*n2+r*n+s]*Rho(s,r)*Rho(r,s)*Eta(p,q)*eta(r)*;\"\>"}, {"\<\"Out(p,q) += 2.*Vi_ss[r*n3+s*n2+q*n+s]*g0(d2,ed2,r,s,q,s,params[\\\ \"CohG0\\\"])*Vi_ss[q*n3+s*n2+r*n+s]*Rho(p,q)*rho(s)*eta(r)*eta(s)*;\"\>"}, {"\<\"Out(p,q) += 2.*Vi ...[TRUNC 6548]
```
**Output [760]**
```
[ TagBox[ {}, Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]], "Output"
```
**Output [761]**
```
[ TagBox[ {}, Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]], "Output"
```
> Toss out all non - oscillating terms :


### [764] Double Secular version.
**Input [765]**
```
[ RowBox[{ RowBox[{ RowBox[{"{", RowBox[{"{", {\<VTs\>2},{\<t\>\<u\>},{\<r\>\<s\>}, "}"}], "}"}], "]"}], ",", RowBox[{"{", RowBox[{"{", {\<VTt\>2},{\<x\>\<y\>},{\<v\>\<w\>}, "}"}], "}"}]}]], "Input", EmphasizeSyntaxErrors->True
```
**Output [766]**
```
CellGroupData[{ Cell[[{ RowBox[{ RowBox[{ RowBox[{"STwo", "[", p_, q_, r_, s_, "]"}], " ", ":=", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {r}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {s}}], "}"}]}], "}"}], ",", " ", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {r}}], "}"}]}], "}"}], ",", " ", RowBox[{"{", RowBox[{ {\<Sign\>-1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {r}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {r}}], "}"}]}], "}"}]}], "}"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"ByHand3", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"DistributeOperatorProduct", "[", RowBox[{ RowBox[{"DistributeOperatorProduct", "[", RowBox[{"FromCode", ",", RowBox[{"STwo", "[", "\<t\>","\<v\>","\<u\>","\<w\>", "]"}]}], "]"}], ",", RowBox[{"STwo", "[", "\<x\>","\<r\>","\<y\>","\<s\>", "]"}]}], "]"}]}], "]"}]}], "]"}]}], "]"}]}], "Input"
```
**Output [768]**
```
[ RowBox[{ RowBox[{ -(1)/(2), " ", "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{ps}}_{{pr}}, " ", "VTt"^{{pr}}_{{ps}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{rs}}, " ", "VTt"^{{pr}}_{{rs}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{ps}}_{{rt}}, " ", "VTt"^{{pr}}_{{st}}, " ", g^{{"t"}}_{{"q"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"t"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"s"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{qr}}, " ", "VTt"^{{qr}}_{{st}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"s"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{qr}}, " ", "VTt"^{{qr}}_{{st}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{qr}}_{{qs}}, " ", "VTt"^{{qs}}_{{qr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{ps}}_{{ps}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{ps}}_{{st}}, " ", "VTt"^{{pr}}_{{rt}}, " ", g^{{"t"}}_{{"q"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"t"}}_{{"p"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qs}}_{{qs}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{rs}}_{{qr}}, " ", "VTt"^{{rs}}_{{qr}}, " ", g^{{"r"}}_{{"r"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"q"}}_{{"s"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{rs}}, " ", "VTt"^{{rs}}_{{pr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{rs}}_{{qr}}, " ", "VTt"^{{qr}}_{{rs}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pt}}_{{rs}}, " ", "VTt"^{{rs}}_{{pt}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rs}}_{{qt}}, " ", "VTt"^{{qt}}_{{rs}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{rs}}_{{qt}}, " ", "VTt"^{{qt}}_{{rs}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"t"}}_{{"q"}}, " ", eta^{{"p"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rs}}_{{qt}}, " ", "VTt"^{{qt}}_{{rs}}, " ", g^{{"r"}}_{{"r"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"t"}}_{{"q"}}, " ", eta^{{"p"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{qs}}, " ", "VTt"^{{st}}_{{qr}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"p"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"q"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{qr}}, " ", "VTt"^{{st}}_{{qs}}, " ", g^{{"r"}}_{{"s"}}, " ", g^{{"p"}}_{{"t"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"q"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{ps}}_{{st}}, " ", "VTt"^{{rt}}_{{pr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"t"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{qr}}, " ", "VTt"^{{qs}}_{{st}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"t"}}_{{"t"}}}]}]], "Output"
```
**Input [770]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"S2terms", "=", RowBox[{"AssociativeSimplify", "[", ByHand3//.NoZerosDenoms, "]"}]}], "]"}]], "Input"
```
**Output [772]**
```
CellGroupData[{ Cell[[ RowBox[{"Map", "[", NumIndicesInTerm,S2terms, "]"}]], "Input"
```
**Input [773]**
```
[ RowBox[{"Map", "[", NumIndicesInTerm,S2terms, "]"}]], "Input"
```
**Output [774]**
```
[ {44555554544455555555}], "Output"
```
**Input [776]**
```
[ RowBox[{ RowBox[{"TermCodeSec", "[", "ATerm_", "]"}], ":=", RowBox[{"(", " ", RowBox[{ RowBox[{ RowBox[{"fac", "=", RowBox[{"ATerm", "[", RowBox[{"[", 1,2, "]"}], "]"}]}], ";", " ", RowBox[{"ndinds", "=", RowBox[{"NumIndicesInTerm", "[", "ATerm", "]"}]}], ";", " ", RowBox[{"Vs", " ", "=", RowBox[{"ATerm", "[", RowBox[{"[", "2", "]"}], "]"}]}], ";", " ", RowBox[{"Vt", " ", "=", RowBox[{"ATerm", "[", RowBox[{"[", "3", "]"}], "]"}]}], ";", " ", RowBox[{"vsinds", "=", RowBox[{"Map", "[", RowBox[{"TlS", ",", RowBox[{"{", RowBox[{ RowBox[{"Vs", "[", RowBox[{"[", 2,1, "]"}], "]"}], ",", RowBox[{"Vs", "[", RowBox[{"[", 2,2, "]"}], "]"}], ",", RowBox[{"Vs", "[", RowBox[{"[", 3,1, "]"}], "]"}], ",", RowBox[{"Vs", "[", RowBox[{"[", 3,2, "]"}], "]"}]}], "}"}]}], "]"}]}], ";", " ", RowBox[{"vtinds", "=", RowBox[{"Map", "[", RowBox[{"TlS", ",", RowBox[{"{", RowBox[{ RowBox[{"Vt", "[", RowBox[{"[", 2,1, "]"}], "]"}], ",", RowBox[{"Vt", "[", RowBox[{"[", 2,2, "]"}], "]"}], ",", RowBox[{"Vt", "[", RowBox[{"[", 3,1, "]"}], "]"}], ",", RowBox[{"Vt", "[", RowBox[{"[", 3,2, "]"}], "]"}]}], "}"}]}], "]"}]}], ";", " ", RowBox[{"Gammas", "=", RowBox[{"ATerm", "[", RowBox[{"[", 4;;, "]"}], "]"}]}], ";", " ", RowBox[{"ContractableGammas", "=", RowBox[{"Select", "[", RowBox[{"Gammas", ",", RowBox[{ RowBox[{"ContractableWithVtsec", "[", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"#", "[", RowBox[{"[", 2,1, "]"}], "]"}], ",", RowBox[{"#", "[", RowBox[{"[", 3,1, "]"}], "]"}]}], "}"}], ",", "vsinds", ",", "vtinds"}], "]"}], "&"}]}], "]"}]}], ";", " ", RowBox[{"UnContractableGammas", "=", RowBox[{"SlowComplement", "[", Gammas,ContractableGammas, "]"}]}], ";", " ", RowBox[{"UGIndices", "=", RowBox[{"{", RowBox[{ RowBox[{"UnContractableGammas", "[", RowBox[{"[", 1,2,1, "]"}], "]"}], ",", RowBox[{"UnContractableGammas", "[", RowBox[{"[", 1,3,1, "]"}], "]"}]}], "}"}]}], ";", " ", (*Get the loop indices.*), " ", RowBox[{"LoopIndices", "=", RowBox[{"Union", "[", RowBox[{"Join", "[", vsinds,vtinds, "]"}], "]"}]}], ";", " ", RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{"Union", "[", "UGIndices", "]"}], "==", {\<p\>\<q\>}}], ",", RowBox[{"(", " ", RowBox[{ RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{"Count", "[", LoopIndices,"\<p\>", "]"}], ">", "0"}], ",", RowBox[{"(", " ", RowBox[{ ExIntermediateIndex1="\<p\>", ";", " ", ExIntermediateIndex2="\<p\>", ";"}], " ", ")"}], ",", RowBox[{"(", " ", RowBox[{ ExIntermediateIndex1="\<q\>", ";", " ", ExIntermediateIndex2="\<q\>", ";"}], " ", ")"}]}], "]"}], ";"}], " ", ")"}], ",", RowBox[{"(", " ", RowBox[{ RowBox[{"ExIntermediateIndex1", "=", RowBox[{ RowBox[{"Complement", "[", {\<p\>\<q\>},UGIndices, "]"}], "[", RowBox[{"[", "1", "]"}], "]"}]}], ";", " ", RowBox[{"ExIntermediateIndex2", "=", RowBox[{ RowBox[{"Complement", "[", UGIndices,{\<p\>\<q\>}, "]"}], "[", RowBox[{"[", "1", "]"}], "]"}]}], ";"}], " ", ")"}]}], "]"}], ";", " ", RowBox[{"(*", RowBox[{ RowBox[{"Print", "[", "\<eii\>",ExIntermediateIndex1,ExIntermediateIndex2, "]"}], ";", " ", RowBox[{"Print", "[", "\<UCG\>",UnContractableGammas,"\<LoopIndices\>",LoopIndices,"\<ExternalIntIndex:\>",ExIntermediateIndex, "]"}], ";"}], "*)"}], " ", RowBox[{"Gstring1", "=", RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{"ContractableGammas", "[", RowBox[{"[", 1,1,1, "]"}], "]"}], "===", "g"}], ",", RowBox[{"\"\<*Rho(\>\"", "<>", RowBox[{"ContractableGammas", "[", RowBox[{"[", 1,2,1, "]"}], "]"}], "<>", "\"\<,\>\"", "<>", RowBox[{"ContractableGammas", "[", RowBox[{"[", 1,3,1, "]"}], "]"}], "<>", "\"\<)\>\""}], ",", RowBox[{"\"\<*Eta(\>\"", "<>", RowBox[{"ContractableGammas", "[", RowBox[{"[", 1,2,1, "]"}], "]"}], "<>", "\"\<,\>\"", "<>", RowBox[{"ContractableGammas", "[", RowBox[{"[", 1,3,1, "]"}], "]"}], "<>", "\"\<)\>\""}]}], "]"}]}], ";", " ", RowBox[{"Gstring2", "=", RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{"ContractableGammas", "[", RowBox[{"[", 2,1,1, "]"}], "]"}], "===", "g"}], ",", RowBox[{"\"\<*Rho(\>\"", "<>", RowBox[{"ContractableGammas", "[", RowBox[{"[", 2,2,1, "]"}], "]"}], "<>", "\"\<,\>\"", "<>", RowBox[{"ContractableGammas", "[", RowBox[{"[", 2,3,1, "]"}], "]"}], "<>", "\"\<)\>\""}], ",", RowBox[{"\"\<*Eta(\>\"", "<>", RowBox[{"ContractableGammas", "[", RowBox[{"[", 2,2,1, "]"}], "]"}], "<>", "\"\<,\>\"", "<>", RowBox[{"ContractableGammas", "[", RowBox[{"[", 2,3,1, "]"}], "]"}], "<>", "\"\<)\>\""}]}], "]"}]}], ";", " ", RowBox[{"Gstring3", "=", RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{"ContractableGammas", "[", RowBox[{"[", 3,1,1, "]"}], "]"}], "===", "g"}], ",", RowBox[{"\"\<*Rho(\>\"", "<>", RowBox[{"ContractableGammas", "[", RowBox[{"[", 3,2,1, "]"}], "]"}], "<>", "\"\<,\>\"", "<>", RowBox[{"ContractableGammas", "[", RowBox[{"[", 3,3,1, "]"}], "]"}], "<>", "\"\<)\>\""}], ",", RowBox[{"\"\<*Eta(\>\"", "<>", RowBox[{"ContractableGammas", "[", RowBox[{"[", 3,2,1, "]"}], "]"}], "<>", "\"\<,\>\"", "<>", RowBox[{"ContractableGammas", "[", RowBox[{"[", 3,3,1, "]"}], "]"}], "<>", "\"\<)\>\""}]}], "]"}]}], ";", " ", RowBox[{"Gstring4", "=", RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{"UnContractableGammas", "[", RowBox[{"[", 1,1,1, "]"}], "]"}], "===", "g"}], ",", RowBox[{"\"\<*Rho(\>\"", "<>", RowBox[{"UnContractableGammas", "[", RowBox[{"[", 1,2,1, "]"}], "]"}], "<>", "\"\<,\>\"", "<>", RowBox[{"UnContractableGammas", "[", RowBox[{"[", 1,3,1, "]"}], "]"}], "<>", "\"\<)\>\""}], ",", RowBox[{"\"\<*Eta(\>\"", "<>", RowBox[{"UnContractableGammas", "[", RowBox[{"[", 1,2,1, "]"}], "]"}], "<>", "\"\<,\>\"", "<>", RowBox[{"UnContractableGammas", "[", RowBox[{"[", 1,3,1, "]"}], "]"}], "<>", "\"\<)\>\""}]}], "]"}]}], ";", " ", cstring="\<\>", ";", " ", RowBox[{"cstring", "=", RowBox[{ "cstring", "<>", "\"\<I(\>\"", "<>", "ExIntermediateIndex1", "<>", "\"\<,\>\"", "<>", "ExIntermediateIndex2", "<>", "\"\<)+=(\>\"", "<>", RowBox[{"ToString", "[", RowBox[{"N", "[", RowBox[{"ATerm", "[", RowBox[{"[", 1,2, "]"}], "]"}], "]"}], "]"}], "<>", "\"\<)*\>\"", "<>", "\"\<Vs[n3*\>\"", "<>", RowBox[{"vsinds", "[", RowBox[{"[", "1", "]"}], "]" ...[TRUNC 6898]
```
**Input [778]**
```
[ RowBox[{"Map", "[", NumIndicesInTerm,ByHand3, "]"}]], "Input"
```
**Output [780]**
```
CellGroupData[{ Cell[[ RowBox[{"RecognizableForm", "[", RowBox[{"AssociativeSimplify", "[", "ByHand3", "]"}], "]"}]], "Input"
```
**Input [781]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"AssociativeSimplify", "[", "ByHand3", "]"}], "]"}]], "Input"
```
**Output [782]**
```
[ RowBox[{ RowBox[{ -(1)/(2), " ", "VTs"^{{pr}}_{{pr}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{ps}}_{{pr}}, " ", "VTt"^{{pr}}_{{ps}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{rs}}, " ", "VTt"^{{pr}}_{{rs}}, " ", g^{{"s"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{ps}}_{{rt}}, " ", "VTt"^{{pr}}_{{st}}, " ", g^{{"t"}}_{{"q"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"t"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"s"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{qr}}, " ", "VTt"^{{qr}}_{{st}}, " ", g^{{"t"}}_{{"s"}}, " ", g^{{"s"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{qr}}, " ", "VTt"^{{qr}}_{{st}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qr}}_{{qr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{qr}}_{{qs}}, " ", "VTt"^{{qs}}_{{qr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{ps}}_{{ps}}, " ", "VTt"^{{pr}}_{{pr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{ps}}_{{st}}, " ", "VTt"^{{pr}}_{{rt}}, " ", g^{{"t"}}_{{"q"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"t"}}_{{"p"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{qr}}_{{qr}}, " ", "VTt"^{{qs}}_{{qs}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"q"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{rs}}_{{qr}}, " ", "VTt"^{{rs}}_{{qr}}, " ", g^{{"r"}}_{{"r"}}, " ", g^{{"p"}}_{{"s"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"q"}}_{{"s"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{rs}}, " ", "VTt"^{{rs}}_{{pr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{rs}}_{{qr}}, " ", "VTt"^{{qr}}_{{rs}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pt}}_{{rs}}, " ", "VTt"^{{rs}}_{{pt}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rs}}_{{qt}}, " ", "VTt"^{{qt}}_{{rs}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{rs}}_{{qt}}, " ", "VTt"^{{qt}}_{{rs}}, " ", g^{{"s"}}_{{"r"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"t"}}_{{"q"}}, " ", eta^{{"p"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rs}}_{{qt}}, " ", "VTt"^{{qt}}_{{rs}}, " ", g^{{"r"}}_{{"r"}}, " ", g^{{"s"}}_{{"s"}}, " ", eta^{{"t"}}_{{"q"}}, " ", eta^{{"p"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{qs}}, " ", "VTt"^{{st}}_{{qr}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"p"}}_{{"t"}}, " ", eta^{{"r"}}_{{"r"}}, " ", eta^{{"q"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{qr}}, " ", "VTt"^{{st}}_{{qs}}, " ", g^{{"r"}}_{{"s"}}, " ", g^{{"p"}}_{{"t"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"q"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{ps}}_{{st}}, " ", "VTt"^{{rt}}_{{pr}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"t"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{qr}}, " ", "VTt"^{{qs}}_{{st}}, " ", g^{{"p"}}_{{"q"}}, " ", g^{{"r"}}_{{"s"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"t"}}_{{"t"}}}]}]], "Output"
```
**Input [783]**
```
[""], "Input"
```


### [830] Direct Loop Codes (Begining with Eta-Variant)
**Input [833]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{"tmp1z", "=", RowBox[{"{", RowBox[{"CollapseDeltas", "[", RowBox[{"ByHand", "[", RowBox[{"[", "1", "]"}], "]"}], "]"}], "}"}]}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"TCLwEtar", "=", RowBox[{ RowBox[{"SpinLabelTerms", "[", "tmp1z", "]"}], "//.", {VsjkVtjkBSpinsBSpint}}]}], "]"}], " ", AllVts={};, " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", TermCode,TCLwEtar, "]"}], "]"}], " ", RowBox[{ RowBox[{"Union", "[", "AllVts", "]"}], " "}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"tmp1z", "=", RowBox[{"{", RowBox[{"CollapseDeltas", "[", RowBox[{"ByHand", "[", RowBox[{"[", "2", "]"}], "]"}], "]"}], "}"}]}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"TCLwEtar", "=", RowBox[{ RowBox[{"SpinLabelTerms", "[", "tmp1z", "]"}], "//.", {VsjkVtjkBSpinsBSpint}}]}], "]"}], " ", AllVts={};, " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", TermCode,TCLwEtar, "]"}], "]"}], " ", RowBox[{ RowBox[{"Union", "[", "AllVts", "]"}], " "}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"tmp1z", "=", RowBox[{"{", RowBox[{"CollapseDeltas", "[", RowBox[{"ByHand", "[", RowBox[{"[", "3", "]"}], "]"}], "]"}], "}"}]}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"TCLwEtar", "=", RowBox[{ RowBox[{"SpinLabelTerms", "[", "tmp1z", "]"}], "//.", {VsjkVtjkBSpinsBSpint}}]}], "]"}], " ", AllVts={};, " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", TermCode,TCLwEtar, "]"}], "]"}], " ", RowBox[{ RowBox[{"Union", "[", "AllVts", "]"}], " "}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"tmp1z", "=", RowBox[{"{", RowBox[{"CollapseDeltas", "[", RowBox[{"ByHand", "[", RowBox[{"[", "4", "]"}], "]"}], "]"}], "}"}]}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"TCLwEtar", "=", RowBox[{ RowBox[{"SpinLabelTerms", "[", "tmp1z", "]"}], "//.", {VsjkVtjkBSpinsBSpint}}]}], "]"}], " ", AllVts={};, " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", TermCode,TCLwEtar, "]"}], "]"}], " ", RowBox[{"Union", "[", "AllVts", "]"}]}], "Input"
```
**Input [975]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"ByHandSec1", "=", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SecularizeATerm,ByHand, "]"}], ",", "1"}], "]"}]}], "]"}]], "Input"
```
**Output [979]**
```
CellGroupData[{ Cell[[ RowBox[{"RecognizableForm", "[", RowBox[{"ByHandSec", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", SecularizeATerm,ByHand, "]"}]}], "]"}]}], "]"}]}], "]"}]], "Input"
```
**Input [980]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"ByHandSec", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", SecularizeATerm,ByHand, "]"}]}], "]"}]}], "]"}]}], "]"}]], "Input"
```
**Output [982]**
```
[ RowBox[{"Total", "[", "ByHand", "]"}]], "Output"
```
**Output [983]**
```
CellGroupData[{ Cell[[{ RowBox[{"RecognizableForm", "[", RowBox[{"ByHandSec", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"DistributeOperatorProduct", "[", RowBox[{"ByHand", ",", RowBox[{"STwo", "[", "\<r\>","\<t\>","\<v\>","\<x\>", "]"}]}], "]"}]}], "]"}]}], "]"}]}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", DiagonalGamma1,ByHandSec, "]"}]}], "]"}]}], "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"TCLwEtar", "=", RowBox[{ RowBox[{"SpinLabelTerms", "[", "ByHandSec", "]"}], "//.", {VsjkVtjkBSpinsBSpint}}]}], "]"}], " ", AllVts={};, " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", TermCodePUp,TCLwEtar, "]"}], "]"}], " ", RowBox[{"Union", "[", "AllVts", "]"}]}], "Input"
```
**Input [984]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{"ByHandSec", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"DistributeOperatorProduct", "[", RowBox[{"ByHand", ",", RowBox[{"STwo", "[", "\<r\>","\<t\>","\<v\>","\<x\>", "]"}]}], "]"}]}], "]"}]}], "]"}]}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", DiagonalGamma1,ByHandSec, "]"}]}], "]"}]}], "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"TCLwEtar", "=", RowBox[{ RowBox[{"SpinLabelTerms", "[", "ByHandSec", "]"}], "//.", {VsjkVtjkBSpinsBSpint}}]}], "]"}], " ", AllVts={};, " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", TermCodePUp,TCLwEtar, "]"}], "]"}], " ", RowBox[{"Union", "[", "AllVts", "]"}]}], "Input"
```
**Output [1022]**
```
[ RowBox[{"Total", "[", "1", "]"}]], "Output"
```
**Output [1044]**
```
["0"], "Output"
```
**Output [1057]**
```
[ RowBox[{"{", RowBox[{ RowBox[{ RowBox[{"\<\"Sign\"\>", " ", TemplateBox[{"\"Sign\"",RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {Sign2 "r"}, ",", RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{"ToUpperCase", "[", RowBox[{ RowBox[{"{", {},{Signr}, "}"}], "[[", 1,2, "]]"}], "]"}]}], "}"}]}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{"ToUpperCase", "[", RowBox[{ RowBox[{"{", {},{Signr}, "}"}], "[[", 1,2, "]]"}], "]"}]}], "}"}]}, "Subsuperscript"]}], "+", RowBox[{"2", " ", "\<\"Sign\"\>", " ", TemplateBox[{"\"Sign\"",RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", "Sign",{}, "}"}], ",", RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{"ToUpperCase", "[", RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"Complement", "[", "]"}], ",", RowBox[{"Complement", "[", "\"r\"", "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{ RowBox[{"Complement", "[", "s","r", "]"}], "[[", 1,2, "]]"}]}], "}"}], "}"}]}], "]"}], "[[", 1,2, "]]"}], "]"}]}], "}"}]}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{"ToUpperCase", "[", RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"Complement", "[", "]"}], ",", RowBox[{"Complement", "[", "\"r\"", "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{ RowBox[{"Complement", "[", "s","r", "]"}], "[[", 1,2, "]]"}]}], "}"}], "}"}]}], "]"}], "[[", 1,2, "]]"}], "]"}]}], "}"}]}, "Subsuperscript"]}]}], ",", RowBox[{ RowBox[{ RowBox[{ RowBox[{"{", {},{\<"Sign"\>\<"r"\>}, "}"}], "[[", 1,2, "]]"}], " ", TemplateBox[{"\"Sign\"",RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {Sign2 "r"}, ",", RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{"ToUpperCase", "[", RowBox[{ RowBox[{"{", {},{Signr}, "}"}], "[[", 1,2, "]]"}], "]"}]}], "}"}]}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{"ToUpperCase", "[", RowBox[{ RowBox[{"{", {},{Signr}, "}"}], "[[", 1,2, "]]"}], "]"}]}], "}"}]}, "Subsuperscript"]}], "+", RowBox[{"2", " ", RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"Complement", "[", "]"}], ",", RowBox[{"Complement", "[", "\<\"r\"\>", "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{ RowBox[{"Complement", "[", \<"s"\>,\<"r"\>, "]"}], "[[", 1,2, "]]"}]}], "}"}], "}"}]}], "]"}], "[[", 1,2, "]]"}], " ", TemplateBox[{"\"Sign\"",RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", "Sign",{}, "}"}], ",", RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{"ToUpperCase", "[", RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"Complement", "[", "]"}], ",", RowBox[{"Complement", "[", "\"r\"", "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{ RowBox[{"Complement", "[", "s","r", "]"}], "[[", 1,2, "]]"}]}], "}"}], "}"}]}], "]"}], "[[", 1,2, "]]"}], "]"}]}], "}"}]}], "}"}], ",", "3"}], "}"}],RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{"ToUpperCase", "[", RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"Complement", "[", "]"}], ",", RowBox[{"Complement", "[", "\"r\"", "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\"Sign\"", ",", RowBox[{ RowBox[{"Complement", "[", "s","r", "]"}], "[[", 1,2, "]]"}]}], "}"}], "}"}]}], "]"}], "[[", 1,2, "]]"}], "]"}]}], "}"}]}, "Subsuperscript"]}]}]}], "}"}]], "Output"
```
**Output [1059]**
```
[ RowBox[{"StringJoin", "[", RowBox[{ RowBox[{"TermCodePUp", "[", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<"Sign"\>2 \<"r"\>}, ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{ RowBox[{"{", {},{\<"Sign"\>\<"r"\>}, "}"}], "[[", 1,2, "]]"}]}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<"Sign"\>2 \<"r"\>}, ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"ToUpperCase", "[", RowBox[{ RowBox[{"{", {},{\<"Sign"\>\<"r"\>}, "}"}], "[[", 1,2, "]]"}], "]"}]}], "}"}]}], "}"}]}], "}"}], "]"}], "[[", RowBox[{ RowBox[{"TermCodePUp", "[", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", \<"Sign"\>,{}, "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"Complement", "[", "]"}], ",", RowBox[{"Complement", "[", "\<\"r\"\>", "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{ RowBox[{"Complement", "[", \<"s"\>,\<"r"\>, "]"}], "[[", 1,2, "]]"}]}], "}"}], "}"}]}], "]"}], "[[", 1,2, "]]"}]}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", \<"Sign"\>,{}, "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"ToUpperCase", "[", RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"Complement", "[", "]"}], ",", RowBox[{"Complement", "[", "\<\"r\"\>", "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{ RowBox[{"Complement", "[", \<"s"\>,\<"r"\>, "]"}], "[[", 1,2, "]]"}]}], "}"}], "}"}]}], "]"}], "[[", 1,2, "]]"}], "]"}]}], "}"}]}], "}"}]}], "}"}], "]"}], ",", RowBox[{"TermCodePUp", "[", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"{", \<"Sign"\>,{}, "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"Complement", "[", "]"}], ",", RowBox[{"Complement", "[", "\<\"r\"\>", "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{ RowBox[{"Complement", "[", \<"s"\>,\<"r"\>, "]"}], "[[", 1,2, "]]"}]}], "}"}], "}"}]}], "]"}], "[[", 1,2, "]]"}]}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", \<"Sign"\>,{}, "}"}], ",", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{"ToUpperCase", "[", RowBox[{ RowBox[{"Join", "[", RowBox[{ RowBox[{"Complement", "[", "]"}], ",", RowBox[{"Complement", "[", "\<\"r\"\>", "]"}], ",", RowBox[{"{", RowBox[{"{", RowBox[{"\<\"Sign\"\>", ",", RowBox[{ RowBox[{"Complement", "[", \<"s"\>,\<"r"\>, "]"}], "[[", 1,2, "]]"}]}], "}"}], "}"}]}], "]"}], "[[", 1,2, "]]"}], "]"}]}], "}"}]}], "}"}]}], "}"}], "]"}]}], "]]"}], "]"}]], "Output"
```
**Output [1060]**
```
[ {}], "Output"
```
**Output [1061]**
```
CellGroupData[{ Cell[[ RowBox[{"Map", "[", Length,TCLwEtar, "]"}]], "Input"
```
**Input [1062]**
```
[ RowBox[{"Map", "[", Length,TCLwEtar, "]"}]], "Input"
```
**Output [1063]**
```
[ RowBox[{"2", "[[", 2,2, "]]"}]], "Output"
```
**Input [1064]**
```
[ (*Now I have to loop over all four secular possibilities.*)], "Input"
```
**Input [1067]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{"dFinalTCL2", "=", RowBox[{"AssociativeSimplify", "[", RowBox[{"Map", "[", RowBox[{"RelabelTerm", ",", RowBox[{"Map", "[", CollapseDeltas,FinalTCL2//.EtaToDelta, "]"}]}], "]"}], "]"}]}], "]"}], " ", RowBox[{ RowBox[{"RecognizableForm", "[", RowBox[{"TCLwEtar", "=", RowBox[{ RowBox[{"SpinLabelTerms", "[", "dFinalTCL2", "]"}], "//.", {VsjkVtjkBSpinsBSpint}}]}], "]"}], ";"}], " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", DirectLoopCodeSec,TCLwEtar, "]"}], "]"}], " "}], "Input"
```


#### [1071] Direct Loop Code


### [1077] Secular Only CC Direct Code
**Input [1079]**
```
[ RowBox[{"RecognizableForm", "[", "sTCLwEta", "]"}]], "Input"
```
**Output [1081]**
```
[ RowBox[{"Total", "[", "sTCLwEta", "]"}]], "Output"
```
**Input [1083]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"Map", "[", RelabelTerm,sTCLwEta, "]"}], "]"}]], "Input"
```
**Input [1087]**
```
[ RowBox[{"Map", "[", RowBox[{ RowBox[{ RowBox[{"Length", "[", RowBox[{"Union", "[", RowBox[{"Flatten", "[", RowBox[{"TermIndices", "[", "#", "]"}], "]"}], "]"}], "]"}], "&"}], ",", "sTCLwEta"}], "]"}]], "Input"
```
> So the cost has gone down, dramatically. Direct loop code?
**Input [1090]**
```
[ Secularization by decoupling populations and coherences.], "Input"
```
**Input [1139]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{"SLT2", "=", RowBox[{"Select", "[", RowBox[{"SLTnopp", ",", RowBox[{ RowBox[{ RowBox[{"NumIndicesInTerm", "[", "#", "]"}], "==", "2"}], "&"}]}], "]"}]}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"SLT3", "=", RowBox[{"Select", "[", RowBox[{"SLTnopp", ",", RowBox[{ RowBox[{ RowBox[{"NumIndicesInTerm", "[", "#", "]"}], "==", "3"}], "&"}]}], "]"}]}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"SLT4", "=", RowBox[{"Select", "[", RowBox[{"SLTnopp", ",", RowBox[{ RowBox[{ RowBox[{"NumIndicesInTerm", "[", "#", "]"}], "==", "4"}], "&"}]}], "]"}]}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"SLT5", "=", RowBox[{"Select", "[", RowBox[{"SLTnopp", ",", RowBox[{ RowBox[{ RowBox[{"NumIndicesInTerm", "[", "#", "]"}], "==", "5"}], "&"}]}], "]"}]}], "]"}]}], "Input"
```
**Output [1148]**
```
CellGroupData[{ Cell[[ RowBox[{"RecognizableForm", "[", SLT5//.{RIsRItBSpinsBSpint}, "]"}]], "Input"
```
**Input [1149]**
```
[ RowBox[{"RecognizableForm", "[", SLT5//.{RIsRItBSpinsBSpint}, "]"}]], "Input"
```
**Output [1150]**
```
["0"], "Output"
```
**Input [1151]**
```
[ (*Figure out what three indices are needed on VTeliminate them to make a 3-index intermediate (always possible in fourth order time, "*)"], "Input"
```
**Input [1156]**
```
[{ RowBox[{"StringJoin", "[", RowBox[{"Map", "[", DirectLoopCodeSep,SLTpp2, "]"}], "]"}], " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", DirectLoopCodeSep,SLTpp3, "]"}], "]"}], " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", DirectLoopCodeSep,SLT2, "]"}], "]"}], " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", DirectLoopCodeSep,SLT3, "]"}], "]"}], " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", DirectLoopCodeSep,SLT4, "]"}], "]"}], " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", DirectLoopCodeSep,SLT5, "]"}], "]"}]}], "Input"
```
**Input [1178]**
```
[{ RowBox[{"StringJoin", "[", RowBox[{"Map", "[", DirectLoopCodeSep2,SLTpp2, "]"}], "]"}], " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", DirectLoopCodeSep2,SLTpp3, "]"}], "]"}], " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", DirectLoopCodeSep2,SLT2, "]"}], "]"}], " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", DirectLoopCodeSep2,SLT3, "]"}], "]"}], " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", DirectLoopCodeSep2,SLT4, "]"}], "]"}], " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", DirectLoopCodeSep2,SLT5, "]"}], "]"}]}], "Input"
```
> <For these terms really I need to eliminate indices in Vs .. It can begin as a cube but eliminate 1 index and it can become a matrix. It' s already three indices in all but a few terms. So: new plan. Use a loop to manually generate a 2-index intermediate for Vs(p, ...) * whatever. and then use the normal plan-multiplications code for the rest. With RI that[CloseCurlyQuote]s nominally O(n^5th) >


#### [1201] MatrixCode
**Output [1202]**
```
CellGroupData[{ Cell[[{ RowBox[{ RowBox[{"CCsl", "=", RowBox[{"AssociativeSimplify", "[", RowBox[{ RowBox[{"(", RowBox[{ RowBox[{"SpinLabelTerms", "[", "OnlyCC", "]"}], "/.", "AddVOscS"}], ")"}], "/.", "AddVOscT"}], "]"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", CCari=CCsl//.{RIsRItBSpin}, "]"}]}], "Input"
```
**Input [1203]**
```
[{ RowBox[{ RowBox[{"CCsl", "=", RowBox[{"AssociativeSimplify", "[", RowBox[{ RowBox[{"(", RowBox[{ RowBox[{"SpinLabelTerms", "[", "OnlyCC", "]"}], "/.", "AddVOscS"}], ")"}], "/.", "AddVOscT"}], "]"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", CCari=CCsl//.{RIsRItBSpin}, "]"}]}], "Input"
```
**Output [1207]**
```
[ RowBox[{ RowBox[{"\<\"Error\"\>", "[[", 1,2, "]]"}], "+", TemplateBox[{RowBox[{ {OnlyCC12}, "[[", 1,1, "]]"}],"2","1"}, "Subsuperscript"]}]], "Output"
```
**Input [1209]**
```
[ RowBox[{"StringJoin", "[", RowBox[{"Map", "[", PlanMultiplications,CCari, "]"}], "]"}]], "Input"
```
**Output [1218]**
```
[ RowBox[{"\<\"{\\n fac=Error[[1,2]];\\nBl1=BpqR.slice(l)%\"\>", "<>", "Null", "<>", "\<\";\\nBl2=BpqR.slice(l)%\"\>", "<>", "Null", "<>", "\<\";\\nBm1=BpqR.slice(m)%\"\>", "<>", "Null", "<>", "\<\";\\nBm2=BpqR.slice(m)%\"\>", "<>", "Null", "<>", "\<\";\\nfac *= accu(\"\>", "<>", "Null", "<>", "\<\");\\nfac *= accu(\"\>", "<>", "Null", "<>", "\<\");\\nfac *= accu(\"\>", "<>", "Null", "<>", "\<\");\\n{\\n fac=1.;\\nBl1=BpqR.slice(l)%\"\>", "<>", RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{ RowBox[{ -2, "+", RowBox[{ {}, "[[", 1,1, "]]"}]}], "==", "1"}], "&&", RowBox[{ RowBox[{ -2, "+", RowBox[{ {}, "[[", 1,1, "]]"}]}], "==", "1"}]}], ",", RowBox[{"Return", "[", "\<\"exp(-1.0*j*timeS*dpdp)\"\>", "]"}], ",", RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{"NegInd1", "==", "0"}], "&&", RowBox[{"NegInd2", "==", "1"}]}], ",", RowBox[{"Return", "[", "\<\"exp(j*timeS*dpdm)\"\>", "]"}], ",", RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{"NegInd1", "==", "1"}], "&&", RowBox[{"NegInd2", "==", "0"}]}], ",", RowBox[{"Return", "[", "\<\"exp(-1.0*j*timeS*dpdm)\"\>", "]"}], ",", RowBox[{"If", "[", RowBox[{ RowBox[{ RowBox[{"NegInd1", "==", "0"}], "&&", RowBox[{"NegInd2", "==", "0"}]}], ",", RowBox[{"Return", "[", "\<\"exp(j*timeS*dpdp)\"\>", "]"}]}], "]"}]}], "]"}]}], "]"}]}], "]"}], "<>", "\<\";\\nBl2=BpqR.slice(l)%\"\>", "<>", "Null", "<>", "\<\";\\nBm1=BpqR.slice(m)%\"\>", "<>", "Null", "<>", "\<\";\\nBm2=BpqR.slice(m)%\"\>", "<>", "Null", "<>", "\<\";\\nfac *= accu(\"\>", "<>", "Null", "<>", "\<\");\\nfac *= accu(\"\>", "<>", "Null", "<>", "\<\");\\nfac *= accu(\"\>", "<>", RowBox[{ {\<"Bl1"\>\<"Bl2"\>\<"Bm1"\>\<"Bm2"\>\<"Rho"\>\<"Rho"\>\<"I"\>}, "[[", RowBox[{"{", RowBox[{"{", RowBox[{ RowBox[{ {}, "[[", 1,1, "]]"}], ",", "1"}], "}"}], "}"}], "]]"}], "<>", "\<\"%Bm2);\\n\"\>"}]], "Output"
```


### [1220] The Oscillating Time Factors
**Input [1226]**
```
[ RowBox[{" ", (*G<*), " ", RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", RowBox[{"-", "oo"}], "t"], RowBox[{ RowBox[{"(", RowBox[{ RowBox[{"I", "*", RowBox[{"HeavisideTheta", "[", RowBox[{"tau", "-", "t"}], "]"}], "*", RowBox[{"Exp", "[", RowBox[{"I", "*", (L), "*", RowBox[{"(", "tau", ")"}]}], "]"}]}], "+", RowBox[{"I", "*", RowBox[{"HeavisideTheta", "[", RowBox[{"t", "-", "tau"}], "]"}], "*", RowBox[{"Exp", "[", I*(-L)*(t), "]"}]}]}], ")"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "FullSimplify"}]}]], "Input"
```
**Input [1227]**
```
[" "], "Input"
```
**Output [1228]**
```
CellGroupData[{ Cell[[ RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", RowBox[{"-", "oo"}], "0"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{"L", "+", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], "*", RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{ -L, "-", RowBox[{"I", "*", "eps"}]}], ")"}], "*", (t)}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "FullSimplify"}]], "Input"
```
**Input [1229]**
```
[ RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", RowBox[{"-", "oo"}], "0"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{"L", "+", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], "*", RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{ -L, "-", RowBox[{"I", "*", "eps"}]}], ")"}], "*", (t)}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "FullSimplify"}]], "Input"
```
**Output [1230]**
```
[ RowBox[{"ConditionalExpression", "[", RowBox[{ FractionBox[ SuperscriptBox["E", RowBox[{"t", " ", RowBox[{"(", RowBox[{ RowBox[{ RowBox[{"-", "I"}], " ", "L"}], "+", "eps"}], ")"}]}]], RowBox[{ RowBox[{"I", " ", "L"}], "-", "eps"}]], ",", RowBox[{"eps", "<", "0"}]}], "]"}]], "Output"
```
**Output [1231]**
```
CellGroupData[{ Cell[[ RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", RowBox[{"-", "oo"}], "t"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{"L", "+", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], "*", RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{ -L, "-", RowBox[{"I", "*", "eps"}]}], ")"}], "*", (t)}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "FullSimplify"}]], "Input"
```
**Input [1232]**
```
[ RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", RowBox[{"-", "oo"}], "t"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{"L", "+", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], "*", RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{ -L, "-", RowBox[{"I", "*", "eps"}]}], ")"}], "*", (t)}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "FullSimplify"}]], "Input"
```
**Output [1233]**
```
[ RowBox[{"ConditionalExpression", "[", RowBox[{ RowBox[{"-", FractionBox["1", RowBox[{ RowBox[{ RowBox[{"-", "I"}], " ", "L"}], "+", "eps"}]]}], ",", RowBox[{"eps", "<", "0"}]}], "]"}]], "Output"
```
**Output [1234]**
```
CellGroupData[{ Cell[[ RowBox[{ RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{"L", "+", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "Simplify"}], "//", "TrigToExp"}]], "Input"
```
**Input [1235]**
```
[ RowBox[{ RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{"L", "+", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "Simplify"}], "//", "TrigToExp"}]], "Input"
```
**Output [1236]**
```
[ RowBox[{"ConditionalExpression", "[", RowBox[{ FractionBox["1", RowBox[{ RowBox[{ RowBox[{"-", "I"}], " ", "L"}], "+", "eps"}]], ",", RowBox[{"eps", ">", "0"}]}], "]"}]], "Output"
```
**Output [1240]**
```
CellGroupData[{ Cell[[ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0, "&&", RowBox[{"eps", ">", "0"}]}], ",", RowBox[{"Simplify", "[", RowBox[{"%", "+", RowBox[{"Conjugate", "[", "%", "]"}]}], "]"}]}], "]"}]], "Input"
```
**Input [1241]**
```
[ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0, "&&", RowBox[{"eps", ">", "0"}]}], ",", RowBox[{"Simplify", "[", RowBox[{"%", "+", RowBox[{"Conjugate", "[", "%", "]"}]}], "]"}]}], "]"}]], "Input"
```
**Output [1242]**
```
[ FractionBox[ RowBox[{ RowBox[{ SuperscriptBox["E", RowBox[{"t", " ", RowBox[{"(", RowBox[{ RowBox[{ RowBox[{"-", "I"}], " ", "L"}], "+", "eps"}], ")"}]}]], " ", RowBox[{"(", RowBox[{"L", "-", RowBox[{"I", " ", "eps"}]}], ")"}]}], "-", RowBox[{ SuperscriptBox["E", RowBox[{"t", " ", RowBox[{"(", RowBox[{ RowBox[{"I", " ", "L"}], "+", "eps"}], ")"}]}]], " ", RowBox[{"(", RowBox[{"L", "+", RowBox[{"I", " ", "eps"}]}], ")"}]}], "+", RowBox[{"2", " ", "I", " ", "eps"}]}], RowBox[{ RowBox[{"(", RowBox[{"L", "-", RowBox[{"I", " ", "eps"}]}], ")"}], " ", RowBox[{"(", RowBox[{ RowBox[{ RowBox[{"-", "I"}], " ", "L"}], "+", "eps"}], ")"}]}]]], "Output"
```
**Output [1246]**
```
CellGroupData[{ Cell[[{ RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", "0", "t"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", (L), "*", RowBox[{"(", "tau", ")"}]}], "]"}], "*", RowBox[{"Exp", "[", I*(-L)*(t), "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "FullSimplify"}], " ", RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", RowBox[{"-", "oo"}], "0"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", (L), "*", RowBox[{"(", "tau", ")"}]}], "]"}], "*", RowBox[{"Exp", "[", I*(-L)*(t), "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "FullSimplify"}]}], "Input"
```
**Input [1247]**
```
[{ RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", "0", "t"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", (L), "*", RowBox[{"(", "tau", ")"}]}], "]"}], "*", RowBox[{"Exp", "[", I*(-L)*(t), "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "FullSimplify"}], " ", RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", RowBox[{"-", "oo"}], "0"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", (L), "*", RowBox[{"(", "tau", ")"}]}], "]"}], "*", RowBox[{"Exp", "[", I*(-L)*(t), "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "FullSimplify"}]}], "Input"
```
**Output [1248]**
```
[ FractionBox[ RowBox[{ RowBox[{"I", " ", RowBox[{"(", RowBox[{ -1, "+", RowBox[{"Cos", "[", L t, "]"}]}], ")"}]}], "+", RowBox[{"Sin", "[", L t, "]"}]}], "L"]], "Output"
```
**Output [1250]**
```
[ RowBox[{ SubsuperscriptBox["Int", RowBox[{"-", "oo"}], "0"], RowBox[{ SuperscriptBox["E", RowBox[{ RowBox[{"-", "I"}], " ", "L", " ", RowBox[{"(", RowBox[{"t", "-", "tau"}], ")"}]}]], RowBox[{"d", "tau"}]}]}]], "Output"
```
**Output [1251]**
```
CellGroupData[{ Cell[[ RowBox[{"FourierTransform", "[", RowBox[{ FractionBox[ RowBox[{"I", " ", RowBox[{"(", RowBox[{ -1, "+", SuperscriptBox["E", RowBox[{ RowBox[{"-", "I"}], " ", "L", " ", "t"}]]}], ")"}]}], "L"], ",", "t", ",", "w"}], "]"}]], "Input"
```
**Input [1252]**
```
[ RowBox[{"FourierTransform", "[", RowBox[{ FractionBox[ RowBox[{"I", " ", RowBox[{"(", RowBox[{ -1, "+", SuperscriptBox["E", RowBox[{ RowBox[{"-", "I"}], " ", "L", " ", "t"}]]}], ")"}]}], "L"], ",", "t", ",", "w"}], "]"}]], "Input"
```
**Output [1253]**
```
[ RowBox[{ RowBox[{"-", FractionBox[ RowBox[{"I", " ", SqrtBox[ RowBox[{"2", " ", "Pi"}]], " ", RowBox[{"DiracDelta", "[", "w", "]"}]}], "L"]}], "+", FractionBox[ RowBox[{"I", " ", SqrtBox[ RowBox[{"2", " ", "Pi"}]], " ", RowBox[{"DiracDelta", "[", RowBox[{ -L, "+", "w"}], "]"}]}], "L"]}]], "Output"
```
**Output [1254]**
```
CellGroupData[{ Cell[[ RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", "0", "t"], RowBox[{ RowBox[{"Exp", "[", RowBox[{ -I, "*", RowBox[{"(", RowBox[{"L", "-", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "FullSimplify"}]], "Input"
```
**Input [1255]**
```
[ RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", "0", "t"], RowBox[{ RowBox[{"Exp", "[", RowBox[{ -I, "*", RowBox[{"(", RowBox[{"L", "-", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "FullSimplify"}]], "Input"
```
**Output [1256]**
```
[ RowBox[{"-", FractionBox[ RowBox[{ -1, "+", SuperscriptBox["E", RowBox[{ RowBox[{ RowBox[{"-", "I"}], " ", "L", " ", "t"}], "-", RowBox[{"t", " ", "eps"}]}]]}], RowBox[{ RowBox[{"I", " ", "L"}], "+", "eps"}]]}]], "Output"
```
**Output [1257]**
```
CellGroupData[{ Cell[[ RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", RowBox[{ -I, "*", RowBox[{"(", RowBox[{"L", "-", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "FullSimplify"}]], "Input"
```
**Input [1258]**
```
[ RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", RowBox[{ -I, "*", RowBox[{"(", RowBox[{"L", "-", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "FullSimplify"}]], "Input"
```
**Output [1259]**
```
[ RowBox[{"ConditionalExpression", "[", RowBox[{ FractionBox["1", RowBox[{ RowBox[{"I", " ", "L"}], "+", "eps"}]], ",", RowBox[{"eps", ">", "0"}]}], "]"}]], "Output"
```
**Output [1260]**
```
CellGroupData[{ Cell[[ RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", RowBox[{ -I, "*", RowBox[{"(", RowBox[{"L", "-", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", RowBox[{"-", "tau"}], ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "FullSimplify"}]], "Input"
```
**Input [1261]**
```
[ RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", RowBox[{ -I, "*", RowBox[{"(", RowBox[{"L", "-", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", RowBox[{"-", "tau"}], ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "//", "FullSimplify"}]], "Input"
```
**Output [1262]**
```
[ RowBox[{"ConditionalExpression", "[", RowBox[{ RowBox[{"-", FractionBox["1", RowBox[{ RowBox[{"I", " ", "L"}], "+", "eps"}]]}], ",", RowBox[{"eps", "<", "0"}]}], "]"}]], "Output"
```
**Output [1263]**
```
CellGroupData[{ Cell[[ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{ RowBox[{"Simplify", "[", RowBox[{"ComplexExpand", "[", RowBox[{ SubsuperscriptBox["Int", RowBox[{"-", "oo"}], "0"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{ -L, "+", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}], "]"}], "+", RowBox[{"Simplify", "[", RowBox[{"ComplexExpand", "[", RowBox[{ SubsuperscriptBox["Int", RowBox[{"-", "oo"}], "0"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{"L", "+", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}], "]"}]}]}], "]"}]], "Input"
```
**Input [1264]**
```
[ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{ RowBox[{"Simplify", "[", RowBox[{"ComplexExpand", "[", RowBox[{ SubsuperscriptBox["Int", RowBox[{"-", "oo"}], "0"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{ -L, "+", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}], "]"}], "+", RowBox[{"Simplify", "[", RowBox[{"ComplexExpand", "[", RowBox[{ SubsuperscriptBox["Int", RowBox[{"-", "oo"}], "0"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{"L", "+", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}], "]"}]}]}], "]"}]], "Input"
```
**Output [1265]**
```
[ RowBox[{"ConditionalExpression", "[", RowBox[{ RowBox[{ FractionBox["I", RowBox[{"L", "-", RowBox[{"I", " ", "eps"}]}]], "-", FractionBox["I", RowBox[{"L", "+", RowBox[{"I", " ", "eps"}]}]]}], ",", RowBox[{"eps", "<", "0"}]}], "]"}]], "Output"
```
`[Output 1266: graphics omitted]`
**Input [1267]**
```
[ RowBox[{"Plot", "[", RowBox[{ RowBox[{ RowBox[{"{", RowBox[{ RowBox[{"Re", "[", RowBox[{ FractionBox["I", RowBox[{"L", "-", RowBox[{"I", " ", "eps"}]}]], "-", FractionBox["I", RowBox[{"L", "+", RowBox[{"I", " ", "eps"}]}]]}], "]"}], ",", RowBox[{"Im", "[", RowBox[{ FractionBox["I", RowBox[{"L", "-", RowBox[{"I", " ", "eps"}]}]], "-", FractionBox["I", RowBox[{"L", "+", RowBox[{"I", " ", "eps"}]}]]}], "]"}]}], "}"}], "/.", RowBox[{"eps", "->", "0.01"}]}], ",", {L-.1.1}, ",", RowBox[{"PlotRange", "->", "All"}]}], "]"}]], "Input"
```
`[Output 1268: graphics omitted]`
> <In the case where the sign of epsilon is opposite you get the sum of a delta and a normal denominator >
**Output [1276]**
```
CellGroupData[{ Cell[[ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{"ComplexExpand", "[", RowBox[{ RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{"L", "+", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}], "+", RowBox[{"Conjugate", "[", RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", RowBox[{ -I, "*", RowBox[{"(", RowBox[{"L", "-", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "]"}]}], "]"}]], "Input"
```
**Input [1277]**
```
[ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ RowBox[{"L", "in", "Reals"}], "&&", RowBox[{"eps", "in", "Reals"}], "&&", L>0, "&&", RowBox[{"t", "in", "Reals"}], "&&", t>0, "&&", L>0}], ",", RowBox[{"Simplify", "[", RowBox[{"ComplexExpand", "[", RowBox[{ RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{"L", "+", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}], "+", RowBox[{"Conjugate", "[", RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", RowBox[{ -I, "*", RowBox[{"(", RowBox[{"L", "-", RowBox[{"I", "*", "eps"}]}], ")"}], "*", RowBox[{"(", "tau", ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}], "]"}]}], "]"}], "]"}]}], "]"}]], "Input"
```
**Output [1278]**
```
[ RowBox[{"ConditionalExpression", "[", RowBox[{ FractionBox[ RowBox[{"2", " ", "I"}], RowBox[{"L", "+", RowBox[{"I", " ", "eps"}]}]], ",", RowBox[{"eps", ">", "0"}]}], "]"}]], "Output"
```
**Output [1279]**
```
CellGroupData[{ Cell[[ RowBox[{"ComplexExpand", "[", FractionBox[ RowBox[{"2", " ", "I", "*", RowBox[{"(", RowBox[{"L", "-", RowBox[{"I", " ", "eps"}]}], ")"}]}], RowBox[{ RowBox[{"(", RowBox[{"L", "+", RowBox[{"I", " ", "eps"}]}], ")"}], "*", RowBox[{"(", RowBox[{"L", "-", RowBox[{"I", " ", "eps"}]}], ")"}]}]], "]"}]], "Input"
```
**Input [1280]**
```
[ RowBox[{"ComplexExpand", "[", FractionBox[ RowBox[{"2", " ", "I", "*", RowBox[{"(", RowBox[{"L", "-", RowBox[{"I", " ", "eps"}]}], ")"}]}], RowBox[{ RowBox[{"(", RowBox[{"L", "+", RowBox[{"I", " ", "eps"}]}], ")"}], "*", RowBox[{"(", RowBox[{"L", "-", RowBox[{"I", " ", "eps"}]}], ")"}]}]], "]"}]], "Input"
```
**Output [1281]**
```
[ RowBox[{ FractionBox[ RowBox[{"2", " ", "I", " ", "L"}], RowBox[{ L^{2}, "+", SuperscriptBox["eps", "2"]}]], "+", FractionBox[ RowBox[{"2", " ", "eps"}], RowBox[{ L^{2}, "+", SuperscriptBox["eps", "2"]}]]}]], "Output"
```
**Output [1282]**
```
CellGroupData[{ Cell[[{ RowBox[{"Limit", "[", RowBox[{ RowBox[{ FractionBox[ RowBox[{"2", " ", "I", " ", "L"}], RowBox[{ L^{2}, "+", SuperscriptBox["eps", "2"]}]], "+", FractionBox[ RowBox[{"2", " ", "eps"}], RowBox[{ L^{2}, "+", SuperscriptBox["eps", "2"]}]]}], ",", RowBox[{"eps", "->", "0"}]}], "]"}], " ", RowBox[{"Limit", "[", RowBox[{ FractionBox[ RowBox[{"2", " ", "I", " ", "L"}], RowBox[{ L^{2}, "+", SuperscriptBox["eps", "2"]}]], ",", RowBox[{"eps", "->", "0"}]}], "]"}], " ", RowBox[{"Limit", "[", RowBox[{ FractionBox[ RowBox[{"2", " ", "eps"}], RowBox[{ L^{2}, "+", SuperscriptBox["eps", "2"]}]], ",", RowBox[{"eps", "->", "0"}]}], "]"}]}], "Input"
```
**Input [1283]**
```
[{ RowBox[{"Limit", "[", RowBox[{ RowBox[{ FractionBox[ RowBox[{"2", " ", "I", " ", "L"}], RowBox[{ L^{2}, "+", SuperscriptBox["eps", "2"]}]], "+", FractionBox[ RowBox[{"2", " ", "eps"}], RowBox[{ L^{2}, "+", SuperscriptBox["eps", "2"]}]]}], ",", RowBox[{"eps", "->", "0"}]}], "]"}], " ", RowBox[{"Limit", "[", RowBox[{ FractionBox[ RowBox[{"2", " ", "I", " ", "L"}], RowBox[{ L^{2}, "+", SuperscriptBox["eps", "2"]}]], ",", RowBox[{"eps", "->", "0"}]}], "]"}], " ", RowBox[{"Limit", "[", RowBox[{ FractionBox[ RowBox[{"2", " ", "eps"}], RowBox[{ L^{2}, "+", SuperscriptBox["eps", "2"]}]], ",", RowBox[{"eps", "->", "0"}]}], "]"}]}], "Input"
```
**Output [1284]**
```
[ FractionBox[ RowBox[{"2", " ", "I"}], "L"]], "Output"
```
**Output [1285]**
```
[ FractionBox[ RowBox[{"2", " ", "I"}], "L"]], "Output"
```
**Output [1286]**
```
["0"], "Output"
```
**Output [1287]**
```
CellGroupData[{ Cell[[ RowBox[{"Int", RowBox[{ RowBox[{"Exp", "[", RowBox[{"L", "*", RowBox[{"(", RowBox[{"t", "-", "tau"}], ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}]], "Input"
```
**Input [1288]**
```
[ RowBox[{"Int", RowBox[{ RowBox[{"Exp", "[", RowBox[{"L", "*", RowBox[{"(", RowBox[{"t", "-", "tau"}], ")"}]}], "]"}], RowBox[{"d", "tau"}]}]}]], "Input"
```
**Output [1289]**
```
[ RowBox[{"-", FractionBox[ SuperscriptBox["E", RowBox[{ L t, "-", RowBox[{"L", " ", "tau"}]}]], "L"]}]], "Output"
```
**Output [1290]**
```
CellGroupData[{ Cell[[{ RowBox[{"Clear", "[", "p", "]"}], " ", RowBox[{"Limit", "[", RowBox[{ RowBox[{"Expand", "[", RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ p>0, "&&", q>0, "&&", RowBox[{"p", "in", "Reals"}], "&&", RowBox[{"q", "in", "Reals"}], "&&", RowBox[{"t", "in", "Reals"}]}], ",", RowBox[{ SubsuperscriptBox["Int", "0", "t"], RowBox[{ RowBox[{"Exp", "[", I*s*(p-q), "]"}], RowBox[{"d", "s"}]}]}]}], "]"}], "*", RowBox[{"Exp", "[", -I*t*(p-q), "]"}]}], "]"}], ",", RowBox[{"t", "->", "oo"}]}], "]"}]}], "Input"
```
**Input [1291]**
```
[{ RowBox[{"Clear", "[", "p", "]"}], " ", RowBox[{"Limit", "[", RowBox[{ RowBox[{"Expand", "[", RowBox[{ RowBox[{"Assuming", "[", RowBox[{ RowBox[{ p>0, "&&", q>0, "&&", RowBox[{"p", "in", "Reals"}], "&&", RowBox[{"q", "in", "Reals"}], "&&", RowBox[{"t", "in", "Reals"}]}], ",", RowBox[{ SubsuperscriptBox["Int", "0", "t"], RowBox[{ RowBox[{"Exp", "[", I*s*(p-q), "]"}], RowBox[{"d", "s"}]}]}]}], "]"}], "*", RowBox[{"Exp", "[", -I*t*(p-q), "]"}]}], "]"}], ",", RowBox[{"t", "->", "oo"}]}], "]"}]}], "Input"
```
**Output [1292]**
```
[ RowBox[{"Limit", "[", RowBox[{ RowBox[{ RowBox[{"-", FractionBox["I", p-q]}], "+", FractionBox[ RowBox[{"I", " ", SuperscriptBox["E", RowBox[{ RowBox[{"-", "I"}], " ", (p-q), " ", "t"}]]}], p-q]}], ",", RowBox[{"t", "->", "oo"}]}], "]"}]], "Output"
```
**Output [1293]**
```
CellGroupData[{ Cell[[ RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", I*(s)*(p-q), "]"}], RowBox[{"d", "s"}]}]}]], "Input"
```
**Input [1294]**
```
[ RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", I*(s)*(p-q), "]"}], RowBox[{"d", "s"}]}]}]], "Input"
```
**Output [1295]**
```
[ RowBox[{"ConditionalExpression", "[", RowBox[{ FractionBox["I", p-q], ",", RowBox[{ RowBox[{"Im", "[", "p", "]"}], ">", RowBox[{"Im", "[", "q", "]"}]}]}], "]"}]], "Output"
```
**Input [1296]**
```
CellGroupData[{ Cell[[ RowBox[{"Integrate", "[", RowBox[{ RowBox[{"Integrate", "[", RowBox[{ RowBox[{"Integrate", "[", r*z,{r01-z}, "]"}], ",", RowBox[{"{", RowBox[{"\[Theta]", ",", "0", ",", RowBox[{"2", "Pi"}]}], "}"}]}], "]"}], ",", {z01}}], "]"}]], "Input"
```
**Input [1297]**
```
[ RowBox[{"Integrate", "[", RowBox[{ RowBox[{"Integrate", "[", RowBox[{ RowBox[{"Integrate", "[", r*z,{r01-z}, "]"}], ",", RowBox[{"{", RowBox[{"\[Theta]", ",", "0", ",", RowBox[{"2", "Pi"}]}], "}"}]}], "]"}], ",", {z01}}], "]"}]], "Input"
```
**Output [1298]**
```
[ FractionBox["Pi", "12"]], "Output"
```
**Input [1299]**
```
CellGroupData[{ Cell[[ RowBox[{"Integrate", "[", RowBox[{ RowBox[{"Integrate", "[", RowBox[{ RowBox[{"Integrate", "[", r,{r01-z}, "]"}], ",", RowBox[{"{", RowBox[{"\[Theta]", ",", "0", ",", RowBox[{"2", "Pi"}]}], "}"}]}], "]"}], ",", {z01}}], "]"}]], "Input"
```
**Input [1300]**
```
[ RowBox[{"Integrate", "[", RowBox[{ RowBox[{"Integrate", "[", RowBox[{ RowBox[{"Integrate", "[", r,{r01-z}, "]"}], ",", RowBox[{"{", RowBox[{"\[Theta]", ",", "0", ",", RowBox[{"2", "Pi"}]}], "}"}]}], "]"}], ",", {z01}}], "]"}]], "Input"
```
**Output [1301]**
```
[ FractionBox["Pi", "3"]], "Output"
```
**Output [1302]**
```
CellGroupData[{ Cell[[ RowBox[{ SubsuperscriptBox["Int", "0", "t"], RowBox[{ RowBox[{"Exp", "[", I*(s-t)*(p-q), "]"}], RowBox[{"d", "s"}]}]}]], "Input"
```
**Input [1303]**
```
[ RowBox[{ SubsuperscriptBox["Int", "0", "t"], RowBox[{ RowBox[{"Exp", "[", I*(s-t)*(p-q), "]"}], RowBox[{"d", "s"}]}]}]], "Input"
```
**Output [1304]**
```
[ RowBox[{"-", FractionBox[ RowBox[{"I", " ", RowBox[{"(", RowBox[{"1", "-", SuperscriptBox["E", RowBox[{ RowBox[{"-", "I"}], " ", (p-q), " ", "t"}]]}], ")"}]}], p-q]}]], "Output"
```
**Output [1305]**
```
CellGroupData[{ Cell[[ RowBox[{ SubsuperscriptBox["Int", "t", 2t], RowBox[{ RowBox[{"Exp", "[", I*(s-t)*(p-q), "]"}], RowBox[{"d", "s"}]}]}]], "Input"
```
**Input [1306]**
```
[ RowBox[{ SubsuperscriptBox["Int", "t", 2t], RowBox[{ RowBox[{"Exp", "[", I*(s-t)*(p-q), "]"}], RowBox[{"d", "s"}]}]}]], "Input"
```
**Output [1307]**
```
[ RowBox[{"-", FractionBox[ RowBox[{"I", " ", RowBox[{"(", RowBox[{ -1, "+", SuperscriptBox["E", RowBox[{"I", " ", (p-q), " ", "t"}]]}], ")"}]}], p-q]}]], "Output"
```
**Output [1308]**
```
CellGroupData[{ Cell[[ RowBox[{ SubsuperscriptBox["Int", RowBox[{"-", "oo"}], "0"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", "tau", ")"}], "*", (p-q)}], "]"}], RowBox[{"d", "tau"}]}]}]], "Input"
```
**Input [1309]**
```
[ RowBox[{ SubsuperscriptBox["Int", RowBox[{"-", "oo"}], "0"], RowBox[{ RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", "tau", ")"}], "*", (p-q)}], "]"}], RowBox[{"d", "tau"}]}]}]], "Input"
```
**Output [1310]**
```
[ RowBox[{"ConditionalExpression", "[", RowBox[{ FractionBox["I", -p+q], ",", RowBox[{ RowBox[{"Im", "[", "p", "]"}], "<", RowBox[{"Im", "[", "q", "]"}]}]}], "]"}]], "Output"
```
**Input [1312]**
```
[ RowBox[{"Limit", "[", RowBox[{ RowBox[{ SubsuperscriptBox["Int", "0", "t"], RowBox[{ RowBox[{"Exp", "[", I*(s-t)*(p-q), "]"}], RowBox[{"d", "s"}]}]}], ",", RowBox[{"t", "->", "oo"}]}], "]"}]], "Input"
```
**Output [1314]**
```
CellGroupData[{ Cell[[ RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", -I*s*(p-q), "]"}], RowBox[{"d", "s"}]}]}]], "Input"
```
**Input [1315]**
```
[ RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", -I*s*(p-q), "]"}], RowBox[{"d", "s"}]}]}]], "Input"
```
**Output [1316]**
```
[ RowBox[{"ConditionalExpression", "[", RowBox[{ FractionBox["I", -p+q], ",", RowBox[{ RowBox[{"Im", "[", "p", "]"}], "<", RowBox[{"Im", "[", "q", "]"}]}]}], "]"}]], "Output"
```
**Input [1318]**
```
[ RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", -I*s*(q-p), "]"}], RowBox[{"d", "s"}]}]}]], "Input"
```
> Trying Out Laplace
**Output [1321]**
```
CellGroupData[{ Cell[[{ RowBox[{ RowBox[{"dtau", "=", "0.5"}], ";"}], " ", p=0.1;, " ", q=0.3;, " ", r=0.005;, " ", s=-0.100;, " ", RowBox[{"dtau", "*", RowBox[{"Total", "[", RowBox[{"Table", "[", RowBox[{ RowBox[{ RowBox[{"Exp", "[", RowBox[{ RowBox[{"-", "tau"}], "*", "p"}], "]"}], RowBox[{"Exp", "[", RowBox[{ RowBox[{"-", "tau"}], "*", "q"}], "]"}], RowBox[{"Exp", "[", RowBox[{"tau", "*", "r"}], "]"}], RowBox[{"Exp", "[", RowBox[{"tau", "*", "s"}], "]"}]}], ",", RowBox[{"{", RowBox[{"tau", ",", "0.0", ",", "20", ",", "dtau"}], "}"}]}], "]"}], "]"}]}], " ", 1/(p+q-r-s), " ", RowBox[{"Clear", "[", RowBox[{"p", ",", "q", ",", "r", ",", "s", ",", "dtau"}], "]"}]}], "Input"
```
**Output [1323]**
```
["2.2804146708936255`"], "Output"
```
**Output [1324]**
```
["2.0202020202020203`"], "Output"
```
**Output [1325]**
```
CellGroupData[{ Cell[[ RowBox[{"Solve", "[", RowBox[{ RowBox[{"x", "==", RowBox[{"Exp", "[", RowBox[{"-", "tau"}], "]"}]}], ",", "tau"}], "]"}]], "Input"
```
**Input [1326]**
```
[ RowBox[{"Solve", "[", RowBox[{ RowBox[{"x", "==", RowBox[{"Exp", "[", RowBox[{"-", "tau"}], "]"}]}], ",", "tau"}], "]"}]], "Input"
```
**Output [1327]**
```
[ RowBox[{"{", RowBox[{"{", RowBox[{"tau", "->", RowBox[{"ConditionalExpression", "[", RowBox[{ RowBox[{ RowBox[{"2", " ", "I", " ", "Pi", " ", RowBox[{"C", "[", "1", "]"}]}], "+", RowBox[{"Log", "[", (1)/(x), "]"}]}], ",", RowBox[{ RowBox[{"C", "[", "1", "]"}], "in", "Integers"}]}], "]"}]}], "}"}], "}"}]], "Output"
```
**Output [1328]**
```
CellGroupData[{ Cell[[ RowBox[{ SubscriptBox["\[PartialD]", "x"], RowBox[{"(", RowBox[{"-", RowBox[{"Log", "[", "x", "]"}]}], ")"}]}]], "Input"
```
**Input [1329]**
```
[ RowBox[{ SubscriptBox["\[PartialD]", "x"], RowBox[{"(", RowBox[{"-", RowBox[{"Log", "[", "x", "]"}]}], ")"}]}]], "Input"
```
**Output [1330]**
```
[ -(1)/(x)], "Output"
```
**Output [1331]**
```
CellGroupData[{ Cell[[ RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", -0.1*x^{2}, "]"}], RowBox[{"d", "x"}]}]}]], "Input"
```
**Input [1332]**
```
[ RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", -0.1*x^{2}, "]"}], RowBox[{"d", "x"}]}]}]], "Input"
```
**Output [1333]**
```
["2.802495608198964`"], "Output"
```
**Output [1334]**
```
CellGroupData[{ Cell[[ RowBox[{"Solve", "[", RowBox[{ RowBox[{ RowBox[{"Exp", "[", RowBox[{"-", "tau"}], "]"}], "==", "0"}], ",", "tau"}], "]"}]], "Input"
```
**Input [1335]**
```
[ RowBox[{"Solve", "[", RowBox[{ RowBox[{ RowBox[{"Exp", "[", RowBox[{"-", "tau"}], "]"}], "==", "0"}], ",", "tau"}], "]"}]], "Input"
```
**Output [1336]**
```
[ {}], "Output"
```
**Output [1346]**
```
CellGroupData[{ Cell[[ RowBox[{ SubsuperscriptBox["Int", "0", "oo"], RowBox[{ RowBox[{"Exp", "[", RowBox[{ -(t-I), "*", RowBox[{"Abs", "[", (p+q-r-s), "]"}]}], "]"}], RowBox[{"d", "t"}]}]}]], "Input"
```
**Output [1348]**
```
[ RowBox[{"ConditionalExpression", "[", RowBox[{ FractionBox[ SuperscriptBox["E", RowBox[{"I", " ", RowBox[{"Abs", "[", p+q-r-s, "]"}]}]], RowBox[{"Abs", "[", p+q-r-s, "]"}]], ",", RowBox[{ p+q, "\[NotEqual]", r+s}]}], "]"}]], "Output"
```
**Input [1357]**
```
[ RowBox[{"ListPlot", "[", RowBox[{ RowBox[{"Re", "[", RowBox[{"Exp", "[", RowBox[{ RowBox[{"(", RowBox[{"Delta", "+", I*0.001}], ")"}], "*", "I", "*", "Wghts"}], "]"}], "]"}], ",", RowBox[{"Joined", "->", "True"}], ",", RowBox[{"PlotRange", "->", "All"}]}], "]"}]], "Input"
```
**Input [1360]**
```
[ RowBox[{"ListPlot", "[", RowBox[{ RowBox[{ RowBox[{"Re", "[", RowBox[{"Exp", "[", (10.0+I*0.001)*I*Wghts, "]"}], "]"}], "*", "Wghts"}], ",", RowBox[{"Joined", "->", "True"}], ",", RowBox[{"PlotRange", "->", "All"}]}], "]"}]], "Input"
```
> <So basically the weights need to be flipped for negative energies. vs. positive. Can a mixed solution work? 4 cases>


## [1487] Secularizer (which actually works)
**Input [1488]**
```
CellGroupData[{ Cell["Derivation of The Secularizer.", "Subsection"
```


### [1489] Derivation of The Secularizer.
> <We want to generate the secularizer for 8 indices. Via a python code I verified that his has the following form written recursively in terms of lower - rank tensors as follows :                                 
>       tt        Cond1 = SecularTwoI (p, q, r, s)*SecularTwoI (t, u, v, w)
>                                   Cond2 = (1.0 - Cond1)*SecularTwoI (p, q, v, w)*SecularTwoI (t, u, r, s)
>                                   Cond3 = (1.0 - Cond1)*(1.0 - Cond2)*SecularTwoI (p, q, v, r)*SecularTwoI (t, u, w, s)
>                                   Cond4 = (1.0 - Cond1)*(1.0 - Cond2)*(1.0 - Cond3)*SecularTwoI (p, q, r, w)*SecularTwoI (t, u, v, s)
>                                   Cond5 = (1.0 - Cond1)*(1.0 - Cond2)*(1.0 - Cond3)*(1.0 - Cond4)*SecularTwoI (p, q, v, s)*SecularTwoI (t, u, r, w)
>                                   Cond  = Cond1 + Cond2 + Cond3 + Cond4 + Cond5 + (1.0 - Cond1)*(1.0 - Cond2)*(1.0 - Cond3)*(1.0 - Cond4)*(1.0 - Cond5)*SecularTwoI (p, q, w, s)*SecularTwoI (t, u, v, r)
>                                   
>   Now I' m using mathematica to simplify the remaining >
**Input [1491]**
```
[ RowBox[{ RowBox[{ RowBox[{"STwo", "[", p_, q_, r_, s_, "]"}], " ", ":=", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {r}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {s}}], "}"}]}], "}"}], ",", " ", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {r}}], "}"}]}], "}"}], ",", " ", RowBox[{"{", RowBox[{ {\<Sign\>-1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {r}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {r}}], "}"}]}], "}"}]}], "}"}]}], ";"}]], "Input"
```
**Input [1494]**
```
[ RowBox[{"RecognizableForm", "[", "C2", "]"}]], "Input"
```
**Input [1497]**
```
[ RowBox[{"RecognizableForm", "[", "C4", "]"}]], "Input"
```
**Input [1499]**
```
[ RowBox[{ RowBox[{"Secular8", " ", "=", " ", RowBox[{"AssociativeSimplify", "[", RowBox[{ RowBox[{ RowBox[{"Join", "[", C1, C2, C3, C4, C5, C6, "]"}], " ", "//.", " ", "CombineMDeltas"}], ",", "2"}], "]"}]}], ";"}]], "Input"
```
**Input [1500]**
```
[ RowBox[{"RecognizableForm", "[", "Secular8", "]"}]], "Input"
```
**Input [1502]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SecularizeATerm,FinalByHand, "]"}], ",", "1"}], "]"}], "]"}]], "Input"
```
**Input [1503]**
```
[ RowBox[{"Map", "[", SecularizeATerm,FinalByHand, "]"}]], "Input"
```
**Input [1505]**
```
[ RowBox[{"RecognizableForm", "[", RowBox[{"Flatten", "[", RowBox[{ RowBox[{"Map", "[", SecularizeATerm,FullySimplifiedTCLNoDelta1s, "]"}], ",", "1"}], "]"}], "]"}]], "Input"
```
**Input [1506]**
```
[ RowBox[{ (*This can be used to test the secularizer in python.*), " ", RowBox[{ RowBox[{"StringJoin", "[", RowBox[{"Map", "[", TermCodeString, Secular8, "]"}], "]"}], ";"}]}]], "Input"
```


### [1508] If I do the diagrams by hand I get:
**Output [1509]**
```
CellGroupData[{ Cell[[ RowBox[{"RecognizableForm", "[", "ByHandFromCode", "]"}]], "Input"
```
**Input [1510]**
```
[ RowBox[{"RecognizableForm", "[", "ByHandFromCode", "]"}]], "Input"
```
**Output [1511]**
```
[ RowBox[{"Total", "[", "ByHandFromCode", "]"}]], "Output"
```
**Output [1512]**
```
CellGroupData[{ Cell[[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{ry}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"y"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{wx}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", d^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{vy}}_{{su}}, " ", "VTt"^{{rt}}_{{wx}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{tv}}_{{sy}}, " ", "VTt"^{{rx}}_{{uw}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"y"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"x"}}}]}]], "Input"
```
**Input [1513]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{ry}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"y"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{wx}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", d^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{vy}}_{{su}}, " ", "VTt"^{{rt}}_{{wx}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{tv}}_{{sy}}, " ", "VTt"^{{rx}}_{{uw}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"y"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"x"}}}]}]], "Input"
```
**Output [1514]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{ry}}_{{uw}}, " ", "VTt"^{{tv}}_{{sx}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"y"}}, " ", eta^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{rt}}_{{wx}}, " ", "VTt"^{{vy}}_{{su}}, " ", g^{{"w"}}_{{"v"}}, " ", g^{{"p"}}_{{"y"}}, " ", d^{{"x"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{vy}}_{{su}}, " ", "VTt"^{{rt}}_{{wx}}, " ", g^{{"x"}}_{{"q"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"p"}}_{{"y"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"u"}}_{{"t"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{tv}}_{{sy}}, " ", "VTt"^{{rx}}_{{uw}}, " ", g^{{"u"}}_{{"t"}}, " ", g^{{"w"}}_{{"v"}}, " ", d^{{"y"}}_{{"q"}}, " ", eta^{{"s"}}_{{"r"}}, " ", eta^{{"p"}}_{{"x"}}}]}]], "Output"
```
**Input [1519]**
```
["ByHand"], "Input"
```
**Input [1522]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"ByHand", "[", RowBox[{"[", 2;;3, "]"}], "]"}]}], "]"}], "]"}], " ", RowBox[{"LatexRecognizableForm", "[", RowBox[{"Map", "[", CollapseDeltas,ByHand, "]"}], "]"}]}], "Input"
```
**Input [1525]**
```
[ (*Thinking now of the off diagonal...*)], "Input"
```
> <This could likely be secularized by only the first two indices : >
**Input [1528]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"ByHand", "[", RowBox[{"[", 2;;3, "]"}], "]"}]}], "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"STwo", "[", "\<r\>","\<t\>","\<v\>","\<x\>", "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"DistributeOperatorProduct", "[", RowBox[{ RowBox[{"ByHand", "[", RowBox[{"[", 2;;2, "]"}], "]"}], ",", RowBox[{"STwo", "[", "\<r\>","\<t\>","\<v\>","\<x\>", "]"}]}], "]"}]}], "]"}]}], "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"DistributeOperatorProduct", "[", RowBox[{ RowBox[{"ByHand", "[", RowBox[{"[", 3;;3, "]"}], "]"}], ",", RowBox[{"STwo", "[", "\<r\>","\<t\>","\<v\>","\<x\>", "]"}]}], "]"}]}], "]"}]}], "]"}], "]"}]}], "Input"
```
**Input [1534]**
```
[{ RowBox[{ RowBox[{ RowBox[{"STwo", "[", p_, q_, r_, s_, "]"}], " ", ":=", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {r}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {s}}], "}"}]}], "}"}], ",", " ", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {r}}], "}"}]}], "}"}], ",", " ", RowBox[{"{", RowBox[{ {\<Sign\>-1}, ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {r}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {p}, ",", " ", {s}}], "}"}], ",", " ", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", " ", "1"}], "}"}], ",", " ", {q}, ",", " ", {r}}], "}"}]}], "}"}]}], "}"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"ByHand3", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"DistributeOperatorProduct", "[", RowBox[{"ByHand", ",", RowBox[{"STwo", "[", "\<r\>","\<t\>","\<v\>","\<x\>", "]"}]}], "]"}]}], "]"}]}], "]"}]}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"TCLSimp4", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"Map", "[", DiagonalGamma1,ByHand3, "]"}]}], "]"}]}], "]"}]}], "]"}]}], "Input"
```
**Input [1539]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{pr}}, " ", "VTt"^{{pr}}_{{st}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}, " ", eta^{{"t"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{pr}}, " ", "VTt"^{{pr}}_{{st}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}, " ", eta^{{"t"}}_{{"t"}}}]}]], "Input"
```
**Output [1540]**
```
[ RowBox[{ RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "+", RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{pr}}, " ", "VTt"^{{pr}}_{{st}}, " ", g^{{"s"}}_{{"s"}}, " ", g^{{"t"}}_{{"t"}}, " ", eta^{{"p"}}_{{"p"}}, " ", eta^{{"r"}}_{{"r"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{pr}}_{{st}}, " ", "VTt"^{{st}}_{{pr}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}, " ", eta^{{"t"}}_{{"t"}}}], "-", RowBox[{ (1)/(2), " ", "VTs"^{{st}}_{{pr}}, " ", "VTt"^{{pr}}_{{st}}, " ", g^{{"p"}}_{{"p"}}, " ", g^{{"r"}}_{{"r"}}, " ", eta^{{"s"}}_{{"s"}}, " ", eta^{{"t"}}_{{"t"}}}]}]], "Output"
```
**Input [1541]**
```
[ RowBox[{"Now", ",", " ", first and second terms are conjugates of each other, ",", " ", RowBox[{ and the blocking factors for s, "->", RowBox[{ p and p, "->", s are correctly opposite in character. So everything looks correct.}]}]}]], "Input"
```
> For this to make sense the "Rates." must all be positive.
**Input [1551]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"DistributeOperatorProduct", "[", RowBox[{ RowBox[{"{", RowBox[{"ByHand", "[", RowBox[{"[", "1", "]"}], "]"}], "}"}], ",", RowBox[{"STwo", "[", "\<r\>","\<t\>","\<v\>","\<x\>", "]"}]}], "]"}]}], "]"}]}], "]"}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"DistributeOperatorProduct", "[", RowBox[{ RowBox[{"{", RowBox[{"ByHand", "[", RowBox[{"[", "4", "]"}], "]"}], "}"}], ",", RowBox[{"STwo", "[", "\<r\>","\<t\>","\<v\>","\<x\>", "]"}]}], "]"}]}], "]"}]}], "]"}], "]"}]}], "Input"
```
**Input [1560]**
```
[""], "Input"
```
**Input [1562]**
```
[{ RowBox[{"RecognizableForm", "[", RowBox[{"ByHand3", "=", RowBox[{"Map", "[", RowBox[{"CollapseMultiDeltas", ",", RowBox[{"Map", "[", RowBox[{"CollapseDeltas", ",", RowBox[{"DistributeOperatorProduct", "[", RowBox[{ RowBox[{"{", RowBox[{"ByHand", "[", RowBox[{"[", "2", "]"}], "]"}], "}"}], ",", RowBox[{"STwo", "[", "\<r\>","\<t\>","\<v\>","\<x\>", "]"}]}], "]"}]}], "]"}]}], "]"}]}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"ByHand5", "=", RowBox[{ RowBox[{"SpinLabelTerms", "[", "ByHand3", "]"}], "//.", {VsjkVtjkBSpinsBSpint}}]}], "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"ByHand6", "=", RowBox[{"AssociativeSimplify", "[", RowBox[{"Map", "[", RelabelTerm,ByHand5, "]"}], "]"}]}], "]"}], " ", AllVts={};, " ", RowBox[{"StringJoin", "[", RowBox[{"Map", "[", TermCodePUp,ByHand6, "]"}], "]"}], " ", RowBox[{"Union", "[", "AllVts", "]"}]}], "Input"
```
**Input [1569]**
```
[ RowBox[{"ByHand5", "[", RowBox[{"[", "1", "]"}], "]"}]], "Input"
```
**Input [1572]**
```
[{ RowBox[{"FormatATerm", "[", RowBox[{"ByHand5", "[", RowBox[{"[", "1", "]"}], "]"}], "]"}], " ", RowBox[{"FormatATerm", "[", RowBox[{"RelabelTerm", "[", RowBox[{"ByHand5", "[", RowBox[{"[", "1", "]"}], "]"}], "]"}], "]"}]}], "Input"
```
**Input [1576]**
```
[ RowBox[{"Length", "[", "FinalTCL2", "]"}]], "Input"
```
`[Output 1578: graphics omitted]`


## [1579] Applications without new Code
**Input [1580]**
```
CellGroupData[{ Cell["Two Particle Products (No functions defined)", "Subsection"
```


### [1581] Two Particle Products (No functions defined)
> <To Work out the two - particle product takes 65 seconds on my i7 (2011) laptop. >
**Input [1583]**
```
[{ RowBox[{ RowBox[{"fac1", "=", RowBox[{"NRta2", "[", {\<b\>\<d\>},{\<c\>\<e\>}, "]"}]}], ";"}], " ", RowBox[{ RowBox[{"fac2", "=", RowBox[{"NRta2", "[", {\<f\>\<h\>},{\<g\>\<i\>}, "]"}]}], ";"}]}], "Input"
```
**Input [1584]**
```
[{ RowBox[{"Timing", "[", RowBox[{ RowBox[{"tpp", "=", RowBox[{"AllStepsOfProductRule", "[", fac1,fac2, "]"}]}], ";"}], "]"}], " ", RowBox[{"RecognizableForm", "[", "tpp", "]"}]}], "Input"
```
**Input [1585]**
```
[ RowBox[{"RecognizableForm", "[", tpp//.{NoCumulantsNormalVaccuumRule}, "]"}]], "Input"
```
`[Output 1586: graphics omitted]`


### [1587] Derivation of BBGKY Heirarchy up to two.
> For now just the Two Particle H with one particle operator.
**Input [1589]**
```
[{ RowBox[{ RowBox[{"fac1", "=", RowBox[{"NRta1", "[", {"\<p\>"},{"\<q\>"}, "]"}]}], ";"}], " ", RowBox[{ RowBox[{"fac2", "=", RowBox[{"NRta2", "[", {\<c\>\<d\>},{\<e\>\<f\>}, "]"}]}], ";"}]}], "Input"
```
**Output [1590]**
```
CellGroupData[{ Cell[[{ RowBox[{ RowBox[{ RowBox[{"a1v2", "=", RowBox[{"AllStepsOfProductRule", "[", fac1,fac2, "]"}]}], ";"}], " ", (*1*)}], " ", RowBox[{ RowBox[{ RowBox[{"v2a1", "=", RowBox[{"AllStepsOfProductRule", "[", fac2,fac1, "]"}]}], ";"}], " ", (*1*)}], " ", RowBox[{ RowBox[{"com", "=", RowBox[{"Join", "[", RowBox[{"a1v2", ",", RowBox[{"MultiplyScalar", "[", v2a1,-1, "]"}]}], "]"}]}], ";"}], " ", RowBox[{ RowBox[{"RecognizableForm", "[", "a1v2", "]"}], ";"}], " ", RowBox[{ RowBox[{"com2", "=", RowBox[{"MultiplyTermsByTensor", "[", RowBox[{"com", ",", RowBox[{"{", RowBox[{"{", {\<VTs\>2},{\<e\>\<f\>},{\<c\>\<d\>}, "}"}], "}"}]}], "]"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", com3=com2//.NormalVaccuumRule, "]"}], " ", RowBox[{"RecognizableForm", "[", RowBox[{"com4", "=", RowBox[{"Map", "[", CollapseDeltas,com3//.EtaToDelta, "]"}]}], "]"}]}], "Input"
```
**Output [1592]**
```
[ RowBox[{ RowBox[{ "VTs"^{{ef}}_{{cd}}, " ", g^{{"d"}}_{{"q"}}, " ", lam^{{cp}}_{{ef}}}], "+", RowBox[{ "VTs"^{{ef}}_{{cd}}, " ", eta^{{"d"}}_{{"q"}}, " ", lam^{{cp}}_{{ef}}}], "-", RowBox[{ "VTs"^{{ef}}_{{cd}}, " ", g^{{"c"}}_{{"q"}}, " ", lam^{{dp}}_{{ef}}}], "-", RowBox[{ "VTs"^{{ef}}_{{cd}}, " ", eta^{{"c"}}_{{"q"}}, " ", lam^{{dp}}_{{ef}}}], "-", RowBox[{ "VTs"^{{ef}}_{{cd}}, " ", g^{{"p"}}_{{"f"}}, " ", lam^{{cd}}_{{eq}}}], "-", RowBox[{ "VTs"^{{ef}}_{{cd}}, " ", eta^{{"p"}}_{{"f"}}, " ", lam^{{cd}}_{{eq}}}], "+", RowBox[{ "VTs"^{{ef}}_{{cd}}, " ", g^{{"p"}}_{{"e"}}, " ", lam^{{cd}}_{{fq}}}], "+", RowBox[{ "VTs"^{{ef}}_{{cd}}, " ", eta^{{"p"}}_{{"e"}}, " ", lam^{{cd}}_{{fq}}}]}]], "Output"
```
**Output [1593]**
```
[ RowBox[{ RowBox[{ "VTs"^{{ef}}_{{cq}}, " ", lam^{{cp}}_{{ef}}}], "+", RowBox[{ "VTs"^{{ef}}_{{dq}}, " ", lam^{{dp}}_{{ef}}}], "-", RowBox[{ "VTs"^{{ep}}_{{cd}}, " ", lam^{{cd}}_{{eq}}}], "-", RowBox[{ "VTs"^{{fp}}_{{cd}}, " ", lam^{{cd}}_{{fq}}}]}]], "Output"
```
**Input [1594]**
```
[ RowBox[{"(*", RowBox[{ RowBox[{ "The", " ", "last", " ", "four", " ", "terms", " ", "can", " ", "be", " ", "significantly", " ", "simplified", " ", "by", " ", "collapsing", " ", "g"}], "+", RowBox[{ "eta", " ", "into", " ", "delta", " ", "yielding", " ", "BBGKY1", " ", "recall", " ", "the", " ", "factor", " ", "of", " ", 1/2, " ", "in", " ", "the", " ", "V", " ", operator. c.f. Tohyama, " ", "and", " ", "Shuck"}]}], "*)"}]], "Input"
```
`[Input 1595: graphics omitted]`
> BBGKY2
**Input [1598]**
```
[{ RowBox[{ RowBox[{ RowBox[{"a1v2", "=", RowBox[{"AllStepsOfProductRule", "[", fac1,fac2, "]"}]}], ";"}], " ", (*1*)}], " ", RowBox[{ RowBox[{ RowBox[{"v2a1", "=", RowBox[{"AllStepsOfProductRule", "[", fac2,fac1, "]"}]}], ";"}], " ", (*1*)}], " ", RowBox[{ RowBox[{"com", "=", RowBox[{"Join", "[", RowBox[{"a1v2", ",", RowBox[{"MultiplyScalar", "[", v2a1,-1, "]"}]}], "]"}]}], ";"}], " ", RowBox[{ RowBox[{"RecognizableForm", "[", "a1v2", "]"}], ";"}]}], "Input"
```
`[Input 1602: graphics omitted]`
**Input [1603]**
```
[ Which looks really imposing, but again comparing with tohyama.], "Input"
```
**Input [1604]**
```
[ The Cumulant decomposition of rho3 is useful:], "Input"
```
**Output [1605]**
```
CellGroupData[{ Cell[[{ RowBox[{ RowBox[{"decompgamma3", "=", RowBox[{"Gd3", "[", {\<p\>\<q\>\<r\>},{\<s\>\<t\>\<u\>}, "]"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", "decompgamma2", "]"}], " ", RowBox[{"RecognizableForm", "[", "decompgamma3", "]"}]}], "Input"
```
**Input [1606]**
```
[{ RowBox[{ RowBox[{"decompgamma3", "=", RowBox[{"Gd3", "[", {\<p\>\<q\>\<r\>},{\<s\>\<t\>\<u\>}, "]"}]}], ";"}], " ", RowBox[{"RecognizableForm", "[", "decompgamma2", "]"}], " ", RowBox[{"RecognizableForm", "[", "decompgamma3", "]"}]}], "Input"
```
**Output [1607]**
```
[ RowBox[{ RowBox[{ RowBox[{"-", g^{{"d"}}_{{"c"}}}], " ", g^{{"b"}}_{{"e"}}}], "+", RowBox[{ g^{{"b"}}_{{"c"}}, " ", g^{{"d"}}_{{"e"}}}], "+", lam^{{bd}}_{{ce}}}]], "Output"
```
**Output [1608]**
```
[ RowBox[{ RowBox[{ RowBox[{"-", g^{{"r"}}_{{"s"}}}], " ", g^{{"q"}}_{{"t"}}, " ", g^{{"p"}}_{{"u"}}}], "+", RowBox[{ g^{{"q"}}_{{"s"}}, " ", g^{{"r"}}_{{"t"}}, " ", g^{{"p"}}_{{"u"}}}], "+", RowBox[{ g^{{"r"}}_{{"s"}}, " ", g^{{"p"}}_{{"t"}}, " ", g^{{"q"}}_{{"u"}}}], "-", RowBox[{ g^{{"p"}}_{{"s"}}, " ", g^{{"r"}}_{{"t"}}, " ", g^{{"q"}}_{{"u"}}}], "-", RowBox[{ g^{{"q"}}_{{"s"}}, " ", g^{{"p"}}_{{"t"}}, " ", g^{{"r"}}_{{"u"}}}], "+", RowBox[{ g^{{"p"}}_{{"s"}}, " ", g^{{"q"}}_{{"t"}}, " ", g^{{"r"}}_{{"u"}}}], "+", RowBox[{ g^{{"r"}}_{{"u"}}, " ", lam^{{pq}}_{{st}}}], "-", RowBox[{ g^{{"q"}}_{{"u"}}, " ", lam^{{pr}}_{{st}}}], "+", RowBox[{ g^{{"p"}}_{{"u"}}, " ", lam^{{qr}}_{{st}}}], "-", RowBox[{ g^{{"r"}}_{{"t"}}, " ", lam^{{pq}}_{{su}}}], "+", RowBox[{ g^{{"q"}}_{{"t"}}, " ", lam^{{pr}}_{{su}}}], "-", RowBox[{ g^{{"p"}}_{{"t"}}, " ", lam^{{qr}}_{{su}}}], "+", RowBox[{ g^{{"r"}}_{{"s"}}, " ", lam^{{pq}}_{{tu}}}], "-", RowBox[{ g^{{"q"}}_{{"s"}}, " ", lam^{{pr}}_{{tu}}}], "+", RowBox[{ g^{{"p"}}_{{"s"}}, " ", lam^{{qr}}_{{tu}}}], "+", lam^{{pqr}}_{{stu}}}]], "Output"
```


### [1610] Triple Commutator of 1-2-2 particle operators required for Correlated dynamics Theory without projectors
**Input [1611]**
```
[{ RowBox[{"an", "[", 2,17, "]"}], " ", RowBox[{"an", "[", 2,21, "]"}]}], "Input"
```
**Input [1613]**
```
[ RowBox[{"RecognizableForm", "[", "temp1", "]"}]], "Input"
```
> <This should require 5 - particle operators for the first term. although such products shouldn' t contribute to the final result.
> This requires ~12 min on my desktop. We need to basically debug these against a diagrammatic version of the same expressions. >
**Input [1615]**
```
[ RowBox[{"Timing", "[", RowBox[{ RowBox[{"rv1v2", "=", RowBox[{"AllStepsOfProductRule1PV", "[", temp1,fac3, "]"}]}], ";"}], "]"}]], "Input"
```
**Input [1616]**
```
[ RowBox[{ RowBox[{"RecognizableForm", "[", rv1v2//.{OneParticleVaccuum1OneParticleVaccuum2}, "]"}], "//", "FullSimplify"}]], "Input"
```
> <Close but there is a sign error in the definition of ta or gamma. These factors of two correspond to terms which should cancel but don' t. Moreover there should be a symmetry between terms that have 3 - etas and terms that have 3 - gammas. >
**Input [1618]**
```
[""], "Input"
```
**Input [1619]**
```
[ RowBox[{ RowBox[{"v1rv2", "=", RowBox[{"AllStepsOfProductRule1PV", "[", RowBox[{ RowBox[{"AllStepsOfProductRule", "[", fac2,fac1, "]"}], ",", "fac3"}], "]"}]}], ";"}]], "Input"
```
**Input [1620]**
```
[ RowBox[{ RowBox[{"v2v1r", "=", RowBox[{"AllStepsOfProductRule1PV", "[", RowBox[{"fac3", ",", RowBox[{"AllStepsOfProductRule", "[", fac2,fac1, "]"}]}], "]"}]}], ";"}]], "Input"
```
**Input [1621]**
```
[ RowBox[{ RowBox[{"v2rv1", "=", RowBox[{"AllStepsOfProductRule1PV", "[", RowBox[{"fac3", ",", RowBox[{"AllStepsOfProductRule", "[", fac1,fac2, "]"}]}], "]"}]}], ";"}]], "Input"
```
**Input [1622]**
```
[{ RowBox[{ RowBox[{"AllTerms", "=", RowBox[{"Join", "[", RowBox[{"rv1v2", ",", "v2v1r", ",", RowBox[{"MultiplyScalar", "[", v1rv2,-1, "]"}], ",", RowBox[{"MultiplyScalar", "[", v2rv1,-1, "]"}]}], "]"}]}], ";"}], " ", RowBox[{"CorrelatedForm", "=", RowBox[{ RowBox[{ RowBox[{"RecognizableForm", "[", AllTerms//.{OneParticleVaccuum1OneParticleVaccuum2}, "]"}], "//", "Simplify"}], "//", "Expand"}]}]}], "Input"
```
**Input [1623]**
```
[ RowBox[{"sAllTerms", "=", RowBox[{"AssociativeSimplify", "[", "AllTerms", "]"}]}]], "Input"
```
**Input [1624]**
```
[ RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{ {\<Sign\>-1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<x\>"}, ",", {"\<u\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<t\>"}, ",", {"\<w\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<r\>"}, ",", {"\<y\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<v\>"}, ",", {"\<q\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<p\>"}, ",", {"\<s\>"}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<x\>"}, ",", {"\<u\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<r\>"}, ",", {"\<w\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<t\>"}, ",", {"\<y\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<v\>"}, ",", {"\<q\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<p\>"}, ",", {"\<s\>"}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<v\>"}, ",", {"\<u\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<t\>"}, ",", {"\<w\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<r\>"}, ",", {"\<y\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<x\>"}, ",", {"\<q\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<p\>"}, ",", {"\<s\>"}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<Sign\>-1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<v\>"}, ",", {"\<u\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<r\>"}, ",", {"\<w\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<t\>"}, ",", {"\<y\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<x\>"}, ",", {"\<q\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<p\>"}, ",", {"\<s\>"}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<x\>"}, ",", {"\<u\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<t\>"}, ",", {"\<w\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<p\>"}, ",", {"\<y\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<r\>"}, ",", {"\<q\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<v\>"}, ",", {"\<s\>"}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<Sign\>-1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<x\>"}, ",", {"\<u\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<p\>"}, ",", {"\<w\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<t\>"}, ",", {"\<y\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<r\>"}, ",", {"\<q\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<v\>"}, ",", {"\<s\>"}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<Sign\>-1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<x\>"}, ",", {"\<u\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<r\>"}, ",", {"\<w\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<p\>"}, ",", {"\<y\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<t\>"}, ",", {"\<q\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<v\>"}, ",", {"\<s\>"}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<x\>"}, ",", {"\<u\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<p\>"}, ",", {"\<w\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<r\>"}, ",", {"\<y\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<t\>"}, ",", {"\<q\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<v\>"}, ",", {"\<s\>"}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<Sign\>-1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<v\>"}, ",", {"\<u\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<t\>"}, ",", {"\<w\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<p\>"}, ",", {"\<y\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<r\>"}, ",", {"\<q\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ",", {"\<x\>"}, ",", {"\<s\>"}}], "}"}]}], "}"}], ",", RowBox[{"{", RowBox[{ {\<Sign\>1}, ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<v\>"}, ",", {"\<u\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<p\>"}, ",", {"\<w\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"g", ",", "1"}], "}"}], ",", {"\<t\>"}, ",", {"\<y\>"}}], "}"}], ",", RowBox[{"{", RowBox[{ RowBox[{"{", RowBox[{"d", ",", "1"}], "}"}], ...[TRUNC 59028]
```
**Input [1625]**
```
[ RowBox[{"Length", "[", RowBox[{"AssociativeSimplify", "[", "AllTerms", "]"}], "]"}]], "Input"
```
> <Wow ... that actually looks correct :) The next step is that we need to multiply all these terms by appropriate interaction picture expressions and then write automatically generated code for these things. >
**Input [1627]**
```
[ RowBox[{"Length", "[", "CorrelatedForm", "]"}]], "Input"
```
`[Output 1628: graphics omitted]`


### [1629] Triple Commutator of 2-1-1 particle operators required for Correlated dynamics Theory
**Input [1630]**
```
[{ RowBox[{"an", "[", 2,1, "]"}], " ", RowBox[{"an", "[", 1,5, "]"}], " ", RowBox[{"an", "[", 1,7, "]"}]}], "Input"
```
**Input [1632]**
```
[ RowBox[{"RecognizableForm", "[", "temp1", "]"}]], "Input"
```
> <This should require 5 - particle operators for the first term. although such products shouldn' t contribute to the final result.
> This requires ~12 min on my desktop. We need to basically debug these against a diagrammatic version of the same expressions. >
**Input [1634]**
```
[ RowBox[{"Timing", "[", RowBox[{ RowBox[{"rv1v2", "=", RowBox[{"AllStepsOfProductRule", "[", temp1,fac3, "]"}]}], ";"}], "]"}]], "Input"
```
**Input [1635]**
```
[ RowBox[{ RowBox[{"RecognizableForm", "[", "rv1v2", "]"}], "//", "FullSimplify"}]], "Input"
```
**Input [1637]**
```
[ RowBox[{ RowBox[{"v2v1r", "=", RowBox[{"AllStepsOfProductRule", "[", RowBox[{"fac3", ",", RowBox[{"AllStepsOfProductRule", "[", fac2,fac1, "]"}]}], "]"}]}], ";"}]], "Input"
```
**Input [1638]**
```
[ RowBox[{ RowBox[{"v2rv1", "=", RowBox[{"AllStepsOfProductRule", "[", RowBox[{"fac3", ",", RowBox[{"AllStepsOfProductRule", "[", fac1,fac2, "]"}]}], "]"}]}], ";"}]], "Input"
```
**Input [1639]**
```
[{ RowBox[{ RowBox[{"AllTerms", "=", RowBox[{"Join", "[", RowBox[{"rv1v2", ",", "v2v1r", ",", RowBox[{"MultiplyScalar", "[", v1rv2,-1, "]"}], ",", RowBox[{"MultiplyScalar", "[", v2rv1,-1, "]"}]}], "]"}]}], ";"}], " ", RowBox[{"CorrelatedForm", "=", RowBox[{ RowBox[{"RecognizableForm", "[", AllTerms//.NormalVaccuumRule, "]"}], "//", "FullSimplify"}]}]}], "Input"
```
> < The next step is that we need to multiply all these terms by appropriate interaction picture expressions and then write automatically generated code for these things. >
**Output [1641]**
```
["$Aborted"], "Output"
```
**Input [1642]**
```
[""], "Input"
```
`[Output 1643: graphics omitted]`
**Input [1644]**
```
[{ RowBox[{ RowBox[{"Cs", "=", RowBox[{"1.401", "*", RowBox[{"Exp", "[", RowBox[{"I", "*", RowBox[{"(", RowBox[{"Pi", "-", "0.00001"}], ")"}]}], "]"}]}]}], ";"}], " ", Z=Cs;, " ", RowBox[{ RowBox[{"Seq", "=", RowBox[{ RowBox[{"Reap", "[", RowBox[{"Do", "[", RowBox[{ RowBox[{"Sow", "[", Z=Z^2+Cs, "]"}], ",", {i11000}}], "]"}], "]"}], "[", RowBox[{"[", 2,1, "]"}], "]"}]}], ";"}], " ", RowBox[{"ListPlot", "[", RowBox[{"Re", "[", "Seq", "]"}], "]"}], " ", RowBox[{"ListPlot", "[", RowBox[{"Im", "[", "Seq", "]"}], "]"}]}], "Input"
```
`[Output 1645: graphics omitted]`
`[Output 1646: graphics omitted]`


