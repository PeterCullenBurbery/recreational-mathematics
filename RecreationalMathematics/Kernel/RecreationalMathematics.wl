(* ::Package:: *)

BeginPackage["PeterBurbery`RecreationalMathematics`"];

(* Declare your packages public symbols here. *)
EvenPermutations;
FivePointConic;
NinePointCubic;
EulerLinePoints;
NinePointQuadric;
BalancedTernary;
CatalanUnrank;
AllBalancedGroupingSymbols;
Derangements;
IntegralNumberQ;
DyckPaths;
DiagonalWalkPlot;
ParenthesizedExpressions;
PermutationGraph;
FullBinaryTrees;
FindRegularPolygonTriangulations;
Multichoose;
Begin["`Private`"];

(* Define your public and private symbols here. *)
FivePointConic//ClearAll
FivePointConic[pts_]/;MatrixQ[pts]&&Dimensions[pts]==={5,2}:=FivePointConic[pts,{\[FormalX],\[FormalY]}];FivePointConic[pts_,{xx_,yy_}]/;MatrixQ[pts]&&Dimensions[pts]==={5,2}:=Times@@Apply[Power,DeleteCases[FactorList[Det[PadLeft[Apply[Function[{x,y},{x,y,x^2,x y,y^2}],Append[pts,{xx,yy}],{1}],{6,6},1]]],{_?NumericQ,_}],{1}]

NinePointCubic//ClearAll
NinePointCubic[pts_]/;MatrixQ[pts]&&Dimensions[pts]==={9,2}:=NinePointCubic[pts,{\[FormalX],\[FormalY]}];NinePointCubic[pts_,{xx_,yy_}]/;MatrixQ[pts]&&Dimensions[pts]==={9,2}:=Times@@Apply[Power,DeleteCases[FactorList[Det[PadLeft[Apply[Function[{x,y},{x,y,x^2,x y,y^2,x^3,x^2 y,x y^2,y^3}],Append[pts,{xx,yy}],{1}],{10,10},1]]],{_?NumericQ,_}],{1}]


EulerLinePoints//ClearAll;

EulerLinePoints[vert_]/;MatrixQ[vert,Internal`RealValuedNumericQ]:=Module[{o,g,n,h},
Which[Dimensions[vert]==={3,2},
o=First[Circumsphere[vert]],
Dimensions[vert]==={3,3},
o=ResourceFunction["Circumcircle3D"][vert,"Center"],
True,Return[$Failed,Module]];
g=Mean[vert];
h=3 g -2 o;
n=Mean[{o,h}];
{o,g,n,h}];

EulerLinePoints[Triangle[vert_]]:=If[MatrixQ[vert],EulerLinePoints[vert],EulerLinePoints/@vert]

NinePointQuadric//ClearAll;
NinePointQuadric[pts_]/;MatrixQ[pts]&&Dimensions[pts]==={9,3}:=NinePointQuadric[pts,{\[FormalX],\[FormalY],\[FormalZ]}];
NinePointQuadric[pts_,{xx_,yy_,zz_}]/;MatrixQ[pts]&&Dimensions[pts]==={9,3}:=Times@@Apply[Power,DeleteCases[FactorList[Det[PadLeft[Apply[Function[{x,y,z},{x,x^2,y,x y,y^2,z,x z,y z,z^2}],Append[pts,{xx,yy,zz}],{1}],{10,10},1]]],{_?NumericQ,_}],{1}]



BalancedTernary//ClearAll

BalancedTernary[n:(_Integer|_Row)]:=Module[{list, len},
If[Head[n]===Row,
list=Normal[n]/.\!\(\*UnderscriptBox[\(1\), \(_\)]\)-> -1;
FromDigits[list,3],
len=Ceiling[Log[3,1+2Abs[n]]];
Row[(IntegerDigits[-n- Quotient[3^len-1,2],3,len]-1)/.{{}->{0},-1->\!\(\*UnderscriptBox[\(1\), \(_\)]\)}]]]

CatalanUnrank//ClearAll
CatalanUnrank[n_,rank_]:= 
Module[{ lo=0,y=0, a=Table[0,2 n],m},
Do[
m=Binomial[2 n -x, n-(x+y+1)/2]-Binomial[2 n -x, n-1-(x+y+1)/2];
(*these terms make the Catalan triangle, or ballot numbers*)
If[rank<= lo+m-1,
y=y+1;
a[[x]]=0,
lo=lo+m;
y=y-1;
a[[x]]=1],
{x,1,2 n}];
a]



Derangements//ClearAll
Derangements[n_?IntegerQ]:=With[
{perms=Permutations[Range@n]},
Pick[perms,Length/@PermutationSupport/@perms,n]
]/;n>=0

AllBalancedGroupingSymbols//ClearAll
AllBalancedGroupingSymbols[{openingsymbol_,closingsymbol_},numberofgroups_]:=StringJoin[ResourceFunction["CatalanUnrank"][numberofgroups,#]/.{0->openingsymbol,1->closingsymbol}]&/@Range[0,CatalanNumber[numberofgroups]-1]

IntegralNumberQ//ClearAll
IntegralNumberQ[n_?NumericQ]:=Floor[n]==Ceiling[n]

DyckPathsUnrankFunction//ClearAll
DyckPaths//ClearAll
DyckPathsUnrankFunction[n_,rank_]:=Module[{ lo=0,y=0, a=Table[0,2 n],m},
Do[
m=Binomial[2 n -x, n-(x+y+1)/2]-Binomial[2 n -x, n-1-(x+y+1)/2];
(*these terms make the Catalan triangle, or ballot numbers*)
If[rank<= lo+m-1,
y=y+1;
a[[x]]=1,
lo=lo+m;
y=y-1;
a[[x]]=-1],
{x,1,2 n}];
a]
DyckPaths[n_,options:OptionsPattern[{Ticks->None,Axes->{True,False},ListLinePlot}]]:=ListLinePlot[Prepend[Accumulate[
DyckPathsUnrankFunction[n,#]],0],PlotRange->n,Ticks->OptionValue[Ticks],Axes->OptionValue[Axes],options]&/@Range[0,CatalanNumber[n]-1]

DiagonalWalkPlot//ClearAll
DiagonalWalkPlot[n_,{horizontaldirectionsign_,
verticaldirectionsign_},
options:OptionsPattern[{Ticks->None,ListLinePlot}]]:=
Flatten[Function[sign,(ListLinePlot[{Prepend[
Accumulate[(sign (CatalanUnrank[n,#1]/.{0->-1}))
/. {-1->{0,verticaldirectionsign},1->{horizontaldirectionsign,0}}],
{0,0}],Transpose[{horizontaldirectionsign Range[0,n],
verticaldirectionsign Range[0,n]}]},Ticks->OptionValue[Ticks],options]&)/@Range[0,CatalanNumber[n]-1]]/@{1,-1}]

ParenthesizedExpressions//ClearAll;
(*ParenthesizedExpressions[n_?IntegerQ]:=Block[{f},SetAttributes[f,{Flat,OneIdentity}];
e:CirclePlus[___,_List,___]:=Distribute[Unevaluated[e],
List];
f[Sequence@@(ToExpression@Array[Subscript["x",##]&,n])]//.{f[x__]:>ReplaceList[f[x],f[u_,v_]:>CirclePlus[u,v]]}//Flatten]*)
(*this doesn't work because it adds the full context like *)
ParenthesizedExpressions[n_,variables_]:=Block[{f},SetAttributes[f,{Flat,OneIdentity}];
e:CirclePlus[___,_List,___]:=Distribute[Unevaluated[e],
List];
f[Sequence@@variables]//.{f[x__]:>ReplaceList[f[x],f[u_,v_]:>CirclePlus[u,v]]}//Flatten]/;Length[variables]==n


PermutationGraph//ClearAll
Options[PermutationGraph]=Options[Graph];

PermutationGraph[p_?PermutationListQ,opts:OptionsPattern[]]:=Module[{q=InversePermutation[p]//PermutationList},RelationGraph[((#1<#2&&q[[#1]]>q[[#2]])||(#1>#2&&q[[#1]]<q[[#2]]))&,Range[Length[q]],opts]
]


EvenPermutations//ClearAll
EvenPermutations[list_List, count_Integer]:= Module[{ag,len},
len=Length[list];
ag=AlternatingGroup[len];
Permute[list,#]&/@RandomPermutation[ag,count]
]
EvenPermutations[list_List]:= Module[{rmax=20,ag,len},
len=Length[list];
ag=AlternatingGroup[len];
If[len<10,
Permute[list,ag],
Permute[list,#]&/@RandomPermutation[ag,rmax]
]
]


FullBinaryTrees[n_]:=UnlabeledTree[ExpressionTree[#]]&/@ParenthesizedExpressions[n,ConstantArray[1,n]]



FindRegularPolygonTriangulations//ClearAll
crosscheck[{{aa_,bb_},{cc_,dd_}}] :=
 If[Length[Union[{aa,bb,cc,dd}]]==3, True,
If[bb<cc || bb>dd,True,False]]
linechecker[zz_List] := If[Union[Map[crosscheck[#]&,Subsets[zz,{2}]]]=={True},True,False]
TriangulationsFunction//ClearAll
TriangulationsFunction[k_]:=First@Table[Last/@Select[Sort[Map[{linechecker[#],#}&,Sort/@Subsets[Sort[Drop[Flatten[Reverse[Table[{a,b},{a,1,sides-2},{b,a+2,sides}]],1],-1]],{sides-3}]]],First[#]==True&],{sides,{k}}];
FindRegularPolygonTriangulations[k_]:=Block[{gg,pts,lines},gg=k;pts=Table[N[{Cos[x],Sin[x]}],{x,2\[Pi]/gg,2Pi,2Pi/gg}];lines=TriangulationsFunction[k];
Table[Graphics[{Line/@Map[pts[[#]]&,lines[[n]]],Line[Append[pts,pts[[1]]]]},ImageSize->Tiny],{n,1,Length[lines]}]]

Multichoose//ClearAll
Multichoose[n_,k_]:=Binomial[n+k-1,k]


End[]; (* End `Private` *)

EndPackage[];
