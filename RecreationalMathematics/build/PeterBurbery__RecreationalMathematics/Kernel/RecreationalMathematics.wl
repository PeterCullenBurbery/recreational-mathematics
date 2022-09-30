(* ::Package:: *)

BeginPackage["PeterBurbery`RecreationalMathematics`"];

(* Declare your packages public symbols here. *)

FivePointConic;

Begin["`Private`"];

(* Define your public and private symbols here. *)
FivePointConic//ClearAll
FivePointConic[pts_]/;MatrixQ[pts]&&Dimensions[pts]==={5,2}:=FivePointConic[pts,{\[FormalX],\[FormalY]}];FivePointConic[pts_,{xx_,yy_}]/;MatrixQ[pts]&&Dimensions[pts]==={5,2}:=Times@@Apply[Power,DeleteCases[FactorList[Det[PadLeft[Apply[Function[{x,y},{x,y,x^2,x y,y^2}],Append[pts,{xx,yy}],{1}],{6,6},1]]],{_?NumericQ,_}],{1}]

End[]; (* End `Private` *)

EndPackage[];
