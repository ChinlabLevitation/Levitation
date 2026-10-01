(* Build and evaluate a notebook from cells.wl (made by mbuild.py).
   Usage: wolframscript -file buildnb.wl cells.wl out.nb out.pdf
   Code cells run in Global`; this builder lives in Builder` so that ClearAll["Global`*"] in the notebook cannot touch it.
   Results of the form VR[ref, claim, status, detail] (or Columns containing them) become formatted verification cells. *)
Begin["Builder`"];
{cellsFile, outNB, outPDF} = Take[$ScriptCommandLine, -3];
Get[cellsFile];

parseBoxes[s_String] := First@FrontEndExecute[FrontEnd`UndocumentedTestFEParserPacket[s, True]];

greekRules = Join[{"Sqrt[pi]" -> "\[Sqrt]\[Pi]", "Sqrt[Pi]" -> "\[Sqrt]\[Pi]", "grad " -> "\[Del]"},
   Map[RegularExpression["\\b" <> #[[1]] <> "(?![A-Za-z])"] -> #[[2]] &,
    {{"alpha", "\[Alpha]"}, {"gamma", "\[Gamma]"}, {"delta", "\[Delta]"}, {"epsilon", "\[CurlyEpsilon]"}, {"eps", "\[CurlyEpsilon]"},
     {"eta", "\[Eta]"}, {"theta", "\[Theta]"}, {"kappa", "\[Kappa]"}, {"lambda", "\[Lambda]"}, {"mu", "\[Mu]"}, {"nu", "\[Nu]"},
     {"xi", "\[Xi]"}, {"pi", "\[Pi]"}, {"rho", "\[Rho]"}, {"sigma", "\[Sigma]"}, {"tau", "\[Tau]"}, {"phi", "\[Phi]"}, {"psi", "\[Psi]"},
     {"omega", "\[Omega]"}, {"Gamma", "\[CapitalGamma]"}, {"Delta", "\[CapitalDelta]"}, {"Lambda", "\[CapitalLambda]"},
     {"Omega", "\[CapitalOmega]"}, {"Pi", "\[CapitalPi]"}}]];
letters = "A-Za-z\\x{0391}-\\x{03C9}";
textRuns[str0_String] := Module[{str = StringReplace[str0, greekRules]},
   StringSplit[str, {
      RegularExpression["([" <> letters <> "]+)_\\{([^}]*)\\}"] :> sub["$1", "$2"],
      RegularExpression["([" <> letters <> "]+)_([A-Za-z0-9\\x{0391}-\\x{03C9}]+)"] :> sub["$1", "$2"],
      RegularExpression["([" <> letters <> "])\\^\\{([^}]*)\\}"] :> sup["$1", "$2"],
      RegularExpression["([" <> letters <> "])\\^([A-Za-z0-9]+)"] :> sup["$1", "$2"]}] /.
    {sub[a_, b_] :> Cell[BoxData[FormBox[SubscriptBox[a, StyleBox[b, FontSlant -> "Plain"]], TraditionalForm]]],
     sup[a_, b_] :> Cell[BoxData[FormBox[SuperscriptBox[a, b], TraditionalForm]]]}];
statusStyle = <|"VERIFIED" -> {"\[Checkmark] VERIFIED", RGBColor[0.1, 0.5, 0.15], RGBColor[0.92, 0.97, 0.92]},
   "DISCREPANCY" -> {"\[WarningSign] DISCREPANCY", RGBColor[0.75, 0.1, 0.1], RGBColor[0.99, 0.92, 0.92]},
   "NOTE" -> {"\[FilledCircle] NOTE", RGBColor[0.15, 0.35, 0.75], RGBColor[0.92, 0.94, 0.99]}|>;
vrCell[vr_] := Module[{st = statusStyle[vr[[3]]]},
   Cell[TextData[Join[{StyleBox[st[[1]], FontWeight -> "Bold", FontColor -> st[[2]]], "    ", StyleBox[vr[[1]], FontWeight -> "Bold"], ":   "},
      textRuns[vr[[2]]], If[vr[[4]] === "", {}, Prepend[Replace[textRuns[vr[[4]]], s_String :> StyleBox[s, FontColor -> GrayLevel[0.3], FontSize -> 11], {1}], "\n"]]]],
    "Text", Background -> st[[3]], CellFrame -> {{4, 0}, {0, 0}}, CellFrameColor -> st[[2]], CellMargins -> {{66, 10}, {4, 4}}]];
outputCells[res_] := Which[
   res === Null, {},
   MatchQ[res, _Global`VR], {vrCell[res]},
   MatchQ[res, Column[_List, ___] | _List] && MemberQ[If[Head[res] === Column, First[res], res], _Global`VR],
   Flatten[If[MatchQ[#, _Global`VR], vrCell[#], Cell[BoxData[ToBoxes[#, StandardForm]], "Output"]] & /@ If[Head[res] === Column, First[res], res]],
   True, {Cell[BoxData[ToBoxes[res, StandardForm]], "Output"]}];

$msgLog = {}; $userPath = {"System`", "Global`"};
$progress = Environment["BUILD_PROGRESS"];
runCode[code_String] := Module[{prints = {}, ed, t0 = AbsoluteTime[], out},
   ed = Block[{Print = Function[Null, AppendTo[prints, Row[{##}]]; Null]},
      Block[{$Context = "Global`", $ContextPath = $userPath}, With[{r = EvaluationData[ToExpression[code, InputForm]]}, $userPath = $ContextPath; r]]];
   If[Lookup[ed, "MessagesText", {}] =!= {}, AppendTo[$msgLog, {StringTake[code, UpTo[90]], ed["MessagesText"]}]];
   If[StringQ[$progress], WriteString[$progress, "eval ", Round[AbsoluteTime[] - t0, 0.01], " s: ", StringTake[StringReplace[code, "\n" -> " "], UpTo[70]], "\n"]];
   out = Join[{Cell[parseBoxes[code], "Input"]},
    Cell[BoxData[ToBoxes[#, StandardForm]], "Print"] & /@ prints,
    Cell[#, "Message", "MSG", FontColor -> Red] & /@ Lookup[ed, "MessagesText", {}],
    outputCells[ed["Result"]]];
   If[StringQ[$progress], WriteString[$progress, "   boxes done ", Round[AbsoluteTime[] - t0, 0.01], " s\n"]];
   out];

UsingFrontEnd[
  Module[{t0 = AbsoluteTime[], flat, nb},
   flat = Flatten[Replace[cells, CodeCell[s_String] :> runCode[s], {1}]];
   nb = NotebookPut[Notebook[flat, WindowSize -> {1150, 950}, StyleDefinitions -> "Default.nb"]];
   If[StringQ[$progress], WriteString[$progress, "saving\n"]];
   NotebookSave[nb, outNB];
   If[StringQ[$progress], WriteString[$progress, "saved\n"]];
   If[outPDF =!= "none", Export[outPDF, nb]];
   Print["built ", Length[flat], " cells in ", Round[AbsoluteTime[] - t0, 0.1], " s"];
   If[ValueQ[Global`$verifications], Print["verification tally: ", Counts[Global`$verifications[[All, 3]]]];
    Scan[If[#[[3]] =!= "VERIFIED", Print["   ", #[[3]], "  ", #[[1]], ": ", #[[2]], "   ", #[[4]]]] &, Global`$verifications]];
   Print["cells with messages: ", Length[$msgLog]];
   Scan[Print["  MSG in: ", #[[1]], "\n     ", #[[2]]] &, $msgLog]]];
End[];
Exit[0];
