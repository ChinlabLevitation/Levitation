(* Export a saved notebook to PDF from a fresh front end (avoids dynamic-evaluation deadlocks in the build session).
   Usage: wolframscript -file exportpdf.wl in.nb out.pdf *)
{inNB, outPDF} = Take[$ScriptCommandLine, -2];
UsingFrontEnd[Module[{nb = NotebookOpen[inNB, Visible -> False]}, Export[outPDF, nb]; NotebookClose[nb]; Print["exported ", outPDF]]];
Exit[0];
