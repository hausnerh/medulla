// Overwrite events/<sample>/Livetime (the single-bin exposure TH1D medulla
// writes per sample) in place. Used by run_cc_sideband.sh --offbeam-exposure to
// fill the off-beam gate count that ReCAF2026 off-beam headers do not carry.
// Prints "SET_LIVETIME_OK <old> -> <new>" on success.
#include <cstdio>
#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TString.h"

void set_livetime(const char * fn, const char * sample, double value)
{
    TFile f(fn, "UPDATE");
    if(f.IsZombie()) { printf("SET_LIVETIME_FAIL cannot open %s\n", fn); return; }
    TDirectory * d = f.GetDirectory(Form("events/%s", sample));
    TH1D * h = d ? (TH1D *) d->Get("Livetime") : nullptr;
    if(!h) { printf("SET_LIVETIME_FAIL no events/%s/Livetime in %s\n", sample, fn); return; }
    double old = h->GetBinContent(1);
    h->SetBinContent(1, value);
    d->cd();
    h->Write("Livetime", TObject::kOverwrite);
    f.Close();
    printf("SET_LIVETIME_OK %.10g -> %.10g\n", old, value);
}
