import ROOT
import math

def compute_s_over_sqrt_b(hist, fit_range=(60, 140), window_width=2):
    fit = ROOT.TF1("gaus_fit", "gaus", *fit_range)
    hist.Fit(fit, "RQ")

    mu = fit.GetParameter(1)
    sigma = fit.GetParameter(2)

    sig_min = mu - window_width * sigma
    sig_max = mu + window_width * sigma
    bin_sig_min = hist.FindBin(sig_min)
    bin_sig_max = hist.FindBin(sig_max)
    signal_total = hist.Integral(bin_sig_min, bin_sig_max)

    sb_left_min = mu - 2 * window_width * sigma
    sb_left_max = mu - window_width * sigma
    sb_right_min = mu + window_width * sigma
    sb_right_max = mu + 2 * window_width * sigma

    bin_sb_left_min = hist.FindBin(sb_left_min)
    bin_sb_left_max = hist.FindBin(sb_left_max)
    bin_sb_right_min = hist.FindBin(sb_right_min)
    bin_sb_right_max = hist.FindBin(sb_right_max)

    bkg_left = hist.Integral(bin_sb_left_min, bin_sb_left_max)
    bkg_right = hist.Integral(bin_sb_right_min, bin_sb_right_max)
    bkg_bins = (bin_sb_left_max - bin_sb_left_min) + (bin_sb_right_max - bin_sb_right_min)
    avg_bkg_per_bin = (bkg_left + bkg_right) / bkg_bins if bkg_bins > 0 else 0

    signal_bins = bin_sig_max - bin_sig_min
    bkg_estimate = avg_bkg_per_bin * signal_bins
    signal_estimate = signal_total - bkg_estimate
    snr = signal_estimate / math.sqrt(bkg_estimate) if bkg_estimate > 0 else float("inf")

    return {
        "mu": mu,
        "sigma": sigma,
        "S": signal_estimate,
        "B": bkg_estimate,
        "SNR": snr,
        "fit": fit
    }

def draw_histogram_with_fit(label, hist, results, output_dir="./"):
    c = ROOT.TCanvas(f"canvas_{label}", f"{label} Canvas", 800, 600)
    hist.SetLineWidth(2)
    hist.SetStats(0)  # Turn off ROOT's stat box
    hist.Draw("HIST")

    # Get the fit object and prepare it for drawing
    fit = results["fit"]
    fit.SetParameter(0, hist.GetMaximum())  # Visually match histogram peak height
    fit.SetLineColor(ROOT.kRed)
    fit.SetNpx(1000)  # Smooth curve
    fit.SetRange(hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())  # Full range
    fit.Draw("SAME")

    # Legend with TeX-style formatting
    legend = ROOT.TLegend(0.60, 0.65, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hist, f"{label} Histogram", "l")
    legend.AddEntry(fit, "Gaussian Fit", "l")
    legend.AddEntry(0, f"#mu = {results['mu']:.2f} GeV", "")
    legend.AddEntry(0, f"#sigma = {results['sigma']:.2f} GeV", "")
    legend.AddEntry(0, f"S/#sqrt{{B}} = {results['SNR']:.2f}", "")
    legend.Draw()

    c.SaveAs(f"{output_dir}{label}_fit.pdf")


# === MAIN SCRIPT ===
filename = "phide100_nobias.root"  # ← Change this as needed
f = ROOT.TFile.Open(filename)

histograms = {
    "Phide": f.Get("phide"),
    "Top": f.Get("top"),
    "W": f.Get("W"),
}

fit_ranges = {
    "Phide": (60, 140),
    "Top": (100, 250),
    "W": (50, 110),
}

print(f"\nS/√B Analysis for: {filename}")
print("-" * 50)

for label, hist in histograms.items():
    if not hist:
        print(f"{label}: Histogram not found!")
        continue

    results = compute_s_over_sqrt_b(hist, fit_range=fit_ranges[label])
    draw_histogram_with_fit(label, hist, results)
    print(f"{label} → μ = {results['mu']:.2f} GeV, σ = {results['sigma']:.2f} GeV, "
          f"S = {results['S']:.1f}, B = {results['B']:.1f}, S/√B = {results['SNR']:.2f}")
