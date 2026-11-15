"""
mz_calculator.py
================
Question this script answers:
    "If a metabolite is in my sample, at which m/z values will the mass
    spectrometer show it?"  and the reverse:
    "I see a peak at this m/z. Which metabolite could it be?"

Background in one paragraph:
    A mass spectrometer cannot see neutral molecules. It only sees ions.
    So each metabolite (neutral mass M) picks up or loses a small charged
    particle in the ion source. The result is called an "adduct".
    Example: glucose + a proton gives the ion written as [M+H]+.
    The instrument reports m/z = mass divided by charge. All ions here have
    a single charge, so m/z is simply the mass of the ion.

No extra packages are needed. Standard Python 3 is enough.

Run it with:
    python mz_calculator.py

Source of the method: https://codanics.com/metabolomics-complete-guide/
"""

# ---------------------------------------------------------------------------
# Masses of the small charged particles that attach to a metabolite (in Da).
# These are ion masses, so the mass of the missing electron is already
# taken into account.
# ---------------------------------------------------------------------------
PROTON = 1.007276   # H+   (a proton)
NA = 22.989218      # Na+  (sodium ion, common from glassware and buffers)
K = 38.963158       # K+   (potassium ion, common in biological samples)
NH4 = 18.033823     # NH4+ (ammonium ion, common when the buffer has ammonium)

# ---------------------------------------------------------------------------
# A tiny "database" of metabolites and their monoisotopic neutral masses (Da).
# "Monoisotopic" means the mass using the lightest, most common isotope of
# every atom (12C, 1H, 16O, 14N). Values are taken from HMDB.
# Add your own metabolites here to extend the calculator.
# ---------------------------------------------------------------------------
METABOLITES = {
    "glucose": 180.063388,
    "phenylalanine": 165.078979,
    "lactate (lactic acid)": 90.031694,
    "citrate (citric acid)": 192.027003,
    "glutamine": 146.069142,
}


def adducts(m):
    """
    Return the common single-charge ions for a neutral mass `m`.

    Positive mode (the molecule GAINS a positive particle):
        [M+H]+, [M+Na]+, [M+K]+, [M+NH4]+
    Negative mode (the molecule LOSES a proton):
        [M-H]-

    The result is a dictionary: ion name -> expected m/z.
    """
    return {
        "[M+H]+": m + PROTON,
        "[M+Na]+": m + NA,
        "[M+K]+": m + K,
        "[M+NH4]+": m + NH4,
        "[M-H]-": m - PROTON,
    }


def ppm_error(observed, theoretical):
    """
    Mass error in parts per million (ppm).

    ppm tells us how far the measured m/z is from the expected m/z,
    relative to the size of the number. Modern high-resolution instruments
    (Orbitrap, Q-TOF) are usually accurate to within 5 ppm.

    A positive value means the measured mass is too high,
    a negative value means it is too low.
    """
    return (observed - theoretical) / theoretical * 1e6


def match(observed_mz, tolerance_ppm=5):
    """
    Find every metabolite + adduct combination that could explain
    an observed m/z, within the given tolerance (default 5 ppm).

    Returns a list of tuples:
        (metabolite name, ion type, expected m/z, error in ppm)
    An empty list means nothing in our small database fits.

    Note: a mass match alone is only a weak identification (MSI level 3-4).
    Real identification also needs retention time and MS/MS spectra.
    """
    hits = []
    # Try every metabolite in the database...
    for name, mass in METABOLITES.items():
        # ...and every ion form that metabolite can take.
        for ion, mz in adducts(mass).items():
            err = ppm_error(observed_mz, mz)
            if abs(err) <= tolerance_ppm:
                hits.append((name, ion, round(mz, 4), round(err, 2)))
    return hits


if __name__ == "__main__":
    # Part 1: forward direction.
    # Print the expected m/z of every adduct for every metabolite.
    for name, mass in METABOLITES.items():
        ions = ", ".join(f"{ion} {mz:.4f}" for ion, mz in adducts(mass).items())
        print(f"{name:22s} M = {mass:.4f} -> {ions}")
    print()

    # Part 2: reverse direction.
    # Pretend we measured these three peaks and ask what they could be.
    for observed in [179.0560, 166.0864, 203.0530]:
        print(f"Observed m/z {observed}: {match(observed) or 'no match within 5 ppm'}")
