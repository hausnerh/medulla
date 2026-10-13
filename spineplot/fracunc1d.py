import re

import numpy as np
import matplotlib.pyplot as plt

from spectra1d import SpineSpectra1D
from systematic import Systematic
from utilities import mark_pot, mark_preliminary

class SpineFracUncertainty1D(SpineSpectra1D):
    """
    Per-bin fractional uncertainty of the stacked prediction, broken down
    by source. Each source is a regex over the names of the individual
    Systematic objects of the samples (the branch names of the multisim /
    multisigma / variation trees, `<sample>_statistical`, ...); the
    covariances of all matching systematics are summed (uncorrelated, as in
    Systematic.combine) over every non-ordinate sample, and the line drawn
    is sqrt(diag(cov)) / prediction. The black total line is the quadrature
    sum of the drawn sources. "_sigma" / "_nsigma" branches never match, as
    for systematic recipes.

    Histogram and systematic bookkeeping is inherited from SpineSpectra1D,
    so the prediction (and its exposure scaling) is exactly the stack drawn
    by the matching SpineSpectra1D figure.

    Attributes
    ----------
    _sources : list[dict]
        The configured sources. Each has a `pattern` (regex) and optionally
        a `label`, a `color` and `each` (true: one line per matching
        systematic instead of one line for their sum).
    _ordinate : str
        The name of the ordinate (data) sample, whose statistical
        uncertainty is not part of the prediction's uncertainty.
    """
    EXCLUDE = re.compile(r'(_nsigma|_sigma)$')

    def __init__(self, variable, categories, colors, category_types, sources,
                 title=None, xrange=None, xtitle=None, yrange=None,
                 ytitle=None, total_label='Total Uncertainty') -> None:
        """
        Initializes the SpineFracUncertainty1D.

        Parameters
        ----------
        variable : Variable
            The Variable object for the spectrum.
        categories : dict
            Map from the category enumeration to its label (as for
            SpineSpectra1D). Only 'histogram' categories enter the
            prediction.
        colors : dict
            Category colors (unused here; kept for the parent class).
        category_types : dict
            Map from the category label to 'histogram' or 'scatter'.
        sources : list[dict]
            The sources to draw (see the class docstring).
        total_label : str, optional
            Legend label of the total line. None disables the total.

        Returns
        -------
        None.
        """
        super().__init__(variable, categories, colors, category_types,
                         title, xrange, xtitle, yrange, ytitle)
        self._sources = sources
        self._total_label = total_label

    def _matching(self, pattern) -> dict:
        """
        Collect the systematics matching `pattern` across all non-ordinate
        samples, grouped by systematic name.

        Parameters
        ----------
        pattern : str
            The regex (re.search) to match against systematic names.

        Returns
        -------
        dict
            Systematic name -> list of Systematic objects (one per sample).
        """
        regxp = re.compile(pattern)
        out = dict()
        for sname, systs in self._systematics.items():
            if sname == self._ordinate:
                continue
            for k, s in systs.items():
                if self.EXCLUDE.search(k) or not regxp.search(k):
                    continue
                out.setdefault(k, []).append(s)
        return out

    def _variance(self, systs) -> np.ndarray:
        """
        Sum of the diagonal covariances of the given systematics for this
        artist's variable.
        """
        return np.sum([np.diag(s.get_covariance(self._variable._key)) for s in systs], axis=0)

    def draw(self, ax, style, logy=False, legend_ncol=2, **kwargs) -> None:
        """
        Draws the fractional uncertainty of each source and the total.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to draw the artist on.
        style : Style
            The style to use when drawing the artist.
        logy : bool, optional
            A flag to indicate if the y-axis should be logarithmic.
        legend_ncol : int, optional
            Number of legend columns. The default is 2.

        Returns
        -------
        None.
        """
        ax.set_xlabel(self._variable._xlabel if self._xtitle is None else self._xtitle)
        ax.set_ylabel('Fractional Uncertainty' if self._ytitle is None else self._ytitle)
        ax.set_xlim(*self._variable._range if self._xrange is None else self._xrange)
        ax.set_title(self._title)

        if self._plotdata is not None:
            labels = [l for l in self._plotdata if self._category_types[l] == 'histogram']
            prediction = np.sum([self._plotdata[l] for l in labels], axis=0)
            edges = self._binedges[labels[0]]
            def frac(var):
                with np.errstate(divide='ignore', invalid='ignore'):
                    return np.where(prediction > 0, np.sqrt(var) / prediction, np.nan)

            lines = list()
            for src in self._sources:
                matches = self._matching(src['pattern'])
                if not matches:
                    print(f"SpineFracUncertainty1D: no systematic matches '{src['pattern']}' — skipped.")
                    continue
                if src.get('each', False):
                    for k, systs in matches.items():
                        lines.append((k, None, self._variance(systs)))
                else:
                    systs = [s for v in matches.values() for s in v]
                    lines.append((src.get('label', src['pattern']), src.get('color', None), self._variance(systs)))

            # The total is drawn first (thick, underneath) so that a source
            # equal to the total in some bins (e.g. the only nonzero one)
            # stays visible on top of it.
            handles = list()
            if self._total_label is not None and lines:
                total = np.sum([var for _, _, var in lines], axis=0)
                total_handle = ax.stairs(frac(total), edges, label=self._total_label, color='black', linewidth=4)
            # Explicit colors: the plot styles override the default color
            # cycle (every line came out black).
            palette = plt.get_cmap('tab10' if len(lines) <= 10 else 'tab20').colors
            for i, (label, color, var) in enumerate(lines):
                handles.append(ax.stairs(frac(var), edges, label=label, linewidth=2,
                               color=color if color is not None else palette[i % len(palette)]))
            if self._total_label is not None and lines:
                handles.append(total_handle)
            ax.legend(handles=handles, ncol=legend_ncol)

        if isinstance(self._yrange, (tuple, list)):
            ax.set_ylim(*self._yrange)
        else:
            ax.set_ylim(0, None)
            if isinstance(self._yrange, (int, float)):
                ax.set_ylim(0, ax.get_ylim()[1] * self._yrange)
        if logy:
            ax.set_yscale('log')

        if style.mark_pot:
            mark_pot(ax, self._exposure, style.mark_pot_horizontal)
        if style.mark_preliminary is not None:
            mark_preliminary(ax, style.mark_preliminary)
