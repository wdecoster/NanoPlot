import os
from base64 import b64encode
from io import BytesIO
from urllib.parse import quote as urlquote
import sys
import logging

import plotly.io as pio
try:
    pio.get_chrome()
except Exception as e:
    logging.warning(
        "Plotly could not fetch or find Chrome automatically. "
        "Static exports may fail unless BROWSER_PATH is set. Details: %s", e
    )

# DPI-aware writer
try:
    from nanoplot.utils import write_static_image
except Exception as e:
    logging.warning("Could not import write_static_image from nanoplot.utils: %s", e)
    write_static_image = None


class Plot(object):
    """A Plot object is defined by a path to the output file and the title of the plot."""

    only_report = False

    def __init__(self, path, title):
        self.path = path
        self.title = title
        self.fig = None      # Plotly fig for HTML/static; Matplotlib Figure in legacy mode
        self.html = None

    def encode(self):
        if self.html:
            return self.html
        elif self.fig:
            return self.encode2()
        else:
            return self.encode1()

    def encode1(self):
        data_uri = b64encode(open(self.path, "rb").read()).decode("utf-8").replace("\n", "")
        return '<img src="data:image/png;base64,{0}">'.format(data_uri)

    def encode2(self):
        # Legacy Matplotlib path only
        buf = BytesIO()
        self.fig.savefig(buf, format="png", bbox_inches="tight", dpi=100)
        buf.seek(0)
        string = b64encode(buf.read())
        return '<img src="data:image/png;base64,{0}">'.format(urlquote(string))

    def save(self, settings):
        if self.only_report:
            return

        if self.html:
            # Save the interactive HTML
            with open(self.path, "w") as html_out:
                html_out.write(self.html)

            # Also save static images unless suppressed
            if not settings.get("no_static", False):
                try:
                    fmts = settings.get("format", ["png"])
                    for fmt in fmts if isinstance(fmts, list) else [fmts]:
                        self.save_static(fmt, settings)  # pass settings
                except (AttributeError, ValueError) as e:
                    p = os.path.splitext(self.path)[0] + ".png"
                    if os.path.exists(p):
                        os.remove(p)
                    logging.warning("No static plots are saved due to an export problem:")
                    logging.warning(e)

        elif self.fig:
            # Legacy Matplotlib path
            fmts = settings.get("format", ["png"])
            dpi = int(settings.get("dpi", 300))
            if isinstance(fmts, list):
                for fmt in fmts:
                    self.fig.savefig(
                        fname=self.path + "." + fmt,
                        format=fmt,
                        bbox_inches="tight",
                        dpi=dpi,
                    )
            else:
                self.fig.savefig(
                    fname=self.path,
                    format=fmts,
                    bbox_inches="tight",
                    dpi=dpi,
                )
        else:
            sys.exit("No method to save plot object: no html or fig defined.")

    def show(self):
        if self.fig:
            return getattr(self.fig, "fig", self.fig)
        else:
            sys.stderr.write(".show not implemented for Plot instance without fig attribute!")

    def save_static(self, figformat, settings):
        """
        Export a Plotly figure as a static image with real DPI.
        Prefers utils.write_static_image; falls back to explicit pixel size.
        """
        output_path = self.path.replace(".html", f".{figformat}")
        dpi = int(settings.get("dpi", 300))

        if self.fig is None:
            logging.warning("No figure attached to Plot; skipping static export for %s", output_path)
            return

        # JSON just dumps the figure spec
        if figformat.lower() == "json":
            try:
                pio.write_json(self.fig, output_path)
                logging.info("Saved %s as JSON", output_path)
            except Exception as e:
                logging.warning("Failed to write JSON for %s: %s", output_path, e)
            return

        # Preferred path: DPI-aware helper
        try:
            if write_static_image is not None:
                write_static_image(self.fig, output_path, dpi=dpi)
                logging.info("Saved %s as %s (dpi=%d)", output_path, figformat, dpi)
                return
        except Exception as e:
            logging.warning("DPI helper failed for %s: %s; falling back to explicit px size", output_path, e)

        # Hard fallback so we don't end up at 700x500 defaults
        width_px  = int(6.4 * dpi)
        height_px = int(4.8 * dpi)
        try:
            pio.write_image(self.fig, output_path, width=width_px, height=height_px, scale=1)
            logging.info("Fallback saved %s at %dx%d px", output_path, width_px, height_px)
        except Exception as e:
            logging.warning("No static plots are saved due to an export problem:")
            logging.warning(e)
